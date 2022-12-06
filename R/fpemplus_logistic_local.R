
#' Fit the fpemplus model to data.
#'
#' @param data a data frame.
#' @param start_year start year of estimates.
#' @param end_year end year of estimates.
#' @param y column name of outcome.
#' @param se column name of outcome standard error.
#' @param year column name of outcome year.
#' @param source column name of data source.
#' @param area column name of the area of each observation
#' @param t_star reference year.
#' @param model which model to fit. Currently only "spline" is supported.
#' @param num_knots number of spline knots.
#' @param spline_degree spline degree. Degree 2 or 3 is supported.
#' @param smoothing whether to include smoothing component.
#' @param hierarchical_asymptote vector specifying hierarchical structure for asymptote (see Details).
#' @param hierarchical_level vector specifying hierarchical structure for the level in reference year (see Details).
#' @param hierarchical_splines vector specifying hierarchical structure for spline coefficients (see Details).
#' @param tau_prior Stan prior definition for the AR tau parameter
#' @param rho_prior Stan prior definition for the AR rho parameter
#' @param held_out binary vector indicating which observations are held out. Set to FALSE to hold out no observations.
#' @param ... additional arguments for CmdStanModel::sample.
#'
#' @return fpemplus object.
#'
#' @details
#' Several area-specific parameters of the fpemplus model have hierarchical priors
#' assigned to them so that information can be shared between areas.
#' The package allows the structure of the hierarchical prior to be configured by the user
#' through the \code{hierarchical_asymptote}, \code{hierarchical_level}, and \code{hierarchical_splines} arguments.
#' These arguments expect a character vector that specifies a nesting hierarchical structure.
#' Each element of the vector must be either "intercept" or a column name in the dataset, where
#' "intercept" will add a global intercept for the parameter.
#' The vector must be in descending order in terms of the hierarchy: that is, it starts with
#' "intercept" and proceeds down the hierarchy.
#'
#' For example, suppose we are fitting country-level data, where the dataset has columns
#' "name_country", "name_sub_region", and "name_region" containing the name of the country,
#' sub-region, and region that each observation belongs to. To specify that the spline coefficients
#' should be fitted with a hierarchical model in which countries are nested within sub-regions within regions within world,
#' we would use the argument
#' \code{hierarchical_splines = c("intercept", "name_region", "name_sub_region", "name_country")}.
#'
#' @importFrom cmdstanr cmdstan_model write_stan_file
#' @importFrom tibble tibble
#' @importFrom splines bs
#' @import dplyr
#' @importFrom readr read_file
#' @importFrom stringr str_replace_all
#'
#' @export
fpemplus_logistic_local <- function(
  # Data and column names
  data,
  y,
  se,
  area,
  year,
  source,

  full_fit,

  start_year = NA,
  end_year = NA,
  t_star = NULL,

  # Model settings
  model = "spline",
  num_knots = 7,
  spline_degree = 2,
  smoothing = TRUE,

  extra_stan_data = list(),

  # Out-of-sample validation
  held_out = FALSE,

  # Stan settings
  ...
) {
  args <- list(...)

  # Save original dataset
  original_data <- data

  ###### Initial argument checks #####
  stopifnot(is.numeric(spline_degree))
  stopifnot(is.numeric(num_knots))
  stopifnot(is.logical(smoothing))

  if(!(spline_degree %in% c(2, 3))) {
    stop("spline_degree must be either 2 or 3.")
  }

  if(num_knots <= 0) {
    stop("num_knots must be greater than zero.")
  }

  if(nrow(data) == 0) {
    stop("Data has no rows.")
  }

  if(length(held_out) > 1) {
    if(length(held_out) != nrow(data)) stop(glue::glue("held_out (length {length(held_out)}) must be same size as dataset ({nrow(data)} rows)."))
  }

  if(start_year > end_year) {
    stop("start_year must be less than end year")
  }

  # Make sure there are no NAs in supplied columns
  check_nas(data, y)
  check_nas(data, se)
  check_nas(data, year)
  check_nas(data, area)
  check_nas(data, source)

  # Make sure SEs are all positive and non-zero
  if(sum(data[[se]] <= 0) > 0) {
    stop(glue::glue("All standard errors must be greater than zero ({sum(data[[se]] <= 0)} observations supplied with zero or negative SEs)."))
  }

  # Initialize start and end year if necessary
  if(is.na(start_year)) start_year <- min(data[[year]])
  if(is.na(end_year)) end_year <- max(data[[year]])

  # Make sure the observed data are within the estimation period
  if(sum(!(data[[year]] %in% start_year:end_year)) > 0) {
    stop(glue::glue("Observations included in dataset that fall outside the estimation period ({start_year} to {end_year})."))
  }

  ###### Load model #####
  include_paths <- system.file("include", package = "BayesTransitionModels")
  stan_file_path <- system.file("stan/fpem_spline_logistic_local.stan", package = "BayesTransitionModels")
  if(model == "spline") {
    stan_model <- cmdstanr::cmdstan_model(
      stan_file_path,
      include_paths = include_paths
    )
  }
  else if(file.exists(model)) {
    stan_model <- cmdstanr::cmdstan_model(
      model,
      include_paths = include_paths
    )
  }
  else {
    stop(glue::glue("Model {model} not supported. Currently \"spline\" is the only supported model."))
  }

  #
  # Setup data for Stan
  #

  # Create year lookup table
  time_index <- tibble::tibble(
    year = seq(start_year, end_year),
    t = 1:length(year)
  )

  area_name <- data[[area]][1]
  country_index = data %>%
    dplyr::distinct(!!! syms(area)) %>%
    dplyr::mutate(c = 1:n())

  source_index <- data %>%
    dplyr::distinct("{source}" := .data[[source]]) %>%
    dplyr::mutate(source = 1:n())

  year_by <- c()
  year_by[year] = "year"
  data <- data %>%
    dplyr::left_join(time_index, by = year_by) %>%
    dplyr::left_join(source_index, by = source)

  if(length(held_out) == 1 && held_out == FALSE) {
    held_out = rep(0, nrow(data))
  }
  else {
    held_out = as.numeric(held_out)
  }

  if(!is.null(t_star)) {
    if(!(t_star %in% time_index$year)) stop(glue::glue("Supplied t_star {t_star} not in range of estimation years ({start_year} to {end_year})."))
    t_stars <- t_star
  }
  else {
    # Otherwise, pick t_star to be the mean of the observation years.
    t_stars <- data %>%
      dplyr::summarize(t_star = round(mean(t))) %>%
      dplyr::pull(t_star)
  }

  t_last <- max(data$t)

  # Set up spline basis
  knots <- c(seq(0, 1, length.out = num_knots), 1000)
  grid <- c(seq(from = 0, to = 1, by = .02), 1000) # generating inputs
  B <- t(splines::bs(grid, knots = knots, degree = spline_degree, intercept = FALSE))
  B <- B[1:(nrow(B) - 1), ]
  num_grid <- length(grid)
  num_basis <- nrow(B)
  ext_knots <- c(rep(knots[1], spline_degree), knots, rep(knots[length(knots)], spline_degree))
  spline_maxima <- find_knot_maxima(knots, num_knots, spline_degree, c(seq(0, 1, 0.001), 1000, 1001))

  a_lower_bound <- 0.01
  a_upper_bound <- 0.3

  #
  # Fixed parameters
  #
  rho <- mean(full_fit$posteriors$ar$est_rho)
  tau <- mean(full_fit$posteriors$ar$est_tau)

  hierarchical_level     <- full_fit$hierarchical_level[length(full_fit$hierarchical_level) - 1]
  hierarchical_asymptote <- full_fit$hierarchical_asymptote[length(full_fit$hierarchical_asymptote) - 1]
  hierarchical_splines   <- full_fit$hierarchical_splines[length(full_fit$hierarchical_splines) - 1]

  Omega_mu    <- full_fit$posteriors$Omega[[hierarchical_level]]$value %>%
    mean()

  Omega_sigma <- full_fit$posteriors$Omega_sigma %>%
    dplyr::ungroup() %>%
    dplyr::filter(i == max(.data$i)) %>%
    dplyr::pull(.data$Omega_sigma) %>%
    mean()

  omega_mu    <- full_fit$posteriors$omega[[hierarchical_level]]$value %>%
    mean()

  omega_sigma <- full_fit$posteriors$omega_sigma %>%
    dplyr::ungroup() %>%
    dplyr::filter(i == max(.data$i)) %>%
    dplyr::pull(.data$omega_sigma) %>%
    mean()

  P_tilde_mu <- full_fit$posteriors$P_tilde[[hierarchical_asymptote]]$value %>%
    mean()

  P_tilde_sigma <- full_fit$posteriors$P_tilde_sigma %>%
    dplyr::ungroup() %>%
    dplyr::filter(i == max(.data$i)) %>%
    dplyr::pull(.data$P_tilde_sigma) %>%
    mean()

  nonse <- source_index %>%
    dplyr::left_join(
      full_fit$posteriors$nonse %>%
        dplyr::group_by(!!sym(fit$source)) %>%
        dplyr::summarize(mean = mean(.data$nonse))
    ) %>%
    pull(.data$mean)

  stan_data <- list(
    T = nrow(time_index),
    N = nrow(data),
    S = nrow(source_index),
    held_out = held_out,

    time = array(data$t),
    source = array(data$source),

    s = array(data[[se]]),
    y = array(data[[y]]),

    t_star = t_stars,
    t_last = t_last,

    # Spline settings
    num_knots = length(knots),
    knots = knots,

    num_grid = num_grid,
    spline_degree = spline_degree,
    grid = grid,
    B = B,
    spline_maxima = spline_maxima,

    a_lower_bound = a_lower_bound,
    a_upper_bound = a_upper_bound,

    Omega_mu = Omega_mu,
    Omega_sigma = Omega_sigma,

    P_tilde_mu = P_tilde_mu,
    P_tilde_sigma = P_tilde_sigma,

    omega_mu = omega_mu,
    omega_sigma = omega_sigma,

    rho = rho,
    tau = tau,

    nonse = nonse,

    # Switches
    smoothing = as.numeric(smoothing)
  )

  stan_data <- c(extra_stan_data, stan_data)

  fit <- stan_model$sample(
    stan_data,
    save_latent_dynamics = TRUE,
    ...
  )

  result <- list(samples = fit,
                 data = original_data,
                 stan_data = stan_data,
                 time_index = time_index,
                 country_index = country_index,
                 source_index = source_index,

                 # Save arguments
                 y = y,
                 se = se,
                 year = year,
                 source = source,
                 area = area,
                 held_out = held_out,
                 model = model)

  cat("Extracting posteriors...\n")

  result$posteriors <- process_logistic_fit_local(result, ifelse(is.null(args$parallel_chains), 1, args$parallel_chains))

  attr(result, "class") <- "fpemplus"

  result
}
