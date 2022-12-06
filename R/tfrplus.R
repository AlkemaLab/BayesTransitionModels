#' Fit TFRplus model
#'
#' @export
tfrplus <- function(
    # Data and column names
  data,
  y,
  year,
  area,

  start_year = NA,
  end_year = NA,
  t_star = NULL,

  # Model settings
  model = "spline",
  num_knots = 7,
  spline_degree = 2,

  hierarchical_splines = c("intercept", area),

  # Priors
  tau_prior_sd = 1,
  rho_prior_sd = 1,

  # Out-of-sample validation
  held_out = FALSE,

  # Stan settings
  ...
) {
  args <- list(...)

  # Save original dataset
  original_data <- data

  ###### Initial argument checks #####
  stopifnot(is.numeric(tau_prior_sd))
  stopifnot(is.numeric(rho_prior_sd))
  stopifnot(is.numeric(spline_degree))
  stopifnot(is.numeric(num_knots))

  if(length(tau_prior_sd) > 1 || tau_prior_sd <= 0) {
    stop("tau_prior_sd must be a single positive number.")
  }

  if(length(rho_prior_sd) > 1 || rho_prior_sd <= 0) {
    stop("tau_prior_sd must be a single positive number.")
  }

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

  if(length(hierarchical_splines) == 0) {
    stop("No hierarchical structure supplied for the spline coefficients. See the hierarchical_splines argument.")
  }

  # Make sure there are no NAs in supplied columns
  check_nas(data, y)
  check_nas(data, year)

  # Initialize start and end year if necessary
  if(is.na(start_year)) start_year <- min(data[[year]])
  if(is.na(end_year)) end_year <- max(data[[year]])

  # Make sure the observed data are within the estimation period
  if(sum(!(data[[year]] %in% start_year:end_year)) > 0) {
    stop(glue::glue("Observations included in dataset that fall outside the estimation period ({start_year} to {end_year})."))
  }

  ###### Load model #####
  include_paths <- system.file("include", package = "BayesTransitionModels")
  stan_file_path <- system.file("stan/tfr_spline.stan", package = "BayesTransitionModels")
  if(model == "spline") {
    stan_model <- cmdstanr::cmdstan_model(
      stan_file_path,
      include_paths = include_paths
    )
  }
  else {
    stop(glue::glue("Model {model} not supported. Currently \"spline\" is the only supported model."))
  }

  #
  # Setup data for Stan
  #

  # Create district index for matching district and district index
  hierarchical_column_names <- unique(c(
    hierarchical_splines
  )) %>%
    setdiff("intercept")

  # Make sure there are no NAs in any of the columns
  for(column in hierarchical_column_names) {
    if(column == "intercept") next
    check_nas(data, column)
  }

  country_index <- data %>%
    dplyr::distinct(!!! syms(hierarchical_column_names)) %>%
    dplyr::mutate(c = 1:n())

  phase2 <- data %>%
    group_by(!!!sym(area)) %>%
    summarize(start_phase2 = min(period),
              end_phase2 = max(period)) %>%
    left_join(country_index)

  # Create year lookup table
  time_index <- tibble(
    period = seq(min(phase2$start_phase2), max(end_year, max(phase2$end_phase2)), by = 5)
  ) %>%
    mutate(t = 1:n())

  phases <- phase2 %>%
    left_join(time_index, by = c("start_phase2" = year)) %>%
    mutate(start_phase2 = t) %>%
    select(-t) %>%
    left_join(time_index, by = c("end_phase2" = year)) %>%
    mutate(end_phase2 = t) %>%
    select(-t) %>%
    arrange(c)

  year_by <- c()
  year_by[year] = year
  data <- data %>%
    dplyr::left_join(time_index, by = year_by) %>%
    dplyr::left_join(country_index, by = hierarchical_column_names)

  if(length(held_out) == 1 && held_out == FALSE) {
    held_out = rep(0, nrow(data))
  }
  else {
    held_out = as.numeric(held_out)
  }

  t_last <- max(data$t)

  # Set up hierarchical structures
  a_data       <- hierarchical_data(country_index, hierarchical_splines)

  # Set up spline basis
  knots <- sort(c(-100, seq(0, 1, length.out = num_knots)))
  grid <- c(-100, seq(from = 0, to = 1, by = .02)) # generating inputs

  B <- t(bs(grid, knots = knots, degree = spline_degree, intercept = FALSE))
  B <- B[1:(nrow(B) - 1), ]
  num_grid <- length(grid)
  num_basis <- nrow(B)
  ext_knots <- c(rep(knots[1], spline_degree), knots, rep(knots[length(knots)], spline_degree))

  a_lower_bound <- -0.01
  a_upper_bound <- -2.5

  stan_data <- list(
    C = nrow(country_index),
    T = nrow(time_index),
    N = nrow(data),
    held_out = held_out,

    time = array(data$t),
    country = array(data$c),

    start_phase2 = phases$start_phase2,
    end_phase2 = phases$end_phase2,

    y = array(data[[y]]),

    a_n_terms = a_data$n_terms,
    a_n_re = a_data$n_re,
    a_re_start = array(a_data$re_start),
    a_re_end = array(a_data$re_end),
    a_model_matrix = a_data$model_matrix$mat,

    # Spline settings
    num_knots = length(knots),
    knots = knots,

    num_grid = num_grid,
    spline_degree = spline_degree,
    grid = grid,
    B = B,

    a_lower_bound = a_lower_bound,
    a_upper_bound = a_upper_bound,


    # AR priors
    rho_prior_sd = rho_prior_sd,
    tau_prior_sd = tau_prior_sd
  )

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
                 phases = phases,

                 # Save arguments
                 y = y,
                 year = year,
                 source = source,
                 area = area,
                 hierarchical_splines = hierarchical_splines,
                 held_out = held_out,
                 model = model)

  cat("Extracting posteriors...\n")

  result$posteriors <- process_tfr_fit(result, ifelse(is.null(args$parallel_chains), 1, args$parallel_chains))

  attr(result, "class") <- "fpemplus"

  result
}
