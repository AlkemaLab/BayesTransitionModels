#' @import dplyr
#' @importFrom tidybayes median_qi spread_draws
#' @importFrom stats plogis
#' @importFrom rlang .data
extract_rate_vs_level_subhierarchical <- function(fit, f, subhierarchy, constrain_zero = "after") {
  hierarchical_a <- hierarchical_data(fit$country_index, f)

  start <- hierarchical_a$model_matrix$index %>%
    dplyr::filter(.data$column == subhierarchy) %>%
    dplyr::pull(i) %>%
    min()

  end <- hierarchical_a$model_matrix$index %>%
    dplyr::filter(.data$column == subhierarchy) %>%
    dplyr::pull(i) %>%
    max()

  B <- fit$stan_data$B

  num_constrained_zero <- fit$stan_data$spline_degree + 1

  a_star <- fit$samples$draws(c("a_star")) %>%
    tidybayes::spread_draws(a_star[j, i]) %>%
    dplyr::group_by(.data$i, .data$.chain, .data$.iteration, .data$.draw) %>%
    tidyr::nest() %>%
    dplyr::mutate(a_star = map(.data$data, `[[`, "a_star")) %>%
    dplyr::select(-.data$data)

  uniq <- unique(hierarchical_a$model_matrix$mat[, 1:end, drop = FALSE])

  titles <- c()
  for(i in 1:nrow(uniq)) {
    index <- rep(0, hierarchical_a$n_terms)
    index[1:end] <- uniq[i, 1:end]
    title <- hierarchical_a$model_matrix$index %>%
      dplyr::filter(i == last(which(index == 1))) %>%
      dplyr::pull(.data$level)
    titles <- c(titles, title)

    a_star[[title]] = map_dbl(a_star$a_star, function(a_star) {
      fit$stan_data$a_lower_bound + (fit$stan_data$a_upper_bound - fit$stan_data$a_lower_bound) * stats::plogis(index %*% a_star)
    })
  }
  a_star <- a_star %>%
    tidyr::pivot_longer(cols = all_of(titles)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$a_star) %>%
    dplyr::group_by(.data$.chain, .data$.iteration, .data$.draw, .data$name) %>%
    tidyr::nest() %>%
    dplyr::mutate(Y = map(.data$data, function(data) {
      if(constrain_zero == "after") {
        tibble::tibble(x = fit$stan_data$grid, Y = t(c(data$value, rep(0, num_constrained_zero)) %*% B)[,1]) %>%
          dplyr::filter(.data$x >= 0 & .data$x <= 1)
      }
      else {
        tibble::tibble(x = fit$stan_data$grid, Y = t(c(rep(0, num_constrained_zero), data$value) %*% B)[,1]) %>%
          dplyr::filter(.data$x >= 0 & .data$x <= 1)
      }
    }))  %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$data, -.data$.draw) %>%
    tidyr::unnest(.data$Y)

  return(a_star)
}

#' @importFrom rlang .data
extract_parameter_subhierarchical <- function(fit, f, subhierarchy, x) {
  hierarchy <- hierarchical_data(fit$country_index, f)

  start <- hierarchy$model_matrix$index %>%
    dplyr::filter(.data$column == subhierarchy) %>%
    dplyr::pull(i) %>%
    min()

  end <- hierarchy$model_matrix$index %>%
    dplyr::filter(.data$column == subhierarchy) %>%
    dplyr::pull(i) %>%
    max()

  pars <- c(glue::glue("{x}_star"))

  star <- fit$samples$draws(pars) %>%
    tidybayes::spread_draws((!!sym(pars[1]))[i]) %>%
    dplyr::group_by(.data$.chain, .data$.iteration, .data$.draw) %>%
    tidyr::nest() %>%
    dplyr::mutate(star = map(.data$data, `[[`, glue::glue("{x}_star"))) %>%
    dplyr::select(-.data$data)

  uniq <- unique(hierarchy$model_matrix$mat[, 1:end, drop = FALSE])

  titles <- c()
  for(i in 1:nrow(uniq)) {
    index <- rep(0, hierarchy$n_terms)
    index[1:end] <- uniq[i, 1:end]
    title <- hierarchy$model_matrix$index %>%
      dplyr::filter(.data$i == last(which(index == 1))) %>%
      dplyr::pull(.data$level)
    titles <- c(titles, title)

    star[[title]] = map_dbl(star[["star"]], function(star) {
      index %*% star
    })
  }

  star <- star %>%
    tidyr::pivot_longer(cols = all_of(titles)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$star)

  return(star)
}

#' @import dplyr tidyr tidybayes stringr
#' @importFrom stats quantile
process_fit <- function(fit, parallel_chains = NULL) {
  if(is.null(parallel_chains)) parallel_chains <- 1

  #
  # eta and epsilon summaries
  #
  if(fit$stan_data$smoothing == 1) {
    temporal_variables <- c("eta", "epsilon")
  }
  else {
    temporal_variables <- c("eta")
  }

  temporal <- fit$samples$summary(temporal_variables, ~stats::quantile(.x, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)), .cores = parallel_chains) %>%
    tidyr::separate(.data$variable, c("variable", "index"), "\\[") %>%
    dplyr::mutate(index = stringr::str_replace_all(.data$index, "\\]", "")) %>%
    tidyr::separate(.data$index, c("c", "t"), ",") %>%
    dplyr::mutate_at(vars(c, t), as.integer) %>%
    dplyr::left_join(fit$country_index, by = "c") %>%
    dplyr::left_join(fit$time_index, by = "t")

  #
  # Transition function summaries
  #
  transition <- list()
  transition_quantiles <- list()
  for(column in fit$hierarchical_splines) {
    transition[[column]] <- extract_rate_vs_level_subhierarchical(fit, fit$hierarchical_splines, column)
    transition_quantiles[[column]] <- transition[[column]] %>%
      dplyr::group_by(.data$name, .data$x) %>%
      tidybayes::median_qi(.data$Y, .width = c(0.5, 0.8, 0.95))
  }


  #
  # Hyperparameters
  #
  P_tilde <- list()
  for(column in fit$hierarchical_asymptote) {
    P_tilde[[column]] <- extract_parameter_subhierarchical(fit, fit$hierarchical_asymptote, column, "P_tilde")
  }

  Omega <- list()
  for(column in fit$hierarchical_level) {
    Omega[[column]] <- extract_parameter_subhierarchical(fit, fit$hierarchical_level, column, "Omega")
  }

  ar <- NULL
  if(fit$stan_data$smoothing == 1) {
    ar <- fit$samples$draws(c("est_rho", "est_tau")) %>%
      tidybayes::spread_draws(est_rho[i], est_tau[i]) %>%
      ungroup() %>%
      select(-.data$i)
  }

  #
  # Hierarchical distributions
  #
  P_tilde_sigma <- fit$samples$draws(c("P_tilde_sigma")) %>%
    tidybayes::spread_draws(P_tilde_sigma[i])
  Omega_sigma <- fit$samples$draws(c("Omega_sigma")) %>%
    tidybayes::spread_draws(Omega_sigma[i])
  a_sigma <- fit$samples$draws(c("a_sigma")) %>%
    tidybayes::spread_draws(a_sigma[i, j])

  #
  # Data model
  #
  nonse <- fit$samples$draws("nonse") %>%
    tidybayes::spread_draws(nonse[source]) %>%
    dplyr::left_join(fit$source_index, by = "source")

  #
  # Generated quantities
  #
  generated_quantities <- fit$samples$draws(c("pit", "resid", "y_pred")) %>%
    tidybayes::spread_draws(pit[i], resid[i], y_pred[i])

  ans <- list(
    temporal = temporal,
    transition = transition,
    transition_quantiles = transition_quantiles,
    ar = ar,
    nonse = nonse,
    P_tilde = P_tilde,
    Omega = Omega,
    generated_quantities = generated_quantities,
    P_tilde_sigma = P_tilde_sigma,
    Omega_sigma = Omega_sigma,
    a_sigma = a_sigma
  )

  ans
}

rate_logistic_logit <- function(P, P_tilde, omega) {
  if(P == 1 && P_tilde == 1) return(0)

  ifelse(P < P_tilde, (P - P_tilde) * omega / (P_tilde * (P - 1)), 0)
}

#' @import dplyr tidyr tidybayes stringr
#' @importFrom stats quantile
process_logistic_fit <- function(fit, parallel_chains = NULL) {
  if(is.null(parallel_chains)) parallel_chains <- 1

  #
  # eta and epsilon summaries
  #
  if(fit$stan_data$smoothing == 1) {
    temporal_variables <- c("eta", "epsilon")
  }
  else {
    temporal_variables <- c("eta")
  }

  temporal <- fit$samples$summary(temporal_variables, ~stats::quantile(.x, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)), .cores = parallel_chains) %>%
    tidyr::separate(.data$variable, c("variable", "index"), "\\[") %>%
    dplyr::mutate(index = stringr::str_replace_all(.data$index, "\\]", "")) %>%
    tidyr::separate(.data$index, c("c", "t"), ",") %>%
    dplyr::mutate_at(vars(c, t), as.integer) %>%
    dplyr::left_join(fit$country_index, by = "c") %>%
    dplyr::left_join(fit$time_index, by = "t")

  #
  # Hyperparameters
  #
  P_tilde <- list()
  for(column in fit$hierarchical_asymptote) {
    P_tilde[[column]] <- extract_parameter_subhierarchical(fit, fit$hierarchical_asymptote, column, "P_tilde")
  }

  omega <- list()
  for(column in fit$hierarchical_rate) {
    omega[[column]] <- extract_parameter_subhierarchical(fit, fit$hierarchical_rate, column, "omega")
  }

  num_constrained_zero <- fit$stan_data$spline_degree + 1
  transition <- list()
  transition_quantiles <- list()
  for(column in intersect(fit$hierarchical_asymptote, fit$hierarchical_rate)) {
    transition[[column]] <- omega[[column]] %>%
      dplyr::rename(omega = .data$value) %>%
      dplyr::left_join(P_tilde[[column]] %>%
      dplyr::rename(P_tilde = .data$value)) %>%
      dplyr::mutate(P_tilde = 0.5 + 0.45 * plogis(P_tilde), omega = 0.5 * plogis(omega)) %>%
      dplyr::mutate(Y = map2(omega, P_tilde, function(omega, P_tilde) {
        coefs <- map_dbl(fit$spline_maxima * P_tilde, rate_logistic_logit, P_tilde = P_tilde, omega = omega)
        coefs <- coefs[1:(length(coefs) - num_constrained_zero + 1)]
        tibble::tibble(x = fit$stan_data$grid, Y = t(c(coefs, rep(0, num_constrained_zero - 1)) %*% fit$stan_data$B)[,1])
      })) %>%
      unnest(.data$Y)

    transition_quantiles[[column]] <- transition[[column]] %>%
      dplyr::group_by(.data$name, .data$x) %>%
      tidybayes::median_qi(.data$Y, .width = c(0.5, 0.8, 0.95))
  }

  Omega <- list()
  for(column in fit$hierarchical_level) {
    Omega[[column]] <- extract_parameter_subhierarchical(fit, fit$hierarchical_level, column, "Omega")
  }

  ar <- NULL
  if(fit$stan_data$smoothing == 1) {
    ar <- fit$samples$draws(c("est_rho", "est_tau")) %>%
      tidybayes::spread_draws(est_rho[i], est_tau[i]) %>%
      ungroup() %>%
      select(-.data$i)
  }

  #
  # Hierarchical distributions
  #
  P_tilde_sigma <- fit$samples$draws(c("P_tilde_sigma")) %>%
    tidybayes::spread_draws(P_tilde_sigma[i])
  omega_sigma <- fit$samples$draws(c("omega_sigma")) %>%
    tidybayes::spread_draws(omega_sigma[i])
  Omega_sigma <- fit$samples$draws(c("Omega_sigma")) %>%
    tidybayes::spread_draws(Omega_sigma[i])

  #
  # Data model
  #
  nonse <- fit$samples$draws("nonse") %>%
    tidybayes::spread_draws(nonse[source]) %>%
    dplyr::left_join(fit$source_index, by = "source")

  #
  # Generated quantities
  #
  generated_quantities <- fit$samples$draws(c("pit", "resid", "y_pred")) %>%
    tidybayes::spread_draws(pit[i], resid[i], y_pred[i])

  ans <- list(
    temporal = temporal,
    transition = transition,
    transition_quantiles = transition_quantiles,
    ar = ar,
    nonse = nonse,
    omega = omega,
    P_tilde = P_tilde,
    Omega = Omega,
    generated_quantities = generated_quantities,
    P_tilde_sigma = P_tilde_sigma,
    Omega_sigma = Omega_sigma,
    omega_sigma = omega_sigma
  )

  ans
}

#' @import dplyr tidyr tidybayes stringr
#' @importFrom stats quantile
process_tfr_fit <- function(fit, parallel_chains = NULL) {
  if(is.null(parallel_chains)) parallel_chains <- 1

  #
  # eta and epsilon summaries
  #
  temporal_variables <- c("eta")

  temporal <- fit$samples$summary(temporal_variables, ~stats::quantile(.x, probs = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)), .cores = parallel_chains) %>%
    tidyr::separate(.data$variable, c("variable", "index"), "\\[") %>%
    dplyr::mutate(index = stringr::str_replace_all(.data$index, "\\]", "")) %>%
    tidyr::separate(.data$index, c("c", "t"), ",") %>%
    dplyr::mutate_at(vars(c, t), as.integer) %>%
    dplyr::left_join(fit$country_index, by = "c") %>%
    dplyr::left_join(fit$time_index, by = "t") %>%
    dplyr::left_join(dplyr::select(fit$phases, .data$start_phase2, .data$end_phase2, .data$c)) %>%
    dplyr::filter(t >= .data$start_phase2)

  #
  # Transition function summaries
  #
  transition <- list()
  transition_quantiles <- list()
  for(column in fit$hierarchical_splines) {
    transition[[column]] <- extract_rate_vs_level_subhierarchical(fit, fit$hierarchical_splines, column, "before")
    transition_quantiles[[column]] <- transition[[column]] %>%
      dplyr::group_by(.data$name, .data$x) %>%
      tidybayes::median_qi(.data$Y, .width = c(0.5, 0.8, 0.95))
  }


  #
  # Hierarchical distributions
  #
  a_sigma <- fit$samples$draws(c("a_sigma")) %>%
    tidybayes::spread_draws(a_sigma[i, j])

  #
  # Data model
  #
  sigma <- fit$samples$draws("sigma")

  #
  # Generated quantities
  #
  generated_quantities <- fit$samples$draws(c("pit", "resid", "y_pred")) %>%
    tidybayes::spread_draws(pit[i], resid[i], y_pred[i])

  ans <- list(
    temporal = temporal,
    transition = transition,
    transition_quantiles = transition_quantiles,
    sigma = sigma,
    generated_quantities = generated_quantities,
    a_sigma = a_sigma
  )

  ans
}
