plot_temporal <- function(x, fit, areas = c(), plot_data = FALSE, color_sources = FALSE) {
  if(length(areas) == 0) areas <- fit$country_index[[fit$area]]

  post <- fit$posteriors$temporal %>%
    filter(!!sym(fit$area) %in% areas) %>%
    filter(.data$variable == x)


  p <- ggplot2::ggplot(post, aes_string(x = fit$year, y = "`50%`")) +
    ggplot2::geom_ribbon(aes(ymin = .data$`2.5%`, ymax = .data$`97.5%`, fill = "95%")) +
    ggplot2::geom_ribbon(aes(ymin = .data$`10%`,  ymax = .data$`90%`, fill = "80%")) +
    ggplot2::geom_ribbon(aes(ymin = .data$`25%`,  ymax = .data$`75%`, fill = "50%")) +
    ggplot2::geom_line() +
    ggplot2::scale_fill_brewer(direction = -1) +
    ggplot2::facet_wrap(vars(!!sym(fit$area))) +
    labs(fill = "Posterior\nQuantile", x = fit$year)

  if(plot_data == TRUE) {
    data_tibble <- tibble::tibble(i = 1:length(fit$held_out), held_out = as.logical(fit$held_out))
    data <- fit$data %>%
      dplyr::mutate(i = 1:n()) %>%
      dplyr::left_join(data_tibble, by = "i")

    data[[fit$source]] <- factor(data[[fit$source]])

    filtered_data <- data %>%
      dplyr::filter(!!sym(fit$area) %in% areas)

    if(!is.null(fit$se)) {
      filtered_data <- filtered_data %>% mutate(
        lower = truncnorm::qtruncnorm(0.025, mean = !!sym(fit$y), sd = !!sym(fit$se), a = 0, b = 1),
        upper = truncnorm::qtruncnorm(0.975, mean = !!sym(fit$y), sd = !!sym(fit$se), a = 0, b = 1)
      )
    }

    some_held_out <- any(fit$held_out == 1)

    if(color_sources == TRUE) {
      if(some_held_out == TRUE) {
        point_aes <- aes_string(y = fit$y, color = fit$source, shape = "held_out")
      }
      else {
        point_aes <- aes_string(y = fit$y, color = fit$source)
      }
      if(!is.null(fit$se)) {
        p <- p + ggplot2::geom_errorbar(aes_string(y = fit$y, ymin = "lower", ymax = "upper", color = fit$source), alpha = 0.3, width = 0, data = filtered_data)
      }
      p <- p + ggplot2::geom_point(point_aes, data = filtered_data, alpha = 0.7)
    }
    else {
      if(some_held_out == TRUE) {
        point_aes <- aes_string(y = fit$y, shape = "held_out")
      }
      else {
        point_aes <- aes_string(y = fit$y)
      }
      if(!is.null(fit$se)) {
        p <- p + ggplot2::geom_errorbar(aes_string(y = fit$y, ymin = "lower", ymax = "upper"), alpha = 0.3, width = 0, data = filtered_data)
      }
      p <- p + ggplot2::geom_point(point_aes, data = filtered_data, alpha = 0.7)
    }
  }

  p
}

#' Plot posterior of the indicator.
#'
#' Plots posterior quantiles of the latent indicator (\eqn{eta_{c,t}}).
#'
#' @param fit fpemplus object
#' @param areas which areas to include in plot.
#' @param color_sources whether to color each data point by data source type
#'
#' @import tidybayes dplyr ggdist
#' @importFrom truncnorm qtruncnorm
#'
#' @export
plot_indicator <- function(fit, areas = c(), color_sources = FALSE) {
  check_is_fpemplus_class(fit)
  plot_temporal("eta", fit, areas, plot_data = TRUE, color_sources = color_sources) +
    labs(y = expression(eta[ct]))
}

#' Plot posterior of the smoothing component.
#'
#' Plots posterior quantiles of the smoothing component (\eqn{epsilon_{c,t}}).
#'
#' @param fit fpemplus object
#' @param areas which areas to include in plot.
#'
#' @import tidybayes dplyr ggdist
#' @importFrom truncnorm qtruncnorm
#'
#' @export
plot_smoother <- function(fit, areas = c()) {
  check_is_fpemplus_class(fit)

  if(fit$stan_data$smoothing == 0) {
    stop("Cannot plot smoothing component because the fpemplus model was fit with smoothing disabled.")
  }

  plot_temporal("epsilon", fit, areas) +
    labs(y = expression(epsilon[ct]))
}

#' Plot transition functions for a level of hierarchy.
#'
#' @param fit fpemplus object
#' @param hierarchical_level which level of hierarchy to plot
#' @param areas which areas to include in the plot
#' @param n how many draws of the transition function to plot
#'
#' @import ggplot2 dplyr
#'
#' @export
plot_transition <- function(fit, hierarchical_level = NA, areas = c(), n = 0) {
  check_is_fpemplus_class(fit)

  if(is.na(hierarchical_level)) {
    hierarchical_level <- names(fit$posteriors$transition)[1]
  }

  check_areas(fit, hierarchical_level, areas)

  data <- fit$posteriors$transition[[hierarchical_level]]
  data_quantiles <- fit$posteriors$transition_quantiles[[hierarchical_level]]
  if(length(areas) > 0) {
    data <- dplyr::filter(data, .data$name %in% areas) %>%
      filter(.data$x < 1000)
    data_quantiles <- dplyr::filter(data_quantiles, .data$name %in% areas) %>%
      filter(.data$x < 1000)
  }

  p <- ggplot2::ggplot(data_quantiles, aes(x = .data$x, y = .data$Y)) +
    ggdist::geom_lineribbon(aes(ymin = .data$.lower, ymax = .data$.upper)) +
    ggplot2::scale_fill_brewer() +
    ggplot2::facet_wrap(~name) +
    labs(x = "x", y = "y")

  if(n > 0) {
    draws <- fit$posteriors$transition[[hierarchical_level]] %>%
      dplyr::distinct(.data$.chain, .data$.iteration) %>%
      dplyr::sample_n(n) %>%
      dplyr::mutate(index = 1:n())

    data <- draws %>%
      dplyr::left_join(data)

    p <- p +
      ggplot2::geom_line(data = data, aes(x = .data$x, y = .data$Y, group = .data$index), alpha = 0.25)
  }

  p
}

#' Plot smoothing model parameters
#'
#' @param fit fpemplus object
#' @param priors list with elements "tau" and "rho", each the density function of the corresponding prior distribution
#'
#' @import ggplot2 dplyr tidyr tibble
#' @importFrom truncnorm dtruncnorm
#'
#' @export
plot_smoothing_hyperparameters <- function(fit, priors = NULL) {
  check_is_fpemplus_class(fit)
  if(!is.null(priors)) {
    stopifnot(is.function(priors$rho))
    stopifnot(is.function(priors$tau))
  }

  if(fit$stan_data$smoothing == 0) {
    stop("Cannot plot smoothing hyperparameters because the fpemplus model was fit with smoothing disabled.")
  }

  data <- fit$posteriors$ar %>%
    tidyr::pivot_longer(c("est_rho", "est_tau")) %>%
    dplyr::mutate(name = ifelse(.data$name == "est_rho", "rho", "tau"))

  p <- ggplot2::ggplot(data, aes(x = .data$value)) +
    ggplot2::geom_density(aes(color = "Posterior")) +
    ggplot2::facet_wrap(~name, scales = "free") +
    labs(x = "value", y = "density")

  if(!is.null(priors)) {
    priors <- data %>%
      dplyr::group_by(.data$name) %>%
      dplyr::summarize(min = min(.data$value), max = max(.data$value)) %>%
      dplyr::mutate(prior = pmap(list(.data$name, .data$min, .data$max), function(name, min, max) {
        tibble::tibble(
          x = seq(min, max, length.out = 100),
          y = ifelse(name == "rho", priors$rho(.data$x), priors$tau(.data$x))
        )
      })) %>%
      tidyr::unnest(.data$prior)

    p <- p +
      ggplot2::geom_line(data = priors, aes(x = .data$x, y = .data$y, color = "Prior"))
  }
  p
}

#' Plot data model parameters
#'
#' Plots the posterior distribution of the non-sampling error parameters for each data source.
#'
#' @param fit fpemplus object
#' @param include_prior whether to include prior distribution in plot
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @export
plot_data_hyperparameters <- function(fit, include_prior = TRUE) {
  check_is_fpemplus_class(fit)
  check_is_boolean(include_prior)

  p <- fit$posteriors$nonse %>%
    ggplot2::ggplot(aes(x = .data$nonse)) +
    ggplot2::geom_density(aes(color = "Posterior")) +
    ggplot2::facet_wrap(vars(!!sym(fit$source))) +
    labs(x = "non-sampling error", y = "density")

  if(include_prior == TRUE) {
    sources <- unique(fit$data[[fit$source]])

    prior <- tibble::tibble(
      x = seq(min(fit$posteriors$nonse$nonse), max(fit$posteriors$nonse$nonse), length.out = 100),
      y = truncnorm::dtruncnorm(.data$x, a = 0, b = Inf, mean = 0, sd = 0.1)
    )

    p <- p + geom_line(data = prior, aes(x = .data$x, y = .data$y, color = "Prior"))
  }

  p
}

#' @importFrom forcats fct_reorder
plot_parameter <- function(fit, hierarchical_level, areas = c(), x, transform, sort = FALSE) {
  if(is.na(hierarchical_level)) {
    hierarchical_level <- names(fit$posteriors[[x]])[1]
  }

  check_is_fpemplus_class(fit)
  check_hierarchical_level(fit, hierarchical_level, "Omega")
  check_is_boolean(sort, "sort")
  check_areas(fit, hierarchical_level, areas)

  data <- fit$posteriors[[x]][[hierarchical_level]] %>%
    dplyr::mutate(value = transform(.data$value))

  if(length(areas) > 0) {
    data <- data %>%
      dplyr::filter(.data$name %in% areas)
  }

  data <- data %>%
    dplyr::group_by(.data$name) %>%
    tidybayes::median_qi(.data$value, .width = c(0.5, 0.8, 0.95))

  if(sort == TRUE) {
    data <- data %>%
      dplyr::mutate(name = forcats::fct_reorder(.data$name, .data$value))
  }

  ggplot(data, aes(x = .data$value, y = .data$name)) +
    ggdist::geom_pointinterval(aes(xmin = .data$.lower, xmax = .data$.upper)) +
    labs(y = hierarchical_level)
}

#' Plot the posterior of the asymptote parameter.
#'
#' Plots the posterior density of the asymptote parameter (\eqn{P}).
#'
#' @param fit fpemplus object
#' @param hierarchical_level which hierarchical level to plot
#' @param areas vector areas to include in plot
#' @param sort boolean, whether to sort the areas by the posterior median asymptote
#'
#' @export
plot_asymptote <- function(fit, hierarchical_level = NA, areas = c(), sort = FALSE) {
  plot_parameter(fit, hierarchical_level, areas, "P_tilde", function(x) 0.5 + 0.45 * plogis(x), sort) +
    labs(x = "asymptote")
}

#' Plot the level in the reference year parameter
#'
#' Plots the posterior density of the level in reference year parameter (\eqn{\Omega}).
#'
#' @param fit fpemplus object
#' @param hierarchical_level which hierarchical level to plot
#' @param areas areas to include in plot
#' @param sort boolean, whether to sort the areas by the posterior median level in reference year
#'
#' @export
plot_level <- function(fit, hierarchical_level = NA, areas = c(), sort = FALSE) {
  plot_parameter(fit, hierarchical_level, areas, "Omega", plogis, sort) +
    labs(x = "level in reference year")
}

#' Plot Probability Integral Transform
#'
#' @param fit fpemplus object
#' @param ... additional arguments to geom_histogram
#'
#' @importFrom rlang .data
#'
#' @export
plot_pit <- function(fit, ...) {
  check_is_fpemplus_class(fit)

  fit$posteriors$generated_quantities %>%
    ggplot(aes(x = .data$pit)) +
    geom_histogram(aes(y = ..density..), color = "white", ...) +
    labs(x = "PIT", y = "density")
}

#' Plot distribution of residuals
#'
#' @param fit fpemplus object
#' @param ... additional arguments to geom_histogram
#'
#' @importFrom rlang .data
#'
#' @export
plot_residuals <- function(fit, ...) {
  check_is_fpemplus_class(fit)

  fit$posteriors$generated_quantities %>%
    ggplot(aes(x = .data$resid)) +
    geom_histogram(aes(y = ..density..), color = "white", ...) +
    labs(x = "residual", y = "density")
}

#' Plot the standard deviation for each level of a hierarchical distribution
#'
#' @param fit fpemplus object
#' @param variable hierarchical distribution to plot, one of "level", "asymptote", or "splines"
#'
#' @export
plot_hierarchical_sd <- function(fit, variable) {
  check_is_fpemplus_class(fit)

  if(!(variable %in% c("level", "asymptote", "splines"))) {
    stop("variable must be one of \"level\", \"asymptote\", or \"splines\"")
  }

  v <- case_when(
    variable == "level"     ~ "Omega",
    variable == "asymptote" ~ "P_tilde",
    variable == "splines"   ~ "a",
  )

  if(v %in% c("Omega", "P_tilde")) {
    data <- fit$posteriors[[paste0(v, "_sigma")]] %>%
      dplyr::filter(.data$i > 1) %>%
      dplyr::mutate(i = fit[[paste0("hierarchical_", variable)]][.data$i])

    data %>%
      ggplot2::ggplot(aes_string(x = paste0(v, "_sigma"))) +
      ggplot2::geom_density() +
      ggplot2::facet_wrap(~i) +
      labs(x = paste0(variable, " hierarchical standard deviation"))
  }
  else if(v == "a") {
    fit$posteriors[["a_sigma"]] %>%
      filter(.data$i > 1) %>%
      dplyr::mutate(i = fit$hierarchical_splines[.data$i]) %>%
      ggplot2::ggplot(aes_string(x = "a_sigma")) +
      ggplot2::geom_density() +
      ggplot2::facet_grid(j ~ i) +
      labs(x = c("spline coefficient hierarchical standard deviation"))
  }
}
