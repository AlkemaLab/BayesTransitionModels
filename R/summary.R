#' Summary of an fpemplus object
#'
#' @param object fpemplus object
#' @param ... additional arguments (not currently used)
#' @import tidyr dplyr tibble
#' @importFrom rlang .data
#'
#' @export
summary.fpemplus <- function(object, ...) {
  res <- list()

  if(!is.null(object$posteriors$ar)) {
    res$ar <- object$posteriors$ar %>%
      tidyr::pivot_longer(c(.data$est_rho, .data$est_tau)) %>%
      dplyr::mutate(name = ifelse(.data$name == "rho", "rho", "tau")) %>%
      dplyr::group_by(.data$name) %>%
      dplyr::summarize(qs = list(quantile(.data$value, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))) %>%
      tidyr::unnest_wider(.data$qs) %>%
      tibble::column_to_rownames("name")
  }

  res$nonse <- object$posteriors$nonse %>%
    dplyr::group_by(.data$data_series_type) %>%
    dplyr::summarize(qs = list(quantile(.data$nonse, c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)))) %>%
    tidyr::unnest_wider(.data$qs) %>%
    tibble::column_to_rownames("data_series_type")

  res$Cn <- object$stan_data$C
  res$Tn <- object$stan_data$T

  res$hierarchical_splines   <- object$hierarchical_splines
  res$hierarchical_level     <- object$hierarchical_level
  res$hierarchical_asymptote <- object$hierarchical_asymptote

  attr(res, "class") <- "summary.fpemplus"
  res
}

#' Print fpemplus summary
#'
#' @param x fpemplus summary object
#' @param ... additional arguments (not currently used)
#'
#' @export
print.summary.fpemplus <- function(x, ...) {
  cat("fpemplus summary object\n")

  cat(glue::glue("{x$Cn} {pluralize(x$Cn, 'area', 'areas')} and {x$Tn} {pluralize(x$Tn, 'time point', 'time points')}"))
  cat("\n\n")


  cat("Hierarchical structures\n")
  cat("=======================\n")

  cat("Spline coefficients: ")
  cat(paste(x$hierarchical_splines, collapse = "/"))
  cat("\n")

  cat("Asymptote: ")
  cat(paste(x$hierarchical_level, collapse = "/"))
  cat("\n")

  cat("Level: ")
  cat(paste(x$hierarchical_asymptote, collapse = "/"))
  cat("\n")

  cat("\nPosterior quantiles\n")
  cat("==================================\n")

  if(!is.null(x$ar)) {
    cat("Smoothing:\n")
    ar_table <- x$ar %>%
      mutate_all(signif, 3)
    print(ar_table)
  }

  cat("\nNon-sampling errors:\n")
  nonse_table <- x$nonse %>%
    mutate_all(signif, 3)
  print(nonse_table)

  invisible("")
}
