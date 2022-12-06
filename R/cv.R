#' Calculate coverage of held-out observations.
#'
#' @param fit fpemplus object
#' @param interval width of the predictive interval
#' @param only_held_out whether to only calculate coverage for held-out observations
#'
#' @importFrom tibble tibble
#' @importFrom rlang .data
#'
#' @export
fpem_coverage <- function(fit, interval = 0.95, only_held_out = TRUE) {
  check_is_fpemplus_class(fit)
  stopifnot(is.numeric(interval))
  stopifnot(is.logical(only_held_out))

  data_tibble <- tibble::tibble(i = 1:length(fit$held_out), held_out = fit$held_out, y = fit$stan_data$y, c = fit$stan_data$country, t = fit$stan_data$time)

  data <- fit$posteriors$generated_quantities %>%
    left_join(data_tibble, by = "i")

  if(only_held_out == TRUE) {
    data <- data %>%
      filter(.data$held_out == 1)
  }

  data %>%
    group_by(.data$i, .data$y) %>%
    tidybayes::median_qi(.data$y_pred, .width = interval) %>%
    mutate(covered = .data$y >= .data$.lower & .data$y <= .data$.upper) %>%
    left_join(data_tibble, by = "i") %>%
    left_join(fit$country_index, by = "c") %>%
    left_join(fit$time_index, by = "t")
}
