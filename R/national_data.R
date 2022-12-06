#' National level dataset.
#'
#' @param countries vector of countries
#' @param start_year first year
#'
#' @import fpemlocal
#' @importFrom stats median
#' @importFrom rlang .data
#'
#' @export
national_data <- function(countries = c(), start_year = 1970) {
  data <- fpemlocal::contraceptive_use %>%
    filter(.data$is_in_union == "Y", .data$age_range == "15-49", .data$group_type_relative_to_baseline == "MW") %>%
    filter(!is.na(.data$contraceptive_use_modern)) %>%
    select(.data$division_numeric_code, .data$start_date, .data$end_date, .data$data_series_type, .data$contraceptive_use_modern, .data$contraceptive_use_traditional, .data$se_modern, .data$contraceptive_use_any, .data$se_traditional, .data$modern_method_bias, .data$has_geographical_region_bias) %>%
    left_join(fpemlocal::divisions, by = "division_numeric_code") %>%
    mutate(year = floor((.data$start_date + .data$end_date) / 2),
           se_modern = ifelse(.data$se_modern == 0, NA, .data$se_modern),
           se_modern_original = .data$se_modern,
           s = .data$se_modern,
           included = 1) %>%
    filter(.data$year >= start_year)

  # Impute missing SEs following strategy of Cahill 2018 et al. appendix p. 16
  max_ses <- data %>%
    filter(!is.na(.data$se_modern)) %>%
    group_by(.data$division_numeric_code) %>%
    summarize(max_se_modern = max(.data$se_modern))

  median_max_se <- stats::median(max_ses$max_se_modern)

  data <- data %>%
    left_join(max_ses, by = "division_numeric_code") %>%
    mutate(se_modern = ifelse(is.na(.data$se_modern) & !is.na(.data$max_se_modern), .data$max_se_modern, .data$se_modern),
           se_modern = ifelse(is.na(.data$se_modern), median_max_se, .data$se_modern))

  if(length(countries) > 0) {
    data <- data %>% filter(.data$name_country %in% countries)
  }

  return(data)
}
