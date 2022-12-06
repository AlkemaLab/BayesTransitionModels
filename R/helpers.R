pluralize <- function(x, singular, plural) {
  ifelse(x > 1, plural, singular)
}

#' @import tibble
hierarchical_model_matrix <- function(columns, data) {
  mat <- matrix(NA, nrow = nrow(data), ncol = 0)
  assign <- c()

  index <- tibble::tibble(i = numeric(0), column = character(0), level = character(0))

  i <- 1
  if("intercept" %in% columns) {
    mat <- cbind(mat, rep(1, nrow(mat)))
    assign <- c(assign, 0)

    index[i, ] <- tibble(i = i, column = "intercept", level = "intercept")

    i <- i + 1
  }

  for(column in columns) {
    for(l in levels(factor(data[[column]]))) {
      mat <- cbind(mat, as.numeric(data[[column]] == l))
      assign <- c(assign, which(column == columns))
      index[i, ] <- tibble::tibble(i = i, column = column, level = l)

      i <- i + 1
    }
  }

  list(
    assign = assign,
    matrix = mat,
    index = index
  )
}

#' @import purrr
hierarchical_data <- function(data, hierarchy) {
  model_matrix <- hierarchical_model_matrix(hierarchy, data)
  n_terms <- ncol(model_matrix$mat)
  re <- unique(model_matrix$assign)
  n_re <- length(re)
  re_start <- map_int(re, function(x) min(which(model_matrix$assign == x)))
  re_end   <- map_int(re, function(x) max(which(model_matrix$assign == x)))

  list(
    model_matrix = model_matrix,
    n_terms = n_terms,
    n_re = n_re,
    re_start = re_start,
    re_end = re_end
  )
}

