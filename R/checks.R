check_is_boolean <- function(x, name) {
  if(!(x %in% c(TRUE, FALSE))) {
    stop(paste0(name, " must be either TRUE or FALSE"))
  }
}

check_is_fpemplus_class <- function(x) {
  if(inherits(x, "fpemplus") == FALSE) {
    stop("fit argument must be an fpemplus object.")
  }
}

check_hierarchical_level <- function(fit, hierarchical_level, x) {
  if(!(hierarchical_level %in% names(fit$posteriors[[x]]))) {
    stop(paste0("hierarchical_level argument must be one of \"", paste0(names(fit$posteriors[[x]]), collapse = "\", \""), "\""))
  }
}

check_areas <- function(fit, hierarchical_level, areas) {
  if(hierarchical_level == "intercept") {
    if(length(areas) > 0) {
      stop("No areas can be specified when hierarchical_level is set to \"intercept\". Did you mean to choose a different hierarchical level?" )
    }
  }
  else {
    available_areas <- unique(fit$data[[hierarchical_level]])

    missing_areas <- areas[!(areas %in% available_areas)]

    if(length(missing_areas) > 0) {
      stop(paste0("The following ", pluralize(missing_areas, "area", "areas"), " could not be found in the dataset: \"", paste0(missing_areas, collapse = "\", \""), "\""))
    }
  }
}
