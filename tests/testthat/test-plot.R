test_that("plot functions reject fpemplus objects", {
  x <- list()
  attr(x, "class") <- "test"

  error_message <- "fit argument must be an fpemplus object."

  expect_error(plot_indicator(x), error_message)
  expect_error(plot_smoother(x), error_message)
  expect_error(plot_transition(x), error_message)
  expect_error(plot_smoothing_hyperparameters(x), error_message)
  expect_error(plot_data_hyperparameters(x), error_message)
  expect_error(plot_asymptote(x), error_message)
  expect_error(plot_level(x), error_message)
  expect_error(plot_pit(x), error_message)
  expect_error(plot_residuals(x), error_message)
})
