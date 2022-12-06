data <- national_data()

test_that("tau_prior checks work", {
  expect_error(fpemplus(data, tau_prior = 1))
})

test_that("rho_prior checks work", {
  expect_error(fpemplus(data, rho_prior = -1))
})


test_that("spline_degree check works", {
  expect_error(fpemplus(data, spline_degree = 1))
})

test_that("num_knots check works", {
  expect_error(fpemplus(data, num_knots = -1))
})

test_that("held_out size check works", {
  expect_error(fpemplus(data, held_out = c(0, 0, 1)))
})

test_that("observations included in estimation period check works", {
  expect_error(fpemplus(
    data,
    y = "contraceptive_use_modern",
    se = "se_modern",
    year = "year",
    source = "data_series_type",
    area = "name_country",

    start_year = 2000,
    end_year = 2030
  ), "Observations included in dataset that fall outside the estimation period")
})

test_that("hierarchical setups have at least one level", {
  f <- function(
    hierarchical_level = c("intercept", "name_region", "name_country"),
    hierarchical_splines = c("intercept", "name_region", "name_country"),
    hierarchical_asymptote = c("intercept", "name_region", "name_country")
  ) {
    fpemplus(
    data,
    y = "contraceptive_use_modern",
    se = "se_modern",
    year = "year",
    source = "data_series_type",
    area = "name_country",

    hierarchical_level     = hierarchical_level,
    hierarchical_splines   = hierarchical_splines,
    hierarchical_asymptote = hierarchical_asymptote,

    start_year = 1970,
    end_year = 2030
    )
  }

  expect_error(f(hierarchical_asymptote = c()), "No hierarchical structure supplied for the asymptote")
  expect_error(f(hierarchical_level = c()), "No hierarchical structure supplied for the level in reference year")
  expect_error(f(hierarchical_splines = c()), "No hierarchical structure supplied for the spline coefficients")
})

test_that("NA checks work", {
  f <- function(data) fpemplus(
    data,
    y = "contraceptive_use_modern",
    se = "se_modern",
    year = "year",
    source = "data_series_type",
    area = "name_country",

    hierarchical_level     = c("intercept", "name_sub_region", "name_region", "name_country"),
    hierarchical_splines   = c("intercept", "name_sub_region", "name_region", "name_country"),
    hierarchical_asymptote = c("intercept", "name_sub_region", "name_region", "name_country"),

    start_year = 1970,
    end_year = 2030
  )

  expect_error(f(data %>% mutate(contraceptive_use_modern = NA)), "NAs in column contraceptive_use_modern")
  expect_error(f(data %>% mutate(se_modern = NA)), "NAs in column se_modern")
  expect_error(f(data %>% mutate(year = NA)), "NAs in column year")
  expect_error(f(data %>% mutate(data_series_type = NA)), "NAs in column data_series_type")
  expect_error(f(data %>% mutate(name_country = NA)), "NAs in column name_country")
  expect_error(f(data %>% mutate(name_region = NA)), "NAs in column name_region")
  expect_error(f(data %>% mutate(name_sub_region = NA)), "NAs in column name_sub_region")
})

test_that("zero or negative SEs are rejected", {
  expect_error(fpemplus(data %>% mutate(se_modern = rep(-1, nrow(data)))))
  expect_error(fpemplus(data %>% mutate(se_modern = rep(0, nrow(data)))))
})

test_that("changing columns names doesn't break model fitting", {
  skip_on_ci()

  data <- national_data(
    countries = c("Somalia", "India", "Palau", "Bangladesh", "Zimbabwe",
                  "Indonesia", "Kenya", "Turkey", "Mexico", "Guatemala"),
    start_year = 1970
  ) %>%
    filter(modern_method_bias == "None", has_geographical_region_bias == "N") %>%
    rename(mcpr = contraceptive_use_modern,
           standarderror = se_modern,
           time = year,
           survey = data_series_type,
           region = name_region,
           place = name_country)

  fit <- fpemplus(
    data,
    y = "mcpr",
    se = "standarderror",
    year = "time",
    source = "survey",
    area = "place",

    start_year = 1970,
    end_year = 2020,

    model = "spline",

    # Spline setup
    spline_degree = 2,
    num_knots = 7,

    # Hierarchical setup
    hierarchical_level     = c("intercept", "region", "place"),
    hierarchical_splines   = c("intercept", "region", "place"),
    hierarchical_asymptote = c("intercept", "region", "place"),

    # Stan sampler settings
    chains = 2,
    adapt_delta = 0.95,
    max_treedepth = 8,
    iter_warmup = 25,
    iter_sampling = 25,
    seed = 5,
    parallel_chains = 2,
    refresh = 50
  )

  expect_equal(fit$samples$return_codes(), c(0, 0))
})

test_that("rearranging order of observations doesn't break model fitting", {
  skip_on_ci()

  data <- national_data(
    countries = c("Somalia", "India", "Palau", "Bangladesh", "Zimbabwe",
                  "Indonesia", "Kenya", "Turkey", "Mexico", "Guatemala"),
    start_year = 1970
  ) %>%
    filter(modern_method_bias == "None", has_geographical_region_bias == "N") %>%
    arrange(year)

  fit <- fpemplus(
    data,
    y = "contraceptive_use_modern",
    se = "se_modern",
    year = "year",
    source = "data_series_type",
    area = "name_country",

    start_year = 1970,
    end_year = 2020,

    model = "spline",

    # Spline setup
    spline_degree = 2,
    num_knots = 7,

    # Hierarchical setup
    hierarchical_level     = c("intercept", "name_region", "name_country"),
    hierarchical_splines   = c("intercept", "name_region", "name_country"),
    hierarchical_asymptote = c("intercept", "name_region", "name_country"),

    # Stan sampler settings
    chains = 2,
    adapt_delta = 0.95,
    max_treedepth = 8,
    iter_warmup = 25,
    iter_sampling = 25,
    seed = 5,
    parallel_chains = 2,
    refresh = 50
  )

  expect_equal(fit$samples$return_codes(), c(0, 0))
})

test_that("end to end test of small national dataset works", {
  skip_on_ci()

  data <- national_data(
    countries = c("Somalia", "India", "Palau", "Bangladesh", "Zimbabwe",
                  "Indonesia", "Kenya", "Turkey", "Mexico", "Guatemala"),
    start_year = 1970
  ) %>%
    filter(modern_method_bias == "None", has_geographical_region_bias == "N")

  fit <- fpemplus(
    data,
    y = "contraceptive_use_modern",
    se = "se_modern",
    year = "year",
    source = "data_series_type",
    area = "name_country",

    start_year = 1970,
    end_year = 2030,

    model = "spline",

    # Spline setup
    spline_degree = 2,
    num_knots = 7,

    # Hierarchical setup
    hierarchical_level     = c("intercept", "name_region", "name_country"),
    hierarchical_splines   = c("intercept", "name_region", "name_country"),
    hierarchical_asymptote = c("intercept", "name_region", "name_country"),

    # Stan sampler settings
    chains = 2,
    adapt_delta = 0.95,
    max_treedepth = 8,
    iter_warmup = 200,
    iter_sampling = 200,
    seed = 5,
    parallel_chains = 2,
    refresh = 50
  )

  expect_equal(fit$samples$return_codes(), c(0, 0))

  test_that("stan_data spot check", {
    expect_equal(fit$stan_data$C, 9)
    expect_equal(fit$stan_data$T, 61)
  })
})
