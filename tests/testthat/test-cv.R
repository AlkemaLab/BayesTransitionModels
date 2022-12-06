expect_invisible({
  skip_on_ci()
  skip("Temporarily skipping")
  data <- national_data(
    countries = c("Somalia", "India", "Palau", "Bangladesh", "Zimbabwe",
                  "Indonesia", "Kenya", "Turkey", "Mexico", "Guatemala"),
    start_year = 1970
  )

  held_out <- rep(FALSE, nrow(data))
  held_out[seq(1, 20, 5)] <- TRUE # Leave out a few data points
  fit <- xfun::cache_rds({
  fpemplus(
    data,
    y = "contraceptive_use_modern",
    se = "se_modern",
    year = "year",
    source = "data_series_type",
    area = "name_country",

    held_out = held_out,

    start_year = 1970,
    end_year = 2030,

    model = "spline",

    # Spline setup
    spline_degree = 2,
    num_knots = 7,

    # Hierarchical setup
    hierarchical_level     = c("intercept", "name_country"),
    hierarchical_splines   = c("intercept", "name_country"),
    hierarchical_asymptote = c("intercept", "name_country"),

    # Stan sampler settings
    adapt_delta = 0.9,
    max_treedepth = 8,
    iter_warmup = 200,
    iter_sampling = 200,
    seed = 5,
    parallel_chains = 2,
    refresh = 50,
    chains = 2
  )
  })
})

test_that("fpemplus object stores which indicators are held out", {
  skip_on_ci()
  expect_equal(fit$held_out, as.numeric(held_out))
})

