test_that("simulated series have the requested shape", {
  sim <- simulate_lagged_series(n_times = 20, n_vars_x = 4, n_vars_y = 3,
                                true_lag = 2, seed = 1)
  expect_equal(dim(sim$X), c(20, 4))
  expect_equal(dim(sim$Y), c(20, 3))
  expect_equal(sim$true_lag, 2)
  expect_null(simulate_lagged_series(seed = 1)$true_lag)
})

test_that("seed is reproducible and leaves the global stream untouched", {
  set.seed(99)
  before <- .Random.seed
  a <- simulate_lagged_series(seed = 5)
  expect_identical(.Random.seed, before)
  b <- simulate_lagged_series(seed = 5)
  expect_identical(a, b)
})

test_that("noise-free lagged series reproduce X shifted by true_lag", {
  sim <- simulate_lagged_series(n_times = 15, true_lag = 3, noise_sd = 0,
                                seed = 1)
  expect_equal(sim$Y[4:15, ], sim$X[1:12, ])
})

test_that("invalid arguments are rejected", {
  expect_error(simulate_lagged_series(true_lag = -1), "non-negative")
  expect_error(simulate_lagged_series(autocorr = 1), "autocorr")
})
