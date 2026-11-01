# Simulation-based power checks. Slow, so skipped on CRAN.

test_that("cyclic lagged PROTEST recovers a known lag", {
  skip_on_cran()
  res <- vapply(seq_len(100), function(i) {
    s <- simulate_lagged_series(n_times = 12, true_lag = 2, autocorr = 0.5,
                                noise_sd = 0.5, seed = i)
    fit <- lagged_protest(s$X, s$Y, lags = 0:3, permutation = "cyclic")
    c(found = fit$optimal_lag == 2,
      significant = fit$profile$p.value[fit$profile$lag == 2] < 0.05)
  }, logical(2))
  expect_gte(mean(res["found", ]), 0.8)
  expect_gte(mean(res["significant", ]), 0.7)
})

test_that("the phase global p-value detects a known lag", {
  skip_on_cran()
  sig <- vapply(seq_len(100), function(i) {
    s <- simulate_lagged_series(n_times = 24, period = 24, true_lag = 2,
                                autocorr = 0.5, noise_sd = 0.5, seed = i)
    fit <- lagged_protest(s$X, s$Y, lags = 0:4, permutation = "phase",
                          nperm = 99, seed = i)
    fit$p.global < 0.05
  }, logical(1))
  expect_gte(mean(sig), 0.8)
})
