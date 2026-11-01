# Simulation-based calibration checks. Slow, so skipped on CRAN.

reject_by_lag <- function(n_rep, permutation, ...) {
  rej <- vapply(seq_len(n_rep), function(i) {
    s <- simulate_lagged_series(n_times = 12, seed = i, ...)
    fit <- lagged_protest(s$X, s$Y, lags = 0:3, permutation = permutation,
                          nperm = 99, seed = i)
    fit$profile$p.value < 0.05
  }, logical(4))
  rowMeans(rej)
}

test_that("cyclic shift keeps Type I error near alpha for seasonal series", {
  skip_on_cran()
  # 300 replicates: binomial 95% upper bound for alpha = 0.05 is ~0.075
  rate <- reject_by_lag(300, "cyclic", autocorr = 0, noise_sd = 0.3)
  expect_true(all(rate <= 0.09))
  rate <- reject_by_lag(300, "cyclic", autocorr = 0.5, noise_sd = 0.5)
  expect_true(all(rate <= 0.09))
})

test_that("free permutation inflates Type I error for seasonal series", {
  skip_on_cran()
  rate <- reject_by_lag(50, "free", autocorr = 0, noise_sd = 0.3)
  expect_true(all(rate > 0.5))
})
