#' @srrstats {G5.6, G5.6a, G5.6b} The known lag of simulated series is
#'   recovered in a defined proportion of 100 replicates, each with its own
#'   seed.
#' @srrstats {G5.7} The proportion of replicates recovering the lag rises
#'   as the series lengthen.
#' @noRd
NULL

# Simulation-based power checks. Slow, so run only as extended tests.

test_that("cyclic lagged PROTEST recovers a known lag", {
  skip_if_not_extended()
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
  skip_if_not_extended()
  sig <- vapply(seq_len(100), function(i) {
    s <- simulate_lagged_series(n_times = 24, period = 24, true_lag = 2,
                                autocorr = 0.5, noise_sd = 0.5, seed = i)
    fit <- lagged_protest(s$X, s$Y, lags = 0:4, permutation = "phase",
                          nperm = 99, seed = i, stat = "raw")
    fit$p.global < 0.05
  }, logical(1))
  expect_gte(mean(sig), 0.8)
})

test_that("lag recovery improves as the series lengthen", {
  skip_if_not_extended()
  found <- function(n) {
    mean(vapply(seq_len(100), function(i) {
      s <- simulate_lagged_series(n_times = n, period = 12, true_lag = 2,
                                  autocorr = 0.5, noise_sd = 1, seed = i)
      fit <- lagged_protest(s$X, s$Y, lags = 0:4, nperm = 19, seed = i,
                            stat = "raw")
      fit$optimal_lag == 2
    }, logical(1)))
  }
  short <- found(12)
  long <- found(48)
  expect_gt(long, short)
  expect_gte(long, 0.9)
})
