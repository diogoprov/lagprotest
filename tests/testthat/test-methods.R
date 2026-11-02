sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
                              seed = 1)
fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 99, seed = 1)

test_that("as.data.frame() returns the profile with null quantile and flag", {
  df <- as.data.frame(fit)
  expect_s3_class(df, "data.frame")
  expect_equal(df[names(fit$profile)], fit$profile)
  expect_equal(df$optimal, df$lag == fit$optimal_lag)
  expect_equal(sum(df$optimal), 1)
  expect_equal(df$null.quantile[1],
               unname(quantile(fit$null[[1]], 0.05)))
  expect_equal(as.data.frame(fit, alpha = 0.1)$null.quantile[3],
               unname(quantile(fit$null[[3]], 0.1)))
})

test_that("autoplot() returns a ggplot with the profile", {
  skip_if_not_installed("ggplot2")
  p <- ggplot2::autoplot(fit)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), nrow(fit$profile))
  expect_match(p$labels$subtitle, "global p")
  built <- ggplot2::ggplot_build(p)
  expect_length(built$data, 5)
})
