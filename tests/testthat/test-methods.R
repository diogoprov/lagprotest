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

test_that("summary() prints the profile and the scan-level measures", {
  s <- summary(fit)
  expect_s3_class(s, c("summary.lagged_protest", "data.frame"))
  expect_identical(attr(s, "p.global"), fit$p.global)
  expect_identical(attr(s, "optimal_lag"), fit$optimal_lag)
  out <- capture.output(print(s))
  expect_true(any(grepl("global p = ", out)))
  expect_true(any(grepl("Effect at lag", out)))
  expect_true(any(grepl("autocorrelation", out)))
  # Subsets drop the attributes and print as a plain profile
  expect_no_error(capture.output(print(s[s$optimal, c("lag", "m2")])))
})

test_that("print() notes when the optimal and effect lags differ", {
  s <- simulate_lagged_series(n_times = 12, seed = 3)
  f <- lagged_protest(s$X, s$Y, lags = 0:3, bidirectional = TRUE,
                      nperm = 199, seed = 1, stat = "raw")
  expect_false(f$optimal_lag == f$effect$lag)
  expect_output(print(f), "smallest m2.*differ")
  expect_false(any(grepl("differ", capture.output(print(fit)))))
})

test_that("print() uses the p-based wording for minP and z scans", {
  g <- fit
  g$stat <- "minP"
  g$effect$lag <- setdiff(g$profile$lag, g$optimal_lag)[1]
  expect_output(print(g), "smallest p\\).*differ\\.")
  expect_false(any(grepl("favors long lags", capture.output(print(g)))))
})
