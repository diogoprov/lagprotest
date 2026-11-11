test_that("effect measures and RV nulls are returned", {
  sim <- simulate_lagged_series(n_times = 24, true_lag = 2, noise_sd = 0.3,
                                seed = 1)
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 199, seed = 1)
  expect_named(fit$effect, c("lag", "Z", "A", "RV"))
  expect_true(all(c("z", "a") %in% names(fit$profile)))
  expect_equal(lengths(fit$null.rv), lengths(fit$null))
  expect_equal(fit$effect$lag, fit$profile$lag[which.max(fit$profile$z)])
  expect_gt(fit$effect$Z, 2)
  expect_gt(fit$effect$A, 0)
  expect_lte(fit$effect$A, 1)
  expect_lte(fit$effect$RV, 1)
  expect_named(fit$ar1, c("series", "raw", "residual"))
  s <- summary(fit)
  expect_identical(attr(s, "effect"), fit$effect)
  expect_output(print(fit), "Effect at lag")
})

test_that("rv_coef matches the definition", {
  set.seed(2)
  X <- matrix(rnorm(40), 10)
  expect_equal(lagprotest:::rv_coef(X, X), 1)
  expect_equal(lagprotest:::rv_coef(X, X %*% diag(c(2, 1, 1, 3))),
               lagprotest:::rv_coef(X %*% diag(c(2, 1, 1, 3)), X))
})

test_that("mean_ar1 removes a pure harmonic", {
  t <- 1:24
  M <- cbind(sin(2 * pi * t / 12), cos(2 * pi * t / 24))
  expect_gt(lagprotest:::mean_ar1(M), 0.8)
  expect_true(is.na(lagprotest:::mean_ar1(M, c(24, 12))))
})

test_that("a message flags an unattainable p < 0.05", {
  sim <- simulate_lagged_series(n_times = 12, true_lag = 2, seed = 1)
  expect_message(lagged_protest(sim$X, sim$Y, lags = 0:2,
                                permutation = "cyclic", mirror = FALSE),
                 "smallest attainable")
  expect_no_message(lagged_protest(sim$X, sim$Y, lags = 0:2,
                                   permutation = "cyclic", mirror = TRUE))
})
