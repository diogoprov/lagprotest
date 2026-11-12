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

test_that("the RV coefficient of pair_stats() matches the definition", {
  set.seed(2)
  X <- matrix(rnorm(40), 10)
  D <- X %*% diag(c(2, 1, 1, 3))
  expect_equal(pair_stats(X, X)[2], 1)
  expect_equal(pair_stats(X, D)[2], pair_stats(D, X)[2])
  Xc <- scale(X, scale = FALSE)
  Dc <- scale(D, scale = FALSE)
  rv <- sum(crossprod(Xc, Dc)^2) /
    sqrt(sum(crossprod(Xc)^2) * sum(crossprod(Dc)^2))
  expect_equal(pair_stats(X, D)[2], rv)
})

test_that("the m2 of pair_stats() matches vegan::procrustes()", {
  set.seed(5)
  for (i in 1:50) {
    n <- sample(4:30, 1)
    A <- matrix(rnorm(n * sample(1:8, 1)), n)
    B <- matrix(rnorm(n * sample(1:8, 1)), n)
    ref <- suppressWarnings(vegan::procrustes(A, B, symmetric = TRUE)$ss)
    expect_equal(pair_stats(A, B)[1], ref, tolerance = 1e-12)
  }
  expect_equal(pair_stats(A, A)[1], 0, tolerance = 1e-12)
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

test_that("the minimum-P message gives advice that applies", {
  s8 <- simulate_lagged_series(n_times = 8, seed = 1)
  msg <- testthat::capture_messages(
    lagged_protest(s8$X, s8$Y, lags = 0:2, permutation = "cyclic")
  )
  expect_match(msg, "consider the phase scheme")
  expect_no_match(msg, "mirror")
  msg <- testthat::capture_messages(
    lagged_protest(s8$X, s8$Y, lags = 0:2, permutation = "cyclic",
                   mirror = FALSE)
  )
  expect_match(msg, "consider `mirror = TRUE` or the phase scheme")
  sim <- simulate_lagged_series(n_times = 12, seed = 1)
  expect_message(lagged_protest(sim$X, sim$Y, lags = 0:2, nperm = 19,
                                seed = 1),
                 "increase `nperm`")
})

test_that("effect_measures() returns NA when the null has no spread", {
  e <- effect_measures(0:1, c(0.4, 0.6), c(0.3, 0.2),
                       list(rep(0.5, 9), rep(0.5, 9)),
                       list(runif(9), runif(9)))
  expect_true(all(is.na(unlist(e$summary))))
  expect_true(all(is.na(e$per_lag$z)))
})

test_that("the corrected RV is NA when the null RV mean is 1", {
  e <- effect_measures(0:1, c(0.4, 0.6), c(0.9, 0.8),
                       list(c(0.5, 0.7, 0.6), c(0.6, 0.8, 0.5)),
                       list(rep(1, 3), rep(1, 3)))
  expect_true(is.na(e$summary$RV))
  expect_false(is.na(e$summary$A))
})
