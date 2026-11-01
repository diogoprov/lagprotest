sim <- simulate_lagged_series(n_times = 12, n_vars_y = 4, true_lag = 2,
                              noise_sd = 0.3, seed = 1)

test_that("m2 at lag 0 matches vegan::protest()", {
  X <- sim$X
  Y <- sim$Y
  fit <- lagged_protest(X, Y, lags = 0, permutation = "free", nperm = 99,
                        scale = FALSE, seed = 1)
  ref <- suppressWarnings(vegan::protest(X, Y, permutations = 99))
  expect_equal(fit$profile$m2, ref$ss)
  expect_equal(fit$profile$r, ref$t0)
})

test_that("lags pair the intended rows", {
  idx <- lag_index(10, 3)
  expect_equal(idx$x, 1:7)
  expect_equal(idx$y, 4:10)
  idx <- lag_index(10, -3)
  expect_equal(idx$x, 4:10)
  expect_equal(idx$y, 1:7)
})

test_that("an exact lagged copy gives m2 = 0 at the true lag", {
  set.seed(2)
  Z <- matrix(rnorm(15 * 3), 15, 3)
  X <- Z[4:15, ]
  Y <- Z[2:13, ]
  # scale = FALSE: z-scoring uses the full series, so truncated windows of
  # an exact copy would differ slightly after scaling
  fit <- lagged_protest(X, Y, lags = 0:4, scale = FALSE)
  expect_equal(fit$optimal_lag, 2L)
  expect_equal(fit$m2, 0, tolerance = 1e-10)
})

test_that("cyclic permutation enumerates 2(n - 1) or n - 1 shifts", {
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:3)
  expect_equal(unique(fit$profile$nperm), 22L)
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:3, mirror = FALSE)
  expect_equal(unique(fit$profile$nperm), 11L)
  expect_true(all(fit$profile$p.value >= 1 / 12))
})

test_that("output is well formed", {
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:4)
  expect_s3_class(fit, "lagged_protest")
  expect_equal(fit$profile$lag, 0:4)
  expect_equal(fit$profile$n, 12 - 0:4)
  expect_equal(fit$m2, min(fit$profile$m2))
  expect_true(all(fit$profile$p.adjusted >= fit$profile$p.value))
  expect_equal(nrow(residuals(fit)), 12 - fit$optimal_lag)
  expect_true(summary(fit)$optimal[summary(fit)$lag == fit$optimal_lag])
  expect_output(print(fit), "Optimal lag")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(fit))
})

test_that("p-values are raw by default and adjusted on request", {
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:4)
  expect_equal(fit$p.adjust.method, "none")
  expect_equal(fit$profile$p.adjusted, fit$profile$p.value)
  best <- fit$profile$lag == fit$optimal_lag
  expect_equal(fit$p.value, fit$profile$p.value[best])
  expect_equal(fit$p.adjusted, fit$p.value)
  expect_equal(fit$n_perm, 22L)
  expect_output(print(fit), "not corrected")

  holm <- lagged_protest(sim$X, sim$Y, lags = 0:4, p.adjust.method = "holm")
  expect_equal(holm$profile$p.adjusted,
               p.adjust(holm$profile$p.value, method = "holm"))
  expect_equal(holm$p.value, fit$p.value)
  expect_equal(holm$p.adjusted, holm$profile$p.adjusted[best])
  expect_output(print(holm), "adjusted p")
})

test_that("the global p-value is the max-statistic over lags", {
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:4)
  null_min <- do.call(pmin, unname(fit$null))
  expect_equal(fit$p.global,
               (sum(null_min <= fit$m2) + 1) / (length(null_min) + 1))
  expect_gte(fit$p.global, fit$p.value)
  expect_output(print(fit), "Global p")

  maxt <- lagged_protest(sim$X, sim$Y, lags = 0:4, p.adjust.method = "maxT")
  expect_equal(maxt$p.adjusted, maxt$p.global)
  expect_true(all(maxt$profile$p.adjusted >= maxt$profile$p.value))
})

test_that("phase surrogates keep the amplitude and cross-spectra of Y", {
  Y <- scale(sim$Y)
  Ys <- with_seed(3, phase_surrogate(Y))
  expect_equal(dim(Ys), dim(Y))
  expect_equal(Mod(mvfft(Ys)), Mod(mvfft(Y)), ignore_attr = TRUE)
  expect_equal(crossprod(Ys), crossprod(Y), ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(Ys, Y, check.attributes = FALSE)))
  # odd n
  Y7 <- Y[1:11, ]
  Ys7 <- with_seed(3, phase_surrogate(Y7))
  expect_equal(Mod(mvfft(Ys7)), Mod(mvfft(Y7)), ignore_attr = TRUE)
})

test_that("the phase scheme is reproducible and uses nperm surrogates", {
  a <- lagged_protest(sim$X, sim$Y, lags = 0:4, permutation = "phase",
                      nperm = 99, seed = 4)
  b <- lagged_protest(sim$X, sim$Y, lags = 0:4, permutation = "phase",
                      nperm = 99, seed = 4)
  expect_identical(a$profile, b$profile)
  expect_identical(a$p.global, b$p.global)
  expect_equal(unique(a$profile$nperm), 99L)
  expect_true(is.na(a$mirror))
  expect_output(print(a), "phase-randomized")
})

test_that("bidirectional adds negative lags", {
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:2, bidirectional = TRUE)
  expect_equal(fit$profile$lag, -2:2)
  expect_equal(fit$profile$n, 12 - abs(-2:2))
})

test_that("free permutation is reproducible with a seed", {
  a <- lagged_protest(sim$X, sim$Y, lags = 0:2, permutation = "free",
                      nperm = 49, seed = 10)
  b <- lagged_protest(sim$X, sim$Y, lags = 0:2, permutation = "free",
                      nperm = 49, seed = 10)
  expect_identical(a$profile, b$profile)
  expect_equal(unique(a$profile$nperm), 49L)
})

test_that("dist inputs are accepted", {
  fit <- lagged_protest(dist(sim$X), dist(sim$Y), lags = 0:3)
  expect_s3_class(fit, "lagged_protest")
  expect_equal(fit$optimal_lag, 2L)
})

test_that("invalid inputs are rejected", {
  expect_error(lagged_protest(sim$X, sim$Y[-1, ]), "same number of rows")
  expect_error(lagged_protest(sim$X, sim$Y, lags = -1), "non-negative")
  expect_error(lagged_protest(sim$X, sim$Y, lags = 10), "at least 3")
  X <- sim$X
  X[1, 1] <- NA
  expect_error(lagged_protest(X, sim$Y), "missing")
})
