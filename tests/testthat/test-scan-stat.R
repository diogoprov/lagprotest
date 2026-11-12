sim_scan <- simulate_lagged_series(n_times = 12, true_lag = 2,
                                   noise_sd = 0.5, seed = 3)
fit_scan <- function(stat, nperm = 199, ...) {
  lagged_protest(sim_scan$X, sim_scan$Y, lags = 0:3, bidirectional = TRUE,
                 nperm = nperm, seed = 1, stat = stat, ...)
}

# Reference: the global statistics computed by hand from the stored
# surrogates, as in the liver example article.
reference_p <- function(fit) {
  null <- do.call(cbind, fit$null)
  B <- nrow(null)
  mu <- colMeans(null)
  s <- apply(null, 2, sd)
  z_obs <- (fit$profile$m2 - mu) / s
  z_null <- sweep(sweep(null, 2, mu), 2, s, "/")
  p_null <- apply(null, 2, rank, ties.method = "max") / (B + 1)
  c(raw = (sum(apply(null, 1, min) <= min(fit$profile$m2)) + 1) / (B + 1),
    z = (sum(apply(z_null, 1, min) <= min(z_obs)) + 1) / (B + 1),
    minP = (sum(apply(p_null, 1, min) <= min(fit$profile$p.value)) + 1) /
      (B + 1))
}

test_that("minP is the default scan statistic", {
  fit <- lagged_protest(sim_scan$X, sim_scan$Y, lags = 0:3, nperm = 99,
                        seed = 1)
  expect_identical(fit$stat, "minP")
})

test_that("the three statistics share permutations and match the reference", {
  fits <- lapply(c(raw = "raw", z = "z", minP = "minP"), fit_scan)
  ref <- reference_p(fits$raw)
  expect_equal(vapply(fits, `[[`, numeric(1), "p.global"), ref)
  expect_identical(fits$raw$null, fits$minP$null)
  expect_identical(fits$raw$profile$p.value, fits$minP$profile$p.value)
})

test_that("the optimal lag follows the statistic", {
  raw <- fit_scan("raw")
  z <- fit_scan("z")
  minp <- fit_scan("minP")
  expect_identical(raw$optimal_lag, raw$profile$lag[which.min(raw$profile$m2)])
  expect_identical(z$optimal_lag, z$effect$lag)
  prof <- minp$profile
  best <- order(prof$p.value, -prof$z, prof$m2)[1]
  expect_identical(minp$optimal_lag, prof$lag[best])
  expect_identical(minp$m2, prof$m2[best])
})

test_that("maxT adjustment uses the chosen statistic", {
  fit <- fit_scan("minP", p.adjust.method = "maxT")
  best <- fit$profile$lag == fit$optimal_lag
  expect_equal(fit$profile$p.adjusted[best], fit$p.global)
  expect_true(all(fit$profile$p.adjusted >= fit$profile$p.value))
})

test_that("stat is validated and printed", {
  expect_error(fit_scan("max"), "should be one of")
  expect_output(print(fit_scan("z")), "z scan")
})

test_that("z is refused when a null has no variance", {
  null <- list(rep(0.5, 9), seq(0.1, 0.9, by = 0.1))
  expect_error(lagprotest:::scan_scores("z", c(0.5, 0.2), c(1, 0.1), null),
               "no variance")
})

test_that("minP warns when too few permutations for the number of lags", {
  expect_message(fit_scan("minP", nperm = 99),
                 "at least 700 permutations")
  expect_no_message(fit_scan("minP", nperm = 999))
  expect_no_message(fit_scan("z", nperm = 99))
})
