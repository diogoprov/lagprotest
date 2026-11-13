sim_scan <- simulate_lagged_series(n_times = 12, true_lag = 2,
                                   noise_sd = 0.5, seed = 3)
fit_scan <- function(stat, nperm = 199, ...) {
  lagged_protest(sim_scan$X, sim_scan$Y, lags = 0:3, bidirectional = TRUE,
                 nperm = nperm, seed = 1, stat = stat, ...)
}

# Reference: the global statistics computed by hand from the stored
# surrogates, with explicit loops. The observed and each surrogate are scored
# in the same way within the B + 1 values of a lag: z against the other B
# values, minP as the rank among all B + 1.
reference_p <- function(fit) {
  null <- do.call(cbind, fit$null)
  obs <- fit$profile$m2
  B <- nrow(null)
  K <- ncol(null)
  z <- p <- matrix(NA_real_, B + 1, K)
  for (k in seq_len(K)) {
    v <- c(obs[k], null[, k])
    for (i in seq_len(B + 1)) {
      z[i, k] <- (v[i] - mean(v[-i])) / sd(v[-i])
      p[i, k] <- sum(v <= v[i]) / (B + 1)
    }
  }
  glob <- function(S) mean(apply(S, 1, min) <= min(S[1, ]))
  c(raw = glob(rbind(obs, null)), z = glob(z), minP = glob(p))
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

test_that("observed and permuted scores are computed alike", {
  # Swapping the observed m2 with a permuted one must swap their scores and
  # leave all others unchanged, for every statistic (exchangeability).
  set.seed(1)
  M <- matrix(runif(4 * 19), 19)
  obs <- runif(4)
  j <- 7
  M2 <- M
  M2[j, ] <- obs
  score_all <- function(stat, o, N) {
    p <- (colSums(sweep(N, 2, o, "<=")) + 1) / (nrow(N) + 1)
    sc <- scan_scores(stat, o, p, split(N, col(N)))
    rbind(sc$observed, sc$null)
  }
  for (stat in c("raw", "z", "minP")) {
    a <- score_all(stat, obs, M)
    b <- score_all(stat, M[j, ], M2)
    expect_equal(b[c(j + 1, 2:j, 1, (j + 2):20), ], a,
                 ignore_attr = TRUE, info = stat)
  }
  e1 <- effect_measures(1:4, obs, obs, split(M, col(M)), split(M, col(M)))
  st <- loo_standardize(obs, M)
  expect_equal(-st$observed, e1$per_lag$z)
})
