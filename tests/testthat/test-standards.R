# Tests written against the rOpenSci statistical software standards
# (general standards G2, G5); the standard is named in each test.

sim <- simulate_lagged_series(n_times = 36, true_lag = 2, noise_sd = 0.3,
                              seed = 1)

test_that("missing and undefined values are errors (G2.13, G2.16)", {
  X <- sim$X
  X[2, 1] <- NaN
  expect_error(lagged_protest(X, sim$Y), "missing values \\(NA or NaN\\)")
  X[2, 1] <- Inf
  expect_error(lagged_protest(X, sim$Y), "infinite values")
  X[2, 1] <- -Inf
  expect_error(lagged_protest(sim$Y, X), "`Y` contains infinite values")
  d <- stats::dist(sim$X)
  d[1] <- Inf
  expect_error(lagged_protest(d, sim$Y), "`X` contains infinite values")
})

test_that("zero-length and too-short input are errors (G5.8a)", {
  expect_error(lagged_protest(sim$X[0, ], sim$Y[0, ]),
               "at least 3 rows \\(time points\\); they have 0")
  expect_error(lagged_protest(sim$X[1:2, ], sim$Y[1:2, ]),
               "they have 2")
})

test_that("classed numeric columns are used as numbers (G2.9, G2.11)", {
  d <- as.data.frame(sim$X)
  d$V1 <- as.difftime(d$V1, units = "days")
  expect_message(fit <- lagged_protest(d, sim$Y, lags = 0:3, nperm = 99,
                                       seed = 1),
                 "Dropping the class of column\\(s\\) V1")
  ref <- lagged_protest(sim$X, sim$Y, lags = 0:3, nperm = 99, seed = 1)
  expect_identical(fit$profile$m2, ref$profile$m2)
})

test_that("list, factor and character columns are errors (G2.12, G5.8b)", {
  d <- as.data.frame(sim$X)
  d$L <- as.list(seq_len(nrow(d)))
  expect_error(lagged_protest(d, sim$Y), "list columns: L")
  d <- as.data.frame(sim$X)
  d$F <- factor(rep(c("a", "b"), length.out = nrow(d)))
  expect_error(lagged_protest(d, sim$Y), "columns: F")
  d <- as.data.frame(sim$X)
  d$C <- "a"
  expect_error(lagged_protest(d, sim$Y), "columns: C")
  M <- matrix(as.character(sim$X), nrow(sim$X))
  expect_error(lagged_protest(M, sim$Y), "`X` must be numeric \\(matrix")
})

test_that("returned values contain no missing or undefined values (G5.3)", {
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, bidirectional = TRUE,
                        nperm = 199, seed = 1)
  num <- vapply(fit$profile, is.numeric, logical(1))
  expect_true(all(is.finite(as.matrix(fit$profile[num]))))
  expect_true(all(is.finite(fit$residuals$residual)))
  expect_true(all(is.finite(c(fit$m2, fit$r, fit$p.value, fit$p.adjusted,
                              fit$p.global))))
})

test_that("m2 matches vegan::protest() on a standard data set (G5.0, G5.4b)", {
  skip_if_not_installed("vegan")
  env <- new.env()
  utils::data("varespec", "varechem", package = "vegan", envir = env)
  fit <- lagged_protest(env$varespec, env$varechem, lags = 0,
                        permutation = "free", nperm = 999, scale = FALSE,
                        seed = 1)
  ref <- suppressWarnings(vegan::protest(env$varespec, env$varechem,
                                         permutations = 999))
  expect_equal(fit$m2, ref$ss, tolerance = 1e-10)
  expect_equal(fit$r, ref$t0, tolerance = 1e-10)
  # Both detect concordance between vegetation and soil chemistry; vegan
  # draws its own permutations, so only the decision is compared.
  expect_lt(fit$p.value, 0.05)
  expect_lt(ref$signif, 0.05)
})

test_that("trivial noise does not change the results (G5.9a)", {
  set.seed(2)
  eps <- .Machine$double.eps
  Xn <- sim$X + matrix(stats::runif(length(sim$X), -eps, eps),
                       nrow(sim$X))
  fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 199, seed = 1)
  fit_n <- lagged_protest(Xn, sim$Y, lags = 0:4, nperm = 199, seed = 1)
  expect_equal(fit_n$profile$m2, fit$profile$m2, tolerance = 1e-8)
  expect_identical(fit_n$optimal_lag, fit$optimal_lag)
  expect_identical(fit_n$p.global, fit$p.global)
})

test_that("different seeds give consistent results (G5.9b)", {
  fits <- lapply(1:3, function(s) {
    lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 999, seed = s)
  })
  lags <- vapply(fits, `[[`, numeric(1), "optimal_lag")
  p <- vapply(fits, `[[`, numeric(1), "p.global")
  expect_identical(unique(lags), 2)
  # The observed statistic does not depend on the seed ...
  expect_identical(fits[[1]]$profile$m2, fits[[2]]$profile$m2)
  # ... and the global p-values agree within Monte Carlo error.
  expect_lt(diff(range(p)), 0.05)
})
