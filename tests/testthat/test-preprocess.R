test_that("remove_cycle_means() subtracts the mean of each cycle position", {
  M <- cbind(a = rep(1:4, 3) + rep(c(0, 10, 20), each = 4), b = 1:12)
  out <- remove_cycle_means(M, 4)
  expect_identical(dimnames(out), dimnames(M))
  expect_equal(unname(out[, "a"]), rep(c(-10, 0, 10), each = 4))
  pos <- rep(1:4, 3)
  expect_equal(as.vector(tapply(out[, "b"], pos, mean)), rep(0, 4))
})

test_that("prais_winsten() keeps n rows and removes AR(1) structure", {
  M <- matrix(c(1, 3, 2, 5, 4, 6), ncol = 1)
  pw <- prais_winsten(M)
  x <- M[, 1] - mean(M[, 1])
  r <- sum(x[-1] * x[-6]) / sum(x^2)
  expect_equal(pw$rho, r)
  expect_equal(pw$series[, 1], c(sqrt(1 - r^2) * x[1], x[-1] - r * x[-6]))
  set.seed(1)
  e <- stats::arima.sim(list(ar = 0.8), n = 2000)
  pw <- prais_winsten(cbind(e, stats::rnorm(2000)))
  expect_equal(nrow(pw$series), 2000)
  expect_lt(max(abs(pw$rho - c(0.8, 0))), 0.06)
  expect_lt(abs(mean_ar1(pw$series[, 1, drop = FALSE])), 0.05)
  # A constant column is left at zero with rho = 0.
  expect_equal(prais_winsten(matrix(2, 5, 1))$rho, 0)
})

test_that("lagged_protest() applies deseason and prewhiten in order", {
  s <- simulate_lagged_series(n_times = 30, true_lag = 1, waveform = "pulse",
                              seed = 3)
  expect_message(
    fit <- lagged_protest(s$X, s$Y, lags = 0:2, period = 12, nperm = 99,
                   deseason = "monthly", prewhiten = "prais", seed = 1),
    "first 24 of 30")
  expect_identical(fit$deseason, "monthly")
  expect_identical(fit$prewhiten, "prais")
  expect_named(fit$rho, c("X", "Y"))
  # Same result by hand: truncate, remove means, prewhiten, scale.
  prep <- function(M) {
    M <- prais_winsten(remove_cycle_means(M[1:24, ], 12))$series
    scale(M)
  }
  hand <- lagged_protest(prep(s$X), prep(s$Y), lags = 0:2, nperm = 99,
                         seed = 1)
  expect_equal(fit$profile[c("m2", "p.value")],
               hand$profile[c("m2", "p.value")])
  expect_equal(fit$rho$X,
               prais_winsten(remove_cycle_means(s$X[1:24, ], 12))$rho)
  expect_output(print(fit), "seasonal means removed, then Prais-Winsten")
})

test_that("defaults leave results unchanged", {
  s <- simulate_lagged_series(n_times = 12, true_lag = 2, seed = 1)
  a <- lagged_protest(s$X, s$Y, lags = 0:3, nperm = 99, seed = 1)
  b <- lagged_protest(s$X, s$Y, lags = 0:3, nperm = 99, seed = 1,
                      deseason = "none", prewhiten = "none")
  a$call <- b$call <- NULL
  expect_identical(a, b)
  expect_null(a$rho)
  expect_false(any(grepl("Preprocessing", utils::capture.output(print(a)))))
})

test_that("deseason works on dist inputs", {
  s <- simulate_lagged_series(n_times = 24, seed = 2)
  fit <- lagged_protest(stats::dist(s$X), s$Y, lags = 0:2, period = 12,
                        deseason = "monthly", nperm = 99, seed = 1)
  expect_s3_class(fit, "lagged_protest")
})

test_that("deseason and prewhiten are validated, with warnings", {
  s <- simulate_lagged_series(n_times = 24, seed = 1)
  expect_error(lagged_protest(s$X, s$Y, deseason = "monthly"), "needs `period`")
  expect_error(suppressMessages(
    lagged_protest(s$X, s$Y, lags = 0:2, period = 12.5, deseason = "monthly")),
    "whole number")
  expect_error(lagged_protest(s$X[1:12, ], s$Y[1:12, ], lags = 0:2,
                              period = 12, deseason = "monthly"),
               "two whole cycles")
  expect_error(lagged_protest(s$X, s$Y, deseason = "yearly"), "should be one")
  expect_warning(
    lagged_protest(s$X, s$Y, lags = 0:2, period = 12, prewhiten = "prais",
                   nperm = 99, seed = 1),
    "without `deseason`")
  expect_message(
    lagged_protest(s$X[1:12, ], s$Y[1:12, ], lags = 0:2, period = 12,
                   nperm = 99, seed = 1),
    "single cycle")
  # A purely seasonal column has nothing left after removing the means.
  X <- cbind(s$X, rep(1:12, 2))
  expect_error(lagged_protest(X, s$Y, lags = 0:2, period = 12,
                              deseason = "monthly", nperm = 99),
               "constant columns after preprocessing")
})
