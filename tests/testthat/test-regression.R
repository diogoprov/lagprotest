# Numerical regression: fixed-seed results recorded with lagprotest 0.6.0
# (commit ec24c79). Refactoring or speeding up the permutation loop must not
# change them beyond floating-point rounding.

reg_sim <- simulate_lagged_series(n_times = 12, n_vars_y = 4, true_lag = 2,
                                  noise_sd = 0.3, seed = 1)

expect_fit <- function(fit, ref) {
  tol <- 1e-9
  expect_equal(fit$profile$m2, ref$m2, tolerance = tol)
  expect_identical(fit$profile$p.value, ref$p)
  expect_identical(fit$p.global, ref$pg)
  expect_equal(fit$profile$z, ref$z, tolerance = tol)
  expect_identical(fit$effect$lag, as.integer(ref$eff$lag))
  expect_equal(fit$effect[c("Z", "A", "RV")], ref$eff[c("Z", "A", "RV")],
               tolerance = tol)
  expect_equal(fit$null[[1]][1:5], ref$null1, tolerance = tol)
}

test_that("phase scheme results are unchanged", {
  fit <- lagged_protest(reg_sim$X, reg_sim$Y, lags = 0:3, nperm = 199,
                        seed = 1)
  expect_fit(fit, list(
    m2 = c(0.141479433217, 0.123578369492, 0.0650119923712, 0.123875165313),
    p = c(0.88, 0.535, 0.005, 0.775), pg = 0.025,
    z = c(-1.25857925054, -0.0877041852856, 3.09676701116, -0.729076174741),
    eff = list(lag = 2, Z = 2.86632111958, A = 0.366846614416,
               RV = 0.386118349814),
    null1 = c(0.102412460157, 0.0985785138123, 0.0965716957053,
              0.117221681105, 0.119602884475)
  ))
})

test_that("cyclic scheme results are unchanged", {
  fit <- lagged_protest(reg_sim$X, reg_sim$Y, lags = 0:3,
                        bidirectional = TRUE, permutation = "cyclic")
  expect_fit(fit, list(
    m2 = c(0.105080220958, 0.0806343041022, 0.0977316638572, 0.141479433217,
           0.123578369492, 0.0650119923712, 0.123875165313),
    p = c(9, 3, 5, 21, 12, 1, 17) / 23, pg = 2 / 23,
    z = c(0.486172166395, 1.65132781646, 1.03349987536, -1.17268658879,
          -0.0597979888649, 3.5933186382, -0.724521824894),
    eff = list(lag = 2, Z = 2.92675368373, A = 0.309068425559,
               RV = 0.304073345968),
    null1 = c(0.0720103685461, 0.0957026133687, 0.151214526656,
              0.123885131106, 0.0796637683079)
  ))
})

test_that("free scheme results are unchanged", {
  fit <- lagged_protest(reg_sim$X, reg_sim$Y, lags = 0:2,
                        permutation = "free", nperm = 99, seed = 3)
  expect_fit(fit, list(
    m2 = c(0.141479433217, 0.123578369492, 0.0650119923712),
    p = c(0.01, 0.01, 0.01), pg = 0.01,
    z = c(8.14937790951, 6.81125061601, 6.93772260082),
    eff = list(lag = 0, Z = 7.54377467753, A = 0.821802671966,
               RV = 0.821071239299),
    null1 = c(0.879567701932, 0.899569378342, 0.940480607743,
              0.731959632645, 0.832017768256)
  ))
})
