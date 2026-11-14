test_that("simulated series have the requested shape", {
  sim <- simulate_lagged_series(n_times = 20, n_vars_x = 4, n_vars_y = 3,
                                true_lag = 2, seed = 1)
  expect_equal(dim(sim$X), c(20, 4))
  expect_equal(dim(sim$Y), c(20, 3))
  expect_equal(sim$true_lag, 2)
  expect_null(simulate_lagged_series(seed = 1)$true_lag)
})

test_that("seed is reproducible and leaves the global stream untouched", {
  set.seed(99)
  before <- .Random.seed
  a <- simulate_lagged_series(seed = 5)
  expect_identical(.Random.seed, before)
  b <- simulate_lagged_series(seed = 5)
  expect_identical(a, b)
})

test_that("noise-free lagged series reproduce X shifted by true_lag", {
  sim <- simulate_lagged_series(n_times = 15, true_lag = 3, noise_sd = 0,
                                seed = 1)
  expect_equal(sim$Y[4:15, ], sim$X[1:12, ])
})

test_that("invalid arguments are rejected", {
  expect_error(simulate_lagged_series(true_lag = -1), "non-negative")
  expect_error(simulate_lagged_series(autocorr = 1), "autocorr")
  expect_error(simulate_lagged_series(true_lag = 1.5), "non-negative")
  expect_error(simulate_lagged_series(n_times = 2), "`n_times`")
  expect_error(simulate_lagged_series(n_vars_x = 0), "`n_vars_x`")
  expect_error(simulate_lagged_series(n_vars_y = 1.5), "`n_vars_y`")
  expect_error(simulate_lagged_series(period = 0), "`period`")
  expect_error(simulate_lagged_series(amplitude = NA), "`amplitude`")
  expect_error(simulate_lagged_series(noise_sd = -1), "`noise_sd`")
  expect_error(simulate_lagged_series(seed = "a"), "`seed`")
})

test_that("the circle waveform reproduces earlier seeds", {
  # Values from version 0.7.0 (sin with evenly spaced offsets).
  expect_identical(simulate_lagged_series(seed = 3),
                   simulate_lagged_series(seed = 3, waveform = "circle"))
  s <- simulate_lagged_series(n_times = 12, n_vars_x = 2, noise_sd = 0,
                              amplitude = 2, seed = 1)
  t <- 1:12
  phi <- with_seed(1, stats::runif(1, 0, 2 * pi))
  expect_equal(s$X[, 2], 2 * sin(2 * pi * t / 12 + phi + pi / 2))
})

test_that("non-circle waveforms are scaled like a sinusoid", {
  for (w in c("sine", "harmonics", "pulse")) {
    s <- simulate_lagged_series(n_times = 1200, period = 12, waveform = w,
                                noise_sd = 0, amplitude = 2, seed = 1)
    sds <- apply(s$X, 2, function(x) sqrt(mean((x - mean(x))^2)))
    expect_equal(sds, rep(sqrt(2), 5), tolerance = 1e-3, info = w)
    expect_equal(unname(colMeans(s$X)), rep(0, 5), tolerance = 1e-3,
                 info = w)
  }
  # The sawtooth sampled 12 times per cycle is close to, not exactly, that.
  s <- simulate_lagged_series(n_times = 1200, waveform = "sawtooth",
                              noise_sd = 0, seed = 1)
  expect_equal(apply(s$X, 2, stats::sd), rep(sqrt(0.5), 5), tolerance = 0.01)
})

test_that("waveforms have the expected shape", {
  saw <- simulate_lagged_series(n_times = 48, waveform = "sawtooth",
                                noise_sd = 0, seed = 2)
  d <- diff(saw$X[, 1])
  # One abrupt drop per cycle, rises otherwise.
  expect_equal(sum(d < 0), 4)
  expect_true(all(abs(d[d > 0] - median(d[d > 0])) < 1e-8))
  pulse <- simulate_lagged_series(n_times = 48, waveform = "pulse",
                                  noise_sd = 0, seed = 2)
  expect_gt(max(pulse$X[, 1]) / abs(min(pulse$X[, 1])), 2)
  sine <- simulate_lagged_series(n_times = 48, waveform = "sine",
                                 noise_sd = 0, seed = 2)
  # A sinusoid is a combination of sin and cos of the period.
  t <- 1:48
  fit <- lm(sine$X ~ cos(2 * pi * t / 12) + sin(2 * pi * t / 12))
  expect_lt(max(abs(residuals(fit))), 1e-10)
})

test_that("spread sets the phase offsets among variables", {
  cor_1 <- function(spread) {
    s <- simulate_lagged_series(n_times = 120, n_vars_x = 4,
                                waveform = "sine", spread = spread,
                                noise_sd = 0, seed = 4)
    stats::cor(s$X)[1, ]
  }
  # Same seed, same offsets: "anti" flips the sign of every second column.
  inphase <- cor_1("inphase")
  expect_true(all(inphase > 0))
  expect_equal(cor_1("anti"), inphase * c(1, -1, 1, -1))
})

test_that("amplitude varies among variables and cycles when asked", {
  s <- simulate_lagged_series(n_times = 36, waveform = "sine",
                              amp_sd_vars = 0.5, noise_sd = 0, seed = 1)
  r <- apply(s$X, 2, function(x) diff(range(x)))
  expect_gt(max(r) / min(r), 1.2)
  s <- simulate_lagged_series(n_times = 36, n_vars_x = 1,
                              amp_sd_cycles = 0.5, noise_sd = 0, seed = 1)
  r <- tapply(s$X[, 1], rep(1:3, each = 12), function(x) diff(range(x)))
  expect_gt(max(r) / min(r), 1.2)
  # Within a cycle the circle keeps its shape: amplitude is a scalar factor.
  base <- simulate_lagged_series(n_times = 36, n_vars_x = 1, noise_sd = 0,
                                 seed = 1)
  ratio <- s$X[, 1] / base$X[, 1]
  expect_lt(max(tapply(ratio, rep(1:3, each = 12), stats::sd)), 1e-8)
})

test_that("lagged series follow X for every waveform", {
  for (w in c("sine", "harmonics", "sawtooth", "pulse")) {
    sim <- simulate_lagged_series(n_times = 15, true_lag = 3, noise_sd = 0,
                                  waveform = w, amp_sd_cycles = 0.3,
                                  seed = 1)
    expect_equal(sim$Y[4:15, ], sim$X[1:12, ], info = w)
  }
})

test_that("new arguments are validated", {
  expect_error(simulate_lagged_series(waveform = "square"), "should be one")
  expect_error(simulate_lagged_series(spread = "wide"), "should be one")
  expect_error(simulate_lagged_series(amp_sd_vars = -1), "`amp_sd_vars`")
  expect_error(simulate_lagged_series(amp_sd_cycles = NA), "`amp_sd_cycles`")
})
