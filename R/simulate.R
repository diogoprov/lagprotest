#' Simulate a pair of multivariate series with a known lag
#'
#' Generates `X` as a multivariate seasonal series plus AR(1) noise. Under
#' the null (`true_lag = NULL`), `Y` is an independent draw from the same
#' generator. Otherwise `Y` at time `t` is a linear map of `X` at time
#' `t - true_lag` plus AR(1) noise; `X` is generated over `n_times + true_lag`
#' points so that every row of `Y` has a source in `X`.
#'
#' The seasonal component is set by `waveform`. The default, `"circle"`,
#' gives every variable a sinusoid with a random phase common to the block
#' and offsets of \eqn{(j - 1)\pi/p} among the \eqn{p} variables. These
#' columns span a single plane with evenly spaced offsets, so the seasonal
#' configuration is a circle and a shift in time is an exact rotation of it,
#' which Procrustes analysis absorbs. Calibration estimated with this
#' generator is therefore a best case. The other waveforms break
#' this symmetry. Each variable follows the same periodic function, shifted
#' by an offset drawn around the block phase as set by `spread`:
#'
#' * `"sine"`: \eqn{\sin\theta};
#' * `"harmonics"`: \eqn{\sin\theta + 0.5\sin(2\theta + \psi_2) +
#'   0.3\sin(3\theta + \psi_3)}, with \eqn{\psi} drawn once per block;
#' * `"sawtooth"`: a linear rise over the cycle followed by an abrupt drop;
#' * `"pulse"`: \eqn{\exp(3\cos\theta)}, a short peak once per cycle.
#'
#' These waveforms are centered and scaled to the standard deviation of a
#' sinusoid of the same `amplitude` (\eqn{amplitude/\sqrt 2} over a full
#' cycle), so that seasonal strength is comparable among them. When `X` and
#' `Y` share a non-sinusoidal waveform (`"sawtooth"` or `"pulse"`), phase
#' surrogates do not preserve the alignment of the harmonics within each
#' cycle, and the phase scheme of [lagged_protest()] can be strongly liberal
#' (see the article "Simulating lagged series").
#'
#' Amplitudes can also vary among variables (`amp_sd_vars`) and among
#' successive cycles (`amp_sd_cycles`): each variable, or each cycle, has its
#' seasonal signal multiplied by \eqn{\exp(N(0, \sigma^2))}. Both apply to
#' every waveform, including `"circle"`.
#'
#' Setting `amplitude = 0` removes the seasonal component and leaves pure
#' AR(1) series; setting `autocorr = 0` gives white noise around the
#' seasonal signal. The article "Simulating lagged series" shows how to use
#' the function for a calibration study.
#'
#' @param n_times Number of time points (at least 3).
#' @param n_vars_x,n_vars_y Number of variables in `X` and `Y`.
#' @param true_lag `NULL` for independent series, or a non-negative integer
#'   lag by which `X` leads `Y`. `0` means a contemporaneous association,
#'   not the absence of one.
#' @param period Length of the seasonal cycle, in time steps.
#' @param amplitude Amplitude of the seasonal signal (that of an equivalent
#'   sinusoid for the non-sinusoidal waveforms).
#' @param waveform Shape of the seasonal signal: `"circle"` (default),
#'   `"sine"`, `"harmonics"`, `"sawtooth"` or `"pulse"`. See Details.
#' @param spread Phase offsets of the variables around the block phase, for
#'   every waveform except `"circle"`: `"inphase"` (default; normal with SD
#'   0.3 rad, so the variables peak at about the same time), `"anti"` (as
#'   `"inphase"`, with every second variable shifted by half a cycle) or
#'   `"random"` (uniform over the cycle).
#' @param amp_sd_vars,amp_sd_cycles Standard deviations of the log-normal
#'   multipliers of the seasonal amplitude among variables and among cycles
#'   of length `period`. `0` (default) keeps amplitudes equal.
#' @param autocorr AR(1) coefficient of the noise, in `(-1, 1)`.
#' @param noise_sd Marginal standard deviation of the AR(1) noise.
#' @param seed Optional integer seed. The global random number stream is
#'   left untouched.
#'
#' @return A list with matrices `X` (`n_times` x `n_vars_x`) and `Y`
#'   (`n_times` x `n_vars_y`), and `true_lag`.
#'
#' @examples
#' null_pair <- simulate_lagged_series(n_times = 12, seed = 42)
#' lagged_pair <- simulate_lagged_series(n_times = 24, true_lag = 3, seed = 42)
#' str(lagged_pair)
#'
#' # Abrupt seasonality that varies in amplitude between years
#' saw <- simulate_lagged_series(n_times = 36, waveform = "sawtooth",
#'                               amp_sd_cycles = 0.3, noise_sd = 0.3,
#'                               seed = 42)
#' matplot(saw$X, type = "l", lty = 1, xlab = "Time", ylab = "X")
#' @export
simulate_lagged_series <- function(n_times = 24, n_vars_x = 5, n_vars_y = 5,
                                   true_lag = NULL, period = 12,
                                   amplitude = 1,
                                   waveform = c("circle", "sine",
                                                "harmonics", "sawtooth",
                                                "pulse"),
                                   spread = c("inphase", "anti", "random"),
                                   amp_sd_vars = 0, amp_sd_cycles = 0,
                                   autocorr = 0.5, noise_sd = 1,
                                   seed = NULL) {
  waveform <- match.arg(waveform)
  spread <- match.arg(spread)
  if (!is.null(true_lag) && !is_whole(true_lag, 0)) {
    stop("`true_lag` must be NULL or a single non-negative integer.",
         call. = FALSE)
  }
  check_count(n_times, "n_times", min = 3)
  check_count(n_vars_x, "n_vars_x")
  check_count(n_vars_y, "n_vars_y")
  check_seed(seed)
  if (!is.numeric(period) || length(period) != 1 || !is.finite(period) ||
      period <= 0) {
    stop("`period` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(amplitude) || length(amplitude) != 1 ||
      !is.finite(amplitude)) {
    stop("`amplitude` must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(noise_sd) || length(noise_sd) != 1 ||
      !is.finite(noise_sd) || noise_sd < 0) {
    stop("`noise_sd` must be a single non-negative number.", call. = FALSE)
  }
  for (a in c("amp_sd_vars", "amp_sd_cycles")) {
    v <- get(a)
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < 0) {
      stop("`", a, "` must be a single non-negative number.", call. = FALSE)
    }
  }
  if (!is.numeric(autocorr) || length(autocorr) != 1 || is.na(autocorr) ||
      abs(autocorr) >= 1) {
    stop("`autocorr` must lie in (-1, 1).", call. = FALSE)
  }

  gen <- function(n, p) {
    sim_seasonal(n, p, period, amplitude, autocorr, noise_sd,
                 waveform = waveform, spread = spread,
                 amp_sd_vars = amp_sd_vars, amp_sd_cycles = amp_sd_cycles)
  }
  with_seed(seed, {
    if (is.null(true_lag)) {
      X <- gen(n_times, n_vars_x)
      Y <- gen(n_times, n_vars_y)
    } else {
      X_full <- gen(n_times + true_lag, n_vars_x)
      B <- if (n_vars_x == n_vars_y) {
        diag(n_vars_x)
      } else {
        matrix(stats::rnorm(n_vars_x * n_vars_y), n_vars_x, n_vars_y) /
          sqrt(n_vars_x)
      }
      X <- X_full[seq.int(true_lag + 1, n_times + true_lag), , drop = FALSE]
      Y <- X_full[seq_len(n_times), , drop = FALSE] %*% B +
        ar1_noise(n_times, n_vars_y, autocorr, noise_sd)
    }
    list(X = X, Y = Y, true_lag = true_lag)
  })
}

# Seasonal signal plus AR(1) noise. With waveform = "circle" the random
# draws are those of versions up to 0.7.0, so seeds reproduce earlier
# results; amplitude multipliers are drawn only when their SD is positive.
sim_seasonal <- function(n, p, period, amplitude, autocorr, noise_sd,
                         waveform = "circle", spread = "inphase",
                         amp_sd_vars = 0, amp_sd_cycles = 0) {
  t <- seq_len(n)
  phase <- stats::runif(1, 0, 2 * pi)
  if (waveform == "circle") {
    base <- vapply(seq_len(p), function(j) {
      sin(2 * pi * t / period + phase + (j - 1) * pi / p)
    }, numeric(n))
  } else {
    g <- wave_function(waveform)
    off <- switch(spread,
      inphase = stats::rnorm(p, 0, 0.3),
      anti = stats::rnorm(p, 0, 0.3) + rep(c(0, pi), length.out = p),
      random = stats::runif(p, 0, 2 * pi))
    # Center and scale over one full cycle to the SD of a unit sinusoid.
    ref <- g(2 * pi * seq_len(1000) / 1000)
    mu <- mean(ref)
    s <- sqrt(mean((ref - mu)^2)) * sqrt(2)
    base <- vapply(seq_len(p), function(j) {
      (g(2 * pi * t / period + phase + off[j]) - mu) / s
    }, numeric(n))
  }
  base <- matrix(base, n, p)
  if (amp_sd_vars > 0) {
    base <- sweep(base, 2, exp(stats::rnorm(p, 0, amp_sd_vars)), "*")
  }
  if (amp_sd_cycles > 0) {
    cyc <- floor((t - 1) / period)
    base <- base * exp(stats::rnorm(max(cyc) + 1, 0, amp_sd_cycles))[cyc + 1]
  }
  amplitude * base + ar1_noise(n, p, autocorr, noise_sd)
}

wave_function <- function(waveform) {
  switch(waveform,
    sine = function(th) sin(th),
    harmonics = {
      psi <- stats::runif(2, 0, 2 * pi)
      function(th) sin(th) + 0.5 * sin(2 * th + psi[1]) +
        0.3 * sin(3 * th + psi[2])
    },
    sawtooth = function(th) (th %% (2 * pi)) / pi - 1,
    pulse = function(th) exp(3 * cos(th))
  )
}

# AR(1) noise with stationary marginal SD noise_sd; started from the
# stationary distribution, so no burn-in is needed.
ar1_noise <- function(n, p, autocorr, noise_sd) {
  innov_sd <- noise_sd * sqrt(1 - autocorr^2)
  out <- matrix(0, n, p)
  out[1, ] <- stats::rnorm(p, 0, noise_sd)
  for (t in seq_len(n)[-1]) {
    out[t, ] <- autocorr * out[t - 1, ] + stats::rnorm(p, 0, innov_sd)
  }
  out
}
