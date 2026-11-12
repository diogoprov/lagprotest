#' Simulate a pair of multivariate series with a known lag
#'
#' Generates `X` as a multivariate seasonal series (sinusoids with a random
#' common phase and staggered offsets among variables) plus AR(1) noise.
#' Under the null (`true_lag = NULL`), `Y` is an independent draw from the
#' same generator. Otherwise `Y` at time `t` is a linear map of `X` at time
#' `t - true_lag` plus AR(1) noise; `X` is generated over `n_times + true_lag`
#' points so that every row of `Y` has a source in `X`.
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
#' @param amplitude Amplitude of the seasonal sinusoids.
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
#' @export
simulate_lagged_series <- function(n_times = 24, n_vars_x = 5, n_vars_y = 5,
                                   true_lag = NULL, period = 12,
                                   amplitude = 1, autocorr = 0.5,
                                   noise_sd = 1, seed = NULL) {
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
  if (!is.numeric(autocorr) || length(autocorr) != 1 || is.na(autocorr) ||
      abs(autocorr) >= 1) {
    stop("`autocorr` must lie in (-1, 1).", call. = FALSE)
  }

  with_seed(seed, {
    if (is.null(true_lag)) {
      X <- sim_seasonal(n_times, n_vars_x, period, amplitude, autocorr,
                        noise_sd)
      Y <- sim_seasonal(n_times, n_vars_y, period, amplitude, autocorr,
                        noise_sd)
    } else {
      X_full <- sim_seasonal(n_times + true_lag, n_vars_x, period, amplitude,
                             autocorr, noise_sd)
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

sim_seasonal <- function(n, p, period, amplitude, autocorr, noise_sd) {
  t <- seq_len(n)
  phase <- stats::runif(1, 0, 2 * pi)
  base <- vapply(seq_len(p), function(j) {
    amplitude * sin(2 * pi * t / period + phase + (j - 1) * pi / p)
  }, numeric(n))
  matrix(base, n, p) + ar1_noise(n, p, autocorr, noise_sd)
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
