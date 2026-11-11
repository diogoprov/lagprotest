# Null distributions of m^2 (and of the RV coefficient) for each lag.
# Each function returns a list with one 2-row matrix per lag (m2, RV), one
# column per permutation.

# Cyclic shift: the full Y series is shifted (and, optionally, time-reversed)
# before lagging, so every permuted series keeps its autocorrelation and
# the same n-row period. Shifting the already-truncated window instead would
# treat it as a cycle of length n - k and break this equivalence.
# The permutation set is enumerated exactly: n - 1 shifts, doubled to
# 2(n - 1) when mirror = TRUE.
null_cyclic <- function(X, Y, lags, mirror) {
  n <- nrow(Y)
  shifted <- lapply(seq_len(n - 1L), function(s) cyclic_shift(Y, s))
  if (mirror) {
    Y_rev <- Y[n:1, , drop = FALSE]
    shifted <- c(shifted,
                 lapply(seq_len(n - 1L), function(s) cyclic_shift(Y_rev, s)))
  }
  lapply(lags, function(k) {
    idx <- lag_index(n, k)
    Xk <- X[idx$x, , drop = FALSE]
    vapply(shifted, function(Ys) pair_stats(Xk, Ys[idx$y, , drop = FALSE]),
           numeric(2))
  })
}

# Free permutation: rows of Y are permuted at random within the
# overlapping window of each lag, as in vegan::protest(). Assumes rows are
# exchangeable under the null, which autocorrelated series violate.
# One permutation of 1..n is drawn per iteration and restricted to each
# window (keeping the order of the window's rows), so the same permutation
# is applied jointly to all lags, as the max-statistic requires.
null_free <- function(X, Y, lags, nperm) {
  n <- nrow(Y)
  perms <- lapply(seq_len(nperm), function(i) sample.int(n))
  lapply(lags, function(k) {
    idx <- lag_index(n, k)
    Xk <- X[idx$x, , drop = FALSE]
    vapply(perms, function(pm) {
      rows <- pm[pm %in% idx$y]
      pair_stats(Xk, Y[rows, , drop = FALSE])
    }, numeric(2))
  })
}

# Phase randomization: each surrogate of Y keeps the amplitude spectrum of
# every column but draws new Fourier phases, shared by all columns so that
# the cross-spectra among the columns of Y are also kept. Surrogates are new
# series rather than shifts of the observed one, so, unlike the cyclic
# scheme, no surrogate reproduces the observed alignment at another lag.
null_phase <- function(X, Y, lags, nperm) {
  n <- nrow(Y)
  surrogates <- lapply(seq_len(nperm), function(i) phase_surrogate(Y))
  lapply(lags, function(k) {
    idx <- lag_index(n, k)
    Xk <- X[idx$x, , drop = FALSE]
    vapply(surrogates,
           function(Ys) pair_stats(Xk, Ys[idx$y, , drop = FALSE]),
           numeric(2))
  })
}

# The mean (frequency 0) and, for even n, the Nyquist component are kept;
# each positive frequency is rotated by a random phase and its conjugate by
# the opposite phase, so the surrogate is real.
phase_surrogate <- function(Y) {
  n <- nrow(Y)
  h <- seq_len(floor((n - 1) / 2)) + 1L
  phases <- stats::runif(length(h), 0, 2 * pi)
  rot <- rep(1 + 0i, n)
  rot[h] <- exp(1i * phases)
  rot[n - h + 2L] <- exp(-1i * phases)
  out <- Re(stats::mvfft(stats::mvfft(Y) * rot, inverse = TRUE)) / n
  dimnames(out) <- dimnames(Y)
  out
}

# m2 (symmetric Procrustes) and RV coefficient of one pair of configurations.
pair_stats <- function(X, Y) c(procrustes_m2(X, Y), rv_coef(X, Y))
