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

# Free permutation: the rows of the full Y series are permuted at random
# and the permuted series is then lagged, as in the other schemes. Assumes
# rows are exchangeable under the null, which autocorrelated series violate.
# The observed pairing (the identity permutation) and every permuted one are
# thus built in the same way, so the joint distribution of m2 across lags is
# the same for all of them and the max-statistic test is exact. Permuting
# within each lag's window instead (keeping the order of the permuted rows)
# made the permuted m2 more correlated across lags than the observed ones,
# and the global test liberal.
null_free <- function(X, Y, lags, nperm) {
  n <- nrow(Y)
  perms <- lapply(seq_len(nperm), function(i) sample.int(n))
  lapply(lags, function(k) {
    idx <- lag_index(n, k)
    Xk <- X[idx$x, , drop = FALSE]
    vapply(perms, function(pm) {
      pair_stats(Xk, Y[pm[idx$y], , drop = FALSE])
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

# m2 (symmetric Procrustes) and RV coefficient of one pair of configurations
# with the same rows. This is the inner loop of every permutation scheme, so
# it avoids the S3 dispatch, warning handling and sweep() of the generic
# functions. With both configurations centered and scaled to unit trace,
# vegan::procrustes(X, Y, symmetric = TRUE)$ss equals 1 - (sum of the
# singular values of X'Y)^2; padding X with zero columns, as vegan does when
# it has fewer columns than Y, leaves the singular values unchanged. The
# result matches vegan to floating-point rounding (tested).
pair_stats <- function(X, Y) {
  X <- X - rep(colMeans(X), each = nrow(X))
  Y <- Y - rep(colMeans(Y), each = nrow(Y))
  sxy <- crossprod(X, Y)
  d <- svd(sxy, nu = 0, nv = 0)$d
  m2 <- 1 - sum(d)^2 / (sum(X^2) * sum(Y^2))
  rv <- sum(sxy^2) / sqrt(sum(crossprod(X)^2) * sum(crossprod(Y)^2))
  c(m2, rv)
}
