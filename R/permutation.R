# Null distributions of m^2 for each lag.
# Both functions return a list with one numeric vector per lag.

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
    vapply(shifted, function(Ys) procrustes_m2(Xk, Ys[idx$y, , drop = FALSE]),
           numeric(1))
  })
}

# Free permutation: rows of Y are permuted at random within the
# overlapping window of each lag, as in vegan::protest(). Assumes rows are
# exchangeable under the null, which autocorrelated series violate.
null_free <- function(X, Y, lags, nperm) {
  n <- nrow(Y)
  lapply(lags, function(k) {
    idx <- lag_index(n, k)
    Xk <- X[idx$x, , drop = FALSE]
    Yk <- Y[idx$y, , drop = FALSE]
    nk <- nrow(Yk)
    vapply(seq_len(nperm),
           function(i) procrustes_m2(Xk, Yk[sample.int(nk), , drop = FALSE]),
           numeric(1))
  })
}
