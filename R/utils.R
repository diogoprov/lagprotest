# Convert input to a numeric matrix of configuration scores.
# dist objects are ordinated by PCoA (all axes with positive eigenvalues);
# matrices and data frames are optionally z-scored column-wise.
as_configuration <- function(M, scale, arg) {
  if (inherits(M, "dist")) {
    n <- attr(M, "Size")
    pcoa <- suppressWarnings(stats::cmdscale(M, k = n - 1, eig = TRUE))
    keep <- pcoa$eig[seq_len(ncol(pcoa$points))] > sqrt(.Machine$double.eps)
    out <- pcoa$points[, keep, drop = FALSE]
    rownames(out) <- attr(M, "Labels")
    return(out)
  }
  M <- as.matrix(M)
  if (!is.numeric(M)) {
    stop("`", arg, "` must be numeric (matrix, data frame or dist).",
         call. = FALSE)
  }
  if (anyNA(M)) stop("`", arg, "` contains missing values.", call. = FALSE)
  if (scale) {
    sds <- apply(M, 2, stats::sd)
    if (any(sds == 0)) {
      stop("`", arg, "` has constant columns; cannot scale.", call. = FALSE)
    }
    M <- scale(M)
    attr(M, "scaled:center") <- attr(M, "scaled:scale") <- NULL
  }
  M
}

n_rows <- function(M) {
  if (inherits(M, "dist")) attr(M, "Size") else NROW(M)
}

# Number of rows spanning whole cycles of length period.
whole_cycle_rows <- function(n, period) {
  if (!is.numeric(period) || length(period) != 1 || is.na(period) ||
      period < 2) {
    stop("`period` must be a single number of at least 2.", call. = FALSE)
  }
  m <- floor(floor(n / period) * period)
  if (m < 3) {
    stop("`period` (", period, ") leaves fewer than 3 rows of whole cycles ",
         "(n = ", n, ").", call. = FALSE)
  }
  m
}

# Keep the first m rows; dist objects are subset as distance matrices.
first_rows <- function(M, m) {
  if (inherits(M, "dist")) {
    stats::as.dist(as.matrix(M)[seq_len(m), seq_len(m)])
  } else {
    M[seq_len(m), , drop = FALSE]
  }
}

# Row indices pairing X and Y at lag k.
# k > 0: X leads Y (X at t paired with Y at t + k).
# k < 0: Y leads X (Y at t paired with X at t + |k|).
lag_index <- function(n, k) {
  if (k >= 0) {
    list(x = seq_len(n - k), y = seq.int(1 + k, n))
  } else {
    list(x = seq.int(1 - k, n), y = seq_len(n + k))
  }
}

# Symmetric Procrustes statistic m^2. Unequal column numbers are padded
# with zeros by vegan::procrustes(), which warns; that is expected here.
procrustes_m2 <- function(X, Y) {
  suppressWarnings(vegan::procrustes(X, Y, symmetric = TRUE)$ss)
}

# Cyclic shift of the rows of M by s positions.
cyclic_shift <- function(M, s) {
  n <- nrow(M)
  M[((seq_len(n) - 1L + s) %% n) + 1L, , drop = FALSE]
}

# Run an expression with a given seed without altering the global RNG state.
with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  code
}
