# Convert input to a numeric matrix of configuration scores.
# dist objects are ordinated by PCoA (all axes with positive eigenvalues).
as_configuration <- function(M, arg) {
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
  M
}

# z-score the columns of a configuration.
scale_columns <- function(M, arg, preprocessed = FALSE) {
  sds <- apply(M, 2, stats::sd)
  if (any(sds <= 1e-12 * max(1, sds))) {
    stop("`", arg, "` has constant columns",
         if (preprocessed) " after preprocessing (`deseason`, `prewhiten`)",
         "; cannot scale.", call. = FALSE)
  }
  M <- scale(M)
  attr(M, "scaled:center") <- attr(M, "scaled:scale") <- NULL
  M
}

# Subtract from every column the mean of its position in the cycle (e.g.,
# each calendar month), over the whole series. Rows are assumed to start at
# position 1 of the cycle and to span whole cycles.
remove_cycle_means <- function(M, period) {
  pos <- factor((seq_len(nrow(M)) - 1) %% period)
  out <- M - apply(M, 2, stats::ave, pos)
  dim(out) <- dim(M)
  dimnames(out) <- dimnames(M)
  out
}

# Prais-Winsten AR(1) prewhitening of every column, which keeps all n rows:
# the first row is scaled by sqrt(1 - rho^2) and the others replaced by
# x_t - rho x_{t-1}, with rho the lag-1 autocorrelation of the centered
# column. Dropping the first row instead breaks the calibration after
# seasonal means are removed.
prais_winsten <- function(M) {
  n <- nrow(M)
  M <- sweep(M, 2, colMeans(M))
  rho <- apply(M, 2, function(x) {
    ss <- sum(x^2)
    if (ss < 1e-12) 0 else sum(x[-1] * x[-n]) / ss
  })
  out <- M
  out[1, ] <- sqrt(1 - rho^2) * M[1, ]
  out[-1, ] <- M[-1, , drop = FALSE] -
    sweep(M[-n, , drop = FALSE], 2, rho, "*")
  list(series = out, rho = rho)
}

# Argument checks with messages naming the argument.
check_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(x)
}

is_whole <- function(x, min) {
  is.numeric(x) && length(x) == 1 && is.finite(x) && x == round(x) &&
    x >= min
}

check_count <- function(x, arg, min = 1) {
  if (!is_whole(x, min)) {
    stop("`", arg, "` must be a single whole number of at least ", min, ".",
         call. = FALSE)
  }
  invisible(x)
}

check_seed <- function(seed) {
  if (!is.null(seed) && !is_whole(seed, -Inf)) {
    stop("`seed` must be NULL or a single whole number.", call. = FALSE)
  }
  invisible(seed)
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

# Mean lag-1 autocorrelation of the columns of a configuration. With
# `periods`, each column is first regressed on the harmonics of those periods
# (in rows), so that the value describes the residual noise rather than the
# seasonal cycle. In short series this residual estimate is strongly biased
# downwards, so it is reported for description only.
mean_ar1 <- function(M, periods = NULL) {
  if (length(periods) > 0) {
    t <- seq_len(nrow(M))
    D <- do.call(cbind, c(list(1), lapply(unique(periods), function(p)
      cbind(cos(2 * pi * t / p), sin(2 * pi * t / p)))))
    D <- D[, qr(D)$pivot[seq_len(qr(D)$rank)], drop = FALSE]
    M <- M - D %*% qr.solve(D, M)
  }
  mean(apply(M, 2, function(x) {
    x <- x - mean(x)
    if (sum(x^2) < 1e-12) return(NA_real_)
    sum(x[-1] * x[-length(x)]) / sum(x^2)
  }), na.rm = TRUE)
}

# Continuous measures of evidence and strength of the scan, from the same
# permutations as the p-values. `obs` and `rv` are the observed m2 and RV per
# lag; `null` and `null_rv` lists (one vector per lag) of permuted values,
# element i of every vector from the same permutation.
#   z_k = (mean* - m2_k) / sd*        standardized concordance at lag k
#   a_k = 1 - m2_k / mean*            concordance relative to its null mean
# The lag is chosen by the largest z; each permutation chooses its own lag in
# the same way. Z is the standardized position of the observed maximum z
# within the permuted maxima; A and the corrected RV are chance-corrected at
# that lag: (obs - mean*) / (1 - mean*), with mean* over the permutations,
# each evaluated at its own chosen lag.
effect_measures <- function(lags, obs, rv, null, null_rv) {
  M <- do.call(cbind, unname(null))
  R <- do.call(cbind, unname(null_rv))
  # Each value is standardized by the other B values at its lag (leave one
  # out), so the observed z uses the permutations only and every permuted z
  # is computed in the same way, with the observed among the others.
  st <- loo_standardize(obs, M)
  mu <- colMeans(M)
  z <- -st$observed
  a <- 1 - obs / mu
  zn <- -st$null
  an <- 1 - sweep(M, 2, mu, "/")
  rows <- seq_len(nrow(M))
  kz <- which.max(z)
  knz <- if (anyNA(zn)) rep(NA_integer_, nrow(M)) else
    max.col(zn, ties.method = "first")
  correct <- function(o, nl) {
    m <- mean(nl)
    if (!is.finite(m) || m >= 1) return(NA_real_)
    (o - m) / (1 - m)
  }
  if (length(kz) == 0 || anyNA(knz)) {
    return(list(per_lag = data.frame(z = z, a = a),
                summary = list(lag = NA_integer_, Z = NA_real_, A = NA_real_,
                               RV = NA_real_)))
  }
  zmax <- zn[cbind(rows, knz)]
  sz <- stats::sd(zmax)
  list(
    per_lag = data.frame(z = z, a = a),
    summary = list(
      lag = lags[kz],
      Z = if (is.finite(sz) && sz > 0) (z[kz] - mean(zmax)) / sz else NA_real_,
      A = correct(a[kz], an[cbind(rows, knz)]),
      RV = correct(rv[kz], R[cbind(rows, knz)])
    )
  )
}

# Standardizes each of the B + 1 values at every lag (the observed m2 and
# its B permuted values) by the mean and SD of the other B values. The
# observed is thus standardized by its permutation distribution, and each
# permuted value by the rest of the permutations plus the observed. Returns
# the standardized observed (vector) and permuted (B x lags) values; NA
# where the other values have no variance.
loo_standardize <- function(observed, M) {
  A <- rbind(observed, M)
  B <- nrow(M)
  if (B < 2) {
    A[] <- NA_real_
    return(list(observed = A[1, ], null = A[-1, , drop = FALSE]))
  }
  A0 <- sweep(A, 2, colMeans(A))
  S <- colSums(A0)
  Q <- colSums(A0^2)
  mu <- sweep(-A0, 2, S, "+") / B
  v <- (sweep(-A0^2, 2, Q, "+") - B * mu^2) / (B - 1)
  # Variances at rounding level relative to the spread of the column are
  # treated as zero (no spread among the other values).
  v[v <= 1e-10 * rep(Q, each = B + 1) / B] <- NA_real_
  s <- sqrt(v)
  Z <- (A0 - mu) / s
  list(observed = Z[1, ], null = Z[-1, , drop = FALSE])
}

# Per-lag scores of the scan statistic, observed and permuted (one row per
# permutation), oriented so that smaller is more extreme. Every score is a
# function of the B + 1 values at that lag (observed and permuted) that
# treats them alike, so under the null the observed and the permuted scores
# are exchangeable and the global test is exact.
#   raw:  m2 itself;
#   z:    m2 centered and scaled by the mean and SD of the other B values
#         at that lag (leave one out);
#   minP: per-lag p-value, the rank (ties counted as "<=") among the B + 1
#         values divided by B + 1 (Westfall and Young 1993; Ge et al. 2003);
#         for the observed, this is (count + 1) / (B + 1).
scan_scores <- function(stat, observed, p_raw, null) {
  M <- do.call(cbind, unname(null))
  B <- nrow(M)
  switch(stat,
    raw = list(observed = observed, null = M),
    z = {
      st <- loo_standardize(observed, M)
      if (anyNA(st$observed) || anyNA(st$null)) {
        stop("The null distribution of m2 has no variance at some lag, so ",
             "`stat = \"z\"` is undefined; use `stat = \"minP\"` or ",
             "`\"raw\"`.", call. = FALSE)
      }
      st
    },
    minP = {
      P <- vapply(seq_len(ncol(M)), function(k) {
        rank(c(observed[k], M[, k]), ties.method = "max")[-1]
      }, numeric(B))
      list(observed = p_raw, null = matrix(P, nrow = B) / (B + 1))
    }
  )
}

# Message when the minP global p-value cannot fall below 0.05: the smallest
# attainable value is reached when the observed p is 1 / (B + 1) at some
# lag, and is then set by the permutations that are the most extreme at any
# lag.
check_min_p_scan <- function(null_min, n_lags, permutation) {
  B <- length(null_min)
  p_floor <- (sum(null_min <= 1 / (B + 1)) + 1) / (B + 1)
  if (p_floor >= 0.05) {
    message("With `stat = \"minP\"`, the smallest attainable global p-value ",
            "is ", format(p_floor, digits = 3), " (", B, " permutations, ",
            n_lags, " lags), so p < 0.05 cannot be reached; ",
            if (permutation == "cyclic") "use `stat = \"z\"`" else
              paste0("use at least ", 100 * n_lags, " permutations or ",
                     "`stat = \"z\"`"), ".")
  }
  invisible(p_floor)
}

# Message when no p-value below 0.05 can be reached with the permutations.
check_min_p <- function(n_null, permutation, mirror = TRUE) {
  p_min <- 1 / (min(n_null) + 1)
  if (p_min >= 0.05) {
    advice <- if (permutation != "cyclic") {
      "increase `nperm`"
    } else if (isTRUE(mirror)) {
      "consider the phase scheme"
    } else {
      "consider `mirror = TRUE` or the phase scheme"
    }
    message("The smallest attainable p-value is ", format(p_min, digits = 3),
            " (", min(n_null), " permutations), so p < 0.05 cannot be ",
            "reached; ", advice, ".")
  }
  invisible(p_min)
}
