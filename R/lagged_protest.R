#' Lagged Procrustean randomization test
#'
#' Scans Procrustes concordance between two multivariate series across
#' candidate temporal lags and tests it by permutation at each lag.
#'
#' Rows of `X` and `Y` are ordered time points (e.g., monthly centroids) and
#' must correspond one to one. At lag `k > 0`, `X` at time `t` is paired with
#' `Y` at time `t + k` (X leads Y), and the series are truncated to the
#' overlapping window of `n - k` rows. Negative lags (Y leads X) are tested
#' when `bidirectional = TRUE`. The statistic at each lag is the symmetric
#' Procrustes sum of squares \eqn{m^2} from [vegan::procrustes()]; the lag of
#' maximal concordance is the one with the smallest \eqn{m^2}.
#'
#' Two permutation schemes are available:
#'
#' * `"cyclic"` (default): the full `Y` series is cyclically shifted by
#'   `1, ..., n - 1` positions before lagging, and optionally also
#'   time-reversed (`mirror = TRUE`). This keeps the temporal autocorrelation
#'   of `Y` under the null. All `n - 1` (or `2(n - 1)`) shifts are enumerated,
#'   so `nperm` is ignored, and the smallest attainable p-value is
#'   `1 / (2(n - 1) + 1)` (0.043 for `n = 12` with `mirror = TRUE`). The scheme
#'   assumes the series spans complete cycles, i.e., that its period matches
#'   the number of rows (e.g., 12 monthly samples of an annual cycle).
#' * `"free"`: rows of `Y` are permuted at random within the overlapping
#'   window, as in [vegan::protest()]. This assumes exchangeable rows and can
#'   strongly inflate Type I error when both series are autocorrelated.
#'
#' P-values are computed as `(b + 1) / (N + 1)`, where `b` is the number of
#' permuted \eqn{m^2} values smaller than or equal to the observed one and
#' `N` the number of permutations. By default (`p.adjust.method = "none"`)
#' they are not corrected for testing several lags, nor for selecting the lag
#' with the smallest \eqn{m^2}; the p-value at the optimal lag should therefore
#' be read together with the whole lag profile. Any method of
#' [stats::p.adjust()] can be requested instead, and both raw and adjusted
#' values are returned. With few permutations, as in the cyclic scheme,
#' corrections such as Holm may make significance unattainable. Note that
#' \eqn{m^2} values at different lags are based on different numbers of rows
#' (`n - |k|`).
#'
#' @param X,Y Numeric matrices or data frames with rows as ordered time
#'   points and columns as variables, or `dist` objects, which are converted
#'   to principal coordinates. Both must have the same number of rows.
#' @param lags Integer vector of non-negative candidate lags.
#' @param permutation Permutation scheme, `"cyclic"` or `"free"`.
#' @param nperm Number of random permutations per lag when
#'   `permutation = "free"`.
#' @param mirror Logical; with `permutation = "cyclic"`, also use the
#'   time-reversed shifts of `Y`, doubling the permutation set.
#' @param bidirectional Logical; if `TRUE`, also test the negative lags
#'   `-lags` (Y leads X).
#' @param scale Logical; z-score the columns of `X` and `Y` (computed on the
#'   full series) before analysis. Ignored for `dist` inputs.
#' @param p.adjust.method Correction for multiple lags, passed to
#'   [stats::p.adjust()]. Defaults to `"none"`.
#' @param seed Optional integer seed for `permutation = "free"`. The global
#'   random number stream is left untouched.
#'
#' @return An object of class `"lagged_protest"`, a list with:
#'   \describe{
#'     \item{`optimal_lag`}{Lag with the smallest \eqn{m^2}.}
#'     \item{`m2`}{\eqn{m^2} at the optimal lag.}
#'     \item{`r`}{Procrustes correlation, \eqn{\sqrt{1 - m^2}}, at the
#'       optimal lag.}
#'     \item{`p.value`}{Raw permutation p-value at the optimal lag.}
#'     \item{`p.adjusted`}{P-value at the optimal lag adjusted with
#'       `p.adjust.method` (equal to `p.value` when `"none"`).}
#'     \item{`n_perm`}{Number of permutations at the optimal lag (the exact
#'       count for the cyclic scheme).}
#'     \item{`profile`}{Data frame with `lag`, `n` (rows used), `m2`, `r`,
#'       `p.value`, `p.adjusted` and `nperm` for every lag.}
#'     \item{`residuals`}{Data frame of per-observation Procrustes residuals
#'       at the optimal lag, with the paired row labels of `X` and `Y`.}
#'     \item{`procrustes`}{The [vegan::procrustes()] fit at the optimal lag.}
#'     \item{`null`}{List of permuted \eqn{m^2} values, one element per lag.}
#'     \item{`permutation`, `p.adjust.method`, `n`, `call`}{Settings used.}
#'   }
#'
#' @seealso [vegan::protest()], [vegan::procrustes()]
#'
#' @examples
#' sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
#'                               seed = 1)
#' fit <- lagged_protest(sim$X, sim$Y, lags = 0:4)
#' fit
#' plot(fit)
#' residuals(fit)
#' @export
lagged_protest <- function(X, Y, lags = 0:6,
                           permutation = c("cyclic", "free"),
                           nperm = 999, mirror = TRUE, bidirectional = FALSE,
                           scale = TRUE, p.adjust.method = "none",
                           seed = NULL) {
  call <- match.call()
  permutation <- match.arg(permutation)
  p.adjust.method <- match.arg(p.adjust.method, stats::p.adjust.methods)

  Xc <- as_configuration(X, scale, "X")
  Yc <- as_configuration(Y, scale, "Y")
  n <- nrow(Xc)
  if (nrow(Yc) != n) {
    stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  }
  if (is.null(rownames(Xc))) rownames(Xc) <- seq_len(n)
  if (is.null(rownames(Yc))) rownames(Yc) <- seq_len(n)

  if (!is.numeric(lags) || anyNA(lags) || any(lags < 0) ||
      any(lags != round(lags))) {
    stop("`lags` must be non-negative integers.", call. = FALSE)
  }
  lags <- as.integer(sort(unique(lags)))
  if (bidirectional) lags <- sort(unique(c(-lags, lags)))
  if (any(n - abs(lags) < 3)) {
    stop("Each lag must leave at least 3 overlapping rows (n = ", n, ").",
         call. = FALSE)
  }
  if (permutation == "free" && (length(nperm) != 1 || nperm < 1)) {
    stop("`nperm` must be a positive integer.", call. = FALSE)
  }

  observed <- vapply(lags, function(k) {
    idx <- lag_index(n, k)
    procrustes_m2(Xc[idx$x, , drop = FALSE], Yc[idx$y, , drop = FALSE])
  }, numeric(1))

  null <- switch(permutation,
    cyclic = null_cyclic(Xc, Yc, lags, mirror),
    free = with_seed(seed, null_free(Xc, Yc, lags, nperm))
  )
  names(null) <- lags

  n_null <- unname(lengths(null))
  p_raw <- (vapply(seq_along(lags),
                   function(i) sum(null[[i]] <= observed[i]),
                   numeric(1)) + 1) / (n_null + 1)

  profile <- data.frame(
    lag = lags,
    n = n - abs(lags),
    m2 = observed,
    r = sqrt(pmax(0, 1 - observed)),
    p.value = p_raw,
    p.adjusted = stats::p.adjust(p_raw, method = p.adjust.method),
    nperm = n_null
  )

  best <- which.min(observed)
  k_best <- lags[best]
  idx <- lag_index(n, k_best)
  fit <- suppressWarnings(
    vegan::procrustes(Xc[idx$x, , drop = FALSE], Yc[idx$y, , drop = FALSE],
                      symmetric = TRUE)
  )
  res <- data.frame(
    x_row = rownames(Xc)[idx$x],
    y_row = rownames(Yc)[idx$y],
    residual = unname(stats::residuals(fit))
  )

  structure(
    list(
      optimal_lag = k_best,
      m2 = observed[best],
      r = profile$r[best],
      p.value = profile$p.value[best],
      p.adjusted = profile$p.adjusted[best],
      n_perm = profile$nperm[best],
      profile = profile,
      residuals = res,
      procrustes = fit,
      null = null,
      permutation = permutation,
      mirror = if (permutation == "cyclic") mirror else NA,
      p.adjust.method = p.adjust.method,
      n = n,
      call = call
    ),
    class = "lagged_protest"
  )
}
