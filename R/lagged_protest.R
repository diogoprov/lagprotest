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
#' Three permutation schemes are available:
#'
#' * `"cyclic"` (default): the full `Y` series is cyclically shifted by
#'   `1, ..., n - 1` positions before lagging, and optionally also
#'   time-reversed (`mirror = TRUE`). This keeps the temporal autocorrelation
#'   of `Y` under the null. All `n - 1` (or `2(n - 1)`) shifts are enumerated,
#'   so `nperm` is ignored, and the smallest attainable p-value is
#'   `1 / (2(n - 1) + 1)` (0.043 for `n = 12` with `mirror = TRUE`). The scheme
#'   assumes the series spans complete cycles, i.e., that its period matches
#'   the number of rows (e.g., 12 monthly samples of an annual cycle).
#'   Because a shift of `Y` combined with lag `k` is equivalent to another
#'   lag, some shifts reproduce the observed alignment at other lags; the
#'   global p-value (see below) is then valid but has low power.
#' * `"phase"`: `nperm` phase-randomized surrogates of `Y` are drawn. Each
#'   surrogate keeps the amplitude spectrum of every column of `Y` (hence its
#'   autocorrelation) and the cross-spectra among its columns, but receives
#'   new Fourier phases, shared by all columns. Surrogates are applied jointly
#'   to all lags and do not reproduce the observed alignment, which makes
#'   this the recommended scheme for the global p-value. It assumes a
#'   stationary `Y`; like the cyclic scheme, it treats the series as
#'   periodic with period `n`.
#' * `"free"`: rows of `Y` are permuted at random within the overlapping
#'   window, as in [vegan::protest()]. This assumes exchangeable rows and can
#'   strongly inflate Type I error when both series are autocorrelated.
#'
#' P-values are computed as `(b + 1) / (N + 1)`, where `b` is the number of
#' permuted \eqn{m^2} values smaller than or equal to the observed one and
#' `N` the number of permutations. Each permutation is applied jointly to all
#' lags.
#'
#' The per-lag p-values (`p.value`) are not corrected for testing several
#' lags, nor for selecting the lag with the smallest \eqn{m^2}. Reporting the
#' raw p-value at the optimal lag therefore inflates Type I error roughly to
#' the rate at which any lag is significant. The global p-value (`p.global`)
#' accounts for both: it compares the observed minimum \eqn{m^2} across lags
#' with the distribution of the minimum across lags under each permutation
#' (a max-statistic test), so it is the p-value of the whole
#' scan-and-select procedure and is the one to report for the optimal lag.
#' With `p.adjust.method = "maxT"`, every lag receives the corresponding
#' single-step max-statistic adjusted p-value, and the adjusted p-value at
#' the optimal lag equals `p.global`. Any method of [stats::p.adjust()] can be
#' requested instead; with few permutations, as in the cyclic scheme,
#' corrections such as Holm may make significance unattainable. Note that
#' \eqn{m^2} values at different lags are based on different numbers of rows
#' (`n - |k|`).
#'
#' @param X,Y Numeric matrices or data frames with rows as ordered time
#'   points and columns as variables, or `dist` objects, which are converted
#'   to principal coordinates. Both must have the same number of rows.
#' @param lags Integer vector of non-negative candidate lags.
#' @param permutation Permutation scheme, `"cyclic"`, `"phase"` or `"free"`.
#' @param nperm Number of surrogates (`permutation = "phase"`) or random
#'   permutations (`permutation = "free"`). Ignored for the cyclic scheme,
#'   whose permutations are enumerated.
#' @param mirror Logical; with `permutation = "cyclic"`, also use the
#'   time-reversed shifts of `Y`, doubling the permutation set.
#' @param bidirectional Logical; if `TRUE`, also test the negative lags
#'   `-lags` (Y leads X).
#' @param scale Logical; z-score the columns of `X` and `Y` (computed on the
#'   full series) before analysis. Ignored for `dist` inputs.
#' @param p.adjust.method Correction of the per-lag p-values for multiple
#'   lags: `"maxT"` for the single-step max-statistic adjustment, or any
#'   method of [stats::p.adjust()]. Defaults to `"none"`.
#' @param seed Optional integer seed for the `"phase"` and `"free"` schemes.
#'   The global random number stream is left untouched.
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
#'     \item{`p.global`}{Max-statistic p-value of the scan across lags,
#'       accounting for testing several lags and selecting the optimal one.}
#'     \item{`n_perm`}{Number of permutations or surrogates at the optimal
#'       lag (the exact count for the cyclic scheme).}
#'     \item{`profile`}{Data frame with `lag`, `n` (rows used), `m2`, `r`,
#'       `p.value`, `p.adjusted` and `nperm` for every lag.}
#'     \item{`residuals`}{Data frame of per-observation Procrustes residuals
#'       at the optimal lag, with the paired row labels of `X` and `Y`.}
#'     \item{`procrustes`}{The [vegan::procrustes()] fit at the optimal lag.}
#'     \item{`null`}{List of permuted \eqn{m^2} values, one element per lag;
#'       element `i` of every vector comes from the same permutation.}
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
#' fit_phase <- lagged_protest(sim$X, sim$Y, lags = 0:4,
#'                             permutation = "phase", seed = 1)
#' fit_phase$p.global
#' plot(fit)
#' residuals(fit)
#' @export
lagged_protest <- function(X, Y, lags = 0:6,
                           permutation = c("cyclic", "phase", "free"),
                           nperm = 999, mirror = TRUE, bidirectional = FALSE,
                           scale = TRUE, p.adjust.method = "none",
                           seed = NULL) {
  call <- match.call()
  permutation <- match.arg(permutation)
  p.adjust.method <- match.arg(p.adjust.method,
                               c(stats::p.adjust.methods, "maxT"))

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
  if (permutation != "cyclic" && (length(nperm) != 1 || nperm < 1)) {
    stop("`nperm` must be a positive integer.", call. = FALSE)
  }

  observed <- vapply(lags, function(k) {
    idx <- lag_index(n, k)
    procrustes_m2(Xc[idx$x, , drop = FALSE], Yc[idx$y, , drop = FALSE])
  }, numeric(1))

  null <- switch(permutation,
    cyclic = null_cyclic(Xc, Yc, lags, mirror),
    phase = with_seed(seed, null_phase(Xc, Yc, lags, nperm)),
    free = with_seed(seed, null_free(Xc, Yc, lags, nperm))
  )
  names(null) <- lags

  n_null <- unname(lengths(null))
  p_raw <- (vapply(seq_along(lags),
                   function(i) sum(null[[i]] <= observed[i]),
                   numeric(1)) + 1) / (n_null + 1)

  # Max-statistic: each permutation is applied jointly to all lags, so the
  # smallest permuted m2 across lags is the null of the scan-and-select
  # procedure itself.
  null_min <- do.call(pmin, unname(null))
  p_max <- (vapply(observed, function(o) sum(null_min <= o), numeric(1)) +
              1) / (length(null_min) + 1)
  p_adj <- if (p.adjust.method == "maxT") {
    p_max
  } else {
    stats::p.adjust(p_raw, method = p.adjust.method)
  }

  profile <- data.frame(
    lag = lags,
    n = n - abs(lags),
    m2 = observed,
    r = sqrt(pmax(0, 1 - observed)),
    p.value = p_raw,
    p.adjusted = p_adj,
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
      p.global = p_max[best],
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
