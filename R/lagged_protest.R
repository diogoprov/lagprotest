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
#' Procrustes sum of squares \eqn{m^2} from [vegan::procrustes()]. The
#' optimal lag is chosen by the scan statistic set by `stat` (by default the
#' lag with the smallest per-lag p-value; see below).
#'
#' Three permutation schemes are available:
#'
#' * `"phase"` (default): `nperm` phase-randomized surrogates of `Y` are
#'   drawn. Each surrogate keeps the amplitude spectrum of every column of
#'   `Y` (hence its autocorrelation) and the cross-spectra among its columns,
#'   but receives new Fourier phases, shared by all columns. Surrogates are
#'   applied jointly to all lags and do not reproduce the observed alignment,
#'   which makes this the recommended scheme for the global p-value. It
#'   assumes a stationary `Y` and treats the series as periodic with period
#'   `n`. In simulations, it was liberal for short, strongly autocorrelated
#'   series (AR(1) of 0.7 or more with `n` of 48 or less).
#' * `"cyclic"`: the full `Y` series is cyclically shifted by
#'   `1, ..., n - 1` positions before lagging, and optionally also
#'   time-reversed (`mirror = TRUE`). This keeps the temporal autocorrelation
#'   of `Y` under the null. All `n - 1` (or `2(n - 1)`) shifts are enumerated,
#'   so `nperm` and `seed` are ignored, and the smallest attainable p-value is
#'   `1 / (2(n - 1) + 1)` (0.043 for `n = 12` with `mirror = TRUE`). The scheme
#'   assumes the series spans complete cycles, i.e., that its period matches
#'   the number of rows (e.g., 12 monthly samples of an annual cycle).
#'   Because a shift of `Y` combined with lag `k` is equivalent to another
#'   lag, some shifts reproduce the observed alignment at other lags; the
#'   global p-value (see below) is then valid but has low power.
#' * `"free"`: rows of `Y` are permuted at random, as in
#'   [vegan::protest()], and the permuted series is then lagged like the
#'   observed one. This assumes exchangeable rows and can
#'   strongly inflate Type I error when both series are autocorrelated
#'   (Manly and Navarro Alberto 2020, ch. 10).
#'
#' Both the phase and the cyclic schemes treat the series as periodic with
#' period `n`. When the series contains a seasonal cycle whose length does
#' not divide `n` (e.g., 18 monthly samples of an annual cycle), no scheme is
#' calibrated, because two series sharing a period are Procrustes-concordant
#' whatever their phase. Supplying `period` truncates both series to their
#' first `floor(n / period) * period` rows, so that the analysis spans whole
#' cycles. In simulations, this restored the Type I error of series that
#' already spanned complete cycles, whereas removing the seasonal component
#' by harmonic regression did not.
#'
#' P-values are computed as `(b + 1) / (N + 1)`, where `b` is the number of
#' permuted \eqn{m^2} values smaller than or equal to the observed one and
#' `N` the number of permutations (Manly and Navarro Alberto 2020, p. 57).
#' Each permutation is applied jointly to all lags. The validity of the test
#' does not depend on `N`, only the precision of the p-value does; Manly and
#' Navarro Alberto (2020, pp. 57--59) recommend at least 1000 randomizations
#' for a test at the 5% level and about 5000 at the 1% level. A message is
#' issued when the smallest attainable p-value, `1 / (N + 1)`, is 0.05 or
#' more.
#'
#' Because the p-value is bounded by `N` and conflates strength with sample
#' size, three continuous measures are computed from the same permutations.
#' At each lag, \eqn{z_k = (\bar m^{2*}_k - m^2_k) / s^*_k} standardizes the
#' observed \eqn{m^2} against its permutation mean and standard deviation,
#' and \eqn{a_k = 1 - m^2_k / \bar m^{2*}_k} expresses the concordance
#' relative to its null mean. The effect lag is the lag of the largest
#' \eqn{z} (which weighs lags of different overlap more evenly than the raw
#' \eqn{m^2}), and each permutation selects its own lag in the same way.
#' `Z` is the standardized position of the observed maximum \eqn{z} within
#' the permuted maxima; `A` and the corrected RV coefficient are
#' chance-corrected at the effect lag as
#' \eqn{(obs - \bar{null}) / (1 - \bar{null})}, with the null mean taken
#' over the permutations at their own selected lags. `A` and the corrected RV
#' are near 0 under the null and 1 for perfect concordance; `Z` is near 0
#' under the null. Each permuted \eqn{z} is computed in the same way as the
#' observed one, against the mean and standard deviation of the other values
#' at that lag (the remaining permutations and the observed), so that the
#' observed and the permuted values are treated alike.
#'
#' The mean lag-1 autocorrelation of the columns of each series is returned
#' in `ar1`, both raw and after removing by regression the harmonics of
#' period `n` (and of `period`, when given). It is descriptive only: the
#' raw value is dominated by any seasonal cycle, and the residual value is
#' biased downwards in short series. The phase scheme was liberal for short,
#' strongly autocorrelated series in simulations, so both values help judge
#' whether that caveat applies.
#'
#' The per-lag p-values (`p.value`) are not corrected for testing several
#' lags, nor for selecting the lag with the smallest \eqn{m^2}. Reporting the
#' raw p-value at the optimal lag therefore inflates Type I error roughly to
#' the rate at which any lag is significant. The global p-value (`p.global`)
#' accounts for both: it compares the most extreme observed score across
#' lags with the distribution of the most extreme score across lags under
#' each permutation (a max-statistic test), so it is the p-value of the whole
#' scan-and-select procedure and is the one to report for the optimal lag.
#' Because \eqn{m^2} at neighboring lags is strongly correlated, a Bonferroni
#' level of \eqn{\alpha / m} for \eqn{m} lags would be conservative: the
#' appropriate level lies between \eqn{\alpha / m} and \eqn{\alpha}, and
#' taking the extreme statistic across tests in each randomization, as here,
#' finds it from the data (Manly and Navarro Alberto 2020, pp. 83--84). The
#' same source recommends at least \eqn{100m} randomizations for \eqn{m}
#' tests.
#' With `p.adjust.method = "maxT"`, every lag receives the corresponding
#' single-step adjusted p-value for the chosen `stat`, and the adjusted p-value at
#' the optimal lag equals `p.global`. Any method of [stats::p.adjust()] can be
#' requested instead; with few permutations, as in the cyclic scheme,
#' corrections such as Holm may make significance unattainable.
#'
#' The score compared across lags is set by `stat`:
#'
#' * `"minP"` (default): the per-lag p-value. The global statistic is the
#'   smallest p-value across lags, and its null distribution is obtained by
#'   converting every permuted \eqn{m^2} to a p-value within its own lag,
#'   as its rank among the observed and permuted values divided by
#'   \eqn{N + 1}, exactly as the observed p-value is computed (the minP
#'   procedure of Westfall and Young 1993; Ge et al. 2003). The
#'   optimal lag is the one with the smallest p-value; ties, common with few
#'   permutations, are broken by \eqn{z} and then by \eqn{m^2}.
#' * `"z"`: \eqn{m^2} centered and scaled by the mean and standard deviation
#'   of its own null distribution, with each permuted value standardized
#'   against the others and the observed (see above). The optimal lag is the
#'   one with the largest \eqn{z} and equals the effect lag.
#' * `"raw"`: \eqn{m^2} itself, the statistic of versions up to 0.6. The
#'   optimal lag is the one with the smallest \eqn{m^2}.
#'
#' All three are computed from the same permutations and give valid global
#' tests. They differ in how power is spread across lags, because
#' \eqn{m^2} values at different lags are based on different numbers of
#' rows (`n - |k|`), and \eqn{m^2} is smaller and more variable when fewer
#' rows overlap. Under `"raw"`, the optimal lag therefore tends to fall at
#' the longest lags scanned, even when the series are independent: in 200
#' pairs of independent series simulated with [simulate_lagged_series()]
#' (`n_times = 12`, lags -3 to 3), it fell at \eqn{\pm 3} in about half of
#' the replicates rather than in two of seven, and the test had less power
#' for a true lag of zero. `"minP"` and `"z"` place every lag on the scale
#' of its own null distribution, so the optimal lag is close to uniform
#' across lags under independence. Because the observed and the permuted
#' scores are computed in the same way, all three give exact tests when the
#' permutations are exchangeable with the observed series (free and cyclic
#' schemes). `"minP"` is conservative when the permutations are few relative
#' to the number of lags, because per-lag p-values are discrete, and `"z"`
#' then has more power. The same reasoning motivates selecting among
#' candidates by p-value in conditional inference trees (Hothorn et al.
#' 2006). Because per-lag p-values are discrete, the smallest global p-value
#' that `"minP"` can reach grows with the number of lags: roughly
#' \eqn{(m + 1) / (N + 1)} for \eqn{m} lags and \eqn{N} permutations, so
#' \eqn{N} should be at least \eqn{100m} (Manly and Navarro Alberto 2020,
#' p. 84), and a message is issued when p < 0.05 cannot be reached. With the
#' cyclic scheme, permutations are few, and `"z"` is then preferable.
#' When the optimal lag and the effect lag differ, `print()` says so;
#' candidate lags should be kept within a range justified by the biology of
#' the system.
#'
#' @param X,Y Numeric matrices or data frames with rows as ordered time
#'   points and columns as variables, or `dist` objects, which are converted
#'   to principal coordinates. Both must have the same number of rows.
#' @param lags Integer vector of non-negative candidate lags (default
#'   `0:6`). Every lag must leave at least 3 overlapping rows; long lags rest
#'   on few rows (6 for lag 6 with `n = 12`) and are favored by the smallest
#'   \eqn{m^2} (see Details).
#' @param permutation Permutation scheme, `"phase"`, `"cyclic"` or `"free"`.
#' @param nperm Number of surrogates (`permutation = "phase"`) or random
#'   permutations (`permutation = "free"`). Ignored for the cyclic scheme,
#'   whose permutations are enumerated. Use 999 or more for a test at the 5%
#'   level and about 4999 when the 1% level matters (Manly and Navarro
#'   Alberto 2020, pp. 57--59). With `stat = "minP"`, use at least 100 per
#'   lag tested (for example, 699 or more for lags -3 to 3), otherwise the
#'   global p-value may be unable to fall below 0.05 (see Details).
#' @param mirror Logical; with `permutation = "cyclic"`, also use the
#'   time-reversed shifts of `Y`, doubling the permutation set.
#' @param bidirectional Logical; if `TRUE`, also test the negative lags
#'   `-lags` (Y leads X).
#' @param scale Logical; z-score the columns of `X` and `Y` (computed on the
#'   full series) before analysis. Ignored for `dist` inputs.
#' @param stat Scan statistic used for the global p-value and the optimal
#'   lag: `"minP"` (default), `"z"` or `"raw"` (see Details).
#' @param p.adjust.method Correction of the per-lag p-values for multiple
#'   lags: `"maxT"` for the single-step max-statistic adjustment, or any
#'   method of [stats::p.adjust()]. Defaults to `"none"`.
#' @param period Optional length of the seasonal cycle, in rows. If given,
#'   `X` and `Y` are truncated to their first `floor(n / period) * period`
#'   rows (whole cycles) before analysis; a message reports the rows dropped.
#'   `NULL` (default) uses all rows.
#' @param seed Optional integer seed for the `"phase"` and `"free"` schemes.
#'   The global random number stream is left untouched.
#'
#' @return An object of class `"lagged_protest"`, a list with:
#'   \describe{
#'     \item{`optimal_lag`}{Lag selected by `stat` (see Details).}
#'     \item{`m2`}{\eqn{m^2} at the optimal lag.}
#'     \item{`r`}{Procrustes correlation, \eqn{\sqrt{1 - m^2}}, at the
#'       optimal lag.}
#'     \item{`p.value`}{Raw permutation p-value at the optimal lag.}
#'     \item{`p.adjusted`}{P-value at the optimal lag adjusted with
#'       `p.adjust.method` (equal to `p.value` when `"none"`).}
#'     \item{`p.global`}{P-value of the scan across lags for `stat`,
#'       accounting for testing several lags and selecting the optimal one.}
#'     \item{`n_perm`}{Number of permutations or surrogates, the same at
#'       every lag (the exact count for the cyclic scheme).}
#'     \item{`profile`}{Data frame with `lag`, `n` (rows used), `m2`, `r`,
#'       `p.value`, `p.adjusted`, `nperm`, `z` and `a` for every lag.}
#'     \item{`residuals`}{Data frame of per-observation Procrustes residuals
#'       at the optimal lag, with the paired row labels of `X` and `Y`.}
#'     \item{`procrustes`}{The [vegan::procrustes()] fit at the optimal lag.}
#'     \item{`effect`}{List with the effect `lag` (largest \eqn{z}), `Z`,
#'       `A` and the chance-corrected `RV` coefficient at that lag (see
#'       Details).}
#'     \item{`ar1`}{Data frame with the mean lag-1 autocorrelation of the
#'       columns of `X` and `Y`, `raw` and `residual` after removing the
#'       harmonics of period `n` and `period`.}
#'     \item{`null`}{List of permuted \eqn{m^2} values, one element per lag;
#'       element `i` of every vector comes from the same permutation.}
#'     \item{`null.rv`}{List of permuted RV coefficients, paired with
#'       `null`.}
#'     \item{`permutation`, `stat`, `mirror`, `p.adjust.method`, `n`, `period`,
#'       `call`}{Settings used; `mirror` is `NA` unless the scheme is cyclic,
#'       and `n` is the number of rows after any truncation to whole cycles.}
#'   }
#'
#' @references
#' Ge, Y., Dudoit, S. and Speed, T. P. (2003). Resampling-based multiple
#' testing for microarray data analysis. *Test* 12, 1--77.
#' \doi{10.1007/BF02595811}
#'
#' Hothorn, T., Hornik, K. and Zeileis, A. (2006). Unbiased recursive
#' partitioning: a conditional inference framework. *Journal of
#' Computational and Graphical Statistics* 15, 651--674.
#' \doi{10.1198/106186006X133933}
#'
#' Manly, B. F. J. and Navarro Alberto, J. A. (2020). *Randomization,
#' Bootstrap and Monte Carlo Methods in Biology*, 4th ed. Chapman and
#' Hall/CRC, New York. \doi{10.1201/9780429329203}
#'
#' Westfall, P. H. and Young, S. S. (1993). *Resampling-Based Multiple
#' Testing: Examples and Methods for p-Value Adjustment*. Wiley, New York.
#'
#' @seealso [vegan::protest()], [vegan::procrustes()]
#'
#' @examples
#' sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
#'                               seed = 1)
#' fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, seed = 1)
#' fit
#' fit$p.global
#' fit_cyclic <- lagged_protest(sim$X, sim$Y, lags = 0:4,
#'                              permutation = "cyclic", stat = "z")
#' fit_cyclic$profile
#' plot(fit)
#' residuals(fit)
#' @export
lagged_protest <- function(X, Y, lags = 0:6,
                           permutation = c("phase", "cyclic", "free"),
                           nperm = 999, mirror = TRUE, bidirectional = FALSE,
                           scale = TRUE, stat = c("minP", "z", "raw"),
                           p.adjust.method = "none",
                           period = NULL, seed = NULL) {
  call <- match.call()
  permutation <- match.arg(permutation)
  stat <- match.arg(stat)
  p.adjust.method <- match.arg(p.adjust.method,
                               c(stats::p.adjust.methods, "maxT"))
  check_flag(mirror, "mirror")
  check_flag(bidirectional, "bidirectional")
  check_flag(scale, "scale")
  check_seed(seed)
  if (permutation != "cyclic") check_count(nperm, "nperm")

  if (n_rows(X) != n_rows(Y)) {
    stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  }
  if (!is.null(period)) {
    n_all <- n_rows(X)
    m <- whole_cycle_rows(n_all, period)
    if (m < n_all) {
      message("Using the first ", m, " of ", n_all, " rows (",
              floor(n_all / period), " whole cycle",
              if (floor(n_all / period) > 1) "s", " of period ", period, ").")
      X <- first_rows(X, m)
      Y <- first_rows(Y, m)
    }
  }
  Xc <- as_configuration(X, scale, "X")
  Yc <- as_configuration(Y, scale, "Y")
  n <- nrow(Xc)
  if (is.null(rownames(Xc))) rownames(Xc) <- seq_len(n)
  if (is.null(rownames(Yc))) rownames(Yc) <- seq_len(n)

  if (!is.numeric(lags) || length(lags) == 0 || anyNA(lags) ||
      any(lags < 0) || any(lags != round(lags))) {
    stop("`lags` must be non-negative integers.", call. = FALSE)
  }
  lags <- as.integer(sort(unique(lags)))
  if (bidirectional) lags <- sort(unique(c(-lags, lags)))
  if (any(n - abs(lags) < 3)) {
    stop("Each lag must leave at least 3 overlapping rows (n = ", n, ").",
         call. = FALSE)
  }

  obs_stats <- vapply(lags, function(k) {
    idx <- lag_index(n, k)
    pair_stats(Xc[idx$x, , drop = FALSE], Yc[idx$y, , drop = FALSE])
  }, numeric(2))
  observed <- obs_stats[1, ]

  null_both <- switch(permutation,
    cyclic = null_cyclic(Xc, Yc, lags, mirror),
    phase = with_seed(seed, null_phase(Xc, Yc, lags, nperm)),
    free = with_seed(seed, null_free(Xc, Yc, lags, nperm))
  )
  null <- lapply(null_both, function(m) m[1, ])
  null_rv <- lapply(null_both, function(m) m[2, ])
  names(null) <- names(null_rv) <- lags

  n_null <- unname(lengths(null))
  p_min_lag <- check_min_p(n_null, permutation, mirror)
  cyc <- c(n, period)
  ar1 <- data.frame(
    series = c("X", "Y"),
    raw = c(mean_ar1(Xc), mean_ar1(Yc)),
    residual = c(mean_ar1(Xc, cyc), mean_ar1(Yc, cyc))
  )
  p_raw <- (vapply(seq_along(lags),
                   function(i) sum(null[[i]] <= observed[i]),
                   numeric(1)) + 1) / (n_null + 1)

  # Scan statistic: each permutation is applied jointly to all lags, so the
  # most extreme permuted score across lags is the null of the
  # scan-and-select procedure itself. Scores are on a common scale across
  # lags for "minP" and "z", and on the m2 scale for "raw".
  eff <- effect_measures(lags, observed, obs_stats[2, ], null, null_rv)
  sc <- scan_scores(stat, observed, p_raw, null)
  null_min <- apply(sc$null, 1, min)
  # Only one message: if no per-lag p reaches 0.05, neither does minP.
  if (stat == "minP" && p_min_lag < 0.05) {
    check_min_p_scan(null_min, length(lags), permutation)
  }
  p_max <- (vapply(sc$observed, function(o) sum(null_min <= o), numeric(1)) +
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
    nperm = n_null,
    z = eff$per_lag$z,
    a = eff$per_lag$a
  )

  # Optimal lag: most extreme score; ties (frequent for minP) are broken by
  # the standardized m2 and then by the raw m2.
  best <- order(sc$observed, -eff$per_lag$z, observed)[1]
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
      effect = eff$summary,
      ar1 = ar1,
      null = null,
      null.rv = null_rv,
      permutation = permutation,
      stat = stat,
      mirror = if (permutation == "cyclic") mirror else NA,
      p.adjust.method = p.adjust.method,
      n = n,
      period = period,
      call = call
    ),
    class = "lagged_protest"
  )
}
