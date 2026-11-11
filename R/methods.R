#' @export
print.lagged_protest <- function(x, digits = 3, ...) {
  cat("\nLagged PROTEST\n\n")
  cat("Call: ", deparse(x$call, width.cutoff = 70L), "\n", sep = "")
  scheme <- switch(x$permutation,
    cyclic = paste0("cyclic shift", if (isTRUE(x$mirror)) " with mirror"),
    phase = "phase-randomized surrogates",
    free = "free row permutation"
  )
  cat("Permutation: ", scheme, "; n = ", x$n, " time points",
      if (!is.null(x$period)) paste0(" (whole cycles of period ", x$period,
                                     ")"),
      "\n", sep = "")
  cat("P-value adjustment across lags: ", x$p.adjust.method, "\n\n", sep = "")
  adjusted <- x$p.adjust.method != "none"
  cat("Optimal lag: ", x$optimal_lag,
      "  (m2 = ", format(x$m2, digits = digits),
      ", r = ", format(x$r, digits = digits),
      ", p = ", format(x$p.value, digits = digits),
      if (adjusted) paste0(", adjusted p = ",
                           format(x$p.adjusted, digits = digits)),
      ")\n", sep = "")
  cat("Global p across lags (max-statistic): ",
      format(x$p.global, digits = digits), "\n", sep = "")
  e <- x$effect
  if (!is.null(e) && !is.na(e$lag)) {
    cat("Effect at lag ", e$lag, " (largest z): Z = ",
        format(e$Z, digits = digits), ", A = ", format(e$A, digits = digits),
        ", corrected RV = ", format(e$RV, digits = digits), "\n", sep = "")
  }
  if (!is.null(x$ar1)) {
    cat("Mean lag-1 autocorrelation, raw (residual after removing the cycle): X ",
        format(x$ar1$raw[1], digits = 2), " (",
        format(x$ar1$residual[1], digits = 2), "), Y ",
        format(x$ar1$raw[2], digits = 2), " (",
        format(x$ar1$residual[2], digits = 2), ")\n", sep = "")
  }
  cat("\n")
  cat("Lag profile:\n")
  profile <- x$profile
  if (!adjusted) profile$p.adjusted <- NULL
  print(profile, digits = digits, row.names = FALSE)
  if (x$permutation == "cyclic") {
    p_min <- 1 / (max(x$profile$nperm) + 1)
    cat("\nSmallest attainable raw p-value: ",
        format(p_min, digits = digits), "\n", sep = "")
  }
  if (!adjusted) {
    cat("Per-lag p-values are not corrected for testing several lags",
        "or for selecting\nthe optimal lag; report the global p for the",
        "optimal lag.\n")
  }
  invisible(x)
}

#' Summarize a lagged PROTEST
#'
#' @param object A `"lagged_protest"` object.
#' @param ... Unused.
#' @return The lag profile (a data frame), with the optimal lag flagged.
#'   The scan-level measures are attached as attributes: `p.global`,
#'   `effect` (lag of the largest z, `Z`, `A` and corrected `RV`) and `ar1`
#'   (raw and residual lag-1 autocorrelation of `X` and `Y`).
#' @export
summary.lagged_protest <- function(object, ...) {
  out <- object$profile
  out$optimal <- out$lag == object$optimal_lag
  attr(out, "p.global") <- object$p.global
  attr(out, "effect") <- object$effect
  attr(out, "ar1") <- object$ar1
  out
}

#' Lag profile of a lagged PROTEST as a data frame
#'
#' Returns the lag profile in long (tidy) form, one row per lag, ready for
#' tables or plotting with other tools.
#'
#' @param x A `"lagged_protest"` object.
#' @param row.names,optional Unused; for compatibility with the generic.
#' @param alpha Level of the permutation null quantile in `null.quantile`.
#' @param ... Unused.
#' @return A data frame with the columns of `x$profile` (`lag`, `n`, `m2`,
#'   `r`, `p.value`, `p.adjusted`, `nperm`, `z`, `a`), plus `null.quantile`, the
#'   `alpha` quantile of the permuted \eqn{m^2} at each lag (unadjusted),
#'   and `optimal`, which flags the optimal lag.
#' @seealso [plot.lagged_protest()], [autoplot.lagged_protest()]
#' @examples
#' sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
#'                               seed = 1)
#' fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 199, seed = 1)
#' as.data.frame(fit)
#' @export
as.data.frame.lagged_protest <- function(x, row.names = NULL,
                                         optional = FALSE, alpha = 0.05,
                                         ...) {
  out <- x$profile
  out$null.quantile <- vapply(x$null, stats::quantile, numeric(1),
                              probs = alpha, names = FALSE)
  out$optimal <- out$lag == x$optimal_lag
  out
}

#' Per-observation Procrustes residuals at the optimal lag
#'
#' @param object A `"lagged_protest"` object.
#' @param ... Unused.
#' @return A data frame with the paired row labels of `X` and `Y` and the
#'   residual of each pair after superimposition.
#' @export
residuals.lagged_protest <- function(object, ...) {
  object$residuals
}

#' Plot the lag profile of a lagged PROTEST
#'
#' Plots \eqn{m^2} against lag. The dashed line is the `alpha` quantile of
#' the permutation null distribution at each lag (unadjusted); observed
#' values below it are significant at that level before correction. Filled
#' points mark lags whose p-value (adjusted, if `p.adjust.method` was not
#' `"none"`) is below `alpha`. With `global = TRUE`, the optimal lag is
#' circled, filled when the global p-value is below `alpha`, and the global
#' p-value is printed above the plot; this is the p-value to report for the
#' optimal lag.
#'
#' @param x A `"lagged_protest"` object.
#' @param alpha Significance level.
#' @param global Logical; mark the optimal lag with its global p-value.
#' @param ... Further arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.lagged_protest <- function(x, alpha = 0.05, global = TRUE, ...) {
  prof <- as.data.frame(x, alpha = alpha)
  null_q <- prof$null.quantile
  args <- utils::modifyList(
    list(x = prof$lag, y = prof$m2, type = "n",
         ylim = range(c(prof$m2, null_q)),
         xlab = "Lag", ylab = expression(italic(m)^2), xaxt = "n"),
    list(...)
  )
  do.call(graphics::plot, args)
  graphics::axis(1, at = prof$lag)
  graphics::lines(prof$lag, null_q, lty = 2, col = "grey50")
  graphics::lines(prof$lag, prof$m2)
  graphics::points(prof$lag, prof$m2,
                   pch = ifelse(prof$p.adjusted < alpha, 19, 1))
  graphics::abline(v = x$optimal_lag, lty = 3, col = "grey50")
  if (global) {
    best <- prof$lag == x$optimal_lag
    graphics::points(x$optimal_lag, prof$m2[best], pch = 21, cex = 2,
                     bg = if (x$p.global < alpha) "grey50" else NA)
    graphics::mtext(paste0("Optimal lag ", x$optimal_lag, ", global p = ",
                           format(x$p.global, digits = 2)),
                    side = 3, line = 0.25, adj = 1, cex = 0.8)
  }
  invisible(x)
}
