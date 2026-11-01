#' @export
print.lagged_protest <- function(x, digits = 3, ...) {
  cat("\nLagged PROTEST\n\n")
  cat("Call: ", deparse(x$call, width.cutoff = 70L), "\n", sep = "")
  scheme <- if (x$permutation == "cyclic") {
    paste0("cyclic shift", if (isTRUE(x$mirror)) " with mirror", "")
  } else {
    "free row permutation"
  }
  cat("Permutation: ", scheme, "; n = ", x$n, " time points\n", sep = "")
  cat("P-value adjustment across lags: ", x$p.adjust.method, "\n\n", sep = "")
  cat("Optimal lag: ", x$optimal_lag,
      "  (m2 = ", format(x$m2, digits = digits),
      ", r = ", format(x$r, digits = digits),
      ", adjusted p = ", format(x$p.value, digits = digits), ")\n\n", sep = "")
  cat("Lag profile:\n")
  print(x$profile, digits = digits, row.names = FALSE)
  if (x$permutation == "cyclic") {
    p_min <- 1 / (max(x$profile$nperm) + 1)
    cat("\nSmallest attainable raw p-value: ",
        format(p_min, digits = digits), "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a lagged PROTEST
#'
#' @param object A `"lagged_protest"` object.
#' @param ... Unused.
#' @return The lag profile (a data frame), with the optimal lag flagged.
#' @export
summary.lagged_protest <- function(object, ...) {
  out <- object$profile
  out$optimal <- out$lag == object$optimal_lag
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
#' points mark lags whose adjusted p-value is below `alpha`.
#'
#' @param x A `"lagged_protest"` object.
#' @param alpha Significance level.
#' @param ... Further arguments passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.lagged_protest <- function(x, alpha = 0.05, ...) {
  prof <- x$profile
  null_q <- vapply(x$null, stats::quantile, numeric(1), probs = alpha,
                   names = FALSE)
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
  invisible(x)
}
