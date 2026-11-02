#' Plot the lag profile of a lagged PROTEST with ggplot2
#'
#' A ggplot2 version of [plot.lagged_protest()]: \eqn{m^2} against lag, with
#' the `alpha` quantile of the permutation null distribution at each lag as
#' a dashed line, filled points for lags whose p-value (adjusted, if
#' `p.adjust.method` was not `"none"`) is below `alpha`, and the optimal lag
#' circled, filled when the global p-value is below `alpha`. The global
#' p-value is shown in the subtitle. Requires the ggplot2 package.
#'
#' @param object A `"lagged_protest"` object.
#' @param alpha Significance level.
#' @param ... Unused.
#' @return A `ggplot` object, which can be modified with further layers.
#' @seealso [plot.lagged_protest()], [as.data.frame.lagged_protest()]
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
#'                               seed = 1)
#' fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 199, seed = 1)
#' library(ggplot2)
#' autoplot(fit)
#' autoplot(fit) + theme_bw()
#' @exportS3Method ggplot2::autoplot
autoplot.lagged_protest <- function(object, alpha = 0.05, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("`autoplot()` requires the ggplot2 package.", call. = FALSE) # nocov
  }
  prof <- as.data.frame(object, alpha = alpha)
  prof$significant <- prof$p.adjusted < alpha
  best <- prof[prof$optimal, ]
  ggplot2::ggplot(prof, ggplot2::aes(x = lag, y = m2)) +
    ggplot2::geom_vline(xintercept = object$optimal_lag, linetype = 3,
                        colour = "grey50") +
    ggplot2::geom_line(ggplot2::aes(y = null.quantile), linetype = 2,
                       colour = "grey50") +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(shape = significant), size = 2,
                        show.legend = FALSE) +
    ggplot2::geom_point(data = best, shape = 21, size = 5,
                        fill = if (object$p.global < alpha) "grey50" else NA) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19)) +
    ggplot2::scale_x_continuous(breaks = prof$lag) +
    ggplot2::labs(x = "Lag", y = expression(italic(m)^2),
                  subtitle = paste0("Optimal lag ", object$optimal_lag,
                                    ", global p = ",
                                    format(object$p.global, digits = 2)))
}

utils::globalVariables(c("lag", "m2", "null.quantile", "significant"))
