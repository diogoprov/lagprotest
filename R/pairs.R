#' Lagged PROTEST for several pairs of data sets
#'
#' Runs [lagged_protest()] on pairs of data sets from a named list, such as
#' the modules of a multi-trait study measured on the same time points, and
#' collects the results in one table.
#'
#' By default every unordered pair is analyzed once with
#' `bidirectional = TRUE`, so both directions are covered by the sign of the
#' lag: a positive optimal lag means that the first data set of the pair
#' (`x`) leads the second (`y`), and a negative one that `y` leads `x`. The
#' global p-values of the pairs are then adjusted for the number of pairs
#' with `p.adjust.method`, because the max-statistic test accounts for the
#' lags scanned within a pair but not for the pairs examined.
#'
#' Every pair uses the same `seed`, so results do not depend on the order or
#' number of pairs.
#'
#' @param data A named list of numeric matrices, data frames or `dist`
#'   objects, all with the same number of rows (time points).
#' @param pairs Pairs to analyze: `NULL` (default) for all unordered pairs,
#'   or a list of character vectors of length 2, or a two-column character
#'   matrix, giving names of `data`. The first element of each pair is
#'   passed as `X`.
#' @param lags,bidirectional Passed to [lagged_protest()]. `bidirectional`
#'   defaults to `TRUE` here.
#' @param p.adjust.method Correction of the global p-values across pairs,
#'   any method of [stats::p.adjust()]. Defaults to `"BH"`.
#' @param seed Passed to [lagged_protest()] for every pair.
#' @param ... Further arguments passed to [lagged_protest()], such as
#'   `permutation`, `nperm` or `period`.
#'
#' @return A data frame with one row per pair and columns `x`, `y`,
#'   `optimal_lag`, `n` (rows at the optimal lag), `m2`, `r`, `p.value` (raw
#'   p at the optimal lag), `p.global`, the effect measures `Z` and `A` at
#'   the effect lag (see `effect` in [lagged_protest()]), and
#'   `p.global.adjusted` (adjusted across pairs). Unlike `m2`, `Z` and `A`
#'   are on a common scale across pairs with different numbers of rows. The [lagged_protest()] fits are stored, named
#'   `"x -> y"`, in the attribute `"fits"`.
#'
#' @seealso [lagged_protest()]
#'
#' @examples
#' sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
#'                               seed = 1)
#' blocks <- list(A = sim$X, B = sim$Y,
#'                C = simulate_lagged_series(n_times = 12, seed = 2)$X)
#' res <- lagged_protest_pairs(blocks, lags = 0:3, nperm = 199, seed = 1)
#' res
#' plot(attr(res, "fits")[["A -> B"]])
#' @export
lagged_protest_pairs <- function(data, pairs = NULL, lags = 0:6,
                                 bidirectional = TRUE,
                                 p.adjust.method = "BH", seed = NULL, ...) {
  if (!is.list(data) || is.data.frame(data) || length(data) < 2 ||
      is.null(names(data)) || any(names(data) == "") ||
      anyDuplicated(names(data))) {
    stop("`data` must be a list of at least two data sets with unique names.",
         call. = FALSE)
  }
  p.adjust.method <- match.arg(p.adjust.method, stats::p.adjust.methods)
  pairs <- as_pairs(pairs, names(data))

  fits <- lapply(seq_along(pairs), function(i) {
    p <- pairs[[i]]
    run <- function() {
      lagged_protest(data[[p[1]]], data[[p[2]]], lags = lags,
                     bidirectional = bidirectional, seed = seed, ...)
    }
    # A truncation message, if any, is the same for every pair
    if (i == 1) run() else suppressMessages(run())
  })
  names(fits) <- vapply(pairs, paste, character(1), collapse = " -> ")

  out <- data.frame(
    x = vapply(pairs, `[`, character(1), 1),
    y = vapply(pairs, `[`, character(1), 2),
    optimal_lag = vapply(fits, `[[`, integer(1), "optimal_lag"),
    n = vapply(fits, function(f) f$profile$n[f$profile$lag == f$optimal_lag],
               numeric(1)),
    m2 = vapply(fits, `[[`, numeric(1), "m2"),
    r = vapply(fits, `[[`, numeric(1), "r"),
    p.value = vapply(fits, `[[`, numeric(1), "p.value"),
    p.global = vapply(fits, `[[`, numeric(1), "p.global"),
    Z = vapply(fits, function(f) f$effect$Z, numeric(1)),
    A = vapply(fits, function(f) f$effect$A, numeric(1)),
    row.names = NULL
  )
  out$p.global.adjusted <- stats::p.adjust(out$p.global,
                                           method = p.adjust.method)
  attr(out, "fits") <- fits
  out
}

as_pairs <- function(pairs, nms) {
  if (is.null(pairs)) return(utils::combn(nms, 2, simplify = FALSE))
  if (is.matrix(pairs)) {
    if (ncol(pairs) != 2) {
      stop("A matrix of `pairs` must have two columns.", call. = FALSE)
    }
    pairs <- split(pairs, row(pairs))
  }
  if (!is.list(pairs) ||
      !all(vapply(pairs, function(p) is.character(p) && length(p) == 2,
                  logical(1)))) {
    stop("`pairs` must be a list of character vectors of length 2.",
         call. = FALSE)
  }
  bad <- setdiff(unlist(pairs), nms)
  if (length(bad)) {
    stop("Unknown data sets in `pairs`: ", paste(bad, collapse = ", "), ".",
         call. = FALSE)
  }
  if (any(vapply(pairs, function(p) p[1] == p[2], logical(1)))) {
    stop("A pair must contain two different data sets.", call. = FALSE)
  }
  unname(lapply(pairs, as.character))
}
