#' Lake Washington plankton and limnology blocks
#'
#' Builds three monthly blocks (environment, phytoplankton and zooplankton)
#' for Lake Washington, Washington, USA, from `lakeWAplanktonRaw` in the
#' \pkg{MARSS} package. The data are not shipped with \pkg{lagprotest};
#' this function reads them from an installed copy of \pkg{MARSS}, which
#' must be installed separately.
#'
#' Rows are restricted to the years `from` to `to`. Missing values (a few
#' isolated months in each column) are filled by linear interpolation over
#' the month index when `interpolate = TRUE`; a gap at either end of the
#' selected period is filled with the nearest observed value. Plankton
#' counts are then transformed with `log1p()`. The environmental variables
#' are left on their original scale, since [lagged_protest()] standardizes
#' columns by default (`scale = TRUE`).
#'
#' @param from,to First and last year to keep (between 1962 and 1994).
#' @param interpolate Logical. Fill missing values by linear interpolation?
#'   If `FALSE`, missing values are kept and must be handled before calling
#'   [lagged_protest()].
#'
#' @return A named list of three numeric matrices with one row per month
#'   (row names `"YYYY-MM"`):
#'   \describe{
#'     \item{`Environment`}{Water temperature (`Temp`, °C), total
#'       phosphorus (`TP`) and `pH`.}
#'     \item{`Phytoplankton`}{`log1p` of the counts (cells per mL) of
#'       `Cryptomonas`, `Diatoms`, `Greens`, `Unicells` and `Other.algae`.}
#'     \item{`Zooplankton`}{`log1p` of the counts (individuals per L) of
#'       `Cyclops`, `Diaptomus`, `Daphnia`, `Epischura`, `Conochilus`,
#'       `Non.daphnid.cladocerans` and `Non.colonial.rotifers`.}
#'   }
#'
#' @source `lakeWAplanktonRaw` in the \pkg{MARSS} package (Holmes et al.,
#'   GPL-2), described there as adapted from the Lake Washington database
#'   of W. T. Edmondson, funded by the Andrew Mellon Foundation, with data
#'   courtesy of D. Schindler, University of Washington. See
#'   `help("lakeWAplankton", package = "MARSS")` for the original
#'   documentation and how to cite it.
#'
#' @examplesIf requireNamespace("MARSS", quietly = TRUE)
#' lw <- lakewa_blocks()
#' str(lw)
#' fit <- lagged_protest(lw$Environment, lw$Zooplankton, lags = 0:3,
#'                       nperm = 199, seed = 1)
#' fit$optimal_lag
#' @export
lakewa_blocks <- function(from = 1962, to = 1994, interpolate = TRUE) {
  if (!has_marss()) {
    stop("`lakewa_blocks()` needs the MARSS package: ",
         "install.packages(\"MARSS\")", call. = FALSE)
  }
  stopifnot(is.numeric(from), length(from) == 1, is.numeric(to),
            length(to) == 1, from <= to, from >= 1962, to <= 1994,
            is.logical(interpolate), length(interpolate) == 1,
            !is.na(interpolate))

  e <- new.env(parent = emptyenv())
  utils::data("lakeWAplankton", package = "MARSS", envir = e)
  d <- as.data.frame(e$lakeWAplanktonRaw)
  d <- d[d$Year >= from & d$Year <= to, , drop = FALSE]

  env <- c("Temp", "TP", "pH")
  phyto <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
  zoo <- c("Cyclops", "Diaptomus", "Daphnia", "Epischura", "Conochilus",
           "Non.daphnid.cladocerans", "Non.colonial.rotifers")
  vars <- c(env, phyto, zoo)

  if (interpolate) {
    d[vars] <- lapply(d[vars], fill_linear)
  }

  rn <- sprintf("%d-%02d", as.integer(d$Year), as.integer(d$Month))
  block <- function(cols, f = identity) {
    m <- f(as.matrix(d[cols]))
    dimnames(m) <- list(rn, cols)
    m
  }
  list(Environment = block(env),
       Phytoplankton = block(phyto, log1p),
       Zooplankton = block(zoo, log1p))
}

has_marss <- function() requireNamespace("MARSS", quietly = TRUE)

# Linear interpolation over the index; end gaps take the nearest value.
fill_linear <- function(x) {
  if (sum(!is.na(x)) < 2) {
    stop("Too few observed values to interpolate in the selected period.",
         call. = FALSE)
  }
  idx <- seq_along(x)
  stats::approx(idx, x, idx, rule = 2)$y
}
