#' srr_stats
#'
#' Package-wide standards of the rOpenSci Statistical Software Peer Review
#' (general and time-series categories). Standards addressed in specific
#' functions or tests are tagged there.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.2} The README has a lifecycle badge (stable) and a
#'   "Lifecycle" section: the API is stable since 0.9.0, and changes will go
#'   through deprecation.
#' @srrstats {G1.4} All exported functions are documented with roxygen2.
#' @srrstats {G1.4a} Internal functions are documented in roxygen2 format
#'   with `@noRd`.
#' @srrstats {G1.5} The calibration and power claims in the documentation
#'   (Type I error of each scheme, recovery of a known lag) are reproduced by
#'   the extended tests (test-type1_error.R, test-power.R) and the article
#'   "Permutation schemes", using `simulate_lagged_series()`.
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G2.4c, G2.4d, G2.4e, G2.5} No input is character or factor
#'   data to be converted; character arguments are options matched with
#'   `match.arg()`, and factor columns are errors.
#' @srrstatsNA {G2.14b, G2.14c, TS2.1b, TS2.1c} Missing values are errors
#'   only. Permutation, cyclic shifts and phase randomization need complete,
#'   regularly spaced series, so ignoring missing values would change the
#'   null distribution, and imputation is a modelling choice left to the
#'   user (the Lake Washington example interpolates with `interpolate =
#'   TRUE`).
#' @srrstatsNA {G3.1, G3.1a} No covariance matrices are computed; the
#'   statistic comes from the singular values of the cross-product of the
#'   two configurations.
#' @srrstatsNA {G4.0} No output is written to local files; `elpa_blocks()`
#'   only downloads its data archive to `tempdir()` or a user path.
#' @srrstatsNA {G5.4c} The method is new, and its correctness is tested
#'   against vegan for the simultaneous case and against exact cases with
#'   known answers rather than against stored published values.
#' @srrstatsNA {G5.11, G5.11a} Extended tests use simulated data and need no
#'   downloads.
#' @srrstatsNA {TS2.5, TS2.6} No autocovariance or autocorrelation matrices
#'   are computed or returned.
#' @srrstatsNA {TS3.0, TS3.1, TS3.2, TS3.3, TS3.3a, TS3.3b} The package does
#'   not forecast.
#' @srrstatsNA {TS4.0a} Results have their own class rather than the class
#'   of the input (TS4.0b).
#' @srrstatsNA {TS4.1} The test statistic is dimensionless, so input units
#'   have no counterpart in the results.
#' @srrstatsNA {TS4.5a, TS4.5b} Residuals describe the fit of the compared
#'   configurations and are reported on that scale (TS4.5c); there are no
#'   forecasts to back-transform.
#' @srrstatsNA {TS4.6, TS4.6a, TS4.6b, TS4.6c, TS4.7, TS4.7a, TS4.7b, TS4.7c}
#'   The package does not forecast.
#' @srrstatsNA {TS5.4} No frequency-domain plots are produced.
#' @srrstatsNA {TS5.5} Data with missing values are not accepted, so plots
#'   never contain them.
#' @srrstatsNA {TS5.6, TS5.7, TS5.8} The package does not forecast.
#' @noRd
NULL
