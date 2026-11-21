#' @srrstats {G5.10} Simulation-based tests of Type I error and power run
#'   only when LAGPROTEST_EXTENDED_TESTS is "true"; GitHub Actions sets it
#'   in the check and coverage workflows.
#' @srrstats {G5.12} CONTRIBUTING.md describes how to run the extended tests
#'   and their runtime.
#' @noRd
NULL

# Extended (slow, simulation-based) tests run only when the environment
# variable LAGPROTEST_EXTENDED_TESTS is "true"; see CONTRIBUTING.md.
skip_if_not_extended <- function() {
  testthat::skip_if_not(
    identical(tolower(Sys.getenv("LAGPROTEST_EXTENDED_TESTS")), "true"),
    "extended tests: set LAGPROTEST_EXTENDED_TESTS=true to run"
  )
}

# autotest traces package functions with typetracer, which draws random
# numbers itself; tests that the global random number stream is untouched
# are skipped while tracing.
skip_if_traced <- function() {
  testthat::skip_if(
    nzchar(Sys.getenv("TYPETRACER_LEAVE_TRACES")),
    "typetracer draws random numbers while tracing"
  )
}
