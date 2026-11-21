# Extended (slow, simulation-based) tests run only when the environment
# variable LAGPROTEST_EXTENDED_TESTS is "true"; see CONTRIBUTING.md.
skip_if_not_extended <- function() {
  testthat::skip_if_not(
    identical(tolower(Sys.getenv("LAGPROTEST_EXTENDED_TESTS")), "true"),
    "extended tests: set LAGPROTEST_EXTENDED_TESTS=true to run"
  )
}
