# Contributing to lagprotest

Thank you for your interest in lagprotest. Bug reports, questions about the
method and suggestions are all welcome.

## Reporting a bug

Open an [issue](https://github.com/diogoprov/lagprotest/issues/new/choose)
with a minimal reproducible example, ideally made with the
[reprex](https://reprex.tidyverse.org) package. Simulated data from
`simulate_lagged_series()` is usually enough to reproduce a problem without
sharing your own data. Please include the output of `sessionInfo()`.

## Methodological questions

Questions about the permutation schemes, the choice of lags or the
interpretation of the lag profile are welcome as issues too. Please describe
your data (number of time points, sampling interval, number of variables in
each matrix) and what you want to test.

## Pull requests

* Open an issue first to discuss substantial changes.
* Work on a branch, and add or update tests in `tests/testthat/` for any
  change in behavior.
* Document functions with roxygen2 and run `devtools::document()`.
* Make sure `devtools::check()` passes without errors, warnings or notes.
* Add a bullet to the top of `NEWS.md` describing the change.
* Follow the existing code style: base R, native pipe `|>` where a pipe is
  needed, and few dependencies.

## Extended tests

The simulation-based checks of Type I error and power in
`tests/testthat/test-type1_error.R` and `tests/testthat/test-power.R` are
extended tests. They run only when the environment variable
`LAGPROTEST_EXTENDED_TESTS` is `"true"`:

```r
Sys.setenv(LAGPROTEST_EXTENDED_TESTS = "true")
devtools::test()
```

They need no data downloads or special hardware, take under a minute on a
laptop, and produce no files to inspect; failures appear as ordinary test
failures. The GitHub Actions workflows (R-CMD-check and test-coverage) set
the variable, so the extended tests run on every push and pull request but
not on CRAN. Larger calibration grids (thousands of replicates per cell)
are not part of the test suite.

## rOpenSci standards (`srr`)

Conformity with the rOpenSci standards for statistical software (general and
time-series categories) is documented with `@srrstats` and `@srrstatsNA`
roxygen tags next to the code or tests that address each standard; the
package-wide ones and the non-applicable ones, with their reasons, are in
`R/srr-stats-standards.R`. `devtools::document()` uses the `srr` roclet, so
install `srr` from the rOpenSci r-universe
(`install.packages("srr", repos = "https://ropensci-review-tools.r-universe.dev")`).
When a change affects a standard, update its tag, and check with
`srr::srr_stats_pre_submit()`.

## Code of conduct

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating you
agree to abide by its terms.
