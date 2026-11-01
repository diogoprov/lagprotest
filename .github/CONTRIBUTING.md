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

## Code of conduct

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating you
agree to abide by its terms.
