## Submission

This is a new submission.

lagprotest extends the Procrustean randomization test (PROTEST) to
multivariate time series: it scans Procrustes concordance across temporal
lags and tests it with phase-randomized surrogates, cyclic shifts or free
permutations, with a global p-value for the scan.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Diogo B. Provete <dbprovete@gmail.com>'
  New submission

## Test environments

* local: macOS 27.2, R 4.6.1 (`R CMD check --as-cran`)
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel,
  release, oldrel-1)
* win-builder: R-release (4.6.2) and R-devel (2026-11-17 r90926), 1 NOTE
  (new submission) on both

## Notes for the reviewers

* `elpa_blocks()` downloads a public data archive (SiBBr IPT, CC BY-NC 4.0)
  to `tempdir()` unless the user gives a `path`. Its example needs an
  internet connection, so it runs only in
  interactive sessions (`@examplesIf interactive()`); its tests use a small
  local archive or mock the download and never access the network.
* `lakewa_blocks()` uses a data set from the suggested package MARSS and
  stops with an informative message when MARSS is not installed; its
  example (`@examplesIf`) and tests are skipped in that case.
* All examples run in under 1 s.

## Downstream dependencies

There are currently no downstream dependencies for this package.
