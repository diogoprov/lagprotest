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

  Possibly misspelled words in DESCRIPTION:
    Neto (10:64)
    Peres (10:58)

  These are not misspellings: they are the surname of Pedro Peres-Neto,
  co-author of a cited reference (Peres-Neto and Jackson 2001).

## Test environments

* local: macOS 27.2, R 4.6.1 (`R CMD check --as-cran`, with remote checks
  and the PDF manual)
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel,
  release, oldrel-1)
* win-builder: R-devel (2026-11-17 r90926) on 0.9.1, 1 NOTE (above);
  R-release (4.6.2) and R-devel on 0.8.0, 1 NOTE (new submission).
* R-hub (on 0.9.1): Linux, Windows, macOS and macOS arm64 (R-devel), plus
  clang ASAN and valgrind: Status OK on all, no valgrind errors

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
