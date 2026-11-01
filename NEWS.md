# lagprotest (development version)

## Breaking changes

* `lagged_protest()` now defaults to `p.adjust.method = "none"`, because
  corrections such as Holm are often unattainable with the few permutations
  of the cyclic scheme. `$p.value` is now the raw p-value at the optimal lag;
  the adjusted value is in the new `$p.adjusted` element, and the number of
  permutations in `$n_perm`. `print()` warns when p-values are unadjusted.

## Infrastructure

* Continuous integration: R CMD check on macOS, Windows and Ubuntu, and test
  coverage reporting, via GitHub Actions.
* Spell check of the documentation, contributing guide, code of conduct and
  issue templates.
* The vignette now uses US English spelling.

# lagprotest 0.1.0

First public release.

* `lagged_protest()` scans Procrustes concordance between two ordered
  multivariate series across candidate lags and tests each lag by
  permutation. Two schemes are available: cyclic shifts of the full
  series, optionally with time reversal (default), and free row
  permutation as in `vegan::protest()`. P-values are adjusted across lags.
  Accepts matrices, data frames or `dist` objects (converted by PCoA).
* `print()`, `summary()`, `plot()` and `residuals()` methods for the
  lag profile, the optimal lag and per-observation Procrustes residuals.
* `simulate_lagged_series()` generates seasonal series with AR(1) noise
  and a known lag, for calibration and power studies.
