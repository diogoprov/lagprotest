# lagprotest 0.5.0

* New data set `dminutus_liver`: monthly means of four hepatic modules and
  of climate for the treefrog *Dendropsophus minutus* over one annual cycle
  (Franco-Belussi et al. 2026, CC BY 4.0), with a new article analyzing it
  with `lagged_protest_pairs()`.

# lagprotest 0.4.0

* New `lagged_protest_pairs()` runs `lagged_protest()` on pairs of data sets
  from a named list (by default all unordered pairs, both directions through
  the sign of the lag) and returns one table, with the global p-values
  adjusted across pairs (`"BH"` by default) and the fits kept in an
  attribute.
* `plot()` now circles the optimal lag, filled when the global p-value is
  below `alpha`, and prints the global p-value above the plot
  (`global = FALSE` restores the previous plot).
* Simulations of the `period` argument with the full design (1000
  replicates per scenario, 12–96 time points) are summarized in the
  article "Choosing a permutation scheme": after truncation the global
  p-value was calibrated up to AR(1) of 0.5 and liberal with strong
  autocorrelation, as for series spanning complete cycles.

# lagprotest 0.3.0

* New `period` argument in `lagged_protest()`: when the seasonal cycle does
  not divide the series length, both series are truncated to their first
  `floor(n / period) * period` rows, so that the phase and cyclic schemes,
  which treat the series as periodic with period `n`, remain calibrated.
  In simulations (24–48 time points, cycle of 0.75 `n`), this brought the
  Type I error of the global p-value from 0.46–0.94 back to 0.04–0.05 for
  AR(1) noise of 0.3. Removing the seasonal component by harmonic
  regression was also tested and left the test liberal, so it was not
  added.

# lagprotest 0.2.0

## Breaking changes

* The default permutation scheme is now `permutation = "phase"`
  (phase-randomized surrogates), which retains power for the global p-value
  across lags. Use `permutation = "cyclic"` for the previous behavior; the
  phase scheme is random, so set `seed` for reproducible results.

* `lagged_protest()` now defaults to `p.adjust.method = "none"`, because
  corrections such as Holm are often unattainable with the few permutations
  of the cyclic scheme. `$p.value` is now the raw p-value at the optimal lag;
  the adjusted value is in the new `$p.adjusted` element, and the number of
  permutations in `$n_perm`. `print()` warns when p-values are unadjusted.

## New features

* `lagged_protest()` returns a global p-value, `$p.global`, for the scan
  across lags. It is a max-statistic test that compares the observed minimum
  *m*² across lags with the minimum across lags under each permutation, and
  so accounts for testing several lags and selecting the optimal one.
  `print()` reports it below the optimal lag. Simulations showed that the
  raw p-value at the optimal lag rejects 10–28% of the time (21–25% on
  average) when the series are independent.
* New scheme `permutation = "phase"`: phase-randomized surrogates of `Y`
  that keep the amplitude spectrum of each column and the cross-spectra
  among columns. Unlike cyclic shifts, surrogates do not reproduce the
  observed alignment at other lags, so the global p-value keeps its power.
* `p.adjust.method = "maxT"` gives single-step max-statistic adjusted
  p-values per lag.
* `permutation = "free"` now applies each permutation jointly to all lags,
  as the max-statistic requires.

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
