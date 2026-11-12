# Changelog

## lagprotest 0.7.0

- **Breaking:**
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  gains `stat = c("minP", "z", "raw")`, the scan statistic used for the
  global p-value and the optimal lag, with `"minP"` as the new default.
  `"minP"` compares lags by their per-lag p-values (Westfall and Young’s
  minP), `"z"` by *m*² standardized by its own null distribution, and
  `"raw"` by *m*² itself, the behavior of earlier versions. All three
  use the same permutations and are valid global tests; `"raw"` favors
  the longest lags, which rest on fewer rows, and has less power for
  short lags. Use `stat = "raw"` to reproduce results from version 0.6
  and earlier. The liver example article keeps `"raw"` to match the
  analysis submitted for publication.

- New
  [`lakewa_blocks()`](https://provetelab.org/lagprotest/reference/lakewa_blocks.md)
  builds monthly environment, phytoplankton and zooplankton blocks for
  Lake Washington (1962–1994) from `lakeWAplanktonRaw` in the MARSS
  package (now in Suggests). The data are read from MARSS rather than
  shipped with lagprotest.

- [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  now returns continuous measures of evidence and strength computed from
  the same permutations as the p-values: per-lag `z` and `a` in the
  profile, and in `effect` the lag of the largest `z`, its standardized
  position `Z` among the permuted maxima, and the chance-corrected `A`
  and RV coefficient. The permuted RV values are kept in `null.rv`.
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) show them, and
  [`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md)
  adds `Z` and `A` as columns.

- [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  now reports the mean lag-1 autocorrelation of each series (`ar1`), raw
  and after removing the seasonal harmonics, as a description of how
  close the data are to the conditions under which the phase scheme was
  liberal in simulations.

- [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  now issues a message when the smallest attainable p-value is 0.05 or
  more (e.g., the cyclic scheme without mirror at `n = 12`).

- The documentation of
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  and the permutation-schemes article cite Manly and Navarro
  Alberto (2020) for the p-value formula and the recommended number of
  permutations, and pp. 83–84 for the max-statistic adjustment of
  correlated tests.

- The liver example article now checks the lag scan with two statistics
  that weigh all lags equally (minP and the standardized m²), computed
  from the surrogates stored in each fit. The raw minimum m² favors long
  lags, and the climate–histochemistry delay is weaker under the
  balanced statistics.

- The permutation loop is about five times faster: *m*² and the RV
  coefficient of each permuted pair are now computed directly from the
  singular values of the cross-product matrix instead of calling
  [`vegan::procrustes()`](https://vegandevs.github.io/vegan/reference/procrustes.html)
  for every permutation. Results with a fixed `seed` are unchanged
  (p-values and lags identical, *m*² equal to floating-point rounding);
  a new regression test guards them.

- [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  and
  [`simulate_lagged_series()`](https://provetelab.org/lagprotest/reference/simulate_lagged_series.md)
  now check their arguments and name the offending one in the error:
  `nperm` and `seed` must be single whole numbers, `lags` must not be
  empty, `mirror`, `bidirectional` and `scale` must be `TRUE` or
  `FALSE`, and the counts, `period`, `amplitude` and `noise_sd` of the
  simulator must be valid. Previously, `nperm = 2.5` was silently
  truncated and an invalid `seed` produced an unrelated warning.

- The minimum-P message no longer suggests `mirror = TRUE` when the
  cyclic scheme already uses it.

- The fits stored by
  [`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md)
  now carry a call that names the data sets and shows the values of
  `lags`, `bidirectional` and `seed`, so
  [`print()`](https://rdrr.io/r/base/print.html) is informative and
  `eval(fit$call)` reproduces the fit.

- [`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md)
  gains the columns `effect_lag` (the lag to which `Z` and `A` refer,
  which can differ from `optimal_lag`) and `RV`.

- [`summary()`](https://rdrr.io/r/base/summary.html) now returns an
  object of class `"summary.lagged_protest"` (still a data frame) whose
  [`print()`](https://rdrr.io/r/base/print.html) method shows the global
  p-value, the effect measures and the autocorrelation below the
  profile; before, these were attributes that were not printed.

- [`print()`](https://rdrr.io/r/base/print.html) notes when the optimal
  lag (smallest *m*²) and the effect lag (largest *z*) differ, and no
  longer runs the lines of a long call together.

- The documentation of
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md),
  the “Get started” vignette and the README now explain that the
  smallest *m*² favors the longest lags scanned, describe the effect
  measures and `ar1`, and point to
  [`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md),
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  and the liver example. Smaller fixes: the `mirror` element is
  documented, `n_perm` is described as the same at every lag, and the
  per-lag adjustment that
  [`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md)
  cannot set is noted.

- New `inst/CITATION`, so `citation("lagprotest")` gives the Zenodo DOI.

## lagprotest 0.6.0

- New [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  method returns the lag profile in long form, with the permutation null
  quantile at each lag and a flag for the optimal lag.
- New
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  method draws the lag profile with ggplot2 (suggested package), marking
  the optimal lag and its global p-value as
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) does.

## lagprotest 0.5.0

- New data set `dminutus_liver`: monthly means of four hepatic modules
  and of climate for the treefrog *Dendropsophus minutus* over one
  annual cycle (Franco-Belussi et al. 2026, CC BY 4.0), with a new
  article analyzing it with
  [`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md).

## lagprotest 0.4.0

- New
  [`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md)
  runs
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  on pairs of data sets from a named list (by default all unordered
  pairs, both directions through the sign of the lag) and returns one
  table, with the global p-values adjusted across pairs (`"BH"` by
  default) and the fits kept in an attribute.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) now circles
  the optimal lag, filled when the global p-value is below `alpha`, and
  prints the global p-value above the plot (`global = FALSE` restores
  the previous plot).
- Simulations of the `period` argument with the full design (1000
  replicates per scenario, 12–96 time points) are summarized in the
  article “Choosing a permutation scheme”: after truncation the global
  p-value was calibrated up to AR(1) of 0.5 and liberal with strong
  autocorrelation, as for series spanning complete cycles.

## lagprotest 0.3.0

- New `period` argument in
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md):
  when the seasonal cycle does not divide the series length, both series
  are truncated to their first `floor(n / period) * period` rows, so
  that the phase and cyclic schemes, which treat the series as periodic
  with period `n`, remain calibrated. In simulations (24–48 time points,
  cycle of 0.75 `n`), this brought the Type I error of the global
  p-value from 0.46–0.94 back to 0.04–0.05 for AR(1) noise of 0.3.
  Removing the seasonal component by harmonic regression was also tested
  and left the test liberal, so it was not added.

## lagprotest 0.2.0

### Breaking changes

- The default permutation scheme is now `permutation = "phase"`
  (phase-randomized surrogates), which retains power for the global
  p-value across lags. Use `permutation = "cyclic"` for the previous
  behavior; the phase scheme is random, so set `seed` for reproducible
  results.

- [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  now defaults to `p.adjust.method = "none"`, because corrections such
  as Holm are often unattainable with the few permutations of the cyclic
  scheme. `$p.value` is now the raw p-value at the optimal lag; the
  adjusted value is in the new `$p.adjusted` element, and the number of
  permutations in `$n_perm`.
  [`print()`](https://rdrr.io/r/base/print.html) warns when p-values are
  unadjusted.

### New features

- [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  returns a global p-value, `$p.global`, for the scan across lags. It is
  a max-statistic test that compares the observed minimum *m*² across
  lags with the minimum across lags under each permutation, and so
  accounts for testing several lags and selecting the optimal one.
  [`print()`](https://rdrr.io/r/base/print.html) reports it below the
  optimal lag. Simulations showed that the raw p-value at the optimal
  lag rejects 10–28% of the time (21–25% on average) when the series are
  independent.
- New scheme `permutation = "phase"`: phase-randomized surrogates of `Y`
  that keep the amplitude spectrum of each column and the cross-spectra
  among columns. Unlike cyclic shifts, surrogates do not reproduce the
  observed alignment at other lags, so the global p-value keeps its
  power.
- `p.adjust.method = "maxT"` gives single-step max-statistic adjusted
  p-values per lag.
- `permutation = "free"` now applies each permutation jointly to all
  lags, as the max-statistic requires.

### Infrastructure

- Continuous integration: R CMD check on macOS, Windows and Ubuntu, and
  test coverage reporting, via GitHub Actions.
- Spell check of the documentation, contributing guide, code of conduct
  and issue templates.
- The vignette now uses US English spelling.

## lagprotest 0.1.0

First public release.

- [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  scans Procrustes concordance between two ordered multivariate series
  across candidate lags and tests each lag by permutation. Two schemes
  are available: cyclic shifts of the full series, optionally with time
  reversal (default), and free row permutation as in
  [`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html).
  P-values are adjusted across lags. Accepts matrices, data frames or
  `dist` objects (converted by PCoA).
- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) methods for
  the lag profile, the optimal lag and per-observation Procrustes
  residuals.
- [`simulate_lagged_series()`](https://provetelab.org/lagprotest/reference/simulate_lagged_series.md)
  generates seasonal series with AR(1) noise and a known lag, for
  calibration and power studies.
