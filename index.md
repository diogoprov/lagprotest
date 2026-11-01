# lagprotest

**lagprotest** extends the Procrustean randomization test (PROTEST) to
ordered multivariate series. It scans Procrustes concordance between two
configurations across candidate temporal lags, tests each lag by
permutation, and returns the lag profile, the lag of maximal concordance
(minimum *m*²) and per-observation Procrustes residuals.

Two permutation schemes are available:

- `permutation = "cyclic"` (default): cyclic shifts of the full series,
  optionally with its time reversal, applied before lagging. This keeps
  temporal autocorrelation under the null and assumes the series spans
  complete cycles (e.g., 12 monthly samples of an annual cycle).
- `permutation = "free"`: free permutation of rows within the
  overlapping window, as in
  [`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html).
  Rows are assumed exchangeable, which autocorrelated series violate;
  Type I error can be severely inflated.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("diogoprov/lagprotest")
```

## Example

``` r

library(lagprotest)

# Two 12-month series in which X leads Y by 2 months
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.5,
                              seed = 1)

fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, seed = 1)
fit             # lag profile, optimal lag, per-lag and global p-values
plot(fit)       # m2 against lag, with the permutation null quantile
residuals(fit)  # per-observation residuals at the optimal lag
```

Report the global p-value (`fit$p.global`) for the optimal lag: it
accounts for scanning several lags and selecting the best one, whereas
the per-lag p-values do not.

Inputs can be matrices or data frames (rows = ordered time points,
columns = variables) or `dist` objects, which are converted to principal
coordinates.

## Status

Under development. The method and its calibration are described in a
manuscript in preparation.

## License

MIT

## Code of Conduct

Please note that the lagprotest project is released with a [Contributor
Code of
Conduct](https://provetelab.org/lagprotest/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.
