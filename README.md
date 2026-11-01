# lagprotest <img src="man/figures/logo.png" align="right" height="139" alt="lagprotest hex logo" />

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23860109.svg)](https://doi.org/10.5281/zenodo.23860109)
<!-- badges: end -->

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
- `permutation = "free"`: free permutation of rows within the overlapping
  window, as in `vegan::protest()`. Rows are assumed exchangeable, which
  autocorrelated series violate; Type I error can be severely inflated.

## Installation

```r
# install.packages("remotes")
remotes::install_github("diogoprov/lagprotest")
```

## Example

```r
library(lagprotest)

# Two 12-month series in which X leads Y by 2 months
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.5,
                              seed = 1)

fit <- lagged_protest(sim$X, sim$Y, lags = 0:4)
fit             # lag profile, optimal lag, adjusted p-value
plot(fit)       # m2 against lag, with the permutation null quantile
residuals(fit)  # per-observation residuals at the optimal lag
```

Inputs can be matrices or data frames (rows = ordered time points,
columns = variables) or `dist` objects, which are converted to principal
coordinates.

## Status

Under development. The method and its calibration are described in a
manuscript in preparation.

## License

MIT
