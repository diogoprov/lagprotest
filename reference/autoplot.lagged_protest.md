# Plot the lag profile of a lagged PROTEST with ggplot2

A ggplot2 version of
[`plot.lagged_protest()`](https://provetelab.org/lagprotest/reference/plot.lagged_protest.md):
\\m^2\\ against lag, with the `alpha` quantile of the permutation null
distribution at each lag as a dashed line, filled points for lags whose
p-value (adjusted, if `p.adjust.method` was not `"none"`) is below
`alpha`, and the optimal lag circled, filled when the global p-value is
below `alpha`. The global p-value is shown in the subtitle. Requires the
ggplot2 package.

## Usage

``` r
# S3 method for class 'lagged_protest'
autoplot(object, alpha = 0.05, ...)
```

## Arguments

- object:

  A `"lagged_protest"` object.

- alpha:

  Significance level.

- ...:

  Unused.

## Value

A `ggplot` object, which can be modified with further layers.

## See also

[`plot.lagged_protest()`](https://provetelab.org/lagprotest/reference/plot.lagged_protest.md),
[`as.data.frame.lagged_protest()`](https://provetelab.org/lagprotest/reference/as.data.frame.lagged_protest.md)

## Examples

``` r
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
                              seed = 1)
fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 199, seed = 1)
library(ggplot2)
autoplot(fit)

autoplot(fit) + theme_bw()
```
