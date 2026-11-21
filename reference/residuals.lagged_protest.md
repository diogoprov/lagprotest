# Per-observation Procrustes residuals at the optimal lag

Per-observation Procrustes residuals at the optimal lag

## Usage

``` r
# S3 method for class 'lagged_protest'
residuals(object, ...)
```

## Arguments

- object:

  A `"lagged_protest"` object.

- ...:

  Unused.

## Value

A data frame with the paired row labels of `X` and `Y` and the residual
of each pair after superimposition. Residuals are on the scale of the
configurations that were compared: after standardization
(`scale = TRUE`) and after any `deseason` or `prewhiten` step, not on
the scale of the original variables. They are not back-transformed.

## Examples

``` r
sim <- simulate_lagged_series(n_times = 36, true_lag = 2, noise_sd = 0.3,
                              seed = 1)
fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 199, seed = 1)
res <- residuals(fit)
head(res)
#>   x_row y_row   residual
#> 1     1     3 0.06270259
#> 2     2     4 0.04671930
#> 3     3     5 0.03868620
#> 4     4     6 0.05368182
#> 5     5     7 0.05948441
#> 6     6     8 0.08208955
# Pairs that fit worst at the optimal lag
res[order(res$residual, decreasing = TRUE)[1:3], ]
#>    x_row y_row   residual
#> 6      6     8 0.08208955
#> 13    13    15 0.07193485
#> 28    28    30 0.06775567
```
