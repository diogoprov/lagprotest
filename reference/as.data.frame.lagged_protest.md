# Lag profile of a lagged PROTEST as a data frame

Returns the lag profile in long (tidy) form, one row per lag, ready for
tables or plotting with other tools.

## Usage

``` r
# S3 method for class 'lagged_protest'
as.data.frame(x, row.names = NULL, optional = FALSE, alpha = 0.05, ...)
```

## Arguments

- x:

  A `"lagged_protest"` object.

- row.names, optional:

  Unused; for compatibility with the generic.

- alpha:

  Level of the permutation null quantile in `null.quantile`.

- ...:

  Unused.

## Value

A data frame with the columns of `x$profile` (`lag`, `n`, `m2`, `r`,
`p.value`, `p.adjusted`, `nperm`, `z`, `a`), plus `null.quantile`, the
`alpha` quantile of the permuted \\m^2\\ at each lag (unadjusted), and
`optimal`, which flags the optimal lag.

## See also

[`plot.lagged_protest()`](https://provetelab.org/lagprotest/reference/plot.lagged_protest.md),
[`autoplot.lagged_protest()`](https://provetelab.org/lagprotest/reference/autoplot.lagged_protest.md)

## Examples

``` r
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
                              seed = 1)
fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 199, seed = 1)
as.data.frame(fit)
#>   lag  n         m2         r p.value p.adjusted nperm          z          a
#> 1   0 12 0.15315445 0.9202421   0.970      0.970   199 -1.7065939 -0.2225066
#> 2   1 11 0.10003719 0.9486637   0.135      0.135   199  1.1300361  0.1812274
#> 3   2 10 0.04963006 0.9748692   0.005      0.005   199  3.4064381  0.5811748
#> 4   3  9 0.13033905 0.9325561   0.705      0.705   199 -0.5527663 -0.1208978
#> 5   4  8 0.15160634 0.9210829   0.815      0.815   199 -0.9442243 -0.2333077
#>   null.quantile optimal
#> 1    0.09891132   FALSE
#> 2    0.08864801   FALSE
#> 3    0.08167769    TRUE
#> 4    0.07840183   FALSE
#> 5    0.07242045   FALSE
```
