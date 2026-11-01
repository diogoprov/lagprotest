# Plot the lag profile of a lagged PROTEST

Plots \\m^2\\ against lag. The dashed line is the `alpha` quantile of
the permutation null distribution at each lag (unadjusted); observed
values below it are significant at that level before correction. Filled
points mark lags whose p-value (adjusted, if `p.adjust.method` was not
`"none"`) is below `alpha`.

## Usage

``` r
# S3 method for class 'lagged_protest'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  A `"lagged_protest"` object.

- alpha:

  Significance level.

- ...:

  Further arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.
