# Summarize a lagged PROTEST

Summarize a lagged PROTEST

## Usage

``` r
# S3 method for class 'lagged_protest'
summary(object, ...)
```

## Arguments

- object:

  A `"lagged_protest"` object.

- ...:

  Unused.

## Value

The lag profile (a data frame), with the optimal lag flagged. The
scan-level measures are attached as attributes: `p.global`, `effect`
(lag of the largest z, `Z`, `A` and corrected `RV`) and `ar1` (raw and
residual lag-1 autocorrelation of `X` and `Y`).
