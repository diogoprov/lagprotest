# Summarize a lagged PROTEST

Summarize a lagged PROTEST

## Usage

``` r
# S3 method for class 'lagged_protest'
summary(object, ...)

# S3 method for class 'summary.lagged_protest'
print(x, digits = 3, ...)
```

## Arguments

- object:

  A `"lagged_protest"` object.

- ...:

  Unused.

- x:

  A `"summary.lagged_protest"` object.

- digits:

  Number of significant digits to print.

## Value

An object of class `"summary.lagged_protest"`: the lag profile (a data
frame), with the optimal lag flagged in `optimal`. The scan-level
measures are attached as attributes: `optimal_lag`, `p.global`, `effect`
(lag of the largest z, `Z`, `A` and corrected `RV`) and `ar1` (raw and
residual lag-1 autocorrelation of `X` and `Y`). Printing it shows the
profile followed by these measures; it can otherwise be used as a data
frame.

## Examples

``` r
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
                              seed = 1)
fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 199, seed = 1)
s <- summary(fit)
s
#>  lag  n     m2     r p.value p.adjusted nperm      z      a optimal
#>    0 12 0.1532 0.920   0.970      0.970   199 -1.707 -0.223   FALSE
#>    1 11 0.1000 0.949   0.135      0.135   199  1.130  0.181   FALSE
#>    2 10 0.0496 0.975   0.005      0.005   199  3.406  0.581    TRUE
#>    3  9 0.1303 0.933   0.705      0.705   199 -0.553 -0.121   FALSE
#>    4  8 0.1516 0.921   0.815      0.815   199 -0.944 -0.233   FALSE
#> 
#> Optimal lag 2, global p = 0.03
#> Effect at lag 2 (largest z): Z = 3.02, A = 0.464, corrected RV = 0.193
#> Mean lag-1 autocorrelation, raw (residual): X 0.76 (-0.018), Y 0.66 (-0.049)
attr(s, "effect")
#> $lag
#> [1] 2
#> 
#> $Z
#> [1] 3.020769
#> 
#> $A
#> [1] 0.4636411
#> 
#> $RV
#> [1] 0.1930574
#> 
```
