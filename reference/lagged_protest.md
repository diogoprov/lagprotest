# Lagged Procrustean randomization test

Scans Procrustes concordance between two multivariate series across
candidate temporal lags and tests it by permutation at each lag.

## Usage

``` r
lagged_protest(
  X,
  Y,
  lags = 0:6,
  permutation = c("cyclic", "free"),
  nperm = 999,
  mirror = TRUE,
  bidirectional = FALSE,
  scale = TRUE,
  p.adjust.method = "holm",
  seed = NULL
)
```

## Arguments

- X, Y:

  Numeric matrices or data frames with rows as ordered time points and
  columns as variables, or `dist` objects, which are converted to
  principal coordinates. Both must have the same number of rows.

- lags:

  Integer vector of non-negative candidate lags.

- permutation:

  Permutation scheme, `"cyclic"` or `"free"`.

- nperm:

  Number of random permutations per lag when `permutation = "free"`.

- mirror:

  Logical; with `permutation = "cyclic"`, also use the time-reversed
  shifts of `Y`, doubling the permutation set.

- bidirectional:

  Logical; if `TRUE`, also test the negative lags `-lags` (Y leads X).

- scale:

  Logical; z-score the columns of `X` and `Y` (computed on the full
  series) before analysis. Ignored for `dist` inputs.

- p.adjust.method:

  Correction for multiple lags, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

- seed:

  Optional integer seed for `permutation = "free"`. The global random
  number stream is left untouched.

## Value

An object of class `"lagged_protest"`, a list with:

- `optimal_lag`:

  Lag with the smallest \\m^2\\.

- `m2`:

  \\m^2\\ at the optimal lag.

- `r`:

  Procrustes correlation, \\\sqrt{1 - m^2}\\, at the optimal lag.

- `p.value`:

  Adjusted p-value at the optimal lag.

- `profile`:

  Data frame with `lag`, `n` (rows used), `m2`, `r`, `p.value`,
  `p.adjusted` and `nperm` for every lag.

- `residuals`:

  Data frame of per-observation Procrustes residuals at the optimal lag,
  with the paired row labels of `X` and `Y`.

- `procrustes`:

  The
  [`vegan::procrustes()`](https://vegandevs.github.io/vegan/reference/procrustes.html)
  fit at the optimal lag.

- `null`:

  List of permuted \\m^2\\ values, one element per lag.

- `permutation`, `p.adjust.method`, `n`, `call`:

  Settings used.

## Details

Rows of `X` and `Y` are ordered time points (e.g., monthly centroids)
and must correspond one to one. At lag `k > 0`, `X` at time `t` is
paired with `Y` at time `t + k` (X leads Y), and the series are
truncated to the overlapping window of `n - k` rows. Negative lags (Y
leads X) are tested when `bidirectional = TRUE`. The statistic at each
lag is the symmetric Procrustes sum of squares \\m^2\\ from
[`vegan::procrustes()`](https://vegandevs.github.io/vegan/reference/procrustes.html);
the lag of maximal concordance is the one with the smallest \\m^2\\.

Two permutation schemes are available:

- `"cyclic"` (default): the full `Y` series is cyclically shifted by
  `1, ..., n - 1` positions before lagging, and optionally also
  time-reversed (`mirror = TRUE`). This keeps the temporal
  autocorrelation of `Y` under the null. All `n - 1` (or `2(n - 1)`)
  shifts are enumerated, so `nperm` is ignored, and the smallest
  attainable p-value is `1 / (2(n - 1) + 1)` (0.043 for `n = 12` with
  `mirror = TRUE`). The scheme assumes the series spans complete cycles,
  i.e., that its period matches the number of rows (e.g., 12 monthly
  samples of an annual cycle).

- `"free"`: rows of `Y` are permuted at random within the overlapping
  window, as in
  [`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html).
  This assumes exchangeable rows and can strongly inflate Type I error
  when both series are autocorrelated.

P-values are computed as `(b + 1) / (N + 1)`, where `b` is the number of
permuted \\m^2\\ values smaller than or equal to the observed one and
`N` the number of permutations, and adjusted across lags with
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Note that
\\m^2\\ values at different lags are based on different numbers of rows
(`n - |k|`).

## See also

[`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html),
[`vegan::procrustes()`](https://vegandevs.github.io/vegan/reference/procrustes.html)

## Examples

``` r
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
                              seed = 1)
fit <- lagged_protest(sim$X, sim$Y, lags = 0:4)
fit
#> 
#> Lagged PROTEST
#> 
#> Call: lagged_protest(X = sim$X, Y = sim$Y, lags = 0:4)
#> Permutation: cyclic shift with mirror; n = 12 time points
#> P-value adjustment across lags: holm
#> 
#> Optimal lag: 2  (m2 = 0.0496, r = 0.975, adjusted p = 0.217)
#> 
#> Lag profile:
#>  lag  n     m2     r p.value p.adjusted nperm
#>    0 12 0.1532 0.920  1.0000      1.000    22
#>    1 11 0.1000 0.949  0.0870      0.348    22
#>    2 10 0.0496 0.975  0.0435      0.217    22
#>    3  9 0.1303 0.933  0.7391      1.000    22
#>    4  8 0.1516 0.921  0.9130      1.000    22
#> 
#> Smallest attainable raw p-value: 0.0435
plot(fit)

residuals(fit)
#>    x_row y_row   residual
#> 1      1     3 0.08219361
#> 2      2     4 0.08428298
#> 3      3     5 0.03545903
#> 4      4     6 0.07590098
#> 5      5     7 0.05382336
#> 6      6     8 0.05169356
#> 7      7     9 0.07466031
#> 8      8    10 0.06927409
#> 9      9    11 0.07813342
#> 10    10    12 0.08188577
```
