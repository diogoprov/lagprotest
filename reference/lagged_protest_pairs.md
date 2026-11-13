# Lagged PROTEST for several pairs of data sets

Runs
[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
on pairs of data sets from a named list, such as the modules of a
multi-trait study measured on the same time points, and collects the
results in one table.

## Usage

``` r
lagged_protest_pairs(
  data,
  pairs = NULL,
  lags = 0:6,
  bidirectional = TRUE,
  p.adjust.method = "BH",
  seed = NULL,
  ...
)
```

## Arguments

- data:

  A named list of numeric matrices, data frames or `dist` objects, all
  with the same number of rows (time points).

- pairs:

  Pairs to analyze: `NULL` (default) for all unordered pairs, or a list
  of character vectors of length 2, or a two-column character matrix,
  giving names of `data`. The first element of each pair is passed as
  `X`.

- lags, bidirectional:

  Passed to
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md).
  `bidirectional` defaults to `TRUE` here.

- p.adjust.method:

  Correction of the global p-values across pairs, any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Defaults
  to `"BH"`. Because this argument is taken by `lagged_protest_pairs()`,
  the per-lag adjustment of
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  cannot be set here and stays `"none"`; `p.global` already accounts for
  the lags scanned within each pair.

- seed:

  Passed to
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
  for every pair.

- ...:

  Further arguments passed to
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md),
  such as `permutation`, `nperm`, `period` or `stat`.

## Value

A data frame with one row per pair and columns `x`, `y`, `optimal_lag`,
`n` (rows at the optimal lag), `m2`, `r`, `p.value` (raw p at the
optimal lag), `p.global`, `effect_lag` (the lag of the largest
standardized concordance \\z\\), the effect measures `Z`, `A` and the
corrected `RV` at `effect_lag` (see `effect` in
[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)),
and `p.global.adjusted` (adjusted across pairs). `effect_lag` can differ
from `optimal_lag`, which favors long lags (see Details in
[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)).
Unlike `m2`, `Z`, `A` and `RV` are on a common scale across pairs with
different numbers of rows. The
[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
fits are stored, named `"x -> y"`, in the attribute `"fits"`.

## Details

By default every unordered pair is analyzed once with
`bidirectional = TRUE`, so both directions are covered by the sign of
the lag: a positive optimal lag means that the first data set of the
pair (`x`) leads the second (`y`), and a negative one that `y` leads
`x`. The global p-values of the pairs are then adjusted for the number
of pairs with `p.adjust.method`, because the max-statistic test accounts
for the lags scanned within a pair but not for the pairs examined.

Every pair uses the same `seed`, so results do not depend on the order
or number of pairs.

## See also

[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)

## Examples

``` r
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
                              seed = 1)
blocks <- list(A = sim$X, B = sim$Y,
               C = simulate_lagged_series(n_times = 12, seed = 2)$X)
res <- lagged_protest_pairs(blocks, lags = 0:3, nperm = 199, seed = 1)
res
#>   x y optimal_lag  n         m2         r p.value p.global effect_lag
#> 1 A B           2 10 0.04963006 0.9748692   0.005    0.035          2
#> 2 A C           3  9 0.21583858 0.8855289   0.095    0.415          3
#> 3 B C           1 11 0.20289455 0.8928076   0.010    0.070          1
#>              Z          A        RV p.global.adjusted
#> 1  2.949030204 0.43844113 0.1275756             0.105
#> 2 -0.007160969 0.09816348 0.1766439             0.415
#> 3  2.093882837 0.19365927 0.2053255             0.105
plot(attr(res, "fits")[["A -> B"]])
```
