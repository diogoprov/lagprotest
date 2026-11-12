# Lagged Procrustean randomization test

Scans Procrustes concordance between two multivariate series across
candidate temporal lags and tests it by permutation at each lag.

## Usage

``` r
lagged_protest(
  X,
  Y,
  lags = 0:6,
  permutation = c("phase", "cyclic", "free"),
  nperm = 999,
  mirror = TRUE,
  bidirectional = FALSE,
  scale = TRUE,
  p.adjust.method = "none",
  period = NULL,
  seed = NULL
)
```

## Arguments

- X, Y:

  Numeric matrices or data frames with rows as ordered time points and
  columns as variables, or `dist` objects, which are converted to
  principal coordinates. Both must have the same number of rows.

- lags:

  Integer vector of non-negative candidate lags (default `0:6`). Every
  lag must leave at least 3 overlapping rows; long lags rest on few rows
  (6 for lag 6 with `n = 12`) and are favored by the smallest \\m^2\\
  (see Details).

- permutation:

  Permutation scheme, `"phase"`, `"cyclic"` or `"free"`.

- nperm:

  Number of surrogates (`permutation = "phase"`) or random permutations
  (`permutation = "free"`). Ignored for the cyclic scheme, whose
  permutations are enumerated. Use 999 or more for a test at the 5%
  level and about 4999 when the 1% level matters (Manly and Navarro
  Alberto 2020, pp. 57–59).

- mirror:

  Logical; with `permutation = "cyclic"`, also use the time-reversed
  shifts of `Y`, doubling the permutation set.

- bidirectional:

  Logical; if `TRUE`, also test the negative lags `-lags` (Y leads X).

- scale:

  Logical; z-score the columns of `X` and `Y` (computed on the full
  series) before analysis. Ignored for `dist` inputs.

- p.adjust.method:

  Correction of the per-lag p-values for multiple lags: `"maxT"` for the
  single-step max-statistic adjustment, or any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Defaults
  to `"none"`.

- period:

  Optional length of the seasonal cycle, in rows. If given, `X` and `Y`
  are truncated to their first `floor(n / period) * period` rows (whole
  cycles) before analysis; a message reports the rows dropped. `NULL`
  (default) uses all rows.

- seed:

  Optional integer seed for the `"phase"` and `"free"` schemes. The
  global random number stream is left untouched.

## Value

An object of class `"lagged_protest"`, a list with:

- `optimal_lag`:

  Lag with the smallest \\m^2\\.

- `m2`:

  \\m^2\\ at the optimal lag.

- `r`:

  Procrustes correlation, \\\sqrt{1 - m^2}\\, at the optimal lag.

- `p.value`:

  Raw permutation p-value at the optimal lag.

- `p.adjusted`:

  P-value at the optimal lag adjusted with `p.adjust.method` (equal to
  `p.value` when `"none"`).

- `p.global`:

  Max-statistic p-value of the scan across lags, accounting for testing
  several lags and selecting the optimal one.

- `n_perm`:

  Number of permutations or surrogates, the same at every lag (the exact
  count for the cyclic scheme).

- `profile`:

  Data frame with `lag`, `n` (rows used), `m2`, `r`, `p.value`,
  `p.adjusted`, `nperm`, `z` and `a` for every lag.

- `residuals`:

  Data frame of per-observation Procrustes residuals at the optimal lag,
  with the paired row labels of `X` and `Y`.

- `procrustes`:

  The
  [`vegan::procrustes()`](https://vegandevs.github.io/vegan/reference/procrustes.html)
  fit at the optimal lag.

- `effect`:

  List with the effect `lag` (largest \\z\\), `Z`, `A` and the
  chance-corrected `RV` coefficient at that lag (see Details).

- `ar1`:

  Data frame with the mean lag-1 autocorrelation of the columns of `X`
  and `Y`, `raw` and `residual` after removing the harmonics of period
  `n` and `period`.

- `null`:

  List of permuted \\m^2\\ values, one element per lag; element `i` of
  every vector comes from the same permutation.

- `null.rv`:

  List of permuted RV coefficients, paired with `null`.

- `permutation`, `mirror`, `p.adjust.method`, `n`, `period`, `call`:

  Settings used; `mirror` is `NA` unless the scheme is cyclic, and `n`
  is the number of rows after any truncation to whole cycles.

## Details

Rows of `X` and `Y` are ordered time points (e.g., monthly centroids)
and must correspond one to one. At lag `k > 0`, `X` at time `t` is
paired with `Y` at time `t + k` (X leads Y), and the series are
truncated to the overlapping window of `n - k` rows. Negative lags (Y
leads X) are tested when `bidirectional = TRUE`. The statistic at each
lag is the symmetric Procrustes sum of squares \\m^2\\ from
[`vegan::procrustes()`](https://vegandevs.github.io/vegan/reference/procrustes.html);
the lag of maximal concordance is the one with the smallest \\m^2\\.

Three permutation schemes are available:

- `"phase"` (default): `nperm` phase-randomized surrogates of `Y` are
  drawn. Each surrogate keeps the amplitude spectrum of every column of
  `Y` (hence its autocorrelation) and the cross-spectra among its
  columns, but receives new Fourier phases, shared by all columns.
  Surrogates are applied jointly to all lags and do not reproduce the
  observed alignment, which makes this the recommended scheme for the
  global p-value. It assumes a stationary `Y` and treats the series as
  periodic with period `n`. In simulations, it was liberal for short,
  strongly autocorrelated series (AR(1) of 0.7 or more with `n` of 48 or
  less).

- `"cyclic"`: the full `Y` series is cyclically shifted by
  `1, ..., n - 1` positions before lagging, and optionally also
  time-reversed (`mirror = TRUE`). This keeps the temporal
  autocorrelation of `Y` under the null. All `n - 1` (or `2(n - 1)`)
  shifts are enumerated, so `nperm` and `seed` are ignored, and the
  smallest attainable p-value is `1 / (2(n - 1) + 1)` (0.043 for
  `n = 12` with `mirror = TRUE`). The scheme assumes the series spans
  complete cycles, i.e., that its period matches the number of rows
  (e.g., 12 monthly samples of an annual cycle). Because a shift of `Y`
  combined with lag `k` is equivalent to another lag, some shifts
  reproduce the observed alignment at other lags; the global p-value
  (see below) is then valid but has low power.

- `"free"`: rows of `Y` are permuted at random within the overlapping
  window, as in
  [`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html).
  This assumes exchangeable rows and can strongly inflate Type I error
  when both series are autocorrelated (Manly and Navarro Alberto 2020,
  ch. 10).

Both the phase and the cyclic schemes treat the series as periodic with
period `n`. When the series contains a seasonal cycle whose length does
not divide `n` (e.g., 18 monthly samples of an annual cycle), no scheme
is calibrated, because two series sharing a period are
Procrustes-concordant whatever their phase. Supplying `period` truncates
both series to their first `floor(n / period) * period` rows, so that
the analysis spans whole cycles. In simulations, this restored the Type
I error of series that already spanned complete cycles, whereas removing
the seasonal component by harmonic regression did not.

P-values are computed as `(b + 1) / (N + 1)`, where `b` is the number of
permuted \\m^2\\ values smaller than or equal to the observed one and
`N` the number of permutations (Manly and Navarro Alberto 2020, p. 57).
Each permutation is applied jointly to all lags. The validity of the
test does not depend on `N`, only the precision of the p-value does;
Manly and Navarro Alberto (2020, pp. 57–59) recommend at least 1000
randomizations for a test at the 5% level and about 5000 at the 1%
level. A message is issued when the smallest attainable p-value,
`1 / (N + 1)`, is 0.05 or more.

Because the p-value is bounded by `N` and conflates strength with sample
size, three continuous measures are computed from the same permutations.
At each lag, \\z_k = (\bar m^{2\*}\_k - m^2_k) / s^\*\_k\\ standardizes
the observed \\m^2\\ against its permutation mean and standard
deviation, and \\a_k = 1 - m^2_k / \bar m^{2\*}\_k\\ expresses the
concordance relative to its null mean. The effect lag is the lag of the
largest \\z\\ (which weighs lags of different overlap more evenly than
the raw \\m^2\\), and each permutation selects its own lag in the same
way. `Z` is the standardized position of the observed maximum \\z\\
within the permuted maxima; `A` and the corrected RV coefficient are
chance-corrected at the effect lag as \\(obs - \bar{null}) / (1 -
\bar{null})\\, with the null mean taken over the permutations at their
own selected lags. `A` and the corrected RV are near 0 under the null
and 1 for perfect concordance; `Z` is near 0 under the null.

The mean lag-1 autocorrelation of the columns of each series is returned
in `ar1`, both raw and after removing by regression the harmonics of
period `n` (and of `period`, when given). It is descriptive only: the
raw value is dominated by any seasonal cycle, and the residual value is
biased downwards in short series. The phase scheme was liberal for
short, strongly autocorrelated series in simulations, so both values
help judge whether that caveat applies.

The per-lag p-values (`p.value`) are not corrected for testing several
lags, nor for selecting the lag with the smallest \\m^2\\. Reporting the
raw p-value at the optimal lag therefore inflates Type I error roughly
to the rate at which any lag is significant. The global p-value
(`p.global`) accounts for both: it compares the observed minimum \\m^2\\
across lags with the distribution of the minimum across lags under each
permutation (a max-statistic test), so it is the p-value of the whole
scan-and-select procedure and is the one to report for the optimal lag.
Because \\m^2\\ at neighboring lags is strongly correlated, a Bonferroni
level of \\\alpha / m\\ for \\m\\ lags would be conservative: the
appropriate level lies between \\\alpha / m\\ and \\\alpha\\, and taking
the extreme statistic across tests in each randomization, as here, finds
it from the data (Manly and Navarro Alberto 2020, pp. 83–84). The same
source recommends at least \\100m\\ randomizations for \\m\\ tests. With
`p.adjust.method = "maxT"`, every lag receives the corresponding
single-step max-statistic adjusted p-value, and the adjusted p-value at
the optimal lag equals `p.global`. Any method of
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) can be
requested instead; with few permutations, as in the cyclic scheme,
corrections such as Holm may make significance unattainable.

\\m^2\\ values at different lags are based on different numbers of rows
(`n - |k|`), and \\m^2\\ is smaller and more variable when fewer rows
overlap. The lag of the smallest \\m^2\\ therefore tends to fall at the
longest lags scanned, even when the series are independent: in 200 pairs
of independent series simulated with
[`simulate_lagged_series()`](https://provetelab.org/lagprotest/reference/simulate_lagged_series.md)
(`n_times = 12`, lags -3 to 3), it fell at \\\pm 3\\ in about half of
the replicates rather than in two of seven. The global p-value remains
valid, but has more power for long lags than for short ones. The effect
lag (largest \\z\\) places every lag on the scale of its own null
distribution and fell at \\\pm 3\\ in about a quarter of the same
replicates. When the two lags differ,
[`print()`](https://rdrr.io/r/base/print.html) says so; an optimal lag
at the end of the scan should be read with caution, and candidate lags
kept within a range justified by the biology of the system.

## References

Manly, B. F. J. and Navarro Alberto, J. A. (2020). *Randomization,
Bootstrap and Monte Carlo Methods in Biology*, 4th ed. Chapman and
Hall/CRC, New York.
[doi:10.1201/9780429329203](https://doi.org/10.1201/9780429329203)

## See also

[`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html),
[`vegan::procrustes()`](https://vegandevs.github.io/vegan/reference/procrustes.html)

## Examples

``` r
sim <- simulate_lagged_series(n_times = 12, true_lag = 2, noise_sd = 0.3,
                              seed = 1)
fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, seed = 1)
fit
#> 
#> Lagged PROTEST
#> 
#> Call: lagged_protest(X = sim$X, Y = sim$Y, lags = 0:4, seed = 1)
#> Permutation: phase-randomized surrogates; n = 12 time points
#> P-value adjustment across lags: none
#> 
#> Optimal lag: 2  (m2 = 0.0496, r = 0.975, p = 0.001)
#> Global p across lags (max-statistic): 0.004
#> Effect at lag 2 (largest z): Z = 3.04, A = 0.464, corrected RV = 0.206
#> Mean lag-1 autocorrelation, raw (residual after removing the cycle): X 0.76 (-0.018), Y 0.66 (-0.049)
#> 
#> Lag profile:
#>  lag  n     m2     r p.value nperm      z      a
#>    0 12 0.1532 0.920   0.969   999 -1.705 -0.227
#>    1 11 0.1000 0.949   0.145   999  1.112  0.179
#>    2 10 0.0496 0.975   0.001   999  3.300  0.578
#>    3  9 0.1303 0.933   0.725   999 -0.592 -0.122
#>    4  8 0.1516 0.921   0.850   999 -1.025 -0.246
#> Per-lag p-values are not corrected for testing several lags or for selecting
#> the optimal lag; report the global p for the optimal lag.
fit$p.global
#> [1] 0.004
fit_cyclic <- lagged_protest(sim$X, sim$Y, lags = 0:4,
                             permutation = "cyclic")
fit_cyclic$profile
#>   lag  n         m2         r    p.value p.adjusted nperm          z          a
#> 1   0 12 0.15315445 0.9202421 1.00000000 1.00000000    22 -1.6979953 -0.2233912
#> 2   1 11 0.10003719 0.9486637 0.08695652 0.08695652    22  1.2470776  0.1865535
#> 3   2 10 0.04963006 0.9748692 0.04347826 0.04347826    22  3.4786381  0.5977914
#> 4   3  9 0.13033905 0.9325561 0.73913043 0.73913043    22 -0.4396758 -0.1086388
#> 5   4  8 0.15160634 0.9210829 0.91304348 0.91304348    22 -0.9413089 -0.2284317
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
