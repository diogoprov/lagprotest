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
  stat = c("minP", "z", "raw"),
  p.adjust.method = "none",
  period = NULL,
  deseason = c("none", "monthly"),
  prewhiten = c("none", "prais"),
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
  Alberto 2020, pp. 57–59). With `stat = "minP"`, use at least 100 per
  lag tested (for example, 699 or more for lags -3 to 3), otherwise the
  global p-value may be unable to fall below 0.05 (see Details).

- mirror:

  Logical; with `permutation = "cyclic"`, also use the time-reversed
  shifts of `Y`, doubling the permutation set.

- bidirectional:

  Logical; if `TRUE`, also test the negative lags `-lags` (Y leads X).

- scale:

  Logical; z-score the columns of `X` and `Y` (computed on the full
  series) before analysis. Ignored for `dist` inputs.

- stat:

  Scan statistic used for the global p-value and the optimal lag:
  `"minP"` (default), `"z"` or `"raw"` (see Details).

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

- deseason:

  `"none"` (default) or `"monthly"`, to subtract from both series the
  mean of each position in the cycle of length `period` (see Details).

- prewhiten:

  `"none"` (default) or `"prais"`, for Prais-Winsten AR(1) prewhitening
  of every column of both series, applied after `deseason` (see
  Details).

- seed:

  Optional integer seed for the `"phase"` and `"free"` schemes. The
  global random number stream is left untouched.

## Value

An object of class `"lagged_protest"`, a list with:

- `optimal_lag`:

  Lag selected by `stat` (see Details).

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

  P-value of the scan across lags for `stat`, accounting for testing
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
  and `Y`, as analyzed (after any preprocessing), `raw` and `residual`
  after removing the harmonics of period `n` and `period`.

- `null`:

  List of permuted \\m^2\\ values, one element per lag; element `i` of
  every vector comes from the same permutation.

- `null.rv`:

  List of permuted RV coefficients, paired with `null`.

- `rho`:

  With `prewhiten = "prais"`, a list with the AR(1) coefficients used
  for the columns of `X` and `Y`; otherwise `NULL`.

- `permutation`, `stat`, `mirror`, `p.adjust.method`, `n`, `period`,
  `deseason`, `prewhiten`, `call`:

  Settings used; `mirror` is `NA` unless the scheme is cyclic, and `n`
  is the number of rows after any truncation to whole cycles.

## Details

Rows of `X` and `Y` are ordered time points (e.g., monthly centroids)
and must correspond one to one. At lag `k > 0`, `X` at time `t` is
paired with `Y` at time `t + k` (X leads Y), and the series are
truncated to the overlapping window of `n - k` rows. Negative lags (Y
leads X) are tested when `bidirectional = TRUE`. The statistic at each
lag is the symmetric Procrustes sum of squares \\m^2\\ from
[`vegan::procrustes()`](https://vegandevs.github.io/vegan/reference/procrustes.html).
The optimal lag is chosen by the scan statistic set by `stat` (by
default the lag with the smallest per-lag p-value; see below).

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
  less). It also assumes that any seasonal cycle shared by `X` and `Y`
  is close to sinusoidal: new phases destroy the alignment among the
  harmonics of an abrupt cycle (a sawtooth or a short pulse), so the
  surrogates no longer share the shape of the cycle, and the test can be
  strongly liberal (Type I error of 0.10 to 0.34 in simulations with
  `simulate_lagged_series(waveform = "sawtooth")` or `"pulse"`, 24 or 48
  monthly samples). When the series span several whole cycles, removing
  from both the mean of each position in the cycle (e.g., the mean of
  each calendar month; `deseason = "monthly"`) restored the Type I error
  in those simulations.

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

- `"free"`: rows of `Y` are permuted at random, as in
  [`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html),
  and the permuted series is then lagged like the observed one. This
  assumes exchangeable rows and can strongly inflate Type I error when
  both series are autocorrelated (Manly and Navarro Alberto 2020, ch.
  10).

Both the phase and the cyclic schemes treat the series as periodic with
period `n`. When the series contains a seasonal cycle whose length does
not divide `n` (e.g., 18 monthly samples of an annual cycle), no scheme
is calibrated, because two series sharing a period are
Procrustes-concordant whatever their phase. Supplying `period` truncates
both series to their first `floor(n / period) * period` rows, so that
the analysis spans whole cycles. In simulations, this restored the Type
I error of series that already spanned complete cycles, whereas removing
the seasonal component by harmonic regression did not.

Two optional preprocessing steps are applied to both series, in this
order, after truncation to whole cycles and before the columns are
standardized (`scale`):

- `deseason = "monthly"` subtracts from every column the mean of each
  position in the cycle (e.g., the mean of each calendar month), so that
  the test addresses the anomalies rather than a seasonal cycle that
  both series may share. It requires `period` and at least two whole
  cycles; the first row is taken as position 1 of the cycle. In
  simulations, it restored the Type I error of the phase scheme for
  shared seasonal cycles of any shape (see the caveat on non-sinusoidal
  cycles above), whereas removing two harmonics by regression did not
  suffice for a sawtooth. The scan then measures concordance among
  anomalies only. With strongly autocorrelated noise (AR(1) of 0.9), it
  remained slightly liberal (up to about 0.08).

- `prewhiten = "prais"` removes AR(1) autocorrelation from every column
  by the Prais-Winsten transformation, \\x_1\sqrt{1 - \rho^2}\\ for the
  first row and \\x_t - \rho x\_{t-1}\\ for the others, with \\\rho\\
  the lag-1 autocorrelation of the column (returned in `rho`). Unlike
  ordinary prewhitening, it keeps all `n` rows: dropping the first row
  after removing seasonal means made the test strongly liberal in
  simulations, and so did prewhitening before removing them. For
  strongly autocorrelated series without seasonality, prewhitening
  brought the Type I error of the phase scheme from 0.10–0.19 back to
  the nominal level, and it increased power when the anomalies were
  strongly autocorrelated. Combined with `deseason`, it can be liberal
  when an abrupt seasonal cycle varies in amplitude among years and the
  noise is strongly autocorrelated (about 0.10 in simulations with a
  sawtooth cycle and AR(1) of 0.9), where `deseason` alone was not. A
  warning is issued when it is requested for a seasonal series (`period`
  given) without `deseason`, a combination that inflated the Type I
  error in simulations.

With a single cycle (`n` equal to `period`), seasonal means cannot be
separated from the signal; a message then recalls that the phase scheme
assumes a close-to-sinusoidal shared cycle.

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
and 1 for perfect concordance; `Z` is near 0 under the null. Each
permuted \\z\\ is computed in the same way as the observed one, against
the mean and standard deviation of the other values at that lag (the
remaining permutations and the observed), so that the observed and the
permuted values are treated alike.

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
(`p.global`) accounts for both: it compares the most extreme observed
score across lags with the distribution of the most extreme score across
lags under each permutation (a max-statistic test), so it is the p-value
of the whole scan-and-select procedure and is the one to report for the
optimal lag. Because \\m^2\\ at neighboring lags is strongly correlated,
a Bonferroni level of \\\alpha / m\\ for \\m\\ lags would be
conservative: the appropriate level lies between \\\alpha / m\\ and
\\\alpha\\, and taking the extreme statistic across tests in each
randomization, as here, finds it from the data (Manly and Navarro
Alberto 2020, pp. 83–84). The same source recommends at least \\100m\\
randomizations for \\m\\ tests. With `p.adjust.method = "maxT"`, every
lag receives the corresponding single-step adjusted p-value for the
chosen `stat`, and the adjusted p-value at the optimal lag equals
`p.global`. Any method of
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) can be
requested instead; with few permutations, as in the cyclic scheme,
corrections such as Holm may make significance unattainable.

The score compared across lags is set by `stat`:

- `"minP"` (default): the per-lag p-value. The global statistic is the
  smallest p-value across lags, and its null distribution is obtained by
  converting every permuted \\m^2\\ to a p-value within its own lag, as
  its rank among the observed and permuted values divided by \\N + 1\\,
  exactly as the observed p-value is computed (the minP procedure of
  Westfall and Young 1993; Ge et al. 2003). The optimal lag is the one
  with the smallest p-value; ties, common with few permutations, are
  broken by \\z\\ and then by \\m^2\\.

- `"z"`: \\m^2\\ centered and scaled by the mean and standard deviation
  of its own null distribution, with each permuted value standardized
  against the others and the observed (see above). The optimal lag is
  the one with the largest \\z\\ and equals the effect lag.

- `"raw"`: \\m^2\\ itself, the statistic of versions up to 0.6. The
  optimal lag is the one with the smallest \\m^2\\.

All three are computed from the same permutations and give valid global
tests. They differ in how power is spread across lags, because \\m^2\\
values at different lags are based on different numbers of rows
(`n - |k|`), and \\m^2\\ is smaller and more variable when fewer rows
overlap. Under `"raw"`, the optimal lag therefore tends to fall at the
longest lags scanned, even when the series are independent: in 200 pairs
of independent series simulated with
[`simulate_lagged_series()`](https://provetelab.org/lagprotest/reference/simulate_lagged_series.md)
(`n_times = 12`, lags -3 to 3), it fell at \\\pm 3\\ in about half of
the replicates rather than in two of seven, and the test had less power
for a true lag of zero. `"minP"` and `"z"` place every lag on the scale
of its own null distribution, so the optimal lag is close to uniform
across lags under independence. Because the observed and the permuted
scores are computed in the same way, all three give exact tests when the
permutations are exchangeable with the observed series (free and cyclic
schemes). `"minP"` is conservative when the permutations are few
relative to the number of lags, because per-lag p-values are discrete,
and `"z"` then has more power. The same reasoning motivates selecting
among candidates by p-value in conditional inference trees (Hothorn et
al. 2006). Because per-lag p-values are discrete, the smallest global
p-value that `"minP"` can reach grows with the number of lags: roughly
\\(m + 1) / (N + 1)\\ for \\m\\ lags and \\N\\ permutations, so \\N\\
should be at least \\100m\\ (Manly and Navarro Alberto 2020, p. 84), and
a message is issued when p \< 0.05 cannot be reached. With the cyclic
scheme, permutations are few, and `"z"` is then preferable. When the
optimal lag and the effect lag differ,
[`print()`](https://rdrr.io/r/base/print.html) says so; candidate lags
should be kept within a range justified by the biology of the system.

## References

Ge, Y., Dudoit, S. and Speed, T. P. (2003). Resampling-based multiple
testing for microarray data analysis. *Test* 12, 1–77.
[doi:10.1007/BF02595811](https://doi.org/10.1007/BF02595811)

Hothorn, T., Hornik, K. and Zeileis, A. (2006). Unbiased recursive
partitioning: a conditional inference framework. *Journal of
Computational and Graphical Statistics* 15, 651–674.
[doi:10.1198/106186006X133933](https://doi.org/10.1198/106186006X133933)

Manly, B. F. J. and Navarro Alberto, J. A. (2020). *Randomization,
Bootstrap and Monte Carlo Methods in Biology*, 4th ed. Chapman and
Hall/CRC, New York.
[doi:10.1201/9780429329203](https://doi.org/10.1201/9780429329203)

Westfall, P. H. and Young, S. S. (1993). *Resampling-Based Multiple
Testing: Examples and Methods for p-Value Adjustment*. Wiley, New York.

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
#> Global p across lags (minP scan): 0.005
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
#> [1] 0.005
fit_cyclic <- lagged_protest(sim$X, sim$Y, lags = 0:4,
                             permutation = "cyclic", stat = "z")
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

# Four years of monthly data sharing an abrupt seasonal cycle: test the
# anomalies, with AR(1) prewhitening
saw <- simulate_lagged_series(n_times = 48, true_lag = 2,
                              waveform = "sawtooth", noise_sd = 1,
                              seed = 1)
lagged_protest(saw$X, saw$Y, lags = 0:4, period = 12,
               deseason = "monthly", prewhiten = "prais", seed = 1)
#> 
#> Lagged PROTEST
#> 
#> Call: lagged_protest(X = saw$X, Y = saw$Y, lags = 0:4, period = 12, deseason = "monthly", 
#>     prewhiten = "prais", seed = 1)
#> Permutation: phase-randomized surrogates; n = 48 time points (whole cycles of period 12)
#> Preprocessing: seasonal means removed, then Prais-Winsten AR(1) prewhitening
#> P-value adjustment across lags: none
#> 
#> Optimal lag: 2  (m2 = 0.507, r = 0.702, p = 0.001)
#> Global p across lags (minP scan): 0.005
#> Effect at lag 2 (largest z): Z = 11.5, A = 0.404, corrected RV = 0.425
#> Mean lag-1 autocorrelation, raw (residual after removing the cycle): X 0.024 (-0.046), Y 0.057 (-0.015)
#> 
#> Lag profile:
#>  lag  n    m2     r p.value nperm        z         a
#>    0 48 0.916 0.290   0.708   999 -0.59714 -0.020735
#>    1 47 0.922 0.280   0.811   999 -0.85029 -0.030585
#>    2 46 0.507 0.702   0.001   999 11.56998  0.431889
#>    3 45 0.892 0.329   0.442   999  0.00522  0.000194
#>    4 44 0.889 0.333   0.454   999  0.01432  0.000521
#> Per-lag p-values are not corrected for testing several lags or for selecting
#> the optimal lag; report the global p for the optimal lag.
```
