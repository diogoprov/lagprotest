# Getting started with lagprotest

``` r

library(lagprotest)
```

## The question

PROTEST (Jackson 1995; Peres-Neto and Jackson 2001) measures how well
two multivariate configurations of the same observations match after
Procrustes superimposition, and tests that match by permutation. When
the observations are ordered time points, for example monthly means of a
set of environmental variables and of a set of phenotypic traits,
concordance may be strongest not between contemporaneous rows but after
a delay: the response may track its driver with a lag of *k* time units.

[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
asks at which lag the two series match best. For each candidate lag *k*,
row *t* of `X` is paired with row *t + k* of `Y`, the series are
truncated to their `n - k` overlapping rows, and the symmetric
Procrustes statistic *m*² is computed. Lower *m*² means stronger
concordance. Each lag is tested by permutation, and the lag with the
strongest evidence of concordance (by default, the smallest per-lag
p-value) is reported together with the full lag profile.

## Input data

Both inputs have one row per time point, in temporal order, and the same
number of rows. Columns are variables and need not be the same in number
or kind. `X` is the configuration expected to lead (the driver); `Y` the
one expected to follow. Inputs can be matrices, data frames or `dist`
objects; dissimilarities are converted to principal coordinates, so any
dissimilarity measure suited to the data can be used upstream.

Here we simulate a 12-month series of five variables (`X`) and a
response (`Y`) that follows `X` with a lag of two months:

``` r

sim <- simulate_lagged_series(n_times = 12, true_lag = 2, autocorr = 0.5,
                              noise_sd = 0.5, seed = 1)
dim(sim$X)
#> [1] 12  5
dim(sim$Y)
#> [1] 12  5
```

[`simulate_lagged_series()`](https://provetelab.org/lagprotest/reference/simulate_lagged_series.md)
builds seasonal sinusoids plus AR(1) noise. With `true_lag = NULL` it
returns two independent series, which is useful for checking Type I
error.

## Running the test

``` r

fit <- lagged_protest(sim$X, sim$Y, lags = 0:4, nperm = 999, seed = 1)
fit
#> 
#> Lagged PROTEST
#> 
#> Call: lagged_protest(X = sim$X, Y = sim$Y, lags = 0:4, nperm = 999, seed = 1)
#> Permutation: phase-randomized surrogates; n = 12 time points
#> P-value adjustment across lags: none
#> 
#> Optimal lag: 2  (m2 = 0.0922, r = 0.953, p = 0.001)
#> Global p across lags (minP scan): 0.006
#> Effect at lag 2 (largest z): Z = 3.22, A = 0.51, corrected RV = 0.318
#> Mean lag-1 autocorrelation, raw (residual after removing the cycle): X 0.7 (-0.018), Y 0.54 (-0.049)
#> 
#> Lag profile:
#>  lag  n     m2     r p.value nperm      z      a
#>    0 12 0.3103 0.830   0.962   999 -1.550 -0.201
#>    1 11 0.2078 0.890   0.138   999  1.107  0.173
#>    2 10 0.0922 0.953   0.001   999  3.533  0.617
#>    3  9 0.2748 0.852   0.818   999 -0.878 -0.171
#>    4  8 0.3046 0.834   0.876   999 -1.181 -0.262
#> Per-lag p-values are not corrected for testing several lags or for selecting
#> the optimal lag; report the global p for the optimal lag.
```

By default the null distribution comes from phase-randomized surrogates
of `Y` (see [Choosing the permutation
scheme](#choosing-the-permutation-scheme)); `seed` makes the result
reproducible. The profile lists, for every lag, the number of rows used
(`n`), *m*², the Procrustes correlation *r* = √(1 − *m*²), the
permutation p-value and the number of surrogates. The optimal lag here
is the true one, 2.

The per-lag p-values are not corrected for testing several lags, nor for
choosing the optimal lag after looking at the profile. Reporting the raw
p-value at the optimal lag therefore inflates Type I error to roughly
the rate at which *any* lag is significant. The p-value to report for
the optimal lag is the global p-value, printed below the optimal lag and
stored in `fit$p.global`. It compares the smallest per-lag p-value with
the smallest per-lag p-value obtained under each permutation (the minP
test of Westfall and Young), so it is the p-value of the whole
scan-and-select procedure:

``` r

fit$p.global
#> [1] 0.006
```

Comparing lags by their p-values, rather than by *m*² itself, matters
because *m*² at lag *k* rests on `n − k` rows and is smaller and more
variable when fewer rows overlap. The smallest *m*² therefore tends to
fall at the longest lags even for independent series, and a test built
on it has less power for short lags. The argument `stat` chooses the
statistic: `"minP"` (default), `"z"` (*m*² standardized by its own null
distribution) or `"raw"` (the smallest *m*², used up to version 0.6).
All three are valid global tests computed from the same permutations:

``` r

sapply(c("minP", "z", "raw"), function(s)
  lagged_protest(sim$X, sim$Y, lags = 0:4, seed = 1, stat = s)$p.global)
#>  minP     z   raw 
#> 0.006 0.002 0.003
```

Per-lag p-values can also be corrected. `p.adjust.method = "maxT"` gives
the single-step adjusted p-value of every lag for the chosen `stat`,
whose value at the optimal lag equals the global p-value; any method of
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) can be used
instead, and raw and adjusted values are then both reported:

``` r

fit_holm <- lagged_protest(sim$X, sim$Y, lags = 0:4, permutation = "cyclic",
                           stat = "z", p.adjust.method = "holm")
fit_holm$profile[, c("lag", "p.value", "p.adjusted")]
#>   lag    p.value p.adjusted
#> 1   0 0.95652174  1.0000000
#> 2   1 0.08695652  0.3478261
#> 3   2 0.04347826  0.2173913
#> 4   3 0.73913043  1.0000000
#> 5   4 0.91304348  1.0000000
```

With the cyclic scheme and `n = 12` there are only 22 permutations, so a
raw p-value at its minimum of 0.043 becomes 0.217 after Holm correction
across five lags, and significance at α = 0.05 is unattainable (see
Caveats).

[`summary()`](https://rdrr.io/r/base/summary.html) returns the profile
as a data frame, with the optimal lag flagged, and prints the global
p-value and the effect measures described next below it:

``` r

summary(fit)
#>  lag  n     m2     r p.value p.adjusted nperm      z      a optimal
#>    0 12 0.3103 0.830   0.962      0.962   999 -1.550 -0.201   FALSE
#>    1 11 0.2078 0.890   0.138      0.138   999  1.107  0.173   FALSE
#>    2 10 0.0922 0.953   0.001      0.001   999  3.533  0.617    TRUE
#>    3  9 0.2748 0.852   0.818      0.818   999 -0.878 -0.171   FALSE
#>    4  8 0.3046 0.834   0.876      0.876   999 -1.181 -0.262   FALSE
#> 
#> Optimal lag 2, global p = 0.006
#> Effect at lag 2 (largest z): Z = 3.22, A = 0.51, corrected RV = 0.318
#> Mean lag-1 autocorrelation, raw (residual): X 0.7 (-0.018), Y 0.54 (-0.049)
```

## Strength of the association

A p-value is bounded by the number of permutations and mixes the
strength of an association with the length of the series.
[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
therefore also returns continuous measures, computed from the same
permutations:

``` r

fit$effect
#> $lag
#> [1] 2
#> 
#> $Z
#> [1] 3.217432
#> 
#> $A
#> [1] 0.5099893
#> 
#> $RV
#> [1] 0.3180133
```

At each lag, `z` in the profile is the observed *m*² standardized by the
mean and standard deviation of its permutation null, with the sign
reversed so that larger values mean stronger concordance, and `a` is the
concordance relative to the null mean. The effect lag is the lag of the
largest `z`. `Z` is the position of that largest `z` among the largest
`z` of each permutation, in standard deviations; `A` and the corrected
RV coefficient are scaled so that 0 is the value expected by chance and
1 is perfect concordance. All three are near 0 for independent series.

The effect lag and the optimal lag usually agree, because both place
every lag on the scale of its own null distribution; with `stat = "z"`
they are the same by definition. With `stat = "raw"` they differ more
often, because the smallest *m*² tends to fall at the longest lags.
[`print()`](https://rdrr.io/r/base/print.html) notes when the two differ
(see Caveats).

`fit$ar1` gives the mean lag-1 autocorrelation of the columns of each
series, raw and after removing the seasonal harmonics. It is
descriptive, and helps judge whether the series are short and strongly
autocorrelated enough for the phase scheme to be liberal (see Caveats):

``` r

fit$ar1
#>   series       raw    residual
#> 1      X 0.6991221 -0.01815912
#> 2      Y 0.5396661 -0.04854199
```

## Plotting the lag profile

``` r

plot(fit)
```

![Lag profile: m2 against lag, with a dashed line for the 5% quantile of
the permutation null distribution at each
lag.](lagprotest_files/figure-html/plot-1.png)

The solid line is the observed *m*² at each lag. The dashed line is the
5% quantile of the permutation null distribution at that lag: observed
values below it are significant at α = 0.05 before correction. Filled
points mark lags with p \< 0.05 (after adjustment, if one was
requested), the dotted vertical line marks the optimal lag, and the
circle around it is filled when the global p-value is below 0.05.
`ggplot2::autoplot(fit)` draws the same plot with ggplot2, and
`as.data.frame(fit)` returns the data behind it.

## Residuals per observation

Procrustes residuals show which time points fit the superimposition well
and which deviate, at the optimal lag:

``` r

residuals(fit)
#>    x_row y_row   residual
#> 1      1     3 0.09796321
#> 2      2     4 0.11137958
#> 3      3     5 0.05870125
#> 4      4     6 0.09649460
#> 5      5     7 0.07735381
#> 6      6     8 0.06387173
#> 7      7     9 0.13230232
#> 8      8    10 0.09243341
#> 9      9    11 0.10137702
#> 10    10    12 0.10519784
```

Each row pairs a time point of `X` (`x_row`) with the time point of `Y`
it is matched to (`y_row`). Row names of the inputs, such as month
labels, are carried through.

## Choosing the permutation scheme

Standard PROTEST permutes rows freely, which assumes that time points
are exchangeable under the null hypothesis. Seasonal or otherwise
autocorrelated series are not: two independent series that both cycle
annually will often look concordant against a null that destroys their
temporal structure, a problem shared with other permutation tests on
autocorrelated data (Legendre et al. 2002).

[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
offers two schemes that keep the temporal structure of `Y`:

- `permutation = "phase"` (default) draws `nperm` phase-randomized
  surrogates of `Y` (Ebisuzaki 1997). Each surrogate keeps the amplitude
  spectrum of every column of `Y`, and hence its autocorrelation, as
  well as the cross-spectra among the columns, but receives new Fourier
  phases.
- `permutation = "cyclic"` shifts the full `Y` series cyclically by 1,
  …, n − 1 positions, and also time-reverses it when `mirror = TRUE`,
  *before* the lag is applied. Each permuted series keeps its
  autocorrelation. All shifts are enumerated, so the test needs no
  random seed.

`permutation = "free"` reproduces the behavior of
[`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html)
(Oksanen et al. 2026).

The cyclic and phase schemes differ in how they handle the scan across
lags. Shifting `Y` by *s* positions and then applying lag *k* pairs the
same rows as lag *k* + *s* without a shift. When several lags are
scanned, some cyclic permutations therefore reproduce the observed
alignment at another lag, and the null distribution of the minimum *m*²
across lags contains the signal itself. The per-lag p-values are
unaffected, but the global p-value loses power. Phase surrogates are new
series rather than shifts of the observed one, so they avoid this
problem; they also allow as many surrogates as desired, instead of the
`2(n − 1)` cyclic permutations.

The choice matters. Under independence, a small simulation with 24
monthly samples spanning one cycle gives the following rejection rates
at α = 0.05: for the raw p-value at the optimal lag, and for the global
p-value under each scheme. A second set of series, in which `Y` follows
`X` with a lag of two, gives the power of the global p-value. The
simulation uses `stat = "raw"` because with only 99 surrogates and five
lags `"minP"` cannot reach p \< 0.05 (it needs at least 100 permutations
per lag):

``` r

rates <- function(true_lag, n_rep) {
  res <- vapply(seq_len(n_rep), function(i) {
    s <- simulate_lagged_series(n_times = 24, period = 24, true_lag = true_lag,
                                autocorr = 0.5, noise_sd = 0.5, seed = i)
    cyc <- lagged_protest(s$X, s$Y, lags = 0:4, permutation = "cyclic",
                          stat = "raw")
    pha <- lagged_protest(s$X, s$Y, lags = 0:4, permutation = "phase",
                          nperm = 99, seed = i, stat = "raw")
    fre <- lagged_protest(s$X, s$Y, lags = 0:4, permutation = "free",
                          nperm = 99, seed = i, stat = "raw")
    c(cyclic_optimal_raw = cyc$p.value, cyclic_global = cyc$p.global,
      phase_global = pha$p.global, free_global = fre$p.global)
  }, numeric(4))
  rowMeans(res < 0.05)
}
rbind(type_I = rates(NULL, 100), power = rates(2, 100))
#>        cyclic_optimal_raw cyclic_global phase_global free_global
#> type_I                0.1          0.03         0.03           1
#> power                 1.0          0.54         1.00           1
```

The raw p-value at the optimal lag is too liberal and free permutation
rejects almost always. Both global p-values keep Type I error near the
nominal level, but only the phase scheme retains power.

## Caveats

- **Complete cycles and shared periodicity.** Both the cyclic and the
  phase schemes treat the series as periodic with period `n`. This is
  appropriate when the series spans whole cycles, such as 12 monthly
  samples of an annual cycle. Cyclic-shift null models can be biased
  when this does not hold (Kalyuzhny 2020). In simulations with this
  package, two independent seasonal series whose period differed from
  the series length were declared concordant far more often than α under
  either scheme: series with the same period trace the same trajectory
  shape, which Procrustes rotation can superimpose whatever the phase
  difference. Supplying the length of the cycle with `period` truncates
  both series to whole cycles (e.g., the first 18 of 24 monthly samples
  for an 18-month cycle), which brought the Type I error back to the
  level of series that span complete cycles, at the cost of the rows
  dropped. Removing the seasonal component by harmonic regression did
  not: it reduced the excess but left the test liberal, and it removes
  the seasonal signal that is often the object of the analysis.
- **Strong autocorrelation.** In simulations with 1000 replicates per
  scenario, the global p-value of the phase scheme kept Type I error at
  0.047–0.053 for AR(1) autocorrelation up to 0.5, but reached 0.063 at
  0.7 and 0.129 at 0.9 (averaged over series length and number of
  variables), and 0.15–0.19 with 12–24 time points or 10 variables. The
  excess faded with longer series (0.063 at AR(1) = 0.9 with 96 time
  points). The cyclic scheme without mirror was less affected (0.081 at
  0.9) but had much lower power. For short, strongly autocorrelated
  series, report the cyclic result alongside the phase one, and estimate
  the autocorrelation of the series.
- **Few permutations.** With `n` time points there are only `2(n − 1)`
  cyclic permutations, so the smallest attainable p-value is
  `1 / (2(n − 1) + 1)`, 0.043 for `n = 12`. Corrections across many lags
  may then be unattainable at α = 0.05. The phase scheme has no such
  floor.
- **Different window sizes.** *m*² at lag *k* is based on `n − k` rows,
  so long lags rest on fewer observations, and the smallest *m*² tends
  to fall at the longest lags scanned even for independent series. The
  default `stat = "minP"` removes this bias from the global test and the
  optimal lag; with `stat = "raw"`, an optimal lag at the end of the
  scan should be compared with the effect lag in `fit$effect`. In all
  cases, keep candidate lags within a range justified by the biology of
  the system.
- **Scaling.** With `scale = TRUE` (default) variables are standardized
  on the full series before truncation.

## Several data sets

[`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md)
runs
[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
on every pair of data sets in a named list (or on chosen pairs),
scanning both directions through the sign of the lag, and adjusts the
global p-values across pairs. The article [Seasonal liver plasticity in
a
treefrog](https://provetelab.org/lagprotest/articles/liver-example.html)
applies it to the `dminutus_liver` data set.

## References

Ebisuzaki, Wesley. 1997. “A Method to Estimate the Statistical
Significance of a Correlation When the Data Are Serially Correlated.”
*Journal of Climate* 10 (9): 2147–53.
<https://doi.org/10.1175/1520-0442(1997)010%3C2147:AMTETS%3E2.0.CO;2>.

Jackson, Donald A. 1995. “PROTEST: A PROcrustean Randomization TEST of
Community Environment Concordance.” *Écoscience* 2 (3): 297–303.
<https://doi.org/10.1080/11956860.1995.11682297>.

Kalyuzhny, Michael. 2020. “Null Models for Community Dynamics: Beware of
the Cyclic Shift Algorithm.” *Global Ecology and Biogeography* 29 (6):
1085–93. <https://doi.org/10.1111/geb.13083>.

Legendre, Pierre, Mark R. T. Dale, Marie-Josée Fortin, Jessica
Gurevitch, Michael Hohn, and Donald Myers. 2002. “The Consequences of
Spatial Structure for the Design and Analysis of Ecological Field
Surveys.” *Ecography* 25 (5): 601–15.
<https://doi.org/10.1034/j.1600-0587.2002.250508.x>.

Oksanen, Jari, Gavin L. Simpson, F. Guillaume Blanchet, et al. 2026.
*Vegan: Community Ecology Package*.
<https://doi.org/10.32614/CRAN.package.vegan>.

Peres-Neto, Pedro R., and Donald A. Jackson. 2001. “How Well Do
Multivariate Data Sets Match? The Advantages of a Procrustean
Superimposition Approach over the Mantel Test.” *Oecologia* 129 (2):
169–78. <https://doi.org/10.1007/s004420100720>.
