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
smallest *m*² is reported together with the full lag profile.

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

fit <- lagged_protest(sim$X, sim$Y, lags = 0:4)
fit
#> 
#> Lagged PROTEST
#> 
#> Call: lagged_protest(X = sim$X, Y = sim$Y, lags = 0:4)
#> Permutation: cyclic shift with mirror; n = 12 time points
#> P-value adjustment across lags: holm
#> 
#> Optimal lag: 2  (m2 = 0.0922, r = 0.953, adjusted p = 0.217)
#> 
#> Lag profile:
#>  lag  n     m2     r p.value p.adjusted nperm
#>    0 12 0.3103 0.830  0.9565      1.000    22
#>    1 11 0.2078 0.890  0.0870      0.348    22
#>    2 10 0.0922 0.953  0.0435      0.217    22
#>    3  9 0.2748 0.852  0.7391      1.000    22
#>    4  8 0.3046 0.834  0.9130      1.000    22
#> 
#> Smallest attainable raw p-value: 0.0435
```

The profile lists, for every lag, the number of rows used (`n`), *m*²,
the Procrustes correlation *r* = √(1 − *m*²), the raw permutation
p-value and the p-value adjusted across lags (Holm by default; see
`p.adjust.method`). The optimal lag here is the true one, 2.

[`summary()`](https://rdrr.io/r/base/summary.html) returns the profile
as a data frame, with the optimal lag flagged:

``` r

summary(fit)
#>   lag  n         m2         r    p.value p.adjusted nperm optimal
#> 1   0 12 0.31032769 0.8304651 0.95652174  1.0000000    22   FALSE
#> 2   1 11 0.20776926 0.8900734 0.08695652  0.3478261    22   FALSE
#> 3   2 10 0.09221418 0.9527779 0.04347826  0.2173913    22    TRUE
#> 4   3  9 0.27479436 0.8515901 0.73913043  1.0000000    22   FALSE
#> 5   4  8 0.30456477 0.8339276 0.91304348  1.0000000    22   FALSE
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
points mark lags that remain significant after adjustment, and the
dotted vertical line marks the optimal lag.

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
therefore defaults to `permutation = "cyclic"`. The full `Y` series is
cyclically shifted by 1, …, n − 1 positions, and also time-reversed when
`mirror = TRUE`, *before* the lag is applied. Each permuted series keeps
its autocorrelation. All shifts are enumerated, so the test is exact and
needs no random seed. `permutation = "free"` reproduces the behaviour of
[`vegan::protest()`](https://vegandevs.github.io/vegan/reference/procrustes.html)
(Oksanen et al. 2026).

The difference matters. Under independence, a small simulation with 12
monthly samples gives:

``` r

reject <- function(permutation, n_rep) {
  rates <- rowMeans(vapply(seq_len(n_rep), function(i) {
    s <- simulate_lagged_series(n_times = 12, autocorr = 0, noise_sd = 0.3,
                                seed = i)
    fit <- lagged_protest(s$X, s$Y, lags = 0:3, permutation = permutation,
                          nperm = 99, seed = i)
    fit$profile$p.value < 0.05
  }, logical(4)))
  setNames(rates, paste("lag", 0:3))
}
rbind(cyclic = reject("cyclic", 100), free = reject("free", 50))
#>        lag 0 lag 1 lag 2 lag 3
#> cyclic  0.06  0.06  0.02  0.04
#> free    1.00  1.00  1.00  1.00
```

Rejection rates stay near the nominal 0.05 with cyclic shifts and
approach 1 with free permutation.

## Caveats

- **Complete cycles.** Cyclic shifting treats the series as one closed
  cycle. This is appropriate when the series spans whole cycles, such as
  12 monthly samples of an annual cycle. Cyclic-shift null models can be
  biased when this does not hold (Kalyuzhny 2020); for series whose
  period differs from their length, a phase-randomisation surrogate
  (Ebisuzaki 1997) is a better reference.
- **Few permutations.** With `n` time points there are only `2(n − 1)`
  cyclic permutations, so the smallest attainable p-value is
  `1 / (2(n − 1) + 1)`, 0.043 for `n = 12`. Corrections across many lags
  may then be unattainable at α = 0.05; report raw and adjusted p-values
  together with the lag profile.
- **Different window sizes.** *m*² at lag *k* is based on `n − k` rows,
  so long lags rest on fewer observations. Keep the candidate lags
  within a range justified by the biology of the system.
- **Scaling.** With `scale = TRUE` (default) variables are standardised
  on the full series before truncation.

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
