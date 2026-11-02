# Example: seasonal liver plasticity in a treefrog

``` r

library(lagprotest)
```

This article analyzes a real data set with
[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
and
[`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md).
The question is whether seasonal changes in the liver of the treefrog
*Dendropsophus minutus* track climate within the same month or with a
delay, and whether hepatic modules lead one another.

## The data

`dminutus_liver` holds monthly means of four hepatic modules and of
climate, from April 2016 to March 2017 \[Franco-Belussi et al. (2026);
see
[`?dminutus_liver`](https://provetelab.org/lagprotest/reference/dminutus_liver.md)
for the variables\]. Each element is a matrix with 12 rows (months) and
one column per variable. Climate combines temperature, rainfall and
humidity from BR-DWGD (Xavier et al. 2022), a drought index from
TerraClimate (Abatzoglou et al. 2018), UV-B irradiance from glUV
(Beckmann et al. 2014) and photoperiod.

``` r

sapply(dminutus_liver, dim)
#>      Somatic Histochemistry Morphometry Density Climate
#> [1,]      12             12          12      12      12
#> [2,]       4              4           5       4       7
```

The series span exactly one annual cycle, so the phase and cyclic
schemes, which treat the series as periodic with period *n*, apply
without truncation (`period` is not needed). With 12 time points, the
power to detect a lag is limited, and strong autocorrelation would make
the phase scheme somewhat liberal (see [Choosing a permutation
scheme](https://provetelab.org/lagprotest/articles/permutation-schemes.md));
the cyclic scheme is used below as a conservative check.

## Contemporaneous concordance

Concordance within the same month is the lag-0 PROTEST for every pair of
data sets. With a single lag, the global p-value equals the raw one, and
[`lagged_protest_pairs()`](https://provetelab.org/lagprotest/reference/lagged_protest_pairs.md)
adjusts it across the ten pairs (Benjamini–Hochberg by default):

``` r

lag0 <- lagged_protest_pairs(dminutus_liver, lags = 0, bidirectional = FALSE,
                             nperm = 9999, seed = 1)
knitr::kable(lag0[, c("x", "y", "m2", "r", "p.value", "p.global.adjusted")],
             digits = 3)
```

| x              | y              |    m2 |     r | p.value | p.global.adjusted |
|:---------------|:---------------|------:|------:|--------:|------------------:|
| Somatic        | Histochemistry | 0.711 | 0.538 |   0.092 |             0.254 |
| Somatic        | Morphometry    | 0.815 | 0.430 |   0.374 |             0.623 |
| Somatic        | Density        | 0.609 | 0.625 |   0.102 |             0.254 |
| Somatic        | Climate        | 0.477 | 0.724 |   0.004 |             0.045 |
| Histochemistry | Morphometry    | 0.806 | 0.440 |   0.609 |             0.751 |
| Histochemistry | Density        | 0.704 | 0.544 |   0.190 |             0.381 |
| Histochemistry | Climate        | 0.791 | 0.457 |   0.743 |             0.751 |
| Morphometry    | Density        | 0.755 | 0.495 |   0.558 |             0.751 |
| Morphometry    | Climate        | 0.772 | 0.478 |   0.751 |             0.751 |
| Density        | Climate        | 0.378 | 0.789 |   0.099 |             0.254 |

Somatic indices and climate were the only pair concordant after the
adjustment. Density and climate had the lowest *m*², yet this
concordance was not distinguishable from that expected between
independent seasonal series, which illustrates why a low *m*² alone is
not evidence of association.

## Lagged concordance

The hypothesized hierarchy has seven directed edges: climate to each of
the four modules, and three links between modules. Each edge is scanned
from −3 to +3 months, so a positive lag means that the first data set
leads the second. The global p-value accounts for the scan within each
edge; the edges are then adjusted within two families of hypotheses,
because the four climate edges ask whether each module responds to the
same driver, whereas the biological edges ask a different question:

``` r

climate_edges <- list(c("Climate", "Histochemistry"),
                      c("Climate", "Morphometry"),
                      c("Climate", "Density"),
                      c("Climate", "Somatic"))
module_edges <- list(c("Histochemistry", "Morphometry"),
                     c("Morphometry", "Density"),
                     c("Density", "Somatic"))
climate <- lagged_protest_pairs(dminutus_liver, pairs = climate_edges,
                                lags = 0:3, nperm = 9999, seed = 1)
modules <- lagged_protest_pairs(dminutus_liver, pairs = module_edges,
                                lags = 0:3, nperm = 9999, seed = 1)
cols <- c("x", "y", "optimal_lag", "n", "m2", "r", "p.value", "p.global",
          "p.global.adjusted")
knitr::kable(climate[, cols], digits = 3)
```

| x | y | optimal_lag | n | m2 | r | p.value | p.global | p.global.adjusted |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| Climate | Histochemistry | 3 | 9 | 0.379 | 0.788 | 0.011 | 0.011 | 0.044 |
| Climate | Morphometry | -3 | 9 | 0.422 | 0.761 | 0.062 | 0.137 | 0.274 |
| Climate | Density | 0 | 12 | 0.378 | 0.789 | 0.102 | 0.828 | 0.828 |
| Climate | Somatic | 0 | 12 | 0.477 | 0.724 | 0.005 | 0.225 | 0.300 |

``` r

knitr::kable(modules[, cols], digits = 3)
```

| x | y | optimal_lag | n | m2 | r | p.value | p.global | p.global.adjusted |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| Histochemistry | Morphometry | 1 | 11 | 0.664 | 0.580 | 0.148 | 0.777 | 0.777 |
| Morphometry | Density | -1 | 11 | 0.515 | 0.696 | 0.036 | 0.459 | 0.688 |
| Density | Somatic | 1 | 11 | 0.436 | 0.751 | 0.006 | 0.100 | 0.301 |

Climate led histochemistry by three months, and this was the only edge
whose global p-value remained below 0.05 after adjustment within its
family. Notice the difference between the raw p-value at the optimal lag
and the global p-value: for climate and somatic indices, the raw p-value
at lag 0 was below 0.01, but the global p-value, which accounts for
having chosen the best of seven lags, was not significant. The
contemporaneous test above, which examines a single lag, is the
appropriate one for that pair.

The lag profile of the climate–histochemistry edge shows the minimum at
lag +3:

``` r

fit <- attr(climate, "fits")[["Climate -> Histochemistry"]]
plot(fit)
```

![Procrustes m2 against lag from -3 to 3 for climate and histochemistry;
the minimum, circled and filled, is at lag 3, with global p =
0.011.](liver-example_files/figure-html/profile-1.png)

``` r

knitr::kable(fit$profile, digits = 3)
```

| lag |   n |    m2 |     r | p.value | p.adjusted | nperm |
|----:|----:|------:|------:|--------:|-----------:|------:|
|  -3 |   9 | 0.717 | 0.532 |   0.496 |      0.496 |  9999 |
|  -2 |  10 | 0.801 | 0.446 |   0.841 |      0.841 |  9999 |
|  -1 |  11 | 0.678 | 0.567 |   0.239 |      0.239 |  9999 |
|   0 |  12 | 0.791 | 0.457 |   0.738 |      0.738 |  9999 |
|   1 |  11 | 0.718 | 0.531 |   0.462 |      0.462 |  9999 |
|   2 |  10 | 0.593 | 0.638 |   0.134 |      0.134 |  9999 |
|   3 |   9 | 0.379 | 0.788 |   0.011 |      0.011 |  9999 |

The cyclic scheme gives a conservative check. With 12 rows and
`mirror = TRUE`, its smallest attainable p-value is 1/23 = 0.043, so a
global p-value at that bound is as strong as the scheme can show:

``` r

cyc <- lagged_protest(dminutus_liver$Climate, dminutus_liver$Histochemistry,
                      lags = 0:3, bidirectional = TRUE,
                      permutation = "cyclic")
cyc$p.global
#> [1] 0.04347826
```

## Reading the results

The analysis supports a delayed response of liver histochemistry to
seasonal climate, of about three months, and a contemporaneous
association between climate and somatic condition. The other edges were
not distinguishable from chance with 12 monthly samples, which does not
show that they are absent: one annual cycle gives little power for weak
or delayed associations, and the optimal lag of a non-significant edge
should not be interpreted. Because the lag is estimated from a single
year, it also describes this cycle rather than a general property of the
species.

## References

Abatzoglou, John T., Solomon Z. Dobrowski, Sean A. Parks, and Katherine
C. Hegewisch. 2018. “TerraClimate, a High-Resolution Global Dataset of
Monthly Climate and Climatic Water Balance from 1958–2015.” *Scientific
Data* 5: 170191. <https://doi.org/10.1038/sdata.2017.191>.

Beckmann, Michael, Tomáš Václavík, Ameur M. Manceur, et al. 2014. “glUV:
A Global UV-B Radiation Data Set for Macroecological Studies.” *Methods
in Ecology and Evolution* 5: 372–83.
<https://doi.org/10.1111/2041-210X.12168>.

Franco-Belussi, Lilian, Luciana Moraes, Classius de Oliveira, Carlos
Eurico Fernandes, and Diogo Borges Provete. 2026. *Data and Code for:
Seasonal Hepatic Plasticity Follows a Temporal Response Hierarchy in a
Neotropical Frog*. Zenodo. <https://doi.org/10.5281/zenodo.21613832>.

Xavier, Alexandre C., Bridget R. Scanlon, Carey W. King, and Aline I.
Alves. 2022. “New Improved Brazilian Daily Weather Gridded Data
(1961–2020).” *International Journal of Climatology* 42: 8390–404.
<https://doi.org/10.1002/joc.7731>.
