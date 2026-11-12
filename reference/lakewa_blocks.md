# Lake Washington plankton and limnology blocks

Builds three monthly blocks (environment, phytoplankton and zooplankton)
for Lake Washington, Washington, USA, from `lakeWAplanktonRaw` in the
MARSS package. The data are not shipped with lagprotest; this function
reads them from an installed copy of MARSS, which must be installed
separately.

## Usage

``` r
lakewa_blocks(from = 1962, to = 1994, interpolate = TRUE)
```

## Source

`lakeWAplanktonRaw` in the MARSS package (Holmes et al., GPL-2),
described there as adapted from the Lake Washington database of W. T.
Edmondson, funded by the Andrew Mellon Foundation, with data courtesy of
D. Schindler, University of Washington. See
[`help("lakeWAplankton", package = "MARSS")`](https://atsa-es.github.io/MARSS/reference/plankton.html)
for the original documentation and how to cite it.

## Arguments

- from, to:

  First and last year to keep (between 1962 and 1994).

- interpolate:

  Logical. Fill missing values by linear interpolation? If `FALSE`,
  missing values are kept and must be handled before calling
  [`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md).

## Value

A named list of three numeric matrices with one row per month (row names
`"YYYY-MM"`):

- `Environment`:

  Water temperature (`Temp`, °C), total phosphorus (`TP`) and `pH`.

- `Phytoplankton`:

  `log1p` of the counts (cells per mL) of `Cryptomonas`, `Diatoms`,
  `Greens`, `Unicells` and `Other.algae`.

- `Zooplankton`:

  `log1p` of the counts (individuals per L) of `Cyclops`, `Diaptomus`,
  `Daphnia`, `Epischura`, `Conochilus`, `Non.daphnid.cladocerans` and
  `Non.colonial.rotifers`.

## Details

Rows are restricted to the years `from` to `to`. Missing values (a few
isolated months in each column) are filled by linear interpolation over
the month index when `interpolate = TRUE`; a gap at either end of the
selected period is filled with the nearest observed value. Plankton
counts are then transformed with
[`log1p()`](https://rdrr.io/r/base/Log.html). The environmental
variables are left on their original scale, since
[`lagged_protest()`](https://provetelab.org/lagprotest/reference/lagged_protest.md)
standardizes columns by default (`scale = TRUE`).

## Examples

``` r
lw <- lakewa_blocks()
str(lw)
#> List of 3
#>  $ Environment  : num [1:396, 1:3] 7.6 7.5 7.67 10.1 12.8 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:396] "1962-01" "1962-02" "1962-03" "1962-04" ...
#>   .. ..$ : chr [1:3] "Temp" "TP" "pH"
#>  $ Phytoplankton: num [1:396, 1:5] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:396] "1962-01" "1962-02" "1962-03" "1962-04" ...
#>   .. ..$ : chr [1:5] "Cryptomonas" "Diatoms" "Greens" "Unicells" ...
#>  $ Zooplankton  : num [1:396, 1:7] 3.08 2.45 1.94 2.05 1.16 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:396] "1962-01" "1962-02" "1962-03" "1962-04" ...
#>   .. ..$ : chr [1:7] "Cyclops" "Diaptomus" "Daphnia" "Epischura" ...
fit <- lagged_protest(lw$Environment, lw$Zooplankton, lags = 0:3,
                      nperm = 199, seed = 1)
fit$optimal_lag
#> [1] 0
```
