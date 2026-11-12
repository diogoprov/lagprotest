# Patos Lagoon juvenile fish and environment blocks (PELD-ELPA)

Builds two monthly blocks, environment and fish, for one habitat of the
Patos Lagoon estuary, southern Brazil, from the PELD-ELPA beach-seine
survey (Darwin Core Archive v1.16). The data are licensed CC BY-NC 4.0
and are not shipped with lagprotest: pass a local copy of the archive in
`path`, or leave `path = NULL` to download it from the SiBBr IPT into
the session's temporary directory.

## Usage

``` r
elpa_blocks(
  habitat = c("estuary", "sandy beach"),
  from = c(1996, 8),
  to = c(2016, 12),
  species = NULL,
  hellinger = TRUE,
  path = NULL
)
```

## Source

Vieira, J. P., Garcia, A. M. and Lemos, V. M. PELD-ELPA: Species
composition and abundance patterns of fish assemblages at shallow waters
of Patos Lagoon estuary. Version 1.16. Sistema de Informação sobre a
Biodiversidade Brasileira (SiBBr), via GBIF,
[doi:10.15468/kci8zb](https://doi.org/10.15468/kci8zb) . Licensed CC
BY-NC 4.0: cite the dataset and do not use it commercially. The program
is described in Lemos et al. (2022), Earth System Science Data 14:
1015–1041,
[doi:10.5194/essd-14-1015-2022](https://doi.org/10.5194/essd-14-1015-2022)
.

## Arguments

- habitat:

  `"estuary"` or `"sandy beach"`.

- from, to:

  Integer vectors `c(year, month)` delimiting the period. The default
  covers 1996-08 to 2016-12, the longest stretch without unsampled
  months.

- species:

  Species (scientific names) to keep in the fish block. `NULL` (default)
  keeps the most abundant species that together make up 90% of the total
  catch over all sampling events in the chosen habitat and period;
  `"all"` keeps every species.

- hellinger:

  Logical. If `TRUE` (default), the fish block holds
  Hellinger-transformed catches, computed on all species before
  `species` are selected; if `FALSE`, mean catch per haul.

- path:

  Path to the archive (`.zip`) or to a directory with its extracted
  files. `NULL` downloads version 1.16 from the SiBBr IPT (about 3 MB),
  once per session.

## Value

A named list of two numeric matrices with one row per month (row names
`"YYYY-MM"`):

- `Environment`:

  Water temperature (`temp`, °C), salinity (`sal`) and Secchi depth
  (`secchi`).

- `Fish`:

  One column per species in `species`.

## Details

Each sampling event is one site in one month, with five seine hauls.
Counts are summed by species within an event and divided by five (catch
per haul), and both catch and environment are averaged over the sites of
the chosen habitat in each month. Months with no record of an
environmental variable are filled by linear interpolation; a month with
no fish sampling at all is an error, since the blocks must be regular
series. No transformation, detrending or deseasonalization is applied to
the environment; see the example for the steps used in the lagprotest
preprint, which also added freshwater discharge from the Brazilian
National Water Agency (not included here).

## Examples

``` r
if (FALSE) { # interactive()
el <- elpa_blocks("estuary")
detrend <- function(m) apply(m, 2, function(x) residuals(lm(x ~ seq_along(x))))
X <- detrend(log1p(el$Environment))
Y <- detrend(el$Fish)
fit <- lagged_protest(X, Y, lags = 0:6, period = 12, nperm = 999,
                      seed = 1)
fit
}
```
