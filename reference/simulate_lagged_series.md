# Simulate a pair of multivariate series with a known lag

Generates `X` as a multivariate seasonal series (sinusoids with a random
common phase and staggered offsets among variables) plus AR(1) noise.
Under the null (`true_lag = NULL`), `Y` is an independent draw from the
same generator. Otherwise `Y` at time `t` is a linear map of `X` at time
`t - true_lag` plus AR(1) noise; `X` is generated over
`n_times + true_lag` points so that every row of `Y` has a source in
`X`.

## Usage

``` r
simulate_lagged_series(
  n_times = 24,
  n_vars_x = 5,
  n_vars_y = 5,
  true_lag = NULL,
  period = 12,
  amplitude = 1,
  autocorr = 0.5,
  noise_sd = 1,
  seed = NULL
)
```

## Arguments

- n_times:

  Number of time points (at least 3).

- n_vars_x, n_vars_y:

  Number of variables in `X` and `Y`.

- true_lag:

  `NULL` for independent series, or a non-negative integer lag by which
  `X` leads `Y`. `0` means a contemporaneous association, not the
  absence of one.

- period:

  Length of the seasonal cycle, in time steps.

- amplitude:

  Amplitude of the seasonal sinusoids.

- autocorr:

  AR(1) coefficient of the noise, in `(-1, 1)`.

- noise_sd:

  Marginal standard deviation of the AR(1) noise.

- seed:

  Optional integer seed. The global random number stream is left
  untouched.

## Value

A list with matrices `X` (`n_times` x `n_vars_x`) and `Y` (`n_times` x
`n_vars_y`), and `true_lag`.

## Details

Setting `amplitude = 0` removes the seasonal component and leaves pure
AR(1) series; setting `autocorr = 0` gives white noise around the
seasonal signal (the design of the original simulation, Appendix S1).

## Examples

``` r
null_pair <- simulate_lagged_series(n_times = 12, seed = 42)
lagged_pair <- simulate_lagged_series(n_times = 24, true_lag = 3, seed = 42)
str(lagged_pair)
#> List of 3
#>  $ X       : num [1:24, 1:5] 1.61 1.037 0.251 1.715 -0.572 ...
#>  $ Y       : num [1:24, 1:5] 2.309 0.727 0.574 1.5 2.748 ...
#>  $ true_lag: num 3
```
