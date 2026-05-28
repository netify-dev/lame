# Cowles MH update for explicit ordinal cutpoints (Z-marginalised)

One Metropolis-Hastings sweep on `alpha` with Z integrated out of the
conditional. The proposal is a symmetric random walk in
`delta = log(diff(alpha))` space; the acceptance ratio uses the marginal
ordinal-probit likelihood. Mixing of `alpha` under this update is
dramatically faster than the data-induced Gibbs convention (Cowles 1996
reports ESS gains of an order of magnitude on standard ordinal-probit
benchmarks).

## Usage

``` r
sample_alpha_cowles(alpha, Y_int, EZ, tau_prop, symmetric = FALSE)
```

## Arguments

- alpha:

  numeric, current cutpoints with `alpha[1] = 0`, length `K-1`.

- Y_int:

  integer matrix (or array) of ordinal categories coded 1..K, with `NA`
  for missing cells.

- EZ:

  linear predictor at the dyad level, same shape as `Y_int`.

- tau_prop:

  numeric, RW proposal SD in delta space (will be adapted in burn-in by
  the caller).

- symmetric:

  logical, whether the network is symmetric (uses upper-triangle cells
  and sqrt(2) precision scaling).

## Value

list with elements `alpha`, `delta`, `accept` (logical).

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
