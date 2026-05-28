# Sample Z under bipartite Poisson (rectangular MH step)

Rectangular drop-in for
[`rZ_pois_fc`](https://netify-dev.github.io/lame/reference/rZ_pois_fc.md).
Each `(i,j)` cell is updated independently with a Metropolis-Hastings
step on the Poisson log-link; there is no upper/lower-triangle coupling
and no diagonal because the bipartite Y is `nA x nB` with disjoint row
and column actor sets.

## Usage

``` r
rZ_pois_bip_fc(Z, EZ, s2, Y, log_exposure = 0)
```

## Arguments

- Z:

  current rectangular latent matrix (nA x nB).

- EZ:

  expected value (regression + additive + multiplicative effects), nA x
  nB.

- s2:

  dyadic variance (proposal scale on the latent log-mean).

- Y:

  observed count matrix (nA x nB); `NA` entries are resampled from the
  prior.

- log_exposure:

  optional per-period log-exposure offset (scalar or matrix conformable
  with Y). Default `0`.

## Value

updated `nA x nB` latent matrix.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
