# Sample Z under bipartite fixed-rank nominations (rectangular)

Rectangular analogue of
[`rZ_frn_fc`](https://netify-dev.github.io/lame/reference/rZ_frn_fc.md).
Row-wise rank constraints: ranked nominees are positive and ordered by
rank; non-nominees fall below all ranked cells in the row (and below 0
when the row is unsaturated). With rho ≡ 0 the unipartite triangle
coupling drops out.

## Usage

``` r
rZ_frn_bip_fc(Z, EZ, Y, YL, odmax, odobs)
```

## Arguments

- Z:

  current rectangular latent matrix (nA x nB).

- EZ:

  expected value (nA x nB).

- Y:

  observed rank matrix (nA x nB); 0 = unranked, 1..ncol(YL) = rank
  order.

- YL:

  per-row ranked-receiver list (nA x max_rank).

- odmax:

  row-wise max nominations.

- odobs:

  row-wise observed outdegree.

## Value

updated nA x nB latent matrix.

## Details

**Monotonicity convention (important for simulation).** The sampler
enforces *higher Y → higher Z*: when `Y[i, j] > Y[i, k]`, the
corresponding latent values satisfy `Z[i, j] > Z[i, k]`. When simulating
data for recovery tests, build Y via `rank(Z[i, ])` (low rank = small
Y), not `rank(-Z[i, ])`. the inverted convention identifies the negative
of the true beta.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
