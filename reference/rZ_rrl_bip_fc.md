# Sample Z under bipartite relative-rank list (rectangular)

Rectangular drop-in for
[`rZ_rrl_fc`](https://netify-dev.github.io/lame/reference/rZ_rrl_fc.md).
Row-wise full ordering: cells with higher rank have higher Z within the
row. No 0-threshold and no `odmax`. With rho ≡ 0 the unipartite triangle
coupling drops out.

## Usage

``` r
rZ_rrl_bip_fc(Z, EZ, Y, YL)
```

## Arguments

- Z:

  current rectangular latent matrix (nA x nB).

- EZ:

  expected value (nA x nB).

- Y:

  observed rank matrix (nA x nB).

- YL:

  per-row ranked-receiver list (nA x max_rank).

## Value

updated nA x nB latent matrix.

## Details

**Monotonicity convention (important for simulation).** The sampler
enforces *higher Y → higher Z*. When simulating data for recovery tests,
build Y via `rank(Z[i, ])`, not `rank(-Z[i, ])`; an inverted simulator
will identify the *negative* of the true beta and look like a sampler
bug.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
