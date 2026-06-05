# Sample Z under bipartite ordinal data (rectangular)

Rectangular analogue of
[`rZ_ord_fc`](https://netify-dev.github.io/lame/reference/rZ_ord_fc.md).
Cutpoints are data-induced (the boundary between rank `w` and rank `w+1`
is the maximum Z in rank `w` and minimum Z in rank `w+1`), matching the
existing unipartite convention. With rho ≡ 0 in bipartite, every cell at
rank `w` is drawn independently from a truncated normal centered at
`EZ`.

## Usage

``` r
rZ_ord_bip_fc(Z, EZ, Y)
```

## Arguments

- Z:

  current rectangular latent matrix (nA x nB).

- EZ:

  expected value (nA x nB).

- Y:

  observed ordinal matrix (nA x nB).

## Value

updated nA x nB latent matrix.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
