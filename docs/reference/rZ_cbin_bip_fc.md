# Sample Z under bipartite censored binary nominations (rectangular)

Rectangular drop-in for
[`rZ_cbin_fc`](https://netify-dev.github.io/lame/reference/rZ_cbin_fc.md).
Row-wise constraints: `Y_ij = 1 => Z_ij > 0`;
`Y_ij = 0 & odobs_i < odmax_i => Z_ij < 0`; nominated alters dominate
non-nominated alters within the row. With rho ≡ 0 there is no
reciprocity term.

## Usage

``` r
rZ_cbin_bip_fc(Z, EZ, Y, odmax, odobs)
```

## Arguments

- Z:

  current rectangular latent matrix (nA x nB).

- EZ:

  expected value (nA x nB).

- Y:

  observed 0/1 matrix (nA x nB); `NA` entries resampled from the
  unconstrained prior.

- odmax:

  row-wise maximum nominations (length nA or scalar).

- odobs:

  row-wise observed outdegree (length nA).

## Value

updated nA x nB latent matrix.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
