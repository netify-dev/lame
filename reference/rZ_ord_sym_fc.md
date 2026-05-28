# Sample Z under symmetric ordinal data (square symmetric matrix)

Symmetric-Z drop-in for
[`rZ_ord_fc`](https://netify-dev.github.io/lame/reference/rZ_ord_fc.md):
each unordered dyad i, j carries a single latent value (Z_ij = Z_ji)
drawn from a truncated normal with mean equal to the average linear
predictor (EZ_ij + EZ_ji) / 2 and variance 1/2 (precision 2). The lower
triangle is mirrored after the upper-triangle draws so the full matrix
is exactly symmetric at every sweep.

## Usage

``` r
rZ_ord_sym_fc(Z, EZ, Y)
```

## Arguments

- Z:

  current symmetric latent matrix (n x n, Z = t(Z)).

- EZ:

  expected value of Z (regression + additive + multiplicative
  contributions); can be asymmetric, the helper averages the two
  directed contributions.

- Y:

  observed symmetric ordinal matrix (n x n, Y = t(Y)); diagonal is
  ignored.

## Value

updated symmetric n x n latent matrix with `NA` on the diagonal.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
