# Symmetric ordinal Z sample given explicit cutpoints (R \>= 0 safe)

Mirror image of
[`rZ_ord_sym_fc`](https://netify-dev.github.io/lame/reference/rZ_ord_sym_fc.md)
(precision-2 / variance-1/2 on the upper triangle, mirrored to the
lower) parameterised by an explicit `alpha` vector instead of
data-induced cutpoints.

## Usage

``` r
rZ_ord_sym_explicit_fc(Z, EZ, Y_int, alpha)
```

## Arguments

- Z:

  current symmetric latent matrix

- EZ:

  expected value of Z

- Y_int:

  integer ordinal Y (symmetric)

- alpha:

  length-(K-1) cutpoints with `alpha[1] = 0`

## Value

updated symmetric Z (diagonal sampled from unconstrained prior)
