# Sample Z given explicit cutpoints alpha (asymmetric ordinal)

Rectangular / square truncated-normal sample on observed cells. This is
the explicit-cutpoint version of
[`rZ_ord_fc`](https://netify-dev.github.io/lame/reference/rZ_ord_fc.md),
parameterised by the explicit `alpha` vector rather than deriving
cutpoints from Z order statistics.

## Usage

``` r
rZ_ord_explicit_fc(Z, EZ, Y_int, alpha)
```

## Arguments

- Z:

  current latent matrix

- EZ:

  expected value of Z (linear predictor)

- Y_int:

  integer-recoded ordinal Y (1..K)

- alpha:

  length-(K-1) cutpoints with `alpha[1] = 0`

## Value

updated Z
