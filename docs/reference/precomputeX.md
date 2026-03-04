# Precompute design matrix statistics

Precomputes summary statistics of the design array X to speed up MCMC
calculations in the AME model

## Usage

``` r
precomputeX(X)
```

## Arguments

- X:

  n x n x p covariate array

## Value

The input array with precomputed statistics as attributes:

- Xr:

  row sums

- Xc:

  column sums

- mX:

  matricized version

- mXt:

  dyad-transposed matricized version

- XX:

  regression sums of squares

- XXt:

  crossproduct sums of squares

## Author

Peter Hoff
