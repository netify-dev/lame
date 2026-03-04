# Gibbs sampling of U and V

A Gibbs sampler for updating the multiplicative effect matrices U and V,
assuming they are the same across replicates.

## Usage

``` r
rUV_rep_fc(E.T,U,V,rho,s2=1,shrink=TRUE)
```

## Arguments

- E.T:

  Array of square residual relational matrix series with additive
  effects and covariates subtracted out. The third dimension of the
  array is for different replicates. Each slice of the array according
  to the third dimension is a square residual relational matrix.

- U:

  current value of U

- V:

  current value of V

- rho:

  dyadic correlation

- s2:

  dyadic variance

- shrink:

  adaptively shrink the factors with a hierarchical prior

## Value

- U:

  a new value of U

- V:

  a new value of V

## Author

Peter Hoff, Yanjun He
