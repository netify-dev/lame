# Gibbs sampling of additive row and column effects and regression coefficient with independent replicate relational data

Simulates from the joint full conditional distribution of (a,b,beta),
assuming same additive row and column effects and regression coefficient
across replicates.

## Usage

``` r
rbeta_ab_rep_fc(Z.T,Sab,rho,X.T,s2=1)
```

## Arguments

- Z.T:

  n x n x T array, with the third dimension for replicates. Each slice
  of the array is a (latent) normal relational matrix, with
  multiplicative effects subtracted out

- Sab:

  row and column covariance

- rho:

  dyadic correlation

- X.T:

  n x n x p x T covariate array

- s2:

  dyadic variance

## Value

- beta:

  regression coefficients

- a:

  additive row effects

- b:

  additive column effects

## Author

Peter Hoff, Yanjun He
