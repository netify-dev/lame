# Compute Xbeta product for bipartite networks

Computes sum_k beta_k \* X_k for a single time slice

## Usage

``` r
Xbeta_bip_cpp(X, beta)
```

## Arguments

- X:

  3D array (nA x nB x p) of covariates for one time period

- beta:

  Coefficient vector of length p

## Value

nA x nB matrix
