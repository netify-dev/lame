# Compute X'X and X'y for bipartite covariate regression

Replaces the O(T x p^2 x n^2) nested R loop for bipartite XtX and Xty
computation with a single C++ call.

## Usage

``` r
compute_XtX_Xty_bip_cpp(Xlist, resid, p)
```

## Arguments

- Xlist:

  List of T arrays, each nA x nB x p

- resid:

  3D array of residuals nA x nB x T

- p:

  Number of covariates

## Value

List with XtX (p x p) and Xty (p vector)
