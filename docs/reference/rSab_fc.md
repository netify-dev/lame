# Gibbs update for additive effects covariance

Gibbs sampling for the covariance matrix of additive row and column
effects in the AME model. This function implements the inverse-Wishart
posterior update for the covariance matrix Sab.

## Usage

``` r
rSab_fc(a, b, Sab0=NULL, eta0=NULL, rvar=TRUE, cvar=TRUE, symmetric=FALSE)
```

## Arguments

- a:

  vector of row random effects (additive sender effects)

- b:

  vector of column random effects (additive receiver effects)

- Sab0:

  prior (inverse) scale matrix for the prior distribution. Default is
  diag(2), which provides a weakly informative prior.

- eta0:

  prior degrees of freedom for the prior distribution. Default is 4,
  which is the minimum for a proper prior with 2x2 matrix.

- rvar:

  logical: should row variance be updated? (default TRUE)

- cvar:

  logical: should column variance be updated? (default TRUE)

- symmetric:

  logical: is this a symmetric network? (default FALSE)

## Value

Updated covariance matrix Sab (2x2 matrix with variances on diagonal and
covariance off-diagonal)

## Details

The function implements different update strategies:

- Full update: When both rvar and cvar are TRUE, updates the full 2x2
  covariance matrix using an inverse-Wishart distribution

- Row variance only: When only rvar is TRUE, updates only Sab\[1,1\]

- Column variance only: When only cvar is TRUE, updates only Sab\[2,2\]

- Symmetric case: When symmetric is TRUE, enforces equal variances and
  high correlation (0.999) between row and column effects

## Author

Peter Hoff, Shahryar Minhas
