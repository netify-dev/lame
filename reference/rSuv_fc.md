# Gibbs update for multiplicative effects covariance

Gibbs sampling for the covariance matrix of multiplicative effects U and
V in the AME model. This function implements the inverse-Wishart
posterior update for the covariance matrix of the stacked UV effects.

## Usage

``` r
rSuv_fc(U, V, Suv0=NULL, kappa0=NULL)
```

## Arguments

- U:

  matrix of multiplicative row effects (n x R matrix where n is the
  number of nodes and R is the dimension of multiplicative effects)

- V:

  matrix of multiplicative column effects (n x R matrix)

- Suv0:

  prior (inverse) scale matrix for the prior distribution. Default is
  identity matrix of dimension 2R x 2R, providing a weakly informative
  prior.

- kappa0:

  prior degrees of freedom for the prior distribution. Default is 2 +
  2R, which is the minimum for a proper prior with a 2R x 2R covariance
  matrix.

## Value

Updated covariance matrix Suv (2R x 2R matrix) for the stacked effects
\[U, V\]. The first R x R block contains covariances for U, the last R x
R block contains covariances for V, and the off-diagonal blocks contain
cross-covariances between U and V.

## Details

The function updates the full covariance matrix for multiplicative
effects using an inverse-Wishart distribution. The posterior
distribution is: \$\$Suv ~ IW(kappa0 \* Suv0 + t(UV) \\\*\\ UV, n +
kappa0)\$\$ where UV = cbind(U, V) is the stacked matrix of effects.

This hierarchical prior allows for adaptive shrinkage of the
multiplicative effects, with the amount of shrinkage determined by the
data through the posterior update.

## Author

Peter Hoff, Shahryar Minhas
