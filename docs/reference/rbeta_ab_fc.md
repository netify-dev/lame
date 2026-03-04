# Conditional simulation of additive effects and regression coefficients

Simulates from the joint full conditional distribution of (beta,a,b) in
a social relations regression model

## Usage

``` r
rbeta_ab_fc(
  Z,
  Sab,
  rho,
  X = NULL,
  s2 = 1,
  offset = 0,
  iV0 = NULL,
  m0 = NULL,
  g = length(Z)
)
```

## Arguments

- Z:

  n X n normal relational matrix

- Sab:

  row and column covariance

- rho:

  dyadic correlation

- X:

  n x n x p covariate array

- s2:

  dyadic variance

- offset:

  a matrix of the same dimension as Z. It is assumed that Z-offset
  follows a SRRM, so the offset should contain any multiplicative
  effects (such as `U%*% t(V) ` )

- iV0:

  prior precision matrix for regression parameters

- m0:

  prior mean vector for regression parameters

- g:

  prior variance scale for g-prior when iV0 is unspecified

## Value

- beta:

  regression coefficients

- a:

  additive row effects

- b:

  additive column effects

## Author

Peter Hoff
