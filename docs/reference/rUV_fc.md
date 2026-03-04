# Gibbs sampling of U and V

A Gibbs sampler for updating the multiplicative effect matrices U and V

## Usage

``` r
rUV_fc(Z, U, V, Suv, rho, s2 = 1, offset = 0)
```

## Arguments

- Z:

  n X n normal relational matrix

- U:

  current value of U

- V:

  current value of V

- Suv:

  covariance of (U V)

- rho:

  dyadic correlation

- s2:

  dyadic variance

- offset:

  a matrix of the same dimension as Z. It is assumed that Z-offset is
  equal to the multiplicative effects plus dyadic noise, so the offset
  should contain any additive effects (such as
  `Xbeta(X,beta+ outer(a,b,"+") ` )

## Value

- U:

  a new value of U

- V:

  a new value of V

## Author

Peter Hoff
