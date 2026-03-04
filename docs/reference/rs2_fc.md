# Gibbs update for dyadic variance

Gibbs update for dyadic variance

## Usage

``` r
rs2_fc(Z, rho,offset=0,nu0=NULL,s20=NULL)
```

## Arguments

- Z:

  n X n normal relational matrix

- rho:

  current value of rho

- offset:

  matrix of the same dimension as Z. It is assumed that Z-offset is
  equal to dyadic noise, so the offset should contain any additive and
  multiplicative effects (such as
  `Xbeta(X,beta+ U%*%t(V) + outer(a,b,"+") ` )

- nu0:

  prior degrees of freedom

- s20:

  prior estimate of s2

## Value

a new value of s2

## Author

Peter Hoff
