# Metropolis update for dyadic correlation

Metropolis update for dyadic correlation

## Usage

``` r
rrho_mh(Z, rho, s2 = 1,offset=0, asp=NULL)
```

## Arguments

- Z:

  n X n normal relational matrix

- rho:

  current value of rho

- s2:

  current value of s2

- offset:

  matrix of the same dimension as Z. It is assumed that Z-offset is
  equal to dyadic noise, so the offset should contain any additive and
  multiplicative effects (such as
  `Xbeta(X,beta+ U%*%t(V) + outer(a,b,"+") ` )

- asp:

  use arc sine prior (TRUE) or uniform prior (FALSE)

## Value

a new value of rho

## Author

Peter Hoff
