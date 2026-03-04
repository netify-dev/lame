# Griddy Gibbs update for dyadic correlation

Simulation of dyadic correlation from its approximate full conditional
distribution using griddy Gibbs sampling

## Usage

``` r
rrho_fc(Z, Sab, s2 = 1, offset = 0, ngp = 100, asp = NULL)
```

## Arguments

- Z:

  n X n normal relational matrix

- Sab:

  covariance of additive effects

- s2:

  residual variance

- offset:

  matrix of the same dimension as Z. It is assumed that Z-offset follows
  an SRM distribution, so the offset should contain any regression terms
  and multiplicative effects (such as `Xbeta(X,beta+ U%*%t(V) ` )

- ngp:

  the number of points for an unevenly-spaced grid on which to
  approximate the full conditional distribution

- asp:

  use arc sine prior (TRUE) or uniform prior (FALSE)

## Value

a value of rho

## Author

Peter Hoff
