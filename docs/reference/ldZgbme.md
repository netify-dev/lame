# log density for GBME models

Calculation of the log conditional density of the latent AMEN matrix `Z`
given observed data `Y`.

## Usage

``` r
ldZgbme(Z, Y, llYZ, EZ, rho, s2 = 1)
```

## Arguments

- Z:

  n X n latent relational matrix following an AMEN model

- Y:

  n X n observed relational matrix

- llYZ:

  a vectorizable function taking two arguments, y and z. See details
  below.

- EZ:

  n X n mean matrix for `Z` based on AMEN model (including additive
  effects)

- rho:

  dyadic correlation in AMEN model for `Z`

- s2:

  residual variance in AMEN model for `Z`

## Value

a symmetric matrix where entry i,j is proportional to the log
conditional bivariate density of `z[i,j],z[j,i]`.

## Details

This function is used for updating dyadic pairs of the latent variable
matrix `Z` based on `Y` and an AMEN model for `Z`. The function `llYZ`
specifies the log likelihood for each single `z[i,j]` based on `y[i,j]`,
that is, `llYZ` gives the log probability density (or mass function) of
`y[i,j]` given `z[i,j]`.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
## For (overdispersed) Poisson regression, use
llYZ<-function(y,z){ dpois(y,z,log=TRUE) }
```
