# Summary of an AME object

Provides a comprehensive summary of a fitted AME (Additive and
Multiplicative Effects) model, including parameter estimates, standard
errors, confidence intervals, and model diagnostics.

## Usage

``` r
# S3 method for class 'ame'
summary(object, ...)
```

## Arguments

- object:

  an object of class "ame", typically the result of fitting an AME model
  using the `ame` function

- ...:

  additional parameters (currently not used)

## Value

A list of class "summary.ame" containing:

- call:

  The original function call

- beta:

  Matrix of regression coefficient estimates and statistics

- variance:

  Matrix of variance component estimates

## Details

The summary includes:

- Regression coefficients:

  Point estimates, standard errors, z-values, p-values, and 95%
  confidence intervals for dyadic, sender, and receiver covariates

- Variance components:

  Estimates and standard errors for:

  va

  :   Variance of additive sender/row effects (asymmetric networks)

  cab

  :   Covariance between sender and receiver effects

  vb

  :   Variance of additive receiver/column effects (asymmetric networks)

  rho

  :   Dyadic correlation (reciprocity in directed networks)

  ve

  :   Residual variance

  For symmetric networks, only va and ve are estimated.

## See also

[`ame`](https://netify-dev.github.io/lame/reference/ame.md),
[`print.summary.ame`](https://netify-dev.github.io/lame/reference/print.summary.ame.md)

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
