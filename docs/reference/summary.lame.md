# Summary of a LAME object

Provides a comprehensive summary of a fitted LAME (Longitudinal Additive
and Multiplicative Effects) model, including parameter estimates,
standard errors, confidence intervals, and model diagnostics.

## Usage

``` r
# S3 method for class 'lame'
summary(object, ...)
```

## Arguments

- object:

  an object of class "lame", typically the result of fitting a
  longitudinal AME model using the `lame` function

- ...:

  additional parameters (currently not used)

## Value

A list of class "summary.lame" containing:

- call:

  The original function call

- beta:

  Matrix of regression coefficient estimates and statistics

- variance:

  Matrix of variance component estimates

- n.periods:

  Number of time periods in the longitudinal data

## Details

The summary includes:

- Regression coefficients:

  Point estimates, standard errors, z-values, p-values, and 95%
  confidence intervals for dyadic, sender, and receiver covariates

- Variance components:

  Estimates and standard errors for:

  va

  :   Variance of additive sender/row effects

  cab

  :   Covariance between sender and receiver effects

  vb

  :   Variance of additive receiver/column effects

  rho

  :   Dyadic correlation (reciprocity)

  ve

  :   Residual variance

## See also

[`lame`](https://netify-dev.github.io/lame/reference/lame.md),
[`print.summary.lame`](https://netify-dev.github.io/lame/reference/print.summary.lame.md)

## Author

Peter Hoff, Cassy Dorff, Shahryar Minhas
