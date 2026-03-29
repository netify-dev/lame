# Summary of a LAME object

Provides a comprehensive summary of a fitted LAME (Longitudinal Additive
and Multiplicative Effects) model, including parameter estimates,
standard errors, credible intervals, and model diagnostics.

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

  Posterior means, posterior standard deviations, z-values, approximate
  p-values, and 95% credible intervals for dyadic, sender, and receiver
  covariates. Note: the z-values are computed as posterior mean /
  posterior SD, and the p-values are derived from a normal
  approximation. These are convenient screening statistics but are not
  formal frequentist test statistics. For rigorous inference, use the
  credible intervals or examine the full posterior via the BETA matrix
  directly.

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
