# Summary of an AME object

Provides a comprehensive summary of a fitted AME (Additive and
Multiplicative Effects) model, including parameter estimates, standard
errors, credible intervals, and model diagnostics.

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
