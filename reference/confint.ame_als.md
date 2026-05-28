# Confidence intervals for a fast AME fit

Wald confidence intervals for the intercept and dyadic-covariate
coefficients of an
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) fit,
built from the conditional sandwich covariance
([`vcov.ame_als`](https://netify-dev.github.io/lame/reference/vcov.ame_als.md)).

## Usage

``` r
# S3 method for class 'ame_als'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  an `ame_als` fit.

- parm:

  character vector of parameter names, or integer indices; if `NULL`
  (default) all available coefficients are returned.

- level:

  confidence level (default 0.95).

- ...:

  passed to
  [`vcov.ame_als`](https://netify-dev.github.io/lame/reference/vcov.ame_als.md)
  (e.g. `cluster`).

## Value

A matrix with one row per coefficient and lower/upper bound columns.

## Details

These intervals are a fast convenience. They are **conditional** (the
additive and multiplicative effects are held fixed) and therefore
anti-conservative, and they cover only the regression coefficients the
sandwich covariance is defined for – not the node-covariate, additive or
multiplicative parameters. For fully-propagated intervals on all
parameters, use
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
and
[`confint.boot_ame`](https://netify-dev.github.io/lame/reference/confint.boot_ame.md).

## See also

[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
for bootstrap intervals on all parameters.
