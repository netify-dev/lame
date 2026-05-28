# Confidence intervals from a fast AME bootstrap

Bootstrap confidence intervals for the intercept and regression
coefficients of a fast AME fit.

## Usage

``` r
# S3 method for class 'boot_ame'
confint(
  object,
  parm = NULL,
  level = 0.95,
  ci_type = c("basic", "percentile"),
  which = c("all", "beta", "vc", "a", "b", "U", "V"),
  ...
)
```

## Arguments

- object:

  a `boot_ame` object.

- parm:

  character vector of parameter names, or integer indices. If `NULL`
  (default), all coefficients are returned.

- level:

  confidence level (default 0.95).

- ci_type:

  interval type: `"basic"` (default) reflects the replicate quantiles
  about the point estimate (\\2\hat\theta - q\\) and so corrects
  first-order bias – essential for the IRLS binary/Poisson point
  estimators, which carry finite-sample bias; `"percentile"` uses the
  replicate quantiles directly and is best when the estimator is
  approximately unbiased.

- which:

  which uncertainty channel to return: `"all"` (default, intercept +
  regression coefficients + variance components + sender effects

  - receiver effects + multiplicative U/V), or a single channel:
    `"beta"`, `"vc"`, `"a"`, `"b"`, `"U"`, `"V"`.

- ...:

  ignored.

## Value

A matrix with one row per parameter and lower/upper bound columns.
