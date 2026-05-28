# Bayesian credible intervals for AME model parameters

Returns posterior **equal-tailed quantile-based** credible intervals
(not highest-posterior-density, HPD). Built directly from
`quantile(object$BETA, c(alpha/2, 1-alpha/2))` and
`quantile(object$VC, ...)`. These are Bayesian credible intervals, not
frequentist confidence intervals.

## Usage

``` r
# S3 method for class 'ame'
confint(object, parm = NULL, level = 0.95, ...)

# S3 method for class 'lame'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  fitted AME / LAME model.

- parm:

  character vector of parameter names, or numeric indices. When `parm`
  is character, both BETA names (e.g. `"intercept"`, `"x1_dyad"`) and
  variance-component names (`"va"`, `"vb"`, `"cab"`, `"rho"`, `"ve"`)
  are accepted. `NULL` (default) returns intervals for all available
  parameters.

- level:

  credible level (default 0.95).

- ...:

  additional arguments (ignored).

## Value

Matrix with one row per parameter and two columns (e.g. `"2.5%"`,
`"97.5%"`).

## Note on interval type

These are equal-tailed quantile intervals, not HPD. For an HPD interval
use e.g.
[`coda::HPDinterval`](https://rdrr.io/pkg/coda/man/HPDinterval.html) on
the columns of `object$BETA` and `object$VC` directly.
