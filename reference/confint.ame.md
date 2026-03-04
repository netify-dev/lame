# Bayesian credible intervals for AME model coefficients

Returns posterior quantile-based credible intervals for regression
coefficients. These are Bayesian credible intervals, not frequentist
confidence intervals.

## Usage

``` r
# S3 method for class 'ame'
confint(object, parm = NULL, level = 0.95, ...)

# S3 method for class 'lame'
confint(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  fitted AME model (class "ame")

- parm:

  character vector of parameter names, or numeric indices. If NULL
  (default), intervals for all parameters are returned.

- level:

  credible level (default 0.95)

- ...:

  additional arguments (ignored)

## Value

matrix with columns for lower and upper bounds
