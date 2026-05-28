# Log-likelihood is not defined for a fast AME fit

The fast estimator minimises a (working-response) least-squares
objective, not a family likelihood, so it has no log-likelihood and
`AIC` / `BIC` are undefined. Calling
[`logLik()`](https://rdrr.io/r/stats/logLik.html) raises an informative
error rather than returning a misleading number.

## Usage

``` r
# S3 method for class 'ame_als'
logLik(object, ...)
```

## Arguments

- object:

  an `ame_als` object.

- ...:

  ignored.

## Value

Never returns; raises an error.
