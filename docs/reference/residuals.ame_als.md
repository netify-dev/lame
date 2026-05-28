# Residuals from a fast AME fit

Returns residuals from an
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) or
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
fit.

## Usage

``` r
# S3 method for class 'ame_als'
residuals(object, type = c("response", "working"), ...)
```

## Arguments

- object:

  an `ame_als` object.

- type:

  `"response"` (default) or `"working"`.

- ...:

  ignored.

## Value

A matrix (cross-sectional) or list of matrices (longitudinal).

## Details

`type = "response"` (default) returns observed `Y` minus the
response-scale fitted values, so
[`residuals()`](https://rdrr.io/r/stats/residuals.html) reconciles with
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html).
`type = "working"` returns the Gaussian working-scale residuals the
estimator actually minimised (identical to `"response"` for the `normal`
family).
