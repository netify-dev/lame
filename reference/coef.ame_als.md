# Extract coefficients from a fast AME fit

Returns the point-estimated regression coefficients (intercept first) of
an `ame_als` fit. There is no posterior here; for uncertainty use
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md).

## Usage

``` r
# S3 method for class 'ame_als'
coef(object, ...)
```

## Arguments

- object:

  an `ame_als` object.

- ...:

  ignored.

## Value

A named numeric vector of coefficients.
