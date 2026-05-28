# Predictions from a fast AME fit

Returns predictions from an
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) or
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
fit: the fitted linear predictor on the link scale (`type = "link"`) or
the response scale (`type = "response"`, the default).

## Usage

``` r
# S3 method for class 'ame_als'
predict(object, newdata = NULL, type = c("response", "link"), ...)
```

## Arguments

- object:

  an `ame_als` object.

- newdata:

  optional dyadic covariate array matching `object$X`
  (`[n_row, n_col, p, T]`); a matrix or 3D array is coerced.

- type:

  `"response"` (default) or `"link"`.

- ...:

  ignored.

## Value

A matrix (cross-sectional) or list of matrices (longitudinal).

## Details

With `newdata = NULL` the training-data fitted values are returned.
Supplying `newdata` – a dyadic covariate array with the *same* actors
and dimensions as the fitted design – substitutes the dyadic covariate
contribution while holding the intercept, additive effects and
multiplicative term fixed; it predicts counterfactual dyadic covariates
for the same network, not out-of-sample actors.
