# Predict method for LAME models

Generate predictions from fitted longitudinal AME models. Returns a list
of matrices (one per time point) on the requested scale.

## Usage

``` r
# S3 method for class 'lame'
predict(object, type = c("response", "link"), ...)
```

## Arguments

- object:

  Fitted LAME model object

- type:

  Character; "response" or "link"

- ...:

  Additional arguments (not used)

## Value

List of prediction matrices (one per time point)
