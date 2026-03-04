# Extract residuals from LAME model

Extract residuals from LAME model

## Usage

``` r
# S3 method for class 'lame'
residuals(object, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  Fitted LAME model

- type:

  Type of residuals ("response" or "pearson")

- ...:

  Additional arguments

## Value

List of residual matrices (one per time point)
