# Bootstrap covariance of the regression coefficients

Returns the sample covariance matrix of the bootstrap replicate
intercept + regression coefficients of a `boot_ame` object – the
covariance underlying the reported standard errors and confidence
intervals.

## Usage

``` r
# S3 method for class 'boot_ame'
vcov(object, ...)
```

## Arguments

- object:

  a `boot_ame` object.

- ...:

  ignored.

## Value

A covariance matrix over the intercept and regression coefficients.
