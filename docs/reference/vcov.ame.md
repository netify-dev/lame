# Posterior covariance of AME model coefficients

Returns the posterior covariance matrix of regression coefficients.

## Usage

``` r
# S3 method for class 'ame'
vcov(object, ...)

# S3 method for class 'lame'
vcov(object, ...)
```

## Arguments

- object:

  fitted AME model (class "ame")

- ...:

  additional arguments (ignored)

## Value

p x p covariance matrix of posterior BETA draws
