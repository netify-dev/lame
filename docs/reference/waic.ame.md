# WAIC for AME / LAME fits

S3 method for [`waic`](https://mc-stan.org/loo/reference/waic.html) that
uses the stored `fit$log_lik`.

## Usage

``` r
# S3 method for class 'ame'
waic(x, ...)

# S3 method for class 'lame'
waic(x, ...)

# S3 method for class 'ame_als'
waic(x, ...)
```

## Arguments

- x:

  A fitted `ame` or `lame` object with `$log_lik`.

- ...:

  Additional arguments forwarded to
  [`loo::waic.matrix`](https://mc-stan.org/loo/reference/waic.html).

## Value

A `waic` object.
