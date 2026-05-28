# Point estimates from a fast AME bootstrap

Returns the intercept and regression coefficient point estimates carried
by a `boot_ame` object (the estimates the bootstrap quantifies). For
their bootstrap standard errors and intervals see
[`confint.boot_ame`](https://netify-dev.github.io/lame/reference/confint.boot_ame.md)
and
[`vcov.boot_ame`](https://netify-dev.github.io/lame/reference/vcov.boot_ame.md).

## Usage

``` r
# S3 method for class 'boot_ame'
coef(object, ...)
```

## Arguments

- object:

  a `boot_ame` object.

- ...:

  ignored.

## Value

A named numeric vector of coefficient point estimates.
