# Number of observed dyads in an AME / LAME fit

Number of finite (non-missing) cells in `Y`, summed over time slices for
longitudinal fits. Useful as a denominator for sample-size reporting and
as the basis for AIC/BIC *if* you also have a log-likelihood (which
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) /
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) do not
directly expose; see
[`logLik.ame`](https://netify-dev.github.io/lame/reference/logLik.ame.md)).

## Usage

``` r
# S3 method for class 'ame'
nobs(object, ...)

# S3 method for class 'lame'
nobs(object, ...)
```

## Arguments

- object:

  an `ame` or `lame` fit.

- ...:

  ignored.

## Value

Integer: number of observed dyads in `Y`.
