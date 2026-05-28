# fitted/residuals are not defined for a bootstrap object

A `boot_ame` stores replicate-level distributional summaries, not
per-cell predictions or residuals. These methods exist only to refuse
the call clearly (instead of falling through to
[`stats::fitted.default`](https://rdrr.io/r/stats/fitted.values.html) /
[`stats::residuals.default`](https://rdrr.io/r/stats/residuals.html),
which would silently return `NULL`). Apply
[`fitted`](https://rdrr.io/r/stats/fitted.values.html) /
[`residuals`](https://rdrr.io/r/stats/residuals.html) to the underlying
`ame_als` fit instead.

## Usage

``` r
# S3 method for class 'boot_ame'
fitted(object, ...)

# S3 method for class 'boot_ame'
residuals(object, ...)
```

## Arguments

- object:

  a `boot_ame` object.

- ...:

  ignored.

## Value

Never returns; raises an error.
