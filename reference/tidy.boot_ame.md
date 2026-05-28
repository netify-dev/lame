# Tidy method for a standalone bootstrap object (`boot_ame`)

[`ame_als_bootstrap()`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
returns an object of class `"boot_ame"` (not `"ame_als"`); this tidy
method exposes the bootstrap estimates as a broom-style data frame so
the standalone object composes with modelsummary the same way an
embedded `ame_als(..., bootstrap = N)` fit does.

## Usage

``` r
# S3 method for class 'boot_ame'
tidy(x, conf.level = 0.95, ...)
```

## Arguments

- x:

  A `boot_ame` object.

- conf.level:

  Confidence level. Default `0.95`.

- ...:

  Ignored.

## Value

Data frame with columns
`term, estimate, std.error, statistic, p.value, conf.low, conf.high, se_source`.
