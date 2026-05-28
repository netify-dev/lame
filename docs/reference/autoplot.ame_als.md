# autoplot method for ALS fits

Coefficient point-and-interval plot for an
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) or
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
fit. The body builds the same frame
[`autoplot.lame`](https://netify-dev.github.io/lame/reference/autoplot.lame.md)
uses, sourcing the point estimate from `coef(fit)` and the interval from
`confint(fit)` (sandwich or bootstrap, depending on which is available).
`which = "beta"` is the only mode supported; `which = "uv"` / `"ab"`
return an informative error.

## Usage

``` r
# S3 method for class 'ame_als'
autoplot(object, which = c("beta", "uv", "ab"), conf.level = 0.95, ...)

# S3 method for class 'lame_als'
autoplot(object, which = c("beta", "uv", "ab"), conf.level = 0.95, ...)
```

## Arguments

- object:

  A fitted `ame_als` / `lame_als` object.

- which:

  One of `"beta"` (currently the only supported value).

- conf.level:

  Confidence level for the interval. Default `0.95`.

- ...:

  Passed to
  [`vcov.ame_als`](https://netify-dev.github.io/lame/reference/vcov.ame_als.md)
  /
  [`confint.ame_als`](https://netify-dev.github.io/lame/reference/confint.ame_als.md)
  (e.g. `cluster`).

## Value

A `ggplot` object.
