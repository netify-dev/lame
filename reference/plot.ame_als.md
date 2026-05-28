# Plot the convergence of a fast AME fit

Plots the block coordinate descent deviance/SSE history – a convergence
diagnostic for an
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) or
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
fit. For substantive plots see
[`uv_plot`](https://netify-dev.github.io/lame/reference/uv_plot.md),
[`ab_plot`](https://netify-dev.github.io/lame/reference/ab_plot.md) and
[`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md).

## Usage

``` r
# S3 method for class 'ame_als'
plot(x, ...)
```

## Arguments

- x:

  an `ame_als` object.

- ...:

  further arguments passed to the underlying plot.

## Value

`x`, invisibly.
