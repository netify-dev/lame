# Additive-effects plot for an ame_als fit

Sender (a) and receiver (b) additive-effect point estimates as a sorted
lollipop chart. Mirrors
[`ab_plot`](https://netify-dev.github.io/lame/reference/ab_plot.md) but
without posterior intervals: this is a point estimator. If the fit was
produced via `ame_als(..., bootstrap = N)`, bootstrap 95\\ error bars.

## Usage

``` r
ab_plot.ame_als(
  fit,
  effect = c("sender", "receiver", "both"),
  top_n = Inf,
  ...
)
```

## Arguments

- fit:

  an `ame_als` fit.

- effect:

  `"sender"` (default), `"receiver"`, or `"both"`.

- top_n:

  integer; show only the top / bottom `top_n` actors by absolute effect
  (default `Inf`, show all).

- ...:

  reserved.

## Value

A `ggplot` object.
