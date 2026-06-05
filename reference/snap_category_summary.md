# Summarize snap indices by actor category

Summarize snap indices by actor category

## Usage

``` r
snap_category_summary(fit, groups, years = NULL, side = c("u", "v"))
```

## Arguments

- fit:

  a [`lame`](https://netify-dev.github.io/lame/reference/lame.md) fit
  with snap-shift MCMC output, or a
  [`lame_snap_als`](https://netify-dev.github.io/lame/reference/lame_snap_als.md)
  fit. MCMC fits should be run with retained snap draws for posterior
  uncertainty; ALS fits contribute their `snap_prob` scores as one score
  row.

- groups:

  a named list mapping group names to actor names/indices, or a named
  vector whose values are group labels and whose names are actors.

- years:

  optional year/period names or indices. Default uses every period.

- side:

  `"u"` for sender/row-side snap indicators or `"v"` for
  receiver/column-side indicators.

## Value

A data frame with one row per group-year.
