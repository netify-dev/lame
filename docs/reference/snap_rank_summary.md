# Summarize posterior rank uncertainty for snap years

Ranks year-level snap indices within each retained draw, with rank 1
assigned to the highest snap index.

## Usage

``` r
snap_rank_summary(fit, actors = NULL, years = NULL, side = c("u", "v"))
```

## Arguments

- fit:

  a [`lame`](https://netify-dev.github.io/lame/reference/lame.md) fit
  with snap-shift MCMC output, or a
  [`lame_snap_als`](https://netify-dev.github.io/lame/reference/lame_snap_als.md)
  fit. MCMC fits should be run with retained snap draws for posterior
  uncertainty; ALS fits contribute their `snap_prob` scores as one score
  row.

- actors:

  optional actor names or indices. Default uses every actor on the
  selected side.

- years:

  optional year/period names or indices. Default uses every period.

- side:

  `"u"` for sender/row-side snap indicators or `"v"` for
  receiver/column-side indicators.

## Value

A data frame with year-level index means and rank probabilities.
