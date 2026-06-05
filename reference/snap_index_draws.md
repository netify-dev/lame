# Extract posterior draws of snap indices

Computes draw-level averages of MCMC snap indicators over selected
actors and years. The result is useful for system-level rupture indices
and focal-actor rupture indices. Fits should be run with
`keep_snap_draws = "draws"` for posterior uncertainty. If only
`snap_prob` is present, the function returns a single score row and
warns that intervals cannot be read as posterior intervals.

## Usage

``` r
snap_index_draws(fit, actors = NULL, years = NULL, side = c("u", "v"))
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

A data frame with columns `draw`, `year`, `side`, `n_actors`, and
`snap_index`.
