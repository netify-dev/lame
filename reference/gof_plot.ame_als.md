# Goodness-of-fit check for an ame_als fit

Bootstrap-style analogue of the MCMC
[`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md):
draws `nsim` simulated networks from the fitted ALS model, computes the
standard network statistics (`sd.rowmean`, `sd.colmean`, `dyad.dep`,
`cycle.dep`, `trans.dep` for unipartite; `sd.rowmean`, `sd.colmean`,
`four.cycles` for bipartite) on each replicate, and overlays the
observed value on a histogram of replicate values.

## Usage

``` r
gof_plot.ame_als(fit, nsim = 100, seed = NULL, ...)
```

## Arguments

- fit:

  an `ame_als` fit.

- nsim:

  integer; number of replicates (default 100).

- seed:

  optional RNG seed.

- ...:

  reserved.

## Value

A `ggplot` (or `patchwork`) object.

## Details

Uses
[`simulate.ame_als`](https://netify-dev.github.io/lame/reference/simulate.ame_als.md)
for the replicates and is therefore subject to the same caveat: the
noise is resampled, but the point estimates of `mu, beta, a, b, U, V`
are held fixed. The MCMC `gof_plot.ame` integrates over the posterior of
those parameters too, so the ALS GOF is a weaker check.
