# Posterior-predictive temporal-trend test

For a fitted `lame` object, computes a network statistic (`density`,
`reciprocity`, or `transitivity`) at each observed period, fits a
least-squares linear trend on period index, and compares the observed
slope to slopes from posterior-predictive replicates. The two-sided
value is `p_pp = 2 * min(p_up, 1 - p_up)`, which runs from 0 (observed
slope in the extreme tail) to 1 (observed slope dead-centre). A static
fit on truly trending data yields `p_pp` near 0; a dynamic fit that
captures the trend yields `p_pp` near 1.

## Usage

``` r
gof_temporal(
  fit,
  stat = c("auto", "density", "mean", "reciprocity", "transitivity"),
  n_rep = 500,
  seed = NULL
)
```

## Arguments

- fit:

  A fitted `lame` object.

- stat:

  One of `"auto"` (default; picks `"reciprocity"` for unipartite
  directed fits with \\T \ge 3\\ and family `normal` / `poisson` /
  `tobit`, where `"density"` is constant and uninformative; otherwise
  `"density"`), `"density"`, `"mean"` (mean of off-diagonal Y; the right
  "density" analogue for continuous outcomes), `"reciprocity"`,
  `"transitivity"`.

- n_rep:

  Number of posterior-predictive replicates to draw (each replicate is a
  full \\T\\-period network from `simulate(fit)`).

- seed:

  Optional RNG seed.

## Value

A list with

- `stat` (chosen statistic name),

- `slope_obs` (observed slope of stat on period index),

- `slope_rep` (length-`n_rep` vector of replicate slopes),

- `p_pp` (two-sided posterior-predictive p-value),

- `stat_obs_by_t`, `stat_rep_by_t` (per-period statistics)

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            dynamic_beta = "dyad",
            nscan = 200, burn = 50, odens = 5, verbose = FALSE)
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
gof_temporal(fit, stat = "density", n_rep = 100)
#> 
#> ── Temporal-trend posterior-predictive check ──
#> 
#> • Statistic: "density"
#> • Observed slope (per-period): -0.00068
#> • Replicates: 100
#> • Posterior-predictive p-value (two-sided): 0.82
#> Observed temporal trend is well covered by the fitted model.
# }
```
