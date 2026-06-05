# Multi-panel lame() with shared coefficients

Fits K independent
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) models
(one per panel) and pools the per-panel beta posteriors into a
precision-weighted shared posterior. Returns a list with the per-panel
fits, the pooled beta posterior, and the panel-specific deviations.

## Usage

``` r
lame_multi(Y_list, Xdyad_list, ..., shared_dynamic_beta = TRUE)
```

## Arguments

- Y_list:

  A list of length K, each element a list (or 3-D array) of T per-panel
  network observations.

- Xdyad_list:

  A list of length K, each element a list of T dyadic covariate arrays.

- ...:

  Arguments forwarded to
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) (e.g.
  family, R, mode, nscan, burn, odens, dynamic_beta, dynamic_beta_kind).

- shared_dynamic_beta:

  Logical; if `TRUE` (default), the shared beta posterior is computed
  (per-period when dynamic_beta is active).

## Value

A list with

- `fits`: list of K per-panel `lame` fits.

- `beta_shared`: pooled posterior mean of beta (per-period when
  dynamic).

- `beta_deviations`: list of K panel-specific deviations from
  `beta_shared`.

- `K`: number of panels.

Class `"lame_multi"`.

## Details

This is an R-level wrapper. A single joint MCMC would share panel state
across iterations, which is a substantial sampler refactor. The wrapper
is exact when the panels are conditionally independent given beta, which
is the standard assumption.

## See also

[`lame_parallel`](https://netify-dev.github.io/lame/reference/lame_parallel.md)
for the unrelated multi-*chain* wrapper that runs K MCMC chains of the
*same* model (used for R-hat / ESS diagnostics and pooled effective
sample size). `lame_multi` is for K *distinct* panels with shared
regression coefficients; `lame_parallel` is for K chains of one model.

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit_multi <- lame_multi(
  Y_list = list(YX_bin_list$Y, YX_bin_list$Y),
  Xdyad_list = list(YX_bin_list$X, YX_bin_list$X),
  family = "binary", R = 0,
  nscan = 100, burn = 25, odens = 5, verbose = FALSE)
dim(fit_multi$beta_shared)
#> NULL
# }
```
