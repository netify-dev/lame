# Post-MCMC per-actor time-varying slopes

Computes a smoothed per-actor time-varying slope coefficient on a nodal
covariate slice. For `kind = "row"`, each row-actor \\i\\ gets a
length-\\T\\ slope path \\\beta\_{i,t}\\ on `Xrow[i, covariate_idx, t]`,
fit by ridge-penalised least squares on the residual after the main MCMC
linear predictor.

## Usage

``` r
per_actor_slopes(fit, kind = c("row", "col"), covariate_idx = 1L, lambda = 1)
```

## Arguments

- fit:

  A fitted `lame` object.

- kind:

  `"row"` (default; per-row-actor slopes) or `"col"` (per-column-actor
  slopes).

- covariate_idx:

  Integer column index into the chosen nodal covariate (`Xrow` or
  `Xcol`). Default `1L`.

- lambda:

  Non-negative smoothing parameter (first-difference penalty across
  periods). Default `1`.

## Value

A list with `slopes` (an \\n\_{actors} \times T\\ matrix), `kind`,
`covariate_idx`, `lambda`, and `label`. Class `"per_actor_slopes"`.

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            nscan = 100, burn = 25, odens = 5, verbose = FALSE)
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
# post-MCMC per-row-actor slopes on the first dyadic covariate
pas <- per_actor_slopes(fit, kind = "row", lambda = 1)
dim(pas$slopes)
#> [1] 50  4
# }
```
