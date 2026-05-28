# Exact rolling-origin leave-future-out cross-validation

For a fitted `lame` object with \\T\\ periods, refits the model on the
first \\t - 1\\ periods (for each \\t\\ in `periods`) and computes the
expected log predictive density (elpd) of period \\t\\ under the refit.
Returns the summed elpd across all leave-out periods, along with
per-period contributions.

## Usage

``` r
lfo(fit, periods = NULL, refit = TRUE, ...)
```

## Arguments

- fit:

  A fitted `lame` object.

- periods:

  Integer vector of leave-out periods to evaluate. Default is the last 3
  periods (`tail(seq_len(T), 3L)`). Each period \\t\\ must satisfy \\t
  \ge 2\\.

- refit:

  Logical; if `TRUE` (default), refits on the training window. If
  `FALSE`, uses the original posterior means (a much rougher
  approximation).

- ...:

  Passed to the refit
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) call
  (typically `nscan`, `burn`, `odens`, `verbose`).

## Value

A list with `elpd_lfo` (total summed elpd), `pointwise` (per-dyad
log-density at each leave-out period; a list of numeric vectors, one per
period — `unlist(pointwise)` gives a flat vector suitable for
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)-style
stacking), `p_lfo` (effective number of parameters), `per_period` (data
frame with `period`, `elpd`, `n_obs`), and `periods` (the periods
evaluated).

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
lfo_res <- lfo(fit, periods = c(3L, 4L), refit = TRUE,
               nscan = 100, burn = 20, odens = 5, verbose = FALSE)
#> Warning: ! Posterior 95% interval for `rho_beta` reaches 0.991.
#> ℹ Under AR(1) the `h`-step forecast variance saturates at the stationary level
#>   `sigma_beta^2 / (1 - rho_beta^2)`, which itself explodes as `rho_beta -> 1`;
#>   horizons beyond "3" become uninformative in this regime.
#> ℹ Consider refitting with `dynamic_beta_kind = "rw1"` if the underlying process
#>   is genuinely unit-root.
#> Warning: ! Posterior 95% interval for `rho_beta` reaches 0.991.
#> ℹ Under AR(1) the `h`-step forecast variance saturates at the stationary level
#>   `sigma_beta^2 / (1 - rho_beta^2)`, which itself explodes as `rho_beta -> 1`;
#>   horizons beyond "3" become uninformative in this regime.
#> ℹ Consider refitting with `dynamic_beta_kind = "rw1"` if the underlying process
#>   is genuinely unit-root.
print(lfo_res)
#> 
#> ── Exact rolling-origin LFO ──
#> 
#> Periods evaluated: 3 and 4
#> Refit per leave-out: TRUE
#> Total elpd_lfo: -17952
#>  period      elpd n_obs
#>       3 -9729.205  2450
#>       4 -8222.901  2450
# }
```
