# Glance method for fitted `ame` / `lame` objects

One-row data frame summarising model-level statistics, in the broom
idiom. Used by
[`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html)
and similar tabling tools to populate the lower goodness-of-fit panel of
a regression table.

## Usage

``` r
# S3 method for class 'ame'
glance(x, ...)

# S3 method for class 'lame'
glance(x, ...)
```

## Arguments

- x:

  A fitted `ame` / `lame` object.

- ...:

  Ignored.

## Value

A one-row data frame with columns:

- `nobs` — number of observed dyads (NA cells excluded).

- `n_actors` — number of distinct actors (for bipartite, row + column
  actors).

- `n_periods` — number of time periods (1 for `ame`).

- `n_stored` — number of stored MCMC draws.

- `family` — outcome family ("normal","binary",...).

- `mode` — "unipartite","bipartite".

- `R` — latent-space dimension (or max of `R_row`, `R_col` for
  bipartite).

- `dynamic_uv`, `dynamic_ab`, `dynamic_beta` — logicals; whether each
  component is time-varying.

- `elpd_loo` — leave-one-out expected log predictive density, if a
  cached `loo` object is on the fit (NA if not fit with
  `save_log_lik = TRUE`).

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            nscan = 100, burn = 20, odens = 5, verbose = FALSE)
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
glance(fit)
#> # A tibble: 1 × 13
#>    nobs n_actors n_row_actors n_col_actors n_periods n_stored family mode      R
#>   <int>    <int>        <int>        <int>     <int>    <int> <chr>  <chr> <int>
#> 1  9800       50           NA           NA         4       20 binary unip…     0
#> # ℹ 4 more variables: dynamic_uv <lgl>, dynamic_ab <lgl>, dynamic_beta <lgl>,
#> #   elpd_loo <dbl>
# }
```
