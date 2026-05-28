# Glance method for fitted `ame_als` / `lame_als` objects

One-row data frame summarising an ALS fit, compatible with broom /
modelsummary. Reports observation count, actor count, latent dimension,
family, mode, ALS convergence flag, and iteration count.

## Usage

``` r
# S3 method for class 'ame_als'
glance(x, ...)

# S3 method for class 'lame_als'
glance(x, ...)
```

## Arguments

- x:

  A fitted `ame_als` / `lame_als` object.

- ...:

  Ignored.

## Value

One-row data frame with columns `nobs`, `n_actors`, `n_periods`,
`family`, `mode`, `R`, `converged`, `iterations`, `se_source`.

## Examples

``` r
# \donttest{
data(YX_bin_list)
Y1 <- 1 * (YX_bin_list$Y[[1]] > 0); diag(Y1) <- NA
fit <- ame_als(Y = Y1, Xdyad = YX_bin_list$X[[1]],
               family = "binary", R = 1, verbose = FALSE)
glance(fit)
#> # A tibble: 1 × 12
#>    nobs n_actors n_row_actors n_col_actors n_periods family mode           R
#>   <int>    <int>        <int>        <int>     <int> <chr>  <chr>      <int>
#> 1  2450       50           NA           NA         1 binary unipartite     1
#> # ℹ 4 more variables: non_normal_method <chr>, converged <lgl>,
#> #   iterations <int>, se_source <chr>
# }
```
