# Tidy method for fitted `ame_als` / `lame_als` objects

Returns a data frame with one row per regression coefficient, compatible
with the broom idiom, so that ALS fits compose with modelsummary /
kableExtra pipelines next to MCMC fits. Standard errors come from the
sandwich covariance
([`vcov.ame_als`](https://netify-dev.github.io/lame/reference/vcov.ame_als.md))
by default, or from the bootstrap object attached to `x$bootstrap` when
present (preferred, fully propagated). `statistic` is
`estimate / std.error`; `p.value` is the Normal-approximation two-sided
tail \\2(1 - \Phi(\|z\|))\\ - the same heuristic the MCMC `tidy.ame`
uses.

## Usage

``` r
# S3 method for class 'ame_als'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)

# S3 method for class 'lame_als'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A fitted `ame_als` / `lame_als` object.

- conf.int:

  Logical; include `conf.low` / `conf.high` columns. Default `TRUE`.

- conf.level:

  Confidence level. Default `0.95`.

- ...:

  Passed to
  [`vcov.ame_als`](https://netify-dev.github.io/lame/reference/vcov.ame_als.md)
  (e.g. `cluster = "dyad"`).

## Value

Data frame with columns `term`, `estimate`, `std.error`, `statistic`,
`p.value`, `conf.low`, `conf.high`, plus a `se_source` column recording
`"bootstrap"` or `"sandwich"`.

## Details

Only the intercept and dyadic-covariate coefficients are returned,
matching `coef(fit)` on the sandwich-covered subset. Additive (`a`,
`b`), multiplicative (`U`, `V`), and node-covariate parameters are not
included; use
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
and inspect the bootstrap object directly if you need them.

## Examples

``` r
# \donttest{
data(YX_bin_list)
Y1 <- 1 * (YX_bin_list$Y[[1]] > 0); diag(Y1) <- NA
fit <- ame_als(Y = Y1, Xdyad = YX_bin_list$X[[1]],
               family = "binary", R = 1, verbose = FALSE)
tidy(fit)
#> Warning: ! `vcov.ame_als()` for "binary" returns the conditional sandwich on the
#>   surrogate working likelihood -- anti-conservative.
#> ℹ Refit with `ame_als(..., bootstrap = 200)` for fully propagated uncertainty.
#> ℹ Or call `confint(fit)` after attaching a bootstrap with
#>   `ame_als_bootstrap()`.
#> # A tibble: 4 × 8
#>   term       estimate std.error statistic  p.value conf.low conf.high se_source
#>   <chr>         <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl> <chr>    
#> 1 intercept     0.236    0.0362      6.51 7.38e-11    0.165     0.307 sandwich 
#> 2 dyad1_dyad    0.823    0.0139     59.4  0           0.796     0.850 sandwich 
#> 3 dyad2_dyad    0.986    0.0156     63.2  0           0.956     1.02  sandwich 
#> 4 dyad3_dyad    1.19     0.0188     63.4  0           1.16      1.23  sandwich 
# }
```
