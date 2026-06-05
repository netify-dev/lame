# Tidy method for fitted `ame` / `lame` objects

Returns a data frame with one row per estimated coefficient, compatible
with the broom idiom. For `dynamic_beta` fits (3-D `BETA`), returns one
row per coefficient *per period* with a `period` column. Standard errors
are posterior standard deviations; `statistic` is
`estimate / std.error`.

## Usage

``` r
# S3 method for class 'ame'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)

# S3 method for class 'lame'
tidy(x, conf.int = TRUE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A fitted `ame` / `lame` object.

- conf.int:

  Logical; include 95\\ (`conf.low`, `conf.high`). Default `TRUE`.

- conf.level:

  Confidence level for the interval. Default `0.95`.

- ...:

  Ignored.

## Value

Data frame with columns `term`, `estimate`, `std.error`, `statistic`,
`p.value`, `conf.low`, `conf.high`, and (for dynamic_beta fits)
`period`.

## Details

**Note on `p.value`.** This column is included for broom compatibility
but is *not* a classical test. It is a two-sided Normal approximation
based on the posterior mean and marginal posterior standard deviation,
matching the calculation in `summary(fit)`. Use it as a compact signal
that the marginal posterior is far from zero, and report it alongside
the `conf.low` / `conf.high` credible interval. When sign certainty
matters, compute it directly from `x$BETA`, for example
`mean(sign(BETA) == sign(mean(BETA)))`.

Loaded as an S3 method against
[`generics::tidy`](https://generics.r-lib.org/reference/tidy.html) when
the generics package is available; works as `tidy(fit)` either way once
broom is loaded.

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            nscan = 100, burn = 20, odens = 5, verbose = FALSE)
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
tidy(fit)
#> # A tibble: 4 × 7
#>   term      estimate std.error statistic  p.value conf.low conf.high
#>   <chr>        <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 intercept   0.0539    0.109      0.493 6.22e- 1   -0.145     0.215
#> 2 X1_dyad     0.510     0.0704     7.25  4.30e-13    0.389     0.610
#> 3 X2_dyad     0.625     0.0810     7.71  1.24e-14    0.488     0.738
#> 4 X3_dyad     0.771     0.102      7.59  3.13e-14    0.596     0.913
# }
```
