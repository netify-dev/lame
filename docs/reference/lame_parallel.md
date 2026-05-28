# Run LAME (longitudinal AME) with multiple parallel chains

Thin wrapper around
[`ame_parallel`](https://netify-dev.github.io/lame/reference/ame_parallel.md)
that forces `fitter = "lame"`. Use this for longitudinal data; in
particular, multi-chain `dynamic_beta` / `dynamic_uv` / `dynamic_ab`
fits are reached through here. The combined fit's `$BETA` is a 3-D array
when `dynamic_beta` is on, and `posterior::rhat(as_draws(fit))` gives
R-hat across chains for every per-period coefficient.

## Usage

``` r
lame_parallel(
  Y,
  n_chains = 4,
  cores = n_chains,
  combine_method = c("pool", "list"),
  ...
)
```

## Arguments

- Y:

  Longitudinal network: list of T relational matrices, or a 3-D array
  `[n, n, T]`.

- n_chains:

  Number of parallel chains (default 4).

- cores:

  CPU cores to use (default `n_chains`; 1 = sequential).

- combine_method:

  `"pool"` (default) or `"list"`.

- ...:

  Additional arguments forwarded to
  [`lame`](https://netify-dev.github.io/lame/reference/lame.md),
  including `dynamic_beta`, `dynamic_uv`, `dynamic_ab`, `family`,
  `nscan`, `burn`, `odens`, etc.

## Value

Same as `ame_parallel`: a combined `lame` fit (when
`combine_method = "pool"`) or a list of fits.

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit_pll <- lame_parallel(YX_bin_list$Y, Xdyad = YX_bin_list$X,
                         family = "binary", n_chains = 2, cores = 1,
                         nscan = 100, burn = 20, odens = 5,
                         dynamic_beta = "dyad", verbose = FALSE)
#> 
#> ── Running 2 chains sequentially ──
#> 
#> Starting chain 1 (`lame()`)
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
#> Completed chain 1
#> Starting chain 2 (`lame()`)
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
#> Completed chain 2
#> Combining chains...
#> 
#> ── MCMC Convergence Diagnostics 
#> Number of chains: 2
#> Samples per chain: "320, 320"
#> Diagnostics cover regression coefficients and variance components.
#> ! 1 parameter have R-hat >= 1.1
#> rho: R-hat = 3.119
dim(fit_pll$BETA)        # [iter * n_chains, p, T]
#> [1] 40  4  4
fit_pll$chain_indicator  # length iter * n_chains
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> [39] 2 2
if (requireNamespace("posterior", quietly = TRUE)) {
  posterior::summarise_draws(posterior::as_draws(fit_pll))
}
#> # A tibble: 21 × 10
#>    variable       mean median     sd    mad     q5   q95  rhat ess_bulk ess_tail
#>    <chr>         <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>    <dbl>    <dbl>
#>  1 intercept[t… 0.0847 0.0832 0.0255 0.0279 0.0438 0.126  1.14    15.1      49.6
#>  2 X1_dyad[t1]  0.557  0.559  0.0804 0.0833 0.427  0.670  1.75     4.62     22.4
#>  3 X2_dyad[t1]  0.657  0.678  0.0999 0.0971 0.480  0.781  1.70     4.71     49.6
#>  4 X3_dyad[t1]  0.804  0.804  0.116  0.118  0.613  0.959  1.83     4.43     17.3
#>  5 intercept[t… 0.0847 0.0832 0.0255 0.0279 0.0438 0.126  1.14    15.1      49.6
#>  6 X1_dyad[t2]  0.547  0.551  0.0795 0.0835 0.414  0.680  1.65     4.85     18.9
#>  7 X2_dyad[t2]  0.680  0.689  0.0967 0.113  0.517  0.802  2.04     4.20     37.7
#>  8 X3_dyad[t2]  0.833  0.840  0.121  0.153  0.638  0.997  1.88     4.35     17.3
#>  9 intercept[t… 0.0847 0.0832 0.0255 0.0279 0.0438 0.126  1.14    15.1      49.6
#> 10 X1_dyad[t3]  0.540  0.530  0.0828 0.0894 0.399  0.666  1.75     4.58     17.3
#> # ℹ 11 more rows
# }
```
