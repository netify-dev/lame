# Multivariate split-R-hat for dynamic_beta coefficient paths

Given one or more multi-chain `lame` fits whose `$BETA` is
`[n_iter, p, T]`, computes the multivariate R-hat per coefficient \\k\\
treating the length-\\T\\ path as one multivariate observation per
iteration.

## Usage

``` r
rhat_dynamic_beta(fit_list, coefs = NULL)
```

## Arguments

- fit_list:

  A list of fitted `lame` objects from `lame_parallel(..., chains = K)`,
  or any `ame_chain_list` produced by re-running
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) with
  different seeds.

- coefs:

  Optional character vector of coefficient names to subset (matches
  `dimnames(fit$BETA)[[2]]`).

## Value

Data frame with one row per coefficient: `coef`, `rhat_mvt`,
`rhat_max_univariate`, `n_chains`, `n_iter_per_chain`, `T`. `rhat_mvt`
is the Brooks-Gelman multivariate statistic; `rhat_max_univariate` is
the max over per-(k,t) split-R-hat values for comparison.

## Details

For chain \\c\\, let \\\beta^{(c)}\_t \in \mathbb{R}^T\\ be the path.
With \\m\\ chains and \\n\\ iterations per chain, define \$\$W_k =
\tfrac{1}{m}\sum_c S_c^{(k)}, \quad B_k = \tfrac{n}{m-1}\sum_c
(\bar\beta_c^{(k)} - \bar\beta^{(k)})(\bar\beta_c^{(k)} -
\bar\beta^{(k)})'\$\$ where \\S_c^{(k)}\\ is the within-chain sample
covariance of path \\k\\ in chain \\c\\. Then \\V_k = ((n-1)/n) W_k +
((m+1)/(mn)) B_k\\ and \\\hat R_k^{mvt} = \sqrt{\lambda\_\max(W_k^{-1}
V_k)}\\. For nearly degenerate covariances we add a tiny ridge to
\\W_k\\.

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit_list <- lapply(c(1L, 2L, 3L, 4L), function(s) {
  set.seed(s)
  lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
       dynamic_beta = "dyad",
       nscan = 200, burn = 50, odens = 5, verbose = FALSE)
})
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
rhat_dynamic_beta(fit_list)
#>        coef  rhat_mvt rhat_max_univariate n_chains n_iter_per_chain T
#> 1 intercept 0.9874209            1.237604        4               40 4
#> 2   X1_dyad 0.9874209            2.235824        4               40 4
#> 3   X2_dyad 0.9874209            1.945086        4               40 4
#> 4   X3_dyad 0.9874209            2.070859        4               40 4
# }
```
