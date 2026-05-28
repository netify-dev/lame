# Detect potential change points in a dynamic_beta posterior path

For each dynamic coefficient in a `dynamic_beta` fit, compares the
posterior distribution of the maximum scaled first-difference \\M =
\max_t \|\beta_t - \beta\_{t-1}\| / \sigma\_\beta\\ to a Monte-Carlo
approximation of the same statistic under the AR(1) prior. A
“change-point-flavored” Bayes factor is reported as \\BF = \Pr(M^{post}
\> m^\*) / 0.05\\, where \\m^\*\\ is the 95\\ quantile of \\M^{prior}\\.
**This is a heuristic, not a real Bayes factor.** Use it to surface
posterior temporal jumps that the AR(1) prior cannot comfortably
accommodate.

## Usage

``` r
detect_change_point(
  fit,
  coefs = NULL,
  n_prior_sims = 2000L,
  threshold_bf = 10,
  seed = NULL
)
```

## Arguments

- fit:

  A fitted `lame` object with `dynamic_beta`.

- coefs:

  Optional character vector of coefficient names to check. Defaults to
  every dynamic coefficient on the fit.

- n_prior_sims:

  Number of prior simulations for the null distribution. Default `2000`.

- threshold_bf:

  Cutoff for the "warn" column. Default `10`; lower to 5 if the
  false-positive rate on AR(1) data is too high.

- seed:

  Optional seed for reproducibility.

## Value

A data frame with one row per dynamic coefficient: `coef`, `bf` (the
heuristic Bayes factor), `m_post_mean` (posterior mean of \\M\\),
`m_prior_q95` (95\\ \\M\\), `t_hat` (period of the largest scaled jump),
`warn` (logical, `bf > threshold_bf`).

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
detect_change_point(fit)
#>      coef bf m_post_mean m_prior_q95 t_hat  warn
#> 1 X1_dyad  0   0.1279427    2.495379     3 FALSE
#> 2 X2_dyad  0   0.1641786    2.582841     4 FALSE
#> 3 X3_dyad  0   0.2464007    2.569418     4 FALSE
# }
```
