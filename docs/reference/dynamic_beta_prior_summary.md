# Summarise the implied prior on a time-varying coefficient path

Draws `ndraws` sample paths of length `T` from the `dynamic_beta` prior
(given a kind, AR(1) hyperparameters, and an inverse-gamma on the
innovation variance) and reports the implied distribution of common
summary statistics: maximum absolute first difference, roughness (sum of
squared first differences), path range, and correlation with time. Use
this before fitting to confirm that your prior is loose / tight enough.

## Usage

``` r
dynamic_beta_prior_summary(
  T = 10L,
  ndraws = 5000L,
  kind = c("ar1", "rw1", "rw2", "matern32"),
  rho_mean = 0.8,
  rho_sd = 0.15,
  rho_lower = 0,
  rho_upper = 0.999,
  sigma_shape = 2,
  sigma_scale = 1,
  matern32_length_scale = 2,
  threshold = 1,
  seed = NULL
)
```

## Arguments

- T:

  Length of the path (default 10).

- ndraws:

  Number of paths to draw (default 5000).

- kind:

  One of `"ar1"`, `"rw1"`, `"rw2"`, `"matern32"`. Default `"ar1"`. For
  `"matern32"` the path is drawn from a continuous-time Matérn-3/2 GP
  with length scale `matern32_length_scale` and marginal variance
  \\\sigma\_\beta^2\\.

- rho_mean, rho_sd:

  Prior mean/SD on \\\rho\\ (AR(1) only). Used to parameterise a Beta
  prior on the standardised \\\rho\\.

- rho_lower, rho_upper:

  Truncation bounds (AR(1) only). Default `[0, 0.999]`.

- sigma_shape, sigma_scale:

  Inverse-Gamma shape / scale on \\\sigma\_\beta^2\\. Defaults
  `shape = 2, scale = 1`.

- matern32_length_scale:

  Length scale for the Matérn-3/2 covariance (`kind = "matern32"` only).
  Default `2`.

- threshold:

  Cutoff for the `prob_max_diff_gt_threshold` summary. Default `1`.

- seed:

  Optional integer seed for reproducibility.

## Value

A list with components:

- `summary`:

  Data frame with quantile rows for max-first-diff, roughness, path
  range, trend correlation.

- `prob_max_diff_gt_threshold`:

  Empirical probability that the maximum absolute first difference
  exceeds `threshold`.

- `rho`:

  Length `ndraws` vector of sampled rho values (`NA` for kinds with
  fixed rho).

- `sigma`:

  Length `ndraws` vector of sampled sigma values.

- `paths`:

  `ndraws x T` matrix of sample paths.

## Details

For `kind = "ar1"`, `rho_mean`/`rho_sd` are translated to a Beta prior
on the standardised \\(\rho - \rho\_{lower}) / (\rho\_{upper} -
\rho\_{lower})\\. For `kind = "rw1"`, \\\rho = 1\\ is fixed and only
\\\sigma\_\beta^2\\ is sampled. For `kind = "rw2"`, the
second-difference variance is sampled and the first two values use a
diffuse \\N(0, 10)\\ initial prior.

## Examples

``` r
# \donttest{
# Default prior: AR(1) with rho_mean = 0.8, rho_sd = 0.15
s <- dynamic_beta_prior_summary(T = 10, ndraws = 2000, seed = 1)
s$summary
#>             quantity        q05        q50        q95
#> 1 max_abs_first_diff  0.7215172 1.46787811  3.7138105
#> 2          roughness  1.4113939 5.56455715 31.4268861
#> 3              range  1.1631870 2.62820813  6.6322426
#> 4          trend_cor -0.8732470 0.01310171  0.8807169
# what's the probability that consecutive beta_t differ by more than 1?
s$prob_max_diff_gt_threshold
#> [1] 0.799

# Tighter prior on innovation variance
s_tight <- dynamic_beta_prior_summary(sigma_scale = 0.1, seed = 1)
s_tight$summary
#>             quantity        q05        q50       q95
#> 1 max_abs_first_diff  0.2176577 0.46198529 1.1278813
#> 2          roughness  0.1372637 0.55602950 3.0626640
#> 3              range  0.3659379 0.82262847 2.0530328
#> 4          trend_cor -0.8787253 0.02435639 0.8788859
# }
```
