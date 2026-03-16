# Simulate posterior distributions from fitted AME model

Simulate posterior distributions from fitted AME model

## Usage

``` r
simulate_posterior(
  fit,
  component = c("UV", "ab", "beta", "Y"),
  n_samples = 100,
  use_full_cov = FALSE,
  seed = NULL
)
```

## Arguments

- fit:

  Fitted ame model object

- component:

  Character; which component to simulate: "UV", "ab", "beta", "Y"

- n_samples:

  Number of posterior samples to generate

- use_full_cov:

  For UV component, whether to use empirical variance scaling (default
  FALSE)

- seed:

  Random seed for reproducibility

## Value

Array or matrix of posterior samples

## Details

This function can simulate posterior distributions even when they
weren't saved during MCMC, by using the posterior means and variance
components.

For more accurate posteriors, use posterior_options() during model
fitting to save the actual MCMC samples.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Fit a model with multiplicative effects
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Get posterior samples of regression coefficients
beta_post <- simulate_posterior(fit, "beta", n_samples = 50)
#> Using saved MCMC samples for beta
# }
```
