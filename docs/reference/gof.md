# Compute GOF statistics from saved posterior samples

Computes goodness-of-fit statistics after model estimation using saved
posterior samples. This avoids the computational overhead of GOF
calculation during MCMC sampling. Requires that the model was run with
appropriate posterior sampling options.

## Usage

``` r
gof(fit, Y = NULL, custom_gof = NULL, nsim = 100, verbose = TRUE)
```

## Arguments

- fit:

  An ame model object that was run with posterior sampling enabled

- Y:

  Original data matrix (if not stored in fit)

- custom_gof:

  Optional custom GOF function(s) - same format as for ame()

- nsim:

  Number of posterior predictive simulations to generate (default 100).
  If NULL, uses all available posterior samples.

- verbose:

  Logical; print progress information

## Value

A matrix of GOF statistics with the same format as if gof=TRUE was used
during model estimation. First row contains observed statistics,
subsequent rows contain posterior predictive statistics.

## Details

This function requires that the model was estimated with posterior
sampling of the parameters needed to generate posterior predictive
datasets. Specifically, it needs:

- BETA: regression coefficients

- VC: variance components

- For models with random effects: samples of a, b

- For models with latent factors: U_samples, V_samples

To enable posterior sampling during model estimation, use:
`posterior_opts = posterior_options(save_UV = TRUE, save_ab = TRUE)`

Computing GOF post-hoc has several advantages:

- Faster MCMC sampling (no GOF overhead)

- Can experiment with different GOF statistics without re-running model

- Can control number of posterior predictive simulations independently

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Run model without GOF during fitting
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2, gof = FALSE,
           nscan = 100, burn = 10, odens = 1, print = FALSE)

# Compute GOF post-hoc
gof_result <- gof(fit)
#> ℹ Computing GOF statistics post-hoc
#> ℹ Using 99 posterior predictive simulations
#> ✔ GOF computation complete
# }
```
