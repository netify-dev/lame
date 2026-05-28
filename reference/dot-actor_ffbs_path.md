# Univariate Carter-Kohn FFBS for one actor's length-T slope path

Given period-wise sufficient statistics
(H[t](https://rdrr.io/r/base/t.html) = sum_j X^2_ij / s2,
h[t](https://rdrr.io/r/base/t.html) = sum_j X_ij R_ij / s2) and AR(1)
hyperparameters, returns one joint draw from the posterior N(m, P) where
the prior is AR(1) and the observation is the Gaussian likelihood
implied by (H, h).

## Usage

``` r
.actor_ffbs_path(H, h, rho_actor, sigma_actor2)
```

## Arguments

- H:

  length-T vector of period-wise observation precisions

- h:

  length-T vector of period-wise observation cross-products

- rho_actor:

  AR(1) coefficient

- sigma_actor2:

  AR(1) innovation variance

## Value

list with `theta` (length-T draw) and `V` (length-T marginal posterior
variance per period – needed for the exact centering projection).
