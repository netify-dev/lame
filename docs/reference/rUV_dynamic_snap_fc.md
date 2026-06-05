# Gibbs sampling of dynamic U and V with snap-shift dynamics

Like `rUV_dynamic_fc` but lets each actor's latent position either drift
under the AR(1) prior or jump under a diffuse N(0, kappa^2) snap prior,
chosen per actor-period by a Gaussian log-marginal model selection.

## Usage

``` r
rUV_dynamic_snap_fc(
  U,
  V,
  ET,
  rho_uv,
  sigma_uv,
  s2,
  kappa,
  pi_snap,
  delta_u_current = NULL,
  delta_v_current = NULL,
  shrink = TRUE,
  symmetric = FALSE
)
```

## Arguments

- U:

  3D array of current U positions (n x R x T)

- V:

  3D array of current V positions (n x R x T)

- ET:

  3D array of residuals (n x n x T)

- rho_uv:

  AR(1) autoregressive parameter for the drift prior

- sigma_uv:

  Innovation standard deviation for the drift prior

- s2:

  dyadic variance

- kappa:

  diffuse snap-prior standard deviation

- pi_snap:

  prior snap probability

- delta_u_current:

  current sender-side snap indicators. If `NULL`, a zero matrix is used.

- delta_v_current:

  current receiver-side snap indicators. If `NULL`, a zero matrix is
  used.

- shrink:

  whether to apply shrinkage (default TRUE)

- symmetric:

  whether the network is symmetric (default FALSE)

## Value

list with updated U, V arrays and delta_u, delta_v snap indicators

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
