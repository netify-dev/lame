# Update dynamic latent positions with snap-shift model selection

Like rUV_dynamic_fc_cpp but, for t \> 0, chooses per actor between an
AR(1) drift prior and a diffuse N(0, kappa^2 I) snap prior via a
Gaussian log-marginal-likelihood model selection, drawing a Bernoulli
snap indicator delta and sampling the latent position from the selected
posterior.

## Usage

``` r
rUV_dynamic_snap_fc_cpp(
  U_current,
  V_current,
  ET,
  rho_uv,
  sigma_uv,
  s2,
  kappa,
  pi_snap,
  delta_u_current,
  delta_v_current,
  shrink,
  symmetric
)
```

## Arguments

- U_current:

  Current 3D array of U positions (n x R x T)

- V_current:

  Current 3D array of V positions (n x R x T)

- ET:

  3D array of residuals (n x n x T)

- rho_uv:

  AR(1) autoregressive parameter for the drift prior

- sigma_uv:

  Innovation standard deviation for the drift prior

- s2:

  Dyadic variance

- kappa:

  Diffuse snap-prior standard deviation (kappa^2 \>\> sigma_uv^2)

- pi_snap:

  Prior snap probability

- delta_u_current:

  Current sender-side snap indicators from the previous sweep

- delta_v_current:

  Current receiver-side snap indicators from the previous sweep

- shrink:

  Whether to apply shrinkage

- symmetric:

  Whether network is symmetric

## Value

List with updated U, V arrays and delta_u, delta_v snap indicators
