# Update dynamic latent positions with heavy-tailed (Student-t) AR(1) innovations

Like rUV_dynamic_fc_cpp but each AR(1) innovation is Student-t rather
than Gaussian, via a scale-mixture: the innovation for u_t,i has
variance sigma^2 / lambda_t,i with lambda_t,i ~ Gamma(nu/2, nu/2).
Provides a continuous heavy-tailed alternative to the discrete
snap-shift model.

## Usage

``` r
rUV_dynamic_t_fc_cpp(
  U_current,
  V_current,
  ET,
  rho_uv,
  sigma_uv,
  s2,
  nu,
  lambda_u,
  lambda_v,
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

  AR(1) autoregressive parameter

- sigma_uv:

  Innovation scale

- s2:

  Dyadic variance

- nu:

  Student-t degrees of freedom

- lambda_u:

  Current local scales for U (n x T)

- lambda_v:

  Current local scales for V (n x T)

- shrink:

  Whether to apply shrinkage

- symmetric:

  Whether network is symmetric

## Value

List with updated U, V arrays and lambda_u, lambda_v local scales
