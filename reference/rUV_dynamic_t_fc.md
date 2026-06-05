# Gibbs sampling of dynamic U and V with heavy-tailed (Student-t) innovations

Like `rUV_dynamic_fc` but the AR(1) innovations are Student-t via a
scale-mixture of normals, a continuous heavy-tailed alternative to
snap-shift.

## Usage

``` r
rUV_dynamic_t_fc(
  U,
  V,
  ET,
  rho_uv,
  sigma_uv,
  s2,
  nu,
  lambda_u = NULL,
  lambda_v = NULL,
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

  AR(1) autoregressive parameter

- sigma_uv:

  Innovation scale

- s2:

  dyadic variance

- nu:

  Student-t degrees of freedom

- lambda_u:

  current local scales for U (n x T)

- lambda_v:

  current local scales for V (n x T)

- shrink:

  whether to apply shrinkage (default TRUE)

- symmetric:

  whether the network is symmetric (default FALSE)

## Value

list with updated U, V arrays and lambda_u, lambda_v local scales

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
