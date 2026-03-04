# Update dynamic latent positions using AR(1) process

Update dynamic latent positions using AR(1) process

## Usage

``` r
rUV_dynamic_fc_cpp(
  U_current,
  V_current,
  ET,
  rho_uv,
  sigma_uv,
  s2,
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

  AR(1) autoregressive parameter for latent positions

- sigma_uv:

  Innovation standard deviation for latent positions

- s2:

  Dyadic variance

- shrink:

  Whether to apply shrinkage

- symmetric:

  Whether network is symmetric

## Value

List with updated U and V arrays
