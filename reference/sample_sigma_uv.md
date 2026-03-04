# Sample innovation variance for dynamic latent factors

Sample innovation variance for dynamic latent factors

## Usage

``` r
sample_sigma_uv(U_cube, V_cube, rho_uv, symmetric)
```

## Arguments

- U_cube:

  3D array of U positions (n x R x T)

- V_cube:

  3D array of V positions (n x R x T)

- rho_uv:

  AR(1) parameter

- symmetric:

  Whether network is symmetric

## Value

Updated sigma_uv value
