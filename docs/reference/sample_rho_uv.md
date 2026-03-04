# Sample AR(1) parameter for dynamic latent factors

Sample AR(1) parameter for dynamic latent factors

## Usage

``` r
sample_rho_uv(U_cube, V_cube, sigma_uv, rho_current, symmetric)
```

## Arguments

- U_cube:

  3D array of U positions (n x R x T)

- V_cube:

  3D array of V positions (n x R x T)

- sigma_uv:

  Innovation standard deviation

- rho_current:

  Current value of rho

- symmetric:

  Whether network is symmetric

## Value

Updated rho value
