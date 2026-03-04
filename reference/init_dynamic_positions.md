# Initialize dynamic latent positions with AR(1) structure

Initialize dynamic latent positions with AR(1) structure

## Usage

``` r
init_dynamic_positions(n, R, T, rho_uv, sigma_uv)
```

## Arguments

- n:

  Number of actors

- R:

  Latent dimension

- T:

  Number of time points

- rho_uv:

  AR(1) parameter

- sigma_uv:

  Innovation standard deviation

## Value

3D array of latent positions (n x R x T)
