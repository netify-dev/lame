# Sample AR(1) parameter for dynamic latent factors

Uses a conjugate Normal(prior_mean, prior_sd^2) prior on rho. Defaults
(prior_mean = 0, prior_sd = 1) preserve the historical N(0,1) behaviour
for backward compatibility; lame::lame() passes the user-set
prior\$rho_uv_mean / prior\$rho_uv_sd explicitly.

## Usage

``` r
sample_rho_uv(
  U_cube,
  V_cube,
  sigma_uv,
  rho_current,
  symmetric,
  prior_mean = 0,
  prior_sd = 1
)
```

## Arguments

- U_cube:

  3D array of U positions (n x R x T)

- V_cube:

  3D array of V positions (n x R x T)

- sigma_uv:

  Innovation standard deviation

- rho_current:

  Current value of rho (ignored; full conditional)

- symmetric:

  Whether network is symmetric

- prior_mean:

  Prior mean for rho (default 0)

- prior_sd:

  Prior SD for rho (default 1)

## Value

Updated rho value
