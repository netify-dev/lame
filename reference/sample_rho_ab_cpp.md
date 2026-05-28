# Sample AR(1) parameter for dynamic additive effects

Sample AR(1) parameter for dynamic additive effects

## Usage

``` r
sample_rho_ab_cpp(
  a_mat,
  b_mat,
  sigma_ab,
  rho_current,
  symmetric,
  prior_mean = 0,
  prior_sd = -1
)
```

## Arguments

- a_mat:

  Matrix of row effects (n x T)

- b_mat:

  Matrix of column effects (n x T)

- sigma_ab:

  Innovation standard deviation

- rho_current:

  Current value of rho

- symmetric:

  Whether the network is symmetric

- prior_mean:

  Prior mean for rho. Used only when `prior_sd > 0`.

- prior_sd:

  Prior SD for rho. `prior_sd < 0` (the default) selects a Jeffreys-like
  flat prior; a positive value switches to a truncated
  Normal(prior_mean, prior_sd^2) prior.

## Value

Updated rho value
