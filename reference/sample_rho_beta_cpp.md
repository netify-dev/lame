# Sample the AR(1) rho for each dynamic block

Truncated-Normal full conditional, one rho per block.

## Usage

``` r
sample_rho_beta_cpp(
  beta_path,
  group_id,
  n_groups,
  Lambda_inv,
  sigma_by_coef,
  rho_current,
  rho_prior_mean,
  rho_prior_sd,
  rho_lower = 0,
  rho_upper = 0.999
)
```

## Arguments

- beta_path:

  (T x p_dyn) matrix of beta draws.

- group_id:

  Length p_dyn integer vector (1-based) of block IDs.

- n_groups:

  Number of distinct block IDs.

- Lambda_inv:

  (p_dyn x p_dyn) inverse of the (full) Lambda scale.

- sigma_by_coef:

  Length p_dyn vector of per-coef sigma (block-shared).

- rho_current:

  Length n_groups vector of current rho values.

- rho_prior_mean:

  Length n_groups vector of prior means.

- rho_prior_sd:

  Length n_groups vector of prior SDs.

- rho_lower:

  Lower truncation bound (typically 0).

- rho_upper:

  Upper truncation bound (typically 0.999).

## Value

Length n_groups vector of new rho values.
