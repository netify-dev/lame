# Sample the AR(1) innovation sigma for each dynamic block

Inverse-Gamma full conditional, one sigma per block.

## Usage

``` r
sample_sigma_beta_cpp(
  beta_path,
  group_id,
  n_groups,
  Lambda_inv,
  rho_by_group,
  prior_shape,
  prior_scale
)
```

## Arguments

- beta_path:

  (T x p_dyn) matrix.

- group_id:

  Length p_dyn integer vector (1-based) of block IDs.

- n_groups:

  Number of distinct block IDs.

- Lambda_inv:

  (p_dyn x p_dyn) inverse of the (full) Lambda scale.

- rho_by_group:

  Length n_groups vector of current rho values.

- prior_shape:

  Length n_groups vector of IG shape parameters.

- prior_scale:

  Length n_groups vector of IG scale parameters.

## Value

Length n_groups vector of new sigma values.
