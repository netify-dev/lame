# Compute EZ when beta is time-varying

Rebuild the (n x n x T) or (nA x nB x T) EZ cube using a per-period beta
vector. Mirrors get_EZ_cpp / get_EZ_bip_cpp but accepts a (T x p) beta
matrix where each row is the beta for that period.

## Usage

``` r
get_EZ_dynamic_beta_cpp(
  Xlist,
  beta_full_path,
  a_mat,
  b_mat,
  U_cube,
  V_cube,
  G,
  bipartite,
  symmetric
)
```

## Arguments

- Xlist:

  T-length list of design arrays (each n x n x p or nA x nB x p).

- beta_full_path:

  (T x p) matrix of per-period beta. For coefficients that are NOT in
  the dynamic block, the rows are identical (the static beta replicated
  across all periods).

- a_mat:

  (n_a x T) row effects (or repmat-ed static effects).

- b_mat:

  (n_b x T) column effects.

- U_cube:

  (n_u x R x T) row latent positions (or replicated-static).

- V_cube:

  (n_v x R x T) column latent positions.

- G:

  (R x R or RA x RB) interaction matrix (bipartite); identity for
  unipartite.

- bipartite:

  TRUE/FALSE.

- symmetric:

  TRUE/FALSE.

## Value

n_a x n_b x T cube of EZ values.
