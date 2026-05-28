# Sample the static-block beta conditional on the dynamic path

Conjugate Gaussian update for the static coefficient block, treating the
dynamic path as known. Uses a flat-ish ridge-style prior `prior_prec` on
the static beta (typically diag(1/g) to match the existing g-prior
path).

## Usage

``` r
sample_beta_static_cpp(
  Xdyn_list,
  Xstat_list,
  Z_list,
  offset_list,
  beta_dyn_path,
  prior_mean,
  prior_prec,
  s2,
  dyad_rho,
  bipartite,
  symmetric,
  use_dyad_rho
)
```

## Arguments

- Xdyn_list:

  T-length list of long-format dynamic design matrices.

- Xstat_list:

  T-length list of long-format static design matrices.

- Z_list:

  T-length list of latent (n x n) matrices.

- offset_list:

  T-length list of (n x n) offsets.

- beta_dyn_path:

  (T x p_dyn) dynamic path matrix.

- prior_mean:

  Length p_static prior mean.

- prior_prec:

  p_static x p_static prior precision.

- s2:

  Dyadic variance.

- dyad_rho:

  Dyadic correlation.

- bipartite:

  TRUE/FALSE.

- symmetric:

  TRUE/FALSE.

- use_dyad_rho:

  TRUE/FALSE.

## Value

List with: beta – length p_static; chol_fail – integer.
