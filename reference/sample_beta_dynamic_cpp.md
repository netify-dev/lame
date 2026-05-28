# Sample the dynamic-block beta path via FFBS

Forward-filter / backward-sample the AR(1) state-space model for the
dynamic-block beta coefficients. Returns the joint draw of beta_dyn at
every time period.

## Usage

``` r
sample_beta_dynamic_cpp(
  Xdyn_list,
  Xstat_list,
  Z_list,
  offset_list,
  beta_static,
  rho_by_coef,
  sigma_by_coef,
  Lambda,
  beta0_mean,
  beta0_cov,
  s2,
  dyad_rho,
  bipartite,
  symmetric,
  use_dyad_rho
)
```

## Arguments

- Xdyn_list:

  T-length list of (n\*n) x p_dyn long-format design matrices for the
  dynamic block (column-major reshape per period).

- Xstat_list:

  T-length list of (n\*n) x p_static long-format design matrices for the
  static block.

- Z_list:

  T-length list of (n x n) latent matrices.

- offset_list:

  T-length list of (n x n) offset matrices (a_i + b_j + U_i'V_j
  contributions; everything that's not in X\*beta).

- beta_static:

  Length p_static current static beta vector.

- rho_by_coef:

  Length p_dyn vector of AR(1) rho values for each dynamic coefficient.
  (Per-block but expanded per-column for vectorised indexing.)

- sigma_by_coef:

  Length p_dyn vector of AR(1) innovation standard deviations for each
  dynamic coefficient.

- Lambda:

  Block-diagonal innovation scale matrix (p_dyn x p_dyn). Combined with
  sigma^2 to give Q = sigma^2 \* Lambda.

- beta0_mean:

  Length p_dyn prior mean for beta at t=0.

- beta0_cov:

  p_dyn x p_dyn prior covariance for beta at t=0.

- s2:

  Dyadic variance.

- dyad_rho:

  Dyadic correlation (ignored when use_dyad_rho=FALSE).

- bipartite:

  Whether the network is bipartite.

- symmetric:

  Whether the network is symmetric.

- use_dyad_rho:

  Whether to use the dyad-corr branch (TRUE only for unipartite,
  asymmetric, with a non-zero rho).

## Value

List with: path – a (T x p_dyn) matrix of beta draws (one row per
period); chol_fail – integer count of Cholesky failures.
