# Exact variance-weighted sum-to-zero projection for the per-actor block

Given an unconstrained sample `theta_star` (n_actors x T) and per-actor
per-period posterior variances `V` (same shape), projects to the
sum-to-zero manifold per period using the variance-weighted Lagrangian
projection (exact under block-diagonal covariance across actors).

## Usage

``` r
.exact_center_per_actor(theta_star, V)
```
