# Fast (MCMC-free) AME estimation for a longitudinal network

Fits an additive and multiplicative effects (AME) model to a
longitudinal (replicated) network by **iterative block coordinate
descent**, producing a fast point estimate with no MCMC and no credible
intervals. This is the longitudinal counterpart of
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md); the
estimated effects (`mu`, `beta`, `a`, `b`, `U`, `V`) are *static*
(pooled across time), as in a non-dynamic
[`lame`](https://netify-dev.github.io/lame/reference/lame.md) fit.

The estimation algorithm adapts the iterative block coordinate descent
estimator of the Social Influence Regression model of Hoff & Minhas
(`sir::sir_alsfit()`) to the AME model; it is a port and adaptation, not
original lame methodology. See
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) for
the algorithm details.

## Usage

``` r
lame_als(
  Y,
  Xdyad = NULL,
  Xrow = NULL,
  Xcol = NULL,
  R = 0,
  family = "normal",
  mode = c("unipartite", "bipartite"),
  symmetric = FALSE,
  max_iter = 200,
  tol = 1e-06,
  lowrank_method = c("mm", "als", "hybrid"),
  non_normal_method = c("irls", "transform"),
  link = c("probit", "logit"),
  linear_solver = c("eigen", "qr", "auto"),
  multistart = c("none", "cheap", "full"),
  bootstrap = 0L,
  bootstrap_type = c("parametric", "block"),
  bootstrap_block_length = 1L,
  bootstrap_seed = NULL,
  verbose = TRUE,
  seed = 6886
)
```

## Arguments

- Y:

  a list of `T` relational matrices, or a 3D array `[n_row, n_col, T]`.
  A named list may have changing actor composition; every slice must
  carry row and column names so actors can be aligned to the union
  panel. Unnamed lists and arrays are treated positionally and must have
  one fixed layout.

- Xdyad:

  a list of `T` dyadic covariate matrices/arrays, or `NULL`.

- Xrow:

  a list of `T` row/sender covariate matrices, or `NULL`.

- Xcol:

  a list of `T` column/receiver covariate matrices, or `NULL`.

- R:

  integer dimension of the multiplicative effects (default `0`). **The
  covariate coefficients are conditional on this choice.** The
  multiplicative term \\u_i'v_j\\ is a flexible high-variance regressor
  that can correlate with the dyadic covariates, so the estimated `beta`
  can shift – and occasionally change sign – as `R` increases. Comparing
  against an `R = 0` fit is a useful check on whether the covariate
  story is being driven by the latent rank.

- family:

  one of `"normal"`, `"binary"`, or `"poisson"`. The rank and censoring
  families are MCMC-only.

- mode:

  `"unipartite"` (square) or `"bipartite"` (rectangular).

- symmetric:

  logical; fit a symmetric (undirected) model. Unipartite only.

- max_iter:

  maximum number of block coordinate descent iterations (default 200).

- tol:

  convergence tolerance on the relative change in residual sum of
  squares (default `1e-6`).

- lowrank_method:

  inner solver for the multiplicative (low-rank) block: `"mm"` (default)
  weighted majorise-minimise; `"als"` alternating least squares;
  `"hybrid"` runs both and keeps the lower-objective result.
  `"als"`/`"hybrid"` converge faster than `"mm"` on strongly unbalanced
  longitudinal panels and are available for directed and bipartite
  models (symmetric fits always use `"mm"`). All three minimise the same
  objective, so the point estimate is unchanged for balanced data.

- non_normal_method:

  for the non-normal ALS families, `"irls"` (default for
  `binary`/`poisson`) runs iteratively reweighted least squares, giving
  a fast approximate GLM AME fit with coefficients on the requested link
  scale (Poisson log, binary logit/probit). `"transform"` fits one fixed
  Gaussian working response (`log(y+1)` for Poisson, rank-normal scores
  for binary); its coefficients are on an uncalibrated working scale and
  are mainly useful for direction/ranking checks. A directed `R > 0`
  IRLS fit uses the hybrid low-rank solver internally because the IRLS
  weights are unbalanced. Uncertainty for either path comes from the
  bootstrap or sandwich covariance.

- link:

  link for `non_normal_method = "irls"` with a `binary` family:
  `"probit"` (default; matches
  [`ame`](https://netify-dev.github.io/lame/reference/ame.md)/[`lame`](https://netify-dev.github.io/lame/reference/lame.md))
  or `"logit"`. `poisson` always uses the log link; ignored otherwise.

- linear_solver:

  solver for the regression block: `"eigen"` (default) eigendecomposes
  the normal equations; `"qr"` uses a QR factorisation of the observed
  design, which is more stable for ill-conditioned covariates; `"auto"`
  picks `"qr"` when the design is ill-conditioned but full rank. All
  give the same answer for well-conditioned designs.

- multistart:

  for `R > 0` (a non-convex objective), `"none"` (default) fits from a
  single deterministic start; `"cheap"` (4 starts) and `"full"` (8
  starts) also try random low-rank starts and keep the lowest-SSE fit,
  warning when the starts reach materially different optima.
  Reproducible given `seed`; the global RNG stream is left unchanged.

- bootstrap:

  integer: if \> 0, additionally run `bootstrap` replicates of the
  parametric or block bootstrap (via
  [`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md))
  after the point fit and attach the result as `fit$bootstrap`. The
  downstream accessors
  ([`confint.ame_als`](https://netify-dev.github.io/lame/reference/confint.ame_als.md),
  `summary`, `print`) then surface bootstrap intervals instead of the
  anti-conservative sandwich Wald intervals. Default `0` (no bootstrap)
  – bootstrap is expensive, so it is not imposed on a user who just
  wants a quick fit.

- bootstrap_type:

  character: `"parametric"` (default) or `"block"` – the bootstrap
  scheme to use when `bootstrap > 0`.

- bootstrap_block_length:

  integer: block length for the block bootstrap; only used when
  `bootstrap_type = "block"`.

- bootstrap_seed:

  optional integer seed for the bootstrap (the point fit uses `seed`).

- verbose:

  logical; print progress (default `TRUE`).

- seed:

  random seed (default `6886`). The block coordinate descent is
  deterministic, so the point estimate is reproducible regardless; the
  argument is retained for API consistency with
  [`ame`](https://netify-dev.github.io/lame/reference/ame.md).

## Value

An object of class `"ame_als"`; see
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md).

## References

Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
Influence Regression. *Political Analysis*. The iterative block
coordinate descent estimator adapted here originates with that work
(implemented in `sir::sir_alsfit()`).

## See also

[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md),
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md),
[`lame`](https://netify-dev.github.io/lame/reference/lame.md) for the
full MCMC estimator,
[`lame_snap_als`](https://netify-dev.github.io/lame/reference/lame_snap_als.md)
for the approximate dynamic snap-shift point estimator,
[`als_dynamic_beta`](https://netify-dev.github.io/lame/reference/als_dynamic_beta.md)
for a regression-only smoother that estimates *only* a time-varying
\\\beta_t\\ (no \\a, b, U, V\\; not a special case of this function).

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
Y <- replicate(4, { m <- matrix(rnorm(400), 20, 20); diag(m) <- NA; m },
               simplify = FALSE)
fit <- lame_als(Y, R = 1, family = "normal", verbose = FALSE)
coef(fit)
#>  intercept 
#> 0.03266262 
```
