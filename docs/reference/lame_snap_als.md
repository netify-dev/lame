# Fast approximate dynamic snap-shift AME estimator

Fits a fast point-estimator approximation to the dynamic snap-shift
latent-factor model for longitudinal normal-valued networks. It targets
the same drift-versus-reset transition estimand as
[`lame`](https://netify-dev.github.io/lame/reference/lame.md) with
`dynamic_uv = TRUE` and `dynamic_uv_kind = "snap"`, but returns ALS snap
scores rather than MCMC draws.

## Usage

``` r
lame_snap_als(
  Y,
  Xdyad = NULL,
  Xrow = NULL,
  Xcol = NULL,
  R = 2L,
  R_row = NULL,
  R_col = NULL,
  family = "normal",
  mode = c("unipartite", "bipartite"),
  symmetric = FALSE,
  max_iter = 200L,
  tol = 1e-06,
  snap_kappa = 2,
  snap_pi_prior = c(a = 1, b = 9),
  snap_update = c("soft", "hard", "annealed"),
  snap_damping = 0.7,
  estimate_rho_uv = TRUE,
  estimate_sigma_uv = TRUE,
  rho_uv = NULL,
  sigma_uv = NULL,
  hyper_update = c("robust", "em"),
  drift_quantile = 0.05,
  drift_min_transitions = 50,
  rho_prior_mean = 0.98,
  rho_prior_weight = 25,
  align = c("sequential", "global", "none"),
  threshold = 0.5,
  min_sigma = 1e-04,
  sigma_floor_fraction = 0.75,
  ridge = 1e-08,
  snap_stability_tol = 0.05,
  snap_convergence = c("quantile", "max", "classification"),
  snap_delta_quantile = 0.95,
  snap_class_change_tol = 0.005,
  unstable_top_n = 10L,
  stability = c("none", "quick", "validation"),
  verbose = TRUE,
  seed = 6886
)
```

## Arguments

- Y:

  a list of relational matrices, or a 3D array. Unipartite fits use
  `[n, n, T]` panels. Bipartite fits use rectangular `[n_row, n_col, T]`
  panels. For named lists, slices may have changing actor composition;
  the estimator pads to the union actor set and treats absent
  actor-periods as unobserved and snap-ineligible.

- Xdyad:

  optional list of dyadic covariate matrices/arrays, or `NULL`.

- Xrow, Xcol:

  not supported by this fast snap-shift estimator; pass `NULL`. Use
  [`lame`](https://netify-dev.github.io/lame/reference/lame.md) with
  `method = "mcmc"` for node-covariate dynamic snap models.

- R:

  positive integer latent rank. For bipartite fits, used as the default
  for `R_row` and `R_col`.

- R_row, R_col:

  positive integer latent ranks for bipartite row and column positions.
  Defaults to `R`.

- family:

  currently only `"normal"`.

- mode:

  `"unipartite"` or `"bipartite"`.

- symmetric:

  logical; if `TRUE`, fit a symmetric latent-factor approximation and
  report one snap-probability matrix.

- max_iter:

  maximum block-coordinate iterations.

- tol:

  convergence tolerance on the approximate penalized objective.

- snap_kappa:

  diffuse snap-prior standard deviation.

- snap_pi_prior:

  length-two vector `c(a, b)` for the beta prior used to regularize the
  snap rate.

- snap_update:

  one of `"soft"`, `"hard"`, or `"annealed"`.

- snap_damping:

  scalar in `(0, 1]`; damping applied to soft snap-score updates. Values
  below 1 mix the new profiled score with the previous score to reduce
  oscillation on large panels. Ignored for `snap_update = "hard"`.

- estimate_rho_uv, estimate_sigma_uv:

  logical flags for updating the AR drift persistence and innovation
  scale.

- rho_uv, sigma_uv:

  optional fixed/initial AR drift parameters.

- hyper_update:

  one of `"robust"` or `"em"`. The default `"robust"` estimates drift
  hyperparameters from the lower tail of transition innovations, which
  prevents broad ruptures from being absorbed into an overly diffuse
  drift process. `"em"` uses the untrimmed soft-classification moment
  update.

- drift_quantile:

  lower-tail transition quantile used by `hyper_update = "robust"` for
  estimating smooth-drift hyperparameters.

- drift_min_transitions:

  minimum effective number of transitions retained by the lower-tail
  update. This keeps the default from overfitting the smooth-drift scale
  on small panels.

- rho_prior_mean, rho_prior_weight:

  weak regularization for the persistence estimate. The defaults encode
  the snap-shift model's intended persistent, low-innovation drift
  baseline.

- align:

  initialization alignment mode; `"global"` aligns each period to the
  pooled static ALS fit, `"sequential"` aligns each period to the
  previous initialized period, and `"none"` leaves per-period factors
  unaligned.

- threshold:

  hard-classification threshold for `snap_class`.

- min_sigma:

  lower bound for the drift scale.

- sigma_floor_fraction:

  for robust hyperparameter updates, the iterative drift scale cannot
  fall below this fraction of the initial `sigma_uv`. This prevents
  broad snap assignments from collapsing the smooth-drift variance to
  `min_sigma`.

- ridge:

  small ridge added to latent-position normal equations.

- snap_stability_tol:

  tolerance for declaring soft snap scores stable. This is separate from
  `tol`, which tracks the penalized objective.

- snap_convergence:

  convergence criterion for snap scores. `"quantile"` uses
  `snap_delta_quantile` of the absolute score changes plus the
  class-change share; `"max"` uses the worst actor-period score change;
  `"classification"` uses only class-change stability. The worst-case
  max is always stored and warned about when high.

- snap_delta_quantile:

  quantile of actor-period score changes used by
  `snap_convergence = "quantile"`.

- snap_class_change_tol:

  maximum share of eligible actor-periods whose `snap_class` may change
  between iterations while still declaring snap convergence.

- unstable_top_n:

  number of largest actor-period snap-score changes to store in
  convergence diagnostics.

- stability:

  optional start-sensitivity preset. `"none"` runs one fit. `"quick"`
  and `"validation"` rerun the same estimator from additional seeds and
  attach `fit$stability` with snap-score, classification, top-period,
  and fitted-surface comparisons.

- verbose:

  logical; print progress.

- seed:

  integer seed for initialization perturbations.

## Value

An object of class `"lame_snap_als"` with dynamic latent positions, snap
scores, fitted values, residuals, and convergence diagnostics.

## Details

On longer panels, the default `align = "sequential"` keeps the
period-to-period movement in view. `align = "global"` can help on very
short panels, but on longer panels it can pull each period back toward
the pooled static fit and make early jumps look too strong. If you have
a fixed drift scale, pass `rho_uv` and `sigma_uv` with
`estimate_rho_uv = FALSE` and `estimate_sigma_uv = FALSE`.

The convergence output separates the usual score summary from the worst
moving actor-period. By default, convergence uses the 95th percentile of
score changes and the share of class changes. `final_max_snap_delta` and
`unstable_transitions` still show the largest local moves. Read
`snap_prob` as an ALS snap score. It is useful for rankings and
heuristic classifications; it is not a Bayesian posterior probability.

In bipartite mode, the fitted multiplicative term is \\U_t G V_t'\\ with
one static interaction matrix \\G\\. The returned `snap_prob` matrix
contains row-actor snap scores and `snap_prob_v` contains column-actor
snap scores. Bipartite snap ALS does not estimate node covariates,
dynamic coefficients, or a dynamic \\G_t\\; those combinations need a
separate model.

## Examples

``` r
set.seed(1)
Y_bip <- lapply(seq_len(3), function(t) {
  m <- matrix(rnorm(6 * 5), 6, 5)
  rownames(m) <- paste0("r", seq_len(6))
  colnames(m) <- paste0("c", seq_len(5))
  m
})
fit_bip <- lame_snap_als(Y_bip, R = 1, mode = "bipartite",
                         max_iter = 3, verbose = FALSE)
dim(fit_bip$snap_prob_v)
#> [1] 5 3
```
