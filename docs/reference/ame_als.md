# Fast (MCMC-free) AME estimation for a cross-sectional network

Fits an additive and multiplicative effects (AME) model to a single
cross-sectional network by **iterative block coordinate descent**,
producing a fast point estimate with no MCMC and no credible intervals.

The estimation algorithm adapts the *iterative block coordinate descent*
estimator of the Social Influence Regression (SIR) model of Hoff &
Minhas (the algorithm implemented in
[`sir::sir_alsfit()`](https://rdrr.io/pkg/sir/man/sir_alsfit.html) and
nicknamed "ALS" in that package) to the AME model. It is a port and
adaptation, not original lame methodology.

## Usage

``` r
ame_als(
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

  a square (unipartite) or rectangular (bipartite) relational matrix.

- Xdyad:

  an `n_row x n_col` matrix or `n_row x n_col x pd` array of dyadic
  covariates, or `NULL`.

- Xrow:

  an `n_row x pr` matrix of row/sender covariates, or `NULL`.

- Xcol:

  an `n_col x pc` matrix of column/receiver covariates, or `NULL`.

- R:

  integer dimension of the multiplicative effects (default `0`). **The
  covariate coefficients are conditional on this choice.** The
  multiplicative term \\u_i'v_j\\ is a flexible high-variance regressor
  that can correlate with the dyadic covariates, so the estimated `beta`
  can shift – and occasionally change sign – as `R` increases. Validate
  against an `R = 0` fit and the MCMC
  [`ame`](https://netify-dev.github.io/lame/reference/ame.md)/[`lame`](https://netify-dev.github.io/lame/reference/lame.md)
  before reporting.

- family:

  one of `"normal"`, `"binary"`, `"poisson"`, `"ordinal"`.

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

  for the non-normal families, `"irls"` (default for `binary`/`poisson`)
  runs iteratively reweighted least squares – a fast approximate GLM AME
  fit with coefficients on the calibrated link scale (Poisson log,
  binary logit/probit). `"transform"` fits a single fixed Gaussian
  working response (`log(y+1)` for Poisson, rank-normal scores for
  binary/ordinal); its coefficients are on an uncalibrated rank scale –
  good for direction/ranking of effects, not for magnitudes. **For
  ordinal data the transform path attenuates \\\beta\\ toward zero
  (classical discretisation bias); prefer `"em"` when \\\beta\\ is an
  effect-size estimand.** `"em"` (for `ordinal`, opt-in; default for
  `tobit`) runs a deterministic EM outer loop that maximises the
  observed ordinal log-likelihood, jointly estimating cutpoints
  (returned on `fit$alpha`, anchored at `alpha[1] = 0`) and the
  regression block; the surrogate is monotone non-decreasing across EM
  iterations. `"irls"` is supported for `binary`/`poisson`; `ordinal`
  accepts `"transform"` (default, back-compat) or `"em"`. A directed
  `R > 0` IRLS fit uses the hybrid low-rank solver internally (the IRLS
  weights are unbalanced). Uncertainty for any of these is the
  bootstrap.

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

An object of class `"ame_als"`: a list with the point estimates `mu`,
`beta`, `a`, `b`, `U`, `V` (`L` for symmetric models), the linear
predictor `EZ`, response-scale `fitted` values, working-scale
`residuals`, convergence information, and the variance-component vector
`VC` with five descriptive entries:

- `va`, `vb`:

  empirical variances of the sender and receiver additive effects.

- `cab`:

  covariance of the sender and receiver effects (`NA` for symmetric or
  bipartite models).

- `rho`:

  dyadic residual reciprocity — the correlation between the residuals of
  \\(i,j)\\ and \\(j,i)\\ — not the sender/receiver correlation
  `cor(a, b)`.

- `ve`:

  residual (working-scale) variance.

These are descriptive summaries of the point estimates, not
random-effect variance components. See
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
for uncertainty and
[`vcov.ame_als`](https://netify-dev.github.io/lame/reference/vcov.ame_als.md)
for a fast analytic covariance of the regression coefficients.

## Details

The AME decomposition \$\$z\_{ij} = \mu + \beta' x\_{ij} + a_i + b_j +
u_i' v_j + \epsilon\_{ij}\$\$ is conditionally linear: it is a linear
regression in \\(\mu, \beta, a, b)\\ for fixed multiplicative factors
\\(U, V)\\, and the optimal rank-`R` \\(U, V)\\ for fixed \\(\mu, \beta,
a, b)\\ is the truncated SVD of the residual matrix. The estimator
therefore cycles, until the residual sum of squares stabilises, through
three blocks — each a monotone (objective-non-increasing) update: a
joint least-squares solve for the intercept \\\mu\\ and regression
coefficients \\\beta\\ (with a
[`ginv`](https://rdrr.io/pkg/MASS/man/ginv.html) pseudoinverse fallback
for rank-deficient designs); Gauss-Seidel sweeps to convergence for the
additive effects \\(a, b)\\; and a weighted low-rank update for the
multiplicative factors \\(U, V)\\ via SVD (eigen-decomposition when
`symmetric = TRUE`). The weighting makes the multiplicative step correct
for unbalanced longitudinal panels (dyads observed at unequal numbers of
time points); severely unbalanced panels may need more iterations to
converge, so raise `max_iter` if the fit reports non-convergence.

For `family = "normal"` this is a least-squares (Gaussian maximum-
likelihood) fit — the exact global solution when `R = 0`, and a local
optimum of the non-convex low-rank objective when `R > 0`. For
`"binary"`, `"poisson"` and `"ordinal"` the same algorithm is run on a
Gaussian working response (the latent transforms
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) uses to
initialise its sampler): a fast *approximation* suitable for exploration
and starting values. For `"tobit"` an EM outer loop replaces censored
cells with their truncated-normal conditional mean and the BCD M-step is
warm-started from the previous iterate (see the “Tobit attenuation”
section below for the documented \\\sigma^2\\ bias and the
moderate-to-heavy censoring regimes where the point estimates of
\\\beta\\ can drift noticeably from the MCMC estimator). For calibrated
family-specific inference, use the MCMC estimator
[`ame`](https://netify-dev.github.io/lame/reference/ame.md). The
remaining censoring/rank families (`"cbin"`, `"frn"`, `"rrl"`) are not
supported and raise an informative error.

## Tobit attenuation

The EM-ALS path for `family = "tobit"` is a deterministic generalised EM
that imputes censored \\Z\_{ij}\\ on the truncated-normal mean
\\\eta\_{ij} -
\sigma\\\phi(-\eta\_{ij}/\sigma)/\Phi(-\eta\_{ij}/\sigma)\\ and refits
the BCD low-rank fit on the imputed working response. The reported
\\\sigma^2\\ carries a downward (attenuation) bias that scales with the
censoring fraction: the additive sender/receiver effects \\a_i, b_j\\
absorb part of the censored signal during the M-step, leaving the
residual scale smaller than the truth. In the package's own internal
simulations (\\n = 80\\, dyadic covariate, true \\\sigma^2 = 1\\) the
reported \\\sigma^2\\ runs roughly 0.2-0.7 at 50\\ 0.1 at 80\\
essentially without bias. The censored-cell EM also drifts the intercept
and the dyadic slope: in the same simulation EM-ALS gives \\(\hat\mu,
\hat\beta\_\mathrm{dyad}) \approx (-1.4, 0.31)\\ at 50\\ censoring and
\\(-0.9, 0.03)\\ at 80\\ \\(0, 0.5)\\ and \\(-1.5, 0.5)\\ respectively),
while MCMC tobit on the same draws gives \\(0.11, 0.47)\\ and \\(-1.28,
0.42)\\. Treat the EM-ALS tobit point estimates as fast exploratory
summaries; use
[`ame`](https://netify-dev.github.io/lame/reference/ame.md)`(..., family = "tobit")`
for calibrated point and interval estimates. The EM iteration count and
\\\sigma^2\\ are surfaced on the fit as `$em_iters` and `$sigma2`.

**Row/column (node) covariates.** A node covariate broadcasts to a
per-actor constant, which is collinear with the additive sender/receiver
effect, so its coefficient is not identified by the objective alone. It
is identified here by an explicit constraint: the additive effects are
taken orthogonal to the node covariates, and `beta_row`/`beta_col` are
the corresponding between-actor regression coefficients (the additive
effects then carry only the residual heterogeneity). This is the
standard estimand under the assumption that the additive effects are
uncorrelated with the node covariates; if that assumption is doubtful
the coefficient absorbs the covariate-correlated part of the additive
heterogeneity. Time-varying node covariates are summarised by their
per-actor mean; if a node covariate varies within actor over time the
discarded within-actor variation triggers a warning.

**Identifiability of the additive and multiplicative terms.** For
`R > 0` the additive term \\a_i + b_j\\ and the multiplicative term
\\u_i' v_j\\ are not separately identified by the objective alone: a
broadcast (row- or column-constant) component can sit in either, since a
pure sender effect \\a 1'\\ is itself rank one. The estimator imposes
the standard AME gauge — the multiplicative term is double-centered
(zero row and column means), so all broadcast structure is carried by
\\a, b\\ — which makes the reported \\a, b, U, V\\ unique given the
fitted values and independent of the optimisation path. For a unipartite
network the double-centering means are taken over the full matrix, which
includes the structurally unobserved self-tie diagonal that the model
fills in by its low-rank completion; the additive/multiplicative split
therefore carries an \\O(1/n)\\ dependence on that completion. On a
disconnected observed-dyad graph the additive effects additionally have
a per-component level shift, which is pinned to a precision-weighted
minimum-norm gauge so that `a`, `b` and the variance components remain
reproducible. A dyadic covariate that is itself (near) low-rank can
still be partially aliased with the multiplicative term, so keep `R`
modest relative to the covariate structure. This residual aliasing is
intrinsic to the AME model — the MCMC estimator resolves it only through
its priors.

**Choosing R.** There is no automatic order-selection criterion (the
working-response objective has no likelihood, so AIC/BIC do not apply).
Fit a few values — e.g.
`lapply(0:4, function(r) ame_als(Y, R = r, ...))` — and inspect
`deviance` (the residual sum of squares): it falls steeply while real
multiplicative signal is being captured and then flattens, so the
“elbow” of that curve is a reasonable choice. With `R > 0` the objective
is non-convex; use `multistart` to guard against local optima.

Uncertainty is obtained separately, by the bootstrap; see
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md).
(The SIR paper's own primary standard errors are Hessian-based,
classical and sandwich/robust. For AME a Hessian-based variance is
awkward on two counts: without an explicit gauge fix the rotational
invariance of the multiplicative factors leaves the joint Hessian
rank-deficient, and even with a gauge fixed the Gaussian
working-response approximation used for the non-normal families leaves a
Hessian-based variance miscalibrated. The bootstrap side-steps both, so
it is preferred here.)

## Limitations vs. [`ame`](https://netify-dev.github.io/lame/reference/ame.md) / [`lame`](https://netify-dev.github.io/lame/reference/lame.md)

The ALS estimator is a fast, frequentist point estimator. It is not a
drop-in replacement for the MCMC estimators – the following features
apply to MCMC only:

- **Families.** ALS supports `normal`, `binary`, `poisson`, `ordinal`
  (transform or EM), and `tobit` (EM only; see the “Tobit attenuation”
  section for the documented \\\sigma^2\\ downward bias and the
  censoring regimes where the point estimates drift from the MCMC
  estimator). It does **not** support `cbin`, `frn`, `rrl`; for those,
  fall back to
  [`ame`](https://netify-dev.github.io/lame/reference/ame.md) /
  [`lame`](https://netify-dev.github.io/lame/reference/lame.md). The
  `bootstrap = N` option is **not** calibrated for `family = "tobit"` –
  the parametric simulator falls through to the binary/probit branch,
  producing degenerate 0/1 resamples and coefficient bootstraps
  collapsed near zero. Use the MCMC tobit for uncertainty
  quantification.

- **Dynamic effects.**
  [`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
  fits a STATIC model pooled across time slices. There is no
  `dynamic_uv`, `dynamic_ab`, `dynamic_G`, or `dynamic_beta`.
  `lame(..., method = "als")` warns and ignores those arguments.

- **Priors.** ALS has no priors. `prior = list(...)` and `g = ...` are
  MCMC-only; the dispatcher `ame(..., method = "als")` warns and ignores
  them.

- **Posterior quantities.** No `$BETA` / `$VC` posterior draws on the
  point fit. Use `bootstrap = N` for a sampling distribution analogue;
  the draws are bootstrap, not Bayesian.

- **Multi-chain.** Not applicable; `n_chains` is dropped.

- **Convergence diagnostics.** No Rhat / ESS /
  [`trace_plot`](https://netify-dev.github.io/lame/reference/trace_plot.md).

Features that work the same on ALS fits: `coef`, `vcov` (sandwich on
regression block), `confint` (auto-routes bootstrap intervals when
present, sandwich Wald otherwise), `predict`, `fitted`, `residuals`,
`summary`, `print`, `nobs`,
[`simulate.ame_als`](https://netify-dev.github.io/lame/reference/simulate.ame_als.md),
[`gof_plot.ame_als`](https://netify-dev.github.io/lame/reference/gof_plot.ame_als.md),
[`ab_plot.ame_als`](https://netify-dev.github.io/lame/reference/ab_plot.ame_als.md),
[`uv_plot`](https://netify-dev.github.io/lame/reference/uv_plot.md),
[`latent_positions`](https://netify-dev.github.io/lame/reference/latent_positions.md).

## References

Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
Influence Regression. *Political Analysis*. The iterative block
coordinate descent estimator adapted here originates with that work
(implemented in
[`sir::sir_alsfit()`](https://rdrr.io/pkg/sir/man/sir_alsfit.html)).

## See also

[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
for longitudinal networks,
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
for bootstrap uncertainty,
[`ame`](https://netify-dev.github.io/lame/reference/ame.md) for the full
MCMC estimator.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
coef(fit)
#>  intercept 
#> 0.05534878 
```
