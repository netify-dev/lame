# AME model fitting routine

An MCMC routine providing a fit to an additive and multiplicative
effects (AME) regression model to cross-sectional relational data of
various types. This function supports both unipartite (square) and
bipartite (rectangular) networks. For longitudinal networks, use the
`lame` function. Original implementation by Peter Hoff.

## Usage

``` r
ame(
  Y,
  Xdyad = NULL,
  Xrow = NULL,
  Xcol = NULL,
  rvar = !(family == "rrl"),
  cvar = TRUE,
  dcor = !symmetric,
  nvar = TRUE,
  R = 0,
  R_row = NULL,
  R_col = NULL,
  mode = c("unipartite", "bipartite"),
  family = "normal",
  intercept = !is.element(family, c("rrl", "ordinal")),
  symmetric = FALSE,
  odmax = rep(max(apply(Y > 0, 1, sum, na.rm = TRUE)), nrow(Y)),
  prior = list(),
  g = NA,
  seed = 6886,
  nscan = 10000,
  burn = 500,
  odens = 25,
  verbose = TRUE,
  gof = TRUE,
  custom_gof = NULL,
  plot = FALSE,
  start_vals = NULL,
  periodic_save = FALSE,
  out_file = NULL,
  save_interval = 0.25,
  model.name = NULL,
  posterior_opts = NULL,
  n_chains = 1,
  cores = 1,
  use_sparse_matrices = FALSE,
  method = c("mcmc", "als"),
  bootstrap = 0L,
  bootstrap_type = c("parametric", "block"),
  bootstrap_block_length = 1L,
  bootstrap_seed = NULL,
  save_log_lik = FALSE,
  ordinal_cutpoints = c("data_induced", "explicit"),
  print,
  ...
)
```

## Arguments

- Y:

  For unipartite: an n x n square relational matrix. For bipartite: an
  nA x nB rectangular relational matrix where nA is the number of row
  nodes and nB is the number of column nodes. See family below for
  various data types.

- Xdyad:

  For unipartite: an n x n x pd array of dyadic covariates (e.g.
  distance, shared group). For bipartite: an nA x nB x pd array. A 2-D
  matrix (single dyadic covariate) must be wrapped as
  `array(x, dim = c(n, n, 1))`; `Inf`/`NaN` entries are rejected.

- Xrow:

  For unipartite: an n x pr matrix of **sender** (row) covariates (e.g.
  sender's age, group). For bipartite: an nA x pr matrix of row-node
  covariates. A `data.frame` is accepted and coerced to numeric matrix
  internally.

- Xcol:

  For unipartite: an n x pc matrix of **receiver** (column) covariates.
  For bipartite: an nB x pc matrix of column-node covariates. A
  `data.frame` is accepted and coerced internally.

- rvar:

  logical: fit row random effects (asymmetric case)?

- cvar:

  logical: fit column random effects (asymmetric case)?

- dcor:

  logical: fit a dyadic correlation (asymmetric case)? Note: not used
  for bipartite networks.

- nvar:

  logical: fit nodal random effects (symmetric case)?

- R:

  integer: dimension of the multiplicative effects (can be zero). For
  bipartite networks, this is used as the default for both R_row and
  R_col if they are not specified.

- R_row:

  integer: for bipartite networks, dimension of row node multiplicative
  effects (defaults to R)

- R_col:

  integer: for bipartite networks, dimension of column node
  multiplicative effects (defaults to R)

- mode:

  character: either "unipartite" (default) for square networks or
  "bipartite" for rectangular networks. Not all combinations of `family`
  and `mode` are supported – see the table under **Supported family x
  mode combinations** below.

- family:

  character: one of
  "normal","tobit","binary","ordinal","cbin","frn","rrl","poisson". See
  the **Supported family x mode combinations** table below for which
  combinations are valid; see **Details** for the model assumptions
  behind each family.

- intercept:

  logical: fit model with an intercept?

- symmetric:

  logical: Is the sociomatrix symmetric by design?

- odmax:

  a scalar integer or vector of length n giving the maximum number of
  nominations that each node may make - used for "frn" and "cbin"
  families

- prior:

  a list containing hyperparameters for the prior distributions.
  Available options and their defaults:

  Sab0

  :   Prior covariance matrix for additive effects (default: diag(2)). A
      2x2 matrix where Sab0\\1,1\\ is the prior variance for row
      effects, Sab0\\2,2\\ is the prior variance for column effects, and
      off-diagonals control correlation between row and column effects.

  eta0

  :   Prior degrees of freedom for covariance of multiplicative effects
      (default: 4 + 3 \\ n/100, where n is the number of actors). Higher
      values impose stronger shrinkage on the latent factors. Common
      values: 4-10 for weak shrinkage, 10-20 for moderate, \>20 for
      strong shrinkage.

  etaab

  :   Prior degrees of freedom for covariance of additive effects
      (default: 4 + 3 \\ n/100). Controls shrinkage of row/column random
      effects. Larger values shrink effects toward zero. *Used by the
      bipartite and longitudinal paths; ignored for a unipartite
      cross-sectional `ame()` fit, which controls the additive prior
      through `Sab0` and `eta0`.*

  s20

  :   Prior variance for regression coefficients (default: 1). Larger
      values allow for larger coefficient values. *Used by the bipartite
      and longitudinal paths; for a unipartite `ame()` fit the
      regression prior is controlled by `g`.*

  s2u0

  :   Prior variance for multiplicative effects (default: 1). *Used by
      the bipartite and longitudinal paths.*

  Suv0

  :   Prior covariance for multiplicative effects (default: identity
      matrix).

  Common usage: prior = list(Sab0 = diag(c(2, 2)), eta0 = 10) for
  moderate shrinkage, or prior = list(Sab0 = diag(c(0.5, 0.5))) for
  tighter control. For a unipartite cross-sectional `ame()` fit the
  prior is controlled by `Sab0`, `eta0`, `Suv0` and `g`.

- g:

  optional **scalar** for the Zellner g-prior on regression coefficients
  (`beta ~ N(0, g * sigma^2 * solve(XtX))` where `XtX` is the design
  cross-product). If not specified, defaults are: for `normal` family,
  `g = n * var(Y)`; for `tobit`, `g = n * var(Y) * 4`; for other
  families, `g = n` (number of non-missing dyads). Per-coefficient
  (vector) `g` is NOT currently supported by the unipartite path – pass
  a scalar. *Note:* `g` is a top-level argument to `ame()`, **not** an
  element of `prior = list(...)`; passing `prior = list(g = 0.1)` is a
  no-op (warned about).

- seed:

  random seed

- nscan:

  number of iterations of the Markov chain (beyond burn-in)

- burn:

  burn in for the Markov chain. Typical use is `burn` far smaller than
  `nscan`; `burn > nscan` is allowed but warns.

- odens:

  output density (thinning interval) for the Markov chain.
  `nscan / odens` samples are stored.

- verbose:

  logical: print progress while running? Default TRUE.

- gof:

  logical: calculate goodness of fit statistics? Setting to TRUE adds
  approximately 2-5% to runtime. For faster sampling without GOF
  overhead, set gof=FALSE and use gof() after model fitting.

- custom_gof:

  optional function or list of named functions for computing custom
  goodness-of-fit statistics. Each function must accept a single matrix
  Y as input and return a numeric vector. If a single function is
  provided, it should return a named vector. If a list of functions is
  provided, each function should return a single value and will be named
  according to the list names. Custom statistics will be computed in
  addition to default statistics. Example: custom_gof = function(Y)
  c(density = mean(Y \> 0, na.rm = TRUE))

- plot:

  accepted for signature parity with
  [`lame`](https://netify-dev.github.io/lame/reference/lame.md);
  ignored. `ame()` runs a single-period MCMC and does not draw live
  trace panels. Default `FALSE`.

- start_vals:

  List from previous model run containing parameter starting values for
  new MCMC

- periodic_save:

  logical: indicating whether to periodically save MCMC results

- out_file:

  character vector indicating name and path in which file should be
  stored if periodic_save is selected. For example, on an Apple OS
  out_file="~/Desktop/ameFit.rda".

- save_interval:

  quantile interval indicating when to save during post burn-in phase.

- model.name:

  optional string for model selection output

- posterior_opts:

  optional list of posterior sampling options

- n_chains:

  integer: number of MCMC chains to run (default: 1)

- cores:

  integer: number of cores for parallel chains (default: 1)

- use_sparse_matrices:

  logical: use sparse matrix storage for large networks? (default:
  FALSE). Recommended only for truly sparse networks (\< 10% non-zero
  entries).

- method:

  character: `"mcmc"` (default, the Bayesian MCMC fit) or `"als"` (the
  fast, MCMC-free iterative block coordinate descent point estimator).
  When `method = "als"`, MCMC-specific arguments (`nscan`, `burn`,
  `odens`, `prior`, ...) are silently ignored and the call forwards to
  [`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md).

- bootstrap:

  integer (only used when `method = "als"`): number of bootstrap
  replicates. `0` (default) skips the bootstrap; `N > 0` runs `N`
  replicates and attaches the result so that
  [`confint`](https://rdrr.io/r/stats/confint.html) returns bootstrap
  intervals.

- bootstrap_type:

  character (only used when `method = "als"`): `"parametric"` (default)
  or `"block"`.

- bootstrap_block_length:

  integer: block length for the block bootstrap.

- bootstrap_seed:

  optional integer seed for the bootstrap.

- save_log_lik:

  logical: when `TRUE`, attach a per-iteration pointwise log-likelihood
  matrix `fit$log_lik` (an `n_iter x n_obs` array) using the augmented-Z
  normal approximation. Required for `loo::loo(fit)` and
  `loo::waic(fit)`. The augmented-Z approximation is the same semantics
  used by
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md); for
  non-normal families it is the predictive criterion on the latent-Z
  scale, not the family-specific Y density. Default `FALSE` (no log_lik
  storage; default fit is byte-identical to previous releases).

- ordinal_cutpoints:

  character: cutpoint convention for `family = "ordinal"`.
  `"data_induced"` (default) uses the data-induced cutpoints;
  `"explicit"` samples explicit cutpoints via a Cowles (1996)
  Metropolis-Hastings update. Ignored for other families.

- print:

  Deprecated. Use `verbose` instead.

- ...:

  reserved for future use. Passing
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md)-only
  arguments (e.g. `dynamic_beta`, `period_exposure`) here triggers a
  clean abort directing you to
  [`lame`](https://netify-dev.github.io/lame/reference/lame.md); passing
  any other unrecognised name warns so typos are visible rather than
  silently dropped.

## Value

**Posterior Samples (full MCMC chains):**

- BETA:

  Regression coefficients (nscan xp matrix)

- VC:

  Variance components (nscan xk matrix)

- GOF:

  Goodness-of-fit statistics (nscan x4 matrix). First row contains
  observed values, remaining rows contain posterior predictive samples.
  See [`gof`](https://netify-dev.github.io/lame/reference/gof.md) for
  post-hoc computation and
  [`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md)
  for visualization.

**Posterior Means (averaged over chain):**

- APM:

  Additive row/sender effects (n-vector)

- BPM:

  Additive column/receiver effects (m-vector); NULL for symmetric
  networks

- U:

  Multiplicative row/sender factors (n xR matrix)

- V:

  Multiplicative column/receiver factors (m xR matrix); NULL for
  symmetric networks

- L:

  Eigenvalue matrix (R xR diagonal); symmetric networks only

- YPM:

  Posterior mean of Y on response scale (for predictions and imputing
  missing values)

**Metadata:**

- family:

  Model family (normal, binary, etc.)

- mode:

  Network mode (unipartite or bipartite)

- symmetric:

  Logical indicating if network is symmetric

- R:

  Dimension of multiplicative effects

**Optional Posterior Samples (if requested via posterior_options):**

- U_samples:

  Samples of U (n xR xiterations array)

- V_samples:

  Samples of V (m xR xiterations array)

- a_samples:

  Samples of row effects (n xiterations matrix)

- b_samples:

  Samples of column effects (m xiterations matrix)

**Note on reconstructing removed matrices:** To save memory, EZ
(expected latent network) and UVPM/ULUPM (multiplicative products) are
not stored but can be reconstructed using:

- `reconstruct_EZ(fit)` - Returns linear predictor (link scale, not
  response scale)

- `reconstruct_UVPM(fit)` - Returns U\\

**Generating posterior distributions:** Use
`simulate_posterior(fit, component="UV")` to generate posterior samples
for components where only means are stored, or use
[`posterior_options()`](https://netify-dev.github.io/lame/reference/posterior_options.md)
during model fitting to save full posterior samples.

- model.name:

  Name of the model (if provided)

## Details

This command provides posterior inference for parameters in AME models
of cross-sectional relational data, assuming one of eight possible data
types/models. The function supports both unipartite networks (square
adjacency matrices) and bipartite networks (rectangular adjacency
matrices with distinct row and column node sets) for single time point
analysis.

**Model.** For a dyad \\(i, j)\\ the AME linear predictor is

\$\$\eta\_{ij} = \beta_0 + x\_{ij}'\beta + a_i + b_j + u_i' D v_j,\$\$

where \\\beta\\ are regression coefficients on dyadic / nodal
covariates, \\a_i\\ is a row (sender) random effect, \\b_j\\ is a column
(receiver) random effect, and \\u_i' D v_j\\ is the multiplicative
latent-factor term (rank `R`). The observation model is \\Y\_{ij} \sim
F(\eta\_{ij}, \theta)\\ with \\F\\ specified by `family` (Gaussian for
`"normal"`, probit for `"binary"` / `"cbin"`, etc.). For unipartite
(`symmetric = FALSE`) the residual error has dyad-level correlation
\\\rho\\ between \\(i, j)\\ and \\(j, i)\\; for bipartite, dyad
correlation is fixed at 0.

**Priors (in brief).** \\\beta\\ has a Zellner-style g-prior (`g`);
\\(a_i, b_i)\\ are jointly Normal with covariance \\\Sigma\_{ab}\\
(Inverse-Wishart prior `Sab0 / eta0`); \\u_i, v_j\\ are independent
Normal with covariance \\\Sigma\_{uv}\\ (Inverse-Wishart prior
`Suv0 / kappa0`); the dyad-correlation \\\rho\\ has a Normal prior on
\\\mathrm{atanh}(\rho)\\. See `prior_summary(fit)` for the priors
actually used.

**Choosing R.** The multiplicative rank `R` controls the dimensionality
of latent homophily / heterogeneity not explained by covariates and
additive effects. `R = 0` fits an additive-only social-relations model;
`R = 1` or `2` is typical for small / medium networks; `R > floor(n/3)`
is rarely identifiable and will issue a warning. Latent factors capture
unobserved structure (clusters, hub patterns, transitive triangles the
covariates miss) and are accessed at `fit$U`, `fit$V`.

**Identifiability.** The latent factor term \\u_i' D v_j\\ is invariant
to rotation and reflection of \\U, V\\; the package canonicalises with
an SVD so successive draws are interpretable. For visual stability
across posterior summaries see
[`procrustes_align`](https://netify-dev.github.io/lame/reference/procrustes_align.md)
and
[`latent_positions`](https://netify-dev.github.io/lame/reference/latent_positions.md).

**Theoretical Foundation:**

The AME model decomposes network structure into several components:
\$\$y\_{ij} = \beta'x\_{ij} + a_i + b_j + u_i'v_j + \epsilon\_{ij}\$\$
where:

- \\\beta'x\_{ij}\\: Fixed effects of dyadic/nodal covariates

- \\a_i\\: Additive sender (row) effect for node i

- \\b_j\\: Additive receiver (column) effect for node j

- \\u_i'v_j\\: Multiplicative interaction between latent factors

- \\\epsilon\_{ij}\\: Dyadic error term (may be correlated)

This specification generalizes the social relations model (Warner et al.
1979) and latent space models (Hoff et al. 2002) within a unified
framework.

**Prior Distributions:**

The model uses conjugate and semi-conjugate priors where possible:

- Regression coefficients: \\\beta \sim N(0, g\sigma^2(X'X)^{-1})\\
  (g-prior)

- Additive effects: \\(a_i, b_i)' \sim N(0, \Sigma\_{ab})\\ jointly

- Covariance: \\\Sigma\_{ab} \sim IW(\eta_0, \eta_0 S\_{ab0})\\
  (inverse-Wishart)

- Multiplicative effects: Hierarchical shrinkage via \\\eta_0\\

- Dyadic correlation: \\\rho \sim Uniform(-1, 1)\\ with Metropolis
  updates

The inverse-Wishart prior on \\\Sigma\_{ab}\\ allows learning
correlation between sender and receiver effects, capturing reciprocity
patterns.

**Multiplicative Effects (Latent Factors):**

When R \> 0, the model includes R-dimensional latent factors:

- Asymmetric case: \\u_i, v_j \in \mathbb{R}^R\\ with \\u_i'v_j\\
  interaction

- Symmetric case: \\u_i = v_i\\ with eigendecomposition \\ULU'\\

- Captures homophily, transitivity, and community structure

- R chosen via model selection or set to 2-3 for visualization

**Estimation Algorithm:**

The model uses a Gibbs sampler with the following updates:

1.  Sample latent Z given parameters (data augmentation for non-normal
    families)

2.  Update regression coefficients \\\beta\\ via g-prior conjugate
    update

3.  Update additive effects (a,b) jointly with \\\beta\\

4.  Update covariance \\\Sigma\_{ab}\\ from inverse-Wishart

5.  Update multiplicative effects U,V via Gibbs or Metropolis-Hastings

6.  Update dyadic correlation \\\rho\\ via Metropolis-Hastings

7.  Update variance \\\sigma^2\\ (for continuous families)

**Standard Model Types:**

The following data types/models are available:

"normal": A normal AME model (identity link: \\E\[Y\] = \eta\\).

"tobit": A tobit AME model for censored continuous data. Values are
censored at zero, appropriate for non-negative continuous relational
data (identity link with censoring).

"binary": A binary probit AME model (probit link: \\P(Y=1) =
\Phi(\eta)\\).

"ordinal": An ordinal probit AME model (cumulative probit link). An
intercept is not identifiable in this model.

"cbin": An AME model for censored binary data (probit link with
censoring). The value of 'odmax' specifies the maximum number of links
each row may have.

"frn": An AME model for fixed rank nomination networks. A higher value
of the rank indicates a stronger relationship. The value of 'odmax'
specifies the maximum number of links each row may have.

"rrl": An AME model based on the row ranks. This is appropriate if the
relationships across rows are not directly comparable in terms of scale.
An intercept, row random effects and row regression effects are not
estimable for this model.

"poisson": An overdispersed Poisson AME model for count data (log link:
\\E\[Y\] = \exp(\eta)\\). The linear predictor \\\eta\\ represents
\\\log(\lambda)\\ where \\\lambda\\ is the expected count.

## Supported family x mode combinations

Every `family` is supported under both `mode`s. The bipartite Z-samplers
live in `R/rZ_bipartite.R` and dispatch per family; see also the inline
comment in `R/lame.R` (“first-class samplers as of the
bipartite-families round”).

|            |                |               |
|------------|----------------|---------------|
| **family** | **unipartite** | **bipartite** |
| normal     | yes            | yes           |
| binary     | yes            | yes           |
| tobit      | yes            | yes           |
| ordinal    | yes            | yes           |
| cbin       | yes            | yes           |
| frn        | yes            | yes           |
| poisson    | yes            | yes           |
| rrl        | yes            | yes           |

Symmetric (`symmetric = TRUE`) fits require a symmetric `Y`.
`family = "ordinal"` with `symmetric = TRUE` is supported via the
dedicated sampler in `R/rZ_ord_sym_fc.R`, which uses the
symmetric-doubled precision and mirrors upper-triangle draws to the
lower triangle so \\Z = t(Z)\\ holds at every sweep. The previous
“ordinal

- symmetric is unsupported

restriction has been lifted. `rrl` + bipartite is first-class: the
rectangular row-rank Z-sampler handles the full per-row permutation
likelihood, and `ame()` recovers regression coefficients comparably to
the unipartite rrl path on similarly-sized data. ”

**Symmetric input with one triangle missing.** When `symmetric = TRUE`
and one triangle of `Y` is fully `NA` (the user stored only the lower or
upper triangle), the symmetry validator `any(is.finite(Y - t(Y)))`
evaluates `FALSE` and the call proceeds: the sampler then treats the
populated triangle as the symmetric data and mirrors it. This is usually
intended, but if the upper / lower triangles were meant to differ, the
model is silently fitting half the data. Audit
`anyNA(Y[upper.tri(Y)]) != anyNA(Y[lower.tri(Y)])` before calling if you
are unsure.

## Converting a network from another R package

`ame()` expects `Y` as a (potentially named) matrix, not a graph object.
The most common conversions:

- igraph:
  `Y <- igraph::as_adjacency_matrix(g, sparse = FALSE); diag(Y) <- NA`

- network / statnet: `Y <- as.matrix(net); diag(Y) <- NA`

- tibble / long-format tidyverse:
  `Y <- tidyr::pivot_wider(df, names_from = receiver, values_from = tie) |> tibble::column_to_rownames("sender") |> as.matrix()`;
  or use the convenience helper
  [`as_lame_y`](https://netify-dev.github.io/lame/reference/as_lame_y.md).

For an undirected/symmetric network, pass `symmetric = TRUE`; for a
rectangular two-mode network (students x courses, donors x candidates),
pass `mode = "bipartite"`. `Xrow` and `Xcol` accept either a numeric
matrix or a data.frame (coerced internally); `Xdyad` must be a 3-D array
`n x n x p` of numeric covariates with no `Inf`/`NaN`.

## When to use AME vs ERGM

AME and ERGM are complementary tools for binary network analysis, not
direct substitutes. **ERGM** is a class of exponential-family models
built around explicit network statistics (counts of edges, mutual ties,
triangles, geometrically-weighted shared partners, ...). You write the
statistics you think matter, ERGM gives you their coefficients. ERGM
excels when you have a substantive theory about which configurations
drive tie formation.

**AME** models latent homophily / heterogeneity directly via sender,
receiver, and multiplicative latent-factor effects. You don't enumerate
triadic terms; the multiplicative-effects rank `R` captures higher-order
structure (clustering, transitivity, hub patterns) implicitly. AME
excels when (a) you have dyadic / nodal covariates whose effects you
want to interpret cleanly without ERGM degeneracy, (b) higher-order
structure is "nuisance" that you want to absorb but not parameterise, or
(c) you need a posterior distribution over predictions for forecasting
or imputation.

Practical guidance: if your research question is "do nodes that share
attribute X tend to form triangles together?", reach for ERGM's `gwesp`.
If your research question is "controlling for unobserved sender /
receiver heterogeneity and latent clustering, what is the effect of
dyadic covariate X?", reach for AME. The `R = 0` additive-only case is
the social relations model (Warner, Kenny, Stoto 1979); `R >= 1` adds
latent space.

## ERGM to AME translation

For users coming from `statnet::ergm`, the rough analogues are:

|  |  |
|----|----|
| **ERGM term** | **AME analogue** |
| `edges` | `intercept` (probit link, not logit) |
| `nodecov("x")` | `Xrow = x` or `Xcol = x` |
| `nodematch("g")` | dyadic covariate via [`nodematch`](https://netify-dev.github.io/lame/reference/nodematch.md)`(g)` into `Xdyad` |
| `nodefactor("g")` | dyadic covariate via [`nodefactor`](https://netify-dev.github.io/lame/reference/nodematch.md)`(g)` into `Xdyad` (drop one level) |
| `absdiff("z")` | dyadic covariate via [`absdiff`](https://netify-dev.github.io/lame/reference/nodematch.md)`(z)` into `Xdyad` |
| `mutual` | `dcor = TRUE` -\> the `rho` parameter (probit-scale, NOT log-odds; not numerically comparable to ERGM's `mutual`) |
| `gwesp` / transitivity | `R >= 1` multiplicative latent factors (not the same statistic) |
| sender activity heterogeneity | `rvar = TRUE`, gives `a_i`, `va` |
| receiver popularity heterog. | `cvar = TRUE`, gives `b_j`, `vb` |

For `family = "binary"` the link is probit, so the intercept is on the
probit scale; do not compare it to an ERGM `edges` estimate by simple
arithmetic.

## Migration from amen

Both `amen` and `lame` export `ame()`; loading both packages fires a
startup warning telling you to call `lame::ame(...)` or `amen::ame(...)`
explicitly.

For default cross-sectional calls, `lame::ame()` is a drop-in
replacement for `amen::ame()`: the `fit$BETA` slot is a 2-D
`[n_stored, p]` matrix in both packages, so scripts that call
`colMeans(fit$BETA)` or `apply(fit$BETA, 2, mean)` continue to work
unchanged. `lame` additionally accepts `family = "binary"` (which `amen`
1.4.5 no longer accepts; `amen` requires `"bin"`). The `print` argument
is deprecated in favour of `verbose`; calls that pass `print = ...`
still work but warn.

The cross-sectional path has no `dynamic_beta` option (an AR(1) prior on
a single-period coefficient is unidentified), so the `BETA` 3-D shape
that [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) can
produce never arises from `ame()`. See
[`lame`](https://netify-dev.github.io/lame/reference/lame.md) for the
longitudinal path and the silent-aggregation hazard with old
`apply(fit$BETA, 2, mean)` scripts under `dynamic_beta = TRUE`.

## Notes on priors

- The regression-coefficient prior is a Zellner g-prior
  (`beta ~ N(0, g * sigma^2 * (XtX)^-1)` where `XtX` is the design
  cross-product). `g` is a *top-level* argument of `ame()`, not an entry
  in `prior = list(...)` (a common slip).

- There is no per-coefficient prior knob. If you need student_t,
  horseshoe, or to centre a slope away from zero, AME does not currently
  expose it; the g-prior structure is the only knob.

- Unknown names in `prior = list(...)` are warned about (a typo like
  `Sab = ...` instead of `Sab0` would otherwise be silently dropped).

## See also

[`lame`](https://netify-dev.github.io/lame/reference/lame.md) for
longitudinal models,
[`gof`](https://netify-dev.github.io/lame/reference/gof.md) for post-hoc
goodness-of-fit computation,
[`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md)
for visualizing GOF results,
[`latent_positions`](https://netify-dev.github.io/lame/reference/latent_positions.md)
for extracting latent positions as a tidy data frame,
[`procrustes_align`](https://netify-dev.github.io/lame/reference/procrustes_align.md)
for Procrustes alignment of latent positions,
[`summary.ame`](https://netify-dev.github.io/lame/reference/summary.ame.md)
for model summaries,
[`coef.ame`](https://netify-dev.github.io/lame/reference/coef.ame.md)
for coefficient extraction

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
data(YX_bin)
fit <- ame(YX_bin$Y, Xdyad = YX_bin$X, burn = 10, nscan = 100, odens = 1,
           family = "binary", verbose = FALSE)
summary(fit)
#> 
#> === AME Model Summary ===
#> 
#> Call:
#> [1] "Y ~ dyad(intercept, rgpa, rsmoke, cgpa, csmoke, igrade, igpa, ismoke) + a[i] + b[j] + rho*e[ji], family = 'binary'"
#> 
#> Regression coefficients:
#> ------------------------
#>                Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept_dyad   -2.543    0.114 -22.296       0   -2.772   -2.319 ***
#> rgpa_dyad         0.178    0.081   2.193   0.028    0.047     0.32   *
#> rsmoke_dyad       0.254    0.112   2.279   0.023    0.054    0.474   *
#> cgpa_dyad         0.188    0.042   4.472       0    0.107    0.273 ***
#> csmoke_dyad       0.167    0.052   3.248   0.001    0.067     0.27  **
#> igrade_dyad        1.13    0.052  21.782       0    1.032    1.229 ***
#> igpa_dyad         0.049    0.038   1.284   0.199   -0.045    0.123    
#> ismoke_dyad       0.038    0.063   0.601   0.548   -0.084     0.15    
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: stars are a visual hint from posterior mean / SD only; for inference use the credible intervals.
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.377    0.088
#> cab    0.045    0.019
#> vb     0.073    0.011
#> rho    0.723    0.085
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
# Note: you should run the Markov chain much longer in practice
# }
 
```
