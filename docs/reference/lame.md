# AME model fitting routine for longitudinal relational data

An MCMC routine providing a fit to an additive and multiplicative
effects (AME) regression model to longitudinal (time-series) relational
data of various types. Supports both unipartite (square) and bipartite
(rectangular) network structures. For cross-sectional (single time
point) networks, use the `ame` function.

## Usage

``` r
lame(
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
  dynamic_uv = FALSE,
  dynamic_ab = FALSE,
  dynamic_G = FALSE,
  dynamic_beta = FALSE,
  dynamic_beta_kind = c("ar1", "rw1", "rw2", "matern32"),
  family = "normal",
  intercept = !is.element(family, c("rrl", "ordinal")),
  symmetric = FALSE,
  odmax = NULL,
  prior = list(),
  g = NA,
  seed = 6886,
  nscan = 10000,
  burn = 500,
  odens = 25,
  plot = FALSE,
  verbose = FALSE,
  gof = TRUE,
  start_vals = NULL,
  periodic_save = FALSE,
  out_file = NULL,
  save_interval = 0.25,
  model.name = NULL,
  save_log_lik = FALSE,
  log_lik_path = NULL,
  log_lik_chunk_size = 10000L,
  freeze_call = FALSE,
  dynamic_beta_pool = c("none", "rho", "sigma", "both"),
  dynamic_beta_per_actor = NULL,
  per_actor_covariate_idx = 1L,
  per_actor_identifiability = c("center", "exact_center", "drop_population"),
  keep_per_actor = c("auto", "draws", "summary", "none"),
  time_index = NULL,
  period_exposure = NULL,
  max_seconds = Inf,
  checkpoint_path = NULL,
  checkpoint_every = 100L,
  log_lik_method = c("observed_exact", "observed_ghk", "augmented"),
  ordinal_cutpoints = c("data_induced", "explicit"),
  method = c("mcmc", "als"),
  bootstrap = 0L,
  bootstrap_type = c("parametric", "block"),
  bootstrap_block_length = 1L,
  bootstrap_seed = NULL,
  resume_from = NULL,
  print
)
```

## Arguments

- Y:

  a T length list of n x n relational matrices, or a 3D array of
  dimensions `[n, n, T]`, where T corresponds to the number of
  replicates (over time, for example). If a 3D array is provided, it is
  automatically converted to list format. See family below for various
  data types.

- Xdyad:

  a T length list of n x n x pd arrays of covariates

- Xrow:

  a T length list of n x pr matrices of nodal row covariates

- Xcol:

  a T length list of n x pc matrices of nodal column covariates

- rvar:

  logical: fit row random effects (asymmetric case)?

- cvar:

  logical: fit column random effects (asymmetric case)?

- dcor:

  logical: fit a dyadic correlation (asymmetric case)?

- nvar:

  logical: fit nodal random effects (symmetric case)?

- R:

  integer: dimension of the multiplicative effects (can be zero)

- R_row:

  integer: for bipartite networks, dimension of row node multiplicative
  effects (defaults to R)

- R_col:

  integer: for bipartite networks, dimension of column node
  multiplicative effects (defaults to R)

- mode:

  character: either "unipartite" (default) for square networks or
  "bipartite" for rectangular networks

- dynamic_uv:

  logical: fit dynamic multiplicative effects (latent factors) that
  evolve over time using AR(1) processes. When TRUE, the latent
  positions U and V become time-varying, following \\U\_{i,t} =
  \rho\_{uv} U\_{i,t-1} + \epsilon\_{i,t}\\, where epsilon follows N(0,
  sigma_uv^2). This allows actors' positions in latent social space to
  drift smoothly over time, capturing evolving network structure and
  community dynamics. Inspired by dynamic latent space models (Sewell
  and Chen 2015, "Latent Space Models for Dynamic Networks", JASA;
  Durante and Dunson 2014, "Nonparametric Bayes Dynamic Modeling of
  Relational Data", Biometrika). The implementation uses efficient
  blocked Gibbs sampling with C++ acceleration for scalability. Default
  FALSE.

- dynamic_ab:

  logical: fit dynamic additive effects (sender/receiver effects) that
  evolve over time using AR(1) processes. When TRUE, the row effects (a)
  and column effects (b) become time-varying, following \\a\_{i,t} =
  \rho\_{ab} a\_{i,t-1} + \epsilon\_{i,t}\\. This captures temporal
  heterogeneity in actors' baseline propensities to send and receive
  ties, allowing for smooth changes in activity levels and popularity
  over time. For example, an actor's tendency to form outgoing ties
  might gradually increase or decrease across observation periods. The
  AR(1) specification ensures temporal smoothness while allowing for
  actor-specific evolution patterns. Implementation uses conjugate
  updates where possible and C++ for computational efficiency. Default
  FALSE.

- dynamic_G:

  logical (bipartite only). When TRUE the bipartite interaction matrix
  \\G\\ is sampled per period and returned as `fit$G_cube`, an
  \\R\_\text{row} \times R\_\text{col} \times T\\ array. The static
  \\G\\ (single matrix) is the default and is the recommended choice;
  `dynamic_G = TRUE` is experimental and emits a one-line note at the
  start of the run. Unipartite fits ignore this argument with a warning.
  Default FALSE.

- dynamic_beta:

  logical, character, integer, or logical-vector flag selecting which
  regression coefficients evolve over time via independent AR(1)
  processes. Default FALSE keeps every coefficient static (the
  historical behaviour). Accepted forms:

  - `FALSE` / `NULL`: no coefficient is dynamic.

  - `TRUE`: every coefficient is dynamic.

  - character vector of block shortcuts (`"intercept"`, `"dyad"`,
    `"row"`, `"col"`) or specific coefficient names from
    `colnames(BETA)` — those become dynamic.

  - integer vector: 1-based column indices into `colnames(BETA)`.

  - logical vector of length `p`: per-coefficient mask.

  Each dynamic-coefficient block (one per distinct intercept / dyad /
  row / col label that has at least one dynamic coefficient) gets its
  own AR(1) parameters \\\rho\_\beta\\ (truncated-Normal prior, default
  mean 0.8, bounded between 0 and 0.999) and \\\sigma\_\beta^2\\
  (inverse-Gamma prior, defaults shape 2 / scale 1). The joint posterior
  over the time path \\\beta_t\\ is drawn by Forward-Filter /
  Backward-Sample (FFBS) inside the MCMC loop, conditional on the
  current static-block beta and on (a, b, U, V). When the intercept or a
  nodal coefficient is dynamic, a sum-to-zero contrast basis is applied
  to the additive-effect sampler to keep \\a_i + intercept_t\\
  identified. Requires at least 2 time periods, ignored for the ALS
  estimator (`method = "als"`), and unavailable for the intercept when
  `family = "ordinal"` or `"rrl"` (or for row-nodal coefficients with
  `family = "rrl"`). Default FALSE.

- dynamic_beta_kind:

  character: state-space prior on the dynamic coefficient block(s).
  `"ar1"` (default) gives mean-reverting AR(1) with a truncated-Normal
  prior on \\\rho\_\beta\\. `"rw1"` gives a random walk (\\\rho\_\beta\\
  pinned at 1, no truncation), appropriate when you expect permanent
  drift with no mean reversion – e.g. trade-gravity coefficients in a
  permanently changing world economy. The alias `"random_walk"` is also
  accepted for `"rw1"`. When the stationarity warning in `summary(fit)`
  fires under the default AR(1), refit with `dynamic_beta_kind = "rw1"`
  for cleaner inference. Decision tree: **mean-reverting?** AR(1).
  **Permanent / unit-root drift?** RW1. **Smooth with curvature?**
  `"rw2"` (second-order random walk). **Smooth with a known
  length-scale?** `"matern32"` (Matern 3/2). Note: `"rw2"` and
  `"matern32"` use an R-level joint Gaussian sampler and run
  substantially slower per iteration than the C++ FFBS used for `"ar1"`
  / `"rw1"` (roughly 2-4x on n = 100, T = 10).

- family:

  character: one of
  "normal","tobit","binary","ordinal","cbin","frn","rrl","poisson" - see
  the details below

- intercept:

  logical: fit model with an intercept?

- symmetric:

  logical: is the sociomatrix symmetric?

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
      values impose stronger shrinkage on the latent factors.

  etaab

  :   Prior degrees of freedom for covariance of additive effects
      (default: 4 + 3 \\ n/100). Controls shrinkage of row/column random
      effects.

  rho_uv_mean

  :   For dynamic_uv=TRUE: Prior mean for UV AR(1) parameter (default:
      0.9). Values close to 1 indicate high temporal persistence.

  rho_uv_sd

  :   For dynamic_uv=TRUE: Prior SD for UV AR(1) parameter (default:
      0.1). Controls uncertainty about temporal dependence.

  sigma_uv_shape

  :   For dynamic_uv=TRUE: Shape parameter for inverse-gamma prior on UV
      innovation variance (default: 2).

  sigma_uv_scale

  :   For dynamic_uv=TRUE: Scale parameter for inverse-gamma prior on UV
      innovation variance (default: 1).

  rho_ab_mean

  :   For dynamic_ab=TRUE: Prior mean for additive effects AR(1)
      parameter (default: 0.8). Controls temporal smoothness of
      sender/receiver effects.

  rho_ab_sd

  :   For dynamic_ab=TRUE: Prior SD for additive effects AR(1) parameter
      (default: 0.15).

  sigma_ab_shape

  :   For dynamic_ab=TRUE: Shape parameter for inverse-gamma prior on
      additive effects innovation variance (default: 2).

  sigma_ab_scale

  :   For dynamic_ab=TRUE: Scale parameter for inverse-gamma prior on
      additive effects innovation variance (default: 1).

  rho_beta_mean

  :   For dynamic_beta: Prior mean for the per-block AR(1) parameter on
      time-varying regression coefficients (default: 0.8). Closer to 1 =
      smoother evolution.

  rho_beta_sd

  :   For dynamic_beta: Prior SD for the per-block AR(1) parameter
      (default: 0.15).

  rho_beta_lower

  :   For dynamic_beta: Lower truncation bound on the AR(1) parameter
      (default: 0). Pass a negative value to allow negative
      autoregression.

  rho_beta_upper

  :   For dynamic_beta: Upper truncation bound (default: 0.999). Closer
      to 1 admits near-unit-root behaviour.

  sigma_beta_shape

  :   For dynamic_beta: Shape parameter for the inverse-Gamma prior on
      the per-block innovation variance (default: 2).

  sigma_beta_scale

  :   For dynamic_beta: Scale parameter for the inverse-Gamma prior
      (default: 1).

  sigma_beta_init

  :   For dynamic_beta: Initial value of the per-block innovation
      standard deviation (default: 0.25). Affects mixing, not the
      stationary distribution.

  beta0_mean

  :   For dynamic_beta: Mean of the Gaussian prior on the dynamic
      coefficients at t = 0 (default: 0).

  beta0_var

  :   For dynamic_beta: Variance of the Gaussian prior on the dynamic
      coefficients at t = 0 (default: 10). Weakly informative.

  Common usage: prior = list(Sab0 = diag(c(1, 1)), eta0 = 10) for
  stronger shrinkage, or prior = list(rho_uv_mean = 0.95) for higher
  temporal persistence, or prior = list(rho_beta_mean = 0.95,
  sigma_beta_scale = 0.1) for very smooth time-varying coefficients with
  tight innovations.

- g:

  optional scalar or vector for g-prior on regression coefficients.
  Default is p^2 where p is the number of regression parameters. The
  g-prior controls the variance of regression coefficients: larger
  values allow for larger coefficient values. Can be a vector of length
  p for parameter-specific control.

- seed:

  random seed

- nscan:

  number of iterations of the Markov chain (beyond burn-in)

- burn:

  burn in for the Markov chain

- odens:

  output density for the Markov chain

- plot:

  logical: plot results while running?

- verbose:

  logical: print progress while running? Default FALSE.

- gof:

  logical: calculate goodness of fit statistics?

- start_vals:

  list of parameter starting values for the MCMC chain

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

- save_log_lik:

  one of `FALSE` (default), `TRUE`, or `"chunked"`. When `TRUE`, stores
  the per-iteration pointwise log-likelihood matrix on `fit$log_lik` (an
  `[n_stored, n_obs]` double matrix). When `"chunked"`, streams the
  log-lik values to per-column-chunk binary files under `log_lik_path`
  so the in-memory cost during MCMC is just one chunk's width; the
  chunks are recovered later via `read_log_lik(fit)`. Required for
  `loo::loo(fit)` and `loo::waic(fit)`.

- log_lik_path:

  directory to write log-lik chunks to when `save_log_lik = "chunked"`.
  Default `NULL` creates a session-scoped subdirectory under
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- log_lik_chunk_size:

  column-width of each on-disk chunk when `save_log_lik = "chunked"`.
  Larger chunks mean fewer files (and one `readBin` per chunk on read)
  but more memory during MCMC. Default `10000L`.

- freeze_call:

  logical: if `TRUE`, store a snapshot of the evaluated `Y`, `Xdyad`,
  `Xrow`, `Xcol` on `fit$data_snapshot` so that a later
  `update(fit, ...)` refits against the same data even if the caller has
  mutated those objects in their workspace. Memory cost equals the size
  of the data; default `FALSE`.

- dynamic_beta_pool:

  one of `"none"` (default), `"rho"`, `"sigma"`, or `"both"`. Stable
  argument name for a hierarchical pooling prior on the per-block
  \\\rho_b\\ and / or \\\sigma_b\\ dynamic hyperparameters. The value is
  recorded on the fit so future-phase samplers and downstream tooling
  can read it uniformly; the hierarchical sampler ships in a later
  release, at which point non-`"none"` values will take effect.

- dynamic_beta_per_actor:

  optional, one of `NULL` (default), `"row"`, or `"col"`. **Currently a
  wired:** the argument is validated and stored on the fit, but the
  in-loop sampler is not yet wired. For per-actor time-varying slopes
  today, fit a standard `lame()` and then call
  [`per_actor_slopes`](https://netify-dev.github.io/lame/reference/per_actor_slopes.md)`(fit, kind, covariate_idx, lambda)`
  for the post-MCMC penalised-LS estimate.

- per_actor_covariate_idx:

  positive integer; index into the dyadic covariate cube to slope on for
  the per-actor extension. Default `1L`.

- per_actor_identifiability:

  one of `"center"` (default) or `"drop_population"`. `"center"`
  preserves the population coefficient and constrains per-actor
  deviations to sum to zero per period via pairwise contrast FFBS.
  Reserved API for the `"drop_population"` mode (not yet wired).

- keep_per_actor:

  one of `"auto"` (default; full draws when memory cost \< 250 MB else
  streaming summary), `"draws"` (full `[n_iter, n_actors, T]` cube),
  `"summary"` (streaming posterior mean + sd only), or `"none"` (only
  hyperparameter chains).

- time_index:

  optional numeric vector of length \\T\\ giving observation times.
  Default `NULL` treats periods as equally spaced. Strictly-increasing
  values are required when supplied. The gap-aware AR(1) update ships in
  a later release; the value is recorded on the fit.

- period_exposure:

  optional non-negative numeric vector of length \\T\\ giving
  period-level exposure offsets. **Wired for Poisson:** when supplied
  with any value != 1, the Poisson observation likelihood becomes
  \\Y\_{ij,t} \sim \text{Poisson}(e_t \cdot \exp(Z\_{ij,t}))\\, while
  the latent \\Z\\ retains its existing semantic of "unexposed
  log-rate". For other families, non-trivial `period_exposure` is an
  error (rescale `Y` or use a covariate offset). `NULL` (default) or
  all-ones uses the unscaled likelihood path.

- max_seconds:

  optional positive scalar; if the MCMC wall-clock time exceeds this
  many seconds, the chain terminates cleanly and `fit$terminated_early`
  is set to `TRUE`.

- checkpoint_path:

  optional file path. When set, the chain periodically writes a snapshot
  of `BETA`, `VC`, the RNG state, and the original call to this file.
  Use
  [`lame_resume`](https://netify-dev.github.io/lame/reference/lame_resume.md)
  or pass the same path back as `resume_from = path` to continue.

- checkpoint_every:

  positive integer; iterations between checkpoint writes. Default
  `100L`.

- log_lik_method:

  one of `"observed_exact"` (default), `"observed_ghk"`, or
  `"augmented"`. Selects the argument for selecting which pointwise
  log-lik is stored on `fit$log_lik`. `"observed_exact"` uses the
  closed-form marginal log-likelihood (currently wired for normal,
  binary, cbin, tobit, poisson, ordinal). `"observed_ghk"` is a stable
  api stub for the GHK Monte Carlo marginal on cbin/frn/rrl/ordinal that
  ships in a later release; it falls back to `"augmented"` now and emits
  a one-line note. `"augmented"` uses the augmented- data
  Gaussian-on-\\Z\\ contribution.

- ordinal_cutpoints:

  character: cutpoint convention for `family = "ordinal"`.
  `"data_induced"` (default) uses the data-induced cutpoints;
  `"explicit"` samples explicit cutpoints via a Cowles (1996)
  Metropolis-Hastings update. Ignored for other families.

- method:

  character: `"mcmc"` (default, the Bayesian MCMC fit) or `"als"` (the
  fast, MCMC-free iterative block coordinate descent point estimator;
  pooled-static over time). When `method = "als"`, MCMC- and
  dynamics-specific arguments (`nscan`, `burn`, `odens`, `prior`,
  `dynamic_uv`, `dynamic_ab`, ...) are silently ignored and the call
  forwards to
  [`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md).

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

- resume_from:

  optional path to a checkpoint file produced by a previous
  `lame(..., checkpoint_path = path)` call. When non-`NULL`, `lame()`
  short-circuits all other input parsing and delegates to
  [`lame_resume`](https://netify-dev.github.io/lame/reference/lame_resume.md)
  with the user-supplied overrides forwarded. Pass `nscan = K` to
  request `K` additional stored draws on the continuation. Note:
  `checkpoint_path` and `max_seconds` cannot currently be overridden on
  a resume (they are stripped before re-evaluating the saved call); open
  both at the original `lame()` call if you need iterative checkpointing
  or a time budget on a chain of resumes. Pass `verbose = FALSE` on
  resume (the sampling progress bar is gated on `burn != 0` and the
  continuation forces `burn = 0`). Default `NULL`.

- print:

  Deprecated. Use `verbose` instead.

## Value

- BETA:

  posterior samples of regression coefficients. A 2-dimensional matrix
  `[n_stored, p]` when all coefficients are static. When `dynamic_beta`
  flags any coefficient as time-varying, `BETA` is a 3-dimensional array
  `[n_stored, p, T]` whose third dimension is the time period. Static
  coefficients are still present with their values replicated across the
  third dimension. `coef(fit)` collapses this to a `[p, T]`
  posterior-mean matrix.

  **Migration from amen.** Under
  [`amen::ame()`](https://rdrr.io/pkg/amen/man/ame.html) the `BETA` slot
  was always 2-D. Scripts that compute `apply(fit$BETA, 2, mean)` on an
  amen fit will silently aggregate across periods when run against a
  `lame()` fit with `dynamic_beta` active (the second margin is the
  coefficient index in both shapes; the new third margin is time). Use
  `length(dim(fit$BETA))` to detect the shape, or call `coef(fit)` which
  returns a `[p, T]` matrix in either case.

- VC:

  posterior samples of the variance parameters

- APM:

  posterior mean of additive row effects a

- BPM:

  posterior mean of additive column effects b

- U:

  posterior mean of multiplicative row effects u. For dynamic_uv=TRUE,
  this is a 3D array (n x R x T)

- V:

  posterior mean of multiplicative column effects v (asymmetric case).
  For dynamic_uv=TRUE, this is a 3D array (n x R x T)

- UVPM:

  posterior mean of UV

- ULUPM:

  posterior mean of ULU (symmetric case)

- L:

  posterior mean of L (symmetric case)

- EZ:

  estimate of expectation of Z matrix. For `mode = "bipartite"`, `EZ` is
  a list of per-period matrices whose row/column dimnames are currently
  `NULL`. Use `names(fit$APM)` / `names(fit$BPM)` (or
  `dimnames(fit$YPM[[t]])`) to recover the row/column actor ordering,
  which is the lexicographic sort of the input actor names (so `"r10"`
  sorts before `"r2"` unless you zero-pad).

- YPM:

  posterior mean of Y (for imputing missing values)

- GOF:

  observed (first row) and posterior predictive (remaining rows) values
  of four goodness-of-fit statistics. See
  [`gof`](https://netify-dev.github.io/lame/reference/gof.md) for
  post-hoc computation and
  [`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md)
  for visualization.

- start_vals:

  Final parameter values from MCMC, can be used as the input for a
  future model run.

- model.name:

  Name of the model (if provided)

## Details

This command provides posterior inference for parameters in AME models
of longitudinal relational data, assuming one of eight possible data
types/models. The model supports both unipartite networks (square
adjacency matrices) and bipartite networks (rectangular adjacency
matrices with distinct row and column node sets) across multiple time
points.

**Dynamic Effects Implementation:**

The dynamic_uv and dynamic_ab parameters enable time-varying latent
representations through autoregressive processes. These extensions are
particularly useful for understanding how network structure evolves over
time.

*Dynamic Multiplicative Effects (dynamic_uv=TRUE):* The latent factors U
and V evolve according to AR(1) processes: \$\$U\_{i,k,t} = \rho\_{uv}
U\_{i,k,t-1} + \epsilon\_{i,k,t}\$\$ where \\\epsilon\_{i,k,t} \sim N(0,
\sigma\_{uv}^2)\\, i indexes actors, k indexes latent dimensions, and t
indexes time. The parameter \\\rho\_{uv}\\ controls temporal persistence
(values near 1 indicate slow evolution). This captures time-varying
homophily, latent community structure, and transitivity dynamics.

Key references:

- Sewell & Chen (2015): Introduced dynamic latent space models with
  actor-specific evolution rates

- Durante & Dunson (2014): Nonparametric Bayesian approach allowing
  flexible evolution of network structure

- Hoff (2011): Hierarchical multilinear models providing theoretical
  foundation for temporal dependencies

*Dynamic Additive Effects (dynamic_ab=TRUE):* The sender (a) and
receiver (b) effects evolve as: \$\$a\_{i,t} = \rho\_{ab} a\_{i,t-1} +
\epsilon\_{i,t}\$\$ \$\$b\_{i,t} = \rho\_{ab} b\_{i,t-1} +
\eta\_{i,t}\$\$ where \\\epsilon\_{i,t}, \eta\_{i,t} \sim N(0,
\sigma\_{ab}^2)\\. This models time-varying individual activity levels
(outdegree) and popularity (indegree).

Applications include:

- Tracking changes in node centrality over time

- Identifying emerging influential actors

- Detecting declining activity patterns

- Modeling life-cycle effects in social networks

*Prior Specification for Dynamic Parameters:*

- \\\rho\_{uv}, \rho\_{ab} \sim TruncNormal(mean, sd, 0, 1)\\: Ensures
  stationarity of AR(1) process

- \\\sigma\_{uv}^2, \sigma\_{ab}^2 \sim InverseGamma(shape, scale)\\:
  Controls innovation variance

- Default priors (\\\rho\_{uv}\\ mean=0.9, \\\rho\_{ab}\\ mean=0.8)
  favor smooth evolution

- Adjust rho\_\*\_mean closer to 1 for slower evolution, closer to 0 for
  more rapid changes

*Computational Considerations:*

- Dynamic effects increase computation by ~30-50\\

- Memory usage scales as O(n*R*T) for dynamic_uv, O(n\*T) for dynamic_ab

- C++ implementation provides ~70\\

- Convergence diagnostics: Monitor rho and sigma parameters carefully

- Effective sample sizes typically lower due to temporal correlation

- Recommend burn \>= 1000 and nscan \>= 20000 for dynamic models

*Model Selection Guidelines:* Use both dynamic_uv and dynamic_ab when:

- Networks show clear temporal trends in density or clustering

- Individual node behavior changes systematically over time

- Community structure evolves (merging, splitting, drift)

Use only dynamic_uv when:

- Latent structure/communities change but individual effects are stable

- Focus is on evolving homophily or clustering patterns

- Network shows structural reconfiguration over time

Use only dynamic_ab when:

- Individual heterogeneity varies but overall structure is stable

- Actors' activity/popularity changes over observation period

- Focus is on individual-level temporal dynamics

**Bipartite Network Models:**

When mode="bipartite", the model handles rectangular adjacency matrices
Y with dimensions n_A x n_B, where n_A and n_B represent the number of
row and column nodes respectively.

*Static Bipartite Case:* The model uses separate latent factor matrices:

- U: n_A x R_row matrix of row node latent positions

- V: n_B x R_col matrix of column node latent positions

- G: R_row x R_col interaction matrix mapping between latent spaces

- Multiplicative term: U G V' captures bipartite community structure

*Dynamic Bipartite Case:* When dynamic_uv=TRUE for bipartite networks:
\$\$U\_{i,k,t} = \rho\_{uv} U\_{i,k,t-1} + \epsilon\_{i,k,t}\$\$
\$\$V\_{j,k,t} = \rho\_{uv} V\_{j,k,t-1} + \eta\_{j,k,t}\$\$ where i
indexes row nodes, j indexes column nodes, k indexes latent dimensions.

When `dynamic_G = TRUE` the bipartite interaction matrix \\G_t\\ is
sampled per period via the bipartite \\G\\-sampler (a separate
\\R\_\text{row} \times R\_\text{col}\\ matrix at every time slice) and
returned as `fit$G_cube`. This is experimental: the canonicalisation
that absorbs rotation and scaling into \\(U, V)\\ is enforced per
period, so the marginal \\U_t G_t V_t'\\ linear predictor is identified
even though the individual \\G_t\\ entries are not. Always sanity-check
against a `dynamic_G = FALSE` fit before relying on per-period \\G_t\\
estimates. `dynamic_G = TRUE` is bipartite only and is ignored (with a
warning) for unipartite networks.

*Key Differences from Unipartite Models:*

- No dyadic correlation (rho): Bipartite edges are inherently directed

- Separate dimensions: R_row and R_col can differ for row/column spaces

- Rectangular structure: Network density patterns differ from square
  matrices

- Community interpretation: Captures affiliation patterns between node
  types

**Standard AME Model Types:**

The following describes the eight standard data types/models available:

"normal": A normal AME model.

"tobit": A tobit AME model for censored continuous data. Values are
censored at zero, appropriate for non-negative continuous relational
data.

"binary": A binary probit AME model.

"ordinal": An ordinal probit AME model. An intercept is not identifiable
in this model.

"cbin": An AME model for censored binary data. The value of 'odmax'
specifies the maximum number of links each row may have.

"frn": An AME model for fixed rank nomination networks. A higher value
of the rank indicates a stronger relationship. The value of 'odmax'
specifies the maximum number of links each row may have.

"rrl": An AME model based on the row ranks. This is appropriate if the
relationships across rows are not directly comparable in terms of scale.
An intercept, row random effects and row regression effects are not
estimable for this model.

"poisson": An overdispersed Poisson AME model for count data. The latent
variable represents the log mean of the Poisson distribution.

## See also

[`ame`](https://netify-dev.github.io/lame/reference/ame.md) for
cross-sectional models,
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
for the fast MCMC-free point estimator,
[`als_dynamic_beta`](https://netify-dev.github.io/lame/reference/als_dynamic_beta.md)
for a regression-only penalised smoother on the time-varying coefficient
path,
[`lame_resume`](https://netify-dev.github.io/lame/reference/lame_resume.md)
for the legacy resume entry point (equivalent to
`lame(resume_from = path)`),
[`gof`](https://netify-dev.github.io/lame/reference/gof.md) for post-hoc
goodness-of-fit computation,
[`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md)
for visualizing GOF results,
[`latent_positions`](https://netify-dev.github.io/lame/reference/latent_positions.md)
for extracting latent positions as a tidy data frame,
[`procrustes_align`](https://netify-dev.github.io/lame/reference/procrustes_align.md)
for Procrustes alignment of latent positions,
[`summary.lame`](https://netify-dev.github.io/lame/reference/summary.lame.md)
for model summaries,
[`coef.lame`](https://netify-dev.github.io/lame/reference/coef.ame.md)
for coefficient extraction,
[`netify_to_lame`](https://netify-dev.github.io/lame/reference/netify_to_lame.md)
for the recommended netify -\> lame bridge (set `lame = TRUE` on
[`netify::to_amen()`](https://netify-dev.github.io/netify/reference/netify_to_amen.html)).

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
data(YX_bin_list)
fit<-lame(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,family="binary")
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
# you should run the Markov chain much longer than this

# \donttest{
## Time-varying regression coefficients (dynamic_beta).
## Make every dyadic coefficient evolve as an AR(1):
fit_dyn <- lame(YX_bin_list$Y, YX_bin_list$X,
                family = "binary", R = 0,
                nscan = 200, burn = 50, odens = 5,
                dynamic_beta = "dyad")
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
dim(fit_dyn$BETA)        # [n_stored, p, T] -- 3-D when dynamic
#> [1] 40  4  4
coef(fit_dyn)            # [p, T] posterior-mean coefficient paths
#>                  t1        t2        t3        t4
#> intercept 0.1199238 0.1199238 0.1199238 0.1199238
#> X1_dyad   0.7118277 0.7078791 0.6908902 0.6654124
#> X2_dyad   0.8650829 0.8748416 0.8417964 0.7950323
#> X3_dyad   1.0582930 1.0777578 1.0792037 0.9840359
confint(fit_dyn)         # per-period 95% credible intervals
#>                     2.5%     97.5%
#> intercept[t1] 0.07566298 0.1711694
#> X1_dyad[t1]   0.53743396 0.8624926
#> X2_dyad[t1]   0.62135516 1.0598551
#> X3_dyad[t1]   0.76479496 1.3208536
#> intercept[t2] 0.07566298 0.1711694
#> X1_dyad[t2]   0.50982957 0.8755414
#> X2_dyad[t2]   0.63062747 1.0549406
#> X3_dyad[t2]   0.77008591 1.2578084
#> intercept[t3] 0.07566298 0.1711694
#> X1_dyad[t3]   0.54407504 0.8382262
#> X2_dyad[t3]   0.65171342 0.9979983
#> X3_dyad[t3]   0.84688357 1.2833791
#> intercept[t4] 0.07566298 0.1711694
#> X1_dyad[t4]   0.51387021 0.7652479
#> X2_dyad[t4]   0.60842819 0.9169838
#> X3_dyad[t4]   0.75128416 1.1401657
summary(fit_dyn)         # prints a "Dynamic coefficients per period" block
#> 
#> === Longitudinal AME Model Summary ===
#> 
#> Call:
#> [1] "Y ~ dyad(X1, X2, X3) + a[i] + b[j] + rho*e[ji], family = 'binary'"
#> 
#> Time periods: 4 
#> Family: binary 
#> Mode: unipartite 
#> Dynamic regression coefficients: enabled, kind = 'ar1' (rho_beta = dyad : 0.93)
#> 
#> Regression coefficients:
#> ------------------------
#>           Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept     0.12    0.031   3.857       0    0.076    0.171 ***
#> X1_dyad      0.694    0.095   7.335       0    0.526    0.822 ***
#> X2_dyad      0.844    0.117   7.244       0    0.631        1 ***
#> X3_dyad       1.05    0.145   7.265       0    0.782    1.242 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: stars are a visual hint from posterior mean / SD only; for inference use the credible intervals.
#> 
#> Dynamic coefficients per period:
#> -------------------------------
#>            Mean   Min   Max Drift_pct Dynamic
#> intercept 0.120 0.120 0.120     0.000       N
#> X1_dyad   0.694 0.665 0.712     6.688       Y
#> X2_dyad   0.844 0.795 0.875     9.454       Y
#> X3_dyad   1.050 0.984 1.079     9.065       Y
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.839    0.321
#> cab    0.437    0.193
#> vb     0.646    0.245
#> rho    0.210    0.059
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
# }
```
