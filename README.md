# **lame** <img src="man/figures/lame_hex.png" align="right" alt="hex" width="200px">

<!-- badges: start -->
[![R-CMD-check](https://github.com/netify-dev/lame/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/netify-dev/lame/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> **L**ongitudinal **A**dditive and **M**ultiplicative **E**ffects Models for Networks

Trade flows, friendship nominations, alliance ties, and sanctions are all examples of *network* data: each observation is a tie between two actors (countries, people, organisations), not an independent unit. That non-independence breaks the assumptions a standard regression makes. Some actors send a lot of ties, some receive a lot, ties tend to be reciprocated, and actors with similar attributes tend to behave similarly. If you ignore those patterns and just run a logistic regression on every (sender, receiver) pair, both your coefficients and your standard errors will be off.

The `lame` package fits **L**ongitudinal **A**dditive and **M**ultiplicative **E**ffects models that account for that structure. Concretely, the model gives every actor two random effects (a *sender* effect $a_i$ for how active they are and a *receiver* effect $b_j$ for how popular they are), plus a position $(u_i, v_j)$ in a low-dimensional *latent space* that captures patterns covariates cannot, such as homophily ("similar actors tie to each other") and transitivity ("friends of friends are friends"). Covariates enter additively, just like in a regression. The "longitudinal" piece lets all of these pieces drift across time periods.

The package builds directly on Peter Hoff's [`amen`](https://CRAN.R-project.org/package=amen) and extends it in four directions: (1) `lame()` fits *panel* networks (the same actors observed at multiple time points) with time-varying effects; (2) it supports *bipartite* (rectangular) networks — such as countries-by-international-organisations (which states belong to which IGOs) — alongside the more familiar unipartite (square) case; (3) the MCMC sampler is written in C++ (via Rcpp / RcppArmadillo) so it scales to the network sizes typical in political-science applications; (4) the fitted object exposes the usual `coef()`, `vcov()`, `confint()`, `predict()`, `fitted()`, `residuals()`, `simulate()` methods, so it plugs into the rest of the R modelling ecosystem the same way `lm()` or `glm()` does.

**Two ways to estimate the same model.** Alongside the Bayesian MCMC sampler, `lame` ships a fast, **MCMC-free point estimator**, `ame_als()` / `lame_als()`, adapted from the **Social Influence Regression (SIR)** estimator of Minhas & Hoff (2025). It fits the same additive-and-multiplicative structure by iterative block coordinate descent (alternating least squares / IRLS) for the normal, binary, and Poisson families, with parametric-bootstrap or sandwich standard errors and the usual S3 methods. Use it for data exploration, model screening, starting values, or large networks; use the MCMC path (`ame()` / `lame()`) when you need a full posterior or the rank/censored families. For longitudinal normal, binary, and Poisson panels, including named panels where actors enter or exit, `lame(..., method = "als")` has a dynamic point-estimation path for time-varying additive effects, selected regression coefficients, AR(1) latent factors in directed, symmetric, and bipartite panels, Student-t latent-factor drift in directed, symmetric, and bipartite panels, and bipartite dynamic `G`. For normal unipartite and bipartite snap-shift panels, `lame_snap_als()` gives a fast ALS snap score using sequential alignment and lower-tail drift updates by default.

If terms like *probit*, *latent space*, or *AR(1)* are new to you, the `vignette("lame")` walk-through introduces each one in context before the first fit.

Five estimators form the core public API. `ame()` and `lame()` fit cross-sectional and longitudinal networks via Bayesian MCMC; `lame()` additionally adds dynamic additive, multiplicative, and (optionally) regression-coefficient effects through autoregressive processes. `ame_als()` and `lame_als()` are their fast, MCMC-free static counterparts (see [Fast (MCMC-free) estimation](#fast-mcmc-free-estimation)), while `lame(..., method = "als")` routes supported dynamic requests to the dynamic ALS point-estimation path. `lame_snap_als()` is a separate fast approximate estimator for dynamic snap-shift latent positions in normal unipartite and bipartite panels.

## Installation

Install from GitHub (requires C++ build tools: Xcode CLI on macOS, Rtools on Windows, or `build-essential` on Linux):

```r
# install.packages("devtools")
devtools::install_github(
	"netify-dev/lame",
	dependencies    = TRUE,
	build_vignettes = TRUE
)
```

`build_vignettes = TRUE` is recommended: the package ships long-form
vignettes (`vignette(package = "lame")` after install), which are
built locally during install. The first install may take a few minutes while
the `Rcpp` / `RcppArmadillo` C++ code compiles.

## Quick Start

```r
library(lame)

# load example data
# `YX_bin` is a synthetic dataset bundled with the package containing a
# 100-actor binary network and 8 dyadic covariates. `data("YX_bin")`
# puts the list `YX_bin` into your environment; run `?YX_bin` for the
# data dictionary, and inspect it before fitting:
data("YX_bin")
str(YX_bin, max.level = 1)        # list of 2: Y (100x100), X (100x100x8)
Y      <- YX_bin$Y                # 100 x 100 binary "sociomatrix":
                                  #   rows = senders, cols = receivers,
                                  #   entries 0/1, diagonal NA (self-ties undefined)
# the bundled YX_bin$X array carries a literal "intercept" slice
# as its first slice, because the package's worked-example tradition is
# to include it explicitly. `ame()` adds its own intercept by default
# (the `intercept = TRUE` argument), so we drop the redundant slice here
# to avoid a rank-deficient design matrix:
Xdyad  <- YX_bin$X[, , -1]        # 100 x 100 x 7 array of dyadic covariates
# the third dim of a 3-D Xdyad array names the covariates (one slice each):
dimnames(Xdyad)[[3]]              # prints: "rgpa" "rsmoke" "cgpa" "csmoke" ...

# fit a cross-sectional ame model on one snapshot
# defaults below are sized for a real fit: nscan = 4000
# with odens = 10 gives ~400 stored draws, which is enough to read trace
# plots and credible intervals. Drop to nscan = 500 only if you just want
# to see the function run.
fit <- ame(
  Y       = Y,            # n x n network matrix
  Xdyad   = Xdyad,        # dyadic covariates (n x n x p array)
  family  = "binary",     # outcome family
  R       = 2,            # latent dimensions (multiplicative effects, see ?ame)
  burn    = 500,          # burn-in iterations, discarded
  nscan   = 4000,         # post-burn iterations
  odens   = 10,           # thinning: keep every odens-th iter, about 400 draws
  verbose = FALSE         # quiet fit; drop this for a live progress bar
)

summary(fit)              # coefficients, variance components, gof
# check convergence on a fresh fit:
trace_plot(fit)           # one panel per parameter
# gof_plot(fit)           # posterior-predictive GOF
```

**What is `R`?** `R` is the dimension of the multiplicative-effects
latent space, a small space where each actor gets coordinates that
say "who I tend to tie to" beyond what the covariates and sender /
receiver effects explain. `R = 0` fits the additive part only
(sender + receiver + dyadic correlation); `R = 2` adds a 2-D latent
space that captures higher-order structure like *homophily* (similar
actors tying to each other), clustering, and transitivity (friends
of friends are friends). The package warns if `R > floor(n/3)`
because high-rank latent factors can absorb structure that belongs
to the additive effects. Start with `R = 2` for a first fit; try
`R = 0, 1, 2, 3` and compare `gof_plot()` panels to pick a value.
See `?ame` for math and guidance.

**A naming gotcha (the one most beginners hit).** `ame()` / `lame()`
align actors across `Y`, `Xdyad`, `Xrow`, and `Xcol` by **sorted**
dimnames. The natural-looking idiom

```r
rownames(Y) <- colnames(Y) <- paste0("N", 1:25)  # "N1","N2",...,"N10",...,"N25"
```

silently breaks this, because the alphabetic sort gives
`"N1","N10","N11",...,"N19","N2","N20",...`, not `"N1","N2",...,"N10",...`.
Any positionally-aligned (but unnamed) covariate array then gets
matched against that re-ordered `Y` and your `X` rows no longer
correspond to your `Y` rows. Two safe patterns:

```r
rownames(Y) <- colnames(Y) <- sprintf("N%02d", 1:25)  # zero-padded: sort matches position
# or just give actors real names:
rownames(Y) <- colnames(Y) <- c("Alice", "Bob", ...)
```

The package will inherit `Y`'s dimnames onto unnamed covariate arrays
when it safely can, but explicitly matching dimnames on `Y`, `Xdyad`,
`Xrow`, and `Xcol` is always the safest option.

### Longitudinal Quick Start

For a panel (sequence of network snapshots), use `lame()` and pass a list of
`Y` matrices (and lists of covariates):

```r
# `data("vignette_data")` brings four objects into your environment:
#   Y, Xdyad, Xrow, Xcol  -- each a length-4 list (one entry per year),
# covering 35 countries observed over 4 years of sanctions data.
# See `?vignette_data`.
data("vignette_data")     # Y is now a list of 4 binary 35x35 matrices

fit_long <- lame(
  Y      = Y,             # List of T matrices
  Xdyad  = Xdyad,         # List of T dyadic-covariate arrays
  Xrow   = Xrow,          # List of T sender-covariate matrices
  Xcol   = Xcol,          # List of T receiver-covariate matrices
  family = "binary", R = 2,
  burn = 500, nscan = 4000, odens = 10,
  verbose = FALSE
)

summary(fit_long)
```

### Coming from another network package?

`lame` expects `Y` as a (possibly named) matrix, not a graph object. The
package ships a one-call `as_lame_y()` helper that handles `igraph`,
`network`, plain matrices, and data.frames with the correct NA diagonal:

```r
# igraph / network / statnet / matrix -> matrix
Y <- as_lame_y(g)                              # works for any of the above

# manual idioms still work
Y <- igraph::as_adjacency_matrix(g, sparse = FALSE); diag(Y) <- NA
Y <- as.matrix(net); diag(Y) <- NA
```

### Coming from a long-format tibble?

If your data lives in a tibble of (sender, receiver, tie), pivot to a matrix
and align row / column order explicitly:

```r
library(tidyr); library(tibble); library(dplyr)
# get a stable actor universe (sorted; or whatever order you want)
actors <- sort(unique(c(df$sender, df$receiver)))
Y <- df |>
  # ensure every (sender, receiver) cell exists, with NA for absent pairs
  complete(sender = actors, receiver = actors) |>
  arrange(factor(sender, levels = actors)) |>
  pivot_wider(names_from = receiver, values_from = tie) |>
  column_to_rownames("sender") |>
  as.matrix()
Y <- Y[actors, actors]    # final guarantee: row / col order matches
diag(Y) <- NA              # unipartite only; skip for bipartite
```

A common mistake: `pivot_wider()` does not preserve row order, so the naive
`pivot_wider() |> column_to_rownames() |> as.matrix()` can produce a `Y`
where `rownames(Y) != colnames(Y)`, and the next `diag(Y) <- NA` writes NA
into the wrong cells. The `arrange()` + `Y[actors, actors]` step above
prevents that.

### Coming from ERGM?

`lame` is not an `ergm` substitute and the two models make
different dependence assumptions. ERGM specifies a joint distribution
over the whole adjacency matrix and uses terms like `gwesp`, `triangle`,
and `kstar` (via *change statistics*) to encode unconditional
higher-order dependence among ties. AME instead assumes **conditional
dyadic independence**: given the sender effects `a_i`, receiver effects
`b_j`, latent positions `u_i, v_j`, dyadic correlation `rho`, and
covariates, the dyads are independent. Higher-order structure
(clustering, transitivity, degree heterogeneity) is captured *indirectly*
through these random effects when you integrate them out. Practical
consequences:

- AME does not suffer the model-degeneracy failure modes of
  triangle-heavy ERGMs.
- AME does not give you change-statistic interpretations of beta --
  coefficients are partial associations on the latent (e.g. probit)
  scale.
- AME typically under-predicts transitivity in friendship-style networks
  (see the GOF sections in `vignette("cross_sec_ame")` and
  `vignette("lame-overview")`). If `gwesp` is your substantive target,
  stay in ERGM; if you want stable estimates with explicit actor
  heterogeneity, dynamic effects, or forecasts, AME is the better tool.
- The GOF section in `vignette("cross_sec_ame")` shows how to score
  ERGM-style raw `triangle`, two-path, and clustering-coefficient
  statistics via `gof(fit, custom_gof = ...)` so you can read AME
  against the same diagnostics you would put under `ergm::gof()`.

### ERGM-style dyadic covariates

`lame` exports `nodematch()`, `absdiff()`, and `nodefactor()` -- direct
analogues of the ERGM terms of the same name -- that drop straight into
`Xdyad`:

```r
Xdyad <- array(
  c(nodematch(group), absdiff(age)),
  dim = c(n, n, 2),
  dimnames = list(NULL, NULL, c("same_group", "age_diff"))
)
```

For a single cross-section use `ame(Y, ...)`; for a panel use
`lame(list(Y_t1, Y_t2, ...), ...)`. For undirected networks pass
`symmetric = TRUE`; for two-mode/rectangular data (countries ×
international organisations, states × treaties) pass
`mode = "bipartite"`. `Xrow` / `Xcol` accept
either a numeric matrix or a `data.frame`; `Xdyad` must be a 3-D
`n × n × p` numeric array.

All eight families are supported under both modes; pick the family
that matches your tie type (continuous, binary, ordered categories,
counts, censored, rank lists). For symmetric (undirected) networks
pass `symmetric = TRUE` and the same family table applies.

| family  | unipartite | bipartite | symmetric |
|---------|------------|-----------|-----------|
| normal  | yes        | yes       | yes       |
| binary  | yes        | yes       | yes       |
| tobit   | yes        | yes       | yes       |
| ordinal | yes        | yes       | yes       |
| poisson | yes        | yes       | yes       |
| cbin    | yes        | yes       | no (row-cone is directed) |
| frn     | yes        | yes       | no (row-cone is directed) |
| rrl     | yes        | yes       | no (row-cone is directed) |

Every family works under both modes, including the rank-and-count
families (`ordinal`, `cbin`, `frn`, `rrl`, `poisson`) on bipartite
networks and symmetric `ordinal`; the
[`vignette("bipartite")`](#documentation) walk-through covers each.

### Dynamic Effects

The key extension in `lame` is time-varying network effects:

```r
# fit model with dynamic effects, building on fit_long and the panel above
fit_dynamic <- lame(
  Y = Y,
  Xdyad = Xdyad,
  Xrow = Xrow,
  Xcol = Xcol,
  family = "binary",
  dynamic_ab = TRUE,        # time-varying sender/receiver effects
  dynamic_uv = TRUE,        # time-varying latent positions
  dynamic_beta = "dyad",    # time-varying regression coefficients on dyadic Xs
  dynamic_beta_kind = "ar1",# "ar1" / "rw1" / "rw2" / "matern32"
  R = 2,
  burn = 500, nscan = 4000, odens = 10,
  prior = list(
    rho_uv_mean = 0.9,
    rho_ab_mean = 0.8
  ),
  verbose = FALSE
)
```

`predict(fit_dynamic, h = 3, type = "response")` then forecasts three
periods ahead by drawing the AR(1) / RW1 hyperparameters from the
posterior; `predict(..., interval = "credible")` returns the per-period
`$lower` / `$median` / `$upper` matrices. See
`vignette("forecasting")` for the full forecasting workflow and
`vignette("dynamic_effects")` for the `dynamic_beta_kind` decision tree.

## Visualization

```r
# additive effects (sender/receiver)
ab_plot(fit, effect = "sender")                    # static effects
ab_plot(fit_dynamic, plot_type = "trajectory")     # dynamic over time

# multiplicative effects (latent factors)
uv_plot(fit)                                       # static latent positions
uv_plot(fit_dynamic, plot_type = "trajectory")     # dynamic trajectories

# diagnostics
trace_plot(fit)                                     # mcmc convergence
gof_plot(fit)                                       # goodness of fit
```

## Key Features

### Network Analysis

- **Cross-sectional**: `ame()` fits AME models for a single network snapshot
- **Longitudinal**: `lame()` fits AME models for networks observed over multiple time periods
- **Bipartite**: full support for two-mode (rectangular) networks in both `ame()` and `lame()`
- **8 families**: normal, binary, ordinal, poisson, tobit, censored binary, fixed rank nomination, row-ranked likelihood

### Dynamic Modeling

- **Time-varying latent positions** (`dynamic_uv`): captures evolving community structure via AR(1) processes
- **Time-varying heterogeneity** (`dynamic_ab`): models changing sender/receiver activity over time
- **Time-varying regression coefficients** (`dynamic_beta`): lets dyadic, nodal, or intercept coefficients evolve via per-block AR(1) processes; selectable per-coefficient or per-block; identifiability is preserved via an automatic sum-to-zero constraint on the additive effects
- **Flexible priors**: customizable temporal persistence and innovation variance

### Inference and Diagnostics

- **S3 methods**: `coef()`, `vcov()`, `confint()`, `residuals()`,
  `predict()`, `summary()`, `update()`, `simulate()`. All work on
  both MCMC fits (`ame` / `lame`) and ALS fits (`ame_als` / `lame_als`).
- **Visualization**: `trace_plot()`, `gof_plot()`, `uv_plot()`,
  `ab_plot()`. `autoplot()` returns the coefficient-with-CI plot for
  either fit class.
- **Broom + modelsummary**: `tidy()` and `glance()` are registered
  against `generics::tidy` / `generics::glance`, so `broom::tidy(fit)`,
  `broom::glance(fit)`, and `modelsummary::modelsummary(list(...))`
  work directly on both MCMC and ALS fits. `tidy.boot_ame()` exposes
  the standalone `ame_als_bootstrap()` return path so it composes
  with `modelsummary` the same way.
- **Posterior tools**: `as_draws(fit)` for the `posterior` /
  `tidybayes` ecosystem; `loo::loo(fit)` for leave-one-out CV when
  `save_log_lik = TRUE` was set at fit time.
- **Simulation and GOF**: `simulate(fit, nsim = K)` and `gof(fit)`
  for posterior predictive checks; `gof_temporal()` adds a temporal
  PP test for longitudinal fits.
- **Multi-chain**: `ame_parallel()` / `lame_parallel()` run multiple
  chains in parallel and pool or return them as a list, with built-in
  R-hat / ESS diagnostics via `compute_mcmc_diagnostics()`.
- **Multi-panel**: `lame_multi()` fits K panels with shared regression
  coefficients (e.g., cooperation + conflict + trade networks on the
  same actors).
- **Checkpoint and resume**: `lame(..., checkpoint_path = "X.rds",
  max_seconds = T)` writes a checkpoint when the wall-clock budget
  expires; resume with either `lame(resume_from = "X.rds", nscan = K)`
  or `lame_resume("X.rds", nscan_more = K)`.

### Fast (MCMC-free) estimation

`ame_als()` and `lame_als()` are a fast, MCMC-free point estimator for
the additive-and-multiplicative model, **adapted from the Social
Influence Regression (SIR) estimator of Minhas & Hoff (2025)** — the
iterative block coordinate descent method introduced in *Decomposing
Network Dynamics: Social Influence Regression* and implemented in the
[`sir`](https://github.com/netify-dev/sir) package. It is a port and
adaptation of that estimator to the AME model, not original `lame`
methodology.

- **`ame_als()` / `lame_als()`**: fast AME point estimates via
  iterative block coordinate descent (alternating least squares /
  IRLS), typically much faster than the MCMC sampler since they return
  point estimates rather than a posterior sample. Useful for data
  exploration, model screening, starting values, and large networks.
  **Family scope**: `normal`, `binary`, `poisson`; for `tobit`,
  `ordinal`, `cbin`, `frn`, `rrl` use the MCMC path.
- **`lame(..., method = "als")` with dynamic effects**: for
  normal, binary, and Poisson panels, the dispatcher
  can fit dynamic additive effects, selected dynamic regression
  coefficients, AR(1) and t-transition dynamic `U`/`V` in directed,
  symmetric, and bipartite panels, and bipartite dynamic `G`. These are
  penalized point estimates. Student-t dynamic UV fits attach final local
  transition-weight matrices as `fit$lambda_u` and `fit$lambda_v`.
  Use `als_max_iter`, `als_tol`, and `als_stability` when a dynamic ALS
  fit needs a longer run or start-sensitivity check.
  Node-covariate coefficients use the same orthogonal additive-effect
  decomposition as `lame_als()`; dynamic node coefficients use
  period-specific node values when selected by `dynamic_beta`, while static
  node coefficients use per-actor means. Named changing-composition panels are
  aligned to the union actor set, with smoothing penalties broken across
  actor-entry gaps. Rank/censored families still use the MCMC path. Static ALS fits still use a
  single latent rank `R`; dynamic bipartite ALS honors separate `R_row` and
  `R_col` values.
- **`lame_snap_als()`**: fast dynamic snap-shift point estimator
  for normal unipartite and bipartite longitudinal networks. It returns dynamic
  `U`/`V`, `rho_uv`, `sigma_uv`, `pi_snap`, `snap_prob`, and transition
  diagnostics. In bipartite fits it uses a static `G` matrix and reports
  separate row- and column-side snap scores. The same path is available from the unified front door
  with `lame(..., method = "als", dynamic_uv = TRUE,
  dynamic_uv_kind = "snap")`. The default uses sequential initialization
  and lower-tail drift updates so broad ruptures are not absorbed
  as high-variance drift. `snap_prob` is an ALS snap score, not a
  Bayesian posterior probability. Read it through rankings, heuristic 0.5
  classifications, and the stability diagnostics; `fit$convergence$max_snap_stable` and
  `fit$convergence$unstable_transitions` show actor-periods that are still
  moving. Snap ALS does not estimate node covariates, dynamic coefficients, or
  dynamic `G`.
- **`ame_als_bootstrap()`**: block and parametric bootstrap standard
  errors and confidence intervals for the fast estimator. The
  one-shot path `ame_als(..., bootstrap = N)` attaches the bootstrap
  to the fit; `tidy(fit)` / `confint(fit)` use the bootstrap
  intervals when present, the sandwich (`vcov.ame_als`) otherwise.

## Documentation

| Vignette | Topic |
| --- | --- |
| `vignette("lame")` | 5-minute getting-started with a small simulated panel |
| `vignette("lame-overview")` | Full tour on the bundled Dutch college friendship data |
| `vignette("cross_sec_ame")` | Cross-sectional AME on the Add Health friendship network |
| `vignette("bipartite")` | Two-mode networks (cross-sectional + longitudinal) |
| `vignette("dynamic_effects")` | `dynamic_uv` / `dynamic_ab` / `dynamic_beta` with the AR(1) / RW1 / RW2 / Matérn 3/2 decision tree |
| `vignette("forecasting")` | `predict(fit, h = K)` h-step-ahead forecasts and counterfactuals |
| `vignette("fast_estimation")` | `ame_als()` / `lame_als()`, the fast MCMC-free point estimator |

## Reproducibility: version pinning and `sessionInfo()`

MCMC results in this package depend on the **RNG state**, the **C++
ABI** of `RcppArmadillo`, and minor internal changes in `lame` itself
across releases. To make a `lame` analysis reproducible for a
collaborator, a reviewer, or your future self, capture both the random
seed and the exact package versions you fit under. The conventional
pattern at the top of any replication script:

```r
set.seed(6886)        # any integer; same seed -> same chain
fit <- lame(...)
sessionInfo()         # paste output into your README / supplementary
```

`sessionInfo()` records the R version, OS, locale, and every loaded
package version, enough information for someone else to recreate the
environment with `renv::restore()` or a manual `remotes::install_version()`
sweep. For long-running fits, save the fitted object plus the session
record together:

```r
saveRDS(list(fit = fit, session = sessionInfo()),
        file = "lame_fit_2026-05-26.rds")
```

For long-term reproducibility, pin `lame` to a specific commit hash
(the package is pre-CRAN as of writing):

```r
remotes::install_github("netify-dev/lame@<commit-sha>")
```

or commit an `renv.lock` to your project repo so reviewers can
recreate the exact toolchain with `renv::restore()`. Random-seed +
package version is the floor; commit hash + `renv.lock` is the ceiling.

## Citation

If you use `lame` in your research, please cite:

```bibtex
@Manual{lame2026,
  title = {lame: Longitudinal Additive and Multiplicative Effects Models for Networks},
  author = {Cassy Dorff and Shahryar Minhas and Tosin Salau},
  year = {2026},
  note = {R package version 1.1.0},
  url = {https://github.com/netify-dev/lame},
}
```

The dynamic effects implementation draws on:

- Sewell, D. K., & Chen, Y. (2015). Latent space models for dynamic networks. *Journal of the American Statistical Association*, 110(512), 1646-1657. [doi:10.1080/01621459.2014.988214](https://doi.org/10.1080/01621459.2014.988214)
- Durante, D., & Dunson, D. B. (2014). Nonparametric Bayes dynamic modeling of relational data. *Biometrika*, 101(4), 883-898. [doi:10.1093/biomet/asu040](https://doi.org/10.1093/biomet/asu040)

The fast `ame_als()` / `lame_als()` estimator and its bootstrap are adapted from:

- Minhas, S., & Hoff, P. D. (2025). Decomposing Network Dynamics: Social Influence Regression. *Political Analysis*. (Implemented in the [`sir`](https://github.com/netify-dev/sir) package.)

## Contributors

- **Cassy Dorff** (Vanderbilt University)
- **Shahryar Minhas** (Michigan State University)
- **Tosin Salau** (Michigan State University)

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Support

- **Bug Reports**: [GitHub Issues](https://github.com/netify-dev/lame/issues)
- **Questions**: [GitHub Discussions](https://github.com/netify-dev/lame/discussions)
- **Contact**: [minhassh@msu.edu](mailto:minhassh@msu.edu)
