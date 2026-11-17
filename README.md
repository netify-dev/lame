# **lame** <img src="man/figures/lame_hex.png" align="right" alt="hex" width="200px">

<!-- badges: start -->
[![build](https://github.com/netify-dev/lame/actions/workflows/build-release-binaries.yaml/badge.svg)](https://github.com/netify-dev/lame/actions/workflows/build-release-binaries.yaml)
<!-- badges: end -->

> **L**ongitudinal **A**dditive and **M**ultiplicative **E**ffects Models for Networks

Network ties -- trade flows, friendships, alliances, sanctions -- are not independent observations. Some actors send a lot of ties, some receive a lot, ties tend to be reciprocated, and similar actors behave similarly. Run a plain regression on every (sender, receiver) pair and both your coefficients and your standard errors come out wrong.

`lame` fits models that account for that structure. Every actor gets a *sender* effect (how active they are), a *receiver* effect (how popular they are), and a position in a low-dimensional *latent space* that soaks up patterns covariates miss -- homophily, clustering, transitivity. Covariates enter additively, just like in a regression, and the "longitudinal" part lets each piece drift across time.

The package implements the additive and multiplicative effects (AME) framework developed by Peter Hoff (and popularised by his [`amen`](https://CRAN.R-project.org/package=amen) package), and extends it in a few directions: `lame()` fits panel networks with time-varying effects; both `ame()` and `lame()` handle bipartite (rectangular) networks alongside the usual square case; the sampler is written in C++ (Rcpp / RcppArmadillo) to scale; and a fit exposes the usual `coef()`, `vcov()`, `confint()`, `predict()`, `fitted()`, `residuals()`, and `simulate()` methods, so it slots into the R modelling ecosystem like `lm()` or `glm()`.

There is also a fast, **MCMC-free empirical-Bayes estimator** that reports a standard error for every estimate -- `ame_als()` / `lame_als()`, adapted from the Social Influence Regression estimator of Minhas & Hoff (2025). Reach for it when exploring, screening models, generating starting values, or fitting large networks; use the MCMC path for a full posterior and for the ordinal, censored-binary and fixed-rank families.

New to *probit*, *latent space*, or *AR(1)*? `vignette("lame")` introduces each in context before the first fit.

## Installation

```r
install.packages("lame")   # CRAN release
```

The development version compiles C++, so it needs build tools (Xcode CLI on macOS, Rtools on Windows, `build-essential` on Linux) and takes a few minutes to install:

```r
# install.packages("devtools")
devtools::install_github("netify-dev/lame", dependencies = TRUE, build_vignettes = TRUE)
```

## Quick Start

```r
library(lame)

# lazegalaw: 71 lawyers in a law firm; who asks whom for advice. ?lazegalaw
data("lazegalaw")
Y     <- lazegalaw$Y[, , "advice"]                 # 71 x 71 0/1 sociomatrix
Xdyad <- lazegalaw$Y[, , "cowork", drop = FALSE]   # dyadic covariate: do i and j work together
Xnode <- lazegalaw$X[, c("seniority", "female", "practice")]

# two chains at the default run length, about two minutes
fit <- ame(
  Y, Xdyad = Xdyad,
  Xrow    = Xnode,         # sender covariates
  Xcol    = Xnode,         # receiver covariates
  family  = "binary",
  R       = 2,             # latent-space dimension (see below)
  verbose = FALSE
)

summary(fit)               # coefficients, variance components, chain agreement
trace_plot(fit)            # check mixing; gof_plot(fit) for posterior-predictive fit
```

**What is `R`?** The dimension of the latent space that captures who-ties-to-whom beyond the covariates and sender/receiver effects. `R = 0` is the additive part only; `R = 2` adds structure like homophily and transitivity. Start at `R = 2`, then compare `gof_plot()` across `R = 0, 1, 2, 3` and keep the smallest rank that reproduces the statistics you care about.

**Actor names.** Actors are matched across `Y`, `Xdyad`, `Xrow`, `Xcol` by name, so give every input the same dimnames. Output is stored in sorted name order (`"N1","N10","N11",...` before `"N2"`); zero-padded names (`sprintf("N%02d", 1:25)`) or real names keep that order readable.

### Longitudinal

Pass a list of `Y` matrices (and lists of covariates) to `lame()`:

```r
data("dutchcollege")       # 32 students, 7 waves of friendship nominations
actors <- sprintf("S%02d", seq_len(nrow(dutchcollege$Y)))
Y <- lapply(2:7, function(t) {          # binarise; wave 1 has almost no ties
  Yt <- (dutchcollege$Y[, , t] > 0) * 1
  diag(Yt) <- NA
  dimnames(Yt) <- list(actors, actors)
  Yt
})
X_node <- dutchcollege$X
rownames(X_node) <- actors
Xrow <- Xcol <- lapply(Y, function(y) X_node)   # sender and receiver covariates

fit_long <- lame(
  Y, Xrow = Xrow, Xcol = Xcol,
  family = "binary", R = 2,
  burn = 500, nscan = 4000, odens = 10, verbose = FALSE
)
summary(fit_long)
```

### Dynamic effects

The headline extension is letting effects vary over time:

```r
fit_dynamic <- lame(
  Y, Xrow = Xrow, Xcol = Xcol,
  family = "binary", R = 2,
  dynamic_ab = TRUE,         # time-varying sender/receiver effects
  dynamic_uv = TRUE,         # time-varying latent positions
  burn = 500, nscan = 4000, odens = 10,
  verbose = FALSE
)
```

`predict(fit_dynamic, h = 3, type = "response")` forecasts three periods ahead; add `interval = "credible"` for per-period `$lower` / `$median` / `$upper`. See `vignette("forecasting")` and, for the `dynamic_beta_kind` decision tree, the [dynamic effects article](https://netify-dev.github.io/lame/articles/dynamic_effects.html).

## Families

Pick the family that matches your tie type. Pass `symmetric = TRUE` for undirected networks and `mode = "bipartite"` for two-mode data. This table is the full support map: every family runs in `ame()` and `lame()`, and three also run in the fast estimator.

| family  | tie type               | networks                       | estimators |
|---------|------------------------|--------------------------------|------------|
| normal  | continuous             | directed, symmetric, bipartite | MCMC, ALS  |
| binary  | 0/1                    | directed, symmetric, bipartite | MCMC, ALS  |
| poisson | counts                 | directed, symmetric, bipartite | MCMC, ALS  |
| ordinal | ordered categories     | directed, symmetric, bipartite | MCMC       |
| cbin    | censored binary        | directed, bipartite            | MCMC       |
| frn     | fixed rank nominations | directed, bipartite            | MCMC       |

`cbin` and `frn` model each actor's own nominations, so they are directed by construction. The fast estimator warns on a binary network whose density (the share of observed cells that are ties) lies outside 0.2 to 0.8, where the MCMC fit is the better choice.

## Coming from another package

**Tidy edgelists, panels, bipartite data:** build the object with [`netify`](https://netify-dev.github.io/netify/) and hand it straight to `ame()` / `lame()`:

```r
net <- netify::netify(edges, actor1 = "sender", actor2 = "receiver", time = "year",
                      weight = "tie", dyad_vars = c("distance", "shared_border"),
                      nodal_vars = c("gdp", "population"), symmetric = FALSE)
fit <- lame(net, family = "binary", R = 2, burn = 500, nscan = 4000, odens = 10, verbose = FALSE)
```

For a lone `igraph` / `network` / adjacency matrix, `as_lame_y(g)` returns a plain `Y`.

**From ERGM:** the two model classes answer different questions. ERGM specifies a joint distribution over the whole graph via change statistics (`gwesp`, `triangle`); AME treats dyads as *conditionally independent* given the actor effects, latent positions, and covariates, and captures higher-order structure through those random effects. Use ERGM when a change statistic such as `gwesp` is the estimand; use AME for actor heterogeneity, dynamics, and forecasts, with estimation that does not degenerate. `vignette("cross_sec_ame")` scores a triangle-over-two-path clustering statistic on a fit through `gof(fit, custom_gof = ...)`. The ERGM-style helpers `nodematch()`, `absdiff()`, and `nodefactor()` build dyadic covariates that drop straight into `Xdyad`.

## Visualization

```r
ab_plot(fit, effect = "sender")                 # additive effects; plot_type = "trajectory" over time
uv_plot(fit)                                     # latent positions;  plot_type = "trajectory" over time
trace_plot(fit); gof_plot(fit)                   # mcmc mixing; goodness of fit
```

## Features

- **Modes:** cross-sectional (`ame()`), longitudinal (`lame()`), unipartite and bipartite, six families.
- **Dynamics:** time-varying latent positions (`dynamic_uv`), sender/receiver effects (`dynamic_ab`), and regression coefficients (`dynamic_beta`, per-block AR(1) / RW1 / RW2 / Matern 3/2), with forecasting via `predict(fit, h = K)` (AR(1) and RW1 coefficient paths).
- **S3 methods** on both MCMC and ALS fits: `coef()`, `vcov()`, `confint()`, `predict()`, `fitted()`, `residuals()`, `summary()`, `update()`, `simulate()`, plus `autoplot()`.
- **Ecosystem:** `broom::tidy()` / `glance()`, `modelsummary`, `as_draws()` for `posterior` / `tidybayes`, and `loo::loo()` when `save_log_lik = TRUE` (directed fits need `dcor = FALSE`, so that each cell has its own likelihood).
- **Fit assessment:** `simulate()` and `gof()` posterior-predictive checks; `gof_temporal()` for panels.
- **Convergence and scale:** `ame()` and `lame()` run two chains by default (`chains = 2`) and compare them by R-hat before they return; `compute_mcmc_diagnostics()` reports R-hat / ESS; `lame_multi()` for several panels sharing coefficients; checkpoint/resume via `checkpoint_path=` and `lame_resume()`.

## Fast (MCMC-free) estimation

`ame_als()` / `lame_als()` fit the same model by iterative block coordinate descent (alternating least squares / IRLS) -- a port of the Social Influence Regression estimator of Minhas & Hoff (2025), implemented in the [`sir`](https://github.com/netify-dev/sir) package. They cover the `normal`, `binary` and `poisson` families, and the fit is an empirical-Bayes / variational estimator of the model the samplers fit (see `vignette("fast_estimation")`). They return point estimates with standard errors rather than posterior draws.

- **Binary density:** a binary network or panel whose density is below 0.2 or above 0.8 is fitted with a warning; there the sampler (`ame()` or `lame()`) is the better choice.
- **Dynamic fits:** `lame(..., method = "als")` fits dynamic additive effects, selected dynamic coefficients, AR(1) latent factors (directed, symmetric, bipartite) and Student-t ones (directed, symmetric). Normal, binary and Poisson panels are fitted by a dynamic empirical-Bayes estimator, with standard errors for the coefficients, the effects and their paths, the variance components and the product cells. Tune with `als_max_iter`, `als_tol`, `als_stability`. The snap-shift kind (`dynamic_uv_kind = "snap"`) is fitted by the MCMC on directed panels.
- **Uncertainty:** `vcov()` / `confint()` report Laplace standard errors for every coefficient, and every fit carries standard errors for its variance components, sender and receiver effects and latent positions. On static fits, `ame_als_bootstrap()`, or `ame_als(..., bootstrap = N)`, adds block/parametric bootstrap intervals, and `confint()` uses them when present.

## Documentation

| Vignette or article | Topic |
| --- | --- |
| `vignette("lame")` | 5-minute getting-started on a small simulated panel |
| [lame overview](https://netify-dev.github.io/lame/articles/lame-overview.html) | Full tour on the bundled Dutch college friendship data |
| `vignette("cross_sec_ame")` | Cross-sectional AME on the Add Health friendship network |
| [Bipartite article](https://netify-dev.github.io/lame/articles/bipartite.html) | Two-mode networks, cross-sectional and longitudinal |
| [Dynamic effects article](https://netify-dev.github.io/lame/articles/dynamic_effects.html) | `dynamic_uv` / `dynamic_ab` / `dynamic_beta` and the transition-kind decision tree |
| `vignette("forecasting")` | `predict(fit, h = K)` forecasts and counterfactuals |
| `vignette("fast_estimation")` | the `ame_als()` / `lame_als()` empirical-Bayes estimator |

## Reproducibility

`ame()` and `lame()` seed their chains from the `seed` argument (default 6886), so the same call returns the same draws whatever `set.seed()` ran before it. Pass another `seed =` for an independent run, and record the environment:

```r
fit <- lame(..., seed = 6886)
sessionInfo()          # or pin the toolchain with an renv.lock / a @<commit-sha> install
```

## Citation

```bibtex
@Manual{lame2026,
  title  = {lame: Longitudinal Additive and Multiplicative Effects Models for Networks},
  author = {Cassy Dorff and Shahryar Minhas and Tosin Salau},
  year   = {2026},
  note   = {R package version 1.4.2},
  url    = {https://github.com/netify-dev/lame},
}
```

The dynamic effects draw on Sewell & Chen (2015, [doi:10.1080/01621459.2014.988214](https://doi.org/10.1080/01621459.2014.988214)) and Durante & Dunson (2014, [doi:10.1093/biomet/asu040](https://doi.org/10.1093/biomet/asu040)); the fast estimator adapts Minhas & Hoff (2025, *Political Analysis*).

## Contributors

Cassy Dorff (Vanderbilt), Shahryar Minhas (Michigan State), Tosin Salau (Michigan State).

## License

MIT. See the `LICENSE` file.

## Support

[Issues](https://github.com/netify-dev/lame/issues) · [minhassh@msu.edu](mailto:minhassh@msu.edu)
