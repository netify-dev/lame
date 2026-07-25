# **lame** <img src="man/figures/lame_hex.png" align="right" alt="hex" width="200px">

<!-- badges: start -->
[![R-CMD-check](https://github.com/netify-dev/lame/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/netify-dev/lame/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> **L**ongitudinal **A**dditive and **M**ultiplicative **E**ffects Models for Networks

Network ties -- trade flows, friendships, alliances, sanctions -- are not independent observations. Some actors send a lot of ties, some receive a lot, ties tend to be reciprocated, and similar actors behave similarly. Run a plain regression on every (sender, receiver) pair and both your coefficients and your standard errors come out wrong.

`lame` fits models that account for that structure. Every actor gets a *sender* effect (how active they are), a *receiver* effect (how popular they are), and a position in a low-dimensional *latent space* that soaks up patterns covariates miss -- homophily, clustering, transitivity. Covariates enter additively, just like in a regression, and the "longitudinal" part lets each piece drift across time.

The package implements the additive and multiplicative effects (AME) framework developed by Peter Hoff (and popularised by his [`amen`](https://CRAN.R-project.org/package=amen) package), and extends it in a few directions: `lame()` fits panel networks with time-varying effects; both `ame()` and `lame()` handle bipartite (rectangular) networks alongside the usual square case; the sampler is written in C++ (Rcpp / RcppArmadillo) to scale; and a fit exposes the usual `coef()`, `vcov()`, `confint()`, `predict()`, `fitted()`, `residuals()`, and `simulate()` methods, so it slots into the R modelling ecosystem like `lm()` or `glm()`.

There is also a fast, **MCMC-free point estimator** -- `ame_als()` / `lame_als()`, adapted from the Social Influence Regression estimator of Minhas & Hoff (2025). Reach for it when exploring, screening models, generating starting values, or fitting large networks; use the MCMC path when you need a full posterior or the rank/censored families.

New to *probit*, *latent space*, or *AR(1)*? `vignette("lame")` introduces each in context before the first fit.

## Installation

Needs C++ build tools (Xcode CLI on macOS, Rtools on Windows, `build-essential` on Linux). The first install takes a few minutes while the C++ compiles.

```r
# install.packages("devtools")
devtools::install_github("netify-dev/lame", dependencies = TRUE, build_vignettes = TRUE)
```

## Quick Start

```r
library(lame)

# YX_bin: a bundled 100-actor binary network (Y) with dyadic covariates (X). ?YX_bin
data("YX_bin")
Y     <- YX_bin$Y          # 100 x 100 0/1 sociomatrix, diagonal NA
Xdyad <- YX_bin$X[, , -1]  # drop the explicit intercept slice; ame() adds its own

fit <- ame(
  Y, Xdyad = Xdyad,
  family  = "binary",
  R       = 2,             # latent-space dimension (see below)
  burn    = 500,
  nscan   = 4000,          # ~400 stored draws with odens = 10
  odens   = 10,
  verbose = FALSE
)

summary(fit)               # coefficients, variance components, gof
trace_plot(fit)            # check mixing; gof_plot(fit) for posterior-predictive fit
```

**What is `R`?** The dimension of the latent space that captures who-ties-to-whom beyond the covariates and sender/receiver effects. `R = 0` is the additive part only; `R = 2` adds structure like homophily and transitivity. `lame` warns above `R > floor(n/3)`, where high-rank factors start stealing signal from the additive effects. Start at `R = 2`, then compare `gof_plot()` across `R = 0, 1, 2, 3`.

**Naming gotcha.** Actors are aligned across `Y`, `Xdyad`, `Xrow`, `Xcol` by **sorted** dimnames, so `paste0("N", 1:25)` misaligns everything (`"N1","N10","N11",...` sorts before `"N2"`). Zero-pad (`sprintf("N%02d", 1:25)`) or use real names, and set matching dimnames on every input.

### Longitudinal

Pass a list of `Y` matrices (and lists of covariates) to `lame()`:

```r
data("vignette_data")      # Y: list of 4 binary 35 x 35 matrices (sanctions, 35 countries x 4 yrs)

fit_long <- lame(
  Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
  family = "binary", R = 2,
  burn = 500, nscan = 4000, odens = 10, verbose = FALSE
)
summary(fit_long)
```

### Dynamic effects

The headline extension is letting effects vary over time:

```r
fit_dynamic <- lame(
  Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
  family = "binary", R = 2,
  dynamic_ab   = TRUE,       # time-varying sender/receiver effects
  dynamic_uv   = TRUE,       # time-varying latent positions
  dynamic_beta = "dyad",     # time-varying dyadic coefficients
  dynamic_beta_kind = "ar1", # "ar1" / "rw1" / "rw2" / "matern32"
  burn = 500, nscan = 4000, odens = 10,
  prior = list(rho_uv_mean = 0.9, rho_ab_mean = 0.8),
  verbose = FALSE
)
```

`predict(fit_dynamic, h = 3, type = "response")` forecasts three periods ahead; add `interval = "credible"` for per-period `$lower` / `$median` / `$upper`. See `vignette("forecasting")` and, for the `dynamic_beta_kind` decision tree, `vignette("dynamic_effects")`.

## Families

Pick the family that matches your tie type. All are supported unipartite, bipartite, and (where the tie type allows) symmetric -- pass `symmetric = TRUE` for undirected networks and `mode = "bipartite"` for two-mode data.

| family  | tie type            | unipartite | bipartite | symmetric |
|---------|---------------------|------------|-----------|-----------|
| normal  | continuous          | yes        | yes       | yes       |
| binary  | 0/1                 | yes        | yes       | yes       |
| ordinal | ordered categories  | yes        | yes       | yes       |
| poisson | counts              | yes        | yes       | yes       |
| cbin    | censored binary     | yes        | yes       | no (directed row-cone) |
| frn     | fixed rank nominations | yes     | yes       | no (directed row-cone) |

## Coming from another package

**Tidy edgelists, panels, bipartite data:** build the object with [`netify`](https://netify-dev.github.io/netify/) and hand it straight to `ame()` / `lame()`:

```r
net <- netify::netify(edges, actor1 = "sender", actor2 = "receiver", time = "year",
                      weight = "tie", dyad_vars = c("distance", "shared_border"),
                      nodal_vars = c("gdp", "population"), symmetric = FALSE)
fit <- lame(net, family = "binary", R = 2, burn = 500, nscan = 4000, odens = 10, verbose = FALSE)
```

For a lone `igraph` / `network` / adjacency matrix, `as_lame_y(g)` returns a plain `Y`.

**From ERGM:** `lame` is not an `ergm` substitute -- ERGM specifies a joint distribution over the whole graph via change statistics (`gwesp`, `triangle`), while AME assumes dyads are *conditionally independent* given the actor effects, latent positions, and covariates, and captures higher-order structure through those random effects. So AME sidesteps ERGM's degeneracy failure modes but does not give change-statistic interpretations, and it tends to under-predict transitivity. If `gwesp` is your target, stay in ERGM; if you want actor heterogeneity, dynamics, or forecasts, AME is the better fit. `vignette("cross_sec_ame")` scores ERGM-style `triangle` / two-path / clustering statistics through `gof(fit, custom_gof = ...)`. The ERGM-style helpers `nodematch()`, `absdiff()`, and `nodefactor()` build dyadic covariates that drop straight into `Xdyad`.

## Visualization

```r
ab_plot(fit, effect = "sender")                 # additive effects; plot_type = "trajectory" over time
uv_plot(fit)                                     # latent positions;  plot_type = "trajectory" over time
trace_plot(fit); gof_plot(fit)                   # mcmc mixing; goodness of fit
```

## Features

- **Modes:** cross-sectional (`ame()`), longitudinal (`lame()`), unipartite and bipartite, six families.
- **Dynamics:** time-varying latent positions (`dynamic_uv`), sender/receiver effects (`dynamic_ab`), and regression coefficients (`dynamic_beta`, per-block AR(1) / RW1 / RW2 / Matern 3/2), with forecasting via `predict(fit, h = K)`.
- **S3 methods** on both MCMC and ALS fits: `coef()`, `vcov()`, `confint()`, `predict()`, `fitted()`, `residuals()`, `summary()`, `update()`, `simulate()`, plus `autoplot()`.
- **Ecosystem:** `broom::tidy()` / `glance()`, `modelsummary`, `as_draws()` for `posterior` / `tidybayes`, and `loo::loo()` when `save_log_lik = TRUE`.
- **Fit assessment:** `simulate()` and `gof()` posterior-predictive checks; `gof_temporal()` for panels.
- **Scale and robustness:** `ame_parallel()` / `lame_parallel()` multi-chain runs with R-hat / ESS via `compute_mcmc_diagnostics()`; `lame_multi()` for several panels sharing coefficients; checkpoint/resume via `checkpoint_path=` and `lame_resume()`.

## Fast (MCMC-free) estimation

`ame_als()` / `lame_als()` fit the same model by iterative block coordinate descent (alternating least squares / IRLS) -- a port of the Social Influence Regression estimator of Minhas & Hoff (2025), implemented in the [`sir`](https://github.com/netify-dev/sir) package. They return point estimates rather than a posterior, so they are much faster.

- **Scope:** `normal`, `binary`, `poisson`; use the MCMC path for `ordinal`, `cbin`, `frn`.
- **Dynamic ALS:** `lame(..., method = "als")` fits dynamic additive effects, selected dynamic coefficients, AR(1) and Student-t latent factors (directed, symmetric, bipartite), and bipartite dynamic `G`, as penalized point estimates. Tune with `als_max_iter`, `als_tol`, `als_stability`.
- **Snap-shift:** `lame_snap_als()` (or `method = "als", dynamic_uv_kind = "snap"`) estimates dynamic positions for normal unipartite/bipartite panels. `snap_prob` is an ALS score, not a posterior probability -- read it through rankings and the stability diagnostics.
- **Uncertainty:** `ame_als_bootstrap()`, or `ame_als(..., bootstrap = N)`, for block/parametric bootstrap intervals; `confint()` uses them when present, the sandwich `vcov()` otherwise.

## Documentation

| Vignette | Topic |
| --- | --- |
| `vignette("lame")` | 5-minute getting-started on a small simulated panel |
| `vignette("lame-overview")` | Full tour on the bundled Dutch college friendship data |
| `vignette("cross_sec_ame")` | Cross-sectional AME on the Add Health friendship network |
| `vignette("bipartite")` | Two-mode networks, cross-sectional and longitudinal |
| `vignette("dynamic_effects")` | `dynamic_uv` / `dynamic_ab` / `dynamic_beta` and the transition-kind decision tree |
| `vignette("forecasting")` | `predict(fit, h = K)` forecasts and counterfactuals |
| `vignette("fast_estimation")` | the `ame_als()` / `lame_als()` point estimator |

## Reproducibility

MCMC results depend on the RNG seed and package versions. Set a seed, and record the environment:

```r
set.seed(6886)
fit <- lame(...)
sessionInfo()          # or pin the toolchain with an renv.lock / a @<commit-sha> install
```

## Citation

```bibtex
@Manual{lame2026,
  title  = {lame: Longitudinal Additive and Multiplicative Effects Models for Networks},
  author = {Cassy Dorff and Shahryar Minhas and Tosin Salau},
  year   = {2026},
  note   = {R package version 1.3.5},
  url    = {https://github.com/netify-dev/lame},
}
```

The dynamic effects draw on Sewell & Chen (2015, [doi:10.1080/01621459.2014.988214](https://doi.org/10.1080/01621459.2014.988214)) and Durante & Dunson (2014, [doi:10.1093/biomet/asu040](https://doi.org/10.1093/biomet/asu040)); the fast estimator adapts Minhas & Hoff (2025, *Political Analysis*).

## Contributors

Cassy Dorff (Vanderbilt), Shahryar Minhas (Michigan State), Tosin Salau (Michigan State).

## License

MIT. See the `LICENSE` file.

## Support

[Issues](https://github.com/netify-dev/lame/issues) · [Discussions](https://github.com/netify-dev/lame/discussions) · [minhassh@msu.edu](mailto:minhassh@msu.edu)
