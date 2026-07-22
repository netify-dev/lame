# lame 1.3.3

* The package is now distributed under GPL-3 rather than MIT, matching the
  license of `amen`, from which several samplers and helpers are derived.
  Peter Hoff is credited as contributor and copyright holder; see
  `inst/COPYRIGHTS`.
* Fitting no longer leaves a `.Random.seed` behind in the global environment
  when the session had none, and `lame()` no longer sets `options(warn)`.
* Fixed Rd markup reported by CRAN's incoming checks.
* Fixed a prior-scaling bug that could silently collapse the latent factors
  (and attenuate coefficients) in asymmetric `ame()` fits with `R > 0`;
  fits now match `amen`.
* Time-varying coefficients (`dynamic_beta`) compose freely with the
  multiplicative latent factors (`R > 0`) and additive sender/receiver effects,
  for every family and network type. Coefficients can follow AR(1),
  random-walk, or Matern dynamics, and `predict()` carries that drift into its
  forecasts.
* `netify` objects can be passed straight to `ame()`, `lame()`, and the ALS
  fitters; a network's `symmetric` attribute is honored automatically, and
  named covariates are aligned to `Y` by actor.
* Set `posterior_opts = list(save_UV = TRUE)` to keep the per-draw latent
  factors (and the bipartite mixing matrix), so `latent_positions()`,
  `uv_plot()`, and the goodness-of-fit tools report posterior uncertainty
  directly. Symmetric fits now store `V_samples` (`U L` per draw) alongside
  `U_samples`, so the per-draw latent similarity `U L U'` is reconstructable
  via `simulate_posterior(fit, "UV")`.
* `summary()` on a `lame` fit exposes the coefficient table under
  `$coefficients`, matching `summary.ame` and the `broom`/`lm` idiom.

## Known limitations

* Bipartite `fit$U` / `fit$V` are a single posterior draw (the rotation is
  unidentified); use `fit$YPM` or `reconstruct_UVPM()` for the stable
  multiplicative structure.
* In the bipartite `dynamic_uv` path, latent persistence `rho_uv` is only
  weakly identified -- read it qualitatively. The unipartite path is unaffected.

# lame 1.2.0

* Cleaned up ALS support and documentation.
* Raised the `lame(method = "als")` iteration cap and exposed dynamic ALS
  convergence component traces.
* Stabilized dynamic-UV MCMC by normalizing the raw `U`/`V` coordinate scale,
  centering additive row/column effects, pooling very sparse additive effects
  toward the prior mean, and carrying ALS fits into MCMC with `als_start_vals()`.

# lame 1.1.0

## New features

* Fast, MCMC-free estimation via `ame_als()` and `lame_als()`: an alternating
  least squares / IRLS point estimator for the normal, binary, and Poisson
  families. Supports parametric-bootstrap and sandwich standard errors and the
  full S3 method set (`coef()`, `vcov()`, `confint()`, `predict()`, `tidy()`,
  `glance()`, and the diagnostic plots).

# lame 1.0.0

## Initial CRAN Submission

### Features

* Cross-sectional network analysis via `ame()` with support for 6 distributional
  families: normal, binary, ordinal, Poisson, censored binary, and fixed rank
  nomination.
* Longitudinal network analysis via `lame()` with dynamic additive and
  multiplicative effects modeled as AR(1) processes.
* Support for both unipartite (square) and bipartite (rectangular) network
  structures.
* Covariate support: dyadic (`Xdyad`), row (`Xrow`), and column (`Xcol`)
  covariates with automatic design matrix construction.
* S3 methods: `print()`, `summary()`, `coef()`, `vcov()`, `confint()`,
  `predict()`, `fitted()`, `residuals()`, `simulate()`, `plot()`.
* Visualization functions: `trace_plot()`, `gof_plot()`, `ab_plot()`,
  `uv_plot()` for MCMC diagnostics and model assessment.
* Goodness-of-fit via posterior predictive checks with `gof()`.
* C++ acceleration via Rcpp and RcppArmadillo for core sampling routines.
* Handles changing actor compositions across time periods.
