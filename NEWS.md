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

* Cross-sectional network analysis via `ame()` with support for 8 distributional
  families: normal, binary, ordinal, tobit, Poisson, censored binary, fixed rank
  nomination, and row-ranked likelihood.
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
