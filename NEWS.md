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
