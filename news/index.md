# Changelog

## lame 1.2.0

- Cleaned up ALS support and documentation.

## lame 1.1.0

### New features

- Fast, MCMC-free estimation via
  [`ame_als()`](https://netify-dev.github.io/lame/reference/ame_als.md)
  and
  [`lame_als()`](https://netify-dev.github.io/lame/reference/lame_als.md):
  an alternating least squares / IRLS point estimator for the normal,
  binary, and Poisson families. Supports parametric-bootstrap and
  sandwich standard errors and the full S3 method set
  ([`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`tidy()`](https://netify-dev.github.io/lame/reference/tidy.md),
  [`glance()`](https://netify-dev.github.io/lame/reference/glance.md),
  and the diagnostic plots).

## lame 1.0.0

### Initial CRAN Submission

#### Features

- Cross-sectional network analysis via
  [`ame()`](https://netify-dev.github.io/lame/reference/ame.md) with
  support for 8 distributional families: normal, binary, ordinal, tobit,
  Poisson, censored binary, fixed rank nomination, and row-ranked
  likelihood.
- Longitudinal network analysis via
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) with
  dynamic additive and multiplicative effects modeled as AR(1)
  processes.
- Support for both unipartite (square) and bipartite (rectangular)
  network structures.
- Covariate support: dyadic (`Xdyad`), row (`Xrow`), and column (`Xcol`)
  covariates with automatic design matrix construction.
- S3 methods: [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).
- Visualization functions:
  [`trace_plot()`](https://netify-dev.github.io/lame/reference/trace_plot.md),
  [`gof_plot()`](https://netify-dev.github.io/lame/reference/gof_plot.md),
  [`ab_plot()`](https://netify-dev.github.io/lame/reference/ab_plot.md),
  [`uv_plot()`](https://netify-dev.github.io/lame/reference/uv_plot.md)
  for MCMC diagnostics and model assessment.
- Goodness-of-fit via posterior predictive checks with
  [`gof()`](https://netify-dev.github.io/lame/reference/gof.md).
- C++ acceleration via Rcpp and RcppArmadillo for core sampling
  routines.
- Handles changing actor compositions across time periods.
