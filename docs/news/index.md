# Changelog

## lame 0.3.0

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
