# Package index

## Model Fitting (MCMC)

Core Bayesian MCMC estimators for AME and LAME models

- [`ame()`](https://netify-dev.github.io/lame/reference/ame.md) : AME
  model fitting routine

- [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) : AME
  model fitting routine for longitudinal relational data

- [`lame-package`](https://netify-dev.github.io/lame/reference/lame-package.md)
  : Longitudinal Additive and Multiplicative Effects Models for Networks

- [`ame_parallel()`](https://netify-dev.github.io/lame/reference/ame_parallel.md)
  : Run AME model with multiple parallel chains

- [`lame_parallel()`](https://netify-dev.github.io/lame/reference/lame_parallel.md)
  : Run LAME (longitudinal AME) with multiple parallel chains

- [`lame_multi()`](https://netify-dev.github.io/lame/reference/lame_multi.md)
  : Multi-panel lame() with shared coefficients

- [`lame_resume()`](https://netify-dev.github.io/lame/reference/lame_resume.md)
  :

  Resume a
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) MCMC
  run from a checkpoint

- [`ame_options()`](https://netify-dev.github.io/lame/reference/ame_options.md)
  : AME model fitting options

## Forecasting and Temporal Diagnostics

Forecasts, counterfactuals, and diagnostics for dynamic fits

- [`gof_temporal()`](https://netify-dev.github.io/lame/reference/gof_temporal.md)
  : Posterior-predictive temporal-trend test
- [`detect_change_point()`](https://netify-dev.github.io/lame/reference/detect_change_point.md)
  : Detect potential change points in a dynamic_beta posterior path
- [`lfo()`](https://netify-dev.github.io/lame/reference/lfo.md) : Exact
  rolling-origin leave-future-out cross-validation
- [`forecast_pit()`](https://netify-dev.github.io/lame/reference/forecast_pit.md)
  : Probability-integral-transform calibration check for h-step
  forecasts
- [`rhat_dynamic_beta()`](https://netify-dev.github.io/lame/reference/rhat_dynamic_beta.md)
  : Multivariate split-R-hat for dynamic_beta coefficient paths
- [`dynamic_beta_prior_summary()`](https://netify-dev.github.io/lame/reference/dynamic_beta_prior_summary.md)
  : Summarise the implied prior on a time-varying coefficient path
- [`prediction_draws_long()`](https://netify-dev.github.io/lame/reference/prediction_draws_long.md)
  : Long-format draws of the linear predictor for marginaleffects-style
  use
- [`per_actor_slopes()`](https://netify-dev.github.io/lame/reference/per_actor_slopes.md)
  : Post-MCMC per-actor time-varying slopes
- [`als_dynamic_beta()`](https://netify-dev.github.io/lame/reference/als_dynamic_beta.md)
  : Penalised ALS time-varying coefficient estimate

## Model Comparison and Posterior Diagnostics

Information criteria, posterior draws, and held-out evaluation

- [`loo()`](https://netify-dev.github.io/lame/reference/loo.md)
  [`waic()`](https://netify-dev.github.io/lame/reference/loo.md) :
  Generic dispatcher for loo / waic on ame / lame fits
- [`loo(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/loo.ame.md)
  [`loo(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/loo.ame.md)
  [`loo(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/loo.ame.md)
  : Approximate leave-one-out cross-validation for AME / LAME fits
- [`waic(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/waic.ame.md)
  [`waic(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/waic.ame.md)
  [`waic(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/waic.ame.md)
  : WAIC for AME / LAME fits
- [`as_draws()`](https://netify-dev.github.io/lame/reference/as_draws.md)
  : Generic dispatcher for posterior::as_draws on lame fits
- [`as_draws(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/as_draws.ame.md)
  [`as_draws(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/as_draws.ame.md)
  [`as_draws(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/as_draws.ame.md)
  : Convert an AME / LAME fit to a posterior draws object
- [`prior_summary()`](https://netify-dev.github.io/lame/reference/prior_summary.md)
  : Print the priors used by an AME / LAME / ame_als fit
- [`evaluate_heldout()`](https://netify-dev.github.io/lame/reference/evaluate_heldout.md)
  : Held-out predictive evaluation for an ame / lame fit
- [`read_log_lik()`](https://netify-dev.github.io/lame/reference/read_log_lik.md)
  : Read the per-iteration log-lik matrix back from on-disk chunks
- [`snap_index_draws()`](https://netify-dev.github.io/lame/reference/snap_index_draws.md)
  : Extract posterior draws of snap indices
- [`snap_index_summary()`](https://netify-dev.github.io/lame/reference/snap_index_summary.md)
  : Summarize posterior snap indices
- [`snap_category_summary()`](https://netify-dev.github.io/lame/reference/snap_category_summary.md)
  : Summarize snap indices by actor category
- [`snap_rank_summary()`](https://netify-dev.github.io/lame/reference/snap_rank_summary.md)
  : Summarize posterior rank uncertainty for snap years

## Tidyverse Integration

broom and ggplot2 methods for fitted models

- [`tidy()`](https://netify-dev.github.io/lame/reference/tidy.md) :

  S3 generic for `tidy`

- [`tidy(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/tidy.ame.md)
  [`tidy(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/tidy.ame.md)
  :

  Tidy method for fitted `ame` / `lame` objects

- [`tidy(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/tidy.ame_als.md)
  [`tidy(`*`<lame_als>`*`)`](https://netify-dev.github.io/lame/reference/tidy.ame_als.md)
  :

  Tidy method for fitted `ame_als` / `lame_als` objects

- [`tidy(`*`<boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/tidy.boot_ame.md)
  :

  Tidy method for a standalone bootstrap object (`boot_ame`)

- [`glance()`](https://netify-dev.github.io/lame/reference/glance.md) :

  S3 generic for `glance`

- [`glance(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/glance.ame.md)
  [`glance(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/glance.ame.md)
  :

  Glance method for fitted `ame` / `lame` objects

- [`glance(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/glance.ame_als.md)
  [`glance(`*`<lame_als>`*`)`](https://netify-dev.github.io/lame/reference/glance.ame_als.md)
  :

  Glance method for fitted `ame_als` / `lame_als` objects

- [`autoplot(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/autoplot.lame.md)
  [`autoplot(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/autoplot.lame.md)
  : Ribbon plot of time-varying coefficients (or coefplot for static
  fits)

- [`autoplot(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/autoplot.ame_als.md)
  [`autoplot(`*`<lame_als>`*`)`](https://netify-dev.github.io/lame/reference/autoplot.ame_als.md)
  : autoplot method for ALS fits

## Fast Estimator (MCMC-free)

Iterative block coordinate descent point estimator with bootstrap
uncertainty

- [`ame_als()`](https://netify-dev.github.io/lame/reference/ame_als.md)
  : Fast (MCMC-free) AME estimation for a cross-sectional network
- [`lame_als()`](https://netify-dev.github.io/lame/reference/lame_als.md)
  : Fast (MCMC-free) AME estimation for a longitudinal network
- [`lame_snap_als()`](https://netify-dev.github.io/lame/reference/lame_snap_als.md)
  : Fast approximate dynamic snap-shift AME estimator
- [`ame_als_bootstrap()`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
  [`boot_ame()`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
  : Bootstrap uncertainty for the fast AME estimator
- [`ame_als_refit()`](https://netify-dev.github.io/lame/reference/ame_als_refit.md)
  : Refit a fast AME model with a warm start
- [`sampler_describe()`](https://netify-dev.github.io/lame/reference/sampler_describe.md)
  : Describe the estimator behind a fitted object

## S3 Methods

Standard R methods for model objects

- [`coef(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/coef.ame.md)
  [`coef(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/coef.ame.md)
  : Extract model coefficients from AME model

- [`vcov(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/vcov.ame.md)
  [`vcov(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/vcov.ame.md)
  : Posterior covariance of AME model coefficients

- [`confint(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/confint.ame.md)
  [`confint(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/confint.ame.md)
  : Bayesian credible intervals for AME model parameters

- [`predict(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/predict.ame.md)
  : Predict method for AME models

- [`predict(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/predict.lame.md)
  : Predict method for LAME models

- [`fitted(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/fitted.ame.md)
  : Extract fitted values from AME model

- [`fitted(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/fitted.lame.md)
  : Extract fitted values from LAME model

- [`residuals(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/residuals.ame.md)
  : Extract residuals from AME model

- [`residuals(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/residuals.lame.md)
  : Extract residuals from LAME model

- [`simulate(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/simulate.ame.md)
  : Simulate networks from a fitted AME model

- [`simulate(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/simulate.lame.md)
  : Simulate longitudinal networks from a fitted LAME model

- [`summary(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/summary.ame.md)
  : Summary of an AME object

- [`summary(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/summary.lame.md)
  : Summary of a LAME object

- [`print(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/print.ame.md)
  : Print method for AME model objects

- [`print(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/print.lame.md)
  : Print method for LAME objects

- [`print(`*`<summary.ame>`*`)`](https://netify-dev.github.io/lame/reference/print.summary.ame.md)
  : Print method for summary.ame objects

- [`print(`*`<summary.lame>`*`)`](https://netify-dev.github.io/lame/reference/print.summary.lame.md)
  : Print method for summary.lame objects

- [`plot(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/plot.ame.md)
  : Simple diagnostic plot for AME model fit

- [`plot(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/plot.lame.md)
  : Plot diagnostics for a LAME model fit

- [`coef(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/coef.ame_als.md)
  : Extract coefficients from a fast AME fit

- [`vcov(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/vcov.ame_als.md)
  : Sandwich covariance for the regression coefficients of a fast AME
  fit

- [`confint(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/confint.ame_als.md)
  : Confidence intervals for a fast AME fit

- [`fitted(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/fitted.ame_als.md)
  : Extract fitted values from a fast AME fit

- [`residuals(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/residuals.ame_als.md)
  : Residuals from a fast AME fit

- [`predict(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/predict.ame_als.md)
  : Predictions from a fast AME fit

- [`logLik(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/logLik.ame_als.md)
  : Log-likelihood is not defined for a fast AME fit

- [`plot(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/plot.ame_als.md)
  : Plot the convergence of a fast AME fit

- [`print(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/print.ame_als.md)
  : Print an ame_als object

- [`summary(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/summary.ame_als.md)
  : Summarize an ame_als object

- [`print(`*`<boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/print.boot_ame.md)
  : Print bootstrap results for a fast AME fit

- [`summary(`*`<boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/summary.boot_ame.md)
  [`print(`*`<summary.boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/summary.boot_ame.md)
  : Summarize bootstrap results for a fast AME fit

- [`confint(`*`<boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/confint.boot_ame.md)
  : Confidence intervals from a fast AME bootstrap

- [`vcov(`*`<boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/vcov.boot_ame.md)
  : Bootstrap covariance of the regression coefficients

- [`coef(`*`<boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/coef.boot_ame.md)
  : Point estimates from a fast AME bootstrap

- [`fitted(`*`<boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/boot_ame-no-fitted.md)
  [`residuals(`*`<boot_ame>`*`)`](https://netify-dev.github.io/lame/reference/boot_ame-no-fitted.md)
  : fitted/residuals are not defined for a bootstrap object

- [`formula(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/formula.ame.md)
  [`formula(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/formula.ame.md)
  : formula() is not defined for an ame() / lame() fit

- [`logLik(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/logLik.ame.md)
  [`logLik(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/logLik.ame.md)
  : Log-likelihood is not directly exposed for ame() / lame() fits

- [`nobs(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/nobs.ame.md)
  [`nobs(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/nobs.ame.md)
  : Number of observed dyads in an AME / LAME fit

- [`nobs(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/nobs.ame_als.md)
  : Number of observed dyads in an ame_als fit

- [`update(`*`<ame>`*`)`](https://netify-dev.github.io/lame/reference/update.ame.md)
  [`update(`*`<lame>`*`)`](https://netify-dev.github.io/lame/reference/update.ame.md)
  : Update an AME / LAME fit

- [`update(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/update.ame_als.md)
  [`update(`*`<lame_als>`*`)`](https://netify-dev.github.io/lame/reference/update.ame_als.md)
  :

  Update an `ame_als` / `lame_als` fit

- [`simulate(`*`<ame_als>`*`)`](https://netify-dev.github.io/lame/reference/simulate.ame_als.md)
  : Simulate networks from a fitted ame_als model

- [`ab_plot.ame_als()`](https://netify-dev.github.io/lame/reference/ab_plot.ame_als.md)
  : Additive-effects plot for an ame_als fit

- [`gof_plot.ame_als()`](https://netify-dev.github.io/lame/reference/gof_plot.ame_als.md)
  : Goodness-of-fit check for an ame_als fit

- [`coef(`*`<als_dynamic_beta>`*`)`](https://netify-dev.github.io/lame/reference/coef.als_dynamic_beta.md)
  : Extract beta path from a penalised-ALS object

- [`print(`*`<als_dynamic_beta>`*`)`](https://netify-dev.github.io/lame/reference/print.als_dynamic_beta.md)
  : Print method for penalised ALS time-varying beta

- [`print(`*`<gof_temporal>`*`)`](https://netify-dev.github.io/lame/reference/print.gof_temporal.md)
  : Print method for gof_temporal output

- [`print(`*`<lame_multi>`*`)`](https://netify-dev.github.io/lame/reference/print.lame_multi.md)
  : Print method for lame_multi

- [`print(`*`<lfo_lame>`*`)`](https://netify-dev.github.io/lame/reference/print.lfo_lame.md)
  : Print method for lfo() results

- [`print(`*`<per_actor_slopes>`*`)`](https://netify-dev.github.io/lame/reference/print.per_actor_slopes.md)
  : Print method for per_actor_slopes

- [`print(`*`<summary.ame_als>`*`)`](https://netify-dev.github.io/lame/reference/print.summary.ame_als.md)
  : Print a fast AME summary

## Visualization

Plotting functions for diagnostics and model exploration

- [`trace_plot()`](https://netify-dev.github.io/lame/reference/trace_plot.md)
  : MCMC trace plots and density plots for AME/LAME model parameters
- [`gof_plot()`](https://netify-dev.github.io/lame/reference/gof_plot.md)
  : Visualize goodness-of-fit statistics for AME and LAME models
- [`ab_plot()`](https://netify-dev.github.io/lame/reference/ab_plot.md)
  : Visualize sender and receiver random effects
- [`uv_plot()`](https://netify-dev.github.io/lame/reference/uv_plot.md)
  : Visualize multiplicative effects (latent factors) from AME models

## Goodness of Fit

Posterior predictive checks and GOF statistics

- [`gof()`](https://netify-dev.github.io/lame/reference/gof.md) :
  Compute GOF statistics from saved posterior samples
- [`gof_stats()`](https://netify-dev.github.io/lame/reference/gof_stats.md)
  : Goodness of fit statistics
- [`gof_stats_unipartite()`](https://netify-dev.github.io/lame/reference/gof_stats_unipartite.md)
  : Goodness of fit statistics for unipartite networks
- [`gof_stats_bipartite()`](https://netify-dev.github.io/lame/reference/gof_stats_bipartite.md)
  : Goodness of fit statistics for bipartite networks
- [`compute_gof_bipartite()`](https://netify-dev.github.io/lame/reference/compute_gof_bipartite.md)
  : Compute GOF statistics for bipartite networks

## Latent Space

Extract and align latent positions

- [`latent_positions()`](https://netify-dev.github.io/lame/reference/latent_positions.md)
  : Extract latent positions as a tidy data frame
- [`procrustes_align()`](https://netify-dev.github.io/lame/reference/procrustes_align.md)
  : Procrustes alignment of latent positions across time

## Posterior Utilities

Tools for working with posterior samples

- [`simulate_posterior()`](https://netify-dev.github.io/lame/reference/simulate_posterior.md)
  : Simulate posterior distributions from fitted AME model
- [`posterior_quantiles()`](https://netify-dev.github.io/lame/reference/posterior_quantiles.md)
  : Extract posterior quantiles for model components
- [`posterior_options()`](https://netify-dev.github.io/lame/reference/posterior_options.md)
  : Options for saving posterior samples during MCMC
- [`reconstruct_EZ()`](https://netify-dev.github.io/lame/reference/reconstruct_EZ.md)
  [`reconstruct_UVPM()`](https://netify-dev.github.io/lame/reference/reconstruct_EZ.md)
  : Reconstruct EZ and UVPM matrices from AME model output
- [`combine_ame_chains()`](https://netify-dev.github.io/lame/reference/combine_ame_chains.md)
  : Combine multiple AME chains
- [`compute_mcmc_diagnostics()`](https://netify-dev.github.io/lame/reference/compute_mcmc_diagnostics.md)
  : Compute MCMC convergence diagnostics for multiple chains

## Simulation

Simulate network data from fitted models

- [`print(`*`<ame.sim>`*`)`](https://netify-dev.github.io/lame/reference/print.ame.sim.md)
  [`print(`*`<lame.sim>`*`)`](https://netify-dev.github.io/lame/reference/print.ame.sim.md)
  : Print methods for AME and LAME simulation objects
- [`summary(`*`<ame.sim>`*`)`](https://netify-dev.github.io/lame/reference/summary.ame.sim.md)
  : Summary method for AME simulations
- [`summary(`*`<lame.sim>`*`)`](https://netify-dev.github.io/lame/reference/summary.lame.sim.md)
  : Summary method for LAME simulations
- [`simY_nrm()`](https://netify-dev.github.io/lame/reference/simY_nrm.md)
  : Simulate a normal relational matrix
- [`simY_bin()`](https://netify-dev.github.io/lame/reference/simY_bin.md)
  : Simulate a network, i.e. a binary relational matrix
- [`simY_ord()`](https://netify-dev.github.io/lame/reference/simY_ord.md)
  : Simulate an ordinal relational matrix
- [`simY_pois()`](https://netify-dev.github.io/lame/reference/simY_pois.md)
  : Simulate a Poisson relational matrix
- [`simY_tob()`](https://netify-dev.github.io/lame/reference/simY_tob.md)
  : Simulate a tobit relational matrix
- [`simY_frn()`](https://netify-dev.github.io/lame/reference/simY_frn.md)
  : Simulate an relational matrix based on a fixed rank nomination
  scheme
- [`simY_rrl()`](https://netify-dev.github.io/lame/reference/simY_rrl.md)
  : Simulate an relational matrix based on a relative rank nomination
  scheme
- [`simZ()`](https://netify-dev.github.io/lame/reference/simZ.md) :
  Simulate Z given its expectation and covariance

## Data Manipulation

Convert between network data formats and build covariates

- [`el2sm()`](https://netify-dev.github.io/lame/reference/el2sm.md) :
  Edgelist to sociomatrix
- [`sm2el()`](https://netify-dev.github.io/lame/reference/sm2el.md) :
  Sociomatrix to edgelist
- [`array_to_list()`](https://netify-dev.github.io/lame/reference/array_to_list.md)
  : Convert array to list.
- [`list_to_array()`](https://netify-dev.github.io/lame/reference/list_to_array.md)
  : Convert list to array
- [`list_to_array_bipartite()`](https://netify-dev.github.io/lame/reference/list_to_array_bipartite.md)
  : Convert bipartite list data to array format
- [`as_lame_y()`](https://netify-dev.github.io/lame/reference/as_lame_y.md)
  : Convert a graph object to a lame-ready adjacency matrix
- [`nodematch()`](https://netify-dev.github.io/lame/reference/nodematch.md)
  [`absdiff()`](https://netify-dev.github.io/lame/reference/nodematch.md)
  [`nodefactor()`](https://netify-dev.github.io/lame/reference/nodematch.md)
  : ERGM-style covariate helpers for ame() / lame()
- [`check_format()`](https://netify-dev.github.io/lame/reference/check_format.md)
  : Validate input data format for lame function
- [`zscores()`](https://netify-dev.github.io/lame/reference/zscores.md)
  : rank-based z-scores

## Memory and Performance

Manage memory usage in large models

- [`compact_ame()`](https://netify-dev.github.io/lame/reference/compact_ame.md)
  : Optimize AME model output for memory efficiency
- [`ame_memory_usage()`](https://netify-dev.github.io/lame/reference/ame_memory_usage.md)
  : Calculate memory usage of AME model components
- [`ame_memory_settings()`](https://netify-dev.github.io/lame/reference/ame_memory_settings.md)
  : Display memory usage information for AME models

## MCMC Samplers

Gibbs sampling functions for latent variables, regression coefficients,
and variance components

- [`rZ_nrm_fc()`](https://netify-dev.github.io/lame/reference/rZ_nrm_fc.md)
  : Simulate missing values in a normal AME model
- [`rZ_bin_fc()`](https://netify-dev.github.io/lame/reference/rZ_bin_fc.md)
  : Simulate Z based on a probit model
- [`rZ_ord_fc()`](https://netify-dev.github.io/lame/reference/rZ_ord_fc.md)
  : Simulate Z given the partial ranks
- [`rZ_pois_fc()`](https://netify-dev.github.io/lame/reference/rZ_pois_fc.md)
  : Gibbs update for latent variable in a Poisson AME model
- [`rZ_tob_fc()`](https://netify-dev.github.io/lame/reference/rZ_tob_fc.md)
  : Simulate Z based on a tobit model
- [`rZ_cbin_fc()`](https://netify-dev.github.io/lame/reference/rZ_cbin_fc.md)
  : Simulate Z given fixed rank nomination data
- [`rZ_frn_fc()`](https://netify-dev.github.io/lame/reference/rZ_frn_fc.md)
  : Simulate Z given fixed rank nomination data
- [`rZ_rrl_fc()`](https://netify-dev.github.io/lame/reference/rZ_rrl_fc.md)
  : Simulate Z given relative rank nomination data
- [`rbeta_ab_fc()`](https://netify-dev.github.io/lame/reference/rbeta_ab_fc.md)
  : Conditional simulation of additive effects and regression
  coefficients
- [`rbeta_ab_rep_fc()`](https://netify-dev.github.io/lame/reference/rbeta_ab_rep_fc.md)
  : Gibbs sampling of additive row and column effects and regression
  coefficient with independent replicate relational data
- [`rSab_fc()`](https://netify-dev.github.io/lame/reference/rSab_fc.md)
  : Gibbs update for additive effects covariance
- [`rSuv_fc()`](https://netify-dev.github.io/lame/reference/rSuv_fc.md)
  : Gibbs update for multiplicative effects covariance
- [`rUV_fc()`](https://netify-dev.github.io/lame/reference/rUV_fc.md) :
  Gibbs sampling of U and V
- [`rUV_rep_fc()`](https://netify-dev.github.io/lame/reference/rUV_rep_fc.md)
  : Gibbs sampling of U and V
- [`rUV_sym_fc()`](https://netify-dev.github.io/lame/reference/rUV_sym_fc.md)
  : Gibbs sampling of U and V
- [`rUV_dynamic_fc()`](https://netify-dev.github.io/lame/reference/rUV_dynamic_fc.md)
  : Gibbs sampling of dynamic U and V with AR(1) evolution
- [`rUV_dynamic_fc_cpp()`](https://netify-dev.github.io/lame/reference/rUV_dynamic_fc_cpp.md)
  : Update dynamic latent positions using AR(1) process
- [`rUV_dynamic_snap_fc()`](https://netify-dev.github.io/lame/reference/rUV_dynamic_snap_fc.md)
  : Gibbs sampling of dynamic U and V with snap-shift dynamics
- [`rUV_dynamic_snap_fc_cpp()`](https://netify-dev.github.io/lame/reference/rUV_dynamic_snap_fc_cpp.md)
  : Update dynamic latent positions with snap-shift model selection
- [`rUV_dynamic_t_fc()`](https://netify-dev.github.io/lame/reference/rUV_dynamic_t_fc.md)
  : Gibbs sampling of dynamic U and V with heavy-tailed (Student-t)
  innovations
- [`rUV_dynamic_t_fc_cpp()`](https://netify-dev.github.io/lame/reference/rUV_dynamic_t_fc_cpp.md)
  : Update dynamic latent positions with heavy-tailed (Student-t) AR(1)
  innovations
- [`raSab_bin_fc()`](https://netify-dev.github.io/lame/reference/raSab_bin_fc.md)
  : Simulate a and Sab from full conditional distributions under bin
  likelihood
- [`raSab_cbin_fc()`](https://netify-dev.github.io/lame/reference/raSab_cbin_fc.md)
  : Simulate a and Sab from full conditional distributions under the
  cbin likelihood
- [`raSab_frn_fc()`](https://netify-dev.github.io/lame/reference/raSab_frn_fc.md)
  : Simulate a and Sab from full conditional distributions under frn
  likelihood
- [`rrho_fc()`](https://netify-dev.github.io/lame/reference/rrho_fc.md)
  : Griddy Gibbs update for dyadic correlation
- [`rrho_mh()`](https://netify-dev.github.io/lame/reference/rrho_mh.md)
  : Metropolis update for dyadic correlation
- [`rrho_mh_rep()`](https://netify-dev.github.io/lame/reference/rrho_mh_rep.md)
  : Metropolis update for dyadic correlation with independent replicate
  data
- [`rs2_fc()`](https://netify-dev.github.io/lame/reference/rs2_fc.md) :
  Gibbs update for dyadic variance
- [`rs2_rep_fc()`](https://netify-dev.github.io/lame/reference/rs2_rep_fc.md)
  : Gibbs update for dyadic variance with independent replicate
  relational data
- [`sample_ab_bipartite()`](https://netify-dev.github.io/lame/reference/sample_ab_bipartite.md)
  : Sample additive effects for bipartite networks
- [`sample_dynamic_ab_cpp()`](https://netify-dev.github.io/lame/reference/sample_dynamic_ab_cpp.md)
  : Sample dynamic additive effects with AR(1) evolution
- [`sample_rho_ab_cpp()`](https://netify-dev.github.io/lame/reference/sample_rho_ab_cpp.md)
  : Sample AR(1) parameter for dynamic additive effects
- [`sample_rho_uv()`](https://netify-dev.github.io/lame/reference/sample_rho_uv.md)
  : Sample AR(1) parameter for dynamic latent factors
- [`sample_sigma_ab_cpp()`](https://netify-dev.github.io/lame/reference/sample_sigma_ab_cpp.md)
  : Sample innovation variance for dynamic additive effects
- [`sample_sigma_uv()`](https://netify-dev.github.io/lame/reference/sample_sigma_uv.md)
  : Sample innovation variance for dynamic latent factors
- [`init_dynamic_ab_cpp()`](https://netify-dev.github.io/lame/reference/init_dynamic_ab_cpp.md)
  : Initialize dynamic additive effects with AR(1) structure
- [`init_dynamic_positions()`](https://netify-dev.github.io/lame/reference/init_dynamic_positions.md)
  : Initialize dynamic latent positions with AR(1) structure
- [`update_variances_bipartite()`](https://netify-dev.github.io/lame/reference/update_variances_bipartite.md)
  : Update variance parameters for bipartite
- [`rUV_dynamic_bip_fc_cpp()`](https://netify-dev.github.io/lame/reference/rUV_dynamic_bip_fc_cpp.md)
  : Bipartite dynamic UV Gibbs update
- [`rZ_nrm_batch_cpp()`](https://netify-dev.github.io/lame/reference/rZ_nrm_batch_cpp.md)
  : Batch normal Z sampling across all time periods
- [`rZ_bin_bip_batch_cpp()`](https://netify-dev.github.io/lame/reference/rZ_bin_bip_batch_cpp.md)
  : Batch binary Z sampling across all time periods (bipartite, rho=0)
- [`rZ_tob_bip_batch_cpp()`](https://netify-dev.github.io/lame/reference/rZ_tob_bip_batch_cpp.md)
  : Batch tobit Z sampling across all time periods (bipartite, rho=0)
- [`compute_XtX_Xty_bip_cpp()`](https://netify-dev.github.io/lame/reference/compute_XtX_Xty_bip_cpp.md)
  : Compute X'X and X'y for bipartite covariate regression
- [`Xbeta_bip_cpp()`](https://netify-dev.github.io/lame/reference/Xbeta_bip_cpp.md)
  : Compute Xbeta product for bipartite networks
- [`rbeta_ab_bip_gibbs_cpp()`](https://netify-dev.github.io/lame/reference/rbeta_ab_bip_gibbs_cpp.md)
  : Full bipartite Gibbs update for beta, a, b
- [`sample_beta_dynamic_cpp()`](https://netify-dev.github.io/lame/reference/sample_beta_dynamic_cpp.md)
  : Sample the dynamic-block beta path via FFBS
- [`sample_beta_static_cpp()`](https://netify-dev.github.io/lame/reference/sample_beta_static_cpp.md)
  : Sample the static-block beta conditional on the dynamic path
- [`sample_rho_beta_cpp()`](https://netify-dev.github.io/lame/reference/sample_rho_beta_cpp.md)
  : Sample the AR(1) rho for each dynamic block
- [`sample_sigma_beta_cpp()`](https://netify-dev.github.io/lame/reference/sample_sigma_beta_cpp.md)
  : Sample the AR(1) innovation sigma for each dynamic block
- [`get_EZ_dynamic_beta_cpp()`](https://netify-dev.github.io/lame/reference/get_EZ_dynamic_beta_cpp.md)
  : Compute EZ when beta is time-varying

## Internal Helpers

Design matrix construction and low-level utilities

- [`init_bipartite_startvals()`](https://netify-dev.github.io/lame/reference/bipartite_helpers.md)
  : Bipartite network helper functions
- [`design_array()`](https://netify-dev.github.io/lame/reference/design_array.md)
  : Computes the design socioarray of covariate values
- [`design_array_listwisedel()`](https://netify-dev.github.io/lame/reference/design_array_listwisedel.md)
  : Computes the design socioarray of covariate values
- [`get_design_rep()`](https://netify-dev.github.io/lame/reference/get_design_rep.md)
  : Create design array for replicate data
- [`get_fit_object()`](https://netify-dev.github.io/lame/reference/get_fit_object.md)
  : Get fitted object from MCMC results
- [`get_start_vals()`](https://netify-dev.github.io/lame/reference/get_start_vals.md)
  : Get fitted object from MCMC results
- [`precomputeX()`](https://netify-dev.github.io/lame/reference/precomputeX.md)
  : Precompute design matrix statistics
- [`Xbeta()`](https://netify-dev.github.io/lame/reference/Xbeta.md) :
  Linear combinations of submatrices of an array
- [`ldZgbme()`](https://netify-dev.github.io/lame/reference/ldZgbme.md)
  : log density for GBME models
- [`llsrmRho()`](https://netify-dev.github.io/lame/reference/llsrmRho.md)
  : SRM log likelihood evaluated on a grid of rho-values
- [`mhalf()`](https://netify-dev.github.io/lame/reference/mhalf.md) :
  Symmetric square root of a matrix
- [`rmvnorm()`](https://netify-dev.github.io/lame/reference/rmvnorm.md)
  : Simulation from a multivariate normal distribution
- [`rwish()`](https://netify-dev.github.io/lame/reference/rwish.md) :
  Simulation from a Wishart distribution

## Datasets

Example network datasets

- [`YX_nrm`](https://netify-dev.github.io/lame/reference/YX_nrm.md) :
  normal relational data and covariates
- [`YX_bin`](https://netify-dev.github.io/lame/reference/YX_bin.md) :
  binary relational data and covariates
- [`YX_bin_list`](https://netify-dev.github.io/lame/reference/YX_bin_list.md)
  : Synthetic longitudinal binary relational data, list-form
  (latent-scale)
- [`YX_bin_long`](https://netify-dev.github.io/lame/reference/YX_bin_long.md)
  : synthetic longitudinal binary relational data (latent-scale)
- [`YX_ord`](https://netify-dev.github.io/lame/reference/YX_ord.md) :
  ordinal relational data and covariates
- [`YX_cbin`](https://netify-dev.github.io/lame/reference/YX_cbin.md) :
  Censored binary nomination data and covariates
- [`YX_frn`](https://netify-dev.github.io/lame/reference/YX_frn.md) :
  Fixed rank nomination data and covariates
- [`YX_rrl`](https://netify-dev.github.io/lame/reference/YX_rrl.md) :
  row-specific ordinal relational data and covariates
- [`IR90s`](https://netify-dev.github.io/lame/reference/IR90s.md) :
  International relations in the 90s
- [`coldwar`](https://netify-dev.github.io/lame/reference/coldwar.md) :
  Cold War data
- [`comtrade`](https://netify-dev.github.io/lame/reference/comtrade.md)
  : Comtrade data
- [`lazegalaw`](https://netify-dev.github.io/lame/reference/lazegalaw.md)
  : Lazega's law firm data
- [`dutchcollege`](https://netify-dev.github.io/lame/reference/dutchcollege.md)
  : Dutch college data
- [`sampsonmonks`](https://netify-dev.github.io/lame/reference/sampsonmonks.md)
  : Sampson's monastery data
- [`sheep`](https://netify-dev.github.io/lame/reference/sheep.md) :
  Sheep dominance data
- [`addhealthc3`](https://netify-dev.github.io/lame/reference/addhealthc3.md)
  : AddHealth community 3 data
- [`addhealthc9`](https://netify-dev.github.io/lame/reference/addhealthc9.md)
  : AddHealth community 9 data
- [`vignette_data`](https://netify-dev.github.io/lame/reference/vignette_data.md)
  : TIES sanctions data for vignettes
- [`Y`](https://netify-dev.github.io/lame/reference/Y.md) : Relational
  matrix
- [`Xdyad`](https://netify-dev.github.io/lame/reference/Xdyad.md) :
  Dyadic covariates
- [`Xrow`](https://netify-dev.github.io/lame/reference/Xrow.md) : Row
  covariates
- [`Xcol`](https://netify-dev.github.io/lame/reference/Xcol.md) : Column
  covariates
