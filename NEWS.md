# lame 1.3.5

* `lfo()` refits no longer carry two intercepts. The stored design
  (`fit$Xlist`) already includes the intercept slice, and `lfo()` passed it
  back through `Xdyad`, so each refit added a second, exactly collinear
  intercept; the stored slice is now stripped before refitting.
* `lame()` no longer errors when called through a wrapper that forwards
  `...` (or with only positional arguments); the `model.name`
  partial-match guard assumed the call had names.
* `uv_plot()` now honors `layout = "circle"` for dynamic snapshot plots; it was
  silently drawing a biplot. Also corrected the latent-space description in the
  help to the inner-product (eigenmodel) interpretation.
* `simulate_posterior()` now returns latent-factor and coefficient draws for
  dynamic fits (`dynamic_uv`, `dynamic_beta`) instead of erroring, and its
  bipartite draws include the interaction matrix `G` (they previously
  reconstructed `U V'` rather than `U G V'`). Credible-interval forecasts keep
  their actor names.
* Fixed the bundled `sampsonmonks` data: the three liking layers (`like_m2`,
  `like_m1`, `like`) had their actors in a different order from the array's
  row/column names (a problem inherited from `amen`), so those layers were
  mislabeled. All ten layers now share the same actor order.
* `ab_plot(effect = "both")` now works for `ame` and `lame` fits, drawing the
  sender and receiver effects as two panels (for dynamic fits, in the
  snapshot views); it previously only worked for `ame_als` fits.
* `ame()` and `lame()` gained a `multistart` argument that is forwarded to
  the static ALS estimators when `method = "als"`; it used to be dropped
  with an unused-argument error (`lame()`) or warning (`ame()`).
* `symmetric = TRUE` with the rank-nomination families (`cbin`, `frn`) is now
  rejected with an explanation; it used to route silently through the
  directed samplers.
* `latent_positions()` and `procrustes_align()` now say plainly that
  alignment is across the periods of one fit, and how to align two fits by
  hand.
* Poisson fits with many zero cells no longer pin the actor effects and latent
  factors at zero and inflate the fitted level (about 3x at 5% density, 1.2-1.4x
  at 20%): their prior scale vanished as zeros dominated. `ame()` also gained a
  level move so the predicted total settles quickly; see the notes in `?ame`.
* Fixed a scale drift in bipartite `dynamic_uv` fits. Only the product
  `U G V'` is identified, and the leftover scale used to wander: `G` shrank
  toward zero while the latent coordinates grew, until they hit
  `prior$uv_max_abs` and the sampler began rejecting 20-40% of its own `U`/`V`
  proposals on realistic panel sizes. The latent coordinates are now kept at
  or below unit scale, with `G` absorbing the scale, which leaves `U G V'`
  untouched. For normal outcomes each margin is held at unit scale, which also
  stops the slow drift of scale into the interaction matrix on long
  `dynamic_G` chains, so the raw coordinates and the persistence parameters
  stay readable. For binary and Poisson outcomes the coordinates are only
  ever shrunk: a saturated likelihood wants the product to grow without
  bound, and pulling the coordinates back up would push all of that growth
  into `G`. Estimates of the identified quantities were not affected, but the
  raw coordinate scale was arbitrary and sampling efficiency suffered. The
  warning for this case now says what to do instead of suggesting a longer
  chain, which made it worse.
* Longitudinal fits estimate the latent-persistence parameters (`rho_uv`,
  `sigma_uv`) more reliably: the unipartite sampler now updates them every
  sweep instead of every tenth, matching the earlier bipartite fix, so they
  no longer sit near their starting value on default-length chains.
* Warm-starting a dynamic MCMC run from an ALS fit now carries the dynamic
  hyperparameters along: `get_start_vals()` (and the bipartite path) pass
  `rho_uv`, `sigma_uv`, `rho_ab`, and `sigma_ab` through instead of dropping
  them.
* ALS fits with `R > 0` expose the fitted multiplicative surface as `UVPM`
  (and `ULUPM` for symmetric fits), matching the MCMC fit objects.
* Dynamic ALS (`lame(method = "als")` with dynamic effects) no longer runs away
  on sparse binary or Poisson panels. It now reweights only once the block
  updates have settled or spent their step budget, keeps a reweighted step
  only when the true penalized deviance falls, and falls back to the best
  state on an overshoot, the same scheme the static route uses; the
  non-normal families also carry a level prior on the dynamic sender and
  receiver paths. A warning still fires if a block update fails to descend,
  and now says so plainly.
* `dynamic_beta` combined with `dynamic_ab` now gives the additive effects a
  real AR(1): the per-period draws pool across neighbouring years and
  `rho_ab` / `sigma_ab` are sampled rather than sitting at their starting
  values.
* Bipartite `dynamic_G = TRUE` fits can no longer be poisoned by a single
  non-finite state draw. The time-varying interaction matrix is accepted only
  when every element is finite, its variance draw keeps the current value when
  the state path is not finite, and the chain falls back to the last finite
  draw otherwise; before, one bad draw made every later update fail and the
  fit returned non-finite values for the rest of the run.
* Tidied up the documentation and test suite.

## Known limitations

* In bipartite `dynamic_G = TRUE` fits with a binary or Poisson outcome, the
  raw latent coordinates can drift toward zero over long chains while `G`
  grows to compensate, and `rho_G` may then sit near its upper bound. The
  fitted product `U G V'` (and so `EZ`, `YPM`, and the forecasts) is
  unaffected; read the raw coordinates and `rho_G` qualitatively in that case.

# lame 1.3.4

* The handful of Gibbs samplers and helpers once adapted from `amen` have been
  rewritten from the underlying published methods, so the package carries no
  derived code and stays under the MIT license.
* Fixed a prior-scaling bug that could silently collapse the latent factors in
  asymmetric `ame()` fits with `R > 0`; fits now match `amen`.
* `dynamic_ab` fits carry per-period posterior SDs, so
  `ab_plot(fit, plot_type = "ribbon")` draws a real credible band.
* A few documented options now do what they say: `log_lik_method = "augmented"`,
  `gof(nsim = NULL)` (uses all draws), and bootstrap error bars in `ab_plot()`
  for `ame_als` fits. Dropped an unused argument and fixed some stale help text.
* Fitting no longer leaves a `.Random.seed` behind or touches `options(warn)`.
* Fixed Rd markup flagged by CRAN.
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
