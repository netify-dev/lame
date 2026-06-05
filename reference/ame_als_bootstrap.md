# Bootstrap uncertainty for the fast AME estimator

Computes bootstrap standard errors and percentile confidence intervals
for a fast AME fit produced by
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) or
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md).
Two strategies are available:

- `parametric`:

  (default) Simulates fresh outcomes from the fitted model and refits.
  For `normal` and IRLS (`binary`, `poisson`) fits the simulation is on
  the calibrated response scale; for a non-normal `transform` fit –
  whose estimator is a Gaussian fit to a fixed transformed response – it
  is on the Gaussian working scale. The simulated residual variance is
  inflated by the mean-model degrees of freedom so the replicates are
  centred on the point estimate. Available for both cross-sectional and
  longitudinal fits.

- `block`:

  Resamples time slices with replacement. With `block_length = 1` slices
  are resampled independently; `block_length > 1` draws contiguous
  blocks (a moving-block bootstrap), which preserves short-run temporal
  dependence. Longitudinal fits only (requires `T > 1`); with few slices
  its intervals are necessarily coarse.

Every replicate is refit by
[`ame_als_refit`](https://netify-dev.github.io/lame/reference/ame_als_refit.md),
warm-started from the original point estimate so that replicates do not
drift to different local optima. Replicates that error or return
non-finite values are dropped and counted (`n_valid` / `n_total`).
Standard errors are replicate column standard deviations; confidence
intervals use the percentile method.

## Usage

``` r
ame_als_bootstrap(
  object,
  R = 200,
  type = c("parametric", "block"),
  block_length = 1,
  seed = NULL,
  verbose = TRUE
)

boot_ame(
  object,
  R = 200,
  type = c("parametric", "block"),
  block_length = 1,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  an `ame_als` object from
  [`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) or
  [`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md).

- R:

  integer number of bootstrap replicates (default 200).

- type:

  `"parametric"` (default) or `"block"`; see Description.

- block_length:

  block length for the block bootstrap: `1` (default) resamples time
  slices independently; an integer `> 1` draws contiguous blocks of that
  many slices (a moving-block bootstrap), appropriate when the series
  has short-run temporal dependence. Capped at `T`; ignored for the
  parametric bootstrap.

- seed:

  optional integer random seed.

- verbose:

  logical; print progress (default `TRUE`).

## Value

An object of class `"boot_ame"` with components including `coefs`
(replicate intercept + regression coefficients), `se`, `ci_lo`, `ci_hi`,
`point_est`, `param_names`; `vc_*` for the variance components;
`a_coefs`, `b_coefs`, `se_a`, `se_b`; `U_aligned`, `V_aligned`, `U_raw`,
`V_raw`, `se_U`, `se_V` (when `R > 0`); and `n_valid`, `n_total`,
`type`, `family`.

## Details

**Choice of inference.** The Social Influence Regression paper of Hoff &
Minhas (2025) derives its primary standard errors from the observed
Hessian (classical \\-H^{-1}\\ and the sandwich/robust estimator
\\H^{-1} S H^{-1}\\). Two features of the AME model make a Hessian-based
variance awkward here. First, without an explicit gauge fix the rank-`R`
multiplicative term is identified only up to a full \\R\times R\\
rotation/reflection (\\U \to U R\\, \\V \to V R^{-\top}\\), leaving the
joint Hessian rank-deficient. Second, even with a gauge fixed, the
non-normal families are fit on a Gaussian working response, so a Hessian
computed from that working objective is not calibrated to the family
likelihood. The bootstrap side-steps both issues and is the recommended
uncertainty tool here, mirroring `sir::boot_sir()`.

**Alignment-sensitive quantities.** The regression coefficients `beta`,
the additive effects `a`, `b` and the variance components are
rotation-invariant and are aggregated directly. The multiplicative
factors `U`, `V` are *not*: each replicate is Procrustes-aligned to the
original fit before its standard errors are computed. The object stores
both the raw and aligned replicate factors so the effect of alignment
can be inspected. The alignment is well-determined only when the
multiplicative singular values are well separated; with near-equal
singular values the per-column `U`/`V` standard errors reflect an
unstable rotation, and the subspace — or, for symmetric models, the
eigenvalues `L` — is the summary of record. For a symmetric model the
eigenvalues `L` of the multiplicative term are aggregated and reported
as the primary multiplicative-uncertainty summary.

**The parametric bootstrap** is a full parametric bootstrap of the AME
model: each replicate draws fresh additive random effects `a`, `b` from
their fitted dispersion (`va`, `vb`, `cab`) and fresh residuals, holding
`mu`, `beta` and the multiplicative term fixed. Regenerating the
additive effects (rather than holding the fitted `a`, `b` fixed) is what
gives the intercept and node-covariate standard errors their actor-level
sampling-variability component.

**Assumptions.** The block bootstrap treats the time slices as
exchangeable replicates of the static-effects model — appropriate for
that model, but not for strongly trended or serially dependent series
(use the dynamic
[`lame`](https://netify-dev.github.io/lame/reference/lame.md) there);
with only a few time slices it is necessarily coarse, its intervals
correspondingly imprecise and somewhat anti-conservative, so
`parametric` is the default. A variance component fit with `R > 0` can
still carry a small residual parametric-bootstrap bias (the low-rank
refit re-absorbs simulated noise);
[`summary()`](https://rdrr.io/r/base/summary.html) flags any point
estimate that falls outside its interval. A binary IRLS fit can carry a
finite-sample (incidental-parameters) bias in the point estimator
itself, which the bootstrap reproduces rather than removes. The MCMC
[`ame`](https://netify-dev.github.io/lame/reference/ame.md) /
[`lame`](https://netify-dev.github.io/lame/reference/lame.md) path gives
posterior summaries when that is the target.

## References

Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
Influence Regression. *Political Analysis*. The block and parametric
bootstrap follow the inference scheme of `sir::boot_sir()` developed for
the SIR estimator.

## See also

[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) and
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md),
both of which accept `bootstrap = N`, `bootstrap_type`,
`bootstrap_block_length`, `bootstrap_seed` to do this work in a single
call when you are about to fit the model anyway. `ame_als_bootstrap()`
is the post-hoc path (bootstrap an *already-fit* object without
refitting from scratch).
[`ame_als_refit`](https://netify-dev.github.io/lame/reference/ame_als_refit.md)
for warm-start refits;
[`confint.boot_ame`](https://netify-dev.github.io/lame/reference/confint.boot_ame.md)
for the bootstrap CI extractor.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
Y <- replicate(6, { m <- matrix(rnorm(225), 15, 15); diag(m) <- NA; m },
               simplify = FALSE)
fit <- lame_als(Y, R = 1, family = "normal", verbose = FALSE)
# post-hoc bootstrap of an existing fit:
bt <- ame_als_bootstrap(fit, R = 50, type = "block",
                            seed = 1, verbose = FALSE)
# equivalent one-shot call:
# fit_b <- lame_als(Y, R = 1, family = "normal", verbose = FALSE,
#                   bootstrap = 50, bootstrap_type = "block",
#                   bootstrap_seed = 1)
print(bt)
#> 
#> Bootstrap results: fast AME estimator
#> Type: block | Family: normal | Replicates: 50/50 valid
#> 
#>           Estimate Boot SE    2.5%  97.5%
#> intercept   0.0092  0.0285 -0.0352 0.0768
#> 
#> beta, a, b and variance components aggregate directly; U/V are
#> Procrustes-aligned.
#> ! Block bootstrap on a static-effects fit: SEs for the intercept,
#> node-covariate coefficients and additive variance components (va, vb, cab) can
#> severely under-cover because slice resampling does not perturb actor-level
#> structure. Use `type = "parametric"` for those parameters.
# }
```
