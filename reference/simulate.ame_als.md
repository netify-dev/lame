# Simulate networks from a fitted ame_als model

Posterior-predictive equivalent for the fast ALS estimator: draws `nsim`
replicates of `Y` from the parametric model implied by the fitted `mu`,
`beta`, `a`, `b`, `U`, `V`, and the family-appropriate noise
distribution. Returned object has class `"ame.sim"` so it is compatible
with `plot_ppc_*`-style consumers (where applicable).

## Usage

``` r
# S3 method for class 'ame_als'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  an `ame_als` fit.

- nsim:

  integer; number of replicates to simulate (default 1).

- seed:

  optional RNG seed.

- ...:

  ignored.

## Value

An object of class `"ame.sim"` with element `Y` – a list of `nsim`
simulated outcome arrays (one matrix for a cross-section, one
list-of-matrices per slice for a longitudinal fit).

## Details

Caveat: ALS is a point estimator. `simulate.ame_als` therefore holds
`mu, beta, a, b, U, V` FIXED at the point estimate and only the noise is
resampled. For uncertainty over the parameters themselves use
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
(whose replicates each carry their own resampled `Y`) and combine those.
