# Log-likelihood is not directly exposed for ame() / lame() fits

[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) and
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) produce
a posterior sample, not a maximum- likelihood point. A pointwise
log-likelihood is computable from the posterior draws but is not stored
on the fit object, so [`logLik()`](https://rdrr.io/r/stats/logLik.html)
(and the AIC / BIC generics that dispatch through it) error out
informatively rather than return a misleading number.

## Usage

``` r
# S3 method for class 'ame'
logLik(object, ...)

# S3 method for class 'lame'
logLik(object, ...)
```

## Arguments

- object:

  an `ame` or `lame` fit.

- ...:

  ignored.

## Value

Never returns; raises an error.

## Details

For Bayesian model comparison use posterior-predictive checks via
[`gof`](https://netify-dev.github.io/lame/reference/gof.md) /
[`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md),
or compute WAIC / LOO yourself from the per-draw log-likelihoods (e.g.
via the loo package on the BETA / VC chains).
