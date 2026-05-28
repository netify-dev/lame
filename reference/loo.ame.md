# Approximate leave-one-out cross-validation for AME / LAME fits

S3 method for [`loo`](https://mc-stan.org/loo/reference/loo.html) that
uses the per-iteration pointwise log-likelihood stored on the fit object
(`fit$log_lik`) when the model was fit with `save_log_lik = TRUE`.
Returns the standard `loo` object with Pareto-k diagnostics.

## Usage

``` r
# S3 method for class 'ame'
loo(x, ...)

# S3 method for class 'lame'
loo(x, ...)

# S3 method for class 'ame_als'
loo(x, ...)
```

## Arguments

- x:

  A fitted `ame` or `lame` object that has `$log_lik`.

- ...:

  Additional arguments forwarded to
  [`loo::loo.matrix`](https://mc-stan.org/loo/reference/loo.html) (e.g.
  `cores`, `r_eff`).

## Value

A `loo` object.

## Details

**What `log_lik` measures.** For `family` in {normal, binary, cbin,
tobit, poisson, ordinal} the stored pointwise log-likelihood is the
exact family-specific Y density on the response scale, so `elpd_loo` is
directly comparable to a
[`loo()`](https://netify-dev.github.io/lame/reference/loo.md) output
from Stan / brms fit to the same family. For the rank likelihoods `frn`
and `rrl` the exact marginal needs GHK Monte Carlo (Halton sequence); on
the longitudinal
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) path you
can opt in with `log_lik_method = "observed_ghk"`, on the
cross-sectional
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) path the
fallback is the augmented-Z normal approximation (with a one-time
warning). Inspect `fit$log_lik_method` on any fit to see which branch
was used.

**Chunked log-lik portability.** When fit with
`save_log_lik = "chunked"`, the on-disk chunk files default to
[`tempdir()`](https://rdrr.io/r/base/tempfile.html), which is cleared at
the end of the R session. If you intend to
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) the fit and reload it
in a fresh session, supply an explicit persistent `log_lik_path` (e.g.
`"./loglik_chunks"`) so the chunks survive the round trip.

## Examples

``` r
# \donttest{
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 0,
           nscan = 200, burn = 50, odens = 5,
           save_log_lik = TRUE, verbose = FALSE)
if (requireNamespace("loo", quietly = TRUE)) {
  loo_res <- loo::loo(fit)
  print(loo_res)
}
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> 
#> Computed from 40 by 9900 log-likelihood matrix.
#> 
#>          Estimate    SE
#> elpd_loo -14947.4  74.9
#> p_loo       704.7  12.9
#> looic     29894.7 149.9
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume independent draws (r_eff=1).
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.38]   (good)     6114  61.8%   12      
#>    (0.38, 1]   (bad)      3700  37.4%   <NA>    
#>     (1, Inf)   (very bad)   86   0.9%   <NA>    
#> See help('pareto-k-diagnostic') for details.
# }
```
