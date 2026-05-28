# Sandwich covariance for the regression coefficients of a fast AME fit

Returns a heteroskedasticity-robust (optionally dyad-clustered) sandwich
covariance matrix for the intercept and *dyadic*-covariate coefficients
of an
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) fit.
This is a fast analytic alternative to the bootstrap for those
coefficients.

## Usage

``` r
# S3 method for class 'ame_als'
vcov(object, cluster = c("dyad", "none"), ...)
```

## Arguments

- object:

  an `ame_als` fit.

- cluster:

  `"dyad"` (default) for a dyad-clustered robust meat, or `"none"` for
  an HC0 (heteroskedasticity-only) meat.

- ...:

  ignored.

## Value

A covariance matrix for `c(intercept, dyadic coefficients)`, with
matching row/column names.

## Details

The estimate is the conditional sandwich \\B^{-} M B^{-}\\ with bread
\\B = D'WD\\ (\\D\\ the observed intercept + dyadic-covariate design,
\\W\\ the fit's observation weights) and meat \\M\\ the
heteroskedasticity-robust (`cluster = "none"`, an HC0 meat) or
dyad-clustered (`cluster = "dyad"`, the default) outer product of the
weighted score contributions \\w\_\ell e\_\ell d\_\ell\\. For a normal
or transform fit the weights are unit, so this reduces to the ordinary
\\D'D\\ sandwich; for an IRLS fit it uses the final IRLS weights,
matching the estimating equation the fit actually solved. Dyad
clustering pools the score across \\(i,j)\\, \\(j,i)\\ and time, so it
reflects dyadic dependence (reciprocity, repeated observation).

It is **conditional**: the additive effects `a`, `b` and the
multiplicative term are held fixed, so it omits their estimation
uncertainty and is anti-conservative. Node-covariate, additive and
multiplicative standard errors are not returned – use
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
for those and for fully-propagated inference.

## See also

[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md)
for bootstrap uncertainty covering all parameters.
