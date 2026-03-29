# Extract model coefficients from AME model

Returns posterior means of regression coefficients from a fitted AME or
LAME model. The coefficients correspond to the columns of the `BETA`
matrix stored in the model object, which records one draw per
post-burn-in MCMC iteration. The posterior mean is computed as
`colMeans(BETA)`.

## Usage

``` r
# S3 method for class 'ame'
coef(object, ...)

# S3 method for class 'lame'
coef(object, ...)
```

## Arguments

- object:

  Fitted AME model (class `"ame"` or `"lame"`).

- ...:

  Additional arguments (ignored).

## Value

Named numeric vector of posterior mean coefficients.

## Details

For binary models, these are on the probit (latent) scale. Use
[`predict.ame`](https://netify-dev.github.io/lame/reference/predict.ame.md)
with `type = "response"` to get predicted probabilities.

## See also

[`vcov.ame`](https://netify-dev.github.io/lame/reference/vcov.ame.md)
for the posterior covariance matrix,
[`confint.ame`](https://netify-dev.github.io/lame/reference/confint.ame.md)
for credible intervals,
[`summary.ame`](https://netify-dev.github.io/lame/reference/summary.ame.md)
for a full summary table
