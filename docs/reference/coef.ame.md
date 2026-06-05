# Extract model coefficients from AME model

Returns posterior means of regression coefficients from a fitted AME or
LAME model.

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

Named numeric vector (static fit) or `p x T` matrix (`dynamic_beta` fit)
of posterior mean coefficients.

## Details

For a static fit (the default, and any model with
`dynamic_beta = FALSE`), coefficients are returned as a named numeric
vector computed as `colMeans(fit$BETA)`.

For a dynamic fit (`lame(..., dynamic_beta = ...)` where some
coefficient is time-varying), `fit$BETA` is a 3-dimensional array
`[n_stored, p, T]` and `coef.lame` returns a `[p, T]` matrix of
per-period posterior means. Rownames are the coefficient names; colnames
are the period labels (from `names(Y)` or `t1, t2, ...`). Static
coefficients in a dynamic fit are constant across the columns.

For binary models, these are on the probit (latent) scale. Use
[`predict.ame`](https://netify-dev.github.io/lame/reference/predict.ame.md)
with `type = "response"` to get predicted probabilities.

**What [`coef()`](https://rdrr.io/r/stats/coef.html) does not return.**
The multiplicative latent positions \\U\\, \\V\\ are not part of the
coefficient vector; they live on `fit$U` and `fit$V` (or as 3-D arrays
`[n, R, T]` when `dynamic_uv` is on). The additive sender / receiver
effects \\a, b\\ are on `fit$APM` and `fit$BPM`. For a tidy frame of
latent positions use
[`latent_positions`](https://netify-dev.github.io/lame/reference/latent_positions.md);
for sender / receiver lollipops use
[`ab_plot`](https://netify-dev.github.io/lame/reference/ab_plot.md).

## See also

[`vcov.ame`](https://netify-dev.github.io/lame/reference/vcov.ame.md)
for the posterior covariance matrix,
[`confint.ame`](https://netify-dev.github.io/lame/reference/confint.ame.md)
for credible intervals,
[`summary.ame`](https://netify-dev.github.io/lame/reference/summary.ame.md)
for a full summary table
