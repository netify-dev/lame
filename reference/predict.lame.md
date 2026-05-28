# Predict method for LAME models

Generate predictions from a fitted longitudinal AME model. Returns a
list of matrices (one per time point) on the requested scale.

## Usage

``` r
# S3 method for class 'lame'
predict(
  object,
  newdata = NULL,
  type = c("response", "link"),
  h = 0L,
  by_draw = FALSE,
  interval = c("none", "credible"),
  probs = c(0.025, 0.975),
  newexposure = NULL,
  ...
)
```

## Arguments

- object:

  Fitted LAME model object.

- newdata:

  Optional list of `T` dyadic covariate arrays (`[n_row, n_col, p]`
  each, with the same actors as the fit) to compute counterfactual
  predictions. When `NULL`, the training-data predictions are returned.

- type:

  Character; `"response"` (default) or `"link"`.

- h:

  Integer \>= 0: forecast horizon. When `h = 0` (default), returns
  in-sample predictions as before. When `h > 0`, propagates the AR(1)
  (or RW1) state-space model forward by `h` periods and returns a list
  of `h` matrices (one per future period). Requires at least one dynamic
  component on the fit. Warns when posterior \\\rho\_\beta\\ is near 1.

- by_draw:

  When `TRUE` and `h > 0`, returns an `n x n x h x n_draws` array of
  per-draw forecasts instead of per-period means.

- interval:

  One of `"none"` (default) or `"credible"`. When `h > 0` and
  `"credible"`, the per-period output is a list of length-3 lists with
  `$lower`, `$median`, `$upper` matrices computed at the `probs`
  quantiles across posterior draws. Ignored for in-sample (`h = 0`)
  predictions.

- probs:

  Length-2 vector of lower / upper quantiles for the credible interval
  when `interval = "credible"`. Default `c(0.025, 0.975)`.

- newexposure:

  Optional length-`h` non-negative numeric vector of future-period
  exposures (Poisson only). When omitted and the fit has
  `period_exposure` stored, defaults to the last observed exposure; when
  both are absent, defaults to 1.

- ...:

  Additional arguments (not used).

## Value

List of prediction matrices (one per time point).
