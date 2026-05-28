# Long-format draws of the linear predictor for marginaleffects-style use

Returns a long-format data frame with one row per `(draw, i, j, period)`
combination, giving the per-draw linear predictor (or response-scale
prediction) at each dyad and period. Intended for marginaleffects- /
tidybayes-style downstream summarisation: column names follow the
`.draw` / `.chain` / `.iteration` / `.value` convention so that
[`tidybayes::spread_draws()`](https://mjskay.github.io/tidybayes/reference/spread_draws.html)
and `marginaleffects::posterior_draws()` auto-dispatch on the returned
data frame. Actor and period names from `fit$Y`'s dimnames are carried
forward into the `actor_i`, `actor_j`, `period_label` columns.

## Usage

``` r
prediction_draws_long(
  object,
  newdata = NULL,
  type = c("link", "response"),
  n_draws = 100L,
  seed = NULL
)
```

## Arguments

- object:

  A fitted `lame` object.

- newdata:

  Optional list of `T` dyadic covariate arrays for counterfactual
  evaluation. Default `NULL` uses in-sample covariates.

- type:

  One of `"link"` (default) or `"response"`.

- n_draws:

  Number of posterior draws to use. Default 100.

- seed:

  Optional RNG seed.

## Value

A long-format data frame with columns `.chain`, `.iteration`, `.draw`,
`period`, `period_label`, `i`, `j`, `actor_i`, `actor_j`, `.value`.
Returned as a tibble when the tibble package is available; otherwise a
plain data frame.
