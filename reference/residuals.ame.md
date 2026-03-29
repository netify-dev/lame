# Extract residuals from AME model

Computes residuals as the difference between observed values and fitted
values. For `type = "response"`, returns `Y - fitted(object)`. For
`type = "pearson"`, returns response residuals scaled by the standard
deviation implied by the family (e.g., `sqrt(p*(1-p))` for binary).

## Usage

``` r
# S3 method for class 'ame'
residuals(object, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  Fitted AME model object (class "ame").

- type:

  Character; `"response"` (default) for raw residuals or `"pearson"` for
  standardized residuals.

- ...:

  Additional arguments (not used).

## Value

An n x n matrix (unipartite) or nA x nB matrix (bipartite) of residuals.
Entries where the original data was `NA` remain `NA`.

## See also

[`fitted.ame`](https://netify-dev.github.io/lame/reference/fitted.ame.md),
[`predict.ame`](https://netify-dev.github.io/lame/reference/predict.ame.md)

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
