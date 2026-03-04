# Extract posterior quantiles for model components

Extract posterior quantiles for model components

## Usage

``` r
posterior_quantiles(
  fit,
  component = c("beta", "UV", "ab"),
  probs = c(0.025, 0.5, 0.975)
)
```

## Arguments

- fit:

  Fitted ame model object

- component:

  Character; which component: "beta", "UV", "ab"

- probs:

  Numeric vector of probabilities for quantiles

## Value

Matrix or array of posterior quantiles

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
