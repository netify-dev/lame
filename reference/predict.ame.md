# Predict method for AME models

Generate predictions from fitted AME models, including point estimates
and predictive distributions.

## Usage

``` r
# S3 method for class 'ame'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "distribution"),
  n_samples = 100,
  include_uncertainty = TRUE,
  ...
)
```

## Arguments

- object:

  Fitted AME model object

- newdata:

  Optional new data (covariates) for prediction

- type:

  Character; type of prediction:

  - "response": predicted values on response scale (default)

  - "link": predicted values on link scale

  - "distribution": full posterior predictive distribution

- n_samples:

  For type="distribution", number of posterior samples

- include_uncertainty:

  Logical; include parameter uncertainty (default TRUE)

- ...:

  Additional arguments (not used)

## Value

Depending on type:

- "response"/"link": Matrix of predictions

- "distribution": Array of posterior predictive samples

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Fit model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Point predictions
Y_pred <- predict(fit)

# Predictions on link scale
Y_link <- predict(fit, type = "link")
# }
```
