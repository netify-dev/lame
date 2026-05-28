# Build the full canonical design for a newdata prediction

Build the full canonical design for a newdata prediction

## Usage

``` r
.build_full_design(newdata, model_names, fitted_design, n, m)
```

## Arguments

- newdata:

  an `n x m x n_dyad` array (or `n x m` matrix) of NEW dyadic
  covariates. If its 3rd-dim names match the model's dyadic coefficient
  names they are mapped by name; otherwise the slices are taken in order
  as the model's dyadic covariates.

- model_names:

  the model's coefficient names (`colnames(fit$BETA)` /
  `dimnames(fit$BETA)[[2]]`), in canonical layout.

- fitted_design:

  the model's fitted design array (`fit$X` or `Xlist[[t]]`); used to
  hold intercept + nodal slices fixed. May be `NULL` only when the model
  has no nodal covariates.

- n, m:

  output dimensions.

## Value

an `n x m x length(model_names)` array whose slices are in `model_names`
order: dyadic slices from `newdata`, all other slices held at their
fitted values (intercept = 1).
