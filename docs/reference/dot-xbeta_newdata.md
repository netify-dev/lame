# Name-safe linear predictor from a newdata prediction

Convenience wrapper: builds the full canonical design via
[`.build_full_design`](https://netify-dev.github.io/lame/reference/dot-build_full_design.md)
and contracts it against a named coefficient vector. Because the design
is rebuilt in the model's slice order, the positional contraction is
correct even with nodal covariates.

## Usage

``` r
.xbeta_newdata(beta, newdata, fitted_design, n, m)
```

## Arguments

- beta:

  named numeric coefficient vector (a single draw or the posterior
  mean).

- newdata:

  dyadic covariate array (see `.build_full_design`).

- fitted_design:

  fitted design array (`fit$X`); held-fixed source for intercept + nodal
  slices.

- n, m:

  output dimensions.

## Value

an `n x m` matrix of the linear predictor.
