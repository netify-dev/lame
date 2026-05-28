# Refit a fast AME model with a warm start

Refits an AME model by iterative block coordinate descent, initialised
(“warm-started”) from an existing `ame_als` fit. This is the workhorse
of
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md):
every bootstrap replicate is refit from the original point estimate
rather than from a cold random start, which prevents replicates from
converging to different local optima or rotations and is essential for
meaningful bootstrap standard errors.

## Usage

``` r
ame_als_refit(
  object,
  Y_new = NULL,
  X_new = NULL,
  Z_new = NULL,
  max_iter = 30,
  tol = 1e-05,
  verbose = FALSE
)
```

## Arguments

- object:

  an `ame_als` object supplying the warm-start values (`mu`, `beta`,
  `a`, `b`, `U`, `V`) and the model configuration (`family`, `mode`,
  `symmetric`, `R`).

- Y_new:

  optional canonical outcome array `[n_row, n_col, T]` to refit on.
  Defaults to `object$Y`.

- X_new:

  optional canonical design array `[n_row, n_col, p, T]`. Defaults to
  `object$X`.

- Z_new:

  optional canonical *working-response* array `[n_row, n_col, T]`. When
  supplied, the model is refit by a Gaussian block coordinate descent
  directly on `Z_new`, bypassing the family transform / IRLS
  reweighting. This is used by the parametric bootstrap of a non-normal
  `transform` fit, whose estimator is a Gaussian fit to a fixed
  transformed response. `Y_new` is still used for the observed-cell
  pattern.

- max_iter:

  maximum block coordinate descent iterations (default 30; fewer are
  needed than for a cold start).

- tol:

  convergence tolerance (default `1e-5`).

- verbose:

  logical; print progress (default `FALSE`).

## Value

An object of class `"ame_als"`; see
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md).

## Details

The estimation algorithm is the iterative block coordinate descent
estimator of the Social Influence Regression model of Hoff & Minhas
([`sir::sir_alsfit()`](https://rdrr.io/pkg/sir/man/sir_alsfit.html)),
adapted to the AME model. Warm-starting the bootstrap replicates is the
adaptation introduced when this routine was ported to sibling packages
and is reproduced here.

## References

Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
Influence Regression. *Political Analysis*. The iterative block
coordinate descent estimator refit here originates with that work
(implemented in
[`sir::sir_alsfit()`](https://rdrr.io/pkg/sir/man/sir_alsfit.html)).

## See also

[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md),
[`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md).

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
refit <- ame_als_refit(fit, verbose = FALSE)
coef(refit)
#> intercept 
#> 0.1243366 
```
