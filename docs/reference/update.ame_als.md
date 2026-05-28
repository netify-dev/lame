# Update an `ame_als` / `lame_als` fit

S3 method for [`update`](https://rdrr.io/r/stats/update.html) that
re-fits an ALS fit with modified arguments. Two routes:

- if `...` contains only the warm-start refit knobs (`Y_new`, `X_new`,
  `Z_new`, `max_iter`, `tol`, `verbose`), it delegates to
  [`ame_als_refit`](https://netify-dev.github.io/lame/reference/ame_als_refit.md)
  for a fast warm-started refit;

- otherwise (e.g. changing `R`, `family`, `lambda`, `bootstrap`) it
  re-evaluates the original `object$call` with the overrides, matching
  the
  [`update.ame`](https://netify-dev.github.io/lame/reference/update.ame.md)
  semantics on the MCMC side.

The split exists because
[`ame_als_refit()`](https://netify-dev.github.io/lame/reference/ame_als_refit.md)
is a genuinely warm-started single-fit refit (faster, narrower argument
set), while changing `R` / `family` requires a cold-start refit through
[`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md).

## Usage

``` r
# S3 method for class 'ame_als'
update(object, ..., evaluate = TRUE)

# S3 method for class 'lame_als'
update(object, ..., evaluate = TRUE)
```

## Arguments

- object:

  A fitted `ame_als` / `lame_als` object.

- ...:

  Named arguments. See
  [`ame_als_refit`](https://netify-dev.github.io/lame/reference/ame_als_refit.md)
  for the warm-start argument list and
  [`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md)
  for the cold-start argument list.

- evaluate:

  Logical: if `TRUE` (default), evaluate the updated call and return the
  new fit; if `FALSE`, return the unevaluated call (cold-start path
  only).

## Value

A new fitted `ame_als` object, or (when `evaluate = FALSE` on the
cold-start path) the modified call.

## Examples

``` r
# \donttest{
Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
# warm-start refit with the same arguments (fast)
fit_w <- update(fit, max_iter = 50)
# cold-start refit with a different R (full restart)
fit_R2 <- update(fit, R = 2)
# }
```
