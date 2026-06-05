# Resume a `lame()` MCMC run from a checkpoint

For a fit started with `lame(..., checkpoint_path = "X.rds")` that
terminated early (either via `max_seconds` or an external interruption),
continues the chain from the most recent checkpoint. The implementation
is pragmatic: it loads the saved RNG state and re-invokes
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) with the
original call arguments. The result is a fresh fit that picks up where
the previous one left off in the random-number stream; the underlying
MCMC counter restarts at 1.

## Usage

``` r
lame_resume(path, nscan_more = NULL, ..., .envir = NULL)
```

## Arguments

- path:

  Checkpoint file path (the one passed as `checkpoint_path` to
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md)).

- nscan_more:

  Optional integer, the new `nscan` for the continuation. If `NULL`, the
  original `nscan` is used.

- ...:

  Additional arguments forwarded to
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md)
  (override the saved values).

- .envir:

  Environment in which to evaluate the saved call's data arguments.
  Defaults to the caller's frame; used internally when
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md)
  forwards a resume.

## Value

A fitted `lame` object.

## Details

**Equivalent consolidated entry point.** `lame(resume_from = path, ...)`
short-circuits to this function with the user-supplied overrides
forwarded; pass `nscan = K` on the resume call to request `K` additional
stored draws. Both call shapes are supported. The
`lame_resume(path, ...)` form is safer for nested calls: the
consolidated form re-evaluates the saved
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) call in
[`parent.frame()`](https://rdrr.io/r/base/sys.parent.html), which is the
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) frame
whose required formal arguments (`Y`, `Xdyad`, etc.) were not supplied
on the resume call. Prefer `lame_resume(path, ...)` when calling from
inside other functions or from non-global scopes.

Use `nscan_more = K` to override the `nscan` value to `K` for the
continuation. Other arguments can be overridden by passing them to
`...`.

## Examples

``` r
# \donttest{
ck <- tempfile(fileext = ".rds")
data(YX_bin_list)
# short run with very-aggressive max_seconds to force early termination
fit1 <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
             nscan = 5000, burn = 50, odens = 5,
             checkpoint_path = ck, checkpoint_every = 50L,
             max_seconds = 0.5, verbose = FALSE)
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
if (isTRUE(fit1$terminated_early)) {
  fit2 <- lame_resume(ck, nscan_more = 200)
  dim(fit2$BETA)
}
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
#> [1] 40  4
# }
```
