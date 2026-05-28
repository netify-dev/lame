# Generic dispatcher for posterior::as_draws on lame fits

Lightweight S3 generic so calls of the form `as_draws(fit)` dispatch
through R's S3 system even when the posterior package is not loaded.
When posterior *is* loaded, its generic of the same name is resolved
first by R's namespace search; this fallback only fires for
bare-namespace use.

## Usage

``` r
as_draws(x, ...)
```

## Arguments

- x:

  a fitted `ame` or `lame` object.

- ...:

  passed to the relevant method.

## Value

An object dispatched by the relevant method (typically a
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
for an `ame` / `lame` fit).
