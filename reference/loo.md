# Generic dispatcher for loo / waic on ame / lame fits

Lightweight S3 generics so calls of the form `loo(fit)` dispatch through
R's S3 system even when the loo package is not loaded. When loo *is*
loaded, its generic resolves first; this fallback only fires for
bare-namespace use.

## Usage

``` r
loo(x, ...)

waic(x, ...)
```

## Arguments

- x:

  a fitted `ame` or `lame` object with `$log_lik`.

- ...:

  passed to the relevant method.

## Value

A `loo` / `waic` object.
