# S3 generic for `tidy`

Light-weight fallback generic so that calls of the form `tidy(fit)`
dispatch through R's S3 system even when the broom or generics packages
are not loaded. When either is loaded, its generic resolves first; this
generic only fires for bare-namespace use.

## Usage

``` r
tidy(x, ...)
```

## Arguments

- x:

  An object to tidy.

- ...:

  Passed to the relevant method.

## Value

A data frame; method-specific schema.
