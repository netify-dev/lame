# Describe the estimator behind a fitted object

Reports which estimation routine produced a fitted object and what
uncertainty information, if any, is available. Methods exist for
`ame_als` fits (the fast block coordinate descent estimator) and
`boot_ame` bootstrap objects.

## Usage

``` r
sampler_describe(object, verbose = TRUE)

# S3 method for class 'ame_als'
sampler_describe(object, verbose = TRUE)

# S3 method for class 'boot_ame'
sampler_describe(object, verbose = TRUE)

# S3 method for class 'ame'
sampler_describe(object, verbose = TRUE)

# S3 method for class 'lame'
sampler_describe(object, verbose = TRUE)
```

## Arguments

- object:

  a fitted object (`ame_als` or `boot_ame`).

- verbose:

  logical; if `TRUE` (default) print the description.

## Value

The input `object`, invisibly.
