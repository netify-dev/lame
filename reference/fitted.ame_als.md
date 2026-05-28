# Extract fitted values from a fast AME fit

Returns response-scale fitted values: a single matrix for a
cross-sectional (`ame_als`) fit, or a list of matrices for a
longitudinal (`lame_als`) fit.

## Usage

``` r
# S3 method for class 'ame_als'
fitted(object, ...)
```

## Arguments

- object:

  an `ame_als` object.

- ...:

  ignored.

## Value

A matrix or list of matrices of fitted values.
