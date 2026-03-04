# Validate input data format for lame function

Internal validation function that checks the format and consistency of
input data for longitudinal AME models. Ensures that network data and
covariates are properly formatted as lists with consistent dimensions
across time periods.

## Usage

``` r
check_format(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL)
```

## Arguments

- Y:

  a list of T network matrices, where T is the number of time periods.
  Each element should be an n x m matrix representing the network at
  time t.

- Xdyad:

  an optional list of T dyadic covariates. Each element can be either an
  n x m matrix (single covariate) or an n x m x p array (p covariates).
  Must have the same length as Y if provided.

- Xrow:

  an optional list of T matrices of row/sender covariates. Each element
  should be an n x pr matrix where pr is the number of row covariates.
  Must have the same length as Y if provided.

- Xcol:

  an optional list of T matrices of column/receiver covariates. Each
  element should be an m x pc matrix where pc is the number of column
  covariates. Must have the same length as Y if provided.

## Value

Invisible TRUE if all checks pass. Throws an error with an informative
message if any validation fails.

## Details

Validates input data for longitudinal network analysis:

- Verifying Y is a non-empty list of matrices

- Checking that all covariates (if provided) are lists of appropriate
  length

- Validating data types for all elements

- Warning about dimension inconsistencies across time periods

The function uses informative error messages via the cli package to help
users identify and correct data formatting issues.

## Note

This is an internal function primarily used by
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) but
exported for advanced users who want to validate their data before model
fitting.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
