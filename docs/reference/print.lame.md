# Print method for LAME objects

Provides a concise print output for fitted LAME models. When the fit has
two or more dynamic effects active simultaneously
(`dynamic_uv + dynamic_ab + dynamic_beta` in any combination), a compact
joint table summarising the AR(1) hyperparameters of each dynamic block
is printed instead of a separate paragraph per component. Disable the
compact mode by passing `compact = FALSE`.

## Usage

``` r
# S3 method for class 'lame'
print(x, compact = TRUE, digits = 3, ...)
```

## Arguments

- x:

  an object of class "lame"

- compact:

  logical: when `TRUE` (default), use a compact joint table for fits
  with 2+ dynamic components; set to `FALSE` to force the per-component
  long-form display.

- digits:

  Number of digits to display in the compact table. Default 3.

- ...:

  additional arguments (not used)

## Value

the lame object invisibly

## Author

Peter Hoff, Cassy Dorff, Shahryar Minhas
