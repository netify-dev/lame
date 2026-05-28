# S3 generic for `glance`

Light-weight fallback so `glance(fit)` dispatches through S3 even when
broom or generics is not loaded.

## Usage

``` r
glance(x, ...)
```

## Arguments

- x:

  An object to glance at.

- ...:

  Passed to the method.

## Value

A one-row data frame.
