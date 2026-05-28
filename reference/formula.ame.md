# formula() is not defined for an ame() / lame() fit

[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) and
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) use
named arguments (`Xrow`, `Xcol`, `Xdyad`, `symmetric`, etc.) rather than
a formula interface. The
[`summary()`](https://rdrr.io/r/base/summary.html) output renders a
pseudo-formula for diagnostic display, but it is not a valid R formula
and cannot be passed back to a fitting function.

## Usage

``` r
# S3 method for class 'ame'
formula(x, ...)

# S3 method for class 'lame'
formula(x, ...)
```

## Arguments

- x:

  an `ame` or `lame` fit.

- ...:

  ignored.

## Value

Never returns; raises an error.
