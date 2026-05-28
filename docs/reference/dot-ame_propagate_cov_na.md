# Propagate covariate missingness into the response

For an AME fit, a dyad whose covariate value is missing cannot
contribute a covariate-coefficient observation. Rather than silently
imputing the missing covariate as 0 (which biases the coefficient when
the covariate is not mean-centred), the dyad is treated as an unobserved
tie and handled by data augmentation. This helper sets the affected `Y`
cells to `NA`;
[`lame`](https://netify-dev.github.io/lame/reference/lame.md) already
does this internally, and this brings
[`ame`](https://netify-dev.github.io/lame/reference/ame.md) into line.

## Usage

``` r
.ame_propagate_cov_na(Y, Xrow = NULL, Xcol = NULL, Xdyad = NULL)
```

## Arguments

- Y:

  an n x n response matrix.

- Xrow, Xcol, Xdyad:

  row, column and dyadic covariates (or `NULL`).

## Value

`Y` with covariate-missing cells set to `NA`.
