# Extract latent positions as a tidy data frame

Extracts multiplicative latent factor positions (U and V) from a fitted
`ame` or `lame` model and returns them as a tidy data frame suitable for
plotting and analysis. Optionally applies Procrustes alignment for
dynamic models and includes posterior standard deviations when posterior
samples are available.

## Usage

``` r
latent_positions(object, ...)

# S3 method for class 'ame'
latent_positions(object, align = FALSE, ...)

# S3 method for class 'lame'
latent_positions(object, align = TRUE, ...)
```

## Arguments

- object:

  A fitted `ame` or `lame` model object with R \> 0.

- ...:

  Additional arguments (currently unused).

- align:

  Logical. For dynamic models (`dynamic_uv = TRUE`), apply Procrustes
  alignment across time to remove rotational indeterminacy. Default is
  `FALSE` for `ame` objects and `TRUE` for `lame` objects.

## Value

A data frame with columns:

- actor:

  Character. Actor name (from rownames of U or V).

- dimension:

  Integer. Latent dimension index (1 to R).

- time:

  Character or NA. Time period label for dynamic models; `NA` for static
  (cross-sectional) models.

- value:

  Numeric. The posterior mean latent position.

- posterior_sd:

  Numeric. Posterior standard deviation of the latent position, or `NA`
  if posterior samples are not available. To enable, fit the model with
  `posterior_opts = posterior_options(save_UV = TRUE)`.

- type:

  Character. `"U"` for sender/row positions, `"V"` for receiver/column
  positions. Symmetric models have only `"U"`.

Returns a zero-row data frame with correct column names if R = 0.

## See also

[`procrustes_align`](https://netify-dev.github.io/lame/reference/procrustes_align.md)
for standalone Procrustes alignment,
[`uv_plot`](https://netify-dev.github.io/lame/reference/uv_plot.md) for
visualizing latent positions,
[`posterior_options`](https://netify-dev.github.io/lame/reference/posterior_options.md)
for enabling posterior sampling of U/V

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           burn = 5, nscan = 5, odens = 1, verbose = FALSE)
lp <- latent_positions(fit)
head(lp)
#>   actor dimension time       value posterior_sd type
#> 1 node1         1 <NA> -0.28119834           NA    U
#> 2 node2         1 <NA>  0.24841529           NA    U
#> 3 node3         1 <NA>  0.45288719           NA    U
#> 4 node4         1 <NA> -0.24169524           NA    U
#> 5 node5         1 <NA>  0.06922696           NA    U
#> 6 node6         1 <NA> -0.28624933           NA    U
# }
```
