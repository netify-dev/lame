# Update an AME / LAME fit

S3 method for [`update`](https://rdrr.io/r/stats/update.html) that
re-fits the model with modified arguments. Reuses the original call
recorded in `fit$call`. If the fit was produced with
`freeze_call = TRUE`, the data snapshot on `fit$data_snapshot` is used
in place of looking up names in the caller's environment.

## Usage

``` r
# S3 method for class 'ame'
update(object, ..., evaluate = TRUE)

# S3 method for class 'lame'
update(object, ..., evaluate = TRUE)
```

## Arguments

- object:

  A fitted `ame` or `lame` object.

- ...:

  Named arguments to overwrite in the original call (e.g.
  `nscan = 5000`, `dynamic_beta = TRUE`, `R = 2`).

- evaluate:

  Logical: if `TRUE` (default), evaluate the updated call and return the
  new fit; if `FALSE`, return the unevaluated call.

## Value

A new fitted object, or (when `evaluate = FALSE`) the modified call.

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary",
            burn = 5, nscan = 20, odens = 1, verbose = FALSE)
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
# toggle dynamic_beta on without rewriting the whole call
fit_dyn <- update(fit, dynamic_beta = "dyad")
#> Warning: `family` = "binary" but `Y` contains values other than 0/1.
#> ℹ `Y` will be thresholded to `1 * (Y > 0)`; if you meant counts, use "poisson",
#>   or "ordinal"/"normal" as appropriate.
dim(fit_dyn$BETA)  # 3-D now
#> [1] 20  4  4
# }
```
