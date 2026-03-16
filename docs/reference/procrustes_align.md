# Procrustes alignment of latent positions across time

Aligns dynamic latent positions across time periods to remove arbitrary
rotational indeterminacy. This is essential for interpreting temporal
trajectories of latent positions, since the latent space is only
identified up to rotation at each time point.

Uses Procrustes rotation: at each time step t, finds the orthogonal
rotation matrix that best aligns U_t to U_t-1 (minimizing Frobenius
norm), then applies it sequentially from t=2 to T.

## Usage

``` r
procrustes_align(
  object = NULL,
  U = NULL,
  V = NULL,
  G = NULL,
  return_fit = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted `ame` or `lame` object, or NULL if raw arrays are provided
  via `U`/`V`/`G`.

- U:

  Optional 3D array `[n, R, T]` of sender latent positions. If `object`
  is provided and `U` is NULL, extracted from `object$U`.

- V:

  Optional 3D array `[n, R, T]` of receiver latent positions. For
  unipartite asymmetric models, aligned independently. For bipartite
  models, aligned jointly with U via the G interaction matrix.

- G:

  Optional `[R_row, R_col]` interaction matrix for bipartite models.
  Updated to maintain the invariant U G V'.

- return_fit:

  Logical. If TRUE and `object` is provided, returns a modified copy of
  the fit object with aligned latent positions. Default FALSE.

- ...:

  Additional arguments (currently unused).

## Value

If `return_fit = FALSE` (default): a list with components `U` (aligned
sender positions), `V` (aligned receiver positions, if applicable), and
`G` (updated interaction matrix, for bipartite).

If `return_fit = TRUE`: a copy of `object` with aligned latent positions
replacing the originals.

## Details

For unipartite networks, U and V are aligned independently using
separate Procrustes rotations. For symmetric networks, only U is present
and aligned.

For bipartite networks, U and V are aligned jointly: separate rotation
matrices are computed for U and V, and the G interaction matrix is
updated as `G_aligned = t(R_U) %*% G %*% R_V` to preserve the product
`U %*% G %*% t(V)`.

If U is a 2D matrix (static model with a single time point), it is
returned unchanged with an informational message.

## See also

[`latent_positions`](https://netify-dev.github.io/lame/reference/latent_positions.md)
for extracting aligned positions as a tidy data frame,
[`uv_plot`](https://netify-dev.github.io/lame/reference/uv_plot.md) for
visualizing latent positions

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, Xdyad = YX_bin_list$X, R = 2,
            family = "binary", dynamic_uv = TRUE,
            burn = 5, nscan = 5, odens = 1,
            verbose = FALSE)
aligned <- procrustes_align(fit)
str(aligned$U)  # aligned 3D array [n, R, T]
#>  num [1:50, 1:2, 1:4] 0.0512 0.083 0.0178 -0.0686 -0.079 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ : chr [1:50] "303" "304" "317" "321" ...
#>   ..$ : NULL
#>   ..$ : NULL
# }
```
