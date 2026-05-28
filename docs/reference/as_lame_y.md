# Convert a graph object to a lame-ready adjacency matrix

Convenience wrapper that takes an igraph or network object and returns
the kind of matrix that
[`ame`](https://netify-dev.github.io/lame/reference/ame.md) /
[`lame`](https://netify-dev.github.io/lame/reference/lame.md) expect: a
numeric adjacency matrix with `NA` on the diagonal (self-ties not
modelled) and (if available) actor names preserved on the row/column
dimnames.

## Usage

``` r
as_lame_y(x, na_diag = TRUE)
```

## Arguments

- x:

  an `igraph`, `network`, `matrix`, or `data.frame` representing an
  adjacency / sociomatrix.

- na_diag:

  logical: replace the diagonal with `NA`? Defaults to `TRUE` for square
  (unipartite) matrices and is ignored for rectangular (bipartite)
  matrices.

## Value

A numeric matrix suitable to pass as `Y` to
[`ame`](https://netify-dev.github.io/lame/reference/ame.md) /
[`lame`](https://netify-dev.github.io/lame/reference/lame.md).

## Details

Plain matrices and data.frames are accepted too; data.frames are coerced
to a numeric matrix and only the diagonal is rewritten to `NA`.

## Examples

``` r
# \donttest{
# mat <- as_lame_y(igraph::sample_smallworld(1, 12, 2, 0.1))
# ame(mat, family = "binary", R = 1, burn = 50, nscan = 200,
#     odens = 5, verbose = FALSE, symmetric = TRUE)
# }
```
