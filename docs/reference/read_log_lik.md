# Read the per-iteration log-lik matrix back from on-disk chunks

For fits run with `save_log_lik = "chunked"`, the per-iteration
pointwise log-likelihood lives in per-column-chunk binary files instead
of `fit$log_lik`. This helper reconstitutes the `[n_stored, n_obs]`
matrix in memory by issuing one `readBin` per chunk.
[`loo.lame()`](https://netify-dev.github.io/lame/reference/loo.ame.md) /
[`waic.lame()`](https://netify-dev.github.io/lame/reference/waic.ame.md)
call this transparently.

## Usage

``` r
read_log_lik(x)
```

## Arguments

- x:

  A fitted `ame`/`lame` object with `$log_lik_chunks` (i.e. fit with
  `save_log_lik = "chunked"`).

## Value

The full `[n_stored, n_obs]` double matrix.
