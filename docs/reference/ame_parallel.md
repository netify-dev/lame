# Run AME model with multiple parallel chains

Run AME model with multiple parallel chains

## Usage

``` r
ame_parallel(
  Y,
  n_chains = 4,
  cores = n_chains,
  combine_method = c("pool", "list"),
  fitter = c("auto", "ame", "lame"),
  ...
)
```

## Arguments

- Y:

  Network data matrix

- n_chains:

  Number of parallel chains to run (default = 4)

- cores:

  Number of CPU cores to use for parallel processing. Default is
  n_chains. Use 1 for sequential processing.

- combine_method:

  Method for combining chains: "pool" (default) or "list"

- fitter:

  Which fitter to use: `"auto"` (default; pick
  [`lame`](https://netify-dev.github.io/lame/reference/lame.md) when `Y`
  is a list / 3-D array with \>1 time slice, else
  [`ame`](https://netify-dev.github.io/lame/reference/ame.md)), `"ame"`
  (force cross-sectional), or `"lame"` (force longitudinal – required
  for `dynamic_*` flags).

- ...:

  Additional arguments passed to
  [`ame()`](https://netify-dev.github.io/lame/reference/ame.md) or
  [`lame()`](https://netify-dev.github.io/lame/reference/lame.md).

## Value

If combine_method = "pool": A single ame object with pooled chains If
combine_method = "list": A list of ame objects, one per chain

## See also

[`lame_multi`](https://netify-dev.github.io/lame/reference/lame_multi.md)
for the unrelated multi-*panel* wrapper that fits K *distinct* networks
with shared regression coefficients. `ame_parallel` / `lame_parallel`
run K MCMC chains of the *same* model (for R-hat / ESS diagnostics);
`lame_multi` runs one chain across K panels with pooled beta.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Run 2 chains sequentially
data(YX_nrm)
fit_parallel <- ame_parallel(YX_nrm$Y, Xdyad = YX_nrm$X,
                             n_chains = 2, cores = 1,
                             nscan = 100, burn = 10, odens = 1,
                             verbose = FALSE)
#> 
#> ── Running 2 chains sequentially ──
#> 
#> Starting chain 1 (`ame()`)
#> Completed chain 1
#> Starting chain 2 (`ame()`)
#> Completed chain 2
#> Combining chains...
#> 
#> ── MCMC Convergence Diagnostics 
#> Number of chains: 2
#> Samples per chain: "100, 100"
#> Diagnostics cover regression coefficients and variance components.
#> ✔ All parameters converged (R-hat < 1.1)
# }
```
