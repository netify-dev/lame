# Run AME model with multiple parallel chains

Run AME model with multiple parallel chains

## Usage

``` r
ame_parallel(
  Y,
  n_chains = 4,
  cores = n_chains,
  combine_method = c("pool", "list"),
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

- ...:

  Additional arguments passed to ame()

## Value

If combine_method = "pool": A single ame object with pooled chains If
combine_method = "list": A list of ame objects, one per chain

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
                             print = FALSE)
#> 
#> ── Running 2 chains sequentially ──
#> 
#> Starting chain 1
#> Completed chain 1
#> Starting chain 2
#> Completed chain 2
#> Combining chains...
#> 
#> ── MCMC Convergence Diagnostics 
#> Number of chains: 2
#> Samples per chain: "100, 100"
#> ✔ All parameters converged (R-hat < 1.1)
# }
```
