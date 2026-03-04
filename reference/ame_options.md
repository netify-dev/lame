# AME model fitting options

Configure options for fitting AME models. Memory efficiency is now
handled automatically based on network size.

## Usage

``` r
ame_options(
  parallel_chains = 1,
  verbose = TRUE,
  odens = 25,
  use_sparse_matrices = FALSE
)
```

## Arguments

- parallel_chains:

  Number of parallel chains to run (default: 1)

- verbose:

  Logical; print progress information (default: TRUE)

- odens:

  Output density - save every odens iterations (default: 25)

- use_sparse_matrices:

  Logical; use sparse matrices for storing results (default: FALSE). Set
  to TRUE if your network is actually sparse (many zero/NA entries) and
  memory is a concern.

## Value

List of options to pass to ame()

## Details

Memory optimization features:

- Redundant matrices (EZ, UVPM) are never stored - they can be
  reconstructed if needed

- use_sparse_matrices = TRUE: Converts large matrices to sparse format

- Posterior samples are always thinned appropriately

When to use sparse matrices:

- Your network has \< 10\\

- Memory usage is a critical concern

- You're willing to trade computational speed for memory efficiency

Note: For dense networks (most edges observed), sparse matrices will be
slower and may use MORE memory than dense storage.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Configure options
opts <- ame_options(verbose = TRUE, odens = 25)
opts
#> $parallel_chains
#> [1] 1
#> 
#> $verbose
#> [1] TRUE
#> 
#> $odens
#> [1] 25
#> 
#> $use_sparse_matrices
#> [1] FALSE
#> 
# }
```
