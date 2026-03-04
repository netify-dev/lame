# Display memory usage information for AME models

Shows estimated memory usage for networks of given size. Memory
optimization is now automatic, so this is informational only.

## Usage

``` r
ame_memory_settings(n_nodes, R = 2)
```

## Arguments

- n_nodes:

  Number of nodes in network

- R:

  Rank of multiplicative effects (default: 2)

## Value

Invisibly returns memory estimates

## Details

Memory optimization is now user-controlled:

- Use remove_redundant=TRUE to remove EZ and UVPM matrices

- Use use_sparse_matrices=TRUE for sparse networks

- Both options can be combined for maximum memory savings

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
