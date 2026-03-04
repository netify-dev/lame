# Optimize AME model output for memory efficiency

Optimizes AME model storage by optionally removing redundant components
and/or using sparse matrices based on user preferences.

## Usage

``` r
compact_ame(fit, use_sparse_matrices = FALSE)
```

## Arguments

- fit:

  Fitted AME model

- use_sparse_matrices:

  Logical; whether to use sparse matrix storage (default: FALSE)

## Value

Memory-optimized AME model object

## Details

This function provides two optimization strategies:

1.  remove_redundant = TRUE: Removes redundant matrices (EZ, UVPM) that
    can be reconstructed if needed, keeping only essential components
    (BETA, VC, YPM). Recommended for networks with \> 100 nodes.

2.  use_sparse_matrices = TRUE: Converts matrices to sparse format. Only
    recommended when your network is actually sparse (\< 10\\ For dense
    networks, sparse storage will be slower and may use more memory.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
