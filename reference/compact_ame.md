# Optimize AME model output for memory efficiency

Optimizes AME model storage by converting matrices to sparse format.

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

When `use_sparse_matrices = TRUE`, converts YPM, APM, and BPM to sparse
matrices via the Matrix package. This is only beneficial when your
network is actually sparse (\< 10\\ networks, sparse storage will be
slower and may use more memory.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
