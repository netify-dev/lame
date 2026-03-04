# Goodness of fit statistics

Calculates goodness of fit statistics for relational data matrices,
evaluating second-order (dyadic) and third-order (triadic) dependence
patterns. Can handle both unipartite and bipartite networks.

## Usage

``` r
gof_stats(Y, mode = NULL, custom_gof = NULL)
```

## Arguments

- Y:

  a relational data matrix. For unipartite networks, a square n x n
  matrix where Y\\i,j\\ represents the relationship from node i to
  node j. For bipartite networks, an nA x nB matrix where Y\\i,j\\
  represents the relationship from node i in set A to node j in set B.
  Missing values (NA) are allowed and will be handled appropriately.

- mode:

  character string specifying the network type: "unipartite" or
  "bipartite". If NULL (default), attempts to infer from matrix
  dimensions (rectangular = bipartite, square = unipartite). Note:
  square bipartite networks must specify mode="bipartite".

- custom_gof:

  optional function or list of functions to compute custom GOF
  statistics. Each function should take Y as input and return a named
  numeric value or vector. Custom statistics will be added to the
  standard statistics.

## Value

A named numeric vector containing goodness-of-fit statistics. For
unipartite networks:

- sd.rowmean:

  Standard deviation of row means. Measures the heterogeneity in
  out-degree centrality (sender effects).

- sd.colmean:

  Standard deviation of column means. Measures the heterogeneity in
  in-degree centrality (receiver effects).

- dyad.dep:

  Dyadic dependence/reciprocity correlation.

- cycle.dep:

  Cyclic/transitive triadic dependence.

- trans.dep:

  Transitive triadic dependence.

For bipartite networks:

- sd.rowmean:

  Standard deviation of row means (sender heterogeneity).

- sd.colmean:

  Standard deviation of column means (receiver heterogeneity).

- four.cycles:

  Count of four-cycles in the bipartite network.

If custom_gof is provided, additional statistics will be included with
their user-specified names.

## Details

The function computes network statistics that capture different aspects
of network structure beyond simple density. These statistics are
particularly useful for:

- Model checking: comparing observed statistics to those from simulated
  networks

- Model selection: choosing between models that better capture network
  dependencies

- Descriptive analysis: summarizing key structural features of the
  network

For bipartite networks with square dimensions (nA = nB), you must
explicitly specify mode="bipartite" to ensure correct statistics are
calculated.

Missing values in Y are handled by pairwise deletion for correlations
and are excluded from matrix products in triadic calculations.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
data(YX_nrm) 

# Auto-detect unipartite
gof_stats(YX_nrm$Y) 
#> sd.rowmean sd.colmean   dyad.dep  cycle.dep  trans.dep 
#> 0.92646818 0.27555881 0.66792884 0.06139376 0.07380099 

# Explicitly specify mode for square bipartite
# Y_bip <- matrix(rnorm(100), 10, 10) # 10x10 bipartite
# gof_stats(Y_bip, mode = "bipartite")

# Custom GOF function
# my_stat <- function(Y) { c(my_measure = sum(Y > 0, na.rm=TRUE)) }
# gof_stats(YX_nrm$Y, custom_gof = my_stat)
```
