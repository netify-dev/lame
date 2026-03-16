# Goodness of fit statistics for unipartite networks

Calculates goodness of fit statistics specifically for unipartite
(square) networks, evaluating second-order (dyadic) and third-order
(triadic) dependence patterns.

## Usage

``` r
gof_stats_unipartite(Y)
```

## Arguments

- Y:

  a square n x n relational data matrix where `Y[i,j]` represents the
  relationship from node i to node j. Missing values (NA) are allowed
  and will be handled appropriately. Diagonal values are typically NA
  for non-self-loop networks.

## Value

A named numeric vector containing five goodness-of-fit statistics:

- sd.rowmean:

  Standard deviation of row means. Measures the heterogeneity in
  out-degree centrality (sender effects).

- sd.colmean:

  Standard deviation of column means. Measures the heterogeneity in
  in-degree centrality (receiver effects).

- dyad.dep:

  Dyadic dependence/reciprocity correlation. Measures the correlation
  between `Y[i,j]` and `Y[j,i]`, capturing reciprocity patterns.

- cycle.dep:

  Cyclic triadic dependence. Measures the tendency for directed cycles
  (i-\>j-\>k-\>i) in the network.

- trans.dep:

  Transitive triadic dependence. Measures the tendency for transitivity
  (if i-\>j and j-\>k, then i-\>k) in the network.

## Details

This function computes network statistics that capture different aspects
of network structure beyond simple density. These statistics are
particularly useful for evaluating how well a model captures the
observed network patterns.

The dyadic dependence statistic captures reciprocity - the tendency for
relationships to be mutual. The triadic statistics capture different
forms of triadic closure that are common in social networks.

Missing values in Y are handled by pairwise deletion for correlations
and are excluded from matrix products in triadic calculations.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# Create a random unipartite network
Y <- matrix(rnorm(100), 10, 10)
diag(Y) <- NA
gof_stats_unipartite(Y)
#>  sd.rowmean  sd.colmean    dyad.dep   cycle.dep   trans.dep 
#>  0.34918446  0.45012473 -0.26811531 -0.06644152  0.01049362 
```
