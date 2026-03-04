# Goodness of fit statistics for bipartite networks

Calculates goodness of fit statistics specifically designed for
bipartite networks, evaluating degree heterogeneity and higher-order
dependencies.

## Usage

``` r
gof_stats_bipartite(Y)
```

## Arguments

- Y:

  a bipartite relational data matrix (nA x nB rectangular matrix) where
  Y\\i,j\\ represents the relationship from node i in set A to node j in
  set B. Missing values (NA) are allowed and will be handled
  appropriately.

## Value

A named numeric vector containing bipartite-specific goodness-of-fit
statistics:

- sd.rowmean:

  Standard deviation of row means. Measures the heterogeneity in
  out-degree from set A nodes (sender effects). Higher values indicate
  more variation in how active A nodes are.

- sd.colmean:

  Standard deviation of column means. Measures the heterogeneity in
  in-degree to set B nodes (receiver effects). Higher values indicate
  more variation in how popular B nodes are.

- four.cycles:

  Count of four-cycles (also called 4-paths or squares) in the bipartite
  network. A four-cycle occurs when two nodes from set A (e.g., i and k)
  both connect to the same two nodes in set B (e.g., j and l), forming a
  closed path: i-\>j-\>k-\>l-\>i. This measures the tendency for pairs
  of A-nodes to share multiple common B-node connections, capturing a
  form of clustering specific to bipartite networks. High four-cycle
  counts indicate that connections are not random but show patterns of
  shared preferences or co-occurrence. For example, in a user-item
  network, many four-cycles suggest that users who like one item tend to
  also like other items that co-occur with it.

## Details

For bipartite networks, reciprocity and triadic closure are not
meaningful concepts since edges only exist between the two node sets.
Instead, this function focuses on:

- Degree heterogeneity in both node sets

- Four-cycles as the simplest higher-order dependence pattern

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Create a random bipartite network
Y <- matrix(rnorm(10*12), 10, 12)

# Calculate GOF statistics
gof_stats_bipartite(Y)
#>   sd.rowmean   sd.colmean  four.cycles 
#>    0.2801100    0.4537885 2970.0000000 
# }
```
