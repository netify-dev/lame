# Gibbs update for dyadic variance with independent replicate relational data

Gibbs update for dyadic variance with independent replicate relational
data

## Usage

``` r
rs2_rep_fc(E.T, rho)
```

## Arguments

- E.T:

  Array of square residual relational matrix series. The third dimension
  of the array is for different replicates. Each slice of the array
  according to the third dimension is a square residual relational
  matrix

- rho:

  current value of rho

## Value

a new value of s2

## Author

Peter Hoff, Yanjun He
