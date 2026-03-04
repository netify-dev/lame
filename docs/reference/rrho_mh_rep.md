# Metropolis update for dyadic correlation with independent replicate data

Metropolis update for dyadic correlation with independent replicate
data.

## Usage

``` r
rrho_mh_rep(E.T, rho, s2 = 1)
```

## Arguments

- E.T:

  Array of square residual relational matrix series. The third dimension
  of the array is for different replicates. Each slice of the array
  according to the third dimension is a square residual relational
  matrix.

- rho:

  current value of rho

- s2:

  current value of s2

## Value

a new value of rho

## Author

Peter Hoff, Yanjun He
