# Simulation from a Wishart distribution

Simulates a random Wishart-distributed matrix

## Usage

``` r
rwish(S0, nu = dim(S0)[1] + 2)
```

## Arguments

- S0:

  a positive definite matrix

- nu:

  a positive integer

## Value

a positive definite matrix

## Author

Peter Hoff

## Examples

``` r
## The expectation is S0*nu

S0<-rwish(diag(3)) 

SS<-matrix(0,3,3) 
for(s in 1:1000) { SS<-SS+rwish(S0,5) }

SS/s 
#>           [,1]      [,2]      [,3]
#> [1,]  28.56606 -16.86061 -28.42098
#> [2,] -16.86061  41.98580  38.50025
#> [3,] -28.42098  38.50025  47.14915

S0*5
#>           [,1]      [,2]      [,3]
#> [1,]  28.86048 -17.20189 -28.90535
#> [2,] -17.20189  43.45814  39.52976
#> [3,] -28.90535  39.52976  47.99669

```
