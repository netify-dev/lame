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
#>           [,1]        [,2]        [,3]
#> [1,] 22.946304 24.40709516  2.45440657
#> [2,] 24.407095 43.17559537 -0.09257386
#> [3,]  2.454407 -0.09257386  0.90743519

S0*5
#>           [,1]       [,2]       [,3]
#> [1,] 22.969883 24.6612057  2.3931178
#> [2,] 24.661206 43.8108527 -0.1529072
#> [3,]  2.393118 -0.1529072  0.8938374

```
