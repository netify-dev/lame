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
#>            [,1]      [,2]      [,3]
#> [1,]  26.216328 -9.074367 -11.97831
#> [2,]  -9.074367 46.079683  32.81260
#> [3,] -11.978313 32.812604  41.56388

S0*5
#>            [,1]      [,2]      [,3]
#> [1,]  26.297314 -8.697926 -12.11344
#> [2,]  -8.697926 46.784113  32.30276
#> [3,] -12.113443 32.302764  40.83354

```
