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
#>            [,1]      [,2]       [,3]
#> [1,] 25.6460215  5.322573  0.8112129
#> [2,]  5.3225730 21.623964 11.5956678
#> [3,]  0.8112129 11.595668 27.6213424

S0*5
#>            [,1]      [,2]       [,3]
#> [1,] 25.9495551  5.008831  0.6264909
#> [2,]  5.0088310 22.073612 11.8918707
#> [3,]  0.6264909 11.891871 27.8402945

```
