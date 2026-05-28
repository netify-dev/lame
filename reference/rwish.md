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
#>           [,1]     [,2]     [,3]
#> [1,] 21.045373 23.53593 2.504484
#> [2,] 23.535935 44.32224 1.296290
#> [3,]  2.504484  1.29629 2.401330

S0*5
#>           [,1]      [,2]     [,3]
#> [1,] 20.829294 23.455492 2.424659
#> [2,] 23.455492 44.581153 1.320969
#> [3,]  2.424659  1.320969 2.356190

```
