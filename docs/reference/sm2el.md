# Sociomatrix to edgelist

Construction of an edgelist from a sociomatrix

## Usage

``` r
sm2el(sm,directed=TRUE)
```

## Arguments

- sm:

  a sociomatrix with possibly valued relations

- directed:

  if TRUE, only use the upper triangular part of the matrix to enumerate
  edges

## Value

an edglist

## Author

Peter Hoff

## Examples

``` r
Y<-matrix(rpois(10*10,.5),10,10) ; diag(Y)<-NA
E<-sm2el(Y) 
el2sm(E) - Y 
#>     1  2  3  4  5  6  7  8  9 10
#> 1  NA  0  0  0  0  0  0  0  0  0
#> 2   0 NA  0  0  0  0  0  0  0  0
#> 3   0  0 NA  0  0  0  0  0  0  0
#> 4   0  0  0 NA  0  0  0  0  0  0
#> 5   0  0  0  0 NA  0  0  0  0  0
#> 6   0  0  0  0  0 NA  0  0  0  0
#> 7   0  0  0  0  0  0 NA  0  0  0
#> 8   0  0  0  0  0  0  0 NA  0  0
#> 9   0  0  0  0  0  0  0  0 NA  0
#> 10  0  0  0  0  0  0  0  0  0 NA
```
