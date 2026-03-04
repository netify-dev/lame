# Edgelist to sociomatrix

Construction of a sociomatrix from an edgelist

## Usage

``` r
el2sm(el,directed=TRUE,nadiag=all(el[,1]!=el[,2]))
```

## Arguments

- el:

  a matrix in which each row contains the indices of an edge and
  possibly the weight for the edge

- directed:

  if FALSE, then a relation is placed in both entry ij and ji of the
  sociomatrix, for each edge ij (or ji)

- nadiag:

  put NAs on the diagonal

## Value

a sociomatrix

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
