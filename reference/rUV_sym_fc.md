# Gibbs sampling of U and V

A Gibbs sampler for updating the multiplicative effect matrices U and V
in the symmetric case. In this case `U%*%t(V)` is symmetric, so this is
parameterized as `V=U%*%L` where `L` is the diagonal matrix of
eigenvalues of `U%*%t(V)`.

## Usage

``` r
rUV_sym_fc(E, U, V, s2 = 1, shrink=TRUE)
```

## Arguments

- E:

  square residual relational matrix

- U:

  current value of U

- V:

  current value of V

- s2:

  dyadic variance

- shrink:

  adaptively shrink the factors with a hierarchical prior

## Value

- U:

  a new value of U

- V:

  a new value of V

## Author

Peter Hoff

## Examples

``` r
U0<-matrix(rnorm(30,2),30,2) ; V0<-U0%*%diag(c(3,-2)) 
E<- U0%*%t(V0) + matrix(rnorm(30^2),30,30) 
rUV_sym_fc 
#> function (E, U, V, s2 = 1, shrink = TRUE) 
#> {
#>     n <- nrow(U)
#>     uLoopIDs <- as.integer(rep(sample(1:n), 4) - 1)
#>     rUV_sym_fc_cpp(E, U, V, s2, shrink, uLoopIDs)
#> }
#> <bytecode: 0x55ad5edd1768>
#> <environment: namespace:lame>
```
