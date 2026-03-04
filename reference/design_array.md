# Computes the design socioarray of covariate values

Computes the design socioarray of covariate values for an AME fit

## Usage

``` r
design_array(Xrow=NULL,Xcol=NULL,Xdyad=NULL,intercept=TRUE,n)
```

## Arguments

- Xrow:

  an n x pr matrix of row covariates

- Xcol:

  an n x pc matrix of column covariates

- Xdyad:

  an n x n x pd array of dyadic covariates

- intercept:

  logical

- n:

  number of rows/columns

## Value

an n x n x (pr+pc+pd+intercept) 3-way array

## Author

Peter Hoff
