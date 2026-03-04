# ordinal relational data and covariates

a synthetic dataset that includes ordinal relational data as well as
information on seven covariates

## Usage

``` r
data(YX_ord)
```

## Format

The format is: List of 2 \$ Y: num \[1:100, 1:100\] NA 0 3 0 3 1 0 1 1 0
... \$ X: num \[1:100, 1:100, 1:7\] 1 1 1 1 1 1 1 1 1 1 ... ..- attr(\*,
"dimnames")=List of 3 .. ..\$ : NULL .. ..\$ : NULL .. ..\$ : chr
\[1:7\] "rgpa" "rsmoke" "cgpa" "csmoke" ...

## Examples

``` r
data(YX_ord)
gof_stats(YX_ord$Y)
#> sd.rowmean sd.colmean   dyad.dep  cycle.dep  trans.dep 
#> 0.64335425 0.27562901 0.68202697 0.03500558 0.04737175 
```
