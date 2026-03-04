# row-specific ordinal relational data and covariates

a synthetic dataset that includes row-specific ordinal relational data
as well as information on five covariates

## Usage

``` r
data(YX_rrl)
```

## Format

The format is: List of 2 \$ Y: num \[1:100, 1:100\] NA 0 3 0 3 1 0 1 1 0
... \$ X: num \[1:100, 1:100, 1:5\] 1 1 1 1 1 1 1 1 1 1 ... ..- attr(\*,
"dimnames")=List of 3 .. ..\$ : NULL .. ..\$ : NULL .. ..\$ : chr
\[1:5\] "cgpa" "csmoke" "igrade" "ismoke" ...

## Examples

``` r
data(YX_rrl)
gof_stats(YX_rrl$Y)
#>  sd.rowmean  sd.colmean    dyad.dep   cycle.dep   trans.dep 
#> 0.214793558 0.151227002 0.544713565 0.005015426 0.010577213 
```
