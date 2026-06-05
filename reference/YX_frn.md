# Fixed rank nomination data and covariates

a synthetic dataset that includes fixed rank nomination data as well as
information on eight covariates

## Usage

``` r
data(YX_frn)
```

## Format

The format is: List of 2 \$ Y: num \[1:100, 1:100\] NA 0 0 0 1 0 0 0 0 3
... \$ X: num \[1:100, 1:100, 1:8\] 1 1 1 1 1 1 1 1 1 1 ... ..- attr(\*,
"dimnames")=List of 3 .. ..\$ : NULL .. ..\$ : NULL .. ..\$ : chr
\[1:8\] "intercept" "rgpa" "rsmoke" "cgpa" ...

## Examples

``` r

data(YX_frn)
gof_stats(YX_frn$Y) 
#>  sd.rowmean  sd.colmean    dyad.dep   cycle.dep   trans.dep 
#> 0.216963190 0.152754547 0.544539588 0.008276425 0.015688679 
```
