# Censored binary nomination data and covariates

a synthetic dataset that includes relational data where the number of
nominations per row is censored at 10, along with information on eight
covariates

## Usage

``` r
data(YX_cbin)
```

## Format

The format is: List of 2 \$ Y: num \[1:100, 1:100\] NA 0 0 0 1 0 0 0 0 3
... \$ X: num \[1:100, 1:100, 1:8\] 1 1 1 1 1 1 1 1 1 1 ... ..- attr(\*,
"dimnames")=List of 3 .. ..\$ : NULL .. ..\$ : NULL .. ..\$ : chr
\[1:8\] "intercept" "rgpa" "rsmoke" "cgpa" ...

## Examples

``` r
data(YX_cbin)
gof_stats(YX_cbin$Y) 
#>  sd.rowmean  sd.colmean    dyad.dep   cycle.dep   trans.dep 
#> 0.038979274 0.028571901 0.539507967 0.004740685 0.009563239 
```
