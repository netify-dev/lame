# binary relational data and covariates

a synthetic dataset that includes longitudinal binary relational data as
well as information on covariates

## Usage

``` r
data(YX_bin_long)
```

## Format

a list

## Examples

``` r
data(YX_bin_long)
gof_stats(YX_bin_long$Y[,,1]) 
#>   sd.rowmean   sd.colmean     dyad.dep    cycle.dep    trans.dep 
#>  1.261836080  1.007759298  0.060561322 -0.001188745 -0.001308371 
```
