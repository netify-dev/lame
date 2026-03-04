# normal relational data and covariates

a synthetic dataset that includes continuous (normal) relational data as
well as information on eight covariates

## Usage

``` r
data(YX_nrm)
```

## Format

The format is: List of 2 \$ Y: num \[1:100, 1:100\] NA -4.05 -0.181
-3.053 -1.579 ... \$ X: num \[1:100, 1:100, 1:8\] 1 1 1 1 1 1 1 1 1 1
... ..- attr(\*, "dimnames")=List of 3 .. ..\$ : NULL .. ..\$ : NULL ..
..\$ : chr \[1:8\] "intercept" "rgpa" "rsmoke" "cgpa" ...

## Examples

``` r
data(YX_nrm)
gof_stats(YX_nrm$Y)
#> sd.rowmean sd.colmean   dyad.dep  cycle.dep  trans.dep 
#> 0.92646818 0.27555881 0.66792884 0.06139376 0.07380099 

```
