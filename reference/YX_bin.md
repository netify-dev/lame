# binary relational data and covariates

a synthetic dataset that includes binary relational data as well as
information on eight covariates

## Usage

``` r
data(YX_bin)
```

## Format

The format is: List of 2 \$ Y: num \[1:100, 1:100\] NA 0 0 0 0 0 0 0 0 1
... \$ X: num \[1:100, 1:100, 1:8\] 1 1 1 1 1 1 1 1 1 1 ... ..- attr(\*,
"dimnames")=List of 3 .. ..\$ : NULL .. ..\$ : NULL .. ..\$ : chr
\[1:8\] "intercept" "rgpa" "rsmoke" "cgpa" ...

## Examples

``` r
data(YX_bin)
gof_stats(YX_bin$Y) 
#> Note: Square matrix assumed to be unipartite. Use mode='bipartite' for square bipartite networks.
#> sd.rowmean sd.colmean   dyad.dep  cycle.dep  trans.dep 
#> 0.04927874 0.03068428 0.52475273 0.01596157 0.01788281 
```
