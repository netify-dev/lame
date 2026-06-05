# synthetic longitudinal binary relational data (latent-scale)

a synthetic 50-actor, 4-period dataset used by the vignettes. The
shipped `Y` array holds the *latent probit-scale predictor* (range
roughly -24 to 20), not 0/1 ties. To recover the binary tie indicator
the dataset name implies, threshold at zero:
`Y_bin <- (YX_bin_long$Y > 0) * 1`. If passed unthresholded to
`lame(..., family = "binary")` the fit will warn and silently apply the
same `Y > 0` threshold.

## Usage

``` r
data(YX_bin_long)
```

## Format

A list with two elements:

- Y:

  Numeric array `[50, 50, 4]`. Latent probit-scale predictor; threshold
  at 0 to get the binary tie indicator.

- X:

  Numeric array `[50, 50, 3, 4]` of dyadic covariates.

## Examples

``` r

data(YX_bin_long)
# threshold latent z to 0/1 before computing binary GOF stats
Yt <- 1 * (YX_bin_long$Y[, , 1] > 0)
diag(Yt) <- NA
gof_stats(Yt)
#>  sd.rowmean  sd.colmean    dyad.dep   cycle.dep   trans.dep 
#> 0.098854028 0.083875884 0.039461682 0.001059680 0.003684353 
```
