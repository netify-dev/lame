# Batch binary Z sampling across all time periods (bipartite, rho=0)

Vectorized probit update for bipartite binary networks without dyadic
correlation.

## Usage

``` r
rZ_bin_bip_batch_cpp(Z, EZ, Y)
```

## Arguments

- Z:

  3D array of latent values (nA x nB x T)

- EZ:

  3D array of expected values (nA x nB x T)

- Y:

  3D array of observed values (nA x nB x T)

## Value

Updated Z array
