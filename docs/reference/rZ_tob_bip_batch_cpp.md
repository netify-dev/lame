# Batch tobit Z sampling across all time periods (bipartite, rho=0)

Batch tobit Z sampling across all time periods (bipartite, rho=0)

## Usage

``` r
rZ_tob_bip_batch_cpp(Z, EZ, s2, Y)
```

## Arguments

- Z:

  3D array of latent values (nA x nB x T)

- EZ:

  3D array of expected values (nA x nB x T)

- s2:

  Dyadic variance

- Y:

  3D array of observed values (nA x nB x T)

## Value

Updated Z array
