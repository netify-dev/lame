# Batch normal Z sampling across all time periods

Replaces the per-time-period R loop with a single C++ call that loops
internally, reducing R-to-C++ transition overhead.

## Usage

``` r
rZ_nrm_batch_cpp(Z, EZ, rho, s2, Y)
```

## Arguments

- Z:

  3D array of latent values (n x n x T)

- EZ:

  3D array of expected values (n x n x T)

- rho:

  Dyadic correlation parameter

- s2:

  Dyadic variance

- Y:

  3D array of observed values (n x n x T)

## Value

List with updated Z and E_nrm (residuals)
