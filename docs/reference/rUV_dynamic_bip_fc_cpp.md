# Bipartite dynamic UV Gibbs update

Replaces the nested R loops for bipartite dynamic_uv in lame.R. Uses
direct 2x2 inverse formula for common R=2 case.

## Usage

``` r
rUV_dynamic_bip_fc_cpp(U_cube, V_cube, E, G, rho_uv, sigma_uv, s2)
```

## Arguments

- U_cube:

  3D array (nA x RA x T)

- V_cube:

  3D array (nB x RB x T)

- E:

  3D array of residuals (nA x nB x T)

- G:

  Interaction matrix (RA x RB)

- rho_uv:

  AR(1) persistence parameter

- sigma_uv:

  Innovation standard deviation

- s2:

  Dyadic variance

## Value

List with updated U_cube, V_cube
