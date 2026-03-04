# Sample dynamic additive effects with AR(1) evolution

Sample dynamic additive effects with AR(1) evolution

## Usage

``` r
sample_dynamic_ab_cpp(
  a_current,
  b_current,
  Z_array,
  EZ_array,
  rho_ab,
  sigma_ab,
  Sab,
  symmetric
)
```

## Arguments

- a_current:

  Current 2D array of row effects (n x T)

- b_current:

  Current 2D array of column effects (n x T)

- Z_array:

  3D array of latent positions (n x n x T)

- EZ_array:

  3D array of expected values without additive effects (n x n x T)

- rho_ab:

  AR(1) parameter for additive effects

- sigma_ab:

  Innovation standard deviation

- Sab:

  Covariance matrix for a and b (2x2)

- symmetric:

  Whether the network is symmetric

## Value

List with updated a and b arrays
