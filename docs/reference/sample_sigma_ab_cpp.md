# Sample innovation variance for dynamic additive effects

Sample innovation variance for dynamic additive effects

## Usage

``` r
sample_sigma_ab_cpp(
  a_mat,
  b_mat,
  rho_ab,
  symmetric,
  prior_shape = 2,
  prior_scale = 1
)
```

## Arguments

- a_mat:

  Matrix of row effects (n x T)

- b_mat:

  Matrix of column effects (n x T)

- rho_ab:

  AR(1) parameter

- symmetric:

  Whether the network is symmetric

- prior_shape:

  Shape parameter for inverse gamma prior

- prior_scale:

  Scale parameter for inverse gamma prior

## Value

Updated sigma_ab value
