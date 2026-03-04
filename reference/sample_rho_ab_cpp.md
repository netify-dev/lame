# Sample AR(1) parameter for dynamic additive effects

Sample AR(1) parameter for dynamic additive effects

## Usage

``` r
sample_rho_ab_cpp(a_mat, b_mat, sigma_ab, rho_current, symmetric)
```

## Arguments

- a_mat:

  Matrix of row effects (n x T)

- b_mat:

  Matrix of column effects (n x T)

- sigma_ab:

  Innovation standard deviation

- rho_current:

  Current value of rho

- symmetric:

  Whether the network is symmetric

## Value

Updated rho value
