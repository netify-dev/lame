# Initialize dynamic additive effects with AR(1) structure

Initialize dynamic additive effects with AR(1) structure

## Usage

``` r
init_dynamic_ab_cpp(n, T, rho_ab, sigma_ab, mean_a = 0, mean_b = 0)
```

## Arguments

- n:

  Number of actors

- T:

  Number of time points

- rho_ab:

  AR(1) parameter

- sigma_ab:

  Innovation standard deviation

- mean_a:

  Mean for row effects

- mean_b:

  Mean for column effects

## Value

List with initialized a and b matrices
