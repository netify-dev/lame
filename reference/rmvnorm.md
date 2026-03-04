# Simulation from a multivariate normal distribution

Simulates a matrix where the rows are i.i.d. samples from a multivariate
normal distribution

## Usage

``` r
rmvnorm(n, mu, Sigma, Sigma.chol = NULL)
```

## Arguments

- n:

  sample size

- mu:

  multivariate mean vector

- Sigma:

  covariance matrix

- Sigma.chol:

  Cholesky factorization of `Sigma`

## Value

a matrix with `n` rows

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
