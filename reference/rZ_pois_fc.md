# Gibbs update for latent variable in a Poisson AME model

Updates the latent variable Z in a Poisson AME model using a
Metropolis-Hastings step. The model assumes y\_{i,j} ~
Poisson(exp(z\_{i,j})) where z\_{i,j} is the latent variable
representing the log mean.

## Usage

``` r
rZ_pois_fc(Z, EZ, rho, s2, Y)
```

## Arguments

- Z:

  a square matrix, the current value of the latent variable

- EZ:

  expected value of Z (regression effects + random effects)

- rho:

  dyadic correlation

- s2:

  dyadic variance (overdispersion parameter)

- Y:

  square relational matrix of observed counts

## Value

updated value of Z

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
