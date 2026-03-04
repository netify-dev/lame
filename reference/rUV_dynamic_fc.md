# Gibbs sampling of dynamic U and V with AR(1) evolution

Updates latent factor positions U and V that evolve over time according
to an AR(1) process: u\_{i,t} = rho \* u\_{i,t-1} + epsilon\_{i,t}

## Usage

``` r
rUV_dynamic_fc(U, V, ET, rho_uv, sigma_uv, s2, shrink=TRUE, symmetric=FALSE)
```

## Arguments

- U:

  3D array of current U positions (n x R x T)

- V:

  3D array of current V positions (n x R x T)

- ET:

  3D array of residuals (n x n x T)

- rho_uv:

  AR(1) autoregressive parameter for latent positions

- sigma_uv:

  Innovation standard deviation for latent positions

- s2:

  dyadic variance

- shrink:

  whether to apply shrinkage (default TRUE)

- symmetric:

  whether the network is symmetric (default FALSE)

## Value

list with updated U and V arrays

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
