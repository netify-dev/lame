# Simulate an relational matrix based on a relative rank nomination scheme

Simulate an relational matrix based on a relative rank nomination scheme

## Usage

``` r
simY_rrl(EZ, rho, odobs, YO)
```

## Arguments

- EZ:

  a square matrix giving the expected value of the latent Z matrix

- rho:

  dyadic correlation

- odobs:

  a scalar or vector giving the observed number of nominations for each
  node

- YO:

  a square matrix identifying where missing values should be maintained

## Value

a square matrix, where higher values represent stronger relationships

## Author

Peter Hoff
