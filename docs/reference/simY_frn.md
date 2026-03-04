# Simulate an relational matrix based on a fixed rank nomination scheme

Simulate an relational matrix based on a fixed rank nomination scheme

## Usage

``` r
simY_frn(EZ, rho, odmax, YO)
```

## Arguments

- EZ:

  a square matrix giving the expected value of the latent Z matrix

- rho:

  dyadic correlation

- odmax:

  a scalar or vector giving the maximum number of nominations for each
  node

- YO:

  a square matrix identifying where missing values should be maintained

## Value

a square matrix, where higher values represent stronger relationships

## Author

Peter Hoff
