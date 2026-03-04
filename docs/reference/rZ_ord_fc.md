# Simulate Z given the partial ranks

Simulates a random latent matrix Z given its expectation, dyadic
correlation and partial rank information provided by W

## Usage

``` r
rZ_ord_fc(Z, EZ, rho, Y)
```

## Arguments

- Z:

  a square matrix, the current value of Z

- EZ:

  expected value of Z

- rho:

  dyadic correlation

- Y:

  matrix of ordinal data

## Value

a square matrix, the new value of Z

## Author

Peter Hoff
