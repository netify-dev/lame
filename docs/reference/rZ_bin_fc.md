# Simulate Z based on a probit model

Simulates a random latent matrix Z given its expectation, dyadic
correlation and a binary relational matrix Y

## Usage

``` r
rZ_bin_fc(Z, EZ, rho, Y)
```

## Arguments

- Z:

  a square matrix, the current value of Z

- EZ:

  expected value of Z

- rho:

  dyadic correlation

- Y:

  square binary relational matrix

## Value

a square matrix , the new value of Z

## Author

Peter Hoff
