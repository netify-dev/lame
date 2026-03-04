# Simulate Z based on a tobit model

Simulates a random latent matrix Z given its expectation, dyadic
correlation and a nonnegative relational matrix Y

## Usage

``` r
rZ_tob_fc(Z, EZ,rho,s2,Y)
```

## Arguments

- Z:

  a square matrix, the current value of Z

- EZ:

  expected value of Z

- rho:

  dyadic correlation

- s2:

  dyadic variance

- Y:

  square relational matrix with nonnegative entries

## Value

a square matrix, the new value of Z

## Author

Peter Hoff
