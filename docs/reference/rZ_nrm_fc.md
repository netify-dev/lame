# Simulate missing values in a normal AME model

Simulates missing values of a sociomatrix under a normal AME model

## Usage

``` r
rZ_nrm_fc(Z, EZ, rho,s2, Y)
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

  square relational matrix

## Value

a square matrix, equal to Y at non-missing values

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
