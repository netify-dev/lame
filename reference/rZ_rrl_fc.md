# Simulate Z given relative rank nomination data

Simulates a random latent matrix Z given its expectation, dyadic
correlation and relative rank nomination data

## Usage

``` r
rZ_rrl_fc(Z, EZ, rho, Y, YL)
```

## Arguments

- Z:

  a square matrix, the current value of Z

- EZ:

  expected value of Z

- rho:

  dyadic correlation

- Y:

  square matrix of ranked nomination data

- YL:

  list of ranked individuals, from least to most preferred in each row

## Value

a square matrix, the new value of Z

## Details

simulates Z under the constraints (1) Y\\i,j\\\>Y\\i,k\\ =\>
Z\\i,j\\\>Z\\i,k\\

## Author

Peter Hoff
