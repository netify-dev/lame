# Simulate Z given fixed rank nomination data

Simulates a random latent matrix Z given its expectation, dyadic
correlation and censored binary nomination data

## Usage

``` r
rZ_cbin_fc(Z, EZ, rho, Y, odmax, odobs)
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

- odmax:

  a scalar or vector giving the maximum number of nominations for each
  individual

- odobs:

  observed outdegree

## Value

a square matrix, the new value of Z

## Details

simulates Z under the constraints (1) Y\\i,j\\=1, Y\\i,k\\=0 =\>
Z\\i,j\\\>Z\\i,k\\ , (2) Y\\i,j\\=1 =\> Z\\i,j\\\>0 , (3) Y\\i,j\\=0 &
odobs\\i\\\<odmax\\i\\ =\> Z\\i,j\\\<0

## Author

Peter Hoff
