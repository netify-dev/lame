# Simulate Z given fixed rank nomination data

Simulates a random latent matrix Z given its expectation, dyadic
correlation and fixed rank nomination data

## Usage

``` r
rZ_frn_fc(Z, EZ, rho, Y, YL, odmax, odobs)
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

- odmax:

  a scalar or vector giving the maximum number of nominations for each
  individual

- odobs:

  observed outdegree

## Value

a square matrix, the new value of Z

## Details

simulates Z under the constraints (1) Y\\i,j\\\>Y\\i,k\\ =\>
Z\\i,j\\\>Z\\i,k\\ , (2) Y\\i,j\\\>0 =\> Z\\i,j\\\>0 , (3) Y\\i,j\\=0 &
odobs\\i\\\<odmax\\i\\ =\> Z\\i,j\\\<0

## Author

Peter Hoff
