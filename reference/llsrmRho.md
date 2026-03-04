# SRM log likelihood evaluated on a grid of rho-values

Calculation of the SRM log-likelihood over a grid of rho-values

## Usage

``` r
llsrmRho(Y, Sab, rhos, s2 = 1)
```

## Arguments

- Y:

  sociomatrix assumed to follow a mean-zero SRM distribution

- Sab:

  covariance matrix of additive effects

- rhos:

  vector of rho-values at which to calculate the log-likelihood

- s2:

  current value of s2

## Value

a vector of log-likelihood values

## Author

Peter Hoff
