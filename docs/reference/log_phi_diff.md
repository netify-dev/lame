# Numerically stable log(Phi(hi) - Phi(lo))

Computes log of the standard-normal CDF difference. Naive
`log(pnorm(hi) - pnorm(lo))` returns `-Inf` whenever both arguments are
above ~7 or below ~-7 (catastrophic cancellation). Uses
`pnorm(., log.p = TRUE)` and a stable log-subtraction.

## Usage

``` r
log_phi_diff(hi, lo)
```

## Arguments

- hi:

  numeric, upper bound

- lo:

  numeric, lower bound (must be \<= hi elementwise)

## Value

numeric, log(Phi(hi) - Phi(lo))
