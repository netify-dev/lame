# Logit-scale MH update on the AR(1) coefficient rho_G

Proposes `rho* = tanh(atanh(rho) + tau * N(0,1))`; under a uniform
`rho ~ U(-1, 1)` prior the proposal is symmetric in the Fisher-z scale
and the Jacobian is `log(1 - rho^2)`. The MH ratio compares the AR(1)
prior density of `vec(G)_path` at the current vs proposed `rho`.

## Usage

``` r
sample_rho_G_mh(rho_G, vecG_path, sigma_G2, tau = 0.3)
```

## Arguments

- rho_G:

  current AR(1) coefficient

- vecG_path:

  p x N matrix of FFBS-sampled states

- sigma_G2:

  state innovation variance

- tau:

  RW proposal SD on the Fisher-z scale

## Value

list with `rho` (updated), `accept` (logical)
