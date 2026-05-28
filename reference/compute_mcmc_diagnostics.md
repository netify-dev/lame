# Compute MCMC convergence diagnostics for multiple chains

Computes Gelman-Rubin R-hat and effective sample size (ESS) across two
or more independently-seeded chains, for both the regression
coefficients (`BETA`) and the variance components (`VC`). ESS is
computed on the pooled draws, so genuinely non-converged chains report a
low ESS rather than an inflated one.

## Usage

``` r
compute_mcmc_diagnostics(chain_list)
```

## Arguments

- chain_list:

  a list of fitted `ame`/`lame` objects, one per chain (each must carry
  a `BETA` matrix).

## Value

A list with `rhat`, `ess`, `param_names` (over `BETA` and `VC`),
`n_chains` and `n_samples`.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
