# Convert an AME / LAME fit to a posterior draws object

Reshapes the BETA + VC posterior draws stored on an `ame` or `lame` fit
into a
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
(or `draws_df`-equivalent) suitable for downstream packages that consume
that representation. Multi-chain fits (`ame(..., n_chains = N)`) are
reshaped with chain as a separate dimension; single-chain fits are
returned as a single-chain `draws_array`.

## Usage

``` r
# S3 method for class 'ame'
as_draws(x, include = c("beta", "vc"), ...)

# S3 method for class 'lame'
as_draws(x, include = c("beta", "vc"), ...)

# S3 method for class 'ame_als'
as_draws(x, ...)
```

## Arguments

- x:

  an `ame` or `lame` fit.

- include:

  character; which parameter blocks to include. Default
  `c("beta", "vc")`; can also include `"rho_uv"`, `"rho_ab"` for dynamic
  LAME fits.

- ...:

  ignored.

## Value

A
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
(when posterior is available) or otherwise a plain 3-D array with named
dimnames.

## Details

The reshape uses `fit$chain_indicator` (set when `n_chains > 1`) to keep
chain identity separate so within-chain Rhat / ESS estimates and
chain-coloured traceplots resolve correctly.
