# Options for saving posterior samples during MCMC

Options for saving posterior samples during MCMC

## Usage

``` r
posterior_options(save_UV = FALSE, save_ab = FALSE, thin_UV = 10, thin_ab = 10)
```

## Arguments

- save_UV:

  Logical; whether to save samples of U and V matrices (default FALSE)

- save_ab:

  Logical; whether to save samples of additive effects a and b (default
  FALSE)

- thin_UV:

  Integer; thinning interval for U/V samples (default 10)

- thin_ab:

  Integer; thinning interval for a/b samples (default 10)

## Value

List of posterior saving options

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
