# Get fitted object from MCMC results

Get fitted object from MCMC results

## Usage

``` r
get_start_vals(start_vals, Y, family, xP, rvar, cvar, R, odmax = NULL)
```

## Arguments

- start_vals:

  List object that is null or contains starting values

- Y:

  dependent variable in array format

- family:

  character vector (e.g. 'bin', 'nrm') specifying family type

- xP:

  number of exogenous covariates

- rvar:

  logical indicating whether to include sender random effects

- cvar:

  logical indicating whether to include receiver random effects

- R:

  Number of dimensions for multiplicative effects

- odmax:

  vector of maximum ranks for cbin/frn families (optional)

## Value

List of starting values for MCMC

## Author

Shahryar Minhas
