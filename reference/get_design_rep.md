# Create design array for replicate data

Create design array for replicate data

## Usage

``` r
get_design_rep(Y, Xdyad, Xrow, Xcol, actorSet, intercept, n, N, pr, pc, pd)
```

## Arguments

- Y:

  dependent variable in array format

- Xdyad:

  dyadic covariates in array format

- Xrow:

  sender covariates in array format

- Xcol:

  receiver covariates in array format

- actorSet:

  vector of actors

- intercept:

  logical indicating whether to include intercept

- n:

  number of actors

- N:

  number of replicates

- pr:

  number of receiver covariates

- pc:

  number of sender covariates

- pd:

  number of dyadic covariates

## Value

returns list of design array values necessary for ame_repL

## Author

Shahryar Minhas
