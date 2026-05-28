# Get fitted object from MCMC results

Get fitted object from MCMC results

## Usage

``` r
get_fit_object(
  APS,
  BPS,
  UVPS,
  YPS,
  BETA,
  VC,
  GOF,
  Xlist,
  actorByYr,
  colActorByYr = NULL,
  start_vals,
  symmetric,
  tryErrorChecks,
  model.name = NULL,
  U = NULL,
  V = NULL,
  dynamic_uv = FALSE,
  dynamic_ab = FALSE,
  bip = FALSE,
  rho_ab = NULL,
  rho_uv = NULL,
  family = NULL,
  odmax = NULL,
  nA = NULL,
  nB = NULL,
  n_time = NULL,
  Y_obs = NULL,
  G = NULL,
  dynamic_beta = FALSE,
  beta_dynamic_mask = NULL,
  beta_dynamic_groups = NULL,
  rho_beta = NULL,
  sigma_beta = NULL,
  RHO_BETA = NULL,
  SIGMA_BETA = NULL
)
```

## Arguments

- APS:

  summed additive sender random effects (or matrix for dynamic)

- BPS:

  summed additive receiver random effects (or matrix for dynamic)

- UVPS:

  summed multiplicative random effects

- YPS:

  summed Y posterior predictive values

- BETA:

  Matrix of draws for regression coefficient estimates

- VC:

  Matrix of draws for variance estimates

- GOF:

  Matrix of draws for goodness of fit calculations

- Xlist:

  List based version of design array

- actorByYr:

  List of actors by time point. In bipartite mode this is the per-year
  list of row actors.

- colActorByYr:

  Bipartite only. List of column actors by time point; defaults to
  `NULL` (unipartite).

- start_vals:

  start_vals for future model run

- symmetric:

  logical indicating whether model is symmetric

- tryErrorChecks:

  list with counts of MCMC errors

- model.name:

  Name of the model (optional)

- U:

  Latent sender positions (optional, for dynamic UV)

- V:

  Latent receiver positions (optional, for dynamic UV)

- dynamic_uv:

  logical indicating whether UV effects are dynamic

- dynamic_ab:

  logical indicating whether additive effects are dynamic

- bip:

  logical indicating whether the network is bipartite

- rho_ab:

  temporal correlation parameter for additive effects (optional)

- rho_uv:

  temporal correlation parameter for multiplicative effects (optional)

- family:

  character string specifying the model family (e.g., "binary",
  "normal", "poisson")

- odmax:

  vector of maximum ranks for ordinal or fixed rank nomination families

- nA:

  number of actors in first mode (for bipartite networks)

- nB:

  number of actors in second mode (for bipartite networks)

- n_time:

  number of time periods (for longitudinal models)

- Y_obs:

  original observed network (stored for residuals computation)

- G:

  bipartite interaction matrix mapping row to column latent spaces

- dynamic_beta:

  logical or scalar; whether the BETA storage is 3-D (dynamic_beta
  path). Default `FALSE`.

- beta_dynamic_mask:

  logical vector marking which coefficients are dynamic.

- beta_dynamic_groups:

  character vector of per-coefficient block labels ("intercept", "dyad",
  "row", "col"); `""` for static coefficients.

- rho_beta:

  named numeric vector of per-block AR(1) rho values (one per dynamic
  block).

- sigma_beta:

  named numeric vector of per-block AR(1) innovation standard
  deviations.

- RHO_BETA:

  matrix of per-iteration rho_beta draws (rows = MCMC draw, cols =
  dynamic block).

- SIGMA_BETA:

  matrix of per-iteration sigma_beta draws.

## Value

Fitted AME object

## Author

Shahryar Minhas
