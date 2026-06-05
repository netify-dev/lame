# Simulate networks from a fitted AME model

Generates multiple network realizations from a fitted AME model. This
function performs conditional posterior predictive simulation: it draws
from stored MCMC samples when they are available and uses posterior
means for latent components that were not retained.

Unlike `print.ame`, which only displays fitted quantities,
`simulate.ame` draws new networks from the posterior predictive
distribution.

## Usage

``` r
# S3 method for class 'ame'
simulate(
  object,
  nsim = 100,
  seed = NULL,
  newdata = NULL,
  burn_in = 0,
  thin = 1,
  return_latent = FALSE,
  ...
)
```

## Arguments

- object:

  fitted model object of class "ame"

- nsim:

  number of networks to simulate (default: 100)

- seed:

  random seed for reproducibility

- newdata:

  optional list containing new covariate data:

  Xdyad

  :   dyadic covariates (n x n x p array or nA x nB x p for bipartite)

  Xrow

  :   row/sender covariates (n x p matrix or nA x p for bipartite)

  Xcol

  :   column/receiver covariates (n x p matrix or nB x p for bipartite)

  If NULL, uses covariates from original model fit

- burn_in:

  number of initial MCMC samples to discard (default: 0, assumes burn-in
  already removed)

- thin:

  thinning interval for MCMC samples (default: 1, use every sample)

- return_latent:

  logical: return latent Z matrices in addition to Y? (default: FALSE)

- ...:

  additional arguments (not currently used)

## Value

A list with components:

- Y:

  list of nsim simulated networks in the same format as the original
  data

- Z:

  if return_latent=TRUE, list of nsim latent Z matrices

- family:

  the family of the model (binary, normal, etc.)

- mode:

  network mode (unipartite or bipartite)

## Details

**Model:**

The AME model represents networks through a latent variable framework:
\$\$Y\_{ij} \sim F(Z\_{ij})\$\$ where F is the observation model (e.g.,
probit for binary) and Z is the latent network: \$\$Z\_{ij} = \beta^T
x\_{ij} + a_i + b_j + u_i^T v_j + \epsilon\_{ij}\$\$

Components:

- \\\beta\\: regression coefficients for dyadic/nodal covariates

- \\a_i, b_j\\: additive sender and receiver random effects

- \\u_i, v_j\\: multiplicative latent factors (dimension R)

- \\\epsilon\_{ij}\\: dyadic random effects with correlation \\\rho\\

**Simulation:**

For each simulated network k = 1, ..., nsim:

1.  **Parameter draw:** Draw parameter set \\\theta^{(k)}\\ from MCMC
    chains:

    - Sample iteration s uniformly from stored MCMC samples

    - Extract \\\beta^{(s)}\\, variance components \\(v_a^{(s)},
      v_b^{(s)}, v_e^{(s)}, \rho^{(s)})\\

2.  **Random effects:** Sample new random effects from posterior
    distributions:

    - \\a_i^{(k)} \sim N(0, v_a^{(s)})\\ for i = 1, ..., n (row effects)

    - \\b_j^{(k)} \sim N(0, v_b^{(s)})\\ for j = 1, ..., m (column
      effects)

    - Fresh draws from the posterior variance carry random-effect
      uncertainty into the simulated networks

3.  **Latent network:** Build expected latent positions:
    \$\$E\[Z\_{ij}^{(k)}\] = \beta^{(s)T} x\_{ij} + a_i^{(k)} +
    b_j^{(k)} + \hat{u}\_i^T \hat{v}\_j\$\$ where \\\hat{u}\_i,
    \hat{v}\_j\\ are posterior mean latent factors

4.  **Dyadic correlation:** Add correlated noise structure:
    \$\$Z\_{ij}^{(k)} = E\[Z\_{ij}^{(k)}\] + \epsilon\_{ij}^{(k)}\$\$
    where \\\epsilon\\ has covariance structure: \$\$Cov(\epsilon\_{ij},
    \epsilon\_{ji}) = \rho^{(s)} v_e^{(s)}\$\$ \$\$Var(\epsilon\_{ij}) =
    v_e^{(s)}\$\$

5.  **Observation model:** Generate the observed network:

    - Binary: \\Y\_{ij}^{(k)} = I(Z\_{ij}^{(k)} \> 0)\\

    - Normal: \\Y\_{ij}^{(k)} = Z\_{ij}^{(k)}\\

    - Poisson: \\Y\_{ij}^{(k)} \sim Poisson(\exp(Z\_{ij}^{(k)}))\\

    - Other families use appropriate link functions

**Sources of uncertainty:**

The simulation captures three types of uncertainty:

1.  **Parameter uncertainty:** Different MCMC samples yield different
    \\\beta, v_a, v_b, v_e, \rho\\

2.  **Random effect uncertainty:** Fresh draws from \\N(0, v_a), N(0,
    v_b)\\ for each simulation

3.  **Dyadic uncertainty:** Correlated random noise \\\epsilon\_{ij}\\

The resulting simulations propagate uncertainty from the stored
parameter draws and from fresh dyadic/random-effect draws. Latent
components that were not stored as MCMC draws are held at their
posterior means.

**Latent-factor draws:**

Multiplicative effects (U, V) use posterior means unless the fit
retained compatible latent-factor draws. Storing full latent-factor
chains can require substantial additional memory.

**Symmetric Networks:**

For symmetric networks, the model enforces \\a_i = b_i\\ and \\u_i =
v_i\\, and the latent matrix Z is symmetrized before generating
observations.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Fit a model
data(YX_bin)
fit <- ame(YX_bin$Y, Xdyad = YX_bin$X, burn = 10, nscan = 100, odens = 1,
           family = "binary", verbose = FALSE)

# Simulate 10 networks from posterior
sims <- simulate(fit, nsim = 10)
# }
```
