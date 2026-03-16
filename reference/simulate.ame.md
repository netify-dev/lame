# Simulate networks from a fitted AME model

Generates multiple network realizations from the posterior distribution
of a fitted AME model. This function performs posterior predictive
simulation by drawing from the full joint posterior distribution of
model parameters, thereby propagating parameter uncertainty into the
simulated networks.

**Key Difference from print():** While `print.ame` displays existing
model results without computation, `simulate.ame` actively generates new
network data by sampling from the posterior predictive distribution.
This is computationally intensive and produces new datasets for
analysis.

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

**Mathematical Framework:**

The AME model represents networks through a latent variable framework:
\$\$Y\_{ij} \sim F(Z\_{ij})\$\$ where F is the observation model (e.g.,
probit for binary) and Z is the latent network: \$\$Z\_{ij} = \beta^T
x\_{ij} + a_i + b_j + u_i^T v_j + \epsilon\_{ij}\$\$

Components:

- \\\beta\\: regression coefficients for dyadic/nodal covariates

- \\a_i, b_j\\: additive sender and receiver random effects

- \\u_i, v_j\\: multiplicative latent factors (dimension R)

- \\\epsilon\_{ij}\\: dyadic random effects with correlation \\\rho\\

**Uncertainty Quantification Process:**

For each simulated network k = 1, ..., nsim:

1.  **Parameter Sampling:** Draw parameter set \\\theta^{(k)}\\ from
    MCMC chains:

    - Sample iteration s uniformly from stored MCMC samples

    - Extract \\\beta^{(s)}\\, variance components \\(v_a^{(s)},
      v_b^{(s)}, v_e^{(s)}, \rho^{(s)})\\

2.  **Random Effects Generation:** Sample new random effects from
    posterior distributions:

    - \\a_i^{(k)} \sim N(0, v_a^{(s)})\\ for i = 1, ..., n (row effects)

    - \\b_j^{(k)} \sim N(0, v_b^{(s)})\\ for j = 1, ..., m (column
      effects)

    - Note: We sample fresh from the posterior variance rather than
      using point estimates to properly propagate uncertainty

3.  **Latent Network Construction:** Build expected latent positions:
    \$\$E\[Z\_{ij}^{(k)}\] = \beta^{(s)T} x\_{ij} + a_i^{(k)} +
    b_j^{(k)} + \hat{u}\_i^T \hat{v}\_j\$\$ where \\\hat{u}\_i,
    \hat{v}\_j\\ are posterior mean latent factors

4.  **Dyadic Correlation:** Add correlated noise structure:
    \$\$Z\_{ij}^{(k)} = E\[Z\_{ij}^{(k)}\] + \epsilon\_{ij}^{(k)}\$\$
    where \\\epsilon\\ has covariance structure: \$\$Cov(\epsilon\_{ij},
    \epsilon\_{ji}) = \rho^{(s)} v_e^{(s)}\$\$ \$\$Var(\epsilon\_{ij}) =
    v_e^{(s)}\$\$

5.  **Observation Model:** Generate observed network based on family:

    - Binary: \\Y\_{ij}^{(k)} = I(Z\_{ij}^{(k)} \> 0)\\

    - Normal: \\Y\_{ij}^{(k)} = Z\_{ij}^{(k)}\\

    - Poisson: \\Y\_{ij}^{(k)} \sim Poisson(\exp(Z\_{ij}^{(k)}))\\

    - Other families use appropriate link functions

**Sources of Uncertainty:**

The simulation captures three types of uncertainty:

1.  **Parameter uncertainty:** Different MCMC samples yield different
    \\\beta, v_a, v_b, v_e, \rho\\

2.  **Random effect uncertainty:** Fresh draws from \\N(0, v_a), N(0,
    v_b)\\ for each simulation

3.  **Dyadic uncertainty:** Correlated random noise \\\epsilon\_{ij}\\

This approach provides proper posterior predictive distributions that
account for all sources of uncertainty in the model. The variation
across simulated networks reflects our posterior uncertainty about the
data generating process.

**Limitations:**

Currently, multiplicative effects (U, V) use posterior means rather than
sampling from their full posterior. For complete uncertainty
quantification, one would need to store and sample from the full MCMC
chains of these latent factors, which would require substantial
additional memory.

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
