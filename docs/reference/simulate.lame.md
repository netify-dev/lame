# Simulate longitudinal networks from a fitted LAME model

Generates multiple longitudinal network realizations from the posterior
distribution of a fitted LAME (Longitudinal AME) model. This function
performs posterior predictive simulation for dynamic networks,
propagating both cross-sectional and temporal uncertainty through the
simulated trajectories.

## Usage

``` r
# S3 method for class 'lame'
simulate(
  object,
  nsim = 100,
  seed = NULL,
  newdata = NULL,
  n_time = NULL,
  burn_in = 0,
  thin = 1,
  return_latent = FALSE,
  start_from = "posterior",
  ...
)
```

## Arguments

- object:

  fitted model object of class "lame"

- nsim:

  number of network trajectories to simulate (default: 100)

- seed:

  random seed for reproducibility

- newdata:

  optional list containing new covariate data:

  Xdyad

  :   list of T dyadic covariate arrays (n x n x p or nA x nB x p)

  Xrow

  :   list of T row/sender covariate matrices (n x p or nA x p)

  Xcol

  :   list of T column/receiver covariate matrices (n x p or nB x p)

  If NULL, uses covariates from original model fit

- n_time:

  number of time periods to simulate. If NULL, uses same as original
  data

- burn_in:

  number of initial MCMC samples to discard (default: 0)

- thin:

  thinning interval for MCMC samples (default: 1, use every sample)

- return_latent:

  logical: return latent Z matrices in addition to Y? (default: FALSE)

- start_from:

  character: how to initialize the simulation

  "posterior"

  :   start from posterior mean (default)

  "random"

  :   random initialization

  "data"

  :   use first time point from original data

- ...:

  additional arguments (not currently used)

## Value

A list with components:

- Y:

  list of nsim simulated longitudinal network trajectories, each element
  is a list of T networks

- Z:

  if return_latent=TRUE, list of nsim latent Z trajectories

- family:

  the family of the model (binary, normal, etc.)

- mode:

  network mode (unipartite or bipartite)

- n_time:

  number of time periods

## Details

**Mathematical Framework for Longitudinal Networks:**

The LAME model extends AME to multiple time periods T with potential
temporal dependencies. For each time t = 1, ..., T: \$\$Y\_{ij,t} \sim
F(Z\_{ij,t})\$\$ where the latent network evolves as: \$\$Z\_{ij,t} =
\beta^T x\_{ij,t} + a\_{i,t} + b\_{j,t} + u\_{i,t}^T v\_{j,t} +
\epsilon\_{ij,t}\$\$

**Temporal Dynamics:**

LAME can incorporate three types of temporal dependencies:

1.  **Static Effects:** Parameters constant over time

    - \\a\_{i,t} = a_i\\, \\b\_{j,t} = b_j\\ for all t

    - \\u\_{i,t} = u_i\\, \\v\_{j,t} = v_j\\ for all t

2.  **Dynamic Additive Effects:** AR(1) process for random effects
    \$\$a\_{i,t} = \rho\_{ab} a\_{i,t-1} + \eta\_{i,t}, \quad
    \eta\_{i,t} \sim N(0, \sigma_a^2(1-\rho\_{ab}^2))\$\$ \$\$b\_{j,t} =
    \rho\_{ab} b\_{j,t-1} + \xi\_{j,t}, \quad \xi\_{j,t} \sim N(0,
    \sigma_b^2(1-\rho\_{ab}^2))\$\$ where \\\rho\_{ab}\\ is the temporal
    correlation parameter

3.  **Dynamic Multiplicative Effects:** AR(1) for latent factors
    \$\$u\_{i,t} = \rho\_{uv} u\_{i,t-1} + \omega\_{i,t}\$\$
    \$\$v\_{j,t} = \rho\_{uv} v\_{j,t-1} + \psi\_{j,t}\$\$

**Uncertainty Quantification Process for Trajectories:**

For each simulated trajectory k = 1, ..., nsim:

**Step 1: Parameter Sampling**

- Draw MCMC iteration s uniformly from stored posterior samples

- Extract static parameters: \\\beta^{(s)}\\, variance components

- Extract temporal parameters if applicable: \\\rho\_{ab}^{(s)}\\,
  \\\rho\_{uv}^{(s)}\\

**Step 2: Initialize at t = 1**

Depending on start_from parameter:

- "posterior": Use posterior means as starting values

- "random": Draw from stationary distribution

- For additive effects: \\a\_{i,1}^{(k)} \sim N(0, \sigma_a^2)\\

- For multiplicative effects: Initialize from prior

**Step 3: Evolve Through Time**

For each t = 2, ..., T:

a\) **Update Dynamic Effects** (if applicable): \$\$a\_{i,t}^{(k)} =
\rho\_{ab}^{(s)} a\_{i,t-1}^{(k)} + \eta\_{i,t}^{(k)}\$\$ where
\\\eta\_{i,t}^{(k)} \sim N(0, \sigma_a^2(1-\[\rho\_{ab}^{(s)}\]^2))\\

The innovation variance \\\sigma_a^2(1-\rho\_{ab}^2)\\ ensures
stationarity

b\) **Construct Latent Network**: \$\$E\[Z\_{ij,t}^{(k)}\] =
\beta^{(s)T} x\_{ij,t} + a\_{i,t}^{(k)} + b\_{j,t}^{(k)} + u\_{i,t}^T
v\_{j,t}\$\$

c\) **Add Dyadic Noise**: \$\$Z\_{ij,t}^{(k)} = E\[Z\_{ij,t}^{(k)}\] +
\epsilon\_{ij,t}^{(k)}\$\$ with correlation structure preserved from AME
model

d\) **Generate Observations**: Apply appropriate link function based on
family

**Sources of Uncertainty in Longitudinal Context:**

1.  **Cross-sectional uncertainty** (as in AME):

    - Parameter uncertainty from MCMC

    - Random effect variability

    - Dyadic noise

2.  **Temporal uncertainty**:

    - Uncertainty in temporal correlation parameters \\\rho\_{ab},
      \rho\_{uv}\\

    - Innovation noise in AR(1) processes

    - Propagation of uncertainty through time (compounds over periods)

3.  **Initial condition uncertainty**:

    - Different starting values lead to different trajectories

    - Captured through start_from options

**Interpretation of Multiple Trajectories:**

Each simulated trajectory represents one possible evolution of the
network consistent with the posterior distribution. Variation across
trajectories captures:

- Model parameter uncertainty

- Stochastic variation in temporal evolution

- Accumulated uncertainty over time periods

The ensemble of trajectories provides prediction intervals that widen
over time, reflecting increasing uncertainty in longer-term forecasts.

**Special Considerations:**

1.  **Temporal Correlation:** Higher \\\rho\\ values create smoother
    trajectories with more persistence

2.  **Stationarity:** The AR(1) innovation variance is scaled to
    maintain stationary marginal distributions

3.  **Missing Time Points:** If simulating beyond observed data (n_time
    \> T_observed), covariates are recycled or set to zero with
    appropriate warnings

**Limitations:**

As with simulate.ame, multiplicative effects currently use posterior
means. Full uncertainty would require storing complete MCMC chains for
\\u\_{i,t}, v\_{j,t}\\ at all time points, which is memory-intensive for
large networks and long time series.

## Author

Shahryar Minhas

## Examples

``` r
# \donttest{
# Create simple longitudinal network data
set.seed(1)
n <- 10
nms <- paste0("n", 1:n)
Y_list <- list(
  matrix(rnorm(n * n), n, n, dimnames = list(nms, nms)),
  matrix(rnorm(n * n), n, n, dimnames = list(nms, nms))
)
diag(Y_list[[1]]) <- diag(Y_list[[2]]) <- NA
fit <- lame(Y_list, family = "normal",
            nscan = 50, burn = 10, odens = 1, verbose = FALSE, plot = FALSE)

# Simulate 10 network trajectories from posterior
sims <- simulate(fit, nsim = 10)
# }
```
