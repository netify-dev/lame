# AME model fitting routine

An MCMC routine providing a fit to an additive and multiplicative
effects (AME) regression model to cross-sectional relational data of
various types. This function supports both unipartite (square) and
bipartite (rectangular) networks. For longitudinal networks, use the
`lame` function. Original implementation by Peter Hoff.

## Usage

``` r
ame(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(family=="rrl") ,
cvar = TRUE,  dcor = !symmetric, nvar=TRUE, R = 0, R_row = NULL, R_col = NULL,
mode = c("unipartite", "bipartite"), family="normal",
intercept=!is.element(family,c("rrl","ordinal")),
symmetric=FALSE,
odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
prior=list(), g=NA,
seed = 6886, nscan = 10000, burn = 500, odens = 25,
verbose = TRUE, gof=TRUE, custom_gof=NULL,
start_vals=NULL, periodic_save=FALSE, out_file=NULL,
save_interval=0.25, model.name=NULL,
posterior_opts = NULL, n_chains = 1, cores = 1,
use_sparse_matrices = FALSE, print)
```

## Arguments

- Y:

  For unipartite: an n x n square relational matrix. For bipartite: an
  nA x nB rectangular relational matrix where nA is the number of row
  nodes and nB is the number of column nodes. See family below for
  various data types.

- Xdyad:

  For unipartite: an n x n x pd array of covariates. For bipartite: an
  nA x nB x pd array of covariates.

- Xrow:

  For unipartite: an n x pr matrix of nodal row covariates. For
  bipartite: an nA x pr matrix of row node covariates.

- Xcol:

  For unipartite: an n x pc matrix of nodal column covariates. For
  bipartite: an nB x pc matrix of column node covariates.

- rvar:

  logical: fit row random effects (asymmetric case)?

- cvar:

  logical: fit column random effects (asymmetric case)?

- dcor:

  logical: fit a dyadic correlation (asymmetric case)? Note: not used
  for bipartite networks.

- nvar:

  logical: fit nodal random effects (symmetric case)?

- R:

  integer: dimension of the multiplicative effects (can be zero). For
  bipartite networks, this is used as the default for both R_row and
  R_col if they are not specified.

- R_row:

  integer: for bipartite networks, dimension of row node multiplicative
  effects (defaults to R)

- R_col:

  integer: for bipartite networks, dimension of column node
  multiplicative effects (defaults to R)

- mode:

  character: either "unipartite" (default) for square networks or
  "bipartite" for rectangular networks

- family:

  character: one of
  "normal","tobit","binary","ordinal","cbin","frn","rrl","poisson" - see
  the details below

- intercept:

  logical: fit model with an intercept?

- symmetric:

  logical: Is the sociomatrix symmetric by design?

- odmax:

  a scalar integer or vector of length n giving the maximum number of
  nominations that each node may make - used for "frn" and "cbin"
  families

- prior:

  a list containing hyperparameters for the prior distributions.
  Available options and their defaults:

  Sab0

  :   Prior covariance matrix for additive effects (default: diag(2)). A
      2x2 matrix where Sab0\\1,1\\ is the prior variance for row
      effects, Sab0\\2,2\\ is the prior variance for column effects, and
      off-diagonals control correlation between row and column effects.

  eta0

  :   Prior degrees of freedom for covariance of multiplicative effects
      (default: 4 + 3 \\ n/100, where n is the number of actors). Higher
      values impose stronger shrinkage on the latent factors. Common
      values: 4-10 for weak shrinkage, 10-20 for moderate, \>20 for
      strong shrinkage.

  etaab

  :   Prior degrees of freedom for covariance of additive effects
      (default: 4 + 3 \\ n/100). Controls shrinkage of row/column random
      effects. Larger values shrink effects toward zero.

  s20

  :   Prior variance for regression coefficients (default: 1). Larger
      values allow for larger coefficient values.

  s2u0

  :   Prior variance for multiplicative effects (default: 1).

  Suv0

  :   Prior covariance for multiplicative effects (default: identity
      matrix).

  Common usage: prior = list(Sab0 = diag(c(2, 2)), eta0 = 10) for
  moderate shrinkage, or prior = list(Sab0 = diag(c(0.5, 0.5))) for
  tighter control.

- g:

  optional scalar or vector of length dim(X)\\3\\ for g-prior on
  regression coefficients. If not specified, defaults are: for normal
  family, g = n*var(Y); for tobit, g = n*var(Y)\*4; for other families,
  g = n, where n is the number of non-missing dyads. The g-prior
  controls the variance of regression coefficients.

- seed:

  random seed

- nscan:

  number of iterations of the Markov chain (beyond burn-in)

- burn:

  burn in for the Markov chain

- odens:

  output density for the Markov chain

- verbose:

  logical: print progress while running? Default TRUE.

- gof:

  logical: calculate goodness of fit statistics? Setting to TRUE adds
  approximately 2-5% to runtime. For faster sampling without GOF
  overhead, set gof=FALSE and use gof() after model fitting.

- custom_gof:

  optional function or list of named functions for computing custom
  goodness-of-fit statistics. Each function must accept a single matrix
  Y as input and return a numeric vector. If a single function is
  provided, it should return a named vector. If a list of functions is
  provided, each function should return a single value and will be named
  according to the list names. Custom statistics will be computed in
  addition to default statistics. Example: custom_gof = function(Y)
  c(density = mean(Y \> 0, na.rm = TRUE))

- start_vals:

  List from previous model run containing parameter starting values for
  new MCMC

- periodic_save:

  logical: indicating whether to periodically save MCMC results

- out_file:

  character vector indicating name and path in which file should be
  stored if periodic_save is selected. For example, on an Apple OS
  out_file="~/Desktop/ameFit.rda".

- save_interval:

  quantile interval indicating when to save during post burn-in phase.

- model.name:

  optional string for model selection output

- posterior_opts:

  optional list of posterior sampling options

- n_chains:

  integer: number of MCMC chains to run (default: 1)

- cores:

  integer: number of cores for parallel chains (default: 1)

- use_sparse_matrices:

  logical: use sparse matrix storage for large networks? (default:
  FALSE). Recommended only for truly sparse networks (\< 10% non-zero
  entries).

- print:

  Deprecated. Use `verbose` instead.

## Value

**Posterior Samples (full MCMC chains):**

- BETA:

  Regression coefficients (nscan xp matrix)

- VC:

  Variance components (nscan xk matrix)

- GOF:

  Goodness-of-fit statistics (nscan x4 matrix). First row contains
  observed values, remaining rows contain posterior predictive samples.
  See [`gof`](https://netify-dev.github.io/lame/reference/gof.md) for
  post-hoc computation and
  [`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md)
  for visualization.

**Posterior Means (averaged over chain):**

- APM:

  Additive row/sender effects (n-vector)

- BPM:

  Additive column/receiver effects (m-vector); NULL for symmetric
  networks

- U:

  Multiplicative row/sender factors (n xR matrix)

- V:

  Multiplicative column/receiver factors (m xR matrix); NULL for
  symmetric networks

- L:

  Eigenvalue matrix (R xR diagonal); symmetric networks only

- YPM:

  Posterior mean of Y on response scale (for predictions and imputing
  missing values)

**Metadata:**

- family:

  Model family (normal, binary, etc.)

- mode:

  Network mode (unipartite or bipartite)

- symmetric:

  Logical indicating if network is symmetric

- R:

  Dimension of multiplicative effects

**Optional Posterior Samples (if requested via posterior_options):**

- U_samples:

  Samples of U (n xR xiterations array)

- V_samples:

  Samples of V (m xR xiterations array)

- a_samples:

  Samples of row effects (n xiterations matrix)

- b_samples:

  Samples of column effects (m xiterations matrix)

**Note on reconstructing removed matrices:** To save memory, EZ
(expected latent network) and UVPM/ULUPM (multiplicative products) are
not stored but can be reconstructed using:

- `reconstruct_EZ(fit)` - Returns linear predictor (link scale, not
  response scale)

- `reconstruct_UVPM(fit)` - Returns U\\

**Generating posterior distributions:** Use
`simulate_posterior(fit, component="UV")` to generate posterior samples
for components where only means are stored, or use
[`posterior_options()`](https://netify-dev.github.io/lame/reference/posterior_options.md)
during model fitting to save full posterior samples.

- model.name:

  Name of the model (if provided)

## Details

This command provides posterior inference for parameters in AME models
of cross-sectional relational data, assuming one of eight possible data
types/models. The function supports both unipartite networks (square
adjacency matrices) and bipartite networks (rectangular adjacency
matrices with distinct row and column node sets) for single time point
analysis.

**Theoretical Foundation:**

The AME model decomposes network structure into several components:
\$\$y\_{ij} = \beta'x\_{ij} + a_i + b_j + u_i'v_j + \epsilon\_{ij}\$\$
where:

- \\\beta'x\_{ij}\\: Fixed effects of dyadic/nodal covariates

- \\a_i\\: Additive sender (row) effect for node i

- \\b_j\\: Additive receiver (column) effect for node j

- \\u_i'v_j\\: Multiplicative interaction between latent factors

- \\\epsilon\_{ij}\\: Dyadic error term (may be correlated)

This specification generalizes the social relations model (Warner et al.
1979) and latent space models (Hoff et al. 2002) within a unified
framework.

**Prior Distributions:**

The model uses conjugate and semi-conjugate priors where possible:

- Regression coefficients: \\\beta \sim N(0, g\sigma^2(X'X)^{-1})\\
  (g-prior)

- Additive effects: \\(a_i, b_i)' \sim N(0, \Sigma\_{ab})\\ jointly

- Covariance: \\\Sigma\_{ab} \sim IW(\eta_0, \eta_0 S\_{ab0})\\
  (inverse-Wishart)

- Multiplicative effects: Hierarchical shrinkage via \\\eta_0\\

- Dyadic correlation: \\\rho \sim Uniform(-1, 1)\\ with Metropolis
  updates

The inverse-Wishart prior on \\\Sigma\_{ab}\\ allows learning
correlation between sender and receiver effects, capturing reciprocity
patterns.

**Multiplicative Effects (Latent Factors):**

When R \> 0, the model includes R-dimensional latent factors:

- Asymmetric case: \\u_i, v_j \in \mathbb{R}^R\\ with \\u_i'v_j\\
  interaction

- Symmetric case: \\u_i = v_i\\ with eigendecomposition \\ULU'\\

- Captures homophily, transitivity, and community structure

- R chosen via model selection or set to 2-3 for visualization

**Estimation Algorithm:**

The model uses a Gibbs sampler with the following updates:

1.  Sample latent Z given parameters (data augmentation for non-normal
    families)

2.  Update regression coefficients \\\beta\\ via g-prior conjugate
    update

3.  Update additive effects (a,b) jointly with \\\beta\\

4.  Update covariance \\\Sigma\_{ab}\\ from inverse-Wishart

5.  Update multiplicative effects U,V via Gibbs or Metropolis-Hastings

6.  Update dyadic correlation \\\rho\\ via Metropolis-Hastings

7.  Update variance \\\sigma^2\\ (for continuous families)

**Standard Model Types:**

The following data types/models are available:

"normal": A normal AME model (identity link: \\E\[Y\] = \eta\\).

"tobit": A tobit AME model for censored continuous data. Values are
censored at zero, appropriate for non-negative continuous relational
data (identity link with censoring).

"binary": A binary probit AME model (probit link: \\P(Y=1) =
\Phi(\eta)\\).

"ordinal": An ordinal probit AME model (cumulative probit link). An
intercept is not identifiable in this model.

"cbin": An AME model for censored binary data (probit link with
censoring). The value of 'odmax' specifies the maximum number of links
each row may have.

"frn": An AME model for fixed rank nomination networks. A higher value
of the rank indicates a stronger relationship. The value of 'odmax'
specifies the maximum number of links each row may have.

"rrl": An AME model based on the row ranks. This is appropriate if the
relationships across rows are not directly comparable in terms of scale.
An intercept, row random effects and row regression effects are not
estimable for this model.

"poisson": An overdispersed Poisson AME model for count data (log link:
\\E\[Y\] = \exp(\eta)\\). The linear predictor \\\eta\\ represents
\\\log(\lambda)\\ where \\\lambda\\ is the expected count.

## See also

[`lame`](https://netify-dev.github.io/lame/reference/lame.md) for
longitudinal models,
[`gof`](https://netify-dev.github.io/lame/reference/gof.md) for post-hoc
goodness-of-fit computation,
[`gof_plot`](https://netify-dev.github.io/lame/reference/gof_plot.md)
for visualizing GOF results,
[`latent_positions`](https://netify-dev.github.io/lame/reference/latent_positions.md)
for extracting latent positions as a tidy data frame,
[`procrustes_align`](https://netify-dev.github.io/lame/reference/procrustes_align.md)
for Procrustes alignment of latent positions,
[`summary.ame`](https://netify-dev.github.io/lame/reference/summary.ame.md)
for model summaries,
[`coef.ame`](https://netify-dev.github.io/lame/reference/coef.ame.md)
for coefficient extraction

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
data(YX_bin)
fit <- ame(YX_bin$Y, Xdyad = YX_bin$X, burn = 10, nscan = 100, odens = 1,
           family = "binary", verbose = FALSE)
summary(fit)
#> 
#> === AME Model Summary ===
#> 
#> Call:
#> [1] "Y ~ dyad(intercept, rgpa, rsmoke, cgpa, csmoke, igrade, igpa, ismoke) + a[i] + b[j] + rho*e[ji], family = 'binary'"
#> 
#> Regression coefficients:
#> ------------------------
#>                Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept_dyad   -2.543    0.114 -22.296       0   -2.772   -2.319 ***
#> rgpa_dyad         0.178    0.081   2.193   0.028    0.047     0.32   *
#> rsmoke_dyad       0.254    0.112   2.279   0.023    0.054    0.474   *
#> cgpa_dyad         0.188    0.042   4.472       0    0.107    0.273 ***
#> csmoke_dyad       0.167    0.052   3.248   0.001    0.067     0.27  **
#> igrade_dyad        1.13    0.052  21.782       0    1.032    1.229 ***
#> igpa_dyad         0.049    0.038   1.284   0.199   -0.045    0.123    
#> ismoke_dyad       0.038    0.063   0.601   0.548   -0.084     0.15    
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: p-values are approximate (posterior mean / SD); use credible intervals for inference.
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.377    0.088
#> cab    0.045    0.019
#> vb     0.073    0.011
#> rho    0.723    0.085
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
# Note: you should run the Markov chain much longer in practice
# }
 
```
