# AME model fitting routine for longitudinal relational data

An MCMC routine providing a fit to an additive and multiplicative
effects (AME) regression model to longitudinal (time-series) relational
data of various types. Supports both unipartite (square) and bipartite
(rectangular) network structures. For cross-sectional (single time
point) networks, use the `ame` function.

## Usage

``` r
lame(Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(family=="rrl"),
  cvar = TRUE, dcor = !symmetric, nvar=TRUE, R = 0, R_row = NULL, R_col = NULL,
  mode = c("unipartite", "bipartite"),
  dynamic_uv = FALSE, dynamic_ab = FALSE, dynamic_G = FALSE, family ="normal",
intercept=!is.element(family,c("rrl","ordinal")),
symmetric=FALSE,
odmax=NULL, prior=list(), g=NA,
seed = 6886, nscan = 10000, burn = 500, odens = 25, plot=FALSE, print = FALSE, gof=TRUE,
start_vals=NULL, periodic_save=FALSE, out_file=NULL, save_interval=0.25, model.name=NULL)
```

## Arguments

- Y:

  a T length list of n x n relational matrices, where T corresponds to
  the number of replicates (over time, for example). See family below
  for various data types.

- Xdyad:

  a T length list of n x n x pd arrays of covariates

- Xrow:

  a T length list of n x pr matrices of nodal row covariates

- Xcol:

  a T length list of n x pc matrices of nodal column covariates

- rvar:

  logical: fit row random effects (asymmetric case)?

- cvar:

  logical: fit column random effects (asymmetric case)?

- dcor:

  logical: fit a dyadic correlation (asymmetric case)?

- nvar:

  logical: fit nodal random effects (symmetric case)?

- R:

  integer: dimension of the multiplicative effects (can be zero)

- R_row:

  integer: for bipartite networks, dimension of row node multiplicative
  effects (defaults to R)

- R_col:

  integer: for bipartite networks, dimension of column node
  multiplicative effects (defaults to R)

- mode:

  character: either "unipartite" (default) for square networks or
  "bipartite" for rectangular networks

- dynamic_uv:

  logical: fit dynamic multiplicative effects (latent factors) that
  evolve over time using AR(1) processes. When TRUE, the latent
  positions U and V become time-varying, following \\U\_{i,t} =
  \rho\_{uv} U\_{i,t-1} + \epsilon\_{i,t}\\, where epsilon follows N(0,
  sigma_uv^2). This allows actors' positions in latent social space to
  drift smoothly over time, capturing evolving network structure and
  community dynamics. Inspired by dynamic latent space models (Sewell
  and Chen 2015, "Latent Space Models for Dynamic Networks", JASA;
  Durante and Dunson 2014, "Nonparametric Bayes Dynamic Modeling of
  Relational Data", Biometrika). The implementation uses efficient
  blocked Gibbs sampling with C++ acceleration for scalability. Default
  FALSE.

- dynamic_ab:

  logical: fit dynamic additive effects (sender/receiver effects) that
  evolve over time using AR(1) processes. When TRUE, the row effects (a)
  and column effects (b) become time-varying, following \\a\_{i,t} =
  \rho\_{ab} a\_{i,t-1} + \epsilon\_{i,t}\\. This captures temporal
  heterogeneity in actors' baseline propensities to send and receive
  ties, allowing for smooth changes in activity levels and popularity
  over time. For example, an actor's tendency to form outgoing ties
  might gradually increase or decrease across observation periods. The
  AR(1) specification ensures temporal smoothness while allowing for
  actor-specific evolution patterns. Implementation uses conjugate
  updates where possible and C++ for computational efficiency. Default
  FALSE.

- dynamic_G:

  logical: for bipartite networks, fit dynamic interaction matrix G that
  evolves over time. When TRUE, the rectangular interaction matrix G
  becomes time-varying, allowing the mapping between row and column
  latent spaces to change over time. Default FALSE.

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
      values impose stronger shrinkage on the latent factors.

  etaab

  :   Prior degrees of freedom for covariance of additive effects
      (default: 4 + 3 \\ n/100). Controls shrinkage of row/column random
      effects.

  rho_uv_mean

  :   For dynamic_uv=TRUE: Prior mean for UV AR(1) parameter (default:
      0.9). Values close to 1 indicate high temporal persistence.

  rho_uv_sd

  :   For dynamic_uv=TRUE: Prior SD for UV AR(1) parameter (default:
      0.1). Controls uncertainty about temporal dependence.

  sigma_uv_shape

  :   For dynamic_uv=TRUE: Shape parameter for inverse-gamma prior on UV
      innovation variance (default: 2).

  sigma_uv_scale

  :   For dynamic_uv=TRUE: Scale parameter for inverse-gamma prior on UV
      innovation variance (default: 1).

  rho_ab_mean

  :   For dynamic_ab=TRUE: Prior mean for additive effects AR(1)
      parameter (default: 0.8). Controls temporal smoothness of
      sender/receiver effects.

  rho_ab_sd

  :   For dynamic_ab=TRUE: Prior SD for additive effects AR(1) parameter
      (default: 0.15).

  sigma_ab_shape

  :   For dynamic_ab=TRUE: Shape parameter for inverse-gamma prior on
      additive effects innovation variance (default: 2).

  sigma_ab_scale

  :   For dynamic_ab=TRUE: Scale parameter for inverse-gamma prior on
      additive effects innovation variance (default: 1).

  Common usage: prior = list(Sab0 = diag(c(1, 1)), eta0 = 10) for
  stronger shrinkage, or prior = list(rho_uv_mean = 0.95) for higher
  temporal persistence.

- g:

  optional scalar or vector for g-prior on regression coefficients.
  Default is p^2 where p is the number of regression parameters. The
  g-prior controls the variance of regression coefficients: larger
  values allow for larger coefficient values. Can be a vector of length
  p for parameter-specific control.

- seed:

  random seed

- nscan:

  number of iterations of the Markov chain (beyond burn-in)

- burn:

  burn in for the Markov chain

- odens:

  output density for the Markov chain

- plot:

  logical: plot results while running?

- print:

  logical: print results while running?

- gof:

  logical: calculate goodness of fit statistics?

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

## Value

- BETA:

  posterior samples of regression coefficients

- VC:

  posterior samples of the variance parameters

- APM:

  posterior mean of additive row effects a

- BPM:

  posterior mean of additive column effects b

- U:

  posterior mean of multiplicative row effects u. For dynamic_uv=TRUE,
  this is a 3D array (n x R x T)

- V:

  posterior mean of multiplicative column effects v (asymmetric case).
  For dynamic_uv=TRUE, this is a 3D array (n x R x T)

- UVPM:

  posterior mean of UV

- ULUPM:

  posterior mean of ULU (symmetric case)

- L:

  posterior mean of L (symmetric case)

- EZ:

  estimate of expectation of Z matrix

- YPM:

  posterior mean of Y (for imputing missing values)

- GOF:

  observed (first row) and posterior predictive (remaining rows) values
  of four goodness-of-fit statistics

- start_vals:

  Final parameter values from MCMC, can be used as the input for a
  future model run.

- model.name:

  Name of the model (if provided)

## Details

This command provides posterior inference for parameters in AME models
of longitudinal relational data, assuming one of eight possible data
types/models. The model supports both unipartite networks (square
adjacency matrices) and bipartite networks (rectangular adjacency
matrices with distinct row and column node sets) across multiple time
points.

**Dynamic Effects Implementation:**

The dynamic_uv and dynamic_ab parameters enable time-varying latent
representations through autoregressive processes. These extensions are
particularly useful for understanding how network structure evolves over
time.

*Dynamic Multiplicative Effects (dynamic_uv=TRUE):* The latent factors U
and V evolve according to AR(1) processes: \$\$U\_{i,k,t} = \rho\_{uv}
U\_{i,k,t-1} + \epsilon\_{i,k,t}\$\$ where \\\epsilon\_{i,k,t} \sim N(0,
\sigma\_{uv}^2)\\, i indexes actors, k indexes latent dimensions, and t
indexes time. The parameter \\\rho\_{uv}\\ controls temporal persistence
(values near 1 indicate slow evolution). This captures time-varying
homophily, latent community structure, and transitivity dynamics.

Key references:

- Sewell & Chen (2015): Introduced dynamic latent space models with
  actor-specific evolution rates

- Durante & Dunson (2014): Nonparametric Bayesian approach allowing
  flexible evolution of network structure

- Hoff (2011): Hierarchical multilinear models providing theoretical
  foundation for temporal dependencies

*Dynamic Additive Effects (dynamic_ab=TRUE):* The sender (a) and
receiver (b) effects evolve as: \$\$a\_{i,t} = \rho\_{ab} a\_{i,t-1} +
\epsilon\_{i,t}\$\$ \$\$b\_{i,t} = \rho\_{ab} b\_{i,t-1} +
\eta\_{i,t}\$\$ where \\\epsilon\_{i,t}, \eta\_{i,t} \sim N(0,
\sigma\_{ab}^2)\\. This models time-varying individual activity levels
(outdegree) and popularity (indegree).

Applications include:

- Tracking changes in node centrality over time

- Identifying emerging influential actors

- Detecting declining activity patterns

- Modeling life-cycle effects in social networks

*Prior Specification for Dynamic Parameters:*

- \\\rho\_{uv}, \rho\_{ab} \sim TruncNormal(mean, sd, 0, 1)\\: Ensures
  stationarity of AR(1) process

- \\\sigma\_{uv}^2, \sigma\_{ab}^2 \sim InverseGamma(shape, scale)\\:
  Controls innovation variance

- Default priors (\\\rho\_{uv}\\ mean=0.9, \\\rho\_{ab}\\ mean=0.8)
  favor smooth evolution

- Adjust rho\_\*\_mean closer to 1 for slower evolution, closer to 0 for
  more rapid changes

*Computational Considerations:*

- Dynamic effects increase computation by ~30-50\\

- Memory usage scales as O(n*R*T) for dynamic_uv, O(n\*T) for dynamic_ab

- C++ implementation provides ~70\\

- Convergence diagnostics: Monitor rho and sigma parameters carefully

- Effective sample sizes typically lower due to temporal correlation

- Recommend burn \>= 1000 and nscan \>= 20000 for dynamic models

*Model Selection Guidelines:* Use both dynamic_uv and dynamic_ab when:

- Networks show clear temporal trends in density or clustering

- Individual node behavior changes systematically over time

- Community structure evolves (merging, splitting, drift)

Use only dynamic_uv when:

- Latent structure/communities change but individual effects are stable

- Focus is on evolving homophily or clustering patterns

- Network shows structural reconfiguration over time

Use only dynamic_ab when:

- Individual heterogeneity varies but overall structure is stable

- Actors' activity/popularity changes over observation period

- Focus is on individual-level temporal dynamics

**Bipartite Network Models:**

When mode="bipartite", the model handles rectangular adjacency matrices
Y with dimensions n_A x n_B, where n_A and n_B represent the number of
row and column nodes respectively.

*Static Bipartite Case:* The model uses separate latent factor matrices:

- U: n_A x R_row matrix of row node latent positions

- V: n_B x R_col matrix of column node latent positions

- G: R_row x R_col interaction matrix mapping between latent spaces

- Multiplicative term: U G V' captures bipartite community structure

*Dynamic Bipartite Case:* When dynamic_uv=TRUE for bipartite networks:
\$\$U\_{i,k,t} = \rho\_{uv} U\_{i,k,t-1} + \epsilon\_{i,k,t}\$\$
\$\$V\_{j,k,t} = \rho\_{uv} V\_{j,k,t-1} + \eta\_{j,k,t}\$\$ where i
indexes row nodes, j indexes column nodes, k indexes latent dimensions.

When dynamic_G=TRUE, the interaction matrix also evolves: \$\$G\_{k,l,t}
= \rho_G G\_{k,l,t-1} + \xi\_{k,l,t}\$\$ allowing the mapping between
row and column latent spaces to change over time.

*Key Differences from Unipartite Models:*

- No dyadic correlation (rho): Bipartite edges are inherently directed

- Separate dimensions: R_row and R_col can differ for row/column spaces

- Rectangular structure: Network density patterns differ from square
  matrices

- Community interpretation: Captures affiliation patterns between node
  types

**Standard AME Model Types:**

The following describes the eight standard data types/models available:

"normal": A normal AME model.

"tobit": A tobit AME model for censored continuous data. Values are
censored at zero, appropriate for non-negative continuous relational
data.

"binary": A binary probit AME model.

"ordinal": An ordinal probit AME model. An intercept is not identifiable
in this model.

"cbin": An AME model for censored binary data. The value of 'odmax'
specifies the maximum number of links each row may have.

"frn": An AME model for fixed rank nomination networks. A higher value
of the rank indicates a stronger relationship. The value of 'odmax'
specifies the maximum number of links each row may have.

"rrl": An AME model based on the row ranks. This is appropriate if the
relationships across rows are not directly comparable in terms of scale.
An intercept, row random effects and row regression effects are not
estimable for this model.

"poisson": An overdispersed Poisson AME model for count data. The latent
variable represents the log mean of the Poisson distribution.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
data(YX_bin_list) 
fit<-lame(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,family="binary")
# you should run the Markov chain much longer than this
```
