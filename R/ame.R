#' AME model fitting routine
#' 
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to cross-sectional relational data of various types. 
#' This function supports both unipartite (square) and bipartite (rectangular) 
#' networks. For longitudinal networks, use the \code{lame} function. Original
#' implementation by Peter Hoff.
#' 
#' @details
#' This command provides posterior inference for parameters in AME models of
#' cross-sectional relational data, assuming one of eight possible data types/models. 
#' The function supports both unipartite networks (square adjacency matrices) and 
#' bipartite networks (rectangular adjacency matrices with distinct row and column 
#' node sets) for single time point analysis.
#' 
#' \strong{Theoretical Foundation:}
#' 
#' The AME model decomposes network structure into several components:
#' \deqn{y_{ij} = \beta'x_{ij} + a_i + b_j + u_i'v_j + \epsilon_{ij}}
#' where:
#' \itemize{
#'   \item \eqn{\beta'x_{ij}}: Fixed effects of dyadic/nodal covariates
#'   \item \eqn{a_i}: Additive sender (row) effect for node i
#'   \item \eqn{b_j}: Additive receiver (column) effect for node j  
#'   \item \eqn{u_i'v_j}: Multiplicative interaction between latent factors
#'   \item \eqn{\epsilon_{ij}}: Dyadic error term (may be correlated)
#' }
#' 
#' This specification generalizes the social relations model (Warner et al. 1979)
#' and latent space models (Hoff et al. 2002) within a unified framework.
#' 
#' \strong{Prior Distributions:}
#' 
#' The model uses conjugate and semi-conjugate priors where possible:
#' \itemize{
#'   \item Regression coefficients: \eqn{\beta \sim N(0, g\sigma^2(X'X)^{-1})} (g-prior)
#'   \item Additive effects: \eqn{(a_i, b_i)' \sim N(0, \Sigma_{ab})} jointly
#'   \item Covariance: \eqn{\Sigma_{ab} \sim IW(\eta_0, \eta_0 S_{ab0})} (inverse-Wishart)
#'   \item Multiplicative effects: Hierarchical shrinkage via \eqn{\eta_0}
#'   \item Dyadic correlation: \eqn{\rho \sim Uniform(-1, 1)} with Metropolis updates
#' }
#' 
#' The inverse-Wishart prior on \eqn{\Sigma_{ab}} allows learning correlation between
#' sender and receiver effects, capturing reciprocity patterns.
#' 
#' \strong{Multiplicative Effects (Latent Factors):}
#' 
#' When R > 0, the model includes R-dimensional latent factors:
#' \itemize{
#'   \item Asymmetric case: \eqn{u_i, v_j \in \mathbb{R}^R} with \eqn{u_i'v_j} interaction
#'   \item Symmetric case: \eqn{u_i = v_i} with eigendecomposition \eqn{ULU'}
#'   \item Captures homophily, transitivity, and community structure
#'   \item R chosen via model selection or set to 2-3 for visualization
#' }
#' 
#' \strong{Estimation Algorithm:}
#' 
#' The model uses a Gibbs sampler with the following updates:
#' 1. Sample latent Z given parameters (data augmentation for non-normal families)
#' 2. Update regression coefficients \eqn{\beta} via g-prior conjugate update
#' 3. Update additive effects (a,b) jointly with \eqn{\beta}
#' 4. Update covariance \eqn{\Sigma_{ab}} from inverse-Wishart
#' 5. Update multiplicative effects U,V via Gibbs or Metropolis-Hastings
#' 6. Update dyadic correlation \eqn{\rho} via Metropolis-Hastings
#' 7. Update variance \eqn{\sigma^2} (for continuous families)
#' 
#' \strong{Standard Model Types:}
#' 
#' The following data types/models are available:
#' 
#' "normal": A normal AME model (identity link: E[Y] = η).
#' 
#' "tobit": A tobit AME model for censored continuous data. Values are censored
#' at zero, appropriate for non-negative continuous relational data (identity link with censoring).
#' 
#' "binary": A binary probit AME model (probit link: P(Y=1) = Φ(η)).
#' 
#' "ordinal": An ordinal probit AME model (cumulative probit link). An intercept is not 
#' identifiable in this model.
#' 
#' "cbin": An AME model for censored binary data (probit link with censoring). The value of 
#' 'odmax' specifies the maximum number of links each row may have.
#' 
#' "frn": An AME model for fixed rank nomination networks. A higher value of
#' the rank indicates a stronger relationship. The value of 'odmax' specifies
#' the maximum number of links each row may have.
#' 
#' "rrl": An AME model based on the row ranks. This is appropriate if the
#' relationships across rows are not directly comparable in terms of scale. An
#' intercept, row random effects and row regression effects are not estimable
#' for this model.
#' 
#' "poisson": An overdispersed Poisson AME model for count data (log link: E[Y] = exp(η)).
#' The linear predictor η represents log(λ) where λ is the expected count.
#' 
#' @usage ame(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(family=="rrl") ,
#' cvar = TRUE,  dcor = !symmetric, nvar=TRUE, R = 0, R_row = NULL, R_col = NULL,
#' mode = c("unipartite", "bipartite"), family="normal",
#' intercept=!is.element(family,c("rrl","ordinal")),
#' symmetric=FALSE,
#' odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)), 
#' prior=list(), g=NA,
#' seed = 6886, nscan = 10000, burn = 500, odens = 25, 
#' print = TRUE, gof=TRUE, 
#' start_vals=NULL, periodic_save=FALSE, out_file=NULL,
#' save_interval=0.25, model.name=NULL)
#' @param Y For unipartite: an n x n square relational matrix. For bipartite: an nA x nB 
#' rectangular relational matrix where nA is the number of row nodes and nB is the 
#' number of column nodes. See family below for various data types.
#' @param Xdyad For unipartite: an n x n x pd array of covariates. For bipartite: 
#' an nA x nB x pd array of covariates.
#' @param Xrow For unipartite: an n x pr matrix of nodal row covariates. For bipartite: 
#' an nA x pr matrix of row node covariates.
#' @param Xcol For unipartite: an n x pc matrix of nodal column covariates. For bipartite: 
#' an nB x pc matrix of column node covariates.
#' @param rvar logical: fit row random effects (asymmetric case)?
#' @param cvar logical: fit column random effects (asymmetric case)?  
#' @param dcor logical: fit a dyadic correlation (asymmetric case)? Note: not used for bipartite networks.
#' @param nvar logical: fit nodal random effects (symmetric case)?
#' @param R integer: dimension of the multiplicative effects (can be zero). For bipartite networks, 
#' this is used as the default for both R_row and R_col if they are not specified.
#' @param R_row integer: for bipartite networks, dimension of row node multiplicative effects (defaults to R)
#' @param R_col integer: for bipartite networks, dimension of column node multiplicative effects (defaults to R)
#' @param mode character: either "unipartite" (default) for square networks or "bipartite" for rectangular networks
#' @param family character: one of "normal","tobit","binary","ordinal","cbin","frn","rrl","poisson" - see
#' the details below
#' @param intercept logical: fit model with an intercept? 
#' @param symmetric logical: Is the sociomatrix symmetric by design?
#' @param odmax a scalar integer or vector of length n giving the maximum
#' number of nominations that each node may make - used for "frn" and "cbin"
#' families
#' @param prior a list containing hyperparameters for the prior distributions.
#' Available options and their defaults:
#' \describe{
#'   \item{Sab0}{Prior covariance matrix for additive effects (default: diag(2)). 
#'         A 2x2 matrix where Sab0\\[1,1\\] is the prior variance for row effects,
#'         Sab0\\[2,2\\] is the prior variance for column effects, and off-diagonals
#'         control correlation between row and column effects.}
#'   \item{eta0}{Prior degrees of freedom for covariance of multiplicative effects 
#'         (default: 4 + 3 \\* n/100, where n is the number of actors). Higher values 
#'         impose stronger shrinkage on the latent factors. Common values: 4-10 for 
#'         weak shrinkage, 10-20 for moderate, >20 for strong shrinkage.}
#'   \item{etaab}{Prior degrees of freedom for covariance of additive effects 
#'         (default: 4 + 3 \\* n/100). Controls shrinkage of row/column random effects.
#'         Larger values shrink effects toward zero.}
#'   \item{s20}{Prior variance for regression coefficients (default: 1). 
#'         Larger values allow for larger coefficient values.}
#'   \item{s2u0}{Prior variance for multiplicative effects (default: 1).}
#'   \item{Suv0}{Prior covariance for multiplicative effects (default: identity matrix).}
#' }
#' Common usage: prior = list(Sab0 = diag(c(2, 2)), eta0 = 10) for moderate 
#' shrinkage, or prior = list(Sab0 = diag(c(0.5, 0.5))) for tighter control.
#' @param g optional scalar or vector of length dim(X)\\[3\\] for g-prior on 
#' regression coefficients. If not specified, defaults are: for normal family, 
#' g = n*var(Y); for tobit, g = n*var(Y)*4; for other families, g = n, where 
#' n is the number of non-missing dyads. The g-prior controls the variance of 
#' regression coefficients.
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param plot logical: plot results while running?
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics? Setting to TRUE 
#'   adds approximately 2-5% to runtime. For faster sampling without GOF overhead,
#'   set gof=FALSE and use gof() after model fitting.
#' @param custom_gof optional function or list of named functions for computing 
#'   custom goodness-of-fit statistics. Each function must accept a single matrix Y 
#'   as input and return a numeric vector. If a single function is provided, it 
#'   should return a named vector. If a list of functions is provided, each function
#'   should return a single value and will be named according to the list names.
#'   Custom statistics will be computed in addition to default statistics.
#'   Example: custom_gof = function(Y) c(density = mean(Y > 0, na.rm = TRUE))
#' @param start_vals List from previous model run containing parameter starting values for new MCMC
#' @param periodic_save logical: indicating whether to periodically save MCMC results
#' @param out_file character vector indicating name and path in which file should be stored if periodic_save is selected. For example, on an Apple OS out_file="~/Desktop/ameFit.rda".
#' @param save_interval quantile interval indicating when to save during post burn-in phase.
#' @param model.name optional string for model selection output
#' @param use_sparse_matrices logical: use sparse matrix storage for large networks? (default: FALSE).
#'   Recommended only for truly sparse networks (< 10% non-zero entries).
#' @return 
#' \strong{Posterior Samples (full MCMC chains):}
#' \item{BETA}{Regression coefficients (nscan × p matrix)}
#' \item{VC}{Variance components (nscan × k matrix)}
#' \item{GOF}{Goodness-of-fit statistics (nscan × 4 matrix)}
#' 
#' \strong{Posterior Means (averaged over chain):}
#' \item{APM}{Additive row/sender effects (n-vector)}
#' \item{BPM}{Additive column/receiver effects (m-vector); NULL for symmetric networks}
#' \item{U}{Multiplicative row/sender factors (n × R matrix)}
#' \item{V}{Multiplicative column/receiver factors (m × R matrix); NULL for symmetric networks}
#' \item{L}{Eigenvalue matrix (R × R diagonal); symmetric networks only}
#' \item{YPM}{Posterior mean of Y on response scale (for predictions and imputing missing values)}
#' 
#' \strong{Metadata:}
#' \item{family}{Model family (normal, binary, etc.)}
#' \item{mode}{Network mode (unipartite or bipartite)}
#' \item{symmetric}{Logical indicating if network is symmetric}
#' \item{R}{Dimension of multiplicative effects}
#' 
#' \strong{Optional Posterior Samples (if requested via posterior_options):}
#' \item{U_samples}{Samples of U (n × R × iterations array)}
#' \item{V_samples}{Samples of V (m × R × iterations array)}
#' \item{a_samples}{Samples of row effects (n × iterations matrix)}
#' \item{b_samples}{Samples of column effects (m × iterations matrix)}
#' 
#' \strong{Note on reconstructing removed matrices:}
#' To save memory, EZ (expected latent network) and UVPM/ULUPM (multiplicative products) 
#' are not stored but can be reconstructed using:
#' \itemize{
#'   \item \code{reconstruct_EZ(fit)} - Returns linear predictor (link scale, not response scale)
#'   \item \code{reconstruct_UVPM(fit)} - Returns U\%*\%t(V) or U\%*\%L\%*\%t(U)
#' }
#' 
#' \strong{Generating posterior distributions:}
#' Use \code{simulate_posterior(fit, component="UV")} to generate posterior samples
#' for components where only means are stored, or use \code{posterior_options()}
#' during model fitting to save full posterior samples.
#' \item{model.name}{Name of the model (if provided)}
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \dontrun{
#' data(YX_bin) 
#' fit<-ame(YX_bin$Y,YX_bin$X,burn=10,nscan=10,odens=1,family="binary")
#' # Note: you should run the Markov chain much longer in practice
#' } 
#'  
#' @export ame
ame<-function (
  Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
  rvar = !(family=="rrl") , cvar = TRUE, dcor = !symmetric, 
  nvar=TRUE, 
  R = 0, R_row = NULL, R_col = NULL,
  mode = c("unipartite", "bipartite"),
  family="normal",
  intercept=!is.element(family,c("rrl","ordinal")), 
  symmetric=FALSE,
  odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
  prior=list(), g=NA,
  seed = 6886, nscan = 10000, burn = 500, odens = 25,
  print = TRUE, gof=TRUE, custom_gof=NULL,
  start_vals=NULL, periodic_save=FALSE, out_file=NULL,
  save_interval=0.25, model.name=NULL,
  posterior_opts = NULL, n_chains = 1, cores = 1,
  use_sparse_matrices = FALSE
  ){
  
  # Process mode argument
  mode <- match.arg(mode)
  
  # Check if running multiple chains
  if(n_chains > 1) {
    # Use parallel chains function
    return(ame_parallel(Y = Y, n_chains = n_chains, cores = cores,
                       combine_method = "pool",
                       Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
                       rvar = rvar, cvar = cvar, dcor = dcor, nvar = nvar,
                       R = R, R_row = R_row, R_col = R_col,
                       mode = mode, family = family, 
                       intercept = intercept, symmetric = symmetric,
                       odmax = odmax, prior = prior, g = g,
                       seed = seed, nscan = nscan, burn = burn, odens = odens,
                       print = print, gof = gof,
                       start_vals = start_vals, periodic_save = periodic_save,
                       out_file = out_file, save_interval = save_interval,
                       model.name = model.name, posterior_opts = posterior_opts))
  }
  
  # Dispatch to appropriate function based on mode
  if(mode == "unipartite") {
    # Check that matrix is square
    if(nrow(Y) != ncol(Y)) {
      cli::cli_abort(c(
        "Matrix Y must be square for unipartite networks.",
        "i" = "Y has dimensions {nrow(Y)} x {ncol(Y)}.",
        "i" = "Use mode='bipartite' for rectangular matrices."
      ))
    }
    
    # Call unipartite-specific function
    fit <- ame_unipartite(
      Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol,
      rvar=rvar, cvar=cvar, dcor=dcor, nvar=nvar, R=R,
      family=family, intercept=intercept, symmetric=symmetric,
      odmax=odmax, prior=prior, g=g,
      seed=seed, nscan=nscan, burn=burn, odens=odens,
      print=print, gof=gof, custom_gof=custom_gof,
      start_vals=start_vals, periodic_save=periodic_save,
      out_file=out_file, save_interval=save_interval,
      model.name=model.name, posterior_opts=posterior_opts,
      use_sparse_matrices=use_sparse_matrices
    )
    
  } else if(mode == "bipartite") {
    # Check that matrix is rectangular (or at least not constrained to be square)
    if(symmetric) {
      cli::cli_warn("Symmetric option ignored for bipartite networks.")
      symmetric <- FALSE
    }
    
    # Use R_row and R_col if specified, otherwise use R for both
    if(is.null(R_row)) R_row <- R
    if(is.null(R_col)) R_col <- R
    
    # Call bipartite-specific function
    fit <- ame_bipartite(
      Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol,
      rvar=rvar, cvar=cvar, R_row=R_row, R_col=R_col,
      family=family, intercept=intercept,
      odmax=odmax, prior=prior, g=g,
      seed=seed, nscan=nscan, burn=burn, odens=odens,
      print=print, gof=gof, custom_gof=custom_gof,
      start_vals=start_vals, periodic_save=periodic_save,
      out_file=out_file, save_interval=save_interval,
      model.name=model.name, posterior_opts=posterior_opts,
      use_sparse_matrices=use_sparse_matrices
    )
  }
  
  # Return the fitted model
  return(fit)
}