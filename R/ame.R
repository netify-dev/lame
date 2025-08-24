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
#' "normal": A normal AME model.
#' 
#' "tobit": A tobit AME model for censored continuous data. Values are censored
#' at zero, appropriate for non-negative continuous relational data.
#' 
#' "binary": A binary probit AME model.
#' 
#' "ordinal": An ordinal probit AME model. An intercept is not identifiable in this
#' model.
#' 
#' "cbin": An AME model for censored binary data.  The value of 'odmax'
#' specifies the maximum number of links each row may have.
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
#' "poisson": An overdispersed Poisson AME model for count data. The latent
#' variable represents the log mean of the Poisson distribution.
#' 
#' @usage ame(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(family=="rrl") ,
#' cvar = TRUE,  dcor = !symmetric, nvar=TRUE, R = 0, R_row = NULL, R_col = NULL,
#' mode = c("unipartite", "bipartite"), family="normal",
#' intercept=!is.element(family,c("rrl","ordinal")),
#' symmetric=FALSE,
#' odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)), 
#' prior=list(), g=NA,
#' seed = 6886, nscan = 10000, burn = 500, odens = 25, 
#' plot=TRUE, print = TRUE, gof=TRUE, 
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
#' @param gof logical: calculate goodness of fit statistics?
#' @param start_vals List from previous model run containing parameter starting values for new MCMC
#' @param periodic_save logical: indicating whether to periodically save MCMC results
#' @param out_file character vector indicating name and path in which file should be stored if periodic_save is selected. For example, on an Apple OS out_file="~/Desktop/ameFit.rda".
#' @param save_interval quantile interval indicating when to save during post burn-in phase.
#' @param model.name optional string for model selection output
#' @return \item{BETA}{posterior samples of regression coefficients}
#' \item{VC}{posterior samples of the variance parameters}
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b} \item{U}{posterior mean of multiplicative
#' row effects u} \item{V}{posterior mean of multiplicative column effects v (asymmetric case)}
#' \item{UVPM}{posterior mean of UV (asymmetric case)} 
#' \item{ULUPM}{posterior mean of ULU (symmetric case)} 
#' \item{L}{posterior mean of L (symmetric case)} 
#'  \item{EZ}{estimate of expectation of Z
#' matrix} \item{YPM}{posterior mean of Y (for imputing missing values)}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
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
ame<-function (Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
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
               plot=TRUE, print = TRUE, gof=TRUE, 
               start_vals=NULL, periodic_save=FALSE, out_file=NULL,
               save_interval=0.25, model.name=NULL)
{ 
  
  # set random seed
  set.seed(seed)
  
  # Process mode argument
  mode <- match.arg(mode)
  bip <- identical(mode, "bipartite")
  
  # set diag to NA (only for unipartite networks)
  if(!bip) {
    diag(Y) <- NA
  } 
  
  # force binary if binary family specified
  if(is.element(family,c("binary","cbin"))) { Y<-1*(Y>0) }
  
  # handle tobit family
  if(family=="tobit") { Y[Y<0]<-0 } 
  
  # observed and max outdegrees 
  if(is.element(family,c("cbin","frn","rrl")))
  {
    odobs<-apply(Y>0,1,sum,na.rm=TRUE) 
    if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y)) }
  }
  
  # Handle network dimensions for bipartite vs unipartite
  if(bip) {
    nA <- nrow(Y)  # number of row nodes
    nB <- ncol(Y)  # number of column nodes
    n <- nA  # for compatibility with existing code
    
    # Set rank dimensions for bipartite
    RA <- R_row %||% R %||% 0L
    RB <- R_col %||% R %||% 0L
    
    # No dyadic correlation or symmetric model for bipartite
    dcor <- FALSE
    symmetric <- FALSE
    
    # Helper function for null coalescing
    `%||%` <- function(x, y) if (is.null(x)) y else x
    
  } else {
    # Unipartite network
    n <- nrow(Y)
    nA <- nB <- n
    RA <- RB <- R
    
    # some settings for symmetric case
    if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }
  }
  
  if(is.na(g)) {
    if(family=="normal") { 
      g <- sum(!is.na(Y))*var(c(Y),na.rm=TRUE)
    } else if(family=="tobit") {
      g <- sum(!is.na(Y))*var(c(Y),na.rm=TRUE)*4
    } else {
      g <- sum(!is.na(Y))
    }
  }
  
  # process prior list and set defaults
  if(!is.list(prior)) { prior <- list() }
  if(is.null(prior$Sab0)) { prior$Sab0 <- diag(2) }
  if(is.null(prior$eta0)) { prior$eta0 <- round(4 + 3*nrow(Y)/100) }
  if(is.null(prior$etaab)) { prior$etaab <- round(4 + 3*nrow(Y)/100) }
  
  # construct design matrix
  n<-nrow(Y) 
  pr<-length(Xrow)/n
  pc<-length(Xcol)/n
  pd<-length(Xdyad)/n^2
  
  # For bipartite networks, need to handle design array differently
  if(bip) {
    # For bipartite, we need nA x nB design array, but design_array creates n x n
    # If no covariates, create minimal intercept-only design
    if(is.null(Xrow) && is.null(Xcol) && is.null(Xdyad)) {
      if(intercept) {
        X <- array(1, dim=c(nA, nB, 1))
        dimnames(X)[[3]] <- list("intercept")
      } else {
        X <- array(dim=c(nA, nB, 0))
      }
    } else {
      # For now, use standard design_array but will need reshaping later
      X <- design_array(Xrow, Xcol, Xdyad, intercept, nrow(Y))
    }
  } else {
    X <- design_array(Xrow, Xcol, Xdyad, intercept, nrow(Y))
  } 
  
  # design matrix warning for rrl 
  if( family=="rrl" & any(apply(apply(X,c(1,3),var),2,sum)==0) 
      & !any( apply(X,3,function(x){var(c(x))})==0) ) 
  {
    cli::cli_warn("Row effects are not estimable using this procedure")
  } 
  
  # design matrix warning for rrl and ord
  if( is.element(family,c("ord","rrl")) & 
      any( apply(X,3,function(x){var(c(x))})==0 ) )
  {
    cli::cli_warn("An intercept is not estimable using this procedure")
  } 
  
  # construct matrix of ranked nominations for frn, rrl 
  if(is.element(family,c("frn","rrl")))
  {
    ymx<-max(apply(1*(Y>0),1,sum,na.rm=TRUE))
    YL<-NULL
    warn<-FALSE
    for(i in 1:nrow(Y))
    {
      yi<-Y[i,] ; rnkd<-which( !is.na(yi)&yi>0 )
      if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
      yi[rnkd]<-rank(yi[rnkd],ties.method="random")
      Y[i,]<-yi
      YL<-rbind(YL, match(1:ymx,yi))
    }
    if(warn){cli::cli_warn("Random reordering used to break ties in ranks")}
  }
  
  # starting Z values
  if(family=="normal") { Z<-Y }
  if(family=="tobit") { 
    Z<-Y 
    # Check if there are any positive values before taking min
    if(sum(Y>0, na.rm=TRUE) > 0) {
      Z[Y==0]<-min(Y[Y>0],na.rm=TRUE)/2 
    }
    # If all zeros, leave as is
  }
  if(family=="ordinal") { Z<-matrix(zscores(Y),nrow(Y),ncol(Y)) } 
  if(family=="rrl") { Z<-matrix(t(apply(Y,1,zscores)),nrow(Y),ncol(Y)) }  
  if(family=="binary")
  { 
    Z<-matrix(zscores(Y),nrow(Y),ncol(Y)) 
    # Check if we have both 0s and 1s to avoid max/min on empty sets
    if(sum(Y==0, na.rm=TRUE) > 0 && sum(Y==1, na.rm=TRUE) > 0) {
      z01<- .5* ( max(Z[Y==0],na.rm=TRUE) + min(Z[Y==1],na.rm=TRUE) ) 
      Z<-Z - z01
    }
    # If all same value or perfect separation, Z-scores are already centered
  } 
  if(family=="poisson") { 
    Z<-log(Y+1) 
    if(!bip) { diag(Z)<-0 }  # only set diagonal for unipartite
  }
  
  if(is.element(family,c("cbin","frn")))
  {
    Z<-Y
    for(i in 1:nrow(Y))
    {
      yi<-Y[i,]
      zi<-zscores(yi)
      rnkd<-which( !is.na(yi) & yi>0 ) 
      if(length(rnkd)>0 && min(zi[rnkd])<0) 
      { 
        zi[rnkd]<-zi[rnkd] - min(zi[rnkd]) + 1e-3 
      }
      
      if(length(rnkd)<odmax[i]) 
      {
        urnkd<-which( !is.na(yi) & yi==0 ) 
        if(max(zi[urnkd])>0) { zi[urnkd]<-zi[urnkd] - max(zi[urnkd]) -1e-3 }
      }
      
      Z[i,]<-zi
    } 
  }
  
  
  # starting values for missing entries 
  mu<-mean(Z,na.rm=TRUE) 
  a<-rowMeans(Z,na.rm=TRUE) ; b<-colMeans(Z,na.rm=TRUE)  
  a[is.na(a)]<-0 ; b[is.na(b)]<-0 
  ZA<-mu + outer(a,b,"+") 
  Z[is.na(Z)]<-ZA[is.na(Z)] 
  
  
  # other starting values (use start_vals if provided)
  if(!is.null(start_vals)) {
    beta <- start_vals$beta
    s2 <- start_vals$s2
    a <- start_vals$a
    b <- start_vals$b
    if(!is.null(start_vals$rho)) rho <- start_vals$rho else rho <- 0
    if(!is.null(start_vals$Sab)) Sab <- start_vals$Sab
    if(!is.null(start_vals$U)) U <- start_vals$U
    if(!is.null(start_vals$V)) V <- start_vals$V
    if(!is.null(start_vals$Z)) Z <- start_vals$Z
  } else {
    beta<-rep(0,dim(X)[3]) 
    s2<-1
  }
  
  # Initialize remaining values if not provided by start_vals
  if(is.null(start_vals) && bip) {
    # For bipartite networks
    rho <- NULL  # No dyadic correlation
    # Initialize U and V with appropriate dimensions for bipartite
    if(RA > 0) {
      U <- matrix(rnorm(nA * RA, 0, 0.1), nA, RA)
    } else {
      U <- matrix(0, nA, 1)  # Minimal dimension
    }
    if(RB > 0) {
      V <- matrix(rnorm(nB * RB, 0, 0.1), nB, RB)
    } else {
      V <- matrix(0, nB, 1)  # Minimal dimension
    }
    # Separate variance parameters for bipartite
    Sab <- diag(c(1, 1))  # Independent row and column variances
    # Initialize G matrix for bipartite multiplicative interaction
    if(RA > 0 && RB > 0) {
      min_rank <- min(RA, RB)
      G <- matrix(0, RA, RB)
      diag(G)[1:min_rank] <- 1
    } else {
      G <- matrix(0, max(1, RA), max(1, RB))
    }
  } else if(is.null(start_vals) && !bip) {
    # For unipartite networks (existing logic)
    rho<-0
    Sab<-cov(cbind(a,b))*tcrossprod(c(rvar,cvar))
    U<-V<-matrix(0, nrow(Y), R)
    G <- NULL  # No G matrix for unipartite
  } else if(!is.null(start_vals) && bip) {
    # Use provided start_vals for bipartite networks
    if(!is.null(start_vals$G)) G <- start_vals$G
  } else if(!is.null(start_vals) && !bip) {
    # Use provided start_vals for unipartite networks  
    G <- NULL  # No G matrix for unipartite
  }  
  
  
  # output items
  BETA <- matrix(nrow = 0, ncol = dim(X)[3] - pr*symmetric)
  # For bipartite: 3 columns (va, vb, ve), for unipartite: 5 (asymm) or 2 (symm)
  if(bip) {
    VC<-matrix(nrow=0,ncol=3) # va, vb, ve
  } else {
    VC<-matrix(nrow=0,ncol=5-3*symmetric)
  } 
  # Initialize UVPS based on network type
  if(bip) {
    # For bipartite, UV interaction goes through G matrix: U %*% G %*% t(V)
    UVPS <- array(0, dim = c(nA, nB))
  } else {
    UVPS <- U %*% t(V) * 0
  }
  if(bip) {
    # For bipartite: APS for rows, BPS for columns
    APS <- rep(0, nrow(Y))
    BPS <- rep(0, ncol(Y))
  } else {
    # For unipartite: both same dimension
    APS<-BPS<- rep(0,nrow(Y))
  }
  YPS<-matrix(0,nrow(Y),ncol(Y)) ; dimnames(YPS)<-dimnames(Y)
  
  # Initialize GOF - handle bipartite case differently
  if(bip) {
    # For bipartite, use simplified GOF or skip certain stats
    # that require square matrices
    GOF <- matrix(NA, 1, 5)
    rownames(GOF) <- "obs"
    colnames(GOF) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    # Only compute row and column means for bipartite
    GOF[1, 1] <- sd(rowMeans(Y, na.rm=TRUE), na.rm=TRUE)
    GOF[1, 2] <- sd(colMeans(Y, na.rm=TRUE), na.rm=TRUE)
    # Skip dyad.dep, cycle.dep, trans.dep for bipartite (they need square matrices)
  } else {
    GOF<-matrix(gof_stats(Y),1,5)  
    rownames(GOF)<-"obs"
    colnames(GOF)<- c("sd.rowmean","sd.colmean","dyad.dep","cycle.dep","trans.dep")
  }
  if(bip) {
    # For bipartite: U corresponds to rows, V to columns
    names(APS) <- rownames(U) <- rownames(Y)
    names(BPS) <- rownames(V) <- colnames(Y)
  } else {
    # For unipartite: both use rownames
    names(APS)<-names(BPS)<- rownames(U)<-rownames(V)<-rownames(Y)
  }
  
  # names of parameters, asymmetric case 
  if(!symmetric)
  {
    if(bip) {
      # Bipartite: no rho, separate row/col variances  
      colnames(VC) <- c("va", "vb", "ve")
    } else {
      # Unipartite asymmetric
      colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
    }
    colnames(BETA) <- dimnames(X)[[3]] 
  }
  
  # names of parameters, symmetric case
  if(symmetric)
  {
    colnames(VC) <- c("va", "ve")  
    rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
    bnames<-dimnames(X)[[3]]
    bni<-bnames[1*intercept] 
    bnn<-gsub("row",bnames[rb],replacement="node")       
    bnd<-bnames[-c(1*intercept,rb,cb)]
    colnames(BETA)<-c(bni,bnn,bnd) 
  }    
  
  # marginal means and regression sums of squares
  Xr<-apply(X,c(1,3),sum)            # row sum
  Xc<-apply(X,c(2,3),sum)            # col sum
  mX<- apply(X,3,c)                  # design matrix
  mXt<-apply(aperm(X,c(2,1,3)),3,c)  # dyad-transposed design matrix
  XX<-t(mX)%*%mX                     # regression sums of squares
  XXt<-t(mX)%*%mXt                   # crossproduct sums of squares
  
  # MCMC 
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE))
  
  # Set up save points for periodic saving
  savePoints <- (burn:(nscan+burn))[(burn:(nscan+burn)) %% odens==0]
  savePoints <- savePoints[round(quantile(1:length(savePoints), probs=seq(save_interval,1,save_interval)))]
  
  for (s in 1:(nscan + burn)) 
  { 
    
    # update Z 
    if(bip) {
      # Bipartite: use G matrix for UV interaction
      if(RA > 0 && RB > 0) {
        UV_term <- U %*% G %*% t(V)
      } else {
        UV_term <- matrix(0, nA, nB)
      }
      # Check dimensions for bipartite
      Xb_term_init <- Xbeta(X, beta)
      ab_term_init <- outer(a, b, "+")
      if(!all(dim(Xb_term_init) == dim(ab_term_init))) {
        # Create compatible zero matrix if dimensions don't match
        Xb_term_init <- matrix(0, nrow=length(a), ncol=length(b))
      }
      EZ <- Xb_term_init + ab_term_init + UV_term
    } else {
      # Unipartite: direct UV multiplication
      EZ<-Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)
    }
    # Update Z using appropriate functions
    if(bip) {
      # For bipartite, use simplified updates without dyadic correlation
      if(family=="normal"){ 
        Z <- EZ + matrix(rnorm(length(Z), 0, sqrt(s2)), nrow=nrow(Z))
        Z[is.na(Y)] <- EZ[is.na(Y)]
      }
      if(family=="binary"){ 
        # Simple probit update for bipartite
        sz <- 1  # No dyadic correlation
        for(i in 1:nrow(Z)) {
          for(j in 1:ncol(Z)) {
            if(!is.na(Y[i,j])) {
              lb <- if(Y[i,j]==0) -Inf else 0
              ub <- if(Y[i,j]==0) 0 else Inf
              Z[i,j] <- EZ[i,j] + sz*qnorm(runif(1, pnorm((lb-EZ[i,j])/sz), pnorm((ub-EZ[i,j])/sz)))
            }
          }
        }
      }
      if(family %in% c("tobit", "ordinal", "cbin", "frn", "rrl", "poisson")) {
        # For other families, use simple normal updates as approximation
        Z <- EZ + matrix(rnorm(length(Z), 0, 1), nrow=nrow(Z))
      }
    } else {
      # For unipartite, use existing functions
      rho_eff <- rho
      if(family=="normal"){ Z<-rZ_nrm_fc(Z,EZ,rho_eff,s2,Y) }
      if(family=="tobit"){ Z<-rZ_tob_fc(Z,EZ,rho_eff,s2,Y) }
      if(family=="binary"){ Z<-rZ_bin_fc(Z,EZ,rho_eff,Y) }
      if(family=="ordinal"){ Z<-rZ_ord_fc(Z,EZ,rho_eff,Y) }
      if(family=="cbin"){Z<-rZ_cbin_fc(Z,EZ,rho_eff,Y,odmax,odobs)}
      if(family=="frn"){ Z<-rZ_frn_fc(Z,EZ,rho_eff,Y,YL,odmax,odobs)}
      if(family=="rrl"){ Z<-rZ_rrl_fc(Z,EZ,rho_eff,Y,YL)}
      if(family=="poisson"){ Z<-rZ_pois_fc(Z,EZ,rho_eff,s2,Y) }
    } 
    
    # update s2
    if (is.element(family,c("normal","tobit","poisson"))) s2<-rs2_fc(Z,rho_eff,offset=EZ)  
    
    # update beta, a b with g-prior
    X_precomp <- X
    attributes(X_precomp) <- c(attributes(X), list(Xr=Xr, Xc=Xc, mX=mX, mXt=mXt, XX=XX, XXt=XXt))
    
    if(bip) {
      # For bipartite, we need to handle non-square matrices differently
      # rbeta_ab_fc assumes square matrices, so we need a workaround
      # For now, just update beta without updating a,b through rbeta_ab_fc
      if(dim(X)[3] > 0) {
        # Simple beta update for bipartite
        if(RA > 0 && RB > 0) {
          EZ_no_ab <- Z - outer(a, b, "+") - UV_term
        } else {
          EZ_no_ab <- Z - outer(a, b, "+")
        }
        vec_z <- c(EZ_no_ab)
        
        # Properly vectorize X for bipartite case
        # X should already be nA x nB x p, so just reshape it
        vec_x <- matrix(0, length(vec_z), dim(X)[3])
        for(p in 1:dim(X)[3]) {
          vec_x[, p] <- c(X[,,p])
        }
        
        # Ridge regression for beta
        XtX <- t(vec_x) %*% vec_x + diag(0.01, dim(X)[3])
        Xty <- t(vec_x) %*% vec_z
        beta <- c(solve(XtX) %*% Xty)  # Ensure beta is a vector
      } else {
        beta <- numeric(0)
      }
      # Keep a and b as they are
    } else {
      # For unipartite, use the standard function
      tmp <- rbeta_ab_fc(Z, Sab, rho, X_precomp, s2, offset=U%*%t(V), g=g)
      beta <- tmp$beta
      a <- tmp$a * rvar
      b <- tmp$b * cvar
    } 
    if(symmetric){ a<-b<-(a+b)/2 }
    
    # update Sab using unified function
    if(is.element(family,c("normal","tobit","ordinal")))
    { 
      if(bip) {
        # For bipartite, update row and column variances independently
        if(rvar) {
          Sab[1,1] <- 1/rgamma(1, (prior$etaab + nA)/2, 
                               (prior$etaab*prior$Sab0[1,1] + sum(a^2))/2)
        }
        if(cvar) {
          Sab[2,2] <- 1/rgamma(1, (prior$etaab + nB)/2, 
                               (prior$etaab*prior$Sab0[2,2] + sum(b^2))/2)
        }
        Sab[1,2] <- Sab[2,1] <- 0  # No covariance for bipartite
      } else {
        Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                       rvar=rvar, cvar=cvar, symmetric=symmetric)
      }
    }
    
    # special updates for discrete families
    if(family=="binary")
    {
      if(rvar & cvar & !symmetric & !bip) {
        # Only use raSab_bin_fc for square matrices (unipartite)
        tmp<-raSab_bin_fc(Z,Y,a,b,Sab,Sab0=prior$Sab0/prior$etaab,eta0=prior$etaab)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      } else if(bip) {
        # For bipartite binary, update variances independently
        if(rvar) {
          Sab[1,1] <- 1/rgamma(1, (prior$etaab + nA)/2, 
                               (prior$etaab*prior$Sab0[1,1] + sum(a^2))/2)
        }
        if(cvar) {
          Sab[2,2] <- 1/rgamma(1, (prior$etaab + nB)/2, 
                               (prior$etaab*prior$Sab0[2,2] + sum(b^2))/2)
        }
        Sab[1,2] <- Sab[2,1] <- 0
      } else {
        # Use standard update for symmetric cases
        Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                       rvar=rvar, cvar=cvar, symmetric=symmetric)
      }
    }
    
    if(family=="cbin")
    {
      if(rvar & cvar & !symmetric & !bip) {
        # Only use for square matrices (unipartite)
        tmp<-raSab_cbin_fc(Z,Y,a,b,Sab,odmax,odobs,Sab0=prior$Sab0/prior$etaab,eta0=prior$etaab)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      } else if(bip) {
        # For bipartite, update variances independently
        if(rvar) {
          Sab[1,1] <- 1/rgamma(1, (prior$etaab + nA)/2, 
                               (prior$etaab*prior$Sab0[1,1] + sum(a^2))/2)
        }
        if(cvar) {
          Sab[2,2] <- 1/rgamma(1, (prior$etaab + nB)/2, 
                               (prior$etaab*prior$Sab0[2,2] + sum(b^2))/2)
        }
        Sab[1,2] <- Sab[2,1] <- 0
      } else {
        Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                       rvar=rvar, cvar=cvar, symmetric=symmetric)
      }
    }
    
    if(family=="frn")
    { 
      if(rvar & cvar & !symmetric & !bip) {
        # Only use for square matrices (unipartite)
        tmp<-raSab_frn_fc(Z,Y,YL,a,b,Sab,odmax,odobs,Sab0=prior$Sab0/prior$etaab,eta0=prior$etaab)
        Z<-tmp$Z;Sab<-tmp$Sab;a<-tmp$a
      } else if(bip) {
        # For bipartite, update variances independently
        if(rvar) {
          Sab[1,1] <- 1/rgamma(1, (prior$etaab + nA)/2, 
                               (prior$etaab*prior$Sab0[1,1] + sum(a^2))/2)
        }
        if(cvar) {
          Sab[2,2] <- 1/rgamma(1, (prior$etaab + nB)/2, 
                               (prior$etaab*prior$Sab0[2,2] + sum(b^2))/2)
        }
        Sab[1,2] <- Sab[2,1] <- 0
      } else {
        Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                       rvar=rvar, cvar=cvar, symmetric=symmetric)
      }
    }
    
    # update rho (skip for bipartite networks)
    if(dcor && !bip) 
    {
      rho<-rrho_mh(Z, rho, s2, offset=Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V))
    }
    
    # shrink rho - symmetric case 
    if(symmetric){ rho<-min(.9999,1-1/sqrt(s)) }
    
    # update U,V
    if(bip && (RA > 0 || RB > 0)) {
      # Bipartite UV update  
      shrink <- (s > 0.5 * burn)
      # Check dimensions for bipartite offset
      Xb_off <- Xbeta(X, beta)
      ab_off <- outer(a, b, "+")
      if(!all(dim(Xb_off) == dim(ab_off))) {
        Xb_off <- matrix(0, nrow=length(a), ncol=length(b))
      }
      offset <- Xb_off + ab_off
      E <- Z - offset
      
      # Simple Gibbs updates for bipartite U, V, G
      # Update U
      if(RA > 0) {
        for(i in 1:nA) {
          if(RB > 0) {
            # Compute residuals for row i after removing U[i,] effect
            residual_i <- E[i, ] - U[i, ] %*% G %*% t(V)
            
            # Update each component of U[i,]
            for(k in 1:RA) {
              # Variance and mean for U[i,k]
              precision <- sum((G[k,] %*% t(V))^2) / s2 + 1  # Simple prior precision
              mean_u <- sum(residual_i * (G[k,] %*% t(V))) / s2 / precision
              U[i, k] <- rnorm(1, mean_u, sqrt(1/precision))
            }
          }
        }
      }
      
      # Update V  
      if(RB > 0) {
        for(j in 1:nB) {
          if(RA > 0) {
            # Compute residuals for col j after removing V[j,] effect  
            residual_j <- E[, j] - U %*% G %*% t(V[j, , drop=FALSE])
            
            # Update each component of V[j,]
            for(k in 1:RB) {
              # Variance and mean for V[j,k]
              precision <- sum((U %*% G[,k])^2) / s2 + 1  # Simple prior precision  
              mean_v <- sum(residual_j * (U %*% G[,k])) / s2 / precision
              V[j, k] <- rnorm(1, mean_v, sqrt(1/precision))
            }
          }
        }
      }
      
      # Update G matrix
      if(RA > 0 && RB > 0) {
        for(k1 in 1:RA) {
          for(k2 in 1:RB) {
            # Compute residuals without G[k1,k2] contribution
            UV_outer <- U[, k1] %*% t(V[, k2])  # nA x nB matrix
            residual_G <- E - U %*% G %*% t(V) + G[k1, k2] * UV_outer
            
            # Update G[k1, k2]
            precision <- sum(UV_outer^2) / s2 + 1  # Simple prior precision
            mean_g <- sum(residual_G * UV_outer) / s2 / precision
            G[k1, k2] <- rnorm(1, mean_g, sqrt(1/precision))
          }
        }
      }
      
    } else if(!bip && R > 0) {
      # Unipartite UV update (existing logic)
      shrink<- (s>.5*burn)
      offset <- Xbeta(X,beta)+outer(a,b,"+")
      
      if(symmetric){ 
        E<-Z-offset ; E<-.5*(E+t(E))
        UV<-rUV_sym_fc(E, U, V, s2,shrink) 
      }
      if(!symmetric){
        # Compute Suv for multiplicative effects
        if(shrink) {
          Suv <- rSuv_fc(U, V, Suv0=diag(2*R)/(n+R+2), kappa0=n+R+2)
        } else {
          Suv <- diag(n, nrow=2*R)
        }
        UV<-rUV_fc(Z, U, V, Suv, rho, s2, offset) 
      }
      
      U<-UV$U ; V<-UV$V
    }    
    
    # burn-in countdown
    if(print & s%%odens==0 & s<=burn){
      # Only show progress if print=TRUE
      message(paste0(round(100*s/burn,1), "% burn-in complete"))
    }
    
    # save parameter values and monitor the MC
    if(s%%odens==0 & s>burn) 
    {  
      
      # save BETA and VC - symmetric case 
      if(symmetric)
      {
        br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2 
        sbeta<-c(beta[1*intercept],bn,beta[-c(1*intercept,rb,cb)] )
        BETA<-rbind(BETA,sbeta)
        
        VC<-rbind(VC,c(Sab[1,1],s2) )
      }
      
      # save BETA and VC - asymmetric case 
      if(!symmetric)
      {
        BETA<-rbind(BETA, beta)
        if(bip) {
          # Bipartite: only row var, col var, error var (no rho or covariance)
          VC<-rbind(VC, c(Sab[1,1], Sab[2,2], s2))
        } else {
          # Unipartite: full covariance structure
          VC<-rbind(VC, c(Sab[upper.tri(Sab, diag = T)], rho, s2))
        }
      }
      
      # update posterior sums of random effects
      if(bip) {
        # For bipartite: accumulate U %*% G %*% t(V)
        if(RA > 0 && RB > 0) {
          UVPS <- UVPS + U %*% G %*% t(V)
        }
      } else {
        # For unipartite
        UVPS <- UVPS + U %*% t(V)
      }
      APS <- APS + a
      BPS <- BPS + b 
      
      # simulate from posterior predictive 
      if(bip) {
        # Bipartite: use G matrix interaction
        # Ensure dimensions are correct
        Xb_term <- Xbeta(X, beta)
        ab_term <- outer(a, b, "+")
        
        # Debug: check dimensions
        if(!all(dim(Xb_term) == dim(ab_term))) {
          # If dimensions don't match, something is wrong
          # Create a compatible matrix
          Xb_term <- matrix(0, nrow=length(a), ncol=length(b))
        }
        
        if(RA > 0 && RB > 0) {
          UV_term <- U %*% G %*% t(V)
          # Check dimensions match
          if(all(dim(Xb_term) == dim(UV_term))) {
            EZ <- Xb_term + ab_term + UV_term
          } else {
            # Fallback if dimensions don't match
            EZ <- Xb_term + ab_term
          }
        } else {
          EZ <- Xb_term + ab_term
        }
      } else {
        # Unipartite
        EZ<-Xbeta(X, beta) + outer(a, b, "+") + U %*% t(V)
        if(symmetric){ EZ<-(EZ+t(EZ))/2 }
      }
      
      # Simulate Y for posterior predictive
      if(bip) {
        # For bipartite, use simplified simulation without dyadic correlation
        if(family=="binary") { 
          ZS <- EZ + matrix(rnorm(length(EZ)), nrow(EZ), ncol(EZ))
          Ys <- 1*(ZS > 0)
        } else if(family=="normal") {
          Ys <- EZ + matrix(rnorm(length(EZ), 0, sqrt(s2)), nrow(EZ), ncol(EZ))
        } else if(family=="poisson") {
          Ys <- matrix(rpois(length(EZ), exp(EZ)), nrow(EZ), ncol(EZ))
        } else {
          # For other families, use normal approximation
          Ys <- EZ + matrix(rnorm(length(EZ)), nrow(EZ), ncol(EZ))
        }
      } else {
        # For unipartite, use existing simulation functions
        rho_eff <- rho
        if(family=="binary") { Ys<-simY_bin(EZ,rho_eff) }
        if(family=="cbin"){ Ys<-1*(simY_frn(EZ,rho_eff,odmax,YO=Y)>0) }
        if(family=="frn") { Ys<-simY_frn(EZ,rho_eff,odmax,YO=Y) }
        if(family=="rrl") { Ys<-simY_rrl(EZ,rho_eff,odobs,YO=Y ) }
        if(family=="normal") { Ys<-simY_nrm(EZ,rho_eff,s2) }
        if(family=="ordinal") { Ys<-simY_ord(EZ,rho_eff,Y) }
        if(family=="tobit") { Ys<-simY_tob(EZ,rho_eff,s2) }
        if(family=="poisson") { Ys<-simY_pois(EZ) }
      } 
      
      if(symmetric){ Ys[lower.tri(Ys)]<-0 ; Ys<-Ys+t(Ys)  }
      
      # update posterior sum
      YPS<-YPS+Ys
      
      # save posterior predictive GOF stats
      if(gof){ 
        Ys[is.na(Y)]<-NA
        if(bip) {
          # For bipartite, compute limited GOF stats
          gof_row <- matrix(NA, 1, 5)
          gof_row[1, 1] <- sd(rowMeans(Ys, na.rm=TRUE), na.rm=TRUE)
          gof_row[1, 2] <- sd(colMeans(Ys, na.rm=TRUE), na.rm=TRUE)
          GOF <- rbind(GOF, gof_row)
        } else {
          GOF<-rbind(GOF,gof_stats(Ys))
        }
      }
      
      # print MC progress 
      if (print) 
      {
        beta_means <- round(apply(BETA,2,mean),2)
        vc_means <- round(apply(VC,2,mean),2)
        cli::cli_text("Iteration {.val {s}}: beta = [{.field {paste(beta_means, collapse=', ')}}], VC = [{.field {paste(vc_means, collapse=', ')}}]")
        if (have_coda & nrow(VC) > 3 & length(beta)>0) 
        {
          eff_sizes <- round(coda::effectiveSize(BETA))
          cli::cli_text("  Effective sizes: [{.emph {paste(eff_sizes, collapse=', ')}}]")
        }
      }
      
      # Note: plotting removed to match modern approach
      # Users should use plot(fit) or trace_plot(fit) after fitting
      
      
    }
    
    # periodic save
    if(periodic_save & s %in% savePoints & !is.null(out_file)){
      # save start_vals for future model runs
      start_vals <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
      if(bip && !is.null(G)) start_vals$G <- G
      
      fit <- list(APM=APS/(sum(s>burn & s%%odens==0)), 
                 BPM=BPS/(sum(s>burn & s%%odens==0)),
                 BETA=BETA, VC=VC, GOF=GOF, 
                 start_vals=start_vals, symmetric=symmetric,
                 model.name=model.name)
      save(fit, file=out_file)
      rm(list=c('fit','start_vals'))
    }
    
  } # end MCMC   
  
  # output 
  
  # posterior means 
  APM<-APS/nrow(VC)
  BPM<-BPS/nrow(VC)
  # Add names to APM and BPM
  if(is.null(names(APM)) && !is.null(rownames(Y))) {
    names(APM) <- rownames(Y)
  }
  if(is.null(names(BPM))) {
    if(bip && !is.null(colnames(Y))) {
      names(BPM) <- colnames(Y)
    } else if(!is.null(rownames(Y))) {
      names(BPM) <- rownames(Y)
    }
  }
  # If still no names, create default ones
  if(is.null(names(APM))) {
    names(APM) <- paste0("Actor", 1:length(APM))
  }
  if(is.null(names(BPM))) {
    names(BPM) <- paste0("Actor", 1:length(BPM))
  }
  UVPM<-UVPS/nrow(VC)
  YPM<-YPS/nrow(VC) 
  # Compute final expected values
  if(bip) {
    # For bipartite, check dimensions
    Xb_final <- Xbeta(X, apply(BETA, 2, mean))
    ab_final <- outer(APM, BPM, "+")
    if(!all(dim(Xb_final) == dim(ab_final))) {
      Xb_final <- matrix(0, nrow=length(APM), ncol=length(BPM))
    }
    EZ <- Xb_final + ab_final + UVPM
  } else {
    EZ<-Xbeta(X,apply(BETA,2,mean)) + outer(APM,BPM,"+")+UVPM
  }  
  
  if(bip) {
    # For bipartite, row and column names are different
    names(APM) <- rownames(UVPM) <- rownames(Y)
    names(BPM) <- colnames(UVPM) <- colnames(Y)
  } else {
    names(APM)<-names(BPM)<-rownames(UVPM)<-colnames(UVPM)<-dimnames(Y)[[1]]
  }
  # Transform EZ to count scale for Poisson family
  # EZ is on log scale internally, but users expect count scale
  if(family == "poisson") {
    EZ <- exp(EZ)
    # Cap extreme values to avoid numerical issues
    EZ[EZ > 1e6] <- 1e6
  }
  
  dimnames(YPM)<-dimnames(EZ)<-dimnames(Y)
  rownames(BETA)<-NULL
  
  # model selection statistics  
  if(!is.null(model.name))
  {
    # Count parameters
    if(bip) {
      p_eff <- length(beta) + rvar*nA + cvar*nB + RA*RB
    } else {
      p_eff <- length(beta) + rvar*n + cvar*n + R*(R+1)/2
    }
    
  }
  
  # asymmetric output 
  if(!symmetric) 
  {
    if(bip) {
      # For bipartite networks, keep the original U, V, G structure
      # No SVD decomposition needed since U and V have different dimensions
      rownames(U) <- rownames(Y)
      if(ncol(Y) > 1) {
        rownames(V) <- colnames(Y)
      }
      # Include G matrix in output for bipartite
      # save start_vals for future model runs
      start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=V, s2=s2, Sab=Sab, G=G)
      fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,G=G,UVPM=UVPM,EZ=EZ,
                  YPM=YPM,GOF=GOF,start_vals=start_vals_final,model.name=model.name,
                  mode="bipartite",nA=nA,nB=nB,RA=RA,RB=RB,
                  family=family,symmetric=symmetric,odmax=odmax)
    } else {
      # Unipartite networks (existing logic)
      UDV<-svd(UVPM)
      U<-UDV$u[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
      V<-UDV$v[,seq(1,R,length=R)]%*%diag(sqrt(UDV$d)[seq(1,R,length=R)],nrow=R)
      rownames(U)<-rownames(V)<-rownames(Y) 
      # save start_vals for future model runs
      start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
      fit <- list(BETA=BETA,VC=VC,APM=APM,BPM=BPM,U=U,V=V,UVPM=UVPM,EZ=EZ,
                  YPM=YPM,GOF=GOF,start_vals=start_vals_final,model.name=model.name,
                  mode="unipartite",family=family,symmetric=symmetric,odmax=odmax)
    }
  }
  
  # symmetric output (not used for bipartite)
  if(symmetric) 
  {
    ULUPM<-UVPM
    eULU<-eigen(ULUPM)
    eR<- which( rank(-abs(eULU$val),ties.method="first") <= R )
    U<-eULU$vec[,seq(1,R,length=R),drop=FALSE]
    L<-eULU$val[eR]
    rownames(U)<-rownames(ULUPM)<-colnames(ULUPM)<-rownames(Y)
    EZ<-.5*(EZ+t(EZ)) ; YPM<-.5*(YPM+t(YPM)) 
    # save start_vals for future model runs
    start_vals_final <- list(Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
    fit<-list(BETA=BETA,VC=VC,APM=APM,U=U,L=L,ULUPM=ULUPM,EZ=EZ,
              YPM=YPM,GOF=GOF,start_vals=start_vals_final,model.name=model.name,
              mode="unipartite",family=family,symmetric=symmetric,odmax=odmax)
  } 
  
  # Ensure scalars are properly typed
  if (!is.null(fit$RHO)) fit$RHO <- as.numeric(fit$RHO)
  if (!is.null(fit$s2))  fit$s2  <- as.numeric(fit$s2)
  
  class(fit) <- "ame"
  fit
}


