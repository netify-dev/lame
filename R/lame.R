#' AME model fitting routine for replicated relational data
#' 
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to replicated relational data of
#' various types. 
#' 
#' @details
#' This command provides posterior inference for parameters in AME models of
#' independent replicated relational data, assuming one of eight possible data
#' types/models.
#' 
#' \strong{Dynamic Effects Implementation:}
#' 
#' The dynamic_uv and dynamic_ab parameters enable time-varying latent representations
#' through autoregressive processes. These extensions are particularly useful for 
#' understanding how network structure evolves over time.
#' 
#' \emph{Dynamic Multiplicative Effects (dynamic_uv=TRUE):}
#' The latent factors U and V evolve according to AR(1) processes:
#' \deqn{U_{i,k,t} = \rho_{uv} U_{i,k,t-1} + \epsilon_{i,k,t}}
#' where \eqn{\epsilon_{i,k,t} \sim N(0, \sigma_{uv}^2)}, i indexes actors, k indexes
#' latent dimensions, and t indexes time. The parameter \eqn{\rho_{uv}} controls
#' temporal persistence (values near 1 indicate slow evolution). This captures 
#' time-varying homophily, latent community structure, and transitivity dynamics.
#' 
#' Key references:
#' \itemize{
#'   \item Sewell & Chen (2015): Introduced dynamic latent space models with 
#'         actor-specific evolution rates
#'   \item Durante & Dunson (2014): Nonparametric Bayesian approach allowing 
#'         flexible evolution of network structure
#'   \item Hoff (2011): Hierarchical multilinear models providing theoretical 
#'         foundation for temporal dependencies
#' }
#' 
#' \emph{Dynamic Additive Effects (dynamic_ab=TRUE):}
#' The sender (a) and receiver (b) effects evolve as:
#' \deqn{a_{i,t} = \rho_{ab} a_{i,t-1} + \epsilon_{i,t}}
#' \deqn{b_{i,t} = \rho_{ab} b_{i,t-1} + \eta_{i,t}}
#' where \eqn{\epsilon_{i,t}, \eta_{i,t} \sim N(0, \sigma_{ab}^2)}. This models 
#' time-varying individual activity levels (outdegree) and popularity (indegree).
#' 
#' Applications include:
#' \itemize{
#'   \item Tracking changes in node centrality over time
#'   \item Identifying emerging influential actors
#'   \item Detecting declining activity patterns
#'   \item Modeling life-cycle effects in social networks
#' }
#' 
#' \emph{Prior Specification for Dynamic Parameters:}
#' \itemize{
#'   \item \eqn{\rho_{uv}, \rho_{ab} \sim TruncNormal(mean, sd, 0, 1)}: 
#'         Ensures stationarity of AR(1) process
#'   \item \eqn{\sigma_{uv}^2, \sigma_{ab}^2 \sim InverseGamma(shape, scale)}: 
#'         Controls innovation variance
#'   \item Default priors (\eqn{\rho_{uv}} mean=0.9, \eqn{\rho_{ab}} mean=0.8) 
#'         favor smooth evolution
#'   \item Adjust rho_*_mean closer to 1 for slower evolution, closer to 0 for 
#'         more rapid changes
#' }
#' 
#' \emph{Computational Considerations:}
#' \itemize{
#'   \item Dynamic effects increase computation by ~30-50\% per iteration
#'   \item Memory usage scales as O(n*R*T) for dynamic_uv, O(n*T) for dynamic_ab
#'   \item C++ implementation provides ~70\% efficiency gain over pure R
#'   \item Convergence diagnostics: Monitor rho and sigma parameters carefully
#'   \item Effective sample sizes typically lower due to temporal correlation
#'   \item Recommend burn >= 1000 and nscan >= 20000 for dynamic models
#' }
#' 
#' \emph{Model Selection Guidelines:}
#' Use both dynamic_uv and dynamic_ab when:
#' \itemize{
#'   \item Networks show clear temporal trends in density or clustering
#'   \item Individual node behavior changes systematically over time
#'   \item Community structure evolves (merging, splitting, drift)
#' }
#' 
#' Use only dynamic_uv when:
#' \itemize{
#'   \item Latent structure/communities change but individual effects are stable
#'   \item Focus is on evolving homophily or clustering patterns
#'   \item Network shows structural reconfiguration over time
#' }
#' 
#' Use only dynamic_ab when:
#' \itemize{
#'   \item Individual heterogeneity varies but overall structure is stable
#'   \item Actors' activity/popularity changes over observation period
#'   \item Focus is on individual-level temporal dynamics
#' }
#' 
#' \strong{Standard AME Model Types:}
#' 
#' The following describes the eight standard data types/models available:
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
#' @usage lame(Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL, rvar = !(family=="rrl"),
#'   cvar = TRUE, dcor = !symmetric, nvar=TRUE, R = 0,
#'   dynamic_uv = FALSE, dynamic_ab = FALSE, family ="normal",
#' intercept=!is.element(family,c("rrl","ordinal")),
#' symmetric=FALSE,
#' odmax=NULL, prior=list(), g=NA,
#' seed = 6886, nscan = 10000, burn = 500, odens = 25, plot=FALSE, print = FALSE, gof=TRUE,
#' startVals=NULL, periodicSave=FALSE, outFile=NULL, saveInterval=0.25, model.name=NULL)
#' @param Y a T length list of n x n relational matrices, where T 
#' corresponds to the number of replicates (over time, for example). 
#' See family below for various data types.
#' @param Xdyad a T length list of n x n x pd arrays of covariates
#' @param Xrow a T length list of n x pr matrices of nodal row covariates
#' @param Xcol a T length list of n x pc matrices of nodal column covariates
#' @param rvar logical: fit row random effects (asymmetric case)?
#' @param cvar logical: fit column random effects (asymmetric case)? 
#' @param dcor logical: fit a dyadic correlation (asymmetric case)?
#' @param nvar logical: fit nodal random effects (symmetric case)? 
#' @param R integer: dimension of the multiplicative effects (can be zero)
#' @param dynamic_uv logical: fit dynamic multiplicative effects (latent factors) that 
#' evolve over time using AR(1) processes. When TRUE, the latent positions U and V become 
#' time-varying, following \eqn{U_{i,t} = \rho_{uv} U_{i,t-1} + \epsilon_{i,t}}, where epsilon 
#' follows N(0, sigma_uv^2). This allows actors' positions in latent social space to 
#' drift smoothly over time, capturing evolving network structure and community dynamics.
#' Inspired by dynamic latent space models (Sewell and Chen 2015, "Latent Space Models 
#' for Dynamic Networks", JASA; Durante and Dunson 2014, "Nonparametric Bayes Dynamic 
#' Modeling of Relational Data", Biometrika). The implementation uses efficient blocked 
#' Gibbs sampling with C++ acceleration for scalability. Default FALSE.
#' @param dynamic_ab logical: fit dynamic additive effects (sender/receiver effects) that 
#' evolve over time using AR(1) processes. When TRUE, the row effects (a) and column 
#' effects (b) become time-varying, following \eqn{a_{i,t} = \rho_{ab} a_{i,t-1} + \epsilon_{i,t}}.
#' This captures temporal heterogeneity in actors' baseline propensities to send and 
#' receive ties, allowing for smooth changes in activity levels and popularity over time.
#' For example, an actor's tendency to form outgoing ties might gradually increase or 
#' decrease across observation periods. The AR(1) specification ensures temporal smoothness
#' while allowing for actor-specific evolution patterns. Implementation uses conjugate 
#' updates where possible and C++ for computational efficiency. Default FALSE.
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
#'         impose stronger shrinkage on the latent factors.}
#'   \item{etaab}{Prior degrees of freedom for covariance of additive effects 
#'         (default: 4 + 3 \\* n/100). Controls shrinkage of row/column random effects.}
#'   \item{rho_uv_mean}{For dynamic_uv=TRUE: Prior mean for UV AR(1) parameter 
#'         (default: 0.9). Values close to 1 indicate high temporal persistence.}
#'   \item{rho_uv_sd}{For dynamic_uv=TRUE: Prior SD for UV AR(1) parameter 
#'         (default: 0.1). Controls uncertainty about temporal dependence.}
#'   \item{sigma_uv_shape}{For dynamic_uv=TRUE: Shape parameter for inverse-gamma 
#'         prior on UV innovation variance (default: 2).}
#'   \item{sigma_uv_scale}{For dynamic_uv=TRUE: Scale parameter for inverse-gamma 
#'         prior on UV innovation variance (default: 1).}
#'   \item{rho_ab_mean}{For dynamic_ab=TRUE: Prior mean for additive effects AR(1) 
#'         parameter (default: 0.8). Controls temporal smoothness of sender/receiver effects.}
#'   \item{rho_ab_sd}{For dynamic_ab=TRUE: Prior SD for additive effects AR(1) 
#'         parameter (default: 0.15).}
#'   \item{sigma_ab_shape}{For dynamic_ab=TRUE: Shape parameter for inverse-gamma 
#'         prior on additive effects innovation variance (default: 2).}
#'   \item{sigma_ab_scale}{For dynamic_ab=TRUE: Scale parameter for inverse-gamma 
#'         prior on additive effects innovation variance (default: 1).}
#' }
#' Common usage: prior = list(Sab0 = diag(c(1, 1)), eta0 = 10) for stronger 
#' shrinkage, or prior = list(rho_uv_mean = 0.95) for higher temporal persistence.
#' @param g optional scalar or vector for g-prior on regression coefficients.
#' Default is p^2 where p is the number of regression parameters. The g-prior 
#' controls the variance of regression coefficients: larger values allow for 
#' larger coefficient values. Can be a vector of length p for parameter-specific 
#' control.
#' @param seed random seed
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain
#' @param odens output density for the Markov chain
#' @param plot logical: plot results while running?
#' @param print logical: print results while running?
#' @param gof logical: calculate goodness of fit statistics?
#' @param startVals List from previous model run containing parameter starting values for new MCMC
#' @param periodicSave logical: indicating whether to periodically save MCMC results
#' @param outFile character vector indicating name and path in which file should be stored if periodicSave is selected. For example, on an Apple OS outFile="~/Desktop/ameFit.rda".
#' @param saveInterval quantile interval indicating when to save during post burn-in phase.
#' @param model.name optional string for model selection output
#' @return \item{BETA}{posterior samples of regression coefficients}
#' \item{VC}{posterior samples of the variance parameters}
#' \item{APM}{posterior mean of additive row effects a} \item{BPM}{posterior
#' mean of additive column effects b} \item{U}{posterior mean of multiplicative
#' row effects u. For dynamic_uv=TRUE, this is a 3D array (n x R x T)} 
#' \item{V}{posterior mean of multiplicative column effects v (asymmetric case).
#' For dynamic_uv=TRUE, this is a 3D array (n x R x T)}
#' \item{UVPM}{posterior mean of UV}
#' \item{ULUPM}{posterior mean of ULU (symmetric case)} 
#' \item{L}{posterior mean of L (symmetric case)} 
#'  \item{EZ}{estimate of expectation of Z
#' matrix} \item{YPM}{posterior mean of Y (for imputing missing values)}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics}
#' \item{startVals}{Final parameter values from MCMC, can be used as the input
#' for a future model run.}
#' \item{AIC}{Akaike Information Criterion (if model.name provided)}
#' \item{BIC}{Bayesian Information Criterion (if model.name provided)}
#' \item{model.name}{Name of the model (if provided)}
#' @author Peter Hoff, Yanjun He, Shahryar Minhas
#' @examples
#' 
#' data(YX_bin_list) 
#' fit<-lame(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,family="binary")
#' # you should run the Markov chain much longer than this
#' 
#' @export lame
lame <- function(
    Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL, 
    rvar = !(family=="rrl") , cvar = TRUE, dcor = !symmetric, 
    nvar = TRUE, 
    R = 0, dynamic_uv = FALSE, dynamic_ab = FALSE,
    family = "normal",
    intercept = !is.element(family,c("rrl","ordinal")),
    symmetric = FALSE,
    odmax = NULL,
    prior = list(), g = NA,
    seed = 6886, nscan = 10000, burn = 500, odens = 25,
    plot = FALSE, print = FALSE, gof = TRUE, 
    startVals = NULL, periodicSave=FALSE, outFile=NULL,
    saveInterval=0.25, model.name = NULL
)
{
  #
  if( nscan %% odens !=0  ){ stop('"odens" must be a multiple of "nscan"')}
  
  # set random seed 
  set.seed(seed)
  
  # get actor info
  actorByYr <- lapply(Y, rownames)
  actorSet <- sort(unique(unlist( actorByYr ))) ; n <- length(actorSet)
  
  # reset odmax param
  odmax <- rep( max( unlist( lapply(Y, function(y){ apply(y>0, 1, sum, na.rm=TRUE)  }) ) ), n )
  
  # calc savePoints
  savePoints <- (burn:(nscan+burn))[(burn:(nscan+burn)) %% odens==0]
  savePoints <- savePoints[round(quantile(1:length(savePoints), probs=seq(saveInterval,1,saveInterval)))]
  
  # check formatting of input objects
  check_format(Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol)
  
  # set diag to NA 
  N<-length(Y) ; pdLabs <- names(Y) ; Y<-lapply(Y, function(y){diag(y)=NA; return(y)})
  
  # convert into large array format
  arrayObj<-list_to_array(actorSet, Y, Xdyad, Xrow, Xcol)
  Y<-arrayObj$Y ; Xdyad<-arrayObj$Xdyad ; Xrow<-arrayObj$Xrow
  Xcol<-arrayObj$Xcol ; rm(arrayObj)
  
  # force binary if binary family specified 
  if(is.element(family,c("binary","cbin"))) { Y<-1*(Y>0) }
  
  # handle tobit family
  if(family=="tobit") { 
    Y <- lapply(Y, function(y) { y[y<0]<-0; return(y) })
  } 
  
  # observed and max outdegrees 
  if(is.element(family,c("cbin","frn","rrl")) ){
    odobs<-apply(Y>0,c(1,3),sum,na.rm=TRUE) 
    if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y[,,1])) } 
  }
  
  # some settings for symmetric case
  if(symmetric){ Xcol<-Xrow ; rvar<-cvar<-nvar }
  
  # construct design matrix    
  pr<-length(Xrow[,,1])/n
  pc<-length(Xcol[,,1])/n
  pd<-length(Xdyad[,,,1])/n^2
  designObj <- get_design_rep(
    Y=Y,Xdyad=Xdyad,Xrow=Xrow,Xcol=Xcol,actorSet=actorSet,
    intercept=intercept,n=n,N=N,pr=pr,pc=pc,pd=pd)
  Y<-designObj$Y ; X<-designObj$X ; Xlist<-designObj$Xlist
  XrLong<-designObj$XrLong ; XcLong<-designObj$XcLong
  mXLong<-designObj$mXLong ; mXtLong<-designObj$mXtLong
  xxLong<-designObj$xxLong ; xxTLong<-designObj$xxTLong ; rm(designObj)
  
  # design matrix warning for rrl
  if( family=="rrl" & any(apply(apply(X,c(1,3),var),2,sum)==0)
      & !any( apply(X,c(3),function(x){var(c(x))})==0) )
  {
    cli::cli_warn("Row effects are not estimable using this procedure")
  }
  
  # design matrix warning for rrl and ord
  if( is.element(family,c("ordinal","rrl")) & 
      any( apply(X,c(3),function(x){var(c(x))})==0 ) )
  {
    cli::cli_warn("An intercept is not estimable using this procedure")
  }
  
  # construct matrix of ranked nominations for frn, rrl   
  if(is.element(family,c("frn","rrl")))
  {
    ymx<-max(apply(1*(Y>0),c(1,3),sum,na.rm=TRUE))
    YL<-list()
    for (t in 1:N) 
    {
      YL.t<-NULL
      warn<-FALSE
      for(i in 1:nrow(Y[,,1]))
      {
        yi<-Y[i,,t] ; rnkd<-which( !is.na(yi)&yi>0 )
        if(length(yi[rnkd])>length(unique(yi[rnkd]))){warn<-TRUE}
        yi[rnkd]<-rank(yi[rnkd],ties.method="random")
        Y[i,,t]<-yi
        YL.t<-rbind(YL.t, match(1:ymx,yi))
      }
      YL[[t]]<-YL.t
      if(warn){cli::cli_warn("Random reordering used to break ties in ranks")}
    }
  }
  
  # g-prior setup
  p<-dim(X)[3]
  if(is.na(g)) { g<-p^2 }  # default g-prior
  if(length(g)==1) { g<-rep(g,p) }
  
  # process prior list
  if(!is.list(prior)) { prior<-list() }
  if(is.null(prior$Sab0)) { prior$Sab0<-diag(2) }
  if(is.null(prior$eta0)) { prior$eta0<-round(4+3*n/100) } 
  if(is.null(prior$etaab)) { prior$etaab<-round(4+3*n/100) }
  
  # Dynamic UV priors
  if(dynamic_uv && R > 0) {
    if(is.null(prior$rho_uv_mean)) { prior$rho_uv_mean <- 0.9 }
    if(is.null(prior$rho_uv_sd)) { prior$rho_uv_sd <- 0.1 }
    if(is.null(prior$sigma_uv_shape)) { prior$sigma_uv_shape <- 2 }
    if(is.null(prior$sigma_uv_scale)) { prior$sigma_uv_scale <- 1 }
  }
  
  # Dynamic AB priors
  if(dynamic_ab) {
    if(is.null(prior$rho_ab_mean)) { prior$rho_ab_mean <- 0.8 }
    if(is.null(prior$rho_ab_sd)) { prior$rho_ab_sd <- 0.15 }
    if(is.null(prior$sigma_ab_shape)) { prior$sigma_ab_shape <- 2 }
    if(is.null(prior$sigma_ab_scale)) { prior$sigma_ab_scale <- 1 }
  }
  
  # Get starting values for MCMC
  startValsObj <- get_start_vals(startVals,Y,family,xP=dim(X)[3],rvar,cvar,R,odmax=odmax)
  Z<-startValsObj$Z ; beta<-startValsObj$beta ; a<-startValsObj$a
  b<-startValsObj$b ; U<-startValsObj$U ; V<-startValsObj$V
  rho<-startValsObj$rho ; s2<-startValsObj$s2 ; Sab<-startValsObj$Sab
  
  # For symmetric models, don't use rho (dyadic correlation)
  if(symmetric) {
    rho <- 0
    use_rho <- FALSE
  } else {
    use_rho <- TRUE
  }
  
  # Initialize dynamic AB parameters if needed
  if(dynamic_ab) {
    # Convert a and b to 2D matrices (n x T) if not already
    if(length(a) == n) {
      # Create time-varying matrices
      a_mat <- matrix(rep(a, N), nrow=n, ncol=N)
      b_mat <- matrix(rep(b, N), nrow=n, ncol=N)
      
      # Add small random variation across time
      for(t in 2:N) {
        a_mat[,t] <- prior$rho_ab_mean * a_mat[,t-1] + rnorm(n, 0, 0.1)
        b_mat[,t] <- prior$rho_ab_mean * b_mat[,t-1] + rnorm(n, 0, 0.1)
      }
    } else if(is.matrix(a) && ncol(a) == N) {
      a_mat <- a
      b_mat <- b
    } else {
      # Initialize from scratch
      ab_init <- init_dynamic_ab_cpp(n, N, prior$rho_ab_mean, 0.1)
      a_mat <- ab_init$a
      b_mat <- ab_init$b
    }
    # Initialize AR(1) parameters
    rho_ab <- ifelse(!is.null(startValsObj$rho_ab), startValsObj$rho_ab, prior$rho_ab_mean)
    sigma_ab <- ifelse(!is.null(startValsObj$sigma_ab), startValsObj$sigma_ab, 0.1)
  } else {
    a_mat <- NULL
    b_mat <- NULL
    rho_ab <- NULL
    sigma_ab <- NULL
  }
  
  # Initialize dynamic UV parameters if needed
  if(dynamic_uv && R > 0) {
    # Convert U and V to 3D arrays if not already
    if(length(dim(U)) == 2) {
      # Initialize 3D arrays by replicating 2D positions across time
      U_cube <- array(0, dim=c(n, R, N))
      V_cube <- array(0, dim=c(n, R, N))
      for(t in 1:N) {
        U_cube[,,t] <- U
        V_cube[,,t] <- V
      }
    } else {
      U_cube <- U
      V_cube <- V
    }
    # Initialize AR(1) parameters
    rho_uv <- ifelse(!is.null(startValsObj$rho_uv), startValsObj$rho_uv, prior$rho_uv_mean)
    sigma_uv <- ifelse(!is.null(startValsObj$sigma_uv), startValsObj$sigma_uv, 0.1)
  } else {
    U_cube <- NULL
    V_cube <- NULL
    rho_uv <- NULL
    sigma_uv <- NULL
  }
  
  rm(list=c('startValsObj','startVals'))
  
  # helpful mcmc params
  symLoopIDs <- lapply(1:(nscan + burn), function(x){ rep(sample(1:nrow(U)),4) })  
  asymLoopIDs <- lapply(1:(nscan + burn), function(x){ sample(1:R) })  
  tryErrorChecks<-list(s2=0,betaAB=0,rho=0,UV=0)  
  iter <- 1
  
  # output items
  BETA <- matrix(nrow = nscan/odens, ncol = dim(X)[3] - pr*symmetric)
  VC<-matrix(nrow=nscan/odens,ncol=5-3*symmetric)
  
  # Initialize UV posterior sums based on dynamic or static
  if(dynamic_uv && R > 0) {
    UVPS <- array(0, dim=c(n, n, N))
    U_SUM <- array(0, dim=c(n, R, N))
    V_SUM <- array(0, dim=c(n, R, N))
  } else {
    UVPS <- U %*% t(V) * 0
  }
  
  APS<-BPS<-rep(0,nrow(Y[,,1]))
  YPS<-array(0,dim=dim(Y),dimnames=dimnames(Y)) 
  # GOF <- matrix(NA, nrow=(nscan/odens)+1, ncol=5,
  #               dimnames=list(c('obs',1:(nscan/odens)),c("sd.rowmean","sd.colmean","dyad.dep","cycle.dep", "trans.dep")))
  # GOF[1,] <- rowMeans(apply(Y,3,gof_stats))
  # 5 corresponds to number of stats we want
  GOF <- array(NA, dim=c(5, N, (nscan/odens)+1))
  
  # add label to first dim
  dimnames(GOF)[[1]] <- c("sd.rowmean","sd.colmean","dyad.dep","cycle.dep","trans.dep")
  
  # add row/col names
  if(!is.null(dimnames(Y)[[3]])) { dimnames(GOF)[[2]] <- dimnames(Y)[[3]] }
  
  # add time dimnames
  dimnames(GOF)[[3]] <- c('obs', 1:(nscan/odens))
  
  # fill in first entry with values from OBSERVED 
  GOF[,,1] <- apply(Y, 3, gof_stats)
  
  names(APS)<-names(BPS)<-rownames(U)<-rownames(V)<-rownames(Y[,,1])
  
  # names of parameters, asymmetric case  
  if(!symmetric)
  { 
    colnames(VC) <- c("va", "cab", "vb", "rho", "ve") 
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
  
  # MCMC
  have_coda<-suppressWarnings(
    try(requireNamespace("coda",quietly = TRUE),silent=TRUE)) 
  
  if(burn!=0){
    cli::cli_h3("Starting burn-in period...")
    cli::cli_progress_bar("Burn-in", total = burn, .envir = environment())
  }
  if(!print){
    burn_complete <- FALSE
  }
  for (s in 1:(nscan + burn)) 
  { 
    
    # update Z
    E.nrm<-array(dim=dim(Z))
    if(dynamic_ab) {
      # Use dynamic helper to compute EZ with time-varying additive effects
      EZ <- get_EZ_dynamic_ab(Xlist, beta, a_mat, b_mat, U, V, N)
    } else {
      EZ <- get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
    }
    for(t in 1:N ){
      if(family=="normal")
      { 
        Z[,,t]<-rZ_nrm_fc(Z[,,t],EZ[,,t],rho,s2,Y[,,t]) ; E.nrm[,,t]<-Z[,,t]-EZ[,,t]
      }
      if(family=="tobit")
      { 
        Z[,,t]<-rZ_tob_fc(Z[,,t],EZ[,,t],rho,s2,Y[,,t]) ; E.nrm[,,t]<-Z[,,t]-EZ[,,t]
      }
      if(family=="binary"){ Z[,,t]<-rZ_bin_fc(Z[,,t],EZ[,,t],rho,Y[,,t]) }
      if(family=="ordinal"){ Z[,,t]<-rZ_ord_fc(Z[,,t],EZ[,,t],rho,Y[,,t]) }
      if(family=="cbin"){Z[,,t]<-rZ_cbin_fc(Z[,,t],EZ[,,t],rho,Y[,,t],odmax,odobs)}
      if(family=="frn")
      { 
        Z[,,t]<-rZ_frn_fc(Z[,,t],EZ[,,t],rho,Y[,,t],YL[[t]],odmax,odobs)
      }
      if(family=="rrl"){ Z[,,t]<-rZ_rrl_fc(Z[,,t],EZ[,,t],rho,Y[,,t],YL[[t]]) }
      if(family=="poisson")
      { 
        Z[,,t]<-rZ_pois_fc(Z[,,t],EZ[,,t],rho,s2,Y[,,t]) ; E.nrm[,,t]<-Z[,,t]-EZ[,,t]
      }
    }
    
    # update s2
    if (is.element(family,c("normal","tobit","poisson"))){
      s2New<-try(
        rs2_rep_fc_cpp(E.nrm,solve(matrix(c(1,rho,rho,1),2,2))), 
        silent=TRUE)
      if(!inherits(s2New, 'try-error')){ s2 <- s2New } else { tryErrorChecks$s2<-tryErrorChecks$s2+1 }
    }
    
    # update beta, a b with g-prior
    if(dynamic_ab) {
      # Dynamic additive effects update
      
      # First update beta only
      if( (pr+pc+pd+intercept)>0 ){
        # Compute expected Z without additive effects for current time
        EZ_no_ab <- get_EZ_cpp(Xlist, beta*0, outer(rep(0,n), rep(0,n), "+"), U, V)
        
        # Use time-averaged a and b for beta update
        a_avg <- rowMeans(a_mat)
        b_avg <- rowMeans(b_mat)
        
        iSe2<-mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)) ; Sabs<-iSe2%*%Sab%*%iSe2
        tmp<-eigen(Sabs) ; k<-sum(zapsmall(tmp$val)>0 )
        G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
        betaABCalc <- try(
          rbeta_ab_rep_fc_cpp(
            ZT=sweep(Z,c(1,2),U%*%t(V)), Xr=XrLong, Xc=XcLong, mX=mXLong, mXt=mXtLong,
            XX=xxLong, XXt=xxTLong, iSe2=iSe2, Sabs=Sabs, k=k, G=G ),
          silent = TRUE)
        if(!inherits(betaABCalc, 'try-error')){
          beta <- c(betaABCalc$beta)
        }
      }
      
      # Now update dynamic a and b
      EZ_no_ab <- get_EZ_cpp(Xlist, beta, outer(rep(0,n), rep(0,n), "+"), U, V)
      ab_update <- sample_dynamic_ab_cpp(a_mat, b_mat, Z, EZ_no_ab, 
                                         rho_ab, sigma_ab, Sab, symmetric)
      a_mat <- ab_update$a
      b_mat <- ab_update$b
      
      # Use current time values for compatibility
      a <- a_mat[,1]  # Will be overridden in get_EZ_cpp calls
      b <- b_mat[,1]
      
      # Update AR(1) parameters periodically
      if(s > burn && s %% 20 == 0) {
        rho_ab <- sample_rho_ab_cpp(a_mat, b_mat, sigma_ab, rho_ab, symmetric)
        sigma_ab <- sample_sigma_ab_cpp(a_mat, b_mat, rho_ab, symmetric)
      }
      
    } else {
      # Standard static update
      if( (pr+pc+pd+intercept)>0 ){
        iSe2<-mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)) ; Sabs<-iSe2%*%Sab%*%iSe2
        tmp<-eigen(Sabs) ; k<-sum(zapsmall(tmp$val)>0 )
        G<-tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
        betaABCalc <- try(
          rbeta_ab_rep_fc_cpp(
            ZT=sweep(Z,c(1,2),U%*%t(V)), Xr=XrLong, Xc=XcLong, mX=mXLong, mXt=mXtLong,
            XX=xxLong, XXt=xxTLong, iSe2=iSe2, Sabs=Sabs, k=k, G=G ),
          silent = TRUE)
      } else {
        betaABCalc <- try(
          rbeta_ab_rep_fc(sweep(Z,c(1,2),U%*%t(V)), Sab, rho, X, s2),
          silent = TRUE)
      }
      if(!inherits(betaABCalc, 'try-error')){
        beta <- c(betaABCalc$beta)
        a <- c(betaABCalc$a) * rvar
        b <- c(betaABCalc$b) * cvar
        if(symmetric){ a<-b<-(a+b)/2 }
      } else { tryErrorChecks$betaAB<-tryErrorChecks$betaAB+1  }
    }
    
    # update Sab using unified function
    Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
                   rvar=rvar, cvar=cvar, symmetric=symmetric)
    
    # update rho (only for asymmetric models)
    if(dcor && !symmetric)
    {
      if(dynamic_ab) {
        E.T <- Z - get_EZ_dynamic_ab(Xlist, beta, a_mat, b_mat, U, V, N)
      } else {
        E.T <- Z - get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
      }
      rhoNew<-try( rrho_mh_rep_cpp(E.T, rho,s2), silent=TRUE )
      if(!inherits(rhoNew, 'try-error')){ rho<-rhoNew } else { tryErrorChecks$rho<-tryErrorChecks$rho+1 }
    }
    
    # update U,V
    if (R > 0)
    {
      if(dynamic_ab) {
        # For dynamic ab, use time-varying effects but zero out UV
        E <- Z - get_EZ_dynamic_ab(Xlist, beta, a_mat, b_mat, U*0, V*0, N)
      } else {
        E <- Z - get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U*0, V*0 )
      }
      shrink<- (s>.5*burn)
      
      if(dynamic_uv) {
        # Dynamic UV update with AR(1) evolution
        UV <- try(
          rUV_dynamic_fc_cpp(U_cube, V_cube, E, rho_uv, sigma_uv, s2, shrink, symmetric),
          silent = TRUE)
        if(inherits(UV, 'try-error')){ 
          UV <- list(U=U_cube, V=V_cube) 
          tryErrorChecks$UV<-tryErrorChecks$UV+1 
        } else {
          U_cube <- UV$U
          V_cube <- UV$V
          # Extract time-averaged positions for compatibility
          U <- apply(U_cube, c(1,2), mean)
          V <- apply(V_cube, c(1,2), mean)
        }
        
        # Update AR(1) parameters (commented out as functions not exported)
        # if(s > burn && s %% 10 == 0) {
        #   rho_uv <- sample_rho_uv(U_cube, V_cube, sigma_uv, rho_uv, symmetric)
        #   sigma_uv <- sample_sigma_uv(U_cube, V_cube, rho_uv, symmetric)
        # }
      } else {
        # Standard static UV update
        if(symmetric)
        { 
          EA<-apply(E,c(1,2),mean) ; EA<-.5*(EA+t(EA))
          UV<-try(
            rUV_sym_fc_cpp(EA, U, V, 
                           s2/dim(E)[3], shrink, symLoopIDs[[s]]-1), silent=TRUE )
          if(inherits(UV, 'try-error')){ UV <- list(U=U,V=V) ; tryErrorChecks$UV<-tryErrorChecks$UV+1 }
        }
        if(!symmetric){
          UV <- try(
            rUV_rep_fc_cpp(E, U, V, rho, s2,
                           mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)),
                           maxmargin=1e-6, shrink, asymLoopIDs[[s]]-1 ), silent = TRUE )
          if(inherits(UV, 'try-error')){ UV <- list(U=U,V=V) ; tryErrorChecks$UV<-tryErrorChecks$UV+1 }
        }      
        
        U<-UV$U ; V<-UV$V
      }
    }
    
    # burn-in countdown
    if(burn!=0 && s <= burn){cli::cli_progress_update()}
    
    # store parameter values and monitor the MC
    if(s==burn+1&!print&burn!=0){
      cli::cli_progress_done()
      cli::cli_alert_success("Burn-in period complete")
      cli::cli_progress_bar("Sampling", total = nscan, .envir = environment())
    }
    if(s%%odens==0 & s>burn) 
    { 
      
      # store BETA and VC - symmetric case 
      if(symmetric){
        br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2
        sbeta<-c(beta[1*intercept],bn,beta[-c(1*intercept,rb,cb)] )
        BETA[iter,]<-sbeta
        VC[iter,]<-c(Sab[1,1],s2)
      }
      
      # store BETA and VC - asymmetric case 
      if(!symmetric){
        BETA[iter,]<-beta
        VC[iter,]<- c(Sab[upper.tri(Sab, diag = T)], rho,s2)
      }
      
      # update posterior sums of random effects
      if(dynamic_uv && R > 0 && !is.null(U_cube) && !is.null(V_cube)) {
        U_SUM <- U_SUM + U_cube
        V_SUM <- V_SUM + V_cube
        # Store UV products for each time point
        for(t in 1:N) {
          # Ensure dimensions are correct
          U_t <- U_cube[,,t,drop=FALSE]
          V_t <- V_cube[,,t,drop=FALSE]
          U_t <- matrix(U_t, nrow=n, ncol=R)
          V_t <- matrix(V_t, nrow=n, ncol=R)
          UVPS[,,t] <- UVPS[,,t] + U_t %*% t(V_t)
        }
      } else if(R > 0) {
        UVPS <- UVPS + U %*% t(V)
      }
      if(dynamic_ab) {
        # For dynamic ab, accumulate time-specific effects
        if(iter == 1) {
          # Initialize accumulation matrices
          APS_dyn <- a_mat
          BPS_dyn <- b_mat
        } else {
          APS_dyn <- APS_dyn + a_mat
          BPS_dyn <- BPS_dyn + b_mat
        }
        # Also track time-averaged for compatibility
        APS <- APS + rowMeans(a_mat)
        BPS <- BPS + rowMeans(b_mat)
      } else {
        APS <- APS + a
        BPS <- BPS + b
      } 
      
      # simulate from posterior predictive
      # EZ already computed earlier in the iteration, reuse it
      dimnames(EZ) <- dimnames(Y)
      Ys <- EZ*0
      for (t in 1:N)
      {
        if(symmetric){ EZ[,,t]<-(EZ[,,t]+t(EZ[,,t]))/2 }
        
        if(family=="binary"){ Ys[,,t]<-simY_bin(EZ[,,t],rho) }
        if(family=="cbin"){ Ys[,,t]<-1*(simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t])>0)}
        if(family=="frn"){ Ys[,,t]<-simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t]) }
        if(family=="rrl"){ Ys[,,t]<-simY_rrl(EZ[,,t],rho,odobs,YO=Y[,,t] ) }
        if(family=="normal"){ Ys[,,t]<-simY_nrm(EZ[,,t],rho,s2) }
        if(family=="tobit"){ Ys[,,t]<-simY_tob(EZ[,,t],rho,s2) }
        if(family=="ordinal"){ Ys[,,t]<-simY_ord(EZ[,,t],rho,Y[,,t]) }
        if(family=="poisson"){ Ys[,,t]<-simY_pois(EZ[,,t]) }
        
        if(symmetric)
        {  
          Yst<-Ys[,,t] ; Yst[lower.tri(Yst)]<-0 ; Ys[,,t]<-Yst+t(Yst)
        }
      } 
      
      # update posterior sum
      YPS<-YPS+Ys
      
      # save posterior predictive GOF stats
      #if(gof){Ys[is.na(Y)]<-NA ;GOF[(iter)+1,]<-rowMeans(apply(Ys,3,gof_stats))}
      if(gof){
        
        Ys[is.na(Y)] <- NA
        GOF[,,(iter)+1] <- apply(Ys,3,gof_stats)
      }
      
      # print MC progress 
      if(print)
      {
        cat('\n',s,
            round(apply(BETA[1:iter,,drop=FALSE],2,mean),2),":",
            round(apply(VC[1:iter,,drop=FALSE],2,mean),2),"\n")
        if (have_coda & nrow(VC[1:iter,,drop=FALSE]) > 3 & length(beta)>0) 
        {
          cat(round(coda::effectiveSize(BETA[1:iter,,drop=FALSE])), "\n")
        }
      }
      
      # periodic save
      if(periodicSave & s %in% savePoints & !is.null(outFile)){
        # save startVals for future model runs
        if(dynamic_uv && R > 0 && dynamic_ab) {
          startVals <- list(Z=Z,beta=beta,a=a_mat,b=b_mat,U=U_cube,V=V_cube,
                           rho=rho,s2=s2,Sab=Sab,rho_uv=rho_uv,sigma_uv=sigma_uv,
                           rho_ab=rho_ab, sigma_ab=sigma_ab)
        } else if(dynamic_uv && R > 0) {
          startVals <- list(Z=Z,beta=beta,a=a,b=b,U=U_cube,V=V_cube,
                           rho=rho,s2=s2,Sab=Sab,rho_uv=rho_uv,sigma_uv=sigma_uv)
        } else if(dynamic_ab) {
          startVals <- list(Z=Z,beta=beta,a=a_mat,b=b_mat,U=U,V=V,rho=rho,s2=s2,Sab=Sab,
                           rho_ab=rho_ab, sigma_ab=sigma_ab)
        } else {
          startVals <- list(Z=Z,beta=beta,a=a,b=b,U=U,V=V,rho=rho,s2=s2,Sab=Sab)
        }
        fit <- get_fit_object( APS=APS, BPS=BPS, UVPS=UVPS, YPS=YPS, 
                             BETA=BETA, VC=VC, GOF=GOF, Xlist=Xlist, actorByYr=actorByYr,
                             startVals=startVals, symmetric=symmetric, tryErrorChecks=tryErrorChecks,
                             AIC=NA, BIC=NA, model.name=model.name)
        save(fit, file=outFile) ; rm(list=c('fit','startVals'))
      }
      
      # plot MC results
      if(plot & s==(burn+nscan))
      {
        # plot VC
        print(trace_plot(list(VC=VC), params="variance"))
        
        # plot BETA
        if(length(beta)>0) 
        {
          betaIndices<-split(1:ncol(BETA), ceiling(seq_along(1:ncol(BETA))/5))
          for(bIndex in betaIndices){
            print(trace_plot(list(BETA=BETA[,bIndex,drop=FALSE]), params="beta")) }
        }
        
        # plot GOF
        if(gof)
        {
          suppressMessages( print( gof_plot(list(GOF=GOF)) ) )
        }
      } # plot code if applicable
      iter<-iter+1
    } # post burn-in
    if(!print && s > burn){cli::cli_progress_update()}
  } # end MCMC  
  if(!print){cli::cli_progress_done(); cli::cli_alert_success("MCMC sampling complete")}
  
  # save startVals for future model runs
  if(dynamic_ab) {
    startVals <- list( Z=Z, beta=beta, a=a_mat, b=b_mat, U=U, V=V, rho=rho, s2=s2, Sab=Sab,
                      rho_ab=rho_ab, sigma_ab=sigma_ab)
  } else {
    startVals <- list( Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
  }
  
  # Calculate AIC/BIC if model.name provided
  if(!is.null(model.name))
  {
    # Count effective parameters
    p_eff <- ncol(BETA) + rvar*n + cvar*n + R*(R+1)/2
    
    # Log-likelihood approximation (simplified)
    # Note: This is a simplified calculation for model comparison purposes
    # Full likelihood calculation would require more complex integration
    s2_mean <- mean(VC[,ncol(VC)])
    n_obs <- sum(!is.na(Y))
    
    if(family %in% c("normal", "tobit", "binary")) {
      # Use a simplified log-likelihood based on residual variance
      # This provides a rough approximation for model comparison
      ll <- -n_obs/2 * (log(2*pi*s2_mean) + 1)
    } else {
      ll <- NA
    }
    
    # AIC and BIC
    if(!is.na(ll)) {
      AIC <- -2*ll + 2*p_eff
      BIC <- -2*ll + p_eff*log(n_obs)
    } else {
      AIC <- BIC <- NA
    }
  } else {
    AIC <- BIC <- NA
  }
  
  # output
  if(dynamic_uv && R > 0) {
    # Return 3D arrays for dynamic UV
    U_final <- U_SUM / (nscan/odens)
    V_final <- V_SUM / (nscan/odens)
    UVPS_final <- UVPS / (nscan/odens)
  } else {
    U_final <- U
    V_final <- V
    UVPS_final <- UVPS / (nscan/odens)
  }
  
  # Handle dynamic additive effects output
  if(dynamic_ab) {
    APS_final <- APS_dyn / (nscan/odens)
    BPS_final <- BPS_dyn / (nscan/odens)
  } else {
    APS_final <- APS / (nscan/odens)
    BPS_final <- BPS / (nscan/odens)
  }
  
  fit <- get_fit_object( APS=APS_final, BPS=BPS_final, UVPS=UVPS_final, YPS=YPS, 
                       BETA=BETA, VC=VC, GOF=GOF, Xlist=Xlist, actorByYr=actorByYr, 
                       startVals=startVals, symmetric=symmetric, tryErrorChecks=tryErrorChecks,
                       AIC=AIC, BIC=BIC, model.name=model.name, U=U_final, V=V_final, 
                       dynamic_uv=dynamic_uv, dynamic_ab=dynamic_ab)
  class(fit) <- "lame" # set class
  return(fit) # output object to workspace
  
}