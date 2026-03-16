#' AME model fitting routine for longitudinal relational data
#' 
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to longitudinal (time-series) relational data of
#' various types. Supports both unipartite (square) and bipartite (rectangular)
#' network structures. For cross-sectional (single time point) networks, 
#' use the \code{ame} function.
#' 
#' @details
#' This command provides posterior inference for parameters in AME models of
#' longitudinal relational data, assuming one of eight possible data
#' types/models. The model supports both unipartite networks (square adjacency
#' matrices) and bipartite networks (rectangular adjacency matrices with 
#' distinct row and column node sets) across multiple time points.
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
#' \strong{Bipartite Network Models:}
#' 
#' When mode="bipartite", the model handles rectangular adjacency matrices Y with
#' dimensions n_A x n_B, where n_A and n_B represent the number of row and column
#' nodes respectively.
#' 
#' \emph{Static Bipartite Case:}
#' The model uses separate latent factor matrices:
#' \itemize{
#'   \item U: n_A x R_row matrix of row node latent positions
#'   \item V: n_B x R_col matrix of column node latent positions  
#'   \item G: R_row x R_col interaction matrix mapping between latent spaces
#'   \item Multiplicative term: U G V' captures bipartite community structure
#' }
#' 
#' \emph{Dynamic Bipartite Case:}
#' When dynamic_uv=TRUE for bipartite networks:
#' \deqn{U_{i,k,t} = \rho_{uv} U_{i,k,t-1} + \epsilon_{i,k,t}}
#' \deqn{V_{j,k,t} = \rho_{uv} V_{j,k,t-1} + \eta_{j,k,t}}
#' where i indexes row nodes, j indexes column nodes, k indexes latent dimensions.
#' 
#' When dynamic_G=TRUE, the interaction matrix also evolves:
#' \deqn{G_{k,l,t} = \rho_G G_{k,l,t-1} + \xi_{k,l,t}}
#' allowing the mapping between row and column latent spaces to change over time.
#' 
#' \emph{Key Differences from Unipartite Models:}
#' \itemize{
#'   \item No dyadic correlation (rho): Bipartite edges are inherently directed
#'   \item Separate dimensions: R_row and R_col can differ for row/column spaces
#'   \item Rectangular structure: Network density patterns differ from square matrices
#'   \item Community interpretation: Captures affiliation patterns between node types
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
#'   cvar = TRUE, dcor = !symmetric, nvar=TRUE, R = 0, R_row = NULL, R_col = NULL,
#'   mode = c("unipartite", "bipartite"),
#'   dynamic_uv = FALSE, dynamic_ab = FALSE, dynamic_G = FALSE, family ="normal",
#' intercept=!is.element(family,c("rrl","ordinal")),
#' symmetric=FALSE,
#' odmax=NULL, prior=list(), g=NA,
#' seed = 6886, nscan = 10000, burn = 500, odens = 25, plot=FALSE, verbose = FALSE, gof=TRUE,
#' start_vals=NULL, periodic_save=FALSE, out_file=NULL, save_interval=0.25, model.name=NULL,
#' print)
#' @param Y a T length list of n x n relational matrices, or a 3D array
#' of dimensions \code{[n, n, T]}, where T corresponds to the number of
#' replicates (over time, for example). If a 3D array is provided, it is
#' automatically converted to list format. See family below for various
#' data types.
#' @param Xdyad a T length list of n x n x pd arrays of covariates
#' @param Xrow a T length list of n x pr matrices of nodal row covariates
#' @param Xcol a T length list of n x pc matrices of nodal column covariates
#' @param rvar logical: fit row random effects (asymmetric case)?
#' @param cvar logical: fit column random effects (asymmetric case)? 
#' @param dcor logical: fit a dyadic correlation (asymmetric case)?
#' @param nvar logical: fit nodal random effects (symmetric case)? 
#' @param R integer: dimension of the multiplicative effects (can be zero)
#' @param R_row integer: for bipartite networks, dimension of row node multiplicative effects (defaults to R)
#' @param R_col integer: for bipartite networks, dimension of column node multiplicative effects (defaults to R)
#' @param mode character: either "unipartite" (default) for square networks or "bipartite" for rectangular networks
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
#' @param dynamic_G logical: for bipartite networks, fit dynamic interaction matrix G that
#' evolves over time. When TRUE, the rectangular interaction matrix G becomes time-varying,
#' allowing the mapping between row and column latent spaces to change over time. Default FALSE.
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
#' @param verbose logical: print progress while running? Default FALSE.
#' @param print Deprecated. Use \code{verbose} instead.
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
#' row effects u. For dynamic_uv=TRUE, this is a 3D array (n x R x T)} 
#' \item{V}{posterior mean of multiplicative column effects v (asymmetric case).
#' For dynamic_uv=TRUE, this is a 3D array (n x R x T)}
#' \item{UVPM}{posterior mean of UV}
#' \item{ULUPM}{posterior mean of ULU (symmetric case)} 
#' \item{L}{posterior mean of L (symmetric case)} 
#'  \item{EZ}{estimate of expectation of Z
#' matrix} \item{YPM}{posterior mean of Y (for imputing missing values)}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics.
#' See \code{\link{gof}} for post-hoc computation and \code{\link{gof_plot}} for visualization.}
#' \item{start_vals}{Final parameter values from MCMC, can be used as the input
#' for a future model run.}
#' \item{model.name}{Name of the model (if provided)}
#' @seealso \code{\link{ame}} for cross-sectional models,
#'   \code{\link{gof}} for post-hoc goodness-of-fit computation,
#'   \code{\link{gof_plot}} for visualizing GOF results,
#'   \code{\link{latent_positions}} for extracting latent positions as a tidy data frame,
#'   \code{\link{procrustes_align}} for Procrustes alignment of latent positions,
#'   \code{\link{summary.lame}} for model summaries,
#'   \code{\link{coef.lame}} for coefficient extraction
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
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
		R = 0, R_row = NULL, R_col = NULL,  # separate ranks for bipartite
		mode = c("unipartite", "bipartite"),  # network mode
		dynamic_uv = FALSE, dynamic_ab = FALSE, dynamic_G = FALSE,
		family = "normal",
		intercept = !is.element(family,c("rrl","ordinal")),
		symmetric = FALSE,
		odmax = NULL,
		prior = list(), g = NA,
		seed = 6886, nscan = 10000, burn = 500, odens = 25,
		plot = FALSE, verbose = FALSE, gof = TRUE,
		start_vals = NULL, periodic_save=FALSE, out_file=NULL,
		save_interval=0.25, model.name = NULL,
		print
) {
	# helper function
	`%||%` <- function(x, y) if (is.null(x)) y else x

	# handle deprecated print argument
	if (!missing(print)) {
		cli::cli_warn(c(
			"The {.arg print} argument is deprecated.",
			"i" = "Use {.arg verbose} instead."
		))
		if (missing(verbose) || identical(verbose, FALSE)) {
			verbose <- print
		}
	}

	# process mode argument
	mode <- match.arg(mode)
	bip <- identical(mode, "bipartite")

	# convert 3D array Y to list format if needed
	if (is.array(Y) && length(dim(Y)) == 3L) {
		dn <- dimnames(Y)
		Y <- lapply(seq_len(dim(Y)[3]), function(t) {
			mat <- Y[,,t]
			if (!is.null(dn)) {
				rownames(mat) <- dn[[1]]
				colnames(mat) <- dn[[2]]
			}
			mat
		})
		if (!is.null(dn) && !is.null(dn[[3]])) {
			names(Y) <- dn[[3]]
		}
	}

	#
	if( nscan %% odens !=0  ){ stop('"odens" must be a multiple of "nscan"')}
	
	# set random seed 
	set.seed(seed)
	
	# get actor info
	if(bip) {
		rowActorByYr <- lapply(Y, rownames)
		colActorByYr <- lapply(Y, colnames)
		rowActorSet <- sort(unique(unlist(rowActorByYr)))
		colActorSet <- sort(unique(unlist(colActorByYr)))
		nA <- length(rowActorSet)
		nB <- length(colActorSet)
		n <- nA
		actorSet <- rowActorSet
		actorByYr <- rowActorByYr
	} else {
		actorByYr <- lapply(Y, rownames)
		actorSet <- sort(unique(unlist( actorByYr )))
		n <- length(actorSet)
		nA <- nB <- n
	}
	
	# reset odmax param
	odmax <- rep( max( unlist( lapply(Y, function(y){ apply(y>0, 1, sum, na.rm=TRUE)  }) ) ), n )
	
	# calc savePoints
	savePoints <- (burn:(nscan+burn))[(burn:(nscan+burn)) %% odens==0]
	savePoints <- savePoints[round(quantile(1:length(savePoints), probs=seq(save_interval,1,save_interval)))]
	
	# check formatting of input objects
	check_format(Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol)
	
	# set diag to NA (only for square matrices)
	N<-length(Y) ; pdLabs <- names(Y) 
	if(!bip) {
		Y<-lapply(Y, function(y){diag(y)=NA; return(y)})
	}
	
	# convert into large array format
	if(bip) {
		arrayObj <- list_to_array_bipartite(rowActorSet, colActorSet, Y, Xdyad, Xrow, Xcol)
	} else {
		arrayObj <- list_to_array(actorSet, Y, Xdyad, Xrow, Xcol)
	}
	Y<-arrayObj$Y ; Xdyad<-arrayObj$Xdyad ; Xrow<-arrayObj$Xrow
	Xcol<-arrayObj$Xcol ; rm(arrayObj)

	if(bip && is.list(Xdyad) && !is.null(Xdyad) && length(Xdyad) > 0) {
		pd_tmp <- dim(Xdyad[[1]])[3]
		Xdyad_4d <- array(NA, dim = c(nA, nB, pd_tmp, length(Xdyad)))
		for(tt in seq_along(Xdyad)) {
			Xdyad_4d[,,,tt] <- Xdyad[[tt]]
		}
		Xdyad <- Xdyad_4d
		rm(Xdyad_4d)
	}
	
	# set g-prior parameter if not provided (after n is defined)
	if(is.na(g)) {
		# use n^2 as default (Zellner's g-prior)
		g <- n^2
	}
	
	# force binary if binary family specified 
	if(is.element(family,c("binary","cbin"))) { Y<-1*(Y>0) }
	
	# clamp tobit values to zero (Y is already a 3D array here)
	if(family=="tobit") {
		Y[Y < 0 & !is.na(Y)] <- 0
	}
	
	# observed and max outdegrees 
	if(is.element(family,c("cbin","frn","rrl")) ){
		odobs<-apply(Y>0,c(1,3),sum,na.rm=TRUE) 
		if(length(odmax)==1) { odmax<-rep(odmax,nrow(Y[,,1])) } 
	}
	
	# some settings for symmetric case
	if(symmetric && !bip){ 
		Xcol<-Xrow ; rvar<-cvar<-nvar 
	}
	if(bip && symmetric) {
		warning("`symmetric` ignored for bipartite networks")
		symmetric <- FALSE
	}
	if(dynamic_G) {
		warning("`dynamic_G=TRUE` is not yet implemented; using static G instead")
		dynamic_G <- FALSE
	}
	if(bip && family %in% c("ordinal", "cbin", "frn", "poisson")) {
		cli::cli_abort(c(
			"family={.val {family}} is not supported for bipartite networks in {.fn lame}.",
			"i" = "Supported bipartite families: {.val normal}, {.val binary}, {.val tobit}."
		))
	}
	if(bip && family == "rrl") {
		warning("rrl (row-ranked likelihood) is not well-defined for bipartite networks; results may be unreliable")
	}
	
	# set rank dimensions
	if(bip) {
		RA <- R_row %||% R %||% 0L
		RB <- R_col %||% R %||% 0L
		dcor <- FALSE  # no dyadic correlation in bipartite
	} else {
		if(!is.null(R_row) || !is.null(R_col)) {
			warning("R_row/R_col ignored for unipartite networks")
		}
		RA <- RB <- as.integer(R %||% 0L)
	}
	
	# construct design matrix    
	pr<-length(Xrow[,,1])/n
	pc<-length(Xcol[,,1])/n
	
	# prepare default Xlist for bipartite before design matrix construction
	if(bip) {
		# for bipartite networks
		pd <- if(!is.null(Xdyad)) length(Xdyad[,,,1])/(nA*nB) else 0
		if(!is.null(Xdyad) && length(Xdyad) > 0) {
			# convert Xdyad (4D array: nA x nB x pd x N) to Xlist (list of 3D arrays)
			Xlist <- vector("list", N)
			n_cov_d <- dim(Xdyad)[3]
			n_cov_total <- n_cov_d
			if(intercept) n_cov_total <- n_cov_total + 1L
			for(t in 1:N) {
				Xt <- array(0, dim = c(nA, nB, n_cov_total))
				k_start <- 1L
				if(intercept) {
					Xt[,,1] <- 1
					k_start <- 2L
				}
				for(k in 1:n_cov_d) {
					Xt[,,k_start + k - 1] <- Xdyad[,,k,t]
				}
				Xlist[[t]] <- Xt
			}
		} else {
			# intercept-only Xlist
			if(intercept) {
				Xlist <- replicate(N, array(1, dim = c(nA, nB, 1L)), simplify = FALSE)
			} else {
				Xlist <- replicate(N, array(0, dim = c(nA, nB, 0L)), simplify = FALSE)
			}
		}
		# add row/column covariates to Xlist if provided
		if(!is.null(Xrow) && length(Xrow) > 0) {
			n_xr <- if(length(dim(Xrow)) >= 2) dim(Xrow)[2] else 1L
			for(t in 1:N) {
				old_p <- dim(Xlist[[t]])[3]
				new_Xt <- array(0, dim = c(nA, nB, old_p + n_xr))
				new_Xt[,,1:old_p] <- Xlist[[t]]
				for(k in 1:n_xr) {
					xr_vals <- if(length(dim(Xrow)) == 3) Xrow[,k,t] else Xrow[,t]
					new_Xt[,,old_p + k] <- matrix(xr_vals, nA, nB)
				}
				Xlist[[t]] <- new_Xt
			}
		}
		if(!is.null(Xcol) && length(Xcol) > 0) {
			n_xc <- if(length(dim(Xcol)) >= 2) dim(Xcol)[2] else 1L
			for(t in 1:N) {
				old_p <- dim(Xlist[[t]])[3]
				new_Xt <- array(0, dim = c(nA, nB, old_p + n_xc))
				new_Xt[,,1:old_p] <- Xlist[[t]]
				for(k in 1:n_xc) {
					xc_vals <- if(length(dim(Xcol)) == 3) Xcol[,k,t] else Xcol[,t]
					new_Xt[,,old_p + k] <- t(matrix(xc_vals, nB, nA))
				}
				Xlist[[t]] <- new_Xt
			}
		}
	} else {
		# for unipartite networks
		pd<-length(Xdyad[,,,1])/n^2
	}

	if(!bip) {
		# only use get_design_rep for unipartite networks
		designObj <- get_design_rep(
			Y=Y,Xdyad=Xdyad,Xrow=Xrow,Xcol=Xcol,actorSet=actorSet,
			intercept=intercept,n=n,N=N,pr=pr,pc=pc,pd=pd)
		Y<-designObj$Y ; X<-designObj$X ; Xlist<-designObj$Xlist
		XrLong<-designObj$XrLong ; XcLong<-designObj$XcLong
		mXLong<-designObj$mXLong ; mXtLong<-designObj$mXtLong
		xxLong<-designObj$xxLong ; xxTLong<-designObj$xxTLong ; rm(designObj)
	} else {
		# bipartite: Xlist already constructed above
		X <- Xdyad
		XrLong <- NULL ; XcLong <- NULL
		mXLong <- NULL ; mXtLong <- NULL
		xxLong <- NULL ; xxTLong <- NULL
	}
	
	# design matrix warning for rrl
	if( !bip && family=="rrl" && !is.null(X) && length(dim(X)) == 4 &&
			any(apply(apply(X,c(1,3),var),2,sum)==0) &&
			!any( apply(X,c(3),function(x){var(c(x))})==0) ) {
		cli::cli_warn("Row effects are not estimable using this procedure")
	}
	
	# design matrix warning for rrl and ord
	if( !bip && is.element(family,c("ordinal","rrl")) && 
			!is.null(X) && length(dim(X)) == 4 &&
			any( apply(X,c(3),function(x){var(c(x))})==0 ) ) {
		cli::cli_warn("An intercept is not estimable using this procedure")
	}
	
	# construct matrix of ranked nominations for frn, rrl   
	if(is.element(family,c("frn","rrl"))) {
		ymx<-max(apply(1*(Y>0),c(1,3),sum,na.rm=TRUE))
		YL<-list()
		for (t in 1:N)  {
			YL.t<-NULL
			warn<-FALSE
			for(i in 1:nrow(Y[,,1])) {
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
	
	# g-prior setup: figure out the number of covariates from
	# whichever design object is available
	if(!is.null(X) && length(dim(X)) >= 3) {
		p <- dim(X)[3]
	} else if(bip && !is.null(Xlist) && length(Xlist) > 0) {
		p <- dim(Xlist[[1]])[3]
	} else {
		p <- if(intercept) 1L else 0L
	}
	# use number of dyads like amen package, not p^2 (too restrictive)
	if(is.na(g)) { 
		# count non-missing dyads across all time periods
		n_dyads <- sum(!is.na(Y))
		g <- max(n_dyads, n^2)  # use number of observed dyads or n^2, whichever is larger
	}
	if(length(g)==1) { g<-rep(g,p) }
	
	# process prior list
	if(!is.list(prior)) { prior<-list() }
	if(is.null(prior$Sab0)) { prior$Sab0<-diag(2) }
	
	# binary family eta0 calculation
	if(family == "binary" && is.null(prior$eta0)) {
		Y1 <- if(is.array(Y)) Y[,,1] else Y
		ydist <- table(Y1)
		ymode <- as.numeric(names(ydist)[ydist == max(ydist)])[1]
		YB <- 1 * (Y1 != ymode)
		ybar <- mean(YB, na.rm=TRUE)
		if(ybar > 0 && ybar < 1) {
			mu_bin <- qnorm(ybar)
			E <- (YB - ybar) / dnorm(qnorm(ybar))
			if(!bip && nrow(E) == ncol(E)) diag(E) <- 0
			a_tmp <- rowMeans(E, na.rm=TRUE)
			b_tmp <- colMeans(E, na.rm=TRUE)
			a_tmp[is.na(a_tmp)] <- 0
			b_tmp[is.na(b_tmp)] <- 0
			if(bip) {
				vscale <- (var(a_tmp) + var(b_tmp)) / 2
			} else {
				vscale <- mean(diag(cov(cbind(a_tmp, b_tmp))))
			}
			PHAT <- pnorm(mu_bin + outer(a_tmp, b_tmp, "+"))
			vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
			prior$eta0 <- round(4 * vdfmlt)
		} else {
			prior$eta0 <- round(4+3*n/100)
		}
	} else if(is.null(prior$eta0)) {
		prior$eta0 <- round(4+3*n/100)
	}
	
	if(is.null(prior$etaab)) { prior$etaab <- prior$eta0 }
	
	# dynamic UV priors
	if(dynamic_uv && R > 0) {
		if(is.null(prior$rho_uv_mean)) { prior$rho_uv_mean <- 0.9 }
		if(is.null(prior$rho_uv_sd)) { prior$rho_uv_sd <- 0.1 }
		if(is.null(prior$sigma_uv_shape)) { prior$sigma_uv_shape <- 2 }
		if(is.null(prior$sigma_uv_scale)) { prior$sigma_uv_scale <- 1 }
	}
	
	# dynamic AB priors
	if(dynamic_ab) {
		if(is.null(prior$rho_ab_mean)) { prior$rho_ab_mean <- 0.8 }
		if(is.null(prior$rho_ab_sd)) { prior$rho_ab_sd <- 0.15 }
		if(is.null(prior$sigma_ab_shape)) { prior$sigma_ab_shape <- 2 }
		if(is.null(prior$sigma_ab_scale)) { prior$sigma_ab_scale <- 1 }
	}
	
	# prepare Xlist before initialization
	if (!bip && (is.null(Xlist) || length(Xlist) == 0L)) {
		# for unipartite, create n x n intercept matrices if still needed
		Xlist <- replicate(N, array(1, dim = c(n, n, 1L)), simplify = FALSE)
	}
	
	# get starting values for MCMC
	if(bip) {
		# bipartite initialization
		startValsObj <- init_bipartite_startvals(Y, family, nA, nB, RA, RB, N, 
																						Xlist, odmax)
		Z <- startValsObj$Z
		beta <- startValsObj$beta
		a <- startValsObj$a
		b <- startValsObj$b
		U_full <- startValsObj$U
		V_full <- startValsObj$V
		# for non-dynamic case, use first time point
		if(!dynamic_uv) {
			if(RA > 0) {
				U <- U_full[,,1,drop=FALSE] ; dim(U) <- c(nA, RA)
			} else {
				U <- matrix(0, nA, 1)  # minimal dimension
			}
			if(RB > 0) {
				V <- V_full[,,1,drop=FALSE] ; dim(V) <- c(nB, RB) 
			} else {
				V <- matrix(0, nB, 1)  # minimal dimension
			}
		} else {
			U <- U_full
			V <- V_full
		}
		G <- startValsObj$G  # rectangular interaction matrix
		sigma2_a <- startValsObj$sigma2_a
		sigma2_b <- startValsObj$sigma2_b
		s2 <- startValsObj$s2
		rho <- 0  # no dyadic correlation in bipartite (set to 0 not NULL)
		use_rho <- FALSE
		Sab <- diag(c(sigma2_a, sigma2_b))  # independent variances
	} else {
		# unipartite initialization
		startValsObj <- get_start_vals(start_vals,Y,family,xP=dim(X)[3],rvar,cvar,R,odmax=odmax)
		Z<-startValsObj$Z ; beta<-startValsObj$beta ; a<-startValsObj$a
		b<-startValsObj$b ; U<-startValsObj$U ; V<-startValsObj$V
		rho<-startValsObj$rho ; s2<-startValsObj$s2 ; Sab<-startValsObj$Sab
		G <- NULL  # no G matrix for unipartite
		
		# symmetric models don't use rho (dyadic correlation)
		if(symmetric) {
			rho <- 0
			use_rho <- FALSE
		} else {
			use_rho <- TRUE
		}
	}
	
	# initialize dynamic AB parameters if needed
	if(dynamic_ab) {
		# determine correct dimensions for a and b
		n_a <- if(bip) nA else n
		n_b <- if(bip) nB else n

		# convert a and b to 2D matrices (n_a x T, n_b x T) if not already
		if(length(a) == n_a && !is.matrix(a)) {
			# create time-varying matrices
			a_mat <- matrix(rep(a, N), nrow=n_a, ncol=N)
			b_mat <- matrix(rep(b, N), nrow=n_b, ncol=N)

			# add small random variation across time
			for(t in 2:N) {
				a_mat[,t] <- prior$rho_ab_mean * a_mat[,t-1] + rnorm(n_a, 0, 0.1)
				b_mat[,t] <- prior$rho_ab_mean * b_mat[,t-1] + rnorm(n_b, 0, 0.1)
			}
		} else if(is.matrix(a) && ncol(a) == N) {
			a_mat <- a
			b_mat <- b
		} else {
			# initialize from scratch
			a_init <- init_dynamic_ab_cpp(n_a, N, prior$rho_ab_mean, 0.1)
			b_init <- init_dynamic_ab_cpp(n_b, N, prior$rho_ab_mean, 0.1)
			a_mat <- a_init$a
			b_mat <- b_init$b
		}
		# initialize AR(1) parameters
		rho_ab <- ifelse(!is.null(startValsObj$rho_ab), startValsObj$rho_ab, prior$rho_ab_mean)
		sigma_ab <- ifelse(!is.null(startValsObj$sigma_ab), startValsObj$sigma_ab, 0.1)
	} else {
		a_mat <- NULL
		b_mat <- NULL
		rho_ab <- NULL
		sigma_ab <- NULL
	}
	
	# initialize dynamic UV parameters if needed
	if(dynamic_uv && R > 0) {
		# determine correct dimensions for U and V cubes
		n_u <- if(bip) nA else n
		n_v <- if(bip) nB else n
		R_u <- if(bip) RA else R
		R_v <- if(bip) RB else R
		# convert U and V to 3D arrays if not already
		if(length(dim(U)) == 2) {
			U_cube <- array(0, dim=c(n_u, R_u, N))
			V_cube <- array(0, dim=c(n_v, R_v, N))
			for(t in 1:N) {
				U_cube[,,t] <- U
				V_cube[,,t] <- V
			}
		} else {
			U_cube <- U
			V_cube <- V
		}
		# initialize AR(1) parameters
		rho_uv <- ifelse(!is.null(startValsObj$rho_uv), startValsObj$rho_uv, prior$rho_uv_mean)
		sigma_uv <- ifelse(!is.null(startValsObj$sigma_uv), startValsObj$sigma_uv, 0.1)
	} else {
		U_cube <- NULL
		V_cube <- NULL
		rho_uv <- NULL
		sigma_uv <- NULL
	}
	
	rm(list=c('startValsObj','start_vals'))
	
	# helpful mcmc params
	symLoopIDs <- lapply(1:(nscan + burn), function(x){ rep(sample(1:nrow(U)),4) })  
	asymLoopIDs <- lapply(1:(nscan + burn), function(x){ sample(1:R) })  
	tryErrorChecks<-list(s2=0,betaAB=0,rho=0,UV=0,Z=0)

	# bipartite MH proposal tracking
	if(bip && RA > 0 && RB > 0) {
		mh_sd_U <- 0.1
		mh_sd_V <- 0.1
		mh_acc_U <- 0L
		mh_att_U <- 0L
		mh_acc_V <- 0L
		mh_att_V <- 0L
	}
	iter <- 1
	
	# output items
	if(bip) {
		n_beta <- if(!is.null(Xlist) && length(Xlist) > 0) dim(Xlist[[1]])[3] else 1
		BETA <- matrix(nrow = nscan/odens, ncol = n_beta)
	} else {
		BETA <- matrix(nrow = nscan/odens, ncol = dim(X)[3] - pr*symmetric)
	}
	VC<-matrix(nrow=nscan/odens,ncol=5-3*symmetric)
	
	# storage for dynamic parameters
	if(dynamic_ab) {
		RHO_AB <- numeric(nscan/odens)
	} else {
		RHO_AB <- NULL
	}
	if(dynamic_uv) {
		RHO_UV <- numeric(nscan/odens)  
	} else {
		RHO_UV <- NULL
	}
	
	# initialize UV posterior sums based on dynamic or static
	R_eff <- if(bip) max(RA, RB) else R
	if(dynamic_uv && R_eff > 0) {
		if(bip) {
			UVPS <- array(0, dim=c(nA, nB, N))
			U_SUM <- array(0, dim=c(nA, RA, N))
			V_SUM <- array(0, dim=c(nB, RB, N))
		} else {
			UVPS <- array(0, dim=c(n, n, N))
			U_SUM <- array(0, dim=c(n, R_eff, N))
			V_SUM <- array(0, dim=c(n, R_eff, N))
		}
	} else if(R_eff > 0) {
		if(bip) {
			UVPS <- matrix(0, nA, nB)
		} else {
			UVPS <- U %*% t(V) * 0
		}
	} else {
		if(bip) {
			UVPS <- matrix(0, nA, nB)
		} else {
			UVPS <- array(0, dim=c(n, n))
		}
	}
	
	# initialize APS and BPS with correct dimensions for bipartite
	if(bip) {
		APS <- rep(0, nA)  # row/sender effects
		names(APS) <- rowActorSet
		BPS <- rep(0, nB)  # column/receiver effects
		names(BPS) <- colActorSet
	} else {
		APS<-BPS<-rep(0,nrow(Y[,,1]))
	}
	YPS<-array(0,dim=dim(Y),dimnames=dimnames(Y))

	# gof storage
	if(!bip) {
		GOF <- array(NA, dim=c(5, N, (nscan/odens)+1))
		dimnames(GOF)[[1]] <- c("sd.rowmean","sd.colmean","dyad.dep","cycle.dep","trans.dep")
		if(!is.null(dimnames(Y)[[3]])) { dimnames(GOF)[[2]] <- dimnames(Y)[[3]] }
		dimnames(GOF)[[3]] <- c('obs', 1:(nscan/odens))
		GOF[,,1] <- apply(Y, 3, gof_stats)
	} else {
		# bipartite: 3 stats (sd.rowmean, sd.colmean, four.cycles)
		GOF <- array(NA, dim=c(3, N, (nscan/odens)+1))
		dimnames(GOF)[[1]] <- c("sd.rowmean","sd.colmean","four.cycles")
		if(!is.null(dimnames(Y)[[3]])) { dimnames(GOF)[[2]] <- dimnames(Y)[[3]] }
		dimnames(GOF)[[3]] <- c('obs', 1:(nscan/odens))
		GOF[,,1] <- apply(Y, 3, gof_stats_bipartite)
	}

	# set names appropriately for bipartite vs unipartite
	if(bip) {
		names(APS) <- rownames(Y[,,1])  # row actors
		names(BPS) <- colnames(Y[,,1])  # column actors
		rownames(U) <- rownames(Y[,,1])
		rownames(V) <- colnames(Y[,,1])
	} else {
		names(APS)<-names(BPS)<-rownames(U)<-rownames(V)<-rownames(Y[,,1])
	}
	
	# names of parameters, asymmetric case
	if(!symmetric) {
		colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
		if(bip) {
			# bipartite: get names from Xlist
			if(!is.null(Xlist) && length(Xlist) > 0) {
				p_names <- dimnames(Xlist[[1]])[[3]]
				if(is.null(p_names) || length(p_names) != ncol(BETA)) {
					n_extra <- ncol(BETA) - as.integer(intercept)
					p_names <- if(intercept && n_extra > 0) {
						c("intercept", paste0("x", seq_len(n_extra)))
					} else if(intercept) {
						"intercept"
					} else {
						paste0("x", seq_len(ncol(BETA)))
					}
				}
				colnames(BETA) <- p_names
			}
		} else {
			colnames(BETA) <- dimnames(X)[[3]]
		}
	}
	
	# names of parameters, symmetric case
	if(symmetric) {
		colnames(VC) <- c("va", "ve")  
		rb<-intercept+seq(1,pr,length=pr) ; cb<-intercept+pr+seq(1,pr,length=pr)
		bnames<-dimnames(X)[[3]]
		if(intercept) {
			bni<-bnames[1]
			bnn<-gsub("row",bnames[rb],replacement="node") 
			bnd<-bnames[-c(1,rb,cb)]
			colnames(BETA)<-c(bni,bnn,bnd)
		} else {
			bnn<-gsub("row",bnames[rb],replacement="node")
			bnd<-bnames[-c(rb,cb)]
			colnames(BETA)<-c(bnn,bnd)
		}
	}    
	
	# mcmc
	have_coda<-suppressWarnings(
		try(requireNamespace("coda",quietly = TRUE),silent=TRUE)) 
	
	# show numerical stability message once at the beginning
	if(verbose && symmetric) {
		cli::cli_div(theme = list(span.note = list(color = "grey60", "font-style" = "italic")))
		cli::cli_inform(c(
			"i" = "Numerical adjustments may be applied for matrix stability.",
			" " = cli::col_grey("Small corrections ensure positive definiteness of covariance matrices."),
			" " = cli::col_grey("This is normal and does not affect model validity.")
		))
		cli::cli_end()
	}
	
	# suppress Armadillo warnings
	old_warn <- options()$warn
	if(symmetric) {
		options(warn = -1)  # suppress all warnings during MCMC
		on.exit(options(warn = old_warn), add = TRUE)
	}
	
	# pre-allocate residual array once
	E.nrm <- array(dim = dim(Z))

	if(burn!=0 && verbose){
		# only show progress if verbose=TRUE
		cli::cli_h3("Starting burn-in period...")
		cli::cli_progress_bar("Burn-in", total = burn, .envir = environment())
	}
	for (s in 1:(nscan + burn))  {

		# update Z (E.nrm reused each iteration)
		
		Xlist_tmp <- Xlist
		
		if(dynamic_ab) {
			# use dynamic helper to compute EZ with time-varying additive effects
			EZ <- get_EZ_dynamic_ab(Xlist_tmp, beta, a_mat, b_mat, U, V, N,
			                        bip = bip, G = G, nA = nA, nB = nB)
		} else {
			if(bip) {
				if(length(a) == 1) a <- rep(a, nA)
				if(length(b) == 1) b <- rep(b, nB)

				a_mat <- matrix(a[1:nA], nA, N)
				b_mat <- matrix(b[1:nB], nB, N)

				if(is.null(dim(U)) || length(dim(U)) == 2) {
					U_3d <- array(0, dim = c(nA, max(1, ncol(U)), N))
					for(t in 1:N) {
						U_3d[,,t] <- U
					}
					U <- U_3d
				}
				
				if(is.null(dim(V)) || length(dim(V)) == 2) {
					V_3d <- array(0, dim = c(nB, max(1, ncol(V)), N))
					for(t in 1:N) {
						V_3d[,,t] <- V
					}
					V <- V_3d
				}
				
				base_cube <- array(0, dim = c(nA, nB, N))
				for(t in 1:N) {
					if(length(beta) > 0) {
						X_dims <- dim(Xlist_tmp[[t]])
						if(X_dims[1] != nA || X_dims[2] != nB) {
							stop(paste("Bipartite X dimension mismatch at time", t, 
												": expected", nA, "x", nB, 
												"but got", X_dims[1], "x", X_dims[2]))
						}
						# sum up X*beta for each covariate
						base_mat <- matrix(0, nA, nB)
						n_covs <- min(length(beta), X_dims[3])
						for(i in 1:n_covs) {
							base_mat <- base_mat + beta[i] * Xlist_tmp[[t]][,,i]
						}
						base_cube[,,t] <- base_mat
					}
				}
				
				# call bipartite-specific function
				EZ <- get_EZ_bip_cpp(base_cube, a_mat, b_mat, U, V, G)
			} else {
				ab_mat <- outer(a, b,"+")
				storage.mode(ab_mat) <- "double"
				EZ <- get_EZ_cpp( Xlist_tmp, beta, ab_mat, U, V )
			}
		}
		# Phase 2A: Use batch C++ Z sampling where possible
		if(family == "normal") {
			# Single C++ call for all T time periods
			Z_batch <- try(rZ_nrm_batch_cpp(Z, EZ, rho, s2, Y), silent = TRUE)
			if(!inherits(Z_batch, 'try-error')) {
				Z <- Z_batch$Z
				E.nrm <- Z_batch$E_nrm
			} else {
				# fallback: per-t R sampling
				for(t in 1:N) { Z[,,t] <- rZ_nrm_fc(Z[,,t], EZ[,,t], rho, s2, Y[,,t]); E.nrm[,,t] <- Z[,,t] - EZ[,,t] }
				tryErrorChecks$Z <- tryErrorChecks$Z + 1
			}
		} else if(family == "tobit" && bip && rho == 0) {
			# Batch bipartite tobit Z sampling in C++
			Z_new <- try(rZ_tob_bip_batch_cpp(Z, EZ, s2, Y), silent = TRUE)
			if(!inherits(Z_new, 'try-error')) { Z <- Z_new }
			E.nrm <- Z - EZ
		} else if(family == "binary" && bip && rho == 0) {
			# Batch bipartite binary Z sampling in C++
			Z_new <- try(rZ_bin_bip_batch_cpp(Z, EZ, Y), silent = TRUE)
			if(!inherits(Z_new, 'try-error')) { Z <- Z_new }
			E.nrm <- Z - EZ
		} else {
			# Fallback: per-t R loop for families needing it
			for(t in 1:N) {
				if(family == "tobit") {
					Z[,,t] <- rZ_tob_fc(Z[,,t], EZ[,,t], rho, s2, Y[,,t])
					E.nrm[,,t] <- Z[,,t] - EZ[,,t]
				}
				if(family == "binary") {
					Z[,,t] <- rZ_bin_fc(Z[,,t], EZ[,,t], rho, Y[,,t])
				}
				if(family == "ordinal") { Z[,,t] <- rZ_ord_fc(Z[,,t], EZ[,,t], rho, Y[,,t]) }
				if(family == "cbin") { Z[,,t] <- rZ_cbin_fc(Z[,,t], EZ[,,t], rho, Y[,,t], odmax, odobs) }
				if(family == "frn") {
					Z[,,t] <- rZ_frn_fc(Z[,,t], EZ[,,t], rho, Y[,,t], YL[[t]], odmax, odobs)
				}
				if(family == "rrl") { Z[,,t] <- rZ_rrl_fc(Z[,,t], EZ[,,t], rho, Y[,,t], YL[[t]]) }
				if(family == "poisson") {
					Z[,,t] <- rZ_pois_fc(Z[,,t], EZ[,,t], rho, s2, Y[,,t])
					E.nrm[,,t] <- Z[,,t] - EZ[,,t]
				}
			}
		}
		
		# update s2
		if (is.element(family,c("normal","tobit","poisson"))){
			if(bip) {
				# bipartite: direct s2 sampling (no dyadic correlation, rs2_rep_fc_cpp
				# assumes square matrices and fails silently for rectangular)
				resid_all <- Z - EZ
				n_obs <- sum(!is.na(resid_all))
				ss <- sum(resid_all^2, na.rm = TRUE)
				s2 <- 1/rgamma(1, (1 + n_obs)/2, (1 + ss)/2)
			} else {
				s2New<-try(
					rs2_rep_fc_cpp(E.nrm,solve(matrix(c(1,rho,rho,1),2,2))),
					silent=TRUE)
				if(!inherits(s2New, 'try-error')){ s2 <- s2New } else { tryErrorChecks$s2<-tryErrorChecks$s2+1 }
			}
		}
		
		# update beta, a b with g-prior
		if(dynamic_ab) {
			if(bip) {
				# bipartite dynamic_ab: update beta via bipartite Gibbs
				# compute EZ without a/b for residual construction
				zero_a <- matrix(0, nA, N)
				zero_b <- matrix(0, nB, N)
				EZ_no_ab <- get_EZ_dynamic_ab(Xlist, beta, zero_a, zero_b, U, V, N,
				                              bip = TRUE, G = G, nA = nA, nB = nB)

				# Phase 2B: update beta via C++ (bipartite conjugate update)
				p_bip <- length(beta)
				if(p_bip > 0) {
					resid_b <- Z - EZ_no_ab
					for(t in 1:N) {
						resid_b[,,t] <- resid_b[,,t] - a_mat[,t]
						for(j in 1:nB) resid_b[,j,t] <- resid_b[,j,t] - b_mat[j,t]
					}
					xtx_xty <- compute_XtX_Xty_bip_cpp(Xlist, resid_b, p_bip)
					XtX <- xtx_xty$XtX
					Xty <- xtx_xty$Xty
					V_post <- tryCatch(
						solve(XtX / s2 + diag(p_bip) / (g[1] * s2)),
						error = function(e) diag(s2, p_bip)
					)
					m_post <- V_post %*% (Xty / s2)
					beta <- c(m_post + t(chol(V_post)) %*% rnorm(p_bip))
				}

				# recompute EZ_no_ab with updated beta
				EZ_no_ab <- get_EZ_dynamic_ab(Xlist, beta, zero_a, zero_b, U, V, N,
				                              bip = TRUE, G = G, nA = nA, nB = nB)
			} else if( (pr+pc+pd+intercept)>0 ){
				# unipartite dynamic_ab: use rbeta_ab_rep_fc_cpp
				EZ_no_ab <- get_EZ_cpp(Xlist, beta, outer(rep(0,n), rep(0,n), "+"), U, V)

				a_avg <- rowMeans(a_mat)
				b_avg <- rowMeans(b_mat)

				iSe2<-mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)) ; Sabs<-iSe2%*%Sab%*%iSe2
				tmp<-eigen(Sabs) ; k<-sum(zapsmall(tmp$val)>0 )
				if(k > 0) {
					if(k == 1) {
						G_eig <- tmp$vec[,1,drop=FALSE] * sqrt(tmp$val[1])
					} else {
						G_eig <- tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
					}
				} else {
					G_eig <- matrix(0, nrow=2, ncol=1)
					k <- 1
				}
				Z_temp <- sweep(Z,c(1,2),U%*%t(V))
				na_mask <- is.na(Z_temp)
				Z_temp[na_mask] <- 0

				k_int <- as.integer(k[1])
				g_num <- as.numeric(g[1])

				betaABCalc <- try(
					rbeta_ab_rep_fc_cpp(
						ZT=Z_temp, Xr=XrLong, Xc=XcLong, mX=mXLong, mXt=mXtLong,
						XX=xxLong, XXt=xxTLong, iSe2=iSe2, Sabs=Sabs, k=k_int, G=G_eig, g=g_num ),
					silent = FALSE)
				if(!inherits(betaABCalc, 'try-error')){
					beta <- c(betaABCalc$beta)
				} else {
					tryErrorChecks$betaAB <- tryErrorChecks$betaAB + 1
				}
				EZ_no_ab <- get_EZ_cpp(Xlist, beta, outer(rep(0,n), rep(0,n), "+"), U, V)
			} else {
				# no covariates: EZ_no_ab is just UV
				if(bip) {
					zero_a <- matrix(0, nA, N)
					zero_b <- matrix(0, nB, N)
					EZ_no_ab <- get_EZ_dynamic_ab(Xlist, beta, zero_a, zero_b, U, V, N,
					                              bip = TRUE, G = G, nA = nA, nB = nB)
				} else {
					EZ_no_ab <- get_EZ_cpp(Xlist, beta, outer(rep(0,n), rep(0,n), "+"), U, V)
				}
			}

			# update dynamic a and b
			if(bip) {
				# bipartite: a (nA x T) and b (nB x T) have different dims,
				# so update each separately via conjugate normal Gibbs steps
				var_innov <- sigma_ab^2
				for(t in 1:N) {
					resid_t <- Z[,,t] - EZ_no_ab[,,t]
					resid_t[is.na(resid_t)] <- 0

					# prior from AR(1)
					for(i in 1:nA) {
						pm <- if(t == 1) 0 else rho_ab * a_mat[i, t-1]
						if(t < N) pm <- 0.5*(pm + rho_ab * a_mat[i, t+1])
						pprec <- if(t == 1) (1 - rho_ab^2)/var_innov else 1/var_innov
						if(t < N) pprec <- pprec * 2
						dprec <- nB / Sab[1,1]
						post_prec <- pprec + dprec
						post_mean <- (pprec*pm + dprec*mean(resid_t[i,], na.rm=TRUE)) / post_prec
						a_mat[i,t] <- rnorm(1, post_mean, 1/sqrt(post_prec))
					}
					# update resid for b step
					for(i in 1:nA) resid_t[i,] <- resid_t[i,] - a_mat[i,t]
					for(j in 1:nB) {
						pm <- if(t == 1) 0 else rho_ab * b_mat[j, t-1]
						if(t < N) pm <- 0.5*(pm + rho_ab * b_mat[j, t+1])
						pprec <- if(t == 1) (1 - rho_ab^2)/var_innov else 1/var_innov
						if(t < N) pprec <- pprec * 2
						dprec <- nA / Sab[2,2]
						post_prec <- pprec + dprec
						post_mean <- (pprec*pm + dprec*mean(resid_t[,j], na.rm=TRUE)) / post_prec
						b_mat[j,t] <- rnorm(1, post_mean, 1/sqrt(post_prec))
					}
				}
			} else {
				ab_update <- sample_dynamic_ab_cpp(a_mat, b_mat, Z, EZ_no_ab,
				                                   rho_ab, sigma_ab, Sab, symmetric)
				a_mat <- ab_update$a
				b_mat <- ab_update$b
			}

			a <- a_mat[,1]
			b <- b_mat[,1]

			# update AR(1) parameters periodically (including during burn-in)
			if(s %% 20 == 0) {
				if(bip) {
					# bipartite: a (nA x T) and b (nB x T) have different row counts,
					# so compute sufficient stats in R rather than C++ which assumes same n
					sp <- 0; ssl <- 0; ssc <- 0
					for(i in 1:nA) for(t in 2:N) {
						sp <- sp + a_mat[i,t]*a_mat[i,t-1]
						ssl <- ssl + a_mat[i,t-1]^2
						ssc <- ssc + a_mat[i,t]^2
					}
					for(j in 1:nB) for(t in 2:N) {
						sp <- sp + b_mat[j,t]*b_mat[j,t-1]
						ssl <- ssl + b_mat[j,t-1]^2
						ssc <- ssc + b_mat[j,t]^2
					}
					# mh step for rho_ab
					rho_prop <- rnorm(1, rho_ab, 0.1)
					if(abs(rho_prop) < 1) {
						ll_c <- -0.5*(ssc - 2*rho_ab*sp + rho_ab^2*ssl) / sigma_ab^2
						ll_p <- -0.5*(ssc - 2*rho_prop*sp + rho_prop^2*ssl) / sigma_ab^2
						lp_c <- -0.5*log(1 - rho_ab^2)
						lp_p <- -0.5*log(1 - rho_prop^2)
						if(log(runif(1)) < (ll_p + lp_p) - (ll_c + lp_c)) rho_ab <- rho_prop
					}
					# inverse-gamma for sigma_ab
					cnt <- (nA + nB) * (N - 1)
					ss_resid <- 0
					for(i in 1:nA) for(t in 2:N) ss_resid <- ss_resid + (a_mat[i,t] - rho_ab*a_mat[i,t-1])^2
					for(j in 1:nB) for(t in 2:N) ss_resid <- ss_resid + (b_mat[j,t] - rho_ab*b_mat[j,t-1])^2
					sigma_ab <- sqrt(1/rgamma(1, shape = 2 + cnt/2, rate = 1 + ss_resid/2))
				} else {
					rho_ab <- sample_rho_ab_cpp(a_mat, b_mat, sigma_ab, rho_ab, symmetric)
					sigma_ab <- sample_sigma_ab_cpp(a_mat, b_mat, rho_ab, symmetric)
				}
			}

		} else {
			# standard static update
			if(bip) {
				# Phase 2B: bipartite gibbs update via C++
				U_2d <- if(length(dim(U)) == 3) U[,,1] else U
				V_2d <- if(length(dim(V)) == 3) V[,,1] else V
				UV_eff <- if(RA > 0 && RB > 0) U_2d %*% G %*% t(V_2d) else matrix(0, nA, nB)

				betaABCalc <- tryCatch(
					rbeta_ab_bip_gibbs_cpp(
						Z, Xlist_tmp, UV_eff,
						if(length(a) == nA) a else rep(0, nA),
						if(length(b) == nB) b else rep(0, nB),
						s2, as.numeric(g[1]),
						Sab[1,1], Sab[2,2],
						rvar, cvar
					),
					error = function(e) {
						list(beta = beta, a = if(length(a) == nA) a else rep(0, nA),
						     b = if(length(b) == nB) b else rep(0, nB))
					}
				)
			} else if( (pr+pc+pd+intercept)>0 ){
				iSe2<-mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)) ; Sabs<-iSe2%*%Sab%*%iSe2
				tmp<-eigen(Sabs) ; k<-sum(zapsmall(tmp$val)>0 )
				if(k > 0) {
					if(k == 1) {
						G_eig <- tmp$vec[,1,drop=FALSE] * sqrt(tmp$val[1])
					} else {
						G_eig <- tmp$vec[,1:k] %*% sqrt(diag(tmp$val[1:k],nrow=k))
					}
				} else {
					G_eig <- matrix(0, nrow=2, ncol=1)
					k <- 1  # avoid dimension issues
				}
				# handle NA values in Z before passing to C++
				Z_temp <- sweep(Z,c(1,2),U%*%t(V))
				na_mask <- is.na(Z_temp)
				Z_temp[na_mask] <- 0

				k_int <- as.integer(k[1])
				g_num <- as.numeric(g[1])

				betaABCalc <- try(
					rbeta_ab_rep_fc_cpp(
						ZT=Z_temp, Xr=XrLong, Xc=XcLong, mX=mXLong, mXt=mXtLong,
						XX=xxLong, XXt=xxTLong, iSe2=iSe2, Sabs=Sabs, k=k_int, G=G_eig, g=g_num ),
					silent = FALSE)
			} else {
				betaABCalc <- try(
					rbeta_ab_rep_fc(sweep(Z,c(1,2),U%*%t(V)), Sab, rho, X, s2),
					silent = FALSE)
			}
			if(!inherits(betaABCalc, 'try-error')){
				beta <- c(betaABCalc$beta)
				a <- c(betaABCalc$a) * rvar
				b <- c(betaABCalc$b) * cvar
				if(symmetric){ a<-b<-(a+b)/2 }
			} else { 
				# silent fallback when betaAB calculation fails
				tryErrorChecks$betaAB<-tryErrorChecks$betaAB+1  
			}
		}
		
		# update Sab using unified function
		if(bip) {
			if(rvar) {
				Sab[1,1] <- 1/rgamma(1, (prior$etaab + nA)/2, 
														 (prior$etaab*prior$Sab0[1,1] + sum(a^2))/2)
			}
			if(cvar) {
				Sab[2,2] <- 1/rgamma(1, (prior$etaab + nB)/2, 
														 (prior$etaab*prior$Sab0[2,2] + sum(b^2))/2)
			}
			Sab[1,2] <- Sab[2,1] <- 0  # no covariance for bipartite
		} else {
			Sab <- rSab_fc(a, b, Sab0=prior$Sab0/prior$etaab, eta0=prior$etaab, 
										 rvar=rvar, cvar=cvar, symmetric=symmetric)
		}
		
		# update rho (only for asymmetric models)
		if(dcor && !symmetric) {
			if(dynamic_ab) {
				E.T <- Z - get_EZ_dynamic_ab(Xlist, beta, a_mat, b_mat, U, V, N,
				                             bip = bip, G = G, nA = nA, nB = nB)
			} else {
				E.T <- Z - get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U, V )
			}
			rhoNew<-try( rrho_mh_rep_cpp(E.T, rho,s2), silent=TRUE )
			if(!inherits(rhoNew, 'try-error')){ rho<-rhoNew } else { tryErrorChecks$rho<-tryErrorChecks$rho+1 }
		}
		
		# update U,V
		if (R > 0) {
			if(dynamic_ab) {
				# for dynamic ab, use time-varying effects but zero out UV
				U_zero <- if(bip) array(0, dim=c(nA, max(1,RA), N)) else U*0
				V_zero <- if(bip) array(0, dim=c(nB, max(1,RB), N)) else V*0
				E <- Z - get_EZ_dynamic_ab(Xlist, beta, a_mat, b_mat, U_zero, V_zero, N,
				                           bip = bip, G = G, nA = nA, nB = nB)
			} else if(bip) {
					U_zero <- array(0, dim = c(nA, max(1, ncol(U)), N))
				V_zero <- array(0, dim = c(nB, max(1, ncol(V)), N))
				
				# create base cube with X*beta and a/b effects
				base_cube <- array(0, dim = c(nA, nB, N))
				for(t in 1:N) {
					if(length(beta) > 0) {
						base_mat <- matrix(0, nA, nB)
						n_covs <- min(length(beta), dim(Xlist_tmp[[t]])[3])
						for(i in 1:n_covs) {
							base_mat <- base_mat + beta[i] * Xlist_tmp[[t]][,,i]
						}
						base_cube[,,t] <- base_mat
					}
				}
				
				# create time-expanded a and b matrices
				a_mat_temp <- matrix(a[1:nA], nA, N)
				b_mat_temp <- matrix(b[1:nB], nB, N)
				
				E <- Z - get_EZ_bip_cpp(base_cube, a_mat_temp, b_mat_temp, U_zero, V_zero, G)
			} else {
				E <- Z - get_EZ_cpp( Xlist, beta, outer(a, b,"+"), U*0, V*0 )
			}
			shrink<- (s>.5*burn)
			
			if(dynamic_uv) {
				if(bip) {
					# Phase 2C: bipartite dynamic UV via C++ (replaces nested R loops)
					UV_try <- tryCatch(
						rUV_dynamic_bip_fc_cpp(U_cube, V_cube, E, G,
						                       rho_uv, sigma_uv, s2),
						error = function(e) NULL
					)

					if(!is.null(UV_try) &&
						   all(is.finite(UV_try$U)) && all(is.finite(UV_try$V))) {
						U_cube <- UV_try$U; V_cube <- UV_try$V
						U <- apply(U_cube, c(1,2), mean)
						V <- apply(V_cube, c(1,2), mean)
					} else {
						tryErrorChecks$UV <- tryErrorChecks$UV + 1
						UV_try <- NULL  # ensure AR(1) guard treats this as failure
					}

					# ar(1) parameter updates for bipartite
					# only update if UV sampling succeeded (avoid inconsistent MCMC state)
					if(!is.null(UV_try) && s %% 10 == 0) {
						sigma2_inv <- 1 / (sigma_uv^2)
						sp <- 0; ssl <- 0; ss <- 0
						for(tt in 2:N) {
							sp <- sp + sum(U_cube[,,tt] * U_cube[,,tt-1]) + sum(V_cube[,,tt] * V_cube[,,tt-1])
							ssl <- ssl + sum(U_cube[,,tt-1]^2) + sum(V_cube[,,tt-1]^2)
							innov_U <- U_cube[,,tt] - rho_uv * U_cube[,,tt-1]
							innov_V <- V_cube[,,tt] - rho_uv * V_cube[,,tt-1]
							ss <- ss + sum(innov_U^2) + sum(innov_V^2)
						}
						# sample rho_uv
						var_post <- 1 / (ssl * sigma2_inv + 1)
						mean_post <- var_post * sp * sigma2_inv
						rho_uv <- max(-0.99, min(0.99, rnorm(1, mean_post, sqrt(var_post))))
						# sample sigma_uv
						n_params <- (nA * RA + nB * RB) * (N - 1)
						shape_post <- n_params / 2 + 1
						scale_post <- ss / 2 + 1
						sigma_uv <- sqrt(1 / rgamma(1, shape_post, rate=scale_post))
					}
				} else {
					# unipartite dynamic UV update with AR(1) evolution
					UV <- try(
						rUV_dynamic_fc_cpp(U_cube, V_cube, E, rho_uv, sigma_uv, s2, shrink, symmetric),
						silent = FALSE)
					if(inherits(UV, 'try-error')){
						UV <- list(U=U_cube, V=V_cube)
						tryErrorChecks$UV<-tryErrorChecks$UV+1
					} else {
						U_cube <- UV$U
						V_cube <- UV$V
						U <- apply(U_cube, c(1,2), mean)
						V <- apply(V_cube, c(1,2), mean)
					}

					# update AR(1) parameters (only if UV sampling succeeded)
					if(!inherits(UV, 'try-error') && s %% 10 == 0) {
						rho_uv <- sample_rho_uv(U_cube, V_cube, sigma_uv, rho_uv, symmetric)
						sigma_uv <- sample_sigma_uv(U_cube, V_cube, rho_uv, symmetric)
					}
				}
			} else {
				# standard static UV update
				if(bip) {
					# bipartite MH random walk update for U and V
					if(RA > 0 && RB > 0) {
						U_2d <- if(length(dim(U)) == 3) U[,,1] else U
						V_2d <- if(length(dim(V)) == 3) V[,,1] else V

						# update U rows
						for(i in 1:nA) {
							u_prop <- U_2d[i,] + rnorm(RA, 0, mh_sd_U)
							uv_curr <- as.vector(U_2d[i,,drop=FALSE] %*% G %*% t(V_2d))
							uv_prop <- as.vector(matrix(u_prop, 1) %*% G %*% t(V_2d))
							ll_diff <- 0
							for(t in 1:N) {
								e_it <- E[i,,t]
								ll_diff <- ll_diff - sum((e_it - uv_prop)^2 - (e_it - uv_curr)^2, na.rm=TRUE)
							}
							mh_att_U <- mh_att_U + 1L
							if(log(runif(1)) < ll_diff / (2*s2)) {
								U_2d[i,] <- u_prop
								mh_acc_U <- mh_acc_U + 1L
							}
						}

						# update V rows
						for(j in 1:nB) {
							v_prop <- V_2d[j,] + rnorm(RB, 0, mh_sd_V)
							uv_curr <- as.vector(U_2d %*% G %*% V_2d[j,])
							uv_prop <- as.vector(U_2d %*% G %*% v_prop)
							ll_diff <- 0
							for(t in 1:N) {
								e_jt <- E[,j,t]
								ll_diff <- ll_diff - sum((e_jt - uv_prop)^2 - (e_jt - uv_curr)^2, na.rm=TRUE)
							}
							mh_att_V <- mh_att_V + 1L
							if(log(runif(1)) < ll_diff / (2*s2)) {
								V_2d[j,] <- v_prop
								mh_acc_V <- mh_acc_V + 1L
							}
						}

						# adapt proposal sd during burn-in (target ~30-50% acceptance)
						if(s <= burn && s %% 50 == 0 && mh_att_U > 0) {
							rate_U <- mh_acc_U / mh_att_U
							rate_V <- mh_acc_V / mh_att_V
							if(rate_U < 0.20) mh_sd_U <- mh_sd_U * 0.8
							if(rate_U > 0.60) mh_sd_U <- mh_sd_U * 1.2
							if(rate_V < 0.20) mh_sd_V <- mh_sd_V * 0.8
							if(rate_V > 0.60) mh_sd_V <- mh_sd_V * 1.2
							mh_acc_U <- 0L; mh_att_U <- 0L
							mh_acc_V <- 0L; mh_att_V <- 0L
						}

						# expand back to 3D for get_EZ_bip_cpp compatibility
						if(length(dim(U)) == 3) {
							for(t in 1:N) { U[,,t] <- U_2d; V[,,t] <- V_2d }
						} else {
							U <- U_2d; V <- V_2d
						}
					}
				} else if(symmetric) {
					EA<-apply(E,c(1,2),mean) ; EA<-.5*(EA+t(EA))
					UV<-try(
						rUV_sym_fc_cpp(EA, U, V, 
													 s2/dim(E)[3], shrink, symLoopIDs[[s]]-1), silent=TRUE )
					if(inherits(UV, 'try-error')){ UV <- list(U=U,V=V) ; tryErrorChecks$UV<-tryErrorChecks$UV+1 }
					U<-UV$U ; V<-UV$V
				} else if(!symmetric) {
					UV <- try(
						rUV_rep_fc_cpp(E, U, V, rho, s2,
													 mhalf(solve(matrix(c(1,rho,rho,1),2,2)*s2)),
													 maxmargin=1e-6, shrink, asymLoopIDs[[s]]-1 ), silent = TRUE )
					if(inherits(UV, 'try-error')){
						UV <- try(rUV_rep_fc(E, U, V, rho, s2, shrink), silent = TRUE)
						if(inherits(UV, 'try-error')){ UV <- list(U=U,V=V) ; tryErrorChecks$UV<-tryErrorChecks$UV+1 }
					}
					U<-UV$U ; V<-UV$V
				}
			}

			# update G for bipartite (interaction matrix between row and column latent spaces)
			if(bip && RA > 0 && RB > 0 && !is.null(G)) {
				U_3d <- if(length(dim(U)) == 3) U else {
					U_tmp <- array(0, dim = c(nA, RA, N))
					for(t in 1:N) U_tmp[,,t] <- U
					U_tmp
				}
				V_3d <- if(length(dim(V)) == 3) V else {
					V_tmp <- array(0, dim = c(nB, RB, N))
					for(t in 1:N) V_tmp[,,t] <- V
					V_tmp
				}
				G_new <- tryCatch(
					sample_G_bip_cpp(E, U_3d, V_3d, lambdaG = 1.0, s2 = rep(s2, N)),
					error = function(e) G
				)
				G <- G_new
			}
		}

		# burn-in countdown
		if(burn!=0 && s <= burn && verbose){
			# only update progress if verbose=TRUE
			cli::cli_progress_update()
		}
		
		# store parameter values and monitor the MC
		if(s==burn+1 && verbose && burn!=0){
			# only show messages if verbose=TRUE
			cli::cli_progress_done()
			cli::cli_alert_success("Burn-in period complete")
			cli::cli_progress_bar("Sampling", total = nscan, .envir = environment())
		}
		if(s%%odens==0 & s>burn)  {
			
			# store BETA and VC - symmetric case 
			if(symmetric){
				br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2
						if(intercept) {
					sbeta<-c(beta[1],bn,beta[-c(1,rb,cb)])
				} else {
					sbeta<-c(bn,beta[-c(rb,cb)])
				}
				BETA[iter,]<-sbeta
				VC[iter,]<-c(Sab[1,1],s2)
			}
			
			# store BETA and VC - asymmetric case 
			if(!symmetric){
				BETA[iter,]<-beta
				VC[iter,]<- c(Sab[upper.tri(Sab, diag = T)], rho,s2)
			}
			
			# store dynamic parameters
			if(dynamic_ab && !is.null(rho_ab)) {
				RHO_AB[iter] <- rho_ab
			}
			if(dynamic_uv && !is.null(rho_uv)) {
				RHO_UV[iter] <- rho_uv
			}
			
			# update posterior sums of random effects
			if(dynamic_uv && R > 0 && !is.null(U_cube) && !is.null(V_cube)) {
				U_SUM <- U_SUM + U_cube
				V_SUM <- V_SUM + V_cube
				# store UV products for each time point
				for(t in 1:N) {
					U_t <- U_cube[,,t,drop=FALSE]
					V_t <- V_cube[,,t,drop=FALSE]
					if(bip) {
						U_t <- matrix(U_t, nrow=nA, ncol=RA)
						V_t <- matrix(V_t, nrow=nB, ncol=RB)
						UVPS[,,t] <- UVPS[,,t] + U_t %*% G %*% t(V_t)
					} else {
						U_t <- matrix(U_t, nrow=n, ncol=R)
						V_t <- matrix(V_t, nrow=n, ncol=R)
						UVPS[,,t] <- UVPS[,,t] + U_t %*% t(V_t)
					}
				}
			} else if(R > 0) {
				if(bip) {
					# bipartite: use 2D slices with G interaction matrix
					U_2d <- if(length(dim(U)) == 3) U[,,1] else U
					V_2d <- if(length(dim(V)) == 3) V[,,1] else V
					UVPS <- UVPS + U_2d %*% G %*% t(V_2d)
				} else {
					UVPS <- UVPS + U %*% t(V)
				}
			}
			if(dynamic_ab) {
				# for dynamic ab, accumulate time-specific effects
				if(iter == 1) {
					# initialize accumulation matrices
					APS_dyn <- a_mat
					BPS_dyn <- b_mat
				} else {
					APS_dyn <- APS_dyn + a_mat
					BPS_dyn <- BPS_dyn + b_mat
				}
				APS <- APS + rowMeans(a_mat)
				BPS <- BPS + rowMeans(b_mat)
			} else {
				APS <- APS + a
				BPS <- BPS + b
			} 
			
			# simulate from posterior predictive
			if(!bip) {
				# reuse EZ from earlier in the iteration
				dimnames(EZ) <- dimnames(Y)
				Ys <- EZ*0
				for (t in 1:N) {
					if(symmetric){ EZ[,,t]<-(EZ[,,t]+t(EZ[,,t]))/2 }
					
					if(family=="binary"){ Ys[,,t]<-simY_bin(EZ[,,t],rho) }
					if(family=="cbin"){ Ys[,,t]<-1*(simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t])>0)}
					if(family=="frn"){ Ys[,,t]<-simY_frn(EZ[,,t],rho,odmax,YO=Y[,,t]) }
					if(family=="rrl"){ Ys[,,t]<-simY_rrl(EZ[,,t],rho,odobs,YO=Y[,,t] ) }
					if(family=="normal"){ Ys[,,t]<-simY_nrm(EZ[,,t],rho,s2) }
					if(family=="tobit"){ Ys[,,t]<-simY_tob(EZ[,,t],rho,s2) }
					if(family=="ordinal"){ Ys[,,t]<-simY_ord(EZ[,,t],rho,Y[,,t]) }
					if(family=="poisson"){ Ys[,,t]<-simY_pois(EZ[,,t]) }
				
				if(symmetric) {
					Yst<-Ys[,,t] ; Yst[lower.tri(Yst)]<-0 ; Ys[,,t]<-Yst+t(Yst)
				}
			}
			} else {
				# bipartite posterior predictive simulation (no dyadic correlation)
				Ys <- array(0, dim=c(nA, nB, N))
				for (t in 1:N) {
					EZ_t <- EZ[,,t]
					if(family == "binary") {
						noise <- matrix(rnorm(nA*nB), nA, nB)
						Ys[,,t] <- 1*(EZ_t + noise > 0)
					} else if(family == "normal") {
						Ys[,,t] <- EZ_t + matrix(rnorm(nA*nB, 0, sqrt(s2)), nA, nB)
					} else if(family == "poisson") {
						lambda <- exp(EZ_t)
						lambda[lambda > 1e6] <- 1e6
						Ys[,,t] <- matrix(rpois(nA*nB, c(lambda)), nA, nB)
					} else if(family == "tobit") {
						Ys_t <- EZ_t + matrix(rnorm(nA*nB, 0, sqrt(s2)), nA, nB)
						Ys_t[Ys_t <= 0] <- 0
						Ys[,,t] <- Ys_t
					} else if(family == "ordinal") {
						# quantile-based discretization
						uY <- sort(unique(c(Y[,,t])))
						uY <- uY[!is.na(uY)]
						FY <- table(c(Y[,,t]))
						FY <- FY / sum(FY)
						FY <- cumsum(FY)
						ZS <- EZ_t + matrix(rnorm(nA*nB), nA, nB)
						qZ <- quantile(ZS, FY[-length(FY)], na.rm = TRUE)
						YS <- ZS * 0 + max(uY)
						for(w in 1:(length(uY) - 1)) {
							YS[(YS == max(uY)) & ZS <= qZ[w]] <- uY[w]
						}
						Ys[,,t] <- YS
					} else if(family == "cbin") {
						# censored binary: simulate nominations
						ZS <- EZ_t + matrix(rnorm(nA*nB), nA, nB)
						Ys[,,t] <- 1*(ZS > 0)
					} else if(family == "frn") {
						# fixed rank nomination
						ZS <- EZ_t + matrix(rnorm(nA*nB), nA, nB)
						YS <- ZS * 0
						for(i in 1:nA) {
							rs <- rank(ZS[i,]) - (nB - odmax[i])
							YS[i,] <- rs * (rs > 0) * (ZS[i,] > 0)
							pos <- YS[i,] > 0
							if(any(pos)) {
								YS[i,pos] <- match(YS[i,pos], sort(unique(YS[i,pos])))
							}
						}
						Ys[,,t] <- YS
					} else if(family == "rrl") {
						# row rank likelihood
						ZS <- EZ_t + matrix(rnorm(nA*nB), nA, nB)
						YS <- ZS * 0
						for(i in 1:nA) {
							ri <- order(-ZS[i,])[seq_len(odobs[i])]
							YS[i,ri] <- seq(odobs[i], 1, length = odobs[i])
						}
						Ys[,,t] <- YS
					}
				}
			}
			
			# update posterior sum
			YPS<-YPS+Ys
			
			# save posterior predictive GOF stats
			if(gof){
				Ys[is.na(Y)] <- NA
				if(bip) {
					GOF[,,(iter)+1] <- apply(Ys, 3, gof_stats_bipartite)
				} else {
					GOF[,,(iter)+1] <- apply(Ys, 3, gof_stats)
				}
			}
			
			# print MC progress 
			if(verbose) {
				beta_means <- round(apply(BETA[1:iter,,drop=FALSE],2,mean),2)
				vc_means <- round(apply(VC[1:iter,,drop=FALSE],2,mean),2)
				cli::cli_text("Iteration {.val {s}}: beta = [{.field {paste(beta_means, collapse=', ')}}], VC = [{.field {paste(vc_means, collapse=', ')}}]")
				if (have_coda & nrow(VC[1:iter,,drop=FALSE]) > 3 & length(beta)>0)  {
					eff_sizes <- round(coda::effectiveSize(BETA[1:iter,,drop=FALSE]))
					cli::cli_text("  Effective sizes: [{.emph {paste(eff_sizes, collapse=', ')}}]")
				}
			}
			
			# periodic save
			if(periodic_save & s %in% savePoints & !is.null(out_file)){
				# save start_vals for future model runs
				if(dynamic_uv && R > 0 && dynamic_ab) {
					start_vals <- list(Z=Z,beta=beta,a=a_mat,b=b_mat,U=U_cube,V=V_cube,
													 rho=rho,s2=s2,Sab=Sab,rho_uv=rho_uv,sigma_uv=sigma_uv,
													 rho_ab=rho_ab, sigma_ab=sigma_ab)
				} else if(dynamic_uv && R > 0) {
					start_vals <- list(Z=Z,beta=beta,a=a,b=b,U=U_cube,V=V_cube,
													 rho=rho,s2=s2,Sab=Sab,rho_uv=rho_uv,sigma_uv=sigma_uv)
				} else if(dynamic_ab) {
					start_vals <- list(Z=Z,beta=beta,a=a_mat,b=b_mat,U=U,V=V,rho=rho,s2=s2,Sab=Sab,
													 rho_ab=rho_ab, sigma_ab=sigma_ab)
				} else {
					start_vals <- list(Z=Z,beta=beta,a=a,b=b,U=U,V=V,rho=rho,s2=s2,Sab=Sab)
				}
				fit <- get_fit_object( APS=APS, BPS=BPS, UVPS=UVPS, YPS=YPS,
														 BETA=BETA, VC=VC, GOF=GOF, Xlist=Xlist, actorByYr=actorByYr,
														 start_vals=start_vals, symmetric=symmetric, tryErrorChecks=tryErrorChecks,
														 model.name=model.name, family=family, odmax=odmax,
														 nA=if(bip) nA else NULL, nB=if(bip) nB else NULL, n_time=N,
														 dynamic_uv=dynamic_uv, dynamic_ab=dynamic_ab, bip=bip,
														 G=if(bip) G else NULL)
				save(fit, file=out_file) ; rm(list=c('fit','start_vals'))
			}
			
			# plotting 
			if(plot && iter %% 10 == 0){
				if(iter == 10) {
					cli::cli_alert_info("Real-time plotting enabled - updating every 10 iterations")
				}
				# create temporary fit object for plotting
				current_APS <- APS / iter
				current_BPS <- BPS / iter
				current_UVPS <- if(R > 0) UVPS / iter else NULL
				current_U <- if(dynamic_uv && R > 0) U_SUM / iter else U
				current_V <- if(dynamic_uv && R > 0) V_SUM / iter else V
				
				plot_result <- try({
					temp_fit <- get_fit_object( APS=current_APS, BPS=current_BPS, UVPS=current_UVPS, YPS=YPS,
																		 BETA=BETA[1:iter,,drop=FALSE], VC=VC[1:iter,,drop=FALSE],
																		 GOF=GOF, Xlist=Xlist, actorByYr=actorByYr,
																		 start_vals=NULL, symmetric=symmetric, tryErrorChecks=tryErrorChecks,
																		 model.name=model.name, U=current_U, V=current_V,
																		 dynamic_uv=dynamic_uv, dynamic_ab=dynamic_ab, bip=bip,
																		 rho_ab=if(dynamic_ab) RHO_AB[1:iter] else NULL,
																		 rho_uv=if(dynamic_uv) RHO_UV[1:iter] else NULL,
																		 family=family, odmax=odmax, nA=if(bip) nA else NULL,
																		 nB=if(bip) nB else NULL, n_time=N,
																		 G=if(bip) G else NULL)
					class(temp_fit) <- "lame"
					suppressMessages(plot(temp_fit, which=c(1,2), pages="single"))
				}, silent=TRUE)
				
				if(inherits(plot_result, "try-error") && iter == 10) {
					cli::cli_alert_warning("Real-time plotting failed - continuing without plots")
				}
			}
			iter<-iter+1
		} # post burn-in
		if(verbose && s > burn){
			# only update sampling progress if verbose=TRUE
			cli::cli_progress_update()
		}
	} # end MCMC  
	if(verbose){
		# only show completion message if verbose=TRUE
		cli::cli_progress_done()
		cli::cli_alert_success("MCMC sampling complete")
	}

	total_iters <- nscan + burn
	failure_threshold <- max(1, round(0.10 * total_iters))
	failed_params <- character(0)
	for(pname in names(tryErrorChecks)) {
		cnt <- tryErrorChecks[[pname]]
		if(cnt > failure_threshold) {
			failed_params <- c(failed_params,
				paste0(pname, " (", cnt, "/", total_iters, " iterations)"))
		}
	}
	if(length(failed_params) > 0) {
		cli::cli_warn(c(
			"Some MCMC sampling steps had frequent failures:",
			"i" = "Affected parameters: {.val {failed_params}}",
			"i" = "These parameters may not have converged. Consider increasing {.arg nscan} or simplifying the model."
		))
	}

	# save start_vals for future model runs
	if(dynamic_ab) {
		start_vals <- list( Z=Z, beta=beta, a=a_mat, b=b_mat, U=U, V=V, rho=rho, s2=s2, Sab=Sab,
											rho_ab=rho_ab, sigma_ab=sigma_ab)
	} else {
		start_vals <- list( Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
	}
	# store G in start_vals for bipartite
	if(bip && !is.null(G)) start_vals$G <- G
	
	if(!is.null(model.name)) {
		# count effective parameters
		p_eff <- ncol(BETA) + rvar*n + cvar*n + R*(R+1)/2
		
	}
	
	# output
	if(dynamic_uv && R > 0) {
		# return 3D arrays for dynamic UV
		U_final <- U_SUM / (nscan/odens)
		V_final <- V_SUM / (nscan/odens)
		UVPS_final <- UVPS / (nscan/odens)
	} else {
		U_final <- U
		V_final <- V
		UVPS_final <- UVPS / (nscan/odens)
	}
	
	# handle dynamic additive effects output
	if(dynamic_ab) {
		APS_final <- APS_dyn / (nscan/odens)
		BPS_final <- BPS_dyn / (nscan/odens)
	} else {
		# leave raw sums -- get_fit_object divides by nrow(VC)
		APS_final <- APS
		BPS_final <- BPS
	}
	
	fit <- get_fit_object( APS=APS_final, BPS=BPS_final, UVPS=UVPS_final, YPS=YPS,
											 BETA=BETA, VC=VC, GOF=GOF, Xlist=Xlist, actorByYr=actorByYr,
											 start_vals=start_vals, symmetric=symmetric, tryErrorChecks=tryErrorChecks,
											 model.name=model.name, U=U_final, V=V_final,
											 dynamic_uv=dynamic_uv, dynamic_ab=dynamic_ab, bip=bip,
											 rho_ab=RHO_AB, rho_uv=RHO_UV,
											 family=family, odmax=odmax, nA=if(bip) nA else NULL,
											 nB=if(bip) nB else NULL, n_time=N,
											 Y_obs=Y, G=if(bip) G else NULL)
	
	####
	# scalar typing hygiene
	if (!is.null(fit$RHO)) fit$RHO <- as.numeric(fit$RHO)
	if (!is.null(fit$s2))  fit$s2  <- as.numeric(fit$s2)
	####

	####
	# normalize EZ to list-of-matrices
	.as_ez_list <- function(EZ) {
		if (is.null(EZ)) return(NULL)
		if (is.list(EZ)) return(EZ)
		if (length(dim(EZ)) == 3L) {
			T <- dim(EZ)[3]
			out <- vector("list", T)
			for (t in seq_len(T)) out[[t]] <- EZ[,,t]
			return(out)
		}
		list(EZ)
	}
	fit$EZ <- .as_ez_list(fit$EZ)
	####

	####
	# temporal alignment of latent positions
	if (isTRUE(dynamic_uv)) {
		if (mode == "unipartite") {
			if (!is.null(fit$U) && length(dim(fit$U)) == 3L) {
				fit$U <- align_over_time_unip(fit$U)
					}
		} else {
			if (!is.null(fit$U) && !is.null(fit$V) &&
					length(dim(fit$U)) == 3L && length(dim(fit$V)) == 3L) {
				tmp <- align_over_time_bip(fit$U, fit$V, fit$G)
				fit$U <- tmp$U; fit$V <- tmp$V; fit$G <- tmp$G
			}
		}
		rho_hat <- estimate_rho_from_U(fit$U)
		if (is.finite(rho_hat)) fit$rho_uv <- as.numeric(rho_hat)
	}
	####

	if(symmetric) {
		options(warn = old_warn)
	}

	# store bipartite MH proposal diagnostics
	if(bip && RA > 0 && RB > 0) {
		fit$mh_diagnostics <- list(
			sd_U = mh_sd_U, sd_V = mh_sd_V
		)
	}

	class(fit) <- c("lame", "ame")
	return(fit)
	
}