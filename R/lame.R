.lame_sample_rho_uv_drift_only <- function(U_cube, V_cube, sigma_uv, rho_current,
		symmetric, delta_u, delta_v = NULL, prior_mean = 0.9, prior_sd = 0.1) {
	n_time <- dim(U_cube)[3]
	if(is.null(n_time) || n_time < 2 || is.null(delta_u)) return(rho_current)
	sigma2_inv <- 1 / (sigma_uv * sigma_uv)
	prior_prec <- 1 / (prior_sd * prior_sd)
	sum_yt_yt1 <- 0
	sum_yt1_yt1 <- 0

	for(t in 2:n_time) {
		keep_u <- delta_u[, t] == 0
		if(any(keep_u)) {
			U_t <- U_cube[keep_u, , t, drop = FALSE]
			U_t1 <- U_cube[keep_u, , t - 1, drop = FALSE]
			sum_yt_yt1 <- sum_yt_yt1 + sum(U_t * U_t1)
			sum_yt1_yt1 <- sum_yt1_yt1 + sum(U_t1 * U_t1)
		}
		if(!symmetric && !is.null(delta_v)) {
			keep_v <- delta_v[, t] == 0
			if(any(keep_v)) {
				V_t <- V_cube[keep_v, , t, drop = FALSE]
				V_t1 <- V_cube[keep_v, , t - 1, drop = FALSE]
				sum_yt_yt1 <- sum_yt_yt1 + sum(V_t * V_t1)
				sum_yt1_yt1 <- sum_yt1_yt1 + sum(V_t1 * V_t1)
			}
		}
	}
	if(!is.finite(sum_yt1_yt1) || sum_yt1_yt1 <= 0) return(rho_current)
	var_post <- 1 / (sum_yt1_yt1 * sigma2_inv + prior_prec)
	mean_post <- var_post * (sum_yt_yt1 * sigma2_inv + prior_mean * prior_prec)
	pmax(-0.99, pmin(0.99, rnorm(1, mean_post, sqrt(var_post))))
}

.lame_sample_sigma_uv_drift_only <- function(U_cube, V_cube, rho_uv, sigma_current,
		symmetric, delta_u, delta_v = NULL, prior_shape = 2, prior_scale = 1) {
	n_time <- dim(U_cube)[3]
	n_rank <- dim(U_cube)[2]
	if(is.null(n_time) || n_time < 2 || is.null(delta_u)) return(sigma_current)
	ss <- 0
	n_trans <- 0

	for(t in 2:n_time) {
		keep_u <- delta_u[, t] == 0
		if(any(keep_u)) {
			innov_u <- U_cube[keep_u, , t, drop = FALSE] -
				rho_uv * U_cube[keep_u, , t - 1, drop = FALSE]
			ss <- ss + sum(innov_u * innov_u)
			n_trans <- n_trans + sum(keep_u)
		}
		if(!symmetric && !is.null(delta_v)) {
			keep_v <- delta_v[, t] == 0
			if(any(keep_v)) {
				innov_v <- V_cube[keep_v, , t, drop = FALSE] -
					rho_uv * V_cube[keep_v, , t - 1, drop = FALSE]
				ss <- ss + sum(innov_v * innov_v)
				n_trans <- n_trans + sum(keep_v)
			}
		}
	}
	count <- n_trans * n_rank
	if(!is.finite(ss) || count <= 0) return(sigma_current)
	shape <- count / 2 + prior_shape
	scale <- ss / 2 + prior_scale
	sqrt(1 / rgamma(1, shape = shape, rate = 1 / scale))
}

.lame_sample_rho_uv_t_weighted <- function(U_cube, V_cube, sigma_uv, rho_current,
		symmetric, lambda_u, lambda_v = NULL, prior_mean = 0.9, prior_sd = 0.1) {
	n_time <- dim(U_cube)[3]
	if(is.null(n_time) || n_time < 2 || is.null(lambda_u)) return(rho_current)
	sigma2_inv <- 1 / (sigma_uv * sigma_uv)
	prior_prec <- 1 / (prior_sd * prior_sd)
	sum_yt_yt1 <- 0
	sum_yt1_yt1 <- 0

	for(t in 2:n_time) {
		w_u <- lambda_u[, t]
		keep_u <- is.finite(w_u) & w_u > 0
		if(any(keep_u)) {
			n_keep <- sum(keep_u)
			n_rank <- dim(U_cube)[2]
			U_t <- matrix(U_cube[keep_u, , t], nrow = n_keep, ncol = n_rank)
			U_t1 <- matrix(U_cube[keep_u, , t - 1], nrow = n_keep, ncol = n_rank)
			w <- matrix(w_u[keep_u], nrow = n_keep, ncol = n_rank)
			sum_yt_yt1 <- sum_yt_yt1 + sum(w * U_t * U_t1)
			sum_yt1_yt1 <- sum_yt1_yt1 + sum(w * U_t1 * U_t1)
		}
		if(!symmetric && !is.null(lambda_v)) {
			w_v <- lambda_v[, t]
			keep_v <- is.finite(w_v) & w_v > 0
			if(any(keep_v)) {
				n_keep <- sum(keep_v)
				n_rank <- dim(V_cube)[2]
				V_t <- matrix(V_cube[keep_v, , t], nrow = n_keep, ncol = n_rank)
				V_t1 <- matrix(V_cube[keep_v, , t - 1], nrow = n_keep, ncol = n_rank)
				w <- matrix(w_v[keep_v], nrow = n_keep, ncol = n_rank)
				sum_yt_yt1 <- sum_yt_yt1 + sum(w * V_t * V_t1)
				sum_yt1_yt1 <- sum_yt1_yt1 + sum(w * V_t1 * V_t1)
			}
		}
	}
	if(!is.finite(sum_yt1_yt1) || sum_yt1_yt1 <= 0) return(rho_current)
	var_post <- 1 / (sum_yt1_yt1 * sigma2_inv + prior_prec)
	mean_post <- var_post * (sum_yt_yt1 * sigma2_inv + prior_mean * prior_prec)
	pmax(-0.99, pmin(0.99, rnorm(1, mean_post, sqrt(var_post))))
}

.lame_sample_sigma_uv_t_weighted <- function(U_cube, V_cube, rho_uv, sigma_current,
		symmetric, lambda_u, lambda_v = NULL, prior_shape = 2, prior_scale = 1) {
	n_time <- dim(U_cube)[3]
	n_rank <- dim(U_cube)[2]
	if(is.null(n_time) || n_time < 2 || is.null(lambda_u)) return(sigma_current)
	ss <- 0
	n_trans <- 0

	for(t in 2:n_time) {
		w_u <- lambda_u[, t]
		keep_u <- is.finite(w_u) & w_u > 0
		if(any(keep_u)) {
			n_keep <- sum(keep_u)
			U_t <- matrix(U_cube[keep_u, , t], nrow = n_keep, ncol = n_rank)
			U_t1 <- matrix(U_cube[keep_u, , t - 1], nrow = n_keep, ncol = n_rank)
			innov_u <- U_t - rho_uv * U_t1
			w <- matrix(w_u[keep_u], nrow = n_keep, ncol = n_rank)
			ss <- ss + sum(w * innov_u * innov_u)
			n_trans <- n_trans + sum(keep_u)
		}
		if(!symmetric && !is.null(lambda_v)) {
			w_v <- lambda_v[, t]
			keep_v <- is.finite(w_v) & w_v > 0
			if(any(keep_v)) {
				n_keep <- sum(keep_v)
				n_rank_v <- dim(V_cube)[2]
				V_t <- matrix(V_cube[keep_v, , t], nrow = n_keep, ncol = n_rank_v)
				V_t1 <- matrix(V_cube[keep_v, , t - 1], nrow = n_keep, ncol = n_rank_v)
				innov_v <- V_t - rho_uv * V_t1
				w <- matrix(w_v[keep_v], nrow = n_keep, ncol = n_rank_v)
				ss <- ss + sum(w * innov_v * innov_v)
				n_trans <- n_trans + sum(keep_v)
			}
		}
	}
	count <- n_trans * n_rank
	if(!is.finite(ss) || count <= 0) return(sigma_current)
	shape <- count / 2 + prior_shape
	scale <- ss / 2 + prior_scale
	sqrt(1 / rgamma(1, shape = shape, rate = 1 / scale))
}

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
#' The latent factors U and V evolve according to the transition selected by
#' \code{dynamic_uv_kind}. The default \code{"ar1"} model uses Gaussian AR(1)
#' drift:
#' \deqn{U_{i,k,t} = \rho_{uv} U_{i,k,t-1} + \epsilon_{i,k,t}}
#' where \eqn{\epsilon_{i,k,t} \sim N(0, \sigma_{uv}^2)}, i indexes actors, k indexes
#' latent dimensions, and t indexes time. The parameter \eqn{\rho_{uv}} controls
#' temporal persistence (values near 1 indicate slow evolution). This captures 
#' time-varying homophily, latent community structure, and transitivity dynamics.
#' The \code{"snap"} transition uses a mixture of AR(1) drift and discontinuous
#' reset transitions, while \code{"t"} uses heavy-tailed Student-t innovations
#' represented by local transition scales. On the MCMC path, snap and t
#' transitions are supported for unipartite directed and symmetric models. With
#' \code{method = "als"}, Student-t dynamic UV is also supported for bipartite
#' normal, binary, and Poisson panels, and snap ALS is supported for normal
#' unipartite and bipartite panels, including named panels where actors enter or
#' exit.
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
#' Dynamic additive effects are useful when sender activity and receiver
#' popularity change over time rather than staying fixed across the whole panel.
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
#'   \item Dynamic effects increase computation by roughly 30 to 50 percent per iteration
#'   \item Memory usage scales as O(n*R*T) for dynamic_uv, O(n*T) for dynamic_ab
#'   \item C++ implementation is substantially faster than a pure R loop
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
#' When \code{dynamic_G = TRUE} the bipartite interaction matrix \eqn{G_t} varies
#' by period (a separate \eqn{R_\text{row} \times R_\text{col}} matrix at every
#' time slice) and is returned as \code{fit$G_cube}. The MCMC estimator samples
#' \eqn{G_t} with a Carter-Kohn/FFBS update; \code{method = "als"} estimates a
#' penalized point path for normal, binary, and Poisson bipartite panels,
#' including named changing-composition panels. The marginal
#' \eqn{U_t G_t V_t'} linear predictor is the
#' identified object; individual \eqn{G_t} entries can shift under equivalent
#' rotations and scalings of \eqn{U_t} and \eqn{V_t}. \code{dynamic_G = TRUE}
#' is bipartite only.
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
#' @param dynamic_uv logical: fit dynamic multiplicative effects (latent
#' factors). The transition model is selected by \code{dynamic_uv_kind}. Default
#' \code{FALSE}.
#' @param dynamic_ab logical: fit dynamic additive effects (sender/receiver effects) that 
#' evolve over time using AR(1) processes. When TRUE, the row effects (a) and column 
#' effects (b) become time-varying, following \eqn{a_{i,t} = \rho_{ab} a_{i,t-1} + \epsilon_{i,t}}.
#' This captures temporal heterogeneity in actors' baseline propensities to send and 
#' receive ties, allowing for smooth changes in activity levels and popularity over time.
#' For example, an actor's tendency to form outgoing ties might gradually increase or 
#' decrease across observation periods. The AR(1) specification ensures temporal smoothness
#' while allowing for actor-specific evolution patterns. Implementation uses conjugate 
#' updates where possible and C++ for computational efficiency. Default FALSE.
#' @param dynamic_G logical (bipartite only). When TRUE the bipartite
#' interaction matrix \eqn{G_t} in \eqn{U_t G_t V_t'} varies by period and is
#' returned as \code{fit$G_cube}, an
#' \eqn{R_\text{row} \times R_\text{col} \times T} array. MCMC fits also return
#' posterior summaries for this cube; ALS fits return the penalized point path
#' using the same storage names. The marginal bilinear predictor is identified
#' under the package's canonicalization conventions, while individual entries
#' of \eqn{G_t} can move under equivalent rotations and scalings of
#' \eqn{U_t} and \eqn{V_t}. Default FALSE.
#' @param dynamic_uv_kind Character string selecting the transition model for
#' dynamic multiplicative latent positions when \code{dynamic_uv = TRUE}.
#' \code{"ar1"} uses Gaussian AR(1) drift, \code{"snap"} uses a mixture of
#' AR(1) drift and discontinuous reset transitions, and \code{"t"} uses
#' heavy-tailed Student-t innovations represented by local transition scales.
#' On the MCMC path, \code{"snap"} and \code{"t"} are currently supported for
#' unipartite directed and symmetric models only. With \code{method = "als"},
#' \code{"snap"} is also supported for normal bipartite panels when it is the
#' only dynamic block and \eqn{G} is static. With \code{method = "als"},
#' \code{"t"} is supported for bipartite normal, binary, and Poisson panels.
#' ALS Student-t dynamic UV fits attach the final local transition-weight
#' matrices as \code{fit$lambda_u} and \code{fit$lambda_v}.
#' @param dynamic_beta logical, character, integer, or logical-vector flag
#' selecting which regression coefficients evolve over time via independent
#' AR(1) processes. Default FALSE keeps every coefficient static (the
#' historical behaviour). Accepted forms:
#' \itemize{
#'   \item \code{FALSE} / \code{NULL}: no coefficient is dynamic.
#'   \item \code{TRUE}: every coefficient is dynamic.
#'   \item character vector of block shortcuts (\code{"intercept"},
#'         \code{"dyad"}, \code{"row"}, \code{"col"}) or specific coefficient
#'         names from \code{colnames(BETA)} — those become dynamic.
#'   \item integer vector: 1-based column indices into \code{colnames(BETA)}.
#'   \item logical vector of length \code{p}: per-coefficient mask.
#' }
#' Each dynamic-coefficient block (one per distinct intercept / dyad / row /
#' col label that has at least one dynamic coefficient) gets its own AR(1)
#' parameters \eqn{\rho_\beta} (truncated-Normal prior, default mean 0.8,
#' bounded between 0 and 0.999) and \eqn{\sigma_\beta^2} (inverse-Gamma
#' prior, defaults shape 2 / scale 1). The joint posterior over the time path \eqn{\beta_t}
#' is drawn by Forward-Filter / Backward-Sample (FFBS) inside the MCMC
#' loop, conditional on the current static-block beta and on (a, b, U, V).
#' When the intercept or a nodal coefficient is dynamic, a sum-to-zero
#' contrast basis is applied to the additive-effect sampler to keep
#' \eqn{a_i + intercept_t} identified. Requires at least 2 time periods.
#' With \code{method = "als"}, selected intercept, dyadic, row-node, and
#' column-node coefficients are supported on normal, binary, and Poisson panels,
#' including named panels where actors enter or exit. Dynamic ALS uses
#' period-specific \code{Xrow} or
#' \code{Xcol} values for selected dynamic node coefficients; static node
#' coefficients still use per-actor means. Dynamic intercepts are
#' unavailable when \code{family = "ordinal"} or \code{"rrl"}
#' (and row-nodal coefficients are unavailable with \code{family = "rrl"}).
#' Default FALSE.
#' @param dynamic_beta_kind character: state-space prior on the dynamic
#'   coefficient block(s). \code{"ar1"} (default) gives mean-reverting
#'   AR(1) with a truncated-Normal prior on \eqn{\rho_\beta}. \code{"rw1"}
#'   gives a random walk (\eqn{\rho_\beta} pinned at 1, no truncation),
#'   appropriate when you expect permanent drift with no mean reversion --
#'   e.g. trade-gravity coefficients in a permanently changing world economy.
#'   The alias \code{"random_walk"} is also accepted for \code{"rw1"}.
#'   When the stationarity warning in \code{summary(fit)} fires under the
#'   default AR(1), refit with \code{dynamic_beta_kind = "rw1"} for cleaner
#'   inference. Decision tree: \strong{mean-reverting?} AR(1).
#'   \strong{Permanent / unit-root drift?} RW1.
#'   \strong{Smooth with curvature?} \code{"rw2"} (second-order random walk).
#'   \strong{Smooth with a known length-scale?} \code{"matern32"} (Matern 3/2).
#'   Note: \code{"rw2"} and \code{"matern32"} use an R-level joint Gaussian
#'   sampler and run substantially slower per iteration than the C++ FFBS
#'   used for \code{"ar1"} / \code{"rw1"} (roughly 2-4x on n = 100, T = 10).
#' @param family character: one of "normal","tobit","binary","ordinal","cbin","frn","rrl","poisson" - see
#' the details below
#' @param intercept logical: fit model with an intercept?
#' @param symmetric logical: is the sociomatrix symmetric?
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
#'   \item{rho_beta_mean}{For dynamic_beta: Prior mean for the per-block AR(1)
#'         parameter on time-varying regression coefficients
#'         (default: 0.8). Closer to 1 = smoother evolution.}
#'   \item{rho_beta_sd}{For dynamic_beta: Prior SD for the per-block AR(1)
#'         parameter (default: 0.15).}
#'   \item{rho_beta_lower}{For dynamic_beta: Lower truncation bound on the
#'         AR(1) parameter (default: 0). Pass a negative value to allow
#'         negative autoregression.}
#'   \item{rho_beta_upper}{For dynamic_beta: Upper truncation bound (default:
#'         0.999). Closer to 1 admits near-unit-root behaviour.}
#'   \item{sigma_beta_shape}{For dynamic_beta: Shape parameter for the
#'         inverse-Gamma prior on the per-block innovation variance
#'         (default: 2).}
#'   \item{sigma_beta_scale}{For dynamic_beta: Scale parameter for the
#'         inverse-Gamma prior (default: 1).}
#'   \item{sigma_beta_init}{For dynamic_beta: Initial value of the per-block
#'         innovation standard deviation (default: 0.25). Affects mixing, not
#'         the stationary distribution.}
#'   \item{beta0_mean}{For dynamic_beta: Mean of the Gaussian prior on the
#'         dynamic coefficients at t = 0 (default: 0).}
#'   \item{beta0_var}{For dynamic_beta: Variance of the Gaussian prior on the
#'         dynamic coefficients at t = 0 (default: 10). Weakly informative.}
#' }
#' Common usage: prior = list(Sab0 = diag(c(1, 1)), eta0 = 10) for stronger
#' shrinkage, or prior = list(rho_uv_mean = 0.95) for higher temporal persistence,
#' or prior = list(rho_beta_mean = 0.95, sigma_beta_scale = 0.1) for very smooth
#' time-varying coefficients with tight innovations.
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
#' @param method character: \code{"mcmc"} (default, the Bayesian MCMC fit) or
#'   \code{"als"} (the fast, MCMC-free point-estimation path). With no dynamic
#'   arguments, \code{method = "als"} uses the pooled-static
#'   \code{\link{lame_als}} estimator. For normal, binary,
#'   and Poisson panels, \code{method = "als"} also supports a dynamic
#'   penalized point-estimation path for \code{dynamic_ab}, selected
#'   \code{dynamic_beta}, smooth AR(1) \code{dynamic_uv} in directed,
#'   symmetric, and bipartite panels, bipartite \code{dynamic_G}, and
#'   Student-t \code{dynamic_uv} in directed, symmetric, and bipartite panels.
#'   The snap-only
#'   \code{dynamic_uv = TRUE, dynamic_uv_kind = "snap"} case routes to
#'   \code{\link{lame_snap_als}} for supported normal unipartite and
#'   bipartite panels.
#'   Node-covariate coefficients use the same orthogonal additive-effect
#'   decomposition as \code{\link{lame_als}}; dynamic node coefficients use
#'   period-specific node values when selected by \code{dynamic_beta}, while
#'   static node coefficients use per-actor means. Named changing-composition
#'   panels are aligned to the union actor set and actor-entry gaps break the
#'   dynamic smoothing penalties. Posterior draws, rank/censored families, and
#'   bipartite snap-shift fits with node covariates, dynamic coefficients, or
#'   dynamic \code{G_t} are not handled by the snap ALS shortcut.
#' @param als_stability character (only used when \code{method = "als"} and a
#'   dynamic ALS path is selected): \code{"none"} (default), \code{"quick"}, or
#'   \code{"validation"}. The non-\code{"none"} presets rerun the dynamic ALS
#'   fit from jittered starts and attach \code{fit$stability}. Smooth dynamic
#'   ALS fits report fitted-surface and coefficient-path differences across
#'   starts; snap ALS fits report snap-score, snap-class, top-period, and
#'   fitted-surface differences.
#' @param als_max_iter positive integer (only used when \code{method = "als"}):
#'   maximum block-coordinate iterations for the ALS fit. Default \code{200}.
#' @param als_tol optional positive scalar (only used when
#'   \code{method = "als"}): convergence tolerance for the ALS objective and
#'   fitted values. \code{NULL} keeps each ALS estimator's built-in default.
#' @param bootstrap integer (only used when \code{method = "als"}): number of
#'   bootstrap replicates. \code{0} (default) skips the bootstrap; \code{N > 0}
#'   runs \code{N} replicates and attaches the result so that \code{\link{confint}}
#'   returns bootstrap intervals.
#' @param bootstrap_type character (only used when \code{method = "als"}):
#'   \code{"parametric"} (default) or \code{"block"}.
#' @param bootstrap_block_length integer: block length for the block bootstrap.
#' @param bootstrap_seed optional integer seed for the bootstrap.
#' @param print Deprecated. Use \code{verbose} instead.
#' @param gof logical: calculate goodness of fit statistics?
#' @param start_vals list of parameter starting values for the MCMC chain
#' @param periodic_save logical: indicating whether to periodically save MCMC results
#' @param out_file character vector indicating name and path in which file should be stored if periodic_save is selected. For example, on an Apple OS out_file="~/Desktop/ameFit.rda".
#' @param save_interval quantile interval indicating when to save during the post-burn-in period.
#' @param log_lik_method one of \code{"observed_exact"} (default),
#'   \code{"observed_ghk"}, or \code{"augmented"}. Selects the
#'   argument for selecting which pointwise log-lik is stored on
#'   \code{fit$log_lik}. \code{"observed_exact"} uses the closed-form
#'   marginal log-likelihood for normal, binary, cbin, tobit, poisson,
#'   and ordinal fits. \code{"observed_ghk"} uses a randomized-Halton GHK
#'   Monte Carlo marginal for rank-family fits (\code{"frn"} and
#'   \code{"rrl"}) and the exact observed likelihood when no GHK step is
#'   needed. \code{"augmented"} uses the augmented-data Gaussian-on-\eqn{Z}
#'   contribution.
#' @param keep_snap_draws controls storage of post-burn-in MCMC snap-shift
#'   indicators when \code{dynamic_uv = TRUE} and
#'   \code{dynamic_uv_kind = "snap"}. \code{"none"} keeps only the posterior
#'   mean snap probabilities. \code{"summary"} stores cell-level count, mean,
#'   and variance summaries. \code{"draws"} stores retained binary indicator
#'   arrays as \code{fit$snap_draws} and, for directed fits,
#'   \code{fit$snap_draws_v}, with dimensions draw x actor x time and first
#'   period set to \code{NA}. \code{"chunked"} currently stores the same arrays
#'   in memory and labels the storage mode for downstream tooling.
#' @param dynamic_beta_per_actor optional, one of \code{NULL} (default),
#'   \code{"row"}, or \code{"col"}. When non-\code{NULL}, fits
#'   actor-specific time-varying deviations for one dyadic covariate on the
#'   selected margin. The default \code{per_actor_identifiability = "center"}
#'   keeps the population coefficient in \code{coef(fit)} and enforces a
#'   per-period sum-to-zero constraint on the actor deviations. Draws are
#'   stored in \code{fit$THETA_ACTOR} when retained; summary mode stores
#'   \code{fit$theta_actor_mean} and \code{fit$theta_actor_sd}.
#' @param per_actor_covariate_idx positive integer; index into the
#'   dyadic covariate cube to slope on for the per-actor extension.
#'   Default \code{1L}.
#' @param per_actor_identifiability one of \code{"center"} (default)
#'   or \code{"drop_population"}. \code{"center"} preserves the
#'   population coefficient and constrains per-actor deviations to
#'   sum to zero per period via pairwise contrast FFBS. \code{"drop_population"}
#'   is reserved and currently falls back to \code{"center"}.
#' @param keep_per_actor one of \code{"auto"} (default; full draws
#'   when memory cost < 250 MB else streaming summary), \code{"draws"}
#'   (full \code{[n_iter, n_actors, T]} cube), \code{"summary"}
#'   (streaming posterior mean + sd only), or \code{"none"} (only
#'   hyperparameter chains).
#' @param dynamic_beta_pool one of \code{"none"} (default), \code{"rho"},
#'   \code{"sigma"}, or \code{"both"}. When two or more coefficient blocks are
#'   dynamic, non-\code{"none"} values update shared hyperparameters across
#'   blocks for the selected persistence and / or innovation-scale component.
#'   With a single dynamic block the value is recorded but there is no other
#'   block to pool against.
#' @param time_index optional numeric vector of length \eqn{T} giving
#'   observation times. Default \code{NULL} treats periods as equally
#'   spaced. Strictly-increasing values are required when supplied. Unequal
#'   gaps scale the AR(1) / RW1 transition variance and are also used by the
#'   RW2 and Matérn 3/2 dynamic-beta precision builders.
#' @param period_exposure optional non-negative numeric vector of length
#'   \eqn{T} giving period-level exposure offsets. **Wired for Poisson:**
#'   when supplied with any value != 1, the Poisson observation
#'   likelihood becomes
#'   \eqn{Y_{ij,t} \sim \text{Poisson}(e_t \cdot \exp(Z_{ij,t}))},
#'   while the latent \eqn{Z} retains its existing semantic of
#'   "unexposed log-rate". For other families, non-trivial
#'   \code{period_exposure} is an error (rescale \code{Y} or use a
#'   covariate offset). \code{NULL} (default) or all-ones uses the
#'   unscaled likelihood path.
#' @param max_seconds optional positive scalar; if the MCMC wall-clock
#'   time exceeds this many seconds, the chain terminates cleanly and
#'   \code{fit$terminated_early} is set to \code{TRUE}.
#' @param checkpoint_path optional file path. When set, the chain
#'   periodically writes a snapshot of \code{BETA}, \code{VC}, the RNG
#'   state, and the original call to this file. Use \code{\link{lame_resume}}
#'   or pass the same path back as \code{resume_from = path} to
#'   continue.
#' @param checkpoint_every positive integer; iterations between
#'   checkpoint writes. Default \code{100L}.
#' @param resume_from optional path to a checkpoint file produced by a
#'   previous \code{lame(..., checkpoint_path = path)} call. When
#'   non-\code{NULL}, \code{lame()} short-circuits all other input
#'   parsing and delegates to \code{\link{lame_resume}} with the
#'   user-supplied overrides forwarded. Pass \code{nscan = K} to
#'   request \code{K} additional stored draws on the continuation.
#'   Note: \code{checkpoint_path} and \code{max_seconds} cannot
#'   currently be overridden on a resume (they are stripped before
#'   re-evaluating the saved call); open both at the original
#'   \code{lame()} call if you need iterative checkpointing or a
#'   time budget on a chain of resumes. Pass \code{verbose = FALSE}
#'   on resume (the sampling progress bar is gated on \code{burn != 0}
#'   and the continuation forces \code{burn = 0}).
#'   Default \code{NULL}.
#' @param freeze_call logical: if \code{TRUE}, store a snapshot of the
#'   evaluated \code{Y}, \code{Xdyad}, \code{Xrow}, \code{Xcol} on
#'   \code{fit$data_snapshot} so that a later \code{update(fit, ...)} refits
#'   against the same data even if the caller has mutated those objects in
#'   their workspace. Memory cost equals the size of the data; default
#'   \code{FALSE}.
#' @param save_log_lik one of \code{FALSE} (default), \code{TRUE}, or
#'   \code{"chunked"}. When \code{TRUE}, stores the per-iteration pointwise
#'   log-likelihood matrix on \code{fit$log_lik} (an
#'   \code{[n_stored, n_obs]} double matrix). When \code{"chunked"},
#'   streams the log-lik values to per-column-chunk binary files under
#'   \code{log_lik_path} so the in-memory cost during MCMC is just one
#'   chunk's width; the chunks are recovered later via
#'   \code{read_log_lik(fit)}. Required for \code{loo::loo(fit)} and
#'   \code{loo::waic(fit)}.
#' @param ordinal_cutpoints character: cutpoint convention for
#'   \code{family = "ordinal"}. \code{"data_induced"} (default) uses the
#'   data-induced cutpoints; \code{"explicit"} samples explicit cutpoints via
#'   a Cowles (1996) Metropolis-Hastings update. Ignored for other families.
#' @param log_lik_path directory to write log-lik chunks to when
#'   \code{save_log_lik = "chunked"}. Default \code{NULL} creates a
#'   session-scoped subdirectory under \code{tempdir()}.
#' @param log_lik_chunk_size column-width of each on-disk chunk when
#'   \code{save_log_lik = "chunked"}. Larger chunks mean fewer files (and
#'   one \code{readBin} per chunk on read) but more memory during MCMC.
#'   Default \code{10000L}.
#' @param model.name optional string for model selection output
#' @return \item{BETA}{posterior samples of regression coefficients. A
#'   2-dimensional matrix \code{[n_stored, p]} when all coefficients are
#'   static. When \code{dynamic_beta} flags any coefficient as time-varying,
#'   \code{BETA} is a 3-dimensional array \code{[n_stored, p, T]} whose third
#'   dimension is the time period. Static coefficients are still present with
#'   their values replicated across the third dimension. \code{coef(fit)}
#'   collapses this to a \code{[p, T]} posterior-mean matrix.
#'
#'   \strong{Migration from amen.} Under \code{amen::ame()} the
#'   \code{BETA} slot was always 2-D. Scripts that compute
#'   \code{apply(fit$BETA, 2, mean)} on an amen fit will silently
#'   aggregate across periods when run against a \code{lame()} fit
#'   with \code{dynamic_beta} active (the second margin is the
#'   coefficient index in both shapes; the new third margin is time).
#'   Use \code{length(dim(fit$BETA))} to detect the shape, or call
#'   \code{coef(fit)} which returns a \code{[p, T]} matrix in either
#'   case.}
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
#' matrix. For \code{mode = "bipartite"}, \code{EZ} is a list of
#' per-period matrices whose row/column dimnames are currently
#' \code{NULL}. Use \code{names(fit$APM)} / \code{names(fit$BPM)} (or
#' \code{dimnames(fit$YPM[[t]])}) to recover the row/column actor
#' ordering, which is the lexicographic sort of the input actor names
#' (so \code{"r10"} sorts before \code{"r2"} unless you zero-pad).}
#' \item{YPM}{posterior mean of Y (for imputing missing values)}
#' \item{GOF}{observed (first row) and posterior predictive (remaining rows)
#' values of four goodness-of-fit statistics.
#' See \code{\link{gof}} for post-hoc computation and \code{\link{gof_plot}} for visualization.}
#' \item{start_vals}{Final parameter values from MCMC, can be used as the input
#' for a future model run.}
#' \item{model.name}{Name of the model (if provided)}
#' @seealso \code{\link{ame}} for cross-sectional models,
#'   \code{\link{lame_als}} for the fast MCMC-free point estimator,
#'   \code{\link{lame_snap_als}} for the approximate dynamic snap-shift
#'   point estimator,
#'   \code{\link{als_dynamic_beta}} for a regression-only penalised
#'   smoother on the time-varying coefficient path,
#'   \code{\link{lame_resume}} for resuming saved checkpoints,
#'   \code{\link{gof}} for post-hoc goodness-of-fit computation,
#'   \code{\link{gof_plot}} for visualizing GOF results,
#'   \code{\link{latent_positions}} for extracting latent positions as a tidy data frame,
#'   \code{\link{procrustes_align}} for Procrustes alignment of latent positions,
#'   \code{\link{summary.lame}} for model summaries,
#'   \code{\link{coef.lame}} for coefficient extraction,
#'   \code{\link{netify_to_lame}} for the recommended netify -> lame
#'   bridge (set \code{lame = TRUE} on \code{netify::to_amen()}).
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#'
#' data(YX_bin_list)
#' fit<-lame(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,family="binary")
#' # you should run the Markov chain much longer than this
#'
#' \donttest{
#' ## Time-varying regression coefficients (dynamic_beta).
#' ## Make every dyadic coefficient evolve as an AR(1):
#' fit_dyn <- lame(YX_bin_list$Y, YX_bin_list$X,
#'                 family = "binary", R = 0,
#'                 nscan = 200, burn = 50, odens = 5,
#'                 dynamic_beta = "dyad")
#' dim(fit_dyn$BETA)        # [n_stored, p, T] -- 3-D when dynamic
#' coef(fit_dyn)            # [p, T] posterior-mean coefficient paths
#' confint(fit_dyn)         # per-period 95% credible intervals
#' summary(fit_dyn)         # prints a "Dynamic coefficients per period" block
#' }
#'
#' @export lame
lame <- function(
		Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL,
		rvar = !(family=="rrl") , cvar = TRUE, dcor = !symmetric,
		nvar = TRUE,
		R = 0, R_row = NULL, R_col = NULL,  # separate ranks for bipartite
		mode = c("unipartite", "bipartite"),  # network mode
		dynamic_uv = FALSE, dynamic_ab = FALSE, dynamic_G = FALSE,
		dynamic_beta = FALSE,
		dynamic_beta_kind = c("ar1", "rw1", "rw2", "matern32"),
		dynamic_uv_kind = c("ar1", "snap", "t"),
		family = "normal",
		intercept = !is.element(family,c("rrl","ordinal")),
		symmetric = FALSE,
		odmax = NULL,
		prior = list(), g = NA,
		seed = 6886, nscan = 10000, burn = 500, odens = 25,
		plot = FALSE, verbose = FALSE, gof = TRUE,
		start_vals = NULL, periodic_save=FALSE, out_file=NULL,
		save_interval=0.25, model.name = NULL,
		save_log_lik = FALSE,
		log_lik_path = NULL,
		log_lik_chunk_size = 10000L,
		keep_snap_draws = c("none", "summary", "draws", "chunked"),
		freeze_call = FALSE,
		dynamic_beta_pool = c("none", "rho", "sigma", "both"),
		dynamic_beta_per_actor = NULL,
		per_actor_covariate_idx = 1L,
		per_actor_identifiability = c("center", "exact_center", "drop_population"),
		keep_per_actor = c("auto", "draws", "summary", "none"),
		time_index = NULL,
		period_exposure = NULL,
		max_seconds = Inf,
		checkpoint_path = NULL,
		checkpoint_every = 100L,
			log_lik_method = c("observed_exact", "observed_ghk", "augmented"),
			ordinal_cutpoints = c("data_induced", "explicit"),
			method = c("mcmc", "als"),
			als_stability = c("none", "quick", "validation"),
			als_max_iter = 200L,
			als_tol = NULL,
			bootstrap = 0L,
		bootstrap_type = c("parametric", "block"),
		bootstrap_block_length = 1L,
		bootstrap_seed = NULL,
		resume_from = NULL,
		print
) {
	# helper function
	`%||%` <- function(x, y) if (is.null(x)) y else x

	# short-circuit to checkpoint-resume when resume_from is supplied. this
	# is the consolidated entry point for lame_resume(); user overrides are
	# forwarded through. the underlying implementation lives in lame_resume()
	# (which reads the saved call from disk and re-evaluates this function).
	# `nscan_more` is sugar for `nscan` on the continuation side.
	#
	# when the user calls `lame(resume_from = path, ...)` without supplying
	# y, the `lame()` frame has y as a missing promise. capture the caller's
	# frame and forward it to lame_resume() via `.envir`, so the saved call's
	# y / xdyad references resolve against the user's bindings, not this
	# slot.
	if (!is.null(resume_from)) {
		caller_env <- parent.frame()
		mc0 <- match.call()
		# user-supplied named args, minus resume_from (resolved here) and
		# any args whose value is still a missing promise (the user did
		# not actually pass them)
		user_named <- setdiff(names(as.list(mc0))[-1L], "resume_from")
		usr <- list()
		for (nm in user_named) {
			# mget with ifnotfound = na prevents a force on a missing-promise
			val <- tryCatch(get(nm, envir = environment(), inherits = FALSE),
			                error = function(e) NULL)
			if (!is.null(val)) usr[[nm]] <- val
		}
		nscan_more <- usr$nscan
		usr$nscan <- NULL
		return(do.call(lame_resume,
			c(list(path = resume_from, nscan_more = nscan_more, .envir = caller_env),
			  usr)))
	}

	# handle deprecated print argument
	if (!missing(print)) {
		# emit the warning only when the user did not also pass verbose;
		# multi-chain wrappers like ame_parallel pass both for back-compat
		# and the duplicate warning per chain is noise.
		if (missing(verbose)) {
			cli::cli_warn(c(
				"The {.arg print} argument is deprecated.",
				"i" = "Use {.arg verbose} instead."
			))
			verbose <- print
		}
	}

	# capture the call so summary.lame can echo it
	mc <- match.call()
	# snapshot the user's data before any internal conversion
	# (list -> 3d array via list_to_array, etc.) so update.lame can re-fit
	# against the original input shape. populated only when
	# freeze_call = true.
	if (isTRUE(freeze_call)) {
		.lame_freeze_snapshot <- list(
			Y     = Y,
			Xdyad = Xdyad,
			Xrow  = Xrow,
			Xcol  = Xcol
		)
	} else {
		.lame_freeze_snapshot <- NULL
	}
	# process mode argument
	mode <- match.arg(mode)
	bip <- identical(mode, "bipartite")
	# validate family early so every front-door route, including als, gets the
	# same amen-era aliases.
	valid_families <- c("normal", "tobit", "binary", "ordinal",
	                    "cbin", "frn", "rrl", "poisson")
	amen_aliases <- c(nrm = "normal", bin = "binary", ord = "ordinal",
	                  pois = "poisson", tob = "tobit")
	if (length(family) == 1L && is.character(family) &&
	    family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		cli::cli_inform(c(
			"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
			"i" = "Use {.code family = \"{new_family}\"} to avoid the alias message."))
		family <- new_family
	}
	if (length(family) != 1L || !is.character(family) ||
	    !family %in% valid_families) {
		cli::cli_abort(c(
			"{.arg family} = {.val {family}} is not a recognised family.",
			"i" = "Choose one of: {.val {valid_families}}."))
	}
	# resolve dynamic_beta_kind. accept "random_walk" as an alias for
	# "rw1". "rw2" and "matern32" use the r-level joint gaussian sampler.
	if (length(dynamic_beta_kind) > 1L) dynamic_beta_kind <- dynamic_beta_kind[[1L]]
	if (identical(dynamic_beta_kind, "random_walk")) dynamic_beta_kind <- "rw1"
	if (!dynamic_beta_kind %in% c("ar1", "rw1", "rw2", "matern32")) {
		cli::cli_abort(c(
			"{.arg dynamic_beta_kind} must be one of {.val ar1}, {.val rw1}, {.val rw2}, {.val matern32}.",
			"i" = "Got {.val {dynamic_beta_kind}}.",
			"i" = "{.val random_walk} is accepted as an alias for {.val rw1}."))
	}
	dynamic_beta_kind_requested <- dynamic_beta_kind
	# resolve dynamic_uv_kind. "ar1" is the default smooth drift; "snap"
	# and "t" are explicit dynamic-uv transition models.
	dynamic_uv_kind <- match.arg(dynamic_uv_kind)
	if (!isTRUE(dynamic_uv) && !identical(dynamic_uv_kind, "ar1")) {
		cli::cli_abort(c(
			"{.arg dynamic_uv_kind = \"{dynamic_uv_kind}\"} requires {.arg dynamic_uv = TRUE}.",
			"i" = "Set {.arg dynamic_uv = TRUE} or leave {.arg dynamic_uv_kind} at {.val ar1}."))
	}
	dynamic_uv_kind_requested <- dynamic_uv_kind
	# rw2 and matern32 use the joint gaussian sampler

	# choose the static, dynamic, or snap point estimator
		method <- match.arg(method)
		als_stability <- match.arg(als_stability)
		if (!is.numeric(als_max_iter) || length(als_max_iter) != 1L ||
		    anyNA(als_max_iter) || als_max_iter <= 0 ||
		    !isTRUE(all.equal(als_max_iter, round(als_max_iter)))) {
			cli::cli_abort("{.arg als_max_iter} must be a positive integer.")
		}
		als_max_iter <- as.integer(als_max_iter)
		if (!is.null(als_tol) &&
		    (!is.numeric(als_tol) || length(als_tol) != 1L ||
		     !is.finite(als_tol) || als_tol <= 0)) {
			cli::cli_abort("{.arg als_tol} must be NULL or a positive finite number.")
		}
	if (identical(method, "als")) {
		if (!isTRUE(intercept)) {
			cli::cli_abort(c(
				"{.fn lame} {.arg method} = {.val als} does not currently support {.code intercept = FALSE}.",
				"i" = "Use {.code method = \"mcmc\"} for a no-intercept longitudinal fit, or keep the ALS intercept and center covariates as needed."))
		}
		bootstrap_type <- match.arg(bootstrap_type)
		if (!is.numeric(bootstrap) || length(bootstrap) != 1L || bootstrap < 0L ||
		    !isTRUE(all.equal(bootstrap, round(bootstrap)))) {
			cli::cli_abort("{.arg bootstrap} must be a non-negative integer (0 = no bootstrap).")
		}
		bootstrap <- as.integer(bootstrap)

		dynamic_beta_requested <- !is.null(dynamic_beta) &&
		    !(is.logical(dynamic_beta) && length(dynamic_beta) == 1L && isFALSE(dynamic_beta))
		als_user_args <- names(as.list(mc))[-1L]
		warn_als_mcmc_controls <- function(route) {
			ignored <- intersect(als_user_args, c(
				"nscan", "burn", "odens", "g", "keep_snap_draws",
				"save_log_lik", "log_lik_path", "log_lik_chunk_size",
				"log_lik_method", "start_vals", "periodic_save", "out_file",
				"save_interval", "freeze_call", "dynamic_beta_pool",
				"dynamic_beta_per_actor", "per_actor_covariate_idx",
				"per_actor_identifiability", "keep_per_actor", "time_index",
				"period_exposure", "max_seconds", "checkpoint_path",
				"checkpoint_every", "ordinal_cutpoints"))
			if (length(ignored) > 0L) {
				ignored_text <- paste(ignored, collapse = ", ")
				cli::cli_warn(c(
					"{.fn lame} {.arg method} = {.val als}: ignored MCMC-only arguments on the {route} route: {ignored_text}.",
					"i" = "Use {.code method = \"mcmc\"} if you need those sampler controls."))
			}
		}

		if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "snap") &&
		    !isTRUE(dynamic_ab) && !isTRUE(dynamic_G) && !dynamic_beta_requested) {
			warn_als_mcmc_controls("snap ALS")
			if (bootstrap > 0L) {
				cli::cli_warn(c(
					"{.arg bootstrap} is not available for snap ALS and is ignored.",
					"i" = "{.arg als_stability} is still available for start-sensitivity checks."))
			}
			if (verbose) {
				cli::cli_inform(c(
					"i" = "Using {.fn lame_snap_als} for snap dynamic UV with {.arg method} = {.val als}.",
					"i" = "Returned snap scores are ALS scores, not MCMC posterior draws."))
			}
			snap_pi_prior_als <- c(prior$snap_pi_a %||% 1,
			                       prior$snap_pi_b %||% 9)
			fit <- lame_snap_als(
				Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
				R = R, R_row = R_row, R_col = R_col,
				family = family, mode = mode, symmetric = symmetric,
				snap_kappa = prior$snap_kappa %||% 2,
				snap_pi_prior = snap_pi_prior_als,
				rho_prior_mean = prior$rho_uv_mean %||% 0.98,
				max_iter = als_max_iter,
				tol = if (is.null(als_tol)) 1e-6 else als_tol,
				stability = als_stability,
				verbose = verbose, seed = seed)
			fit$call <- mc
			return(fit)
		}

		dyn_requested <- character()
		if (isTRUE(dynamic_uv)) dyn_requested <- c(dyn_requested, "dynamic_uv")
		if (isTRUE(dynamic_ab)) dyn_requested <- c(dyn_requested, "dynamic_ab")
		if (isTRUE(dynamic_G)) dyn_requested <- c(dyn_requested, "dynamic_G")
		if (dynamic_beta_requested) {
			dyn_requested <- c(dyn_requested, "dynamic_beta")
		}
		# static lame_als() uses one latent rank
		R_als <- R
		R_row_als <- R_row
		R_col_als <- R_col
		if (!is.null(R_row) || !is.null(R_col)) {
			rr <- R_row %||% R; rc <- R_col %||% R
			if (rr != rc && !(length(dyn_requested) > 0L &&
			    identical(mode, "bipartite"))) {
				cli::cli_warn(c(
					"{.arg R_row} ({rr}) != {.arg R_col} ({rc}), but the ALS estimator uses a single latent dimension.",
					"i" = "Using {.code R = {max(rr, rc)}} for the ALS fit (the larger of the two)."))
			}
			R_als <- max(rr, rc, R)
			R_row_als <- rr
			R_col_als <- rc
		}
		if (length(dyn_requested) > 0L) {
			warn_als_mcmc_controls("dynamic ALS")
			fit <- lame_dynamic_als(
				Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
				R = R_als, R_row = R_row_als, R_col = R_col_als,
				family = family, mode = mode, symmetric = symmetric,
				dynamic_ab = dynamic_ab,
				dynamic_uv = dynamic_uv,
				dynamic_uv_kind = dynamic_uv_kind,
				dynamic_beta = dynamic_beta,
				dynamic_beta_kind = dynamic_beta_kind,
				dynamic_G = dynamic_G,
				prior = prior,
				bootstrap = bootstrap,
				bootstrap_type = bootstrap_type,
				bootstrap_block_length = bootstrap_block_length,
				bootstrap_seed = bootstrap_seed,
				stability = als_stability,
				max_iter = als_max_iter,
				tol = if (is.null(als_tol)) 5e-5 else als_tol,
				verbose = verbose, seed = seed, call = mc)
			fit$call <- mc
			return(fit)
		}
		user_args <- names(as.list(mc))[-1L]
		forwarded <- c("Y", "Xdyad", "Xrow", "Xcol", "R", "R_row", "R_col",
		               "family", "mode", "symmetric",
		               "bootstrap", "bootstrap_type",
		               "bootstrap_block_length", "bootstrap_seed",
		               "als_stability", "als_max_iter", "als_tol",
		               "verbose", "seed", "method", "intercept", "odmax",
		               # absorb cosmetic-only mcmc args that have no als
		               # equivalent and don't affect the point fit
		               "plot", "model.name")
		dropped <- setdiff(user_args, forwarded)
		if (length(dropped) > 0L) {
			cli::cli_warn(c(
				"{.fn lame} {.arg method} = {.val als}: ignored {length(dropped)} arguments that apply only to the MCMC path: {.arg {dropped}}.",
				"i" = "Use {.code method = \"mcmc\"} (the default) if you need custom MCMC-only features."))
		}
		fit <- lame_als(
			Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
			R = R_als, family = family, mode = mode, symmetric = symmetric,
			bootstrap = bootstrap, bootstrap_type = bootstrap_type,
			bootstrap_block_length = bootstrap_block_length,
			bootstrap_seed = bootstrap_seed,
			max_iter = als_max_iter,
			tol = if (is.null(als_tol)) 1e-6 else als_tol,
			verbose = verbose, seed = seed)
		fit$call <- mc
		return(fit)
	}

	# log-lik method labeling. observed_exact uses the closed-form
	# marginal; observed_ghk wires the ghk monte carlo with a downward-
	# bias note; augmented falls back to the gaussian-on-z fallback.
	# observed_exact  - closed-form pointwise lik for {normal, binary, cbin,
	#                   tobit, poisson, ordinal}
	# observed_ghk    - ghk monte carlo for rank-family marginals; log of an
	#                   unbiased probability estimate is biased downward.
	# augmented       - the existing gaussian-on-z augmented-data fallback
	#                   for cbin/frn/rrl/ordinal
	log_lik_method <- match.arg(log_lik_method)
	keep_snap_draws <- match.arg(keep_snap_draws)
	if (identical(log_lik_method, "observed_ghk")) {
		# the ghk marginal log-likelihood is a monte carlo estimate and
		# log of an unbiased probability estimate is biased downward; surface
		# the caveat once per fit.
		cli::cli_inform(c(
			"i" = "{.arg log_lik_method} = {.val observed_ghk}: GHK Monte Carlo marginal for rank/ordinal families.",
			"!" = "The log of an unbiased probability estimate is biased downward; treat this as a Monte Carlo estimate."))
	}

	# ordinal cutpoints use either data-induced thresholds or explicit mh updates
	ordinal_cutpoints <- match.arg(ordinal_cutpoints)
	if (identical(ordinal_cutpoints, "explicit") && family != "ordinal") {
		cli::cli_warn(c(
			"{.arg ordinal_cutpoints} = {.val explicit} only takes effect when {.code family = \"ordinal\"}.",
			"i" = "Argument ignored for {.code family = \"{family}\"}."))
		ordinal_cutpoints <- "data_induced"
	}

	# dynamic_beta_pool, time_index, and period_exposure handle dynamic
	# coefficient pooling and exposure-aware paths.
	dynamic_beta_pool <- match.arg(dynamic_beta_pool)
	# default ("none") runs per-block ig / truncated-normal updates;
	# "rho" / "sigma" / "both" add the hierarchical pool across blocks.
	if (!is.null(time_index)) {
		if (!is.numeric(time_index)) {
			cli::cli_abort("{.arg time_index} must be numeric.")
		}
		if (any(diff(time_index) <= 0)) {
			cli::cli_abort("{.arg time_index} must be strictly increasing.")
		}
		# unequal gaps use the joint gaussian sampler
	}
	if (!is.null(period_exposure)) {
		if (!is.numeric(period_exposure) || any(!is.finite(period_exposure)) ||
		    any(period_exposure < 0)) {
			cli::cli_abort("{.arg period_exposure} must be a non-negative finite numeric vector.")
		}
		if (family != "poisson" && any(period_exposure != 1)) {
			cli::cli_abort(c(
				"{.arg period_exposure} != 1 is currently meaningful only for {.code family = \"poisson\"}.",
				"i" = "Got {.val {family}} family with non-trivial exposure.",
				"i" = "For other families, scale {.arg Y} manually or use a covariate offset."))
		}
		# poisson exposure path: no note when exposures are all-ones
		# (the unscaled path runs), only when any value differs.
	}
	# precompute the activation flag so the inner loop can dispatch
	# without re-checking on every iteration. null or all-ones routes
	# through the unscaled likelihood path.
	poisson_exposure_active <- !is.null(period_exposure) &&
		identical(family, "poisson") &&
		any(period_exposure != 1)
	# dynamic_beta_per_actor: pairwise contrast ffbs with sum-to-zero
	# per period. "center" keeps the population coefficient (default);
	# "drop_population" would remove the selected covariate from the
	# static beta block.
	per_actor_identifiability <- match.arg(per_actor_identifiability)
	keep_per_actor <- match.arg(keep_per_actor)
	per_actor_active <- !is.null(dynamic_beta_per_actor)
	if (per_actor_active) {
		if (!is.character(dynamic_beta_per_actor) ||
		    !dynamic_beta_per_actor %in% c("row", "col")) {
			cli::cli_abort("{.arg dynamic_beta_per_actor} must be one of {.val NULL}, {.val row}, or {.val col}.")
		}
		# drop_population is reserved; the fitted path stays centered
		if (identical(per_actor_identifiability, "drop_population")) {
			cli::cli_warn(c(
				"!" = "{.arg per_actor_identifiability} = {.val drop_population} is reserved but not available yet.",
				"i" = "Using {.val center}: zero-sum per-actor deviations on top of the population coefficient.",
				"i" = "Pass {.val center} to silence this warning."))
		}
		if (!is.numeric(per_actor_covariate_idx) ||
		    per_actor_covariate_idx < 1L) {
			cli::cli_abort("{.arg per_actor_covariate_idx} must be a positive integer.")
		}
	}
	# max_seconds and checkpointing
	if (!is.null(max_seconds) && (!is.numeric(max_seconds) || max_seconds <= 0)) {
		cli::cli_abort("{.arg max_seconds} must be a positive number or {.val Inf}.")
	}
	if (!is.null(checkpoint_path)) {
		ck_dir <- dirname(checkpoint_path)
		if (!dir.exists(ck_dir)) dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
	}
	if (!is.numeric(checkpoint_every) || checkpoint_every <= 0) {
		cli::cli_abort("{.arg checkpoint_every} must be a positive integer.")
	}
	.lame_start_time <- Sys.time()
	.lame_terminated_early <- FALSE
	# bipartite ordinal / cbin / frn / poisson / rrl use rectangular
	# samplers; see r/rz_bipartite.r. each family is dispatched in the
	# rectangular z-sample branch later in this function.

	# convert a 3d y array to list format if needed
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

	# validate y and the key mcmc arguments up front
	if (!is.list(Y) || length(Y) == 0L) {
		cli::cli_abort("{.arg Y} must be a non-empty list of relational matrices (or a 3D array).")
	}
	obs_y <- unlist(lapply(Y, function(y) y[is.finite(y)]), use.names = FALSE)
	if (length(obs_y) == 0L) {
		cli::cli_abort("{.arg Y} has no observed (non-missing, finite) values.")
	}
	if (length(nscan) != 1L || anyNA(nscan) || nscan < 1L) {
		cli::cli_abort("{.arg nscan} must be a positive integer (post-burn-in iterations).")
	}
	if (length(burn) != 1L || anyNA(burn) || burn < 0L) {
		cli::cli_abort("{.arg burn} must be a non-negative integer.")
	}
	if (length(odens) != 1L || anyNA(odens) || odens < 1L || odens > nscan) {
		cli::cli_abort(c(
			"{.arg odens} must be a positive integer no greater than {.arg nscan} = {nscan}.",
			"i" = "Otherwise no posterior draws would be stored."))
	}
	if (family == "poisson" && any(obs_y < 0)) {
		cli::cli_abort(c(
			"{.arg family} = {.val poisson} requires non-negative counts.",
			"i" = "{.arg Y} has negative observed values."))
	}
	# all-na row / col guard: a slice with an entirely na row or column
	# would otherwise abort the bipartite-binary z sampler deep in the
	# stack. surface as a clean error.
	.lame_check_all_na <- function(Y_per, t_label) {
		if (is.null(Y_per)) return(invisible())
		drow <- apply(Y_per, 1, function(r) all(is.na(r) | !is.finite(r)))
		dcol <- apply(Y_per, 2, function(c) all(is.na(c) | !is.finite(c)))
		if (bip) {
			if (any(drow) || any(dcol)) {
				if (family %in% c("binary", "tobit")) {
					cli::cli_abort(c(
						"{.arg Y[[{t_label}]]} (bipartite, {.val {family}}) has {sum(drow)} fully-unobserved row{?s} and {sum(dcol)} fully-unobserved column{?s}.",
						"x" = "The bipartite {.val {family}} sampler cannot proceed with rows/columns that have no observed data.",
						"i" = "Drop the empty actor(s) from {.arg Y} (and matching covariates) before fitting, or use {.val normal} family."))
				}
				cli::cli_inform(c(
					"i" = "{.arg Y[[{t_label}]]} has {sum(drow)} fully-unobserved row{?s} and {sum(dcol)} fully-unobserved column{?s}.",
					"i" = "Their additive effects will be pooled with the prior for that period."))
			}
		} else {
			# unipartite: drop the diagonal na from the all-na judgement
			drow_ex <- vapply(seq_len(nrow(Y_per)), function(i)
				all(is.na(Y_per[i, -i]) | !is.finite(Y_per[i, -i])), logical(1))
			dcol_ex <- vapply(seq_len(ncol(Y_per)), function(j)
				all(is.na(Y_per[-j, j]) | !is.finite(Y_per[-j, j])), logical(1))
			if (any(drow_ex) || any(dcol_ex)) {
				cli::cli_inform(c(
					"i" = "{.arg Y[[{t_label}]]} has {sum(drow_ex)} fully-unobserved row{?s} and {sum(dcol_ex)} fully-unobserved column{?s}.",
					"i" = "Their additive effects will be pooled with the prior for that period."))
			}
		}
	}
	for (.t in seq_along(Y)) {
		.lame_check_all_na(Y[[.t]], if (!is.null(names(Y))) names(Y)[.t] else as.character(.t))
	}
	# binary/cbin threshold via 1*(y>0); warn when counts are collapsed
	if (family %in% c("binary", "cbin")) {
		uy <- unique(obs_y)
		if (any(!uy %in% c(0, 1))) {
			cli::cli_warn(c(
				"{.arg family} = {.val {family}} but {.arg Y} contains values other than 0/1.",
				"i" = "{.arg Y} will be thresholded to {.code 1 * (Y > 0)}; if you meant counts, use {.val poisson}, or {.val ordinal}/{.val normal} as appropriate."))
		}
	}
	if (family == "normal") {
		uy <- unique(obs_y)
		if (length(uy) <= 2L && all(uy %in% c(0, 1))) {
			cli::cli_inform(c(
				"i" = "{.arg Y} looks binary (values in {.val {{0, 1}}}) but {.arg family} = {.val normal}.",
				"i" = "If this is a 0/1 network, pass {.code family = \"binary\"}."))
		}
	}
	# warn when burn dwarfs nscan
	if (burn > nscan) {
		cli::cli_warn(c(
			"{.arg burn} = {burn} is greater than {.arg nscan} = {nscan}.",
			"i" = "Most iterations will be discarded; consider {.arg burn} <= {.arg nscan}."))
	}
	# ar(1) dynamics need at least two time slices to be identified
	if ((isTRUE(dynamic_uv) || isTRUE(dynamic_ab)) && length(Y) < 2L) {
		cli::cli_abort(c(
			"{.arg dynamic_uv}/{.arg dynamic_ab} = TRUE require at least 2 time periods.",
			"i" = "Got {length(Y)} slice{?s}; use {.fn ame} for a single cross-section, or supply more periods."))
	}
	# dynamic_beta needs at least two periods, like the other dynamics
	if (!is.null(dynamic_beta) &&
	    !(is.logical(dynamic_beta) && length(dynamic_beta) == 1L && isFALSE(dynamic_beta)) &&
	    length(Y) < 2L) {
		cli::cli_abort(c(
			"{.arg dynamic_beta} requires at least 2 time periods.",
			"i" = "Got {length(Y)} slice{?s}; use {.fn ame} for a single cross-section."))
	}
	# unknown prior-list keys usually mean a typo
	if (length(prior) > 0L && !is.null(names(prior))) {
		known_prior <- c("Sab0", "eta0", "etaab", "Suv0", "kappa0",
		                 "s20", "s2u0",
		                 "rho_uv_mean", "rho_uv_sd",
		                 "rho_ab_mean", "rho_ab_sd",
		                 "sigma_uv_shape", "sigma_uv_scale",
		                 "sigma_ab_shape", "sigma_ab_scale",
		                 "rho_beta_mean", "rho_beta_sd",
		                 "sigma_beta_shape", "sigma_beta_scale",
		                 "sigma_beta_init",
		                 "rho_beta_lower", "rho_beta_upper",
		                 "beta0_mean", "beta0_var",
		                 "matern32_length_scale",
		                 "snap_kappa", "snap_pi_a", "snap_pi_b", "t_nu")
		bad <- setdiff(names(prior), known_prior)
		if (length(bad) > 0L) {
			cli::cli_warn(c(
				"Unknown {.arg prior} list element{?s}: {.val {bad}}.",
				"i" = "Recognised names: {.val {known_prior}}.",
				"i" = "Note: {.code g} is a top-level argument to {.fn lame}, not a {.arg prior} entry."))
		}
	}

	if (isTRUE(periodic_save) && is.null(out_file)) {
		cli::cli_abort(c(
			"{.arg periodic_save} = TRUE requires an {.arg out_file} path.",
			"i" = "Supply {.arg out_file} or set {.arg periodic_save} = FALSE."))
	}

	#
	if( nscan %% odens !=0  ){ stop('"odens" must be a multiple of "nscan"')}

	# set random seed locally: restore the global rng stream on exit so a
	# downstream random draw is not silently perturbed by fitting a model
	if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
		.old_seed <- get(".Random.seed", envir = globalenv())
		on.exit(assign(".Random.seed", .old_seed, envir = globalenv()),
		        add = TRUE)
	}
	set.seed(seed)

	# lame() identifies actors by the row/column names of the y matrices.
	# auto-name matrices that carry no names (a plain simulated network would
		# otherwise collapse the actor set to length zero); assumes a fixed
	# actor composition across slices.
	y_named <- vapply(Y, function(y) {
		!is.null(rownames(y)) || !is.null(colnames(y))
	}, logical(1))
	if (!any(y_named)) {
			# give y and every covariate array the same default actor names
		nr0 <- nrow(Y[[1]]); nc0 <- ncol(Y[[1]])
		row_names_def <- if (bip) paste0("row", seq_len(nr0)) else paste0("a", seq_len(nr0))
		col_names_def <- if (bip) paste0("col", seq_len(nc0)) else row_names_def
		Y <- lapply(Y, function(y) {
			rownames(y) <- row_names_def; colnames(y) <- col_names_def; y
		})
		.apply_default_dimnames <- function(arr, row_nm, col_nm) {
			if (is.null(arr)) return(arr)
			d <- dim(arr)
			if (length(d) == 2L) {
				if (d[1] == length(row_nm)) rownames(arr) <- row_nm
				if (d[2] == length(col_nm)) colnames(arr) <- col_nm
			} else if (length(d) == 3L) {
				dn <- dimnames(arr); if (is.null(dn)) dn <- list(NULL, NULL, NULL)
				if (d[1] == length(row_nm) && is.null(dn[[1]])) dn[[1]] <- row_nm
				if (d[2] == length(col_nm) && is.null(dn[[2]])) dn[[2]] <- col_nm
				dimnames(arr) <- dn
			}
			arr
		}
		if (!is.null(Xdyad)) {
			Xdyad <- lapply(Xdyad, .apply_default_dimnames, row_names_def, col_names_def)
		}
		if (!is.null(Xrow)) {
			Xrow <- lapply(Xrow, .apply_default_dimnames, row_names_def,
			               colnames(Xrow[[1]]) %||% character(0))
		}
		if (!is.null(Xcol)) {
			Xcol <- lapply(Xcol, .apply_default_dimnames, col_names_def,
			               colnames(Xcol[[1]]) %||% character(0))
		}
	} else if (!all(y_named)) {
		cli::cli_abort(c(
			"Some {.arg Y} slices carry actor names and others do not.",
			"i" = "Give every slice row/column names, or leave them all unnamed (names are then generated automatically)."))
	} else {
		# y is named on every slice. if covariate arrays have no dimnames,
		# attach y's per-slice dimnames so list_to_array reorders x with the
		# same actor sort as y. without this, alphabetic-sorted y names that
		# do not match positional order (e.g. "a1","a10","a11",...,"a2") would
		# silently scramble the y/x correspondence and collapse coefficients.
		.inherit_dimnames <- function(X, Y, axes) {
			if (is.null(X)) return(X)
			lapply(seq_along(X), function(t) {
				arr <- X[[t]]
				if (is.null(arr)) return(arr)
				d <- dim(arr)
				if (length(d) == 2L) {
					if (is.null(rownames(arr)) && d[1] == nrow(Y[[t]]))
						rownames(arr) <- rownames(Y[[t]])
					if (is.null(colnames(arr)) && d[2] == ncol(Y[[t]]) &&
					    "col" %in% axes)
						colnames(arr) <- colnames(Y[[t]])
				} else if (length(d) == 3L) {
					dn <- dimnames(arr); if (is.null(dn)) dn <- list(NULL, NULL, NULL)
					if (is.null(dn[[1L]]) && d[1] == nrow(Y[[t]]))
						dn[[1L]] <- rownames(Y[[t]])
					if (is.null(dn[[2L]]) && d[2] == ncol(Y[[t]]) &&
					    "col" %in% axes)
						dn[[2L]] <- colnames(Y[[t]])
					dimnames(arr) <- dn
				}
				arr
			})
		}
		if (!is.null(Xdyad)) Xdyad <- .inherit_dimnames(Xdyad, Y, c("row", "col"))
		if (!is.null(Xrow))  Xrow  <- .inherit_dimnames(Xrow,  Y, "row")
		if (!is.null(Xcol))  Xcol  <- .inherit_dimnames(Xcol,  Y, "row")
	}

	# get actor info
	if(bip) {
		rowActorByYr <- lapply(Y, rownames)
		colActorByYr <- lapply(Y, colnames)
		rowActorSet <- sort(unique(unlist(rowActorByYr)))
		colActorSet <- sort(unique(unlist(colActorByYr)))
		nA <- length(rowActorSet)
		nB <- length(colActorSet)
		# bipartite row and column actor sets should be disjoint.
		# when they overlap, the row-side `apm["alice"]` and column-side
		# `bpm["alice"]` refer to two *different* actors and get assigned
		# distinct estimates under the same key, which a downstream user
		# reasonably interprets as the same actor's two roles.
		overlap <- intersect(rowActorSet, colActorSet)
		if (length(overlap) > 0L) {
			cli::cli_warn(c(
				"{.arg mode} = {.val bipartite} but {length(overlap)} actor label{?s} appear{?s/} on both the row and column sides: {.val {head(overlap, 5)}}{.if {length(overlap) > 5} (and {length(overlap) - 5} more)}.",
				"i" = "In bipartite networks the row and column actor sets should be disjoint; APM and BPM will carry the same label but refer to different actors.",
				"i" = "Rename to disambiguate (e.g. prefix row actors with {.val r_} and column actors with {.val c_}) or refit with {.code mode = \"unipartite\"} if the data is square."))
		}
		n <- nA
		actorSet <- rowActorSet
		actorByYr <- rowActorByYr
		} else {
			actorByYr <- lapply(Y, rownames)
			actorSet <- sort(unique(unlist( actorByYr )))
			n <- length(actorSet)
			nA <- nB <- n
		}
		.lame_presence_matrix <- function(actor_set, actors_by_period) {
			out <- matrix(FALSE, length(actor_set), length(actors_by_period),
			              dimnames = list(actor_set, names(Y)))
			for (tt in seq_along(actors_by_period)) {
				out[actors_by_period[[tt]], tt] <- TRUE
			}
			out
		}
		.lame_transition_eligible <- function(presence) {
			out <- matrix(FALSE, nrow(presence), ncol(presence),
			              dimnames = dimnames(presence))
			if (ncol(presence) > 1L) {
				out[, -1L] <- presence[, -1L, drop = FALSE] &
					presence[, -ncol(presence), drop = FALSE]
			}
			out
		}
		if (bip) {
			snap_eligible_u <- .lame_transition_eligible(
				.lame_presence_matrix(rowActorSet, rowActorByYr))
			snap_eligible_v <- .lame_transition_eligible(
				.lame_presence_matrix(colActorSet, colActorByYr))
		} else {
			snap_presence <- .lame_presence_matrix(actorSet, actorByYr)
			snap_eligible_u <- .lame_transition_eligible(snap_presence)
			snap_eligible_v <- snap_eligible_u
		}
		
		# reset odmax param
		odmax <- rep( max( unlist( lapply(Y, function(y){ apply(y>0, 1, sum, na.rm=TRUE)  }) ) ), n )
	
	# calc savepoints
	savePoints <- (burn:(nscan+burn))[(burn:(nscan+burn)) %% odens==0]
	savePoints <- savePoints[round(quantile(1:length(savePoints), probs=seq(save_interval,1,save_interval)))]
	
	# check formatting of input objects
	check_format(Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol)
	
	# set diag to na (only for square matrices)
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
		# preserve the per-slice covariate names from the first period;
		# list_to_array_bipartite has already aligned the row/col indices
		slice_names <- dimnames(Xdyad[[1]])[[3]]
		Xdyad_4d <- array(NA, dim = c(nA, nB, pd_tmp, length(Xdyad)),
		                  dimnames = list(NULL, NULL, slice_names, NULL))
		for(tt in seq_along(Xdyad)) {
			Xdyad_4d[,,,tt] <- Xdyad[[tt]]
		}
		Xdyad <- Xdyad_4d
		rm(Xdyad_4d)
	}
	
	# set g-prior parameter if not provided. `g` may be a scalar or a
	# per-coefficient vector; guard against the vector path.
	if(length(g) == 1L && is.na(g)) {
		# zellner's g-prior default. for a continuous outcome (normal / tobit)
		# the prior on beta must scale with the response variance, otherwise
		# the default over-shrinks the coefficients toward zero on any
		# unstandardised, large-magnitude y (a normal lame fit on y ~ 1e3
		# would return beta ~ 0 instead of the true slope). `max(1, var(y))`
		# keeps g = n^2 byte-identical for standardised / unit-scale outcomes
		# (var <= 1) while scaling it up for large-variance outcomes. this
		# mirrors the variance-aware default already used by ame() and makes
		# the two entry points agree. discrete / rank families (binary, cbin,
		# frn, rrl, poisson, ordinal) live on a fixed o(1) latent scale, so
		# the count-based n^2 default is left untouched for them.
		if (is.element(family, c("normal", "tobit"))) {
			vY <- stats::var(c(Y), na.rm = TRUE)
			if (!is.finite(vY) || vY <= 0) vY <- 1
			g <- n^2 * max(1, vY)
		} else {
			g <- n^2
		}
	}
	
	# force binary if binary family specified 
	if(is.element(family,c("binary","cbin"))) { Y<-1*(Y>0) }
	
	# clamp tobit values to zero (y is already a 3d array here)
	if(family=="tobit") {
		Y[Y < 0 & !is.na(Y)] <- 0
	}
	
	# observed and max outdegrees
	if(is.element(family,c("cbin","frn","rrl")) ){
		odobs<-apply(Y>0,c(1,3),sum,na.rm=TRUE)
		# reject odmax = 0 with observed nominations (inconsistent with the
		# data) and odmax = na (no meaningful constraint)
		if (!is.null(odmax)) {
			if (anyNA(odmax)) {
				cli::cli_abort(c(
					"{.arg odmax} contains {.val NA}.",
					"i" = "Supply a positive integer (scalar or per-row vector). Pass {.code odmax = NULL} to let {.fn lame} infer it from {.arg Y}."))
			}
			if (any(odmax <= 0L) && any(Y > 0, na.rm = TRUE)) {
				cli::cli_abort(c(
					"{.arg odmax} <= 0 but {.arg Y} contains nominations.",
					"i" = "{.code odmax = 0} prohibits any nominations; this is inconsistent with the observed data."))
			}
		}
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

	.lame_has_uv_rank <- function(bip, R, RA, RB) {
		if (isTRUE(bip)) return(isTRUE(RA > 0L) && isTRUE(RB > 0L))
		isTRUE(R > 0L)
	}
	if (isTRUE(dynamic_uv) && !.lame_has_uv_rank(bip, R, RA, RB)) {
		cli::cli_abort(c(
			"{.arg dynamic_uv = TRUE} requires a positive multiplicative latent rank.",
			"i" = if (bip) {
				"Set positive {.arg R_row} and {.arg R_col} values, or set {.arg R} > 0 for both."
			} else {
				"Set {.arg R} > 0 or disable {.arg dynamic_uv}."
			}))
	}
	if (isTRUE(dynamic_uv) && isTRUE(bip) &&
	    dynamic_uv_kind %in% c("snap", "t")) {
		cli::cli_abort(c(
			"{.arg dynamic_uv_kind = \"{dynamic_uv_kind}\"} is not available for this bipartite MCMC path.",
			"i" = "Use {.arg dynamic_uv_kind = \"ar1\"} here, or use {.code method = \"als\"} for supported bipartite snap/t point estimates."))
	}
	if(isTRUE(dynamic_G)) {
		if (!isTRUE(bip)) {
			cli::cli_abort(c(
				"{.arg dynamic_G = TRUE} only applies to bipartite networks.",
				"i" = "For unipartite networks the latent interaction matrix is the identity, so there is no {.code G_t} to estimate."))
		}
		if (!isTRUE(RA > 0L) || !isTRUE(RB > 0L)) {
			cli::cli_abort(c(
				"{.arg dynamic_G = TRUE} requires positive bipartite latent ranks.",
				"i" = "Set positive {.arg R_row} and {.arg R_col} values, or set {.arg R} > 0 for both."))
		}
		if (verbose) {
			cli::cli_inform(c(
				"i" = "{.code dynamic_G = TRUE}: running Carter-Kohn FFBS on {.code vec(G_t)} under an AR(1) state-space prior.",
				"i" = "Interpret {.code U_t %*% G_t %*% t(V_t)}, fitted values, and predictive diagnostics; raw {.code G_t} entries are gauge-dependent."))
		}
		dynamic_G_was_requested <- TRUE
	}
	if(dynamic_uv && verbose) {
		msg <- switch(dynamic_uv_kind,
			ar1 = c(
				"i" = "Fitting AR(1) time-varying multiplicative latent effects.",
				"i" = "Latent coordinates are gauge-dependent; interpret aligned trajectories, bilinear contributions, and predictive diagnostics."),
			snap = c(
				"i" = "Fitting snap-shift dynamic multiplicative latent effects.",
				"i" = "First-period snap probabilities are transition-ineligible and reported as NA."),
			t = c(
				"i" = "Fitting heavy-tailed Student-t dynamic multiplicative latent effects.",
				"i" = "Local transition scales diagnose unusually large movement; they are not discrete snap probabilities."))
		cli::cli_inform(msg)
	}

	# warn when the requested latent rank is too large for the network size
	{
		R_used <- max(RA, RB)
		n_use <- if (bip) min(nA, nB) else n
		n_third <- max(1L, floor(n_use / 3))
		if (R_used > n_third) {
			cli::cli_warn(c(
				"!" = "R = {.val {R_used}} is large relative to n = {.val {n_use}} (recommended cap ~ n/3 = {.val {n_third}}).",
				"i" = "High-rank latent factors can absorb structure that belongs to the additive effects, hurting interpretability.",
				"i" = "Consider {.code R = 2} or {.code R = 3} unless you have a specific reason."))
		}
	}
	
	# construct design matrix    
	pr<-length(Xrow[,,1])/n
	pc<-length(Xcol[,,1])/n
	
	# prepare default xlist for bipartite before design matrix construction
	if(bip) {
		# for bipartite networks
		pd <- if(!is.null(Xdyad)) length(Xdyad[,,,1])/(nA*nB) else 0
		if(!is.null(Xdyad) && length(Xdyad) > 0) {
			# convert xdyad (4d array: na x nb x pd x n) to xlist (list of
			# 3d arrays); preserve the dyadic slice names so downstream
			# coefficient labels carry the user's covariate names.
			Xlist <- vector("list", N)
			n_cov_d <- dim(Xdyad)[3]
			n_cov_total <- n_cov_d
			if(intercept) n_cov_total <- n_cov_total + 1L
			dyad_names <- dimnames(Xdyad)[[3]]
			if (is.null(dyad_names) || length(dyad_names) == 0L) {
				dyad_names <- paste0("x", seq_len(n_cov_d))
			}
			slice_names <- if (intercept) c("intercept", dyad_names) else dyad_names
			for(t in 1:N) {
				Xt <- array(0, dim = c(nA, nB, n_cov_total),
				            dimnames = list(NULL, NULL, slice_names))
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
			# intercept-only xlist
			if(intercept) {
				Xlist <- replicate(N, array(1, dim = c(nA, nB, 1L),
				                            dimnames = list(NULL, NULL, "intercept")),
				                   simplify = FALSE)
			} else {
				Xlist <- replicate(N, array(0, dim = c(nA, nB, 0L)), simplify = FALSE)
			}
		}
		# add row/column covariates to xlist if provided
		if(!is.null(Xrow) && length(Xrow) > 0) {
			n_xr <- if(length(dim(Xrow)) >= 2) dim(Xrow)[2] else 1L
				# carry the user's column names and append `_row` only when needed
			xr_names <- if (length(dim(Xrow)) >= 2) dimnames(Xrow)[[2]] else NULL
			if (is.null(xr_names) || length(xr_names) == 0L)
				xr_names <- paste0("xrow", seq_len(n_xr))
			xr_names <- .lame_apply_suffix(xr_names, "row")
			for(t in 1:N) {
				old_p <- dim(Xlist[[t]])[3]
				old_nm <- dimnames(Xlist[[t]])[[3]]
				new_Xt <- array(0, dim = c(nA, nB, old_p + n_xr),
				                dimnames = list(NULL, NULL, c(old_nm, xr_names)))
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
			xc_names <- if (length(dim(Xcol)) >= 2) dimnames(Xcol)[[2]] else NULL
			if (is.null(xc_names) || length(xc_names) == 0L)
				xc_names <- paste0("xcol", seq_len(n_xc))
			xc_names <- .lame_apply_suffix(xc_names, "col")
			for(t in 1:N) {
				old_p <- dim(Xlist[[t]])[3]
				old_nm <- dimnames(Xlist[[t]])[[3]]
				new_Xt <- array(0, dim = c(nA, nB, old_p + n_xc),
				                dimnames = list(NULL, NULL, c(old_nm, xc_names)))
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
		# bipartite: xlist already constructed above
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
	# use number of dyads like amen package, not p^2 (too restrictive). guard
	# against a per-coefficient vector for g (anyna on a vector returns false).
	if(length(g) == 1L && is.na(g)) {
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
	
	# dynamic uv needs the relevant multiplicative ranks
	.dyn_uv_has_rank <- isTRUE(dynamic_uv) && .lame_has_uv_rank(bip, R, RA, RB)

	# dynamic uv priors
	if(.dyn_uv_has_rank) {
		if(is.null(prior$rho_uv_mean)) { prior$rho_uv_mean <- 0.9 }
		if(is.null(prior$rho_uv_sd)) { prior$rho_uv_sd <- 0.1 }
		if(is.null(prior$sigma_uv_shape)) { prior$sigma_uv_shape <- 2 }
		if(is.null(prior$sigma_uv_scale)) { prior$sigma_uv_scale <- 1 }
		# snap-shift jump scale and beta prior on the snap rate
		if(identical(dynamic_uv_kind, "snap")) {
			if(is.null(prior$snap_kappa)) { prior$snap_kappa <- 2 }
			if(is.null(prior$snap_pi_a)) { prior$snap_pi_a <- 1 }
			if(is.null(prior$snap_pi_b)) { prior$snap_pi_b <- 9 }
		}
		if(identical(dynamic_uv_kind, "t")) {
			if(is.null(prior$t_nu)) { prior$t_nu <- 4 }
		}
	}
	
	# dynamic ab priors
	if(dynamic_ab) {
		if(is.null(prior$rho_ab_mean)) { prior$rho_ab_mean <- 0.8 }
		if(is.null(prior$rho_ab_sd)) { prior$rho_ab_sd <- 0.15 }
		if(is.null(prior$sigma_ab_shape)) { prior$sigma_ab_shape <- 2 }
		if(is.null(prior$sigma_ab_scale)) { prior$sigma_ab_scale <- 1 }
	}

	# dynamic beta priors. mirrors the dynamic_ab block above.
	if(is.null(prior$rho_beta_mean))    { prior$rho_beta_mean    <- 0.8  }
	if(is.null(prior$rho_beta_sd))      { prior$rho_beta_sd      <- 0.15 }
	if(is.null(prior$sigma_beta_shape)) { prior$sigma_beta_shape <- 2    }
	if(is.null(prior$sigma_beta_scale)) { prior$sigma_beta_scale <- 1    }
	if(is.null(prior$sigma_beta_init))  { prior$sigma_beta_init  <- 0.25 }
	if(is.null(prior$rho_beta_lower))   { prior$rho_beta_lower   <- 0    }
	if(is.null(prior$rho_beta_upper))   { prior$rho_beta_upper   <- 0.999}
	
	# prepare xlist before initialization
	if (!bip && (is.null(Xlist) || length(Xlist) == 0L)) {
		# for unipartite, create n x n intercept matrices if still needed
		Xlist <- replicate(N, array(1, dim = c(n, n, 1L)), simplify = FALSE)
	}
	
	# get starting values for mcmc
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
		# dynamic_g starts with one g per
		# period, seeded from the static g.
		# rho_g + sigma_g2 are the ar(1) hyperparameters for the
		# vec(g_t) ffbs state-space prior. defaults: rho_g = 0.7 (mild
		# persistence), sigma_g2 = 0.1 (modest innovation variance). both
		# are sampled inside the mcmc loop (logit mh for rho, ig conjugate
		# for sigma_g2).
		G_cube <- NULL
		rho_G_state <- 0.7
		sigma_G2_state <- 0.1
		if (isTRUE(dynamic_G) && !is.null(G) && length(G) > 0L) {
			G_cube <- array(0, dim = c(nrow(G), ncol(G), N))
			for (t in seq_len(N)) G_cube[, , t] <- G
		}
		sigma2_a <- startValsObj$sigma2_a
		sigma2_b <- startValsObj$sigma2_b
		s2 <- startValsObj$s2
		rho <- 0  # no dyadic correlation in bipartite (set to 0 not null)
		use_rho <- FALSE
		Sab <- diag(c(sigma2_a, sigma2_b))  # independent variances
	} else {
		# unipartite initialization
		startValsObj <- get_start_vals(start_vals,Y,family,xP=dim(X)[3],rvar,cvar,R,odmax=odmax)
		Z<-startValsObj$Z ; beta<-startValsObj$beta ; a<-startValsObj$a
		b<-startValsObj$b ; U<-startValsObj$U ; V<-startValsObj$V
		rho<-startValsObj$rho ; s2<-startValsObj$s2 ; Sab<-startValsObj$Sab
		G <- NULL  # no g matrix for unipartite
		
		# symmetric models don't use rho (dyadic correlation)
		if(symmetric) {
			rho <- 0
			use_rho <- FALSE
		} else {
			use_rho <- TRUE
		}
	}

	# ordinal_cutpoints == "explicit" state. y_int is the integer-recoded
	# ordinal array (1..k), alpha is the length-(k-1) free-cutpoint vector
	# anchored at alpha[1] = 0, and tau_prop is the robbins-monro rw
	# proposal sd in delta space (adapted in burn-in only). alpha storage
	# is allocated below alongside the beta matrix.
	use_explicit_cutpoints <- (family == "ordinal" &&
	                            identical(ordinal_cutpoints, "explicit"))
	if (use_explicit_cutpoints) {
		# explicit cutpoints are not available for bipartite fits; data-induced
		# cutpoints are the only ordinal sampler with a bipartite branch.
		if (bip) {
			cli::cli_warn(c(
				"{.arg ordinal_cutpoints} = {.val explicit} is not available for bipartite ordinal fits.",
				"i" = "Using {.val data_induced} instead."))
			use_explicit_cutpoints <- FALSE
		}
	}
	if (use_explicit_cutpoints) {
		ord_lvls <- sort(unique(c(Y[is.finite(Y)])))
		K_ord <- length(ord_lvls)
		if (K_ord < 2L) {
			cli::cli_abort(c(
				"{.arg ordinal_cutpoints} = {.val explicit} requires at least 2 distinct observed categories.",
				"i" = "Found {K_ord}; check {.arg Y} or use a different family."))
		}
		Y_int <- array(match(Y, ord_lvls), dim = dim(Y), dimnames = dimnames(Y))
		alpha <- .init_alpha_from_data(Y)
		# rw proposal sd in delta = log(diff(alpha)) space. robbins-monro
		# adaptation in burn-in only; capped [0.01, 2.0] to avoid pathologies.
		tau_prop <- 0.5
		tau_target <- 0.30  # mh acceptance target for the cowles update
		alpha_accept_n <- 0L
		alpha_accept_total <- 0L
		# alpha storage is allocated after beta storage below
	} else {
		alpha <- numeric(0)
		Y_int <- NULL
		tau_prop <- NA_real_
		K_ord <- 0L
	}

	# initialize dynamic ab parameters if needed
	if(dynamic_ab) {
		# determine correct dimensions for a and b
		n_a <- if(bip) nA else n
		n_b <- if(bip) nB else n

		# convert a and b to 2d matrices (n_a x t, n_b x t) if not already
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
		# initialize ar(1) parameters
		rho_ab <- if (!is.null(startValsObj$rho_ab)) startValsObj$rho_ab
		          else if (!is.null(prior$rho_ab_mean)) prior$rho_ab_mean
		          else 0.8
		sigma_ab <- if (!is.null(startValsObj$sigma_ab)) startValsObj$sigma_ab
		            else 0.1
	} else {
		a_mat <- NULL
		b_mat <- NULL
		rho_ab <- NULL
		sigma_ab <- NULL
	}
	
	# initialize dynamic uv parameters if needed
	if(.dyn_uv_has_rank) {
		# determine correct dimensions for u and v cubes
		n_u <- if(bip) nA else n
		n_v <- if(bip) nB else n
		R_u <- if(bip) RA else R
		R_v <- if(bip) RB else R
		# convert u and v to 3d arrays if not already
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
		# initialize ar(1) parameters. when neither start_vals nor a user prior
		# pins rho_uv, use 0.9 (the documented default prior mean) so the first
		# uv sweep enforces persistence; otherwise an unset prior collapses
		# rho_uv to null and the chain gets stuck at zero on the bipartite path.
		rho_uv <- if (!is.null(startValsObj$rho_uv)) startValsObj$rho_uv
		          else if (!is.null(prior$rho_uv_mean)) prior$rho_uv_mean
		          else 0.9
		sigma_uv <- if (!is.null(startValsObj$sigma_uv)) startValsObj$sigma_uv
		            else 0.1
		# snap rate initialized at the beta prior mean
		pi_snap <- if(identical(dynamic_uv_kind, "snap")) {
			prior$snap_pi_a / (prior$snap_pi_a + prior$snap_pi_b)
		} else { NULL }
	} else {
		U_cube <- NULL
		V_cube <- NULL
		rho_uv <- NULL
		sigma_uv <- NULL
		pi_snap <- NULL
	}
	
	rm(list=c('startValsObj','start_vals'))
	
	# helpful mcmc params
	symLoopIDs <- lapply(1:(nscan + burn), function(x){ rep(sample(1:nrow(U)),4) })  
	asymLoopIDs <- lapply(1:(nscan + burn), function(x){ sample(1:R) })  
	# track per-iteration sampler failures so the post-mcmc warning loop can
	# surface them. each slot is bumped from the relevant trycatch handler.
	tryErrorChecks<-list(
		s2=0, betaAB=0, rho=0, UV=0, Z=0, beta=0,
		rho_sigma_beta=0,   # dynamic_beta hyperparameter updates
		beta_static=0,      # static-block solve in dynamic_ab bipartite path
		pool=0,             # hierarchical pooling mh step
		G=0                 # g update
	)

	# bipartite mh proposal tracking
	if(bip && RA > 0 && RB > 0) {
		mh_sd_U <- 0.1
		mh_sd_V <- 0.1
		mh_acc_U <- 0L
		mh_att_U <- 0L
		mh_acc_V <- 0L
		mh_att_V <- 0L
	}
	iter <- 1

	# initialise per-actor slope state when active.
	# the default path (per_actor_active = false) skips this
	if (per_actor_active) {
		n_actors_init <- if (identical(dynamic_beta_per_actor, "row")) {
			if (bip) nA else n
		} else {
			if (bip) nB else n
		}
		theta_actor <- matrix(0, n_actors_init, N)
		theta_actor_prev <- theta_actor
		rho_actor <- 0.8
		sigma_actor_sq <- 0.1
		# storage for per-actor draws / summary
		n_stored_iters_pa <- nscan / odens
		if (keep_per_actor %in% c("draws", "auto")) {
			# memory check for auto mode
			pa_bytes <- 8 * n_stored_iters_pa * n_actors_init * N
			if (keep_per_actor == "auto" && pa_bytes > 250 * 1024^2) {
				cli::cli_warn(c(
					"!" = "{.code keep_per_actor = \"auto\"}: per-actor draws would require {.val {round(pa_bytes/1024^2, 1)}} MB.",
					"i" = "Falling back to summary mode. Pass {.code keep_per_actor = \"draws\"} to force full storage."))
				keep_per_actor <- "summary"
			} else {
				THETA_ACTOR <- array(NA_real_, dim = c(n_stored_iters_pa, n_actors_init, N))
			}
		}
		if (keep_per_actor == "summary") {
			theta_actor_mean <- matrix(0, n_actors_init, N)
			theta_actor_m2   <- matrix(0, n_actors_init, N)  # streaming second moment
			theta_actor_n_obs <- 0L
		}
		RHO_ACTOR <- numeric(n_stored_iters_pa)
		SIGMA_ACTOR <- numeric(n_stored_iters_pa)
	}

	# output items
	if(bip) {
		n_beta <- if(!is.null(Xlist) && length(Xlist) > 0) dim(Xlist[[1]])[3] else 1
		BETA <- matrix(nrow = nscan/odens, ncol = n_beta)
	} else {
		BETA <- matrix(nrow = nscan/odens, ncol = dim(X)[3] - pr*symmetric)
	}
	VC<-matrix(nrow=nscan/odens,ncol=5-3*symmetric)

	# storage for the dynamic_g hyperparameters. rho_g and sigma_g2 are
	# scalar per stored draw; g_cube draws are accumulated in the
	# g_cube_running holds the running sum below
	if (isTRUE(dynamic_G)) {
		RHO_G <- numeric(nscan/odens)
		SIGMA_G2 <- numeric(nscan/odens)
		# running sum for g_cube posterior mean; same shape as g_cube
		G_cube_sum <- if (!is.null(G_cube)) array(0, dim = dim(G_cube)) else NULL
		G_cube_sumsq <- if (!is.null(G_cube)) array(0, dim = dim(G_cube)) else NULL
		G_cube_n <- 0L
	} else {
		RHO_G <- NULL; SIGMA_G2 <- NULL
		G_cube_sum <- NULL; G_cube_sumsq <- NULL; G_cube_n <- 0L
	}

	# storage for the cowles explicit-cutpoint posterior. alpha has length
	# for k-1 with alpha[1] = 0 (identification anchor; not stored). free
	# cutpoints stored are alpha_2..alpha_{k-1}, i.e. k-2 columns. for
	# for k = 2, binary-equivalent ordinal has zero free cutpoints, so
	# alpha storage is suppressed.
	if (use_explicit_cutpoints && K_ord >= 3L) {
		ALPHA <- matrix(NA_real_, nrow = nscan/odens, ncol = K_ord - 2L,
		                dimnames = list(NULL, paste0("alpha_", seq.int(2L, K_ord - 1L))))
	} else {
		ALPHA <- NULL
	}

	# storage for dynamic parameters. store the innovation-sd chains
	# (sigma_ab / sigma_uv) alongside the persistence chains (rho_ab /
	# alongside rho_uv so the returned fit carries fit$sigma_ab / fit$sigma_uv. the
	# fields print.lame()'s "sigma_innov" column and the dynamic-effects
	# vignette both reference (square them to recover the innovation
	# variance; the stationary marginal variance is sigma^2 / (1 - rho^2)).
	if(dynamic_ab) {
		RHO_AB <- numeric(nscan/odens)
		SIGMA_AB <- numeric(nscan/odens)
	} else {
		RHO_AB <- NULL
		SIGMA_AB <- NULL
	}
	if(dynamic_uv) {
		RHO_UV <- numeric(nscan/odens)
		SIGMA_UV <- numeric(nscan/odens)
	} else {
		RHO_UV <- NULL
		SIGMA_UV <- NULL
	}
	if(dynamic_uv && identical(dynamic_uv_kind, "snap")) {
		PI_SNAP <- numeric(nscan/odens)
	} else {
		PI_SNAP <- NULL
	}
	
	# initialize uv posterior sums based on dynamic or static
	R_eff <- if(bip) max(RA, RB) else R
	# snap-indicator accumulators (sender and receiver sides)
	DELTA_SUM <- NULL
	DELTA_SUM_V <- NULL
	DELTA_COUNT <- NULL
	DELTA_COUNT_V <- NULL
	SNAP_DRAWS <- NULL
	SNAP_DRAWS_V <- NULL
	SNAP_DRAWS_MODE <- keep_snap_draws
	if(dynamic_uv && identical(dynamic_uv_kind, "snap") && R_eff > 0) {
		DELTA_SUM <- matrix(0, if(bip) nA else n, N)
		DELTA_SUM_V <- matrix(0, if(bip) nB else n, N)
		DELTA_COUNT <- matrix(0L, if(bip) nA else n, N)
		DELTA_COUNT_V <- matrix(0L, if(bip) nB else n, N)
		if (SNAP_DRAWS_MODE %in% c("draws", "chunked")) {
			snap_bytes <- 8 * (nscan / odens) * N *
				((if(bip) nA else n) + if (!symmetric) (if(bip) nB else n) else 0)
			if (identical(SNAP_DRAWS_MODE, "chunked")) {
				cli::cli_warn(c(
					"{.code keep_snap_draws = \"chunked\"} currently stores snap indicators as in-memory arrays.",
					"i" = "Estimated storage for this fit is {.val {round(snap_bytes / 1024^2, 1)}} MB.",
					"i" = "Use {.code keep_snap_draws = \"summary\"} for cell-level summaries without retained draw arrays."))
			}
		}
		if (SNAP_DRAWS_MODE %in% c("draws", "chunked")) {
			SNAP_DRAWS <- array(NA_real_, dim = c(nscan / odens, if(bip) nA else n, N))
			if (!symmetric) {
				SNAP_DRAWS_V <- array(NA_real_, dim = c(nscan / odens, if(bip) nB else n, N))
			}
		}
	} else if (!identical(SNAP_DRAWS_MODE, "none")) {
		cli::cli_warn(c(
			"{.arg keep_snap_draws} only applies to MCMC fits with {.code dynamic_uv = TRUE} and {.code dynamic_uv_kind = \"snap\"}.",
			"i" = "No snap-indicator draws will be stored for this fit."))
	}
	if(dynamic_uv && R_eff > 0) {
		if(bip) {
			UVPS <- array(0, dim=c(nA, nB, N))
			U_SUM <- array(0, dim=c(nA, RA, N))
			V_SUM <- array(0, dim=c(nB, RB, N))
			U_SUM_SQ <- array(0, dim=c(nA, RA, N))
		} else {
			UVPS <- array(0, dim=c(n, n, N))
			U_SUM <- array(0, dim=c(n, R_eff, N))
			V_SUM <- array(0, dim=c(n, R_eff, N))
			U_SUM_SQ <- array(0, dim=c(n, R_eff, N))
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
	
	# initialize aps and bps with correct dimensions for bipartite
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
		GOF[,,1] <- apply(Y, 3, gof_stats_bipartite, warn_square = FALSE)
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
			# bipartite: get names from xlist
			if(!is.null(Xlist) && length(Xlist) > 0) {
				p_names <- dimnames(Xlist[[1]])[[3]]
				if(is.null(p_names) || length(p_names) != ncol(BETA)) {
					nc <- ncol(BETA)
					n_extra <- nc - as.integer(intercept)
					p_names <- if(nc == 0L) {
						character(0)
					} else if(intercept && n_extra > 0) {
						c("intercept", paste0("x", seq_len(n_extra)))
					} else if(intercept) {
						"intercept"
					} else {
						paste0("x", seq_len(nc))
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
	
	# dynamic_beta: parse, initialise state, and build the design caches that
	# the ffbs sampler needs. nothing below this block is touched in the
	# default `dynamic_beta = false` branch -- beta_dyn$any is false and every
	# downstream branch falls through to the historical static path.
	beta_dyn <- parse_dynamic_beta(
		dynamic_beta = dynamic_beta,
		coef_names   = colnames(BETA),
		coef_block   = .lame_coef_blocks(colnames(BETA), Xrow, Xcol, Xdyad, intercept, symmetric, bip),
		intercept    = intercept,
		T            = N,
		family       = family,
		mode         = mode)
	# pre-build per-period long-format design matrices (n*n rows x p cols)
	# and the per-period xtx_per_t / n_dyads_per_t caches the ffbs needs.
	# only do this when dynamic_beta is active; the default path doesn't pay
	# for these allocations.
	if (beta_dyn$any) {
		p_beta <- ncol(BETA)
		Xdyn_list_per_t  <- vector("list", N)
		Xstat_list_per_t <- vector("list", N)
		XtX_per_t        <- vector("list", N)
		n_dyads_per_t    <- vector("list", N)
		# unipartite uses x[,,,t] -> n*n x p long; bipartite uses xlist[[t]] -> na*nb x p.
		# the actor-ordering matches what the rest of the loop assumes because
		# we vectorise the same arrays the loop is already using.
		for (t in seq_len(N)) {
			if (bip) {
				Xt <- Xlist[[t]]                      # na x nb x p
				n_long <- nA * nB
			} else {
				Xt <- X[,,,t]                          # n x n x p
				if (length(dim(Xt)) == 2L) Xt <- array(Xt, c(n, n, 1))
				n_long <- n * n
			}
			p_t <- dim(Xt)[3]
			Xmat <- matrix(0, n_long, p_t)
			for (k in seq_len(p_t)) Xmat[, k] <- as.numeric(Xt[,,k])
			if (length(beta_dyn$dynamic_idx) > 0L) {
				Xdyn_list_per_t[[t]]  <- Xmat[, beta_dyn$dynamic_idx, drop = FALSE]
			} else {
				Xdyn_list_per_t[[t]]  <- matrix(0, n_long, 0)
			}
			if (length(beta_dyn$static_idx) > 0L) {
				Xstat_list_per_t[[t]] <- Xmat[, beta_dyn$static_idx, drop = FALSE]
			} else {
				Xstat_list_per_t[[t]] <- matrix(0, n_long, 0)
			}
			XtX_per_t[[t]]     <- crossprod(Xmat)
			n_dyads_per_t[[t]] <- n_long
		}
		# innovation scales: one diag(p_dyn) matrix per period (currently constant)
		Lambda_list <- build_beta_state_scales(beta_dyn, XtX_per_t, n_dyads_per_t,
		                                       g = g, beta0_var = prior$beta0_var %||% 10)
		Lambda_full <- Lambda_list[[1]]
		# per-coef rho / sigma expansion (one entry per dynamic column, broadcast
		# from the per-block parameter via beta_dyn$groups[dynamic_idx])
		rho_beta_by_group   <- init_rho_beta(beta_dyn$group_names, prior)
		# when dynamic_beta_kind = "rw1", pin rho to 1 for every block and
		# skip its sampling in the mcmc loop. the ffbs forward filter
		# handles rho = 1 because the stationary-variance step clamps
		# `1 - rho^2` away from zero.
		if (identical(dynamic_beta_kind, "rw1") && length(rho_beta_by_group) > 0L) {
			rho_beta_by_group[] <- 1
		}
		sigma_beta_by_group <- init_sigma_beta(beta_dyn$group_names, prior)
		# counts how many sweeps hit the sigma_beta robustness cap or the
		# beta-path/static finite clip (a post-fit warning fires if these are
		# non-trivial -- the signal that the dynamic-beta chain was diverging,
		# typically near-separable sparse binary data or an unstandardised
		# large-magnitude outcome).
		beta_sigma_capped_hits <- 0L
		beta_path_clipped_hits <- 0L
		# integer group_id for the c++ ffbs / rho / sigma helpers
		dyn_groups <- beta_dyn$groups[beta_dyn$dynamic_idx]
		group_id   <- as.integer(match(dyn_groups, beta_dyn$group_names))
		rho_by_coef   <- rho_beta_by_group[group_id]
		sigma_by_coef <- sigma_beta_by_group[group_id]
		# ar(1) state prior at t = 0
		p_dyn <- length(beta_dyn$dynamic_idx)
		beta0_mean <- rep(prior$beta0_mean %||% 0, p_dyn)
		beta0_var  <- prior$beta0_var %||% 10
		beta0_cov  <- diag(beta0_var, p_dyn)
		# initial beta path: replicate the static draw (or zero) across periods,
		# then take a single ffbs-style pass to seed each period (cheap, but a
		# zero start works fine too because burn-in absorbs it)
		beta_dyn_path <- matrix(0, N, p_dyn)
		for (t in seq_len(N)) beta_dyn_path[t, ] <- beta[beta_dyn$dynamic_idx]
		# static piece: keep the existing draw
		beta_static <- beta[beta_dyn$static_idx]
		# 3d beta storage for dynamic beta; iter slot is [iter, p, t]
		BETA <- array(NA_real_, dim = c(nscan/odens, p_beta, N),
		              dimnames = list(NULL, colnames(BETA),
		                              if (!is.null(pdLabs)) pdLabs else paste0("t", seq_len(N))))
		# storage for the per-block ar(1) hyperparameters
		RHO_BETA   <- matrix(NA_real_, nrow = nscan/odens,
		                     ncol = beta_dyn$n_groups,
		                     dimnames = list(NULL, beta_dyn$group_names))
		SIGMA_BETA <- matrix(NA_real_, nrow = nscan/odens,
		                     ncol = beta_dyn$n_groups,
		                     dimnames = list(NULL, beta_dyn$group_names))
		# memory cap warning for large 3d beta storage
		est_bytes <- prod(dim(BETA)) * 8
		if (est_bytes > 1024^3) {
			cli::cli_warn(c(
				"{.arg dynamic_beta} will allocate {.val {round(est_bytes/1024^3, 2)}} GB for 3-D BETA storage ({nscan/odens} draws x {p_beta} coefs x {N} periods).",
				"i" = "Increase {.arg odens}, decrease {.arg nscan}, or pass {.arg dynamic_beta} as a subset of coefficients to bound this."))
		}
	} else {
		# defaults for variables that downstream branches may still reference;
		# none of these are touched when beta_dyn$any is false
		Xdyn_list_per_t <- NULL; Xstat_list_per_t <- NULL
		XtX_per_t <- NULL; n_dyads_per_t <- NULL
		Lambda_full <- NULL; Lambda_list <- NULL
		rho_beta_by_group <- numeric(0)
		sigma_beta_by_group <- numeric(0)
		rho_by_coef <- numeric(0); sigma_by_coef <- numeric(0)
		group_id <- integer(0)
		beta_dyn_path <- NULL; beta_static <- NULL
		beta0_mean <- NULL; beta0_cov <- NULL
		RHO_BETA <- NULL; SIGMA_BETA <- NULL
	}

	# log_lik storage. with save_log_lik = true the pointwise log-likelihood
	# is recorded at each stored mcmc iteration so that loo::loo() / waic()
	# can be computed afterwards. the in-memory matrix is
	# [nscan/odens x n_obs] where n_obs is the number of observed dyads
	# (counting per-period dyads for a longitudinal fit). the warning below
	# fires when the in-memory cap is hit.
	#
	# save_log_lik = "chunked" streams the values to per-column-chunk binary
	# files instead of allocating one big matrix. each chunk file holds, in
	# row-major within-chunk order, the iteration-by-column log-lik values
	# for the dyads in that column range. on readback, read_log_lik(fit)
	# does one readbin per chunk and stitches them.
	log_lik_mode <- if (identical(save_log_lik, "chunked")) "chunked"
	                else if (isTRUE(save_log_lik)) "memory"
	                else "off"
	if (log_lik_mode != "off") {
		n_obs_per_t <- vapply(seq_len(N), function(t_) sum(!is.na(Y[,,t_])), integer(1))
		n_obs_total <- sum(n_obs_per_t)
		n_stored    <- nscan / odens
		est_bytes   <- as.numeric(n_stored) * n_obs_total * 8
		if (log_lik_mode == "memory" && est_bytes > 1024^3) {
			cli::cli_warn(c(
				"{.arg save_log_lik} = TRUE will allocate {.val {round(est_bytes/1024^3, 2)}} GB for the log-lik matrix ({n_stored} draws x {n_obs_total} obs).",
				"i" = "Increase {.arg odens}, switch to {.code save_log_lik = \"chunked\"}, or set {.arg save_log_lik} = FALSE."))
		}
		obs_idx <- which(!is.na(Y), arr.ind = TRUE)  # rows = (i, j, t)
		if (log_lik_mode == "memory") {
			LOG_LIK <- matrix(NA_real_, nrow = n_stored, ncol = n_obs_total)
			LOG_LIK_CONNS <- NULL
			LOG_LIK_META  <- NULL
		} else {
			# chunked: open one binary connection per chunk
			LOG_LIK <- NULL
			ll_dir <- if (is.null(log_lik_path)) {
				tempfile(pattern = "lame_loglik_", tmpdir = tempdir())
			} else log_lik_path
			dir.create(ll_dir, showWarnings = FALSE, recursive = TRUE)
			chunk_size <- max(1L, as.integer(log_lik_chunk_size))
			n_chunks <- as.integer(ceiling(n_obs_total / chunk_size))
			chunk_starts <- seq.int(1L, n_obs_total, by = chunk_size)
			chunk_ends   <- pmin(chunk_starts + chunk_size - 1L, n_obs_total)
			chunk_files  <- file.path(ll_dir,
				sprintf("loglik_chunk_%03d.bin", seq_len(n_chunks)))
			LOG_LIK_CONNS <- lapply(chunk_files, function(p)
				file(p, open = "wb"))
			LOG_LIK_META <- list(
				path        = ll_dir,
				files       = chunk_files,
				chunk_size  = chunk_size,
				n_obs_total = n_obs_total,
				n_stored    = n_stored,
				chunk_starts = chunk_starts,
				chunk_ends   = chunk_ends,
				n_chunks    = n_chunks
			)
		}
	} else {
		LOG_LIK <- NULL
		LOG_LIK_CONNS <- NULL
		LOG_LIK_META  <- NULL
		obs_idx <- NULL
	}

	# randomised-halton shift for the ghk pointwise log-lik. deterministic
	# given the fit `seed` so re-runs produce identical loo summaries. the
	# shift is in [0, 1) and added per-digit to the halton base sequence per
	# the .halton_seq_randomised helper. (%||% is defined at the top of
	# lame() and is already in scope here.)
	.lame_log_lik_halton_shift <- if (identical(log_lik_method, "observed_ghk")) {
		# derive a per-seed shift via a small rng kick (does not perturb the
	# main mcmc stream because .random.seed is restored)
		old_seed <- if (exists(".Random.seed", envir = globalenv(),
		                       inherits = FALSE)) {
			get(".Random.seed", envir = globalenv())
		} else NULL
		set.seed(seed + 7919L)  # offset so it is not collinear with main seed
		shift <- stats::runif(1)
		if (!is.null(old_seed)) {
			assign(".Random.seed", old_seed, envir = globalenv())
		}
		shift
	} else 0

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
	
	# suppress armadillo warnings
	old_warn <- options()$warn
	if(symmetric) {
		options(warn = -1)  # suppress all warnings during mcmc
		on.exit(options(warn = old_warn), add = TRUE)
	}
	
	# pre-allocate residual array once
	E.nrm <- array(dim = dim(Z))

	if(burn!=0 && verbose){
		# only show progress if verbose=true
		cli::cli_h3("Starting burn-in period...")
		cli::cli_progress_bar("Burn-in", total = burn, .envir = environment())
	}
	# helper: per-period x*beta cube when dynamic_beta is on. this is added
	# to the ez returned by the static helpers, which are called with beta = 0.
	# when beta_dyn$any is false, the beta override and ez offset are unchanged.
	.lame_loop_beta <- function() {
		if (beta_dyn$any) numeric(ncol(BETA_static_view())) else beta
	}
	BETA_static_view <- function() {
		# returns a layout with length = full p for the zero-beta override
		# we wrap so beta_eff is freshly recomputed every iteration as the
		# dynamic-block path / static-block draw evolves.
		matrix(0, 1, length(beta))
	}
	.lame_beta_full_path <- function() {
		# returns a t x p matrix of per-period beta when dynamic_beta is on,
		# or null otherwise.
		if (!beta_dyn$any) return(NULL)
		out <- matrix(0, N, length(beta))
		for (t_ in seq_len(N)) {
			out[t_, ] <- assemble_beta_path(
				beta_dyn_path[t_, ], beta_static,
				beta_dyn$dynamic_idx, beta_dyn$static_idx,
				length(beta))
		}
		out
	}
	.lame_xbeta_cube <- function(beta_full_path) {
		if (is.null(beta_full_path)) return(NULL)
		n_a_ <- if (bip) nA else n
		n_b_ <- if (bip) nB else n
		out <- array(0, dim = c(n_a_, n_b_, N))
		for (t_ in seq_len(N)) {
			if (bip) Xt_ <- Xlist[[t_]] else Xt_ <- X[,,,t_]
			if (length(dim(Xt_)) == 2L) Xt_ <- array(Xt_, c(n_a_, n_b_, 1))
			p_t_ <- dim(Xt_)[3]
			b_t_ <- beta_full_path[t_, ]
			if (p_t_ > 0 && length(b_t_) >= p_t_) {
				for (k_ in seq_len(p_t_)) {
					if (b_t_[k_] != 0) out[,,t_] <- out[,,t_] + b_t_[k_] * Xt_[,,k_]
				}
			}
		}
		out
	}

	# snap indicators carried across scans for the storage gate
	delta_uv <- NULL
	delta_uv_v <- NULL
	# local scales carried across scans for the heavy-tailed dynamic uv kind
	if(dynamic_uv && identical(dynamic_uv_kind, "t")) {
		lambda_u <- matrix(1, n, N)
		lambda_v <- matrix(1, n, N)
	}

	for (s in 1:(nscan + burn))  {

		# max_seconds: break out cleanly if the wall clock is
		# past the user's budget. when terminating early we also force one
		# final checkpoint so a resumer always has a state to restart from
		# (otherwise a very-aggressive max_seconds could prevent the first
		# scheduled checkpoint write from ever firing).
		if (is.finite(max_seconds)) {
			elapsed <- as.numeric(difftime(Sys.time(), .lame_start_time, units = "secs"))
			if (elapsed > max_seconds) {
				if (verbose) cli::cli_inform(c(
					"i" = "Reached {.arg max_seconds} = {.val {max_seconds}}s at iteration {s} (of {nscan + burn}); terminating early."))
				.lame_terminated_early <- TRUE
				if (!is.null(checkpoint_path)) {
					.ckpt <- list(
						iter        = s,
						BETA        = BETA,
						VC          = VC,
						rng_state   = if (exists(".Random.seed")) .Random.seed else NULL,
						call_args   = mc,
						timestamp   = Sys.time()
					)
					try(saveRDS(.ckpt, file = checkpoint_path), silent = TRUE)
				}
				break
			}
		}
		# periodic checkpoint of the in-memory state so a long
		# run can be resumed via lame_resume(path). we snapshot the
		# minimal set of state variables needed for restart.
		if (!is.null(checkpoint_path) && s > 1L && (s %% checkpoint_every == 0L)) {
			.ckpt <- list(
				iter        = s,
				BETA        = BETA,
				VC          = VC,
				rng_state   = if (exists(".Random.seed")) .Random.seed else NULL,
				call_args   = mc,
				timestamp   = Sys.time()
			)
			# best-effort write; if it fails, continue the mcmc
			try(saveRDS(.ckpt, file = checkpoint_path), silent = TRUE)
		}

		# update z (e.nrm reused each iteration)

		Xlist_tmp <- Xlist

		# precompute beta_full_path and the xbeta-cube offset for the dynamic
		# branch. when dynamic_beta is off, both are null.
		beta_full_path_iter <- .lame_beta_full_path()
		Xbeta_cube_iter <- .lame_xbeta_cube(beta_full_path_iter)
		# `beta_for_ez` is what we pass to the existing ez helpers: when
		# dynamic_beta is on we zero it out so the helper sees no design
		# contribution (we add xbeta_cube_iter instead).
		beta_for_EZ <- if (beta_dyn$any) numeric(length(beta)) else beta

		if(dynamic_ab) {
			# use dynamic helper to compute ez with time-varying additive effects
			EZ <- get_EZ_dynamic_ab(Xlist_tmp, beta_for_EZ, a_mat, b_mat, U, V, N,
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
				if (!beta_dyn$any) {
					for(t in 1:N) {
						if(length(beta) > 0) {
							X_dims <- dim(Xlist_tmp[[t]])
							if(X_dims[1] != nA || X_dims[2] != nB) {
								stop(paste("Bipartite X dimension mismatch at time", t,
													": expected", nA, "x", nB,
													"but got", X_dims[1], "x", X_dims[2]))
							}
							# sum up x*beta for each covariate
							base_mat <- matrix(0, nA, nB)
							n_covs <- min(length(beta), X_dims[3])
							for(i in 1:n_covs) {
								base_mat <- base_mat + beta[i] * Xlist_tmp[[t]][,,i]
							}
							base_cube[,,t] <- base_mat
						}
					}
				}
				# dynamic_beta on: x*beta is precomputed in xbeta_cube_iter
				# (the existing static base_cube stays zero so we can add later)

				# call bipartite-specific function
				EZ <- get_EZ_bip_cpp(base_cube, a_mat, b_mat, U, V, G)
			} else {
				ab_mat <- outer(a, b,"+")
				storage.mode(ab_mat) <- "double"
				EZ <- get_EZ_cpp( Xlist_tmp, beta_for_EZ, ab_mat, U, V )
			}
		}
		# add per-period x*beta contribution when dynamic_beta is on (no-op
		# when beta_dyn$any is false because xbeta_cube_iter is null)
		if (beta_dyn$any && !is.null(Xbeta_cube_iter)) {
			EZ <- EZ + Xbeta_cube_iter
		}
		# batch c++ z sampling where possible
		if(family == "normal") {
			# single c++ call for all t time periods
			Z_batch <- try(rZ_nrm_batch_cpp(Z, EZ, rho, s2, Y), silent = TRUE)
			if(!inherits(Z_batch, 'try-error')) {
				Z <- Z_batch$Z
				E.nrm <- Z_batch$E_nrm
			} else {
				# fallback: per-t r sampling
				for(t in 1:N) { Z[,,t] <- rZ_nrm_fc(Z[,,t], EZ[,,t], rho, s2, Y[,,t]); E.nrm[,,t] <- Z[,,t] - EZ[,,t] }
				tryErrorChecks$Z <- tryErrorChecks$Z + 1
			}
		} else if(family == "tobit" && bip && rho == 0) {
				# bipartite tobit z sampling uses the batch c++ path when available
			Z_new <- tryCatch(rZ_tob_bip_batch_cpp(Z, EZ, s2, Y),
			                  error = function(e) { tryErrorChecks$Z <<- tryErrorChecks$Z + 1; NULL })
			if (!is.null(Z_new)) {
				Z <- Z_new
			} else {
				for(t in 1:N) Z[,,t] <- rZ_tob_fc(Z[,,t], EZ[,,t], 0, s2, Y[,,t])
			}
			E.nrm <- Z - EZ
		} else if(family == "binary" && bip && rho == 0) {
			# bipartite binary z sampling, same pattern as tobit above
			Z_new <- tryCatch(rZ_bin_bip_batch_cpp(Z, EZ, Y),
			                  error = function(e) { tryErrorChecks$Z <<- tryErrorChecks$Z + 1; NULL })
			if (!is.null(Z_new)) {
				Z <- Z_new
			} else {
				for(t in 1:N) Z[,,t] <- rZ_bin_fc(Z[,,t], EZ[,,t], 0, Y[,,t])
			}
			E.nrm <- Z - EZ
		} else {
			# fallback: per-t r loop for remaining families. when bip is true,
			# dispatch to the rectangular `_bip_fc` samplers in r/rz_bipartite.r,
			# which drop the upper/lower-triangle / diagonal / reciprocity
			# machinery and treat each cell as independent given ez. rho is
			# forced to 0 in the bipartite branches because the bipartite
			# linear predictor has no reciprocal-cell contribution.
			for(t in 1:N) {
				if(family == "tobit") {
					Z[,,t] <- rZ_tob_fc(Z[,,t], EZ[,,t], rho, s2, Y[,,t])
					E.nrm[,,t] <- Z[,,t] - EZ[,,t]
				}
				if(family == "binary") {
					Z[,,t] <- rZ_bin_fc(Z[,,t], EZ[,,t], rho, Y[,,t])
				}
				if(family == "ordinal") {
					if (bip) {
						Z[,,t] <- rZ_ord_bip_fc(Z[,,t], EZ[,,t], Y[,,t])
					} else if (use_explicit_cutpoints) {
					# cowles (1996) explicit-cutpoint path: truncated-normal
					# z conditional on the current alpha. symmetric branch
						# uses precision-2 / variance-1/2 on the upper triangle
						# and mirrors to the lower; asymmetric branch is rect.
						if (symmetric) {
							Z[,,t] <- rZ_ord_sym_explicit_fc(
								Z[,,t], EZ[,,t], Y_int[,,t], alpha)
						} else {
							Z[,,t] <- rZ_ord_explicit_fc(
								Z[,,t], EZ[,,t], Y_int[,,t], alpha)
						}
					} else {
						Z[,,t] <- rZ_ord_fc(Z[,,t], EZ[,,t], rho, Y[,,t])
					}
				}
				if(family == "cbin") {
					if (bip) {
						# odobs is [na, t]; pass the per-period row outdegree
						Z[,,t] <- rZ_cbin_bip_fc(Z[,,t], EZ[,,t], Y[,,t], odmax, odobs[, t])
					} else {
						Z[,,t] <- rZ_cbin_fc(Z[,,t], EZ[,,t], rho, Y[,,t], odmax, odobs)
					}
				}
				if(family == "frn") {
					if (bip) {
						Z[,,t] <- rZ_frn_bip_fc(Z[,,t], EZ[,,t], Y[,,t], YL[[t]], odmax, odobs[, t])
					} else {
						Z[,,t] <- rZ_frn_fc(Z[,,t], EZ[,,t], rho, Y[,,t], YL[[t]], odmax, odobs)
					}
				}
				if(family == "rrl") {
					if (bip) {
						Z[,,t] <- rZ_rrl_bip_fc(Z[,,t], EZ[,,t], Y[,,t], YL[[t]])
					} else {
						Z[,,t] <- rZ_rrl_fc(Z[,,t], EZ[,,t], rho, Y[,,t], YL[[t]])
					}
				}
				if(family == "poisson") {
					if (bip) {
						# rectangular mh on z; rho ≡ 0 in bipartite. carries the
						# per-period exposure offset on the bipartite path too.
						log_exp <- if (poisson_exposure_active) log(period_exposure[t]) else 0
						Z[,,t] <- rZ_pois_bip_fc(Z[,,t], EZ[,,t], s2, Y[,,t],
						                          log_exposure = log_exp)
					} else if (poisson_exposure_active) {
						# unipartite poisson with a per-period exposure offset
						Z[,,t] <- suppressWarnings(rZ_pois_fc_exposure_cpp(
							Z[,,t], EZ[,,t], rho, s2, Y[,,t],
							log_exposure = log(period_exposure[t])))
					} else {
						Z[,,t] <- rZ_pois_fc(Z[,,t], EZ[,,t], rho, s2, Y[,,t])
					}
					E.nrm[,,t] <- Z[,,t] - EZ[,,t]
				}
			}
		}
		
		# cowles (1996) mh update on the free ordinal cutpoints alpha,
		# z-marginalised. run after the per-t z sweep so ez is current.
		# robbins-monro adaptation on tau_prop fires in burn-in only.
		if (use_explicit_cutpoints) {
			alpha_step <- sample_alpha_cowles(
				alpha = alpha,
				Y_int = Y_int,
				EZ = EZ,
				tau_prop = tau_prop,
				symmetric = symmetric)
			alpha <- alpha_step$alpha
			alpha_accept_total <- alpha_accept_total + isTRUE(alpha_step$accept)
			alpha_accept_n <- alpha_accept_n + 1L
			# robbins-monro: log(tau) drift toward target acceptance, capped.
			# only adapt during burn-in; freeze afterwards for a stationary chain.
			# `s` is the outer mcmc iteration; freezing post-burn keeps the
			# chain a valid markov sampler with stationary distribution.
			if (s <= burn && alpha_accept_n >= 25L) {
				acc_rate <- alpha_accept_total / alpha_accept_n
				log_tau <- log(tau_prop) + (acc_rate - tau_target) / sqrt(s)
				tau_prop <- min(max(exp(log_tau), 0.01), 2.0)
				alpha_accept_n <- 0L
				alpha_accept_total <- 0L
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
		if (beta_dyn$any) {
			# === dynamic_beta path ====================================
			# 1) sample beta_dyn_path via ffbs (conditional on a, b, u, v)
			# 2) sample beta_static via conjugate gaussian (conditional on path)
			# 3) sample a, b via a conditional gibbs given the per-period beta
			# build the per-period offset for the ffbs update
			n_a_ <- if (bip) nA else n
			n_b_ <- if (bip) nB else n
			offset_list_ffbs <- vector("list", N)
			Z_list_ffbs      <- vector("list", N)
			# resolve per-period a, b, u, v matrices (build them from current
			# state so the ffbs sees the same a/b that the rest of the loop will)
			a_mat_eff <- if (dynamic_ab) a_mat else matrix(rep(a, N), n_a_, N)
			b_mat_eff <- if (dynamic_ab) b_mat else matrix(rep(b, N), n_b_, N)
			# dynamic_uv uses cubes; static uv is replicated across periods
			U_cube_eff <- if (.dyn_uv_has_rank) U_cube else {
				if (R > 0 || (bip && RA > 0)) {
					Uc <- array(0, dim = c(nrow(U), max(1L, ncol(U)), N))
					for (t in seq_len(N)) Uc[,,t] <- U
					Uc
				} else array(0, dim = c(n_a_, 1L, N))
			}
			V_cube_eff <- if (.dyn_uv_has_rank) V_cube else {
				if (R > 0 || (bip && RB > 0)) {
					Vc <- array(0, dim = c(nrow(V), max(1L, ncol(V)), N))
					for (t in seq_len(N)) Vc[,,t] <- V
					Vc
				} else array(0, dim = c(n_b_, 1L, N))
			}
			G_eff <- if (bip) G else diag(1, max(1L, dim(U_cube_eff)[2]))
			# offset = a + b + uv (no x*beta)
			for (t in seq_len(N)) {
				if (bip) {
					off_t <- matrix(0, nA, nB)
					off_t <- off_t + matrix(a_mat_eff[, t], nA, nB)
					off_t <- off_t + matrix(b_mat_eff[, t], nB, nA, byrow = FALSE) |> t() |> (function(m) matrix(rep(b_mat_eff[, t], each = nA), nA, nB))()
					# (the convoluted line above is awkward; do it cleanly)
					off_t <- outer(a_mat_eff[, t], b_mat_eff[, t], "+")
					if (R > 0 && nrow(G_eff) > 0 && ncol(G_eff) > 0) {
						Ut <- U_cube_eff[,,t]; Vt <- V_cube_eff[,,t]
						if (!is.matrix(Ut)) Ut <- matrix(Ut, n_a_, 1)
						if (!is.matrix(Vt)) Vt <- matrix(Vt, n_b_, 1)
						off_t <- off_t + Ut %*% G_eff %*% t(Vt)
					}
				} else {
					off_t <- outer(a_mat_eff[, t], b_mat_eff[, t], "+")
					if (R > 0) off_t <- off_t + U %*% t(V)
				}
				offset_list_ffbs[[t]] <- off_t
				Z_list_ffbs[[t]]      <- Z[,,t]
			}
			use_dyad_rho <- (!bip) && (!symmetric) && isTRUE(dcor) && (abs(rho) > 1e-12)
			# expand rho_by_coef / sigma_by_coef from per-block values
			rho_by_coef   <- rho_beta_by_group[group_id]
			sigma_by_coef <- sigma_beta_by_group[group_id]
			# sample the dynamic coefficient block
			.unequal_gaps <- !is.null(time_index) &&
				length(unique(diff(as.numeric(time_index)))) > 1L
			ffbs_out <- if (dynamic_beta_kind %in% c("rw2", "matern32") || .unequal_gaps) {
				tryCatch(
					.sample_beta_dynamic_alt(
						Xdyn_list   = Xdyn_list_per_t,
						Xstat_list  = Xstat_list_per_t,
						Z_list      = Z_list_ffbs,
						offset_list = offset_list_ffbs,
						beta_static = if (length(beta_static) > 0L) beta_static else numeric(0),
						sigma_by_coef = sigma_by_coef,
						Lambda      = Lambda_full,
						s2          = s2,
						kind        = dynamic_beta_kind,
						length_scale = prior$matern32_length_scale %||% NULL,
						time_positions = if (!is.null(time_index)) as.numeric(time_index) else NULL,
						rho_by_coef = rho_by_coef,
						dyad_rho    = if (is.null(rho)) 0 else rho,
						use_dyad_rho = use_dyad_rho,
						bipartite   = bip,
						symmetric   = symmetric),
					error = function(e) NULL)
			} else {
				tryCatch(
					sample_beta_dynamic_cpp(
						Xdyn_list  = Xdyn_list_per_t,
						Xstat_list = Xstat_list_per_t,
						Z_list     = Z_list_ffbs,
						offset_list = offset_list_ffbs,
						beta_static = if (length(beta_static) > 0L) beta_static else numeric(0),
						rho_by_coef = rho_by_coef,
						sigma_by_coef = sigma_by_coef,
						Lambda     = Lambda_full,
						beta0_mean = beta0_mean,
						beta0_cov  = beta0_cov,
						s2         = s2,
						dyad_rho   = if (is.null(rho)) 0 else rho,
						bipartite  = bip,
						symmetric  = symmetric,
						use_dyad_rho = use_dyad_rho),
					error = function(e) NULL)
			}
			if (!is.null(ffbs_out) && all(is.finite(ffbs_out$path))) {
				beta_dyn_path <- ffbs_out$path
				tryErrorChecks$beta <- tryErrorChecks$beta + (ffbs_out$chol_fail %||% 0L)
			} else {
				# fallback: keep the prior draw and bump the failure counter
				tryErrorChecks$beta <- tryErrorChecks$beta + 1L
			}
			# scale-aware finite clip on the sampled dynamic path. the
			# per-coefficient prior marginal sd is ~ sqrt(s2 * lambda_kk)
			# (lambda absorbs the covariate cross-product, s2 the response
			# scale), so a path excursion beyond 50 prior-sds is never real
			# signal; it is a diverging chain (ar(1) variance feedback, or
			# the probit latent-z <-> beta feedback on near-separable sparse
			# binary data). clipping keeps coef()/predict() finite instead of
			# returning +/-1e25; the bound is loose enough that a
			# well-identified fit is unchanged because it never binds. a
			# post-mcmc warning fires if the clip bound is hit.
			if (!is.null(Lambda_full) && ncol(beta_dyn_path) > 0L) {
				lam_diag <- pmax(diag(as.matrix(Lambda_full)), 1e-12)
				path_bound <- 50 * sqrt(max(s2, 1e-8) * lam_diag)
				for (j_ in seq_len(ncol(beta_dyn_path))) {
					bj <- path_bound[min(j_, length(path_bound))]
					if (any(abs(beta_dyn_path[, j_]) > bj)) {
						beta_path_clipped_hits <- beta_path_clipped_hits + 1L
						beta_dyn_path[, j_] <- pmin(pmax(beta_dyn_path[, j_], -bj), bj)
					}
				}
			}
			# 2) static-block conjugate update (if any static coefs)
			if (length(beta_dyn$static_idx) > 0L) {
				# prior precision = diag(1 / (g * s2)) to mirror the g-prior idiom
				p_static <- length(beta_dyn$static_idx)
				g_diag <- if (length(g) >= length(beta_dyn$static_idx))
					g[beta_dyn$static_idx] else rep(g[1], p_static)
				prior_prec <- diag(1 / (g_diag * s2), p_static)
				ffbs_stat <- tryCatch(
					sample_beta_static_cpp(
						Xdyn_list  = Xdyn_list_per_t,
						Xstat_list = Xstat_list_per_t,
						Z_list     = Z_list_ffbs,
						offset_list = offset_list_ffbs,
						beta_dyn_path = beta_dyn_path,
						prior_mean = rep(0, p_static),
						prior_prec = prior_prec,
						s2         = s2,
						dyad_rho   = if (is.null(rho)) 0 else rho,
						bipartite  = bip,
						symmetric  = symmetric,
						use_dyad_rho = use_dyad_rho),
					error = function(e) NULL)
				if (!is.null(ffbs_stat) && all(is.finite(ffbs_stat$beta))) {
					beta_static <- as.numeric(ffbs_stat$beta)
					tryErrorChecks$beta <- tryErrorChecks$beta + (ffbs_stat$chol_fail %||% 0L)
				} else {
					tryErrorChecks$beta <- tryErrorChecks$beta + 1L
				}
			# scale-aware finite clip on the static coefficients. prior
			# g-prior marginal sd for static coef k is sqrt(g_k * s2); a
				# value beyond 50 of those is a diverging chain (e.g. the static
				# intercept compensating for a near-separable binary latent z).
			# same rationale as the dynamic-path clip: keeps coef() finite
				# without binding on a real fit.
				stat_bound <- 50 * sqrt(g_diag * max(s2, 1e-8))
				for (j_ in seq_along(beta_static)) {
					bj <- stat_bound[min(j_, length(stat_bound))]
					if (is.finite(bj) && abs(beta_static[j_]) > bj) {
						beta_path_clipped_hits <- beta_path_clipped_hits + 1L
						beta_static[j_] <- min(max(beta_static[j_], -bj), bj)
					}
				}
			}
			# assemble the canonical p-length beta vector (time-mean of the path)
			# for downstream summaries; the per-period beta_full_path is what
			# actually feeds back into ez at the start of the next iteration.
			beta_path_now <- matrix(0, N, length(beta))
			for (t in seq_len(N)) {
				beta_path_now[t, ] <- assemble_beta_path(
					beta_dyn_path[t, ], beta_static,
					beta_dyn$dynamic_idx, beta_dyn$static_idx,
					length(beta))
			}
			beta <- colMeans(beta_path_now)
			# re-update a and b given the per-period beta path
			# build residual cube r[t] = z[t] - x*beta_full_path[t,] - uv[t]
			n_a__ <- if (bip) nA else n
			n_b__ <- if (bip) nB else n
			resid_cube <- array(0, dim = c(n_a__, n_b__, N))
			for (t in seq_len(N)) {
				if (bip) Xt_ <- Xlist[[t]] else Xt_ <- X[,,,t]
				if (length(dim(Xt_)) == 2L) Xt_ <- array(Xt_, c(n_a__, n_b__, 1))
				p_t_ <- dim(Xt_)[3]
				xb <- matrix(0, n_a__, n_b__)
				if (p_t_ > 0) {
					b_t_ <- beta_path_now[t, ]
					for (k in seq_len(p_t_)) if (b_t_[k] != 0) xb <- xb + b_t_[k] * Xt_[,,k]
				}
				uv_ <- if (R > 0) {
					Ut <- if (dynamic_uv) U_cube_eff[,,t] else U
					Vt <- if (dynamic_uv) V_cube_eff[,,t] else V
					if (bip) Ut %*% G_eff %*% t(Vt) else Ut %*% t(Vt)
				} else 0
				resid_cube[,,t] <- Z[,,t] - xb - uv_
			}
			# sample a and b from the residual conditional. resid_cube[t]
			# equals (a_i + b_j + noise_ij) for unipartite (with na diagonal)
			# or (a_i + b_j + noise_ij) for bipartite. use a standard
			# coordinate-gibbs over a then b (or a single joint draw pooled
			# across periods for the static case).
			# pass the coef_block labels (not beta_dyn$groups, which is masked by
			# the dynamic_idx) so the helper knows which block "intercept" is.
			ab_bases <- build_a_b_constraint_bases(
				beta_dyn,
				.lame_coef_blocks(colnames(BETA) %||% dimnames(BETA)[[2]],
				                  NULL, NULL, NULL, intercept, symmetric, bip),
				n_a__, n_b__, bip, symmetric,
				intercept = intercept)
			s2_floor <- max(s2, 1e-8)
			sab_a_floor <- max(Sab[1,1], 1e-8)
			sab_b_floor <- max(Sab[2,2], 1e-8)
			# helper to apply constraint projection to a/b
			.apply_H <- function(v, H) if (is.null(H)) v else as.numeric(H %*% (t(H) %*% v))
			if (!bip) {
				if (dynamic_ab) {
					for (t in seq_len(N)) {
						rt <- resid_cube[,,t]
						# unipartite: na diagonal. for row sums we exclude the diagonal.
						if (!symmetric) {
							# vectorised a_i update
							row_sums <- rowSums(rt, na.rm = TRUE)
							b_now <- b_mat[, t]
							# sum_{j != i}(b_j) = sum(b) - b_i
							b_sum  <- sum(b_now)
							row_clean <- row_sums - (b_sum - b_now)
							prec_a <- (n - 1) / s2_floor + 1 / sab_a_floor
							mean_a <- (row_clean / s2_floor) / prec_a
							a_mat[, t] <- rnorm(n, mean_a, 1 / sqrt(prec_a)) * rvar
							a_mat[, t] <- .apply_H(a_mat[, t], ab_bases$H_a)
							# b_j update
							col_sums <- colSums(rt, na.rm = TRUE)
							a_now <- a_mat[, t]
							a_sum <- sum(a_now)
							col_clean <- col_sums - (a_sum - a_now)
							prec_b <- (n - 1) / s2_floor + 1 / sab_b_floor
							mean_b <- (col_clean / s2_floor) / prec_b
							b_mat[, t] <- rnorm(n, mean_b, 1 / sqrt(prec_b)) * cvar
							b_mat[, t] <- .apply_H(b_mat[, t], ab_bases$H_b)
						} else {
							# symmetric: a == b; use both row and column information
							row_sums <- rowSums(rt, na.rm = TRUE)
							col_sums <- colSums(rt, na.rm = TRUE)
							a_now <- a_mat[, t]
							a_sum <- sum(a_now)
							# residual for actor i sums (n-1) row entries + (n-1) col entries
							clean <- row_sums + col_sums - 2 * (a_sum - a_now)
							prec  <- 2 * (n - 1) / s2_floor + 1 / sab_a_floor
							a_mat[, t] <- rnorm(n, (clean / s2_floor) / prec, 1 / sqrt(prec)) * nvar
							a_mat[, t] <- .apply_H(a_mat[, t], ab_bases$H_a)
							b_mat[, t] <- a_mat[, t]
						}
					}
					a <- a_mat[, 1]; b <- b_mat[, 1]
				} else {
					# static a, b pooled over periods
					row_clean_total <- numeric(n)
					col_clean_total <- numeric(n)
					for (t in seq_len(N)) {
						rt <- resid_cube[,,t]
						row_sums <- rowSums(rt, na.rm = TRUE)
						col_sums <- colSums(rt, na.rm = TRUE)
						b_sum <- sum(b); a_sum <- sum(a)
						row_clean_total <- row_clean_total + row_sums - (b_sum - b)
						col_clean_total <- col_clean_total + col_sums - (a_sum - a)
					}
					if (!symmetric) {
						prec_a <- (N * (n - 1)) / s2_floor + 1 / sab_a_floor
						prec_b <- (N * (n - 1)) / s2_floor + 1 / sab_b_floor
						a <- rnorm(n, (row_clean_total / s2_floor) / prec_a, 1 / sqrt(prec_a)) * rvar
						a <- .apply_H(a, ab_bases$H_a)
						b <- rnorm(n, (col_clean_total / s2_floor) / prec_b, 1 / sqrt(prec_b)) * cvar
						b <- .apply_H(b, ab_bases$H_b)
					} else {
						# pool row+col evidence
						clean <- row_clean_total + col_clean_total
						prec  <- (2 * N * (n - 1)) / s2_floor + 1 / sab_a_floor
						a <- rnorm(n, (clean / s2_floor) / prec, 1 / sqrt(prec)) * nvar
						a <- .apply_H(a, ab_bases$H_a)
						b <- a
					}
				}
			} else {
				# bipartite case (asymmetric by construction)
				if (dynamic_ab) {
					for (t in seq_len(N)) {
						rt <- resid_cube[,,t]
						# a_mat[,t] update
						row_sums <- rowSums(rt, na.rm = TRUE)
						b_now <- b_mat[, t]
						row_clean <- row_sums - sum(b_now) + 0  # this drops the b's
						# sufficient stat for a_i is sum_j (rt[i,j] - b_j)
						row_clean <- row_sums - sum(b_now)
						prec_a <- nB / s2_floor + 1 / sab_a_floor
						a_mat[, t] <- rnorm(nA, (row_clean / s2_floor) / prec_a, 1 / sqrt(prec_a)) * rvar
						a_mat[, t] <- .apply_H(a_mat[, t], ab_bases$H_a)
						# b_mat[,t] update
						col_sums <- colSums(rt, na.rm = TRUE)
						a_now <- a_mat[, t]
						col_clean <- col_sums - sum(a_now)
						prec_b <- nA / s2_floor + 1 / sab_b_floor
						b_mat[, t] <- rnorm(nB, (col_clean / s2_floor) / prec_b, 1 / sqrt(prec_b)) * cvar
						b_mat[, t] <- .apply_H(b_mat[, t], ab_bases$H_b)
					}
					a <- a_mat[, 1]; b <- b_mat[, 1]
				} else {
					# static a, b (bipartite): pool residuals
					row_total <- numeric(nA); col_total <- numeric(nB)
					for (t in seq_len(N)) {
						rt <- resid_cube[,,t]
						row_total <- row_total + rowSums(rt, na.rm = TRUE) - sum(b)
						col_total <- col_total + colSums(rt, na.rm = TRUE) - sum(a)
					}
					prec_a <- (nB * N) / s2 + 1 / max(Sab[1,1], 1e-8)
					prec_b <- (nA * N) / s2 + 1 / max(Sab[2,2], 1e-8)
					a <- rvar * rnorm(nA, (row_total / s2) / prec_a, 1 / sqrt(prec_a))
					b <- cvar * rnorm(nB, (col_total / s2) / prec_b, 1 / sqrt(prec_b))
					if (!is.null(ab_bases$H_a)) a <- ab_bases$H_a %*% (t(ab_bases$H_a) %*% a) |> as.numeric()
					if (!is.null(ab_bases$H_b)) b <- ab_bases$H_b %*% (t(ab_bases$H_b) %*% b) |> as.numeric()
				}
			}
			# update the per-block rho_beta and sigma_beta periodically. matches
			# the dynamic_ab cadence (every 20 iterations) to keep this cheap.
			# skip the rho update entirely when dynamic_beta_kind = "rw1"
			# (rho is pinned at 1 by definition; the ar(1) reduces to a random walk).
			# sigma_beta still updates conditional on the rw1 increments.
			is_rw1 <- identical(dynamic_beta_kind, "rw1")
			is_alt_kind <- dynamic_beta_kind %in% c("rw2", "matern32")
			if (s %% 20 == 0 && N >= 2 && length(beta_dyn$dynamic_idx) > 0L) {
				# lambda inverse for the c++ helpers
				Lambda_inv_full <- tryCatch(solve(Lambda_full),
				                            error = function(e) diag(1, nrow(Lambda_full)))
				# hierarchical pooling updates the shared per-block prior
				pool_sigma_scale  <- prior$sigma_beta_scale
				pool_rho_mean_vec <- rep(prior$rho_beta_mean, beta_dyn$n_groups)
				pool_rho_sd_vec   <- rep(prior$rho_beta_sd,   beta_dyn$n_groups)
				if (!identical(dynamic_beta_pool, "none") &&
				    beta_dyn$n_groups >= 2L) {
					if (dynamic_beta_pool %in% c("sigma", "both")) {
						pool_sigma_scale <- .pool_update_sigma_scale(
							sigma_beta_by_group,
							alpha   = prior$sigma_beta_shape,
							hyper_a = prior$sigma_beta_shape,
							hyper_b = prior$sigma_beta_scale)
					}
					if (dynamic_beta_pool %in% c("rho", "both")) {
						hyper <- tryCatch(
							.pool_update_rho_hyper(
								rho_beta_by_group,
								prior_mu_mean = 0,
								prior_mu_sd   = 10),
							error = function(e) {
								tryErrorChecks$pool <<- tryErrorChecks$pool + 1L
								NULL
							})
						if (!is.null(hyper)) {
							# transform pooled (mu, tau) on logit-rho back to
							# rho-prior-mean / sd on the (-1, 1) scale and use
							# as the per-block prior in the truncated-normal mh
							mu_r <- tanh(hyper$mu)
							sd_r <- max(hyper$tau, 1e-3)
							pool_rho_mean_vec <- rep(mu_r, beta_dyn$n_groups)
							pool_rho_sd_vec   <- rep(sd_r, beta_dyn$n_groups)
						}
					}
				}
				rho_new <- if (is_rw1 || is_alt_kind) {
					# rho is not meaningful for rw1/rw2/matern; report pinned 1
					rep(1, beta_dyn$n_groups)
				} else tryCatch(
					sample_rho_beta_cpp(
						beta_path  = beta_dyn_path,
						group_id   = as.integer(group_id),
						n_groups   = as.integer(beta_dyn$n_groups),
						Lambda_inv = Lambda_inv_full,
						sigma_by_coef = sigma_by_coef,
						rho_current   = rho_beta_by_group,
						rho_prior_mean = pool_rho_mean_vec,
						rho_prior_sd   = pool_rho_sd_vec,
						rho_lower = prior$rho_beta_lower,
						rho_upper = prior$rho_beta_upper),
					error = function(e) {
						tryErrorChecks$rho_sigma_beta <<- tryErrorChecks$rho_sigma_beta + 1L
						NULL
					})
				if (!is.null(rho_new) && all(is.finite(rho_new))) {
					rho_beta_by_group <- as.numeric(rho_new)
					names(rho_beta_by_group) <- beta_dyn$group_names
				}
				sigma_new <- if (is_alt_kind) {
					# rw2 and matern32 use kind-specific innovation residuals
					# quadratic form. aggregate to per-group by mean across the
					# coefs in each group (matches the existing per-group storage
					# format).
					sig_per_coef <- tryCatch(
						.sample_sigma_beta_alt(
							beta_path    = beta_dyn_path,
							kind         = dynamic_beta_kind,
							length_scale = prior$matern32_length_scale %||% NULL,
							prior_shape  = prior$sigma_beta_shape,
							prior_scale  = prior$sigma_beta_scale),
						error = function(e) {
							tryErrorChecks$rho_sigma_beta <<- tryErrorChecks$rho_sigma_beta + 1L
							NULL
						})
					if (is.null(sig_per_coef)) NULL else {
						# group-average: per-coef sqrt(var) collapsed to per-group
						per_group <- numeric(beta_dyn$n_groups)
						for (g_ in seq_len(beta_dyn$n_groups)) {
							idx_g <- which(group_id == g_)
							if (length(idx_g) > 0L)
								per_group[g_] <- sqrt(mean(sig_per_coef[idx_g]^2))
						}
						per_group
					}
				} else tryCatch(
					sample_sigma_beta_cpp(
						beta_path  = beta_dyn_path,
						group_id   = as.integer(group_id),
						n_groups   = as.integer(beta_dyn$n_groups),
						Lambda_inv = Lambda_inv_full,
						rho_by_group = rho_beta_by_group,
						prior_shape  = rep(prior$sigma_beta_shape, beta_dyn$n_groups),
						prior_scale  = rep(pool_sigma_scale, beta_dyn$n_groups)),
					error = function(e) {
						tryErrorChecks$rho_sigma_beta <<- tryErrorChecks$rho_sigma_beta + 1L
						NULL
					})
				if (!is.null(sigma_new) && all(is.finite(sigma_new))) {
					# scale-aware robustness cap on the ar(1) innovation sd.
					# the per-coefficient lambda scaling absorbs the covariate
					# cross-product but not the response scale s2 (the ffbs data
					# precision carries the 1/s2, the prior innovation precision
					# does not), so the natural magnitude of sigma_beta is
					# o(sqrt(s2)), not o(1). a fixed cap in raw units would
					# silently attenuate a legitimate dynamic coefficient on an
					# unstandardised large-magnitude outcome (e.g. normal y in
					# the thousands). capping at 10*sqrt(s2) is scale-invariant
					# in both x (via lambda) and y (via s2): it never binds on a
					# well-identified fit yet still catches the ar(1)
					# variance-feedback runaway on weakly-identified data. for
					# binary/probit s2 = 1 so the cap is 10 on the probit scale.
					sigma_beta_cap <- 10 * sqrt(max(s2, 1e-8))
					if (any(sigma_new > sigma_beta_cap)) {
						beta_sigma_capped_hits <- beta_sigma_capped_hits + 1L
					}
					sigma_new <- pmin(sigma_new, sigma_beta_cap)
					sigma_beta_by_group <- as.numeric(sigma_new)
					names(sigma_beta_by_group) <- beta_dyn$group_names
				}
			}
			# === end dynamic_beta path =================================
		} else if(dynamic_ab) {
			if(bip) {
				# bipartite dynamic_ab: update beta via bipartite gibbs
				# compute ez without a/b for residual construction
				zero_a <- matrix(0, nA, N)
				zero_b <- matrix(0, nB, N)
				EZ_no_ab <- get_EZ_dynamic_ab(Xlist, beta, zero_a, zero_b, U, V, N,
				                              bip = TRUE, G = G, nA = nA, nB = nB)

				# update beta via c++ (bipartite conjugate update)
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
						error = function(e) {
							tryErrorChecks$beta_static <<- tryErrorChecks$beta_static + 1L
							diag(s2, p_bip)
						}
					)
					m_post <- V_post %*% (Xty / s2)
					beta <- c(m_post + t(chol(V_post)) %*% rnorm(p_bip))
				}

				# recompute ez_no_ab after beta changes
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
				# no covariates: ez_no_ab is just uv
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
				# bipartite: a (na x t) and b (nb x t) have different dims,
				# so update each separately via conjugate normal gibbs steps
				var_innov <- sigma_ab^2
				for(t in 1:N) {
					resid_t <- Z[,,t] - EZ_no_ab[,,t]
					resid_t[is.na(resid_t)] <- 0

					# prior from ar(1)
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

			# update ar(1) parameters periodically (including during burn-in)
			if(s %% 20 == 0) {
				if(bip) {
					# bipartite: a (na x t) and b (nb x t) have different row counts,
					# so compute sufficient stats in r rather than c++ which assumes same n
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
					# bipartite dynamic_ab: same ig(sigma_ab_shape, sigma_ab_scale) prior
					# the unipartite c++ path uses. defaults match the historical "+1/+1".
					ig_shape <- prior$sigma_ab_shape %||% 2
					ig_scale <- prior$sigma_ab_scale %||% 1
					cnt <- (nA + nB) * (N - 1)
					ss_resid <- 0
					for(i in 1:nA) for(t in 2:N) ss_resid <- ss_resid + (a_mat[i,t] - rho_ab*a_mat[i,t-1])^2
					for(j in 1:nB) for(t in 2:N) ss_resid <- ss_resid + (b_mat[j,t] - rho_ab*b_mat[j,t-1])^2
					sigma_ab <- sqrt(1/rgamma(1, shape = ig_shape + cnt/2,
					                            rate  = ig_scale + ss_resid/2))
					# bipartite rho_ab still uses a conjugate normal full-conditional
					# (matching the rho_uv structure above)
					pmean_ab <- prior$rho_ab_mean %||% 0.8
					psd_ab   <- prior$rho_ab_sd   %||% 0.15
					prior_prec_ab <- 1 / (psd_ab * psd_ab)
					sigma2_inv_ab <- 1 / (sigma_ab * sigma_ab)
					sp_ab <- 0; ssl_ab <- 0
					for(i in 1:nA) for(t in 2:N) {
						sp_ab  <- sp_ab  + a_mat[i,t]   * a_mat[i,t-1]
						ssl_ab <- ssl_ab + a_mat[i,t-1] * a_mat[i,t-1]
					}
					for(j in 1:nB) for(t in 2:N) {
						sp_ab  <- sp_ab  + b_mat[j,t]   * b_mat[j,t-1]
						ssl_ab <- ssl_ab + b_mat[j,t-1] * b_mat[j,t-1]
					}
					var_post_ab  <- 1 / (ssl_ab * sigma2_inv_ab + prior_prec_ab)
					mean_post_ab <- var_post_ab * (sp_ab * sigma2_inv_ab +
					                               pmean_ab * prior_prec_ab)
					rho_ab <- max(-0.99, min(0.99,
					                          rnorm(1, mean_post_ab, sqrt(var_post_ab))))
				} else {
					# unipartite c++ path: forward user-set priors on rho_ab and sigma_ab
					pmean_ab <- prior$rho_ab_mean %||% 0.8
					psd_ab   <- prior$rho_ab_sd   %||% 0.15
					rho_ab <- sample_rho_ab_cpp(a_mat, b_mat, sigma_ab, rho_ab,
					                            symmetric,
					                            prior_mean = pmean_ab,
					                            prior_sd   = psd_ab)
					sigma_ab <- sample_sigma_ab_cpp(a_mat, b_mat, rho_ab, symmetric,
					                                prior_shape = prior$sigma_ab_shape %||% 2,
					                                prior_scale = prior$sigma_ab_scale %||% 1)
				}
			}

		} else {
			# standard static update
			if(bip) {
				# bipartite gibbs update via c++
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
				# handle na values in z before passing to c++
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
				# silent fallback when betaab calculation fails
				tryErrorChecks$betaAB<-tryErrorChecks$betaAB+1  
			}
		}
		
		# update sab using unified function
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
			# refresh per-period x*beta when dynamic_beta is on (beta or
			# beta_dyn_path may have moved since the iteration start)
			if (beta_dyn$any) {
				beta_full_path_iter <- .lame_beta_full_path()
				Xbeta_cube_iter     <- .lame_xbeta_cube(beta_full_path_iter)
			}
			beta_for_EZ <- if (beta_dyn$any) numeric(length(beta)) else beta
			if(dynamic_ab) {
				E.T <- Z - get_EZ_dynamic_ab(Xlist, beta_for_EZ, a_mat, b_mat, U, V, N,
				                             bip = bip, G = G, nA = nA, nB = nB)
			} else {
				E.T <- Z - get_EZ_cpp( Xlist, beta_for_EZ, outer(a, b,"+"), U, V )
			}
			if (beta_dyn$any && !is.null(Xbeta_cube_iter)) E.T <- E.T - Xbeta_cube_iter
			rhoNew<-try( rrho_mh_rep_cpp(E.T, rho,s2), silent=TRUE )
			if(!inherits(rhoNew, 'try-error')){ rho<-rhoNew } else { tryErrorChecks$rho<-tryErrorChecks$rho+1 }
		}
		
		# update u,v — bipartite uses ra / rb instead of r
		if ((!bip && R > 0) || (bip && (RA > 0 || RB > 0))) {
			# refresh per-period x*beta cube when dynamic_beta is on
			if (beta_dyn$any) {
				beta_full_path_iter <- .lame_beta_full_path()
				Xbeta_cube_iter     <- .lame_xbeta_cube(beta_full_path_iter)
			}
			beta_for_EZ <- if (beta_dyn$any) numeric(length(beta)) else beta
			if(dynamic_ab) {
				# for dynamic ab, use time-varying effects but zero out uv
				U_zero <- if(bip) array(0, dim=c(nA, max(1,RA), N)) else U*0
				V_zero <- if(bip) array(0, dim=c(nB, max(1,RB), N)) else V*0
				E <- Z - get_EZ_dynamic_ab(Xlist, beta_for_EZ, a_mat, b_mat, U_zero, V_zero, N,
				                           bip = bip, G = G, nA = nA, nB = nB)
			} else if(bip) {
					U_zero <- array(0, dim = c(nA, max(1, ncol(U)), N))
				V_zero <- array(0, dim = c(nB, max(1, ncol(V)), N))

				# create base cube with x*beta and a/b effects (zero when
				# dynamic_beta is on; the per-period contribution is added below)
				base_cube <- array(0, dim = c(nA, nB, N))
				if (!beta_dyn$any) {
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
				}

				# create time-expanded a and b matrices
				a_mat_temp <- matrix(a[1:nA], nA, N)
				b_mat_temp <- matrix(b[1:nB], nB, N)

				E <- Z - get_EZ_bip_cpp(base_cube, a_mat_temp, b_mat_temp, U_zero, V_zero, G)
			} else {
				E <- Z - get_EZ_cpp( Xlist, beta_for_EZ, outer(a, b,"+"), U*0, V*0 )
			}
			if (beta_dyn$any && !is.null(Xbeta_cube_iter)) E <- E - Xbeta_cube_iter
			shrink<- (s>.5*burn)
			
			if(dynamic_uv) {
				if(bip) {
					if(!identical(dynamic_uv_kind, "ar1")) {
						cli::cli_abort(c(
							"Internal dynamic-UV dispatch error.",
							"x" = "Bipartite dynamic UV can only use {.val ar1}; got {.val {dynamic_uv_kind}}."))
					}
					# bipartite dynamic uv via c++
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
						UV_try <- NULL  # ensure ar(1) guard treats this as failure
					}

					# ar(1) hyperparameter updates for bipartite: route through
					# the same c++ samplers used by the unipartite path so the
					# bipartite case picks up the user-set prior.
					if(!is.null(UV_try) && s %% 10 == 0) {
						pmean_uv <- prior$rho_uv_mean %||% 0.9
						psd_uv   <- prior$rho_uv_sd   %||% 0.1
						rho_uv <- sample_rho_uv(U_cube, V_cube, sigma_uv, rho_uv,
						                        symmetric = FALSE,
						                        prior_mean = pmean_uv,
						                        prior_sd   = psd_uv)
						sigma_uv <- sample_sigma_uv(U_cube, V_cube, rho_uv,
						                            symmetric = FALSE)
					}
					} else {
						# unipartite dynamic uv update: ar(1) drift or snap-shift
						if(identical(dynamic_uv_kind, "snap")) {
							if(is.null(delta_uv)) { delta_uv <- matrix(0, n, N) }
							if(is.null(delta_uv_v)) { delta_uv_v <- matrix(0, n, N) }
							UV <- try(
								rUV_dynamic_snap_fc_cpp(U_cube, V_cube, E, rho_uv, sigma_uv, s2,
									prior$snap_kappa, pi_snap, delta_uv, delta_uv_v,
									shrink, symmetric),
								silent = FALSE)
						} else if(identical(dynamic_uv_kind, "t")) {
						UV <- try(
							rUV_dynamic_t_fc_cpp(U_cube, V_cube, E, rho_uv, sigma_uv, s2,
								prior$t_nu, lambda_u, lambda_v, shrink, symmetric),
							silent = FALSE)
					} else {
						UV <- try(
							rUV_dynamic_fc_cpp(U_cube, V_cube, E, rho_uv, sigma_uv, s2, shrink, symmetric),
							silent = FALSE)
					}
					if(inherits(UV, 'try-error')){
						UV <- list(U=U_cube, V=V_cube)
						tryErrorChecks$UV<-tryErrorChecks$UV+1
					} else {
						U_cube <- UV$U
						V_cube <- UV$V
						U <- apply(U_cube, c(1,2), mean)
						V <- apply(V_cube, c(1,2), mean)
						# snap indicators and beta-conjugate snap-rate update
						if(identical(dynamic_uv_kind, "snap")) {
							delta_uv <- UV$delta_u
							delta_uv_v <- UV$delta_v
							n_snap <- sum(delta_uv[,-1]) +
								if(symmetric) 0 else sum(delta_uv_v[,-1])
							n_elig <- (if(symmetric) 1 else 2) * nrow(delta_uv) * (N - 1)
							pi_snap <- stats::rbeta(1, prior$snap_pi_a + n_snap,
								prior$snap_pi_b + (n_elig - n_snap))
						}
						if(identical(dynamic_uv_kind, "t")) {
							lambda_u <- UV$lambda_u
							lambda_v <- UV$lambda_v
						}
					}

					# update ar(1) parameters (only if uv sampling succeeded). forward
					# the user-set normal(rho_uv_mean, rho_uv_sd^2) prior to the c++
					# sampler instead of the historical hard-coded n(0,1).
						if(!inherits(UV, 'try-error') && s %% 10 == 0) {
							pmean_uv <- prior$rho_uv_mean %||% 0.9
							psd_uv   <- prior$rho_uv_sd   %||% 0.1
							if(identical(dynamic_uv_kind, "snap") && !is.null(delta_uv)) {
								rho_uv <- .lame_sample_rho_uv_drift_only(
									U_cube, V_cube, sigma_uv, rho_uv, symmetric,
									delta_uv, delta_uv_v,
									prior_mean = pmean_uv, prior_sd = psd_uv
								)
								sigma_uv <- .lame_sample_sigma_uv_drift_only(
									U_cube, V_cube, rho_uv, sigma_uv, symmetric,
									delta_uv, delta_uv_v,
									prior$sigma_uv_shape, prior$sigma_uv_scale
								)
							} else if(identical(dynamic_uv_kind, "t") && !is.null(lambda_u)) {
								rho_uv <- .lame_sample_rho_uv_t_weighted(
									U_cube, V_cube, sigma_uv, rho_uv, symmetric,
									lambda_u, lambda_v,
									prior_mean = pmean_uv, prior_sd = psd_uv
								)
								sigma_uv <- .lame_sample_sigma_uv_t_weighted(
									U_cube, V_cube, rho_uv, sigma_uv, symmetric,
									lambda_u, lambda_v,
									prior$sigma_uv_shape, prior$sigma_uv_scale
								)
							} else {
								rho_uv <- sample_rho_uv(U_cube, V_cube, sigma_uv, rho_uv,
								                        symmetric,
								                        prior_mean = pmean_uv,
								                        prior_sd   = psd_uv)
								sigma_uv <- sample_sigma_uv(U_cube, V_cube, rho_uv, symmetric)
							}
						}
				}
			} else {
				# standard static uv update
				if(bip) {
					# bipartite mh random walk update for u and v
					if(RA > 0 && RB > 0) {
						# drop=false: `u[,,1]` with ra=1 collapses a 3-d
						# `[na, 1, t]` cube to a length-na vector, then
						# `u_2d[i,]` errors. force a matrix.
						U_2d <- if(length(dim(U)) == 3)
							U[, , 1, drop = FALSE] else U
						V_2d <- if(length(dim(V)) == 3)
							V[, , 1, drop = FALSE] else V
						if (length(dim(U_2d)) == 3) dim(U_2d) <- dim(U_2d)[1:2]
						if (length(dim(V_2d)) == 3) dim(V_2d) <- dim(V_2d)[1:2]

						# update u rows
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

						# update v rows
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

						# expand back to 3d for get_ez_bip_cpp compatibility
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

			# update g for bipartite (interaction matrix between row and column latent spaces)
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
				if (isTRUE(dynamic_G)) {
					# carter-kohn ffbs on vec(g_t) under an ar(1) state-space
					# prior with hyperparameters (rho_g, sigma_g2). this
					# couples g_t across time and gives a coherent state-space
					# estimate, with rho_g + sigma_g2 sampled inside the loop
					# (logit mh and ig conjugate, resp.).
					if (!exists("G_cube", inherits = FALSE) || is.null(G_cube)) {
						G_cube <- array(0, dim = c(RA, RB, N))
						for (t in seq_len(N)) G_cube[, , t] <- G
					}
					ffbs_out <- tryCatch(
						ffbs_vecG(E_cube = E, U_cube = U_3d, V_cube = V_3d,
							s2 = s2, rho_G = rho_G_state,
							sigma_G2 = sigma_G2_state),
						error = function(e) {
							tryErrorChecks$G <<- tryErrorChecks$G + 1L
							NULL
						})
					if (!is.null(ffbs_out)) {
						G_cube <- ffbs_out$G_cube
						# sample sigma_g2 conjugate ig, with the stationary-
						# variance scale-identification cap (tied to the
						# observation variance s2; for binary s2 = 1). this
						# stops the u-g-v scale degeneracy that otherwise
						# inflates g_cube by orders of magnitude.
						sigma_G2_state <- sample_sigma_G2(
							vecG_path = ffbs_out$vecG_path,
							rho_G = rho_G_state,
							s2_obs = s2)
						# sample rho_g via fisher-z logit mh
						rho_step <- sample_rho_G_mh(
							rho_G = rho_G_state,
							vecG_path = ffbs_out$vecG_path,
							sigma_G2 = sigma_G2_state,
							tau = 0.3)
						rho_G_state <- rho_step$rho
					}
					# working g uses the posterior-mean g_t across time for the
					# next iteration's ez rebuild. the full cube is stored for
					# reporting.
					G <- apply(G_cube, c(1, 2), mean)
				} else {
					G_new <- tryCatch(
						sample_G_bip_cpp(E, U_3d, V_3d, lambdaG = 1.0, s2 = rep(s2, N)),
						error = function(e) {
							tryErrorChecks$G <<- tryErrorChecks$G + 1L
							G
						}
					)
					G <- G_new
				}
			}
		}

		# per-actor pairwise contrast ffbs sweep. only runs when
		# dynamic_beta_per_actor is non-null.
		#
		# ez at this point reflects beta + a + b + uv' + g but does not
		# include the per-actor contribution; the actor term is added on
		# top by this block and overwritten by ez's rebuild at the top of
		# the next iter, so r_t = z[,,t] - ez[,,t] is the correct residual
		# for the per-actor update.
		if (per_actor_active) {
			n_a_pa <- if (bip) nA else n
			n_b_pa <- if (bip) nB else n
			n_actors_pa <- if (identical(dynamic_beta_per_actor, "row"))
				n_a_pa else n_b_pa
			# the internal xlist prepends the intercept (column of 1's) as
			# slot 1 when intercept = true; per_actor_covariate_idx is
			# 1-based into the user's xdyad. offset by +1 when intercept
			# is present so the index picks the intended covariate.
			cov_offset <- if (isTRUE(intercept)) 1L else 0L
			cov_idx <- min(per_actor_covariate_idx + cov_offset,
				if (length(Xlist) > 0L && length(dim(Xlist[[1L]])) >= 3L)
					dim(Xlist[[1L]])[3L] else 1L)
			# per-period sufficient stats
			H_mat <- matrix(0, n_actors_pa, N)
			h_mat <- matrix(0, n_actors_pa, N)
			for (t_pa in seq_len(N)) {
				Xt <- Xlist[[t_pa]]
				X_slice_t <- if (length(dim(Xt)) >= 3L && dim(Xt)[3L] >= cov_idx)
					Xt[, , cov_idx] else matrix(0, n_a_pa, n_b_pa)
				R_t <- Z[, , t_pa] - EZ[, , t_pa]
				suff <- .per_actor_suffstats_period(R_t, X_slice_t,
					dynamic_beta_per_actor, s2)
				H_mat[, t_pa] <- suff$H
				h_mat[, t_pa] <- suff$h
			}
			if (identical(per_actor_identifiability, "exact_center")) {
				# exact-center path: sample each actor's path from its
				# conditional ar(1) posterior via a 1-d carter-kohn ffbs
				# (n_actors independent univariate sweeps), then project to
				# the sum-to-zero manifold per period using variance-weighted
				# lagrangian conditioning. exact under the block-diagonal-
				# across-actors covariance of the lame per-actor model.
				theta_actor <- .sweep_per_actor_exact(
					H_mat = H_mat, h_mat = h_mat,
					rho_actor = rho_actor,
					sigma_actor2 = sigma_actor_sq)
			} else {
				# centered path: build ar(1) precision and run a pairwise
				# contrast sweep. preserves zero-sum exactly per draw.
				Q_AR1 <- .ar1_path_precision(N, rho_actor, sigma_actor_sq)
				theta_actor <- sweep(theta_actor, 2L, colMeans(theta_actor), "-")
				theta_actor <- .sweep_per_actor(theta_actor, H_mat, h_mat, Q_AR1)
			}
			# hyperparameter updates every 20 iterations (skip during the
			# first burn-in period to avoid the feedback loop where a noisy
			# first theta inflates sigma_actor and lets theta grow unbounded)
			if (s > burn / 2 && s %% 20 == 0) {
				sigma_actor_sq <- .sample_sigma_actor(theta_actor, rho_actor)
				rho_actor <- .sample_rho_actor(theta_actor, sigma_actor_sq, rho_actor)
			}
		}

		# burn-in countdown
		if(burn!=0 && s <= burn && verbose){
			# only update progress if verbose=true
			cli::cli_progress_update()
		}
		
		# store parameter values and monitor the mc
		if(s==burn+1 && verbose){
			# only show messages if verbose=true. when burn = 0 (typical for
			# the resume path), the burn progress bar was never created, so
			# the cli_progress_done() / "burn-in period complete" message is
			# also gated on burn != 0; the sampling bar must still come up so
			# the cli_progress_update() at line ~3650 has a current bar.
			if (burn != 0) {
				cli::cli_progress_done()
				cli::cli_alert_success("Burn-in period complete")
			}
			cli::cli_progress_bar("Sampling", total = nscan, .envir = environment())
		}
		if(s%%odens==0 & s>burn)  {
			
			# assemble the per-period full beta path when dynamic_beta is on
			# (used for both 3-d beta storage and downstream summaries)
			beta_path_now_store <- if (beta_dyn$any) {
				bp <- matrix(0, N, length(beta))
				for (t in seq_len(N)) {
					bp[t, ] <- assemble_beta_path(
						beta_dyn_path[t, ], beta_static,
						beta_dyn$dynamic_idx, beta_dyn$static_idx,
						length(beta))
				}
				bp
			} else NULL

			# store beta and vc - symmetric case
			if(symmetric){
				br<-beta[rb] ; bc<-beta[cb] ; bn<-(br+bc)/2
						if(intercept) {
					sbeta<-c(beta[1],bn,beta[-c(1,rb,cb)])
				} else {
					sbeta<-c(bn,beta[-c(rb,cb)])
				}
				if (beta_dyn$any) {
					# 3-d beta storage: per-period symmetric collapse. apply the
					# same rb/cb/bn merge to each period's beta_full row.
					for (t in seq_len(N)) {
						b_t <- beta_path_now_store[t, ]
						brt <- b_t[rb]; bct <- b_t[cb]; bnt <- (brt + bct)/2
						if (intercept) {
							st <- c(b_t[1], bnt, b_t[-c(1, rb, cb)])
						} else {
							st <- c(bnt, b_t[-c(rb, cb)])
						}
						BETA[iter, , t] <- st
					}
				} else {
					BETA[iter,]<-sbeta
				}
				VC[iter,]<-c(Sab[1,1],s2)
			}

			# store beta and vc - asymmetric case
			if(!symmetric){
				if (beta_dyn$any) {
					BETA[iter, , ] <- t(beta_path_now_store)  # p x t slot
				} else {
					BETA[iter,]<-beta
				}
				VC[iter,]<- c(Sab[upper.tri(Sab, diag = T)], rho,s2)
			}

			# store dynamic parameters
			if(dynamic_ab && !is.null(rho_ab)) {
				RHO_AB[iter] <- rho_ab
				if (!is.null(sigma_ab)) SIGMA_AB[iter] <- sigma_ab
			}
			if(dynamic_uv && !is.null(rho_uv)) {
				RHO_UV[iter] <- rho_uv
				if (!is.null(sigma_uv)) SIGMA_UV[iter] <- sigma_uv
			}
			if(dynamic_uv && identical(dynamic_uv_kind, "snap") && !is.null(pi_snap)) {
				PI_SNAP[iter] <- pi_snap
			}
			# store dynamic_beta hyper-parameters
			if (beta_dyn$any) {
				RHO_BETA[iter, ]   <- rho_beta_by_group
				SIGMA_BETA[iter, ] <- sigma_beta_by_group
			}
			# store explicit ordinal cutpoints alpha_2..alpha_{k-1}. the k=2
			# case has zero free cutpoints (alpha null); skip the store.
			if (use_explicit_cutpoints && !is.null(ALPHA)) {
				ALPHA[iter, ] <- alpha[-1L]
			}
			# store dynamic_g hyperparameters + accumulate g_cube draws
			if (isTRUE(dynamic_G)) {
				RHO_G[iter] <- rho_G_state
				SIGMA_G2[iter] <- sigma_G2_state
				if (!is.null(G_cube_sum) && !is.null(G_cube)) {
					G_cube_sum <- G_cube_sum + G_cube
					G_cube_sumsq <- G_cube_sumsq + G_cube^2
					G_cube_n <- G_cube_n + 1L
				}
			}
			# log-likelihood at this draw (pointwise per observed dyad).
			# ez is the linear predictor from the top of the iteration; pair
			# with the current s2 / rho to compute the per-cell density that
			# loo::loo() / waic() consume via loo.lame() / loo.ame().
			if (!is.null(LOG_LIK) || !is.null(LOG_LIK_CONNS)) {
				n_obs_iter <- if (!is.null(LOG_LIK)) ncol(LOG_LIK) else LOG_LIK_META$n_obs_total
				ll_row <- numeric(n_obs_iter)
				# extract observed y values and matching ez values in the same order
				y_obs <- Y[obs_idx]
				ez_obs <- EZ[obs_idx]
				# observed_ghk routes to the ghk helper, which uses
				# closed-form for {normal, binary, cbin, tobit, poisson},
				# closed-form cumulative-probit for ordinal, and halton-ghk
				# for {frn, rrl}.
				if (identical(log_lik_method, "observed_ghk")) {
					# route to the ghk + gauss-hermite + dim-capped pipeline.
					# halton_shift is captured below from
					# fit$log_lik_meta if it exists; otherwise default 0.
					ll_row <- .pointwise_loglik_observed_ghk_rank(
						y_obs = y_obs, ez_obs = ez_obs,
						z_obs = Z[obs_idx], obs_idx = obs_idx,
						family = family, s2 = s2,
						dyad_rho = if (is.null(rho)) 0 else rho,
						halton_shift = .lame_log_lik_halton_shift %||% 0)
				} else if (family == "normal") {
					ll_row <- dnorm(y_obs, mean = ez_obs, sd = sqrt(s2), log = TRUE)
				} else if (family == "binary" || family == "cbin") {
					# probit link: p(y=1) = pnorm(ez)
					p_hat <- pnorm(ez_obs)
					p_hat <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)  # numerical safety
					ll_row <- y_obs * log(p_hat) + (1 - y_obs) * log(1 - p_hat)
				} else if (family == "tobit") {
					# censored gaussian at 0: density above 0, mass at 0 below
					sigma_t <- sqrt(s2)
					is_zero <- y_obs <= 0
					ll_row[is_zero]  <- pnorm(0, mean = ez_obs[is_zero], sd = sigma_t, log.p = TRUE)
					ll_row[!is_zero] <- dnorm(y_obs[!is_zero],
					                          mean = ez_obs[!is_zero],
					                          sd   = sigma_t, log = TRUE)
				} else if (family == "poisson") {
					# overdispersed poisson: rate = exposure_t * exp(ez).
					# when poisson_exposure_active, build the exposure
					# vector that aligns with obs_idx[, 3] (the period
					# index of each observed dyad).
					if (poisson_exposure_active) {
						t_per_obs <- obs_idx[, 3L]
						exposure_per_obs <- period_exposure[t_per_obs]
						lam <- pmin(exposure_per_obs * exp(ez_obs), 1e8)
					} else {
						lam <- pmin(exp(ez_obs), 1e8)
					}
					ll_row <- dpois(round(y_obs), lambda = lam, log = TRUE)
				} else if (family == "ordinal") {
					# closed-form cumulative-probit log-lik using empirical
					# cutpoints inferred from observed y.
					ll_row <- .ordinal_pointwise_loglik(y_obs, ez_obs)
				} else {
						# frn/rrl ranked likelihoods use the augmented-data normal
						# contribution unless observed_ghk is requested
					if (iter == 1L) {
						cli::cli_warn(c(
							"{.arg save_log_lik} = TRUE: family {.val {family}} pointwise log-lik uses the augmented-data normal approximation.",
							"i" = "Use with caution for {.fn loo::loo} / {.fn waic}; the exact marginal log-lik for {.val {family}} is not implemented."))
					}
					ll_row <- dnorm(Z[obs_idx], mean = ez_obs, sd = sqrt(s2), log = TRUE)
				}
				if (!is.null(LOG_LIK)) {
					LOG_LIK[iter, ] <- ll_row
				} else {
					# chunked: append this iteration's chunk-sized vectors to
					# each connection in column-range order. each chunk file
					# ends up holding ll values in [iter, col] row-major order
					# within its column range.
					for (cc in seq_len(LOG_LIK_META$n_chunks)) {
						writeBin(
							ll_row[LOG_LIK_META$chunk_starts[cc]:LOG_LIK_META$chunk_ends[cc]],
							LOG_LIK_CONNS[[cc]], size = 8L)
					}
				}
			}

			# update posterior sums of random effects
			if(.dyn_uv_has_rank && !is.null(U_cube) && !is.null(V_cube)) {
				U_SUM <- U_SUM + U_cube
				V_SUM <- V_SUM + V_cube
				if(exists("U_SUM_SQ")) { U_SUM_SQ <- U_SUM_SQ + U_cube^2 }
				if(identical(dynamic_uv_kind, "snap") && !is.null(DELTA_SUM) && !is.null(delta_uv)) {
					du_store <- delta_uv
					du_store[!snap_eligible_u] <- NA_real_
					DELTA_SUM <- DELTA_SUM + ifelse(is.finite(du_store), du_store, 0)
					DELTA_COUNT <- DELTA_COUNT + ifelse(is.finite(du_store), 1L, 0L)
					if (!is.null(DELTA_SUM_V) && !is.null(delta_uv_v)) {
						dv_store <- delta_uv_v
						dv_store[!snap_eligible_v] <- NA_real_
						DELTA_SUM_V <- DELTA_SUM_V + ifelse(is.finite(dv_store), dv_store, 0)
						DELTA_COUNT_V <- DELTA_COUNT_V + ifelse(is.finite(dv_store), 1L, 0L)
					}
					if (!is.null(SNAP_DRAWS)) {
						SNAP_DRAWS[iter, , ] <- du_store
						if (!symmetric && !is.null(SNAP_DRAWS_V) && !is.null(delta_uv_v)) {
							SNAP_DRAWS_V[iter, , ] <- dv_store
						}
					}
				}
				# store uv products for each time point
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
					# bipartite: use 2d slices with g interaction matrix
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
				# reuse ez from earlier in the iteration
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
					if(family=="poisson"){
						# route through the exposure-aware simulate when
						# exposures are non-trivial; otherwise the unscaled
						# simy_pois path.
						if (poisson_exposure_active) {
							Ys[,,t] <- simY_pois_exposure(EZ[,,t], period_exposure[t])
						} else {
							Ys[,,t] <- simY_pois(EZ[,,t])
						}
					}
				
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
						# fixed rank nomination. guard na / inf in zs so the
						# rank() does not propagate na into ys, which downstream
						# code (match() on sort(unique())) cannot accept.
						ZS <- EZ_t + matrix(rnorm(nA*nB), nA, nB)
						ZS[!is.finite(ZS)] <- 0
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
						# row rank likelihood. same na guard as frn above.
						ZS <- EZ_t + matrix(rnorm(nA*nB), nA, nB)
						ZS[!is.finite(ZS)] <- 0
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
			
			# save posterior predictive gof stats
			if(gof){
				Ys[is.na(Y)] <- NA
				if(bip) {
					GOF[,,(iter)+1] <- apply(Ys, 3, gof_stats_bipartite, warn_square = FALSE)
				} else {
					GOF[,,(iter)+1] <- apply(Ys, 3, gof_stats)
				}
			}
			
			# print mc progress 
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
				if(.dyn_uv_has_rank && dynamic_ab) {
					start_vals <- list(Z=Z,beta=beta,a=a_mat,b=b_mat,U=U_cube,V=V_cube,
													 rho=rho,s2=s2,Sab=Sab,rho_uv=rho_uv,sigma_uv=sigma_uv,
													 rho_ab=rho_ab, sigma_ab=sigma_ab)
				} else if(.dyn_uv_has_rank) {
					start_vals <- list(Z=Z,beta=beta,a=a,b=b,U=U_cube,V=V_cube,
													 rho=rho,s2=s2,Sab=Sab,rho_uv=rho_uv,sigma_uv=sigma_uv)
				} else if(dynamic_ab) {
					start_vals <- list(Z=Z,beta=beta,a=a_mat,b=b_mat,U=U,V=V,rho=rho,s2=s2,Sab=Sab,
													 rho_ab=rho_ab, sigma_ab=sigma_ab)
				} else {
					start_vals <- list(Z=Z,beta=beta,a=a,b=b,U=U,V=V,rho=rho,s2=s2,Sab=Sab)
				}
				fit <- get_fit_object( APS=APS, BPS=BPS, UVPS=UVPS, YPS=YPS,
														 BETA=BETA, VC=VC, GOF=GOF, Xlist=Xlist, actorByYr=actorByYr, colActorByYr=if(bip) colActorByYr else NULL,
														 start_vals=start_vals, symmetric=symmetric, tryErrorChecks=tryErrorChecks,
														 model.name=model.name, family=family, odmax=odmax,
														 nA=if(bip) nA else NULL, nB=if(bip) nB else NULL, n_time=N,
														 dynamic_uv=dynamic_uv, dynamic_ab=dynamic_ab, bip=bip,
														 G=if(bip) G else NULL,
														 dynamic_beta = beta_dyn$any,
														 beta_dynamic_mask = beta_dyn$mask,
														 beta_dynamic_groups = beta_dyn$groups,
														 rho_beta   = rho_beta_by_group,
														 sigma_beta = sigma_beta_by_group,
														 RHO_BETA   = RHO_BETA,
														 SIGMA_BETA = SIGMA_BETA)
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
				current_U <- if(.dyn_uv_has_rank) U_SUM / iter else U
				current_V <- if(.dyn_uv_has_rank) V_SUM / iter else V
				
				plot_result <- try({
					temp_fit <- get_fit_object( APS=current_APS, BPS=current_BPS, UVPS=current_UVPS, YPS=YPS,
																		 BETA=BETA[1:iter,,drop=FALSE], VC=VC[1:iter,,drop=FALSE],
																		 GOF=GOF, Xlist=Xlist, actorByYr=actorByYr, colActorByYr=if(bip) colActorByYr else NULL,
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
			# store per-actor draws / streaming summaries
			if (per_actor_active) {
				if (keep_per_actor == "draws" || keep_per_actor == "auto") {
					if (exists("THETA_ACTOR", inherits = FALSE)) {
						THETA_ACTOR[iter, , ] <- theta_actor
					}
				}
				if (keep_per_actor == "summary") {
					theta_actor_n_obs <- theta_actor_n_obs + 1L
					# welford streaming mean / second moment
					delta_run <- theta_actor - theta_actor_mean
					theta_actor_mean <- theta_actor_mean + delta_run / theta_actor_n_obs
					theta_actor_m2   <- theta_actor_m2 + delta_run *
						(theta_actor - theta_actor_mean)
				}
				RHO_ACTOR[iter]   <- rho_actor
				SIGMA_ACTOR[iter] <- sqrt(sigma_actor_sq)
			}
			iter<-iter+1
		} # post burn-in
		if(verbose && s > burn){
			# only update sampling progress if verbose=true
			cli::cli_progress_update()
		}
	} # end mcmc  
	if(verbose){
		# only show completion message if verbose=true
		cli::cli_progress_done()
		cli::cli_alert_success("MCMC sampling complete")
	}

	total_iters <- nscan + burn
	# tightened from 10% to 5% (at 10% the posterior is already unreliable),
	# and split into a hard "warn loudly" tier (>5%) and a soft "fyi" tier
	# (any failures at all). the user can also inspect fit$tryerrorchecks
	# programmatically for the raw counts.
	loud_threshold <- max(1, round(0.05 * total_iters))
	loud_params <- character(0)
	soft_params <- character(0)
	for(pname in names(tryErrorChecks)) {
		cnt <- tryErrorChecks[[pname]]
		if(cnt > loud_threshold) {
			loud_params <- c(loud_params,
				paste0(pname, " (", cnt, "/", total_iters, " iterations, ",
				       round(100 * cnt / total_iters, 1), "%)"))
		} else if (cnt > 0) {
			soft_params <- c(soft_params,
				paste0(pname, " (", cnt, "/", total_iters, ")"))
		}
	}
	if(length(loud_params) > 0) {
		cli::cli_warn(c(
			"MCMC sampling failures exceeded 5% of iterations:",
			"x" = "Affected parameters: {.val {loud_params}}",
			"i" = "These parameters may not have converged. Consider increasing {.arg nscan}, simplifying the model, or inspecting {.code fit$tryErrorChecks}."))
	}
	if(length(soft_params) > 0) {
		cli::cli_inform(c(
			"i" = "Some MCMC steps had occasional sampling failures (< 5%): {.val {soft_params}}.",
			"i" = "These were caught by the per-step fallback and are unlikely to bias the posterior. Inspect {.code fit$tryErrorChecks} for the raw counts."))
	}

	# warn when dynamic beta paths hit clipping often enough to affect the fit
	if (exists("beta_path_clipped_hits", inherits = FALSE)) {
		# each sweep can register several clip events (one per dynamic /
		# static coefficient), so normalise by the per-sweep opportunity
		# count and cap at 1 to report a stable fraction-of-sweeps figure.
		n_dyn_sweeps <- max(1L, nscan + burn)
		n_clip_slots <- max(1L, length(beta_dyn$dynamic_idx %||% integer(0)) +
		                          length(beta_dyn$static_idx %||% integer(0)))
		clip_rate <- min(1, (beta_path_clipped_hits + beta_sigma_capped_hits) /
		                     (n_dyn_sweeps * n_clip_slots))
		if (clip_rate > 0.05) {
			cli::cli_warn(c(
				"!" = "{.arg dynamic_beta} chain hit the divergence safeguard on roughly {.val {round(100 * clip_rate)}}% of sweeps; the dynamic coefficients were clipped to stay finite.",
				"i" = "Usual causes: near-separable or very sparse binary data (too few ties per period to identify a per-period coefficient), or an unstandardised large-magnitude outcome.",
				"i" = "Refit after standardising a continuous outcome, simplifying {.arg dynamic_beta}, using a static coefficient, or switching to {.code dynamic_beta_kind = \"rw1\"}."))
		}
	}

	# save start_vals for future model runs
	if(dynamic_ab) {
		start_vals <- list( Z=Z, beta=beta, a=a_mat, b=b_mat, U=U, V=V, rho=rho, s2=s2, Sab=Sab,
											rho_ab=rho_ab, sigma_ab=sigma_ab)
	} else {
		start_vals <- list( Z=Z, beta=beta, a=a, b=b, U=U, V=V, rho=rho, s2=s2, Sab=Sab)
	}
	# store g in start_vals for bipartite
	if(bip && !is.null(G)) start_vals$G <- G
	
	if(!is.null(model.name)) {
		# count effective parameters
		p_eff <- ncol(BETA) + rvar*n + cvar*n + R*(R+1)/2
		
	}
	
	# output
	if(.dyn_uv_has_rank) {
		# return 3d arrays for dynamic uv
		n_draws_uv <- nscan/odens
		U_final <- U_SUM / n_draws_uv
		V_final <- V_SUM / n_draws_uv
		UVPS_final <- UVPS / n_draws_uv
		# posterior sd of latent positions from the running sum of squares
		if(exists("U_SUM_SQ")) {
			U_var <- U_SUM_SQ / n_draws_uv - U_final^2
			U_var[U_var < 0] <- 0
			U_SD <- sqrt(U_var)
		} else { U_SD <- NULL }
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
		# leave raw sums -- get_fit_object divides by nrow(vc)
		APS_final <- APS
		BPS_final <- BPS
	}
	
	fit <- get_fit_object( APS=APS_final, BPS=BPS_final, UVPS=UVPS_final, YPS=YPS,
											 BETA=BETA, VC=VC, GOF=GOF, Xlist=Xlist, actorByYr=actorByYr, colActorByYr=if(bip) colActorByYr else NULL,
											 start_vals=start_vals, symmetric=symmetric, tryErrorChecks=tryErrorChecks,
											 model.name=model.name, U=U_final, V=V_final,
											 dynamic_uv=dynamic_uv, dynamic_ab=dynamic_ab, bip=bip,
											 rho_ab=RHO_AB, rho_uv=RHO_UV,
											 family=family, odmax=odmax, nA=if(bip) nA else NULL,
											 nB=if(bip) nB else NULL, n_time=N,
											 Y_obs=Y, G=if(bip) G else NULL,
											 dynamic_beta = beta_dyn$any,
											 beta_dynamic_mask = beta_dyn$mask,
											 beta_dynamic_groups = beta_dyn$groups,
											 rho_beta   = rho_beta_by_group,
											 sigma_beta = sigma_beta_by_group,
											 RHO_BETA   = RHO_BETA,
											 SIGMA_BETA = SIGMA_BETA)
	# attach the log-likelihood matrix so loo.lame() can find it
	if (!is.null(LOG_LIK)) {
		fit$log_lik <- LOG_LIK
		fit$log_lik_index <- obs_idx
	}
	# chunked mode: close all connections and attach metadata
	if (!is.null(LOG_LIK_CONNS)) {
		for (cc in LOG_LIK_CONNS) close(cc)
		fit$log_lik_index <- obs_idx
		fit$log_lik_chunks <- LOG_LIK_META
	}
		# record dynamic argument values on the fit for downstream tools
	fit$dynamic_beta_pool <- dynamic_beta_pool
	fit$log_lik_method <- log_lik_method
	fit$dynamic <- list(
		uv = isTRUE(dynamic_uv),
		uv_kind = if (isTRUE(dynamic_uv)) dynamic_uv_kind else NULL,
		ab = isTRUE(dynamic_ab),
		beta = isTRUE(beta_dyn$any),
		beta_kind = if (isTRUE(beta_dyn$any)) dynamic_beta_kind else NULL,
		G = isTRUE(dynamic_G),
		snap_model_estimated = isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "snap") &&
			isTRUE(.dyn_uv_has_rank),
		t_model_estimated = isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t") &&
			isTRUE(.dyn_uv_has_rank),
		mode = mode)
	# innovation-sd chains for the dynamic additive / multiplicative blocks.
	# get_fit_object() does not take these, so attach them here so the fit
	# carries fit$sigma_ab / fit$sigma_uv (used by print.lame()'s
	# sigma_innov column and documented in the dynamic-effects vignette).
	if (dynamic_ab && !is.null(SIGMA_AB)) fit$sigma_ab <- SIGMA_AB
	if (dynamic_uv && !is.null(SIGMA_UV)) fit$sigma_uv <- SIGMA_UV
	# posterior sd of the dynamic latent positions (same dims as fit$u)
	if (.dyn_uv_has_rank && exists("U_SD") && !is.null(U_SD)) {
		dimnames(U_SD) <- dimnames(fit$U)
		fit$U_sd <- U_SD
	}
	# attach snap-shift outputs: posterior snap probabilities, optional
	# retained indicator draws, and snap-rate chain.
		if (isTRUE(fit$dynamic$snap_model_estimated)) {
			fit$pi_snap <- PI_SNAP
			if (!is.null(DELTA_SUM)) {
				.snap_prob_from_counts <- function(sum_mat, count_mat) {
					out <- sum_mat / count_mat
					out[!is.finite(out) | count_mat <= 0L] <- NA_real_
					out
				}
				fit$snap_prob <- .snap_prob_from_counts(DELTA_SUM, DELTA_COUNT)
				if (!is.null(fit$U) && length(dim(fit$U)) == 3) {
					dimnames(fit$snap_prob) <- list(dimnames(fit$U)[[1]], dimnames(fit$U)[[3]])
					names(dimnames(fit$snap_prob)) <- c("actor", "time")
					dimnames(DELTA_COUNT) <- dimnames(fit$snap_prob)
				}
				if (!symmetric) {
					fit$snap_prob_v <- .snap_prob_from_counts(DELTA_SUM_V, DELTA_COUNT_V)
					if (!is.null(fit$V) && length(dim(fit$V)) == 3) {
						dimnames(fit$snap_prob_v) <- list(dimnames(fit$V)[[1]], dimnames(fit$V)[[3]])
						names(dimnames(fit$snap_prob_v)) <- c("actor", "time")
						dimnames(DELTA_COUNT_V) <- dimnames(fit$snap_prob_v)
					}
				}
				.snap_summary_from_prob <- function(prob, count_mat) {
					n_mat <- count_mat
					dimnames(n_mat) <- dimnames(prob)
					n_mat[!is.finite(prob)] <- 0
					var_mat <- matrix(NA_real_, nrow(prob), ncol(prob), dimnames = dimnames(prob))
					ok <- is.finite(prob) & n_mat > 1L
					var_mat[ok] <- prob[ok] * (1 - prob[ok]) * n_mat[ok] / (n_mat[ok] - 1L)
					list(n = n_mat, mean = prob, variance = var_mat)
				}
				if (SNAP_DRAWS_MODE %in% c("summary", "draws", "chunked")) {
					fit$snap_draws_summary <- .snap_summary_from_prob(fit$snap_prob, DELTA_COUNT)
					if (!symmetric && !is.null(fit$snap_prob_v)) {
						fit$snap_draws_summary_v <- .snap_summary_from_prob(fit$snap_prob_v, DELTA_COUNT_V)
					}
				}
			if (!is.null(SNAP_DRAWS)) {
				dimnames(SNAP_DRAWS) <- list(
					draw = paste0("draw", seq_len(dim(SNAP_DRAWS)[1L])),
					actor = rownames(fit$snap_prob),
					time = colnames(fit$snap_prob))
				fit$snap_draws <- SNAP_DRAWS
				if (!symmetric && !is.null(SNAP_DRAWS_V) && !is.null(fit$snap_prob_v)) {
					dimnames(SNAP_DRAWS_V) <- list(
						draw = paste0("draw", seq_len(dim(SNAP_DRAWS_V)[1L])),
						actor = rownames(fit$snap_prob_v),
						time = colnames(fit$snap_prob_v))
					fit$snap_draws_v <- SNAP_DRAWS_V
				}
				fit$snap_draws_meta <- list(
					mode = SNAP_DRAWS_MODE,
						storage = "in_memory_array",
						dim_order = c("draw", "actor", "time"),
						n_draws = dim(SNAP_DRAWS)[1L],
						first_period = "ineligible_na",
						transition_eligibility = "actor must be present in both adjacent periods")
				}
			}
		}
	# attach the per-fit randomised-halton metadata so re-runs are
	# reproducible and downstream tools can record what mc sequence was
	# used. seed is the fit's main seed; shift is the [0, 1) digital
	# shift; n_mc is the ghk chain length (currently hard-coded 64).
	if (identical(log_lik_method, "observed_ghk")) {
		fit$log_lik_meta <- list(
			seed = seed,
			halton_shift = .lame_log_lik_halton_shift,
			n_mc = 64L,
			method = "GHK + Halton (randomised) + Gauss-Hermite (cbin)",
			dim_cap = 15L)
	}
	# attach explicit ordinal cutpoint posterior. the posterior is on the
	# (k-1)-vector alpha with alpha_1 = 0 (identification anchor, not stored).
	# `ordinal_levels` records the unique observed levels (positionally
	# matched to category indices 1..k) so downstream code can map alpha
	# back to user-facing labels.
	fit$ordinal_cutpoints <- ordinal_cutpoints
	if (use_explicit_cutpoints) {
		fit$ordinal_levels <- ord_lvls
		if (!is.null(ALPHA)) {
			fit$ALPHA <- ALPHA
			fit$alpha_post_mean <- c(0, colMeans(ALPHA, na.rm = TRUE))
		} else {
			# for k = 2, the posterior is point mass at alpha_1 = 0
			fit$alpha_post_mean <- 0
		}
	}
	# attach the per-period g cube when dynamic_g is active, surfacing the
	# posterior-mean cube
	# (g_cube_post_mean), an sd cube, the ar(1) hyperparameter draws
	# (rho_g, sigma_g2), and a rotation-drift diagnostic that flags when
	# apparent g_t variation is dominated by u/v rotation rather than real
	# temporal change.
	if (isTRUE(dynamic_G) && exists("G_cube", inherits = FALSE) &&
	    !is.null(G_cube)) {
		fit$G_cube <- G_cube
		fit$dynamic_G <- TRUE
		fit$RHO_G <- RHO_G
		fit$SIGMA_G2 <- SIGMA_G2
		if (!is.null(G_cube_sum) && G_cube_n > 0L) {
			fit$G_cube_post_mean <- G_cube_sum / G_cube_n
			# unbiased posterior sd per cell
			var_est <- (G_cube_sumsq / G_cube_n) - (fit$G_cube_post_mean)^2
			var_est[var_est < 0] <- 0
			fit$G_cube_post_sd <- sqrt(var_est)
		}
		# rotation-drift diagnostic using the working u/v at the end of
		# mcmc (a reasonable summary; per-draw canonicalisation would be
		# more expensive). uses fit$g_cube as the raw input.
		U_for_diag <- tryCatch({
			if (!is.null(fit$U) && length(dim(fit$U)) == 3L) fit$U
			else if (!is.null(fit$U)) {
				Uc <- array(0, dim = c(nrow(fit$U), ncol(fit$U), N))
				for (t in seq_len(N)) Uc[, , t] <- fit$U
				Uc
			} else NULL
		}, error = function(e) NULL)
		V_for_diag <- tryCatch({
			if (!is.null(fit$V) && length(dim(fit$V)) == 3L) fit$V
			else if (!is.null(fit$V)) {
				Vc <- array(0, dim = c(nrow(fit$V), ncol(fit$V), N))
				for (t in seq_len(N)) Vc[, , t] <- fit$V
				Vc
			} else NULL
		}, error = function(e) NULL)
		if (!is.null(U_for_diag) && !is.null(V_for_diag)) {
			fit$G_rotation_drift <- tryCatch(
				.rotation_drift_diagnostic(fit$G_cube, U_for_diag, V_for_diag),
				error = function(e) NULL)
			if (!is.null(fit$G_rotation_drift) &&
			    isTRUE(fit$G_rotation_drift$flag)) {
				cli::cli_warn(c(
					"{.code dynamic_G} rotation-drift diagnostic flagged the fit (ratio = {.val {round(fit$G_rotation_drift$ratio, 2)}}).",
					"i" = "Apparent G_t variation is dominated by U/V rotation drift, not real temporal change.",
					"i" = "Use {.code fit$G_cube_post_mean} (canonical reporting) rather than per-draw {.code fit$G_cube}."))
			}
		}
	}
	# attach per-actor slopes (draws or streaming summary)
	if (isTRUE(per_actor_active)) {
		fit$dynamic_beta_per_actor <- dynamic_beta_per_actor
		fit$per_actor_covariate_idx <- per_actor_covariate_idx
		fit$per_actor_identifiability <- per_actor_identifiability
		fit$RHO_ACTOR <- RHO_ACTOR
		fit$SIGMA_ACTOR <- SIGMA_ACTOR
		if (keep_per_actor %in% c("draws", "auto") &&
		    exists("THETA_ACTOR", inherits = FALSE)) {
			fit$THETA_ACTOR <- THETA_ACTOR
			# posterior-mean summary
			fit$theta_actor_mean <- apply(THETA_ACTOR, c(2, 3), mean)
		} else if (keep_per_actor == "summary") {
			fit$theta_actor_mean <- theta_actor_mean
			fit$theta_actor_sd   <- sqrt(theta_actor_m2 /
				max(1L, theta_actor_n_obs - 1L))
		}
	}
	# reporting-only uv recentering  the
	# canonical $u / $v / $beta stay in the *sampled* parameterisation so
	# downstream eta-reconstruction is unambiguous; centred copies live in
	# $u_centered / $v_centered / $beta_centered for plotting.
	if (!is.null(fit$U) && !is.null(fit$V) &&
	    NCOL(fit$U) > 0L && NCOL(fit$V) > 0L) {
		.uv_centered <- tryCatch(.compute_uv_centered(fit),
		                          error = function(e) NULL)
		if (!is.null(.uv_centered)) {
			fit$U_centered    <- .uv_centered$U_centered
			fit$V_centered    <- .uv_centered$V_centered
			fit$BETA_centered <- .uv_centered$BETA_centered
		}
	}
	if (!is.null(time_index))       fit$time_index <- time_index
	if (!is.null(period_exposure))  fit$period_exposure <- period_exposure
	if (is.finite(max_seconds))     fit$max_seconds <- max_seconds
	fit$terminated_early <- isTRUE(.lame_terminated_early)
	if (!is.null(checkpoint_path))  fit$checkpoint_path <- checkpoint_path

	####
	# scalar typing hygiene
	if (!is.null(fit$RHO)) fit$RHO <- as.numeric(fit$RHO)
	if (!is.null(fit$s2))  fit$s2  <- as.numeric(fit$s2)
	####

	####
	# normalize ez to list-of-matrices
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
		if (is.finite(rho_hat)) fit$rho_uv_aligned <- as.numeric(rho_hat)
	}
	####

	if(symmetric) {
		options(warn = old_warn)
	}

	# store bipartite mh proposal diagnostics
	if(bip && RA > 0 && RB > 0) {
		fit$mh_diagnostics <- list(
			sd_U = mh_sd_U, sd_V = mh_sd_V
		)
	}

	fit$call <- mc
	# record which dynamic_beta_kind was used so summary/print and
	# any downstream forecasting (predict h>0) can propagate the right state-
	# space model. only meaningful when beta_dyn$any is true.
	# also record the *requested* kind so downstream code can tell when
	# the effective kind fell back to "ar1".
	if (isTRUE(beta_dyn$any)) {
		fit$dynamic_beta_kind <- dynamic_beta_kind
		if (!identical(dynamic_beta_kind, dynamic_beta_kind_requested)) {
			fit$dynamic_beta_kind_requested <- dynamic_beta_kind_requested
		}
	}
	if (isTRUE(dynamic_uv)) {
		fit$dynamic_uv_kind <- dynamic_uv_kind
		if (!identical(dynamic_uv_kind, dynamic_uv_kind_requested)) {
			fit$dynamic_uv_kind_requested <- dynamic_uv_kind_requested
		}
	}
	# attach the pre-conversion data snapshot if freeze_call was set.
	# the snapshot was taken at the top of lame() before list_to_array etc.,
	# so update(fit, ...) can re-fit against the original user-facing shape.
	if (isTRUE(freeze_call)) {
		fit$data_snapshot <- .lame_freeze_snapshot
		fit$freeze_call   <- TRUE
	}
	# expose the resolved g and the prior list (with sampler-side defaults
	# filled in) so prior_summary() can show what was actually used
	fit$g     <- g
	fit$prior <- prior
	# expose the raw sampling-failure counts so users can audit failure
	# masking themselves; the post-mcmc warning above only fires above a
	# threshold. also alias the same counter list as `fit$mh_counters`.
	fit$tryErrorChecks <- tryErrorChecks
	fit$mh_counters    <- tryErrorChecks
	# expose the structural flags / latent dims so format_ame_call() can
	# render the model formula on summary.lame
	fit$rvar   <- isTRUE(rvar)
	fit$cvar   <- isTRUE(cvar)
	fit$dcor   <- isTRUE(dcor)
	fit$R      <- R
	if (bip) {
		fit$R_row <- if (!is.null(R_row)) R_row else R
		fit$R_col <- if (!is.null(R_col)) R_col else R
	}
	class(fit) <- c("lame", "ame")
	return(fit)
	
}
