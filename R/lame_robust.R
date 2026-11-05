# validated inputs for the continuous, directed robust model
.lame_robust_input <- function(Y, Xdyad, R) {
	if (!is.list(Y) || length(Y) < 2L || !all(vapply(Y, is.matrix, logical(1))))
		cli::cli_abort("Y must be a list of at least two square numeric matrices.", call = NULL)
	n <- nrow(Y[[1L]])
	if (n < 2L || !all(vapply(
		Y,
		function(y) is.numeric(y) && !is.complex(y) &&
			identical(dim(y), c(n, n)) && !any(is.infinite(y)),
		logical(1)
	)))
		cli::cli_abort("Y matrices must be numeric, square, equally sized, and free of infinite values.", call = NULL)
	if (length(R) != 1L || !is.numeric(R) || is.complex(R) || !is.finite(R) || R != floor(R) || R < 0L || R >= n)
		cli::cli_abort("R must be an integer between zero and the number of actors minus one.", call = NULL)
	actors <- rownames(Y[[1L]])
	if (is.null(actors)) actors <- colnames(Y[[1L]])
	if (is.null(actors)) actors <- as.character(seq_len(n))
	if (anyNA(actors) || any(!nzchar(actors)) || anyDuplicated(actors))
		cli::cli_abort("Actor names must be nonmissing and unique.", call = NULL)
	align <- function(y) {
		for (margin in 1:2) {
			nm <- dimnames(y)[[margin]]
			if (!is.null(nm) && (anyDuplicated(nm) || !setequal(nm, actors)))
				cli::cli_abort("Every matrix must contain the same actors on both margins.", call = NULL)
		}
		ri <- if (is.null(rownames(y))) seq_len(n) else match(actors, rownames(y))
		ci <- if (is.null(colnames(y))) seq_len(n) else match(actors, colnames(y))
		y <- y[ri, ci, drop = FALSE]
		dimnames(y) <- list(actors, actors)
		diag(y) <- NA_real_
		y
	}
	Y <- lapply(Y, align)
	if (!any(vapply(Y, function(y) any(is.finite(y)), logical(1))))
		cli::cli_abort("Y must contain at least one observed off-diagonal value.", call = NULL)
	years <- names(Y)
	if (is.null(years)) years <- as.character(seq_along(Y))
	if (anyNA(years) || any(!nzchar(years)) || anyDuplicated(years))
		cli::cli_abort("Period names must be nonmissing and unique.", call = NULL)
	names(Y) <- years
	if (!is.null(Xdyad)) {
		if (!is.list(Xdyad) || length(Xdyad) != length(Y))
			cli::cli_abort("Xdyad must be NULL or a list with one array per period.", call = NULL)
		if (!is.null(names(Xdyad))) {
			if (anyDuplicated(names(Xdyad)) || !setequal(names(Xdyad), years))
				cli::cli_abort("Xdyad period names must match Y.", call = NULL)
			Xdyad <- Xdyad[years]
		}
		p <- if (is.matrix(Xdyad[[1L]])) 1L else dim(Xdyad[[1L]])[3L]
		if (length(p) != 1L || is.na(p) || p < 1L)
			cli::cli_abort("Each Xdyad element must have dimensions actor by actor by covariate.", call = NULL)
		cn <- if (length(dim(Xdyad[[1L]])) == 3L) dimnames(Xdyad[[1L]])[[3L]] else NULL
		if (is.null(cn)) cn <- paste0("x", seq_len(p))
		if (anyNA(cn) || any(!nzchar(cn)) || anyDuplicated(cn))
			cli::cli_abort("Covariate names must be nonmissing and unique.", call = NULL)
		Xdyad <- lapply(Xdyad, function(x) {
			if (is.matrix(x)) x <- array(x, c(dim(x), 1L), dimnames = c(dimnames(x), list(NULL)))
			if (!is.numeric(x) || is.complex(x) || !identical(dim(x), c(n, n, p)))
				cli::cli_abort("All Xdyad arrays must have the same actor and covariate dimensions.", call = NULL)
			dn <- dimnames(x)
			if (is.null(dn)) dn <- vector("list", 3L)
			wanted <- list(actors, actors, cn)
			index <- lapply(1:3, function(k) {
				if (is.null(dn[[k]])) return(seq_along(wanted[[k]]))
				if (anyDuplicated(dn[[k]]) || !setequal(dn[[k]], wanted[[k]]))
					cli::cli_abort("Xdyad actor and covariate names must match across periods.", call = NULL)
				match(wanted[[k]], dn[[k]])
			})
			x <- x[index[[1L]], index[[2L]], index[[3L]], drop = FALSE]
			for (k in seq_len(p)) {
				if (any(!is.finite(x[, , k][row(x[, , k]) != col(x[, , k])])))
					cli::cli_abort("Xdyad must be finite off the diagonal, including at missing outcomes.", call = NULL)
				diag(x[, , k]) <- 0
			}
			dimnames(x) <- wanted
			x
		})
		names(Xdyad) <- years
	}
	list(Y = Y, Xdyad = Xdyad, R = as.integer(R))
}

.lame_robust_seed <- function(seed, fun) {
	if (is.null(seed)) return(fun())
	if (length(seed) != 1L || !is.numeric(seed) || is.complex(seed) || !is.finite(seed) ||
			seed != floor(seed) || seed < 0 || seed > .Machine$integer.max)
		cli::cli_abort("seed must be NULL or a nonnegative integer.", call = NULL)
	had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	if (had) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	on.exit({
		if (had) assign(".Random.seed", old, envir = .GlobalEnv)
		else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
			rm(".Random.seed", envir = .GlobalEnv)
	}, add = TRUE)
	set.seed(as.integer(seed))
	fun()
}

.lame_robust_prior_check <- function(prior, p) {
	allowed <- c(
		"ve_shape", "ve_rate", "ar_shape", "mu_rate", "ab_rate",
		"uv_rate", "level_df", "level_scale", "beta_precision"
	)
	if (!is.list(prior) || (length(prior) && (is.null(names(prior)) ||
			any(!names(prior) %in% allowed) || anyDuplicated(names(prior)) ||
			any(vapply(prior, is.null, logical(1))))))
		cli::cli_abort("prior must be a named list of documented robust-model prior settings.", call = NULL)
	for (nm in setdiff(names(prior), c("level_scale", "beta_precision"))) {
		x <- prior[[nm]]
		lower <- if (nm %in% c("ve_shape", "ar_shape")) 1 else if (nm == "level_df") 3 else 0
		if (length(x) != 1L || !is.numeric(x) || is.complex(x) || !is.finite(x) || x <= lower)
			cli::cli_abort("Invalid robust-model prior setting: {.field {nm}}.", call = NULL)
	}
	if (!is.null(prior$level_scale)) {
		x <- prior$level_scale
		if (!is.matrix(x) || !identical(dim(x), c(2L, 2L)) || !is.numeric(x) || is.complex(x) ||
				any(!is.finite(x)) || !isTRUE(all.equal(x, t(x))) ||
				inherits(try(chol(x), silent = TRUE), "try-error"))
			cli::cli_abort("prior$level_scale must be a symmetric positive-definite 2 by 2 matrix.", call = NULL)
	}
	if (!is.null(prior$beta_precision)) {
		x <- prior$beta_precision
		if (!is.numeric(x) || is.complex(x) || !length(x) %in% unique(c(1L, p)) || any(!is.finite(x) | x <= 0))
			cli::cli_abort(
				"prior$beta_precision must be positive, with length one or the number of covariates.", call = NULL
			)
		prior$beta_precision <- rep(x, length.out = p)
	}
	prior
}

#' Robust longitudinal models for directed continuous networks
#'
#' Fit a longitudinal additive and multiplicative effects model with a residual
#' precision shared across all periods of each directed pair. This permits
#' persistent differences in volatility between pairs. All mean components and
#' variance parameters are estimated jointly by MCMC.
#'
#' @param Y List of at least two equally sized square numeric matrices. Rows
#'   are senders, columns are receivers, and list order is temporal order with
#'   equally spaced periods. Missing outcomes are allowed; diagonals are ignored.
#'   The actor set must be constant. When supplied, unique actor and period names
#'   are used to align inputs.
#' @param Xdyad Optional list of sender-by-receiver-by-covariate numeric arrays,
#'   one per period. A matrix denotes a single covariate. Covariates must be
#'   finite off the diagonal, including for missing outcomes. No intercept is
#'   inserted into this design; the model contains a separate time intercept.
#' @param R Nonnegative integer factor rank, smaller than the number of actors.
#' @param df Degrees of freedom greater than two for the persistent pair
#'   precision prior. The default is four. Use `Inf` to remove persistent
#'   pair heterogeneity. With `observation_df = Inf`, this gives a homogeneous
#'   Gaussian residual model with the same mean specification and other priors.
#' @param actor_df Degrees of freedom greater than two for independent
#'   sender-period and receiver-period innovation precisions. `Inf` (default)
#'   retains Gaussian actor innovations. A finite value permits periods with
#'   unusually large changes in sender or receiver activity.
#' @param observation_df Degrees of freedom greater than two for independent
#'   observation precisions. `Inf` (default) adds no temporary scale; finite
#'   values permit isolated residual shocks within otherwise stable pairs.
#' @param burn Number of initial iterations discarded.
#' @param nscan Number of iterations after burn-in; must be a multiple of `odens`.
#' @param odens Retain every `odens`-th iteration. Every retained iteration
#'   includes the complete joint state needed for prediction and simulation.
#' @param seed Optional nonnegative integer sampler seed. An explicit seed
#'   preserves the caller's RNG state. With `NULL`, a seed is drawn from the
#'   caller's stream and recorded in the result.
#' @param prior Named list overriding the proper priors described below.
#' @param verbose Logical; print progress every 5,000 iterations.
#'
#' @details
#' For off-diagonal observations, the conditional model is
#' \deqn{y_{ijt} = x_{ijt}^{T}\beta + \mu_t + a_{it} + b_{jt} +
#'   u_{it}^{T}v_{jt} + e_{ijt}, \qquad
#'   e_{ijt} \mid v_e,\lambda_{ij},\omega_{ijt} \sim N(0,v_e/(\lambda_{ij}\omega_{ijt})).}
#' For finite `df`, independent directed-pair precisions have
#' \eqn{\lambda_{ij}\sim\mathrm{Gamma}(df/2,df/2)} in shape-rate
#' parameterization and persist through observed and forecast periods.
#' With `observation_df = Inf`, each residual has a Student distribution
#' after integrating out the pair precision; residuals for the same pair
#' share that random scale across time. Opposite directions have independent
#' precisions. For `df = Inf`, all pair precisions equal one.
#'
#' For finite `observation_df`, independent precisions
#' \eqn{\omega_{ijt}\sim\mathrm{Gamma}(observation\_df/2,observation\_df/2)}
#' add a temporary residual scale to each observation. Otherwise these
#' precisions equal one. Conditional on the pair precision, integrating out
#' the temporary precision gives a Student residual. When both degrees of
#' freedom are finite, the fully marginal residual is a mixture over both
#' scales and is not a single Student distribution. Replication and forecasting
#' draw fresh temporary precisions; the persistent pair precision is retained.
#'
#' The time intercept follows a zero-mean stationary AR(1). Sender and receiver
#' effects each equal a persistent actor level plus an AR(1)
#' deviation (stationary Gaussian when `actor_df = Inf`). Levels and deviations sum to zero across actors in every period.
#' Sender and receiver levels have an estimated bivariate covariance; their
#' deviations have independent innovations with a shared variance and AR
#' coefficient. Factor coordinates follow independent zero-mean stationary
#' AR(1) processes, sharing an innovation variance and AR coefficient.
#' With finite `actor_df`, each actor role and period has an independent
#' precision \eqn{\tau_{rt}\sim\mathrm{Gamma}(actor\_df/2,actor\_df/2)}.
#' Conditional on this precision, centered actor innovations have covariance
#' \eqn{(v_{ab}/\tau_{rt})(I-11^T/n)}. The initial deviation uses covariance
#' \eqn{v_{ab}(I-11^T/n)/[\tau_{r1}(1-\rho_{ab}^2)]}. This specifies the
#' initial distribution for the finite panel; it is not asserted to be the
#' stationary distribution of the Student-innovation process. Future periods
#' draw new role precisions, while actor levels persist.
#'
#' There are three AR coefficients: for the intercept, actor deviations,
#' and factors. Each has a normal prior with mean 0.8 and standard deviation
#' 0.3, truncated to (-0.995, 0.995).
#'
#' Let `v` be the variance of observed off-diagonal responses, resolved once
#' before sampling (one if that variance is undefined or nonpositive).
#' The prior overrides are:
#' \describe{
#'   \item{`ve_shape`, `ve_rate`}{Residual variance inverse-gamma shape and
#'     rate, defaulting to 2 and `v`. The density is proportional to
#'     \eqn{x^{-a-1}\exp(-b/x)}. The shape must exceed one.}
#'   \item{`ar_shape`, `mu_rate`, `ab_rate`, `uv_rate`}{Inverse-gamma priors
#'     for intercept, actor-deviation, and factor innovation variances.
#'     The common shape defaults to 2 and must exceed one. Rates default to
#'     `0.1 * v`, `0.1 * v`, and `0.1 * sqrt(v)`, respectively.}
#'   \item{`level_df`, `level_scale`}{Inverse-Wishart prior for the bivariate
#'     actor-level covariance. Defaults are 4 degrees of freedom and `v * I`.
#'     Degrees of freedom must exceed three; scale must be positive definite.
#'     The density is proportional to
#'     \eqn{|S|^{-(\nu+3)/2}\exp\{-\mathrm{tr}(\Psi S^{-1})/2\}}.}
#'   \item{`beta_precision`}{Positive scalar or one value per covariate for
#'     independent zero-mean Gaussian coefficient priors. Default `0.01 / v`
#'     gives variance `100 * v`, in the unit of the response and independently of
#'     the residual variance.}
#' }
#' All rates must be positive. The defaults are written in the unit of the
#' response (a fit of `c * Y` is the fit of `Y` in other units), so only the
#' units of the covariates matter for the priors; inspect `fit$prior` when
#' comparing specifications.
#'
#' This estimator supports directed continuous networks with stable actors.
#' It does not estimate residual reciprocity, changing membership, bipartite
#' structure, binary responses, dynamic covariate coefficients, offsets, or
#' alternative latent laws of motion. It has its own `lame_robust` class and
#' does not use `lame()`'s priors or inherit its fit class.
#'
#' Run at least two independently seeded fits and pass them to
#' [lame::compute_mcmc_diagnostics()]. Diagnostics cover
#' all retained identified parameters, including every off-diagonal latent
#' product. Raw factor coordinates are rotation dependent.
#' [lame::simulate.lame_robust()] compares replicated networks with the
#' observed statistics relevant to the analysis.
#'
#' @return An object of class `lame_robust`, with `draws`, aligned `Y` and
#'   `Xdyad`, resolved `prior`, `seed`, and sampling settings. In `draws`,
#'   `beta`, `mu`, and `hyper` have iterations in rows; `A` and `B` are
#'   actor-by-period-by-iteration arrays; `a0` and `b0` are actor-by-iteration
#'   matrices; `U` and `V` are actor-by-rank-by-period-by-iteration arrays;
#'   `lambda` is sender-by-receiver-by-iteration. With finite `actor_df`,
#'   `tau_ab` is role-by-period-by-iteration (sender, receiver). With finite
#'   `observation_df`, `omega` is sender-by-receiver-by-period-by-iteration.
#'   All components share the
#'   same iteration order. Methods include prediction, simulation, coefficient
#'   summaries, posterior covariance, credible intervals, and observation count.
#' @seealso [lame::predict.lame_robust()], [lame::compute_mcmc_diagnostics()], [lame::lame()]
#' @examples
#' set.seed(21)
#' Y <- lapply(1:3, function(t) {
#'   y <- matrix(rnorm(25), 5, 5)
#'   diag(y) <- NA
#'   y
#' })
#' # Short execution example; use longer independent chains for inference.
#' fit <- lame_robust(Y, R = 0, burn = 10, nscan = 20, odens = 2,
#'                    seed = 1, verbose = FALSE)
#' dim(predict(fit, h = 1))
#' dim(simulate(fit, nsim = 2, seed = 2))
#' @export
lame_robust <- function(
	Y,
	Xdyad = NULL,
	R = 2L,
	df = 4,
	burn = 5000L,
	nscan = 20000L,
	odens = 20L,
	seed = NULL,
	prior = list(),
	verbose = interactive(),
	actor_df = Inf,
	observation_df = Inf
) {
	call <- match.call()
	input <- .lame_robust_input(Y, Xdyad, R)
	for (nm in c("burn", "nscan", "odens")) {
		x <- get(nm)
		if (length(x) != 1L || !is.numeric(x) || is.complex(x) || !is.finite(x) ||
				x != floor(x) || x > .Machine$integer.max || x < if (nm == "burn") 0 else 1)
			cli::cli_abort("{.arg {nm}} must be a valid integer sampling count.", call = NULL)
	}
	if (as.double(burn) + as.double(nscan) > .Machine$integer.max)
		cli::cli_abort("burn plus nscan must not exceed the integer sampling limit.", call = NULL)
	if (nscan %% odens != 0L) cli::cli_abort("nscan must be a multiple of odens.", call = NULL)
	if (length(df) != 1L || !is.numeric(df) || is.complex(df) || is.na(df) || df <= 2)
		cli::cli_abort("df must exceed two; use Inf for homogeneous Gaussian errors.", call = NULL)
	if (length(actor_df) != 1L || !is.numeric(actor_df) || is.complex(actor_df) ||
			is.na(actor_df) || actor_df <= 2)
		cli::cli_abort("actor_df must exceed two; use Inf for Gaussian actor innovations.", call = NULL)
	if (length(observation_df) != 1L || !is.numeric(observation_df) || is.complex(observation_df) ||
			is.na(observation_df) || observation_df <= 2)
		cli::cli_abort("observation_df must exceed two; use Inf for no transient residual scale.", call = NULL)
	if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose))
		cli::cli_abort("verbose must be TRUE or FALSE.", call = NULL)
	D <- .lame_robust_data(input$Y, input$Xdyad, input$R, persistent = is.finite(df))
	prior <- .lame_robust_prior_check(prior, D$p)
	prior$nu <- df
	prior$actor_df <- actor_df
	prior$observation_df <- observation_df
	kernel <- list(
		.lame_robust_weighted_beta = robust_weighted_beta,
		.lame_robust_weighted_mu = robust_weighted_mu,
		.lame_robust_weighted_actor = robust_weighted_actor,
		.lame_robust_weighted_uv = robust_weighted_uv
	)
	fit <- .lame_robust_seed(seed, function() {
		actual_seed <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else as.integer(seed)
		.lame_robust_fit(
			D,
			kernel,
			actual_seed,
			burn = as.integer(burn),
			nscan = as.integer(nscan),
			thin = as.integer(odens),
			prior_override = prior,
			verbose = verbose
		)
	})
	fit$call <- call
	fit$Y <- input$Y
	fit$Xdyad <- input$Xdyad
	fit$R <- input$R
	fit$df <- df
	fit$actor_df <- actor_df
	fit$observation_df <- observation_df
	fit$family <- if (is.finite(observation_df)) {
		if (is.finite(df)) "persistent and transient scale mixtures" else "Student"
	} else if (is.finite(df)) "persistent Student" else "normal"
	fit$method <- "mcmc"
	fit$scope <- "Directed continuous networks, stable actors, AR(1) mean/actor/factor paths with the specified actor innovation law, independent directed residuals conditional on the specified pair and observation precisions."
	colnames(fit$draws$beta) <- colnames(D$X)
	colnames(fit$draws$mu) <- D$years
	class(fit) <- "lame_robust"
	fit
}
