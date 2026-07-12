# penalised als time-varying beta point estimator.
# minimises sum_t ||y_t - x_t beta_t||^2 + lambda sum_t ||beta_t - beta_{t-1}||^2.
# lambda = 0 reduces to per-period ols; large lambda collapses to a
# single pooled beta. returns a point estimate, not a posterior.

#' Penalised ALS time-varying coefficient estimate
#'
#' Computes a fast point estimate of a time-varying coefficient vector
#' \eqn{\beta_t} for a longitudinal network model by solving the
#' first-difference-penalised least-squares problem. No MCMC, no
#' multiplicative effects, no random effects -- purely a regression
#' point estimate with a smoothing penalty on the coefficient path.
#' Useful for rapid exploration and for generating starting values for
#' a full \code{lame(dynamic_beta = ...)} fit.
#'
#' \deqn{\min_{\beta_{1:T}} \sum_t \|y_t - X_t \beta_t\|^2 + \lambda
#'   \sum_{t=2}^{T} \|\beta_t - \beta_{t-1}\|^2}
#'
#' At \eqn{\lambda = 0} this is \eqn{T} independent per-period OLS
#' fits. As \eqn{\lambda \to \infty} the solution converges to a
#' single pooled \eqn{\beta} (constant over time).
#'
#' @param Y A list of \eqn{T} response matrices (or a 3-D array with
#'   third dimension time).
#' @param Xdyad A list of \eqn{T} dyadic covariate arrays (each
#'   \eqn{n \times n \times p}, or a single \eqn{n \times n} matrix if
#'   \eqn{p=1}).
#' @param lambda Non-negative smoothing parameter. Default \code{0}
#'   (no penalty, per-period OLS).
#' @param intercept Logical; include an intercept (default
#'   \code{TRUE}).
#'
#' @return A list with \code{beta} (a \eqn{p \times T} matrix of
#'   estimates), \code{lambda}, \code{residual_ss} (sum of squared
#'   residuals across periods), \code{intercept}, and
#'   \code{call}. Class \code{"als_dynamic_beta"}.
#'
#' @section Relation to other entry points:
#' \code{als_dynamic_beta} is a \emph{regression-only smoother}; the
#' returned object carries time-varying \eqn{\beta_t} alone, no
#' additive (\eqn{a, b}) or multiplicative (\eqn{U, V}) AME
#' components. It is therefore not a special case of
#' \code{\link{lame_als}} (which always estimates the full AME
#' decomposition) and not a special case of
#' \code{\link{lame}(dynamic_beta = TRUE)} (which is Bayesian with
#' an AR(1) / RW1 / RW2 / Matern 3/2 state-space prior). Use it when
#' you want a fast, point-only, deterministic smoother for the
#' coefficient path; use \code{\link{lame}(dynamic_beta = ...)} when
#' you need full posterior uncertainty and additive / multiplicative
#' effects.
#'
#' @seealso \code{\link{lame}} for the full Bayesian dynamic-\eqn{\beta}
#'   fit (AR(1) / RW1 / RW2 / Matern 3/2);
#'   \code{\link{lame_als}} for the longitudinal AME point estimator
#'   (static \eqn{\beta}, full \eqn{a, b, U, V}).
#'
#' @examples
#' set.seed(1)
#' n <- 10; T <- 4
#' X <- replicate(T, array(rnorm(n*n*2), c(n, n, 2)), simplify = FALSE)
#' beta_true <- rbind(seq(-1, 1, length.out = T), seq(0.5, -0.5, length.out = T))
#' Y <- vector("list", T)
#' for (t in seq_len(T)) {
#'   Yt <- X[[t]][, , 1] * beta_true[1, t] + X[[t]][, , 2] * beta_true[2, t] +
#'         matrix(rnorm(n*n, 0, 0.2), n, n)
#'   diag(Yt) <- NA
#'   Y[[t]] <- Yt
#' }
#' fit_als <- als_dynamic_beta(Y, X, lambda = 0)     # per-period LS
#' fit_smooth <- als_dynamic_beta(Y, X, lambda = 10) # smoother path
#'
#' @export
als_dynamic_beta <- function(Y, Xdyad, lambda = 0, intercept = TRUE) {
	if (lambda < 0) cli::cli_abort("{.arg lambda} must be non-negative.")

	# normalise y to list-of-matrices
	if (is.array(Y) && length(dim(Y)) == 3L) {
		Y <- lapply(seq_len(dim(Y)[3L]), function(t) Y[, , t])
	}
	if (!is.list(Y)) cli::cli_abort("{.arg Y} must be a list of T matrices.")
	T_per <- length(Y)
	if (T_per < 2L) cli::cli_abort("Need at least two periods.")

	# normalise xdyad to list-of-arrays
	if (is.array(Xdyad) && length(dim(Xdyad)) == 4L) {
		Xdyad <- lapply(seq_len(dim(Xdyad)[4L]), function(t) Xdyad[, , , t])
	}
	if (length(Xdyad) != T_per) {
		cli::cli_abort("{.arg Xdyad} must have one element per time period.")
	}
	# per-period (y_t, x_t) on observed dyads
	y_list <- list(); X_list <- list()
	for (t in seq_len(T_per)) {
		Yt <- Y[[t]]; Xt <- Xdyad[[t]]
		if (length(dim(Xt)) == 2L) Xt <- array(Xt, c(dim(Xt), 1L))
		ok <- !is.na(Yt)
		y_obs <- as.numeric(Yt[ok])
		p_in <- dim(Xt)[3L]
		X_obs <- matrix(NA_real_, nrow = length(y_obs), ncol = p_in)
		for (k in seq_len(p_in)) X_obs[, k] <- Xt[, , k][ok]
		if (intercept) X_obs <- cbind(intercept = 1, X_obs)
		y_list[[t]] <- y_obs
		X_list[[t]] <- X_obs
	}
	p <- ncol(X_list[[1L]])
	cov_names <- colnames(X_list[[1L]]) %||% paste0("v", seq_len(p))

	# build the block-diagonal normal equations:
	#   h = blockdiag(x_t' x_t) + lambda * d' d
	#   g = blockdiag(x_t' y_t)
	# d is the first-difference operator on the stacked beta = (beta_1, ..., beta_t)
	H <- matrix(0, p * T_per, p * T_per)
	g <- numeric(p * T_per)
	for (t in seq_len(T_per)) {
		idx <- ((t - 1L) * p + 1L):(t * p)
		Xt <- X_list[[t]]; yt <- y_list[[t]]
		H[idx, idx] <- crossprod(Xt)
		g[idx] <- crossprod(Xt, yt)
	}
	if (lambda > 0) {
		# add lambda * d'd to h where d differences adjacent periods
		Dblock <- diag(p)
		for (t in seq_len(T_per - 1L)) {
			i1 <- ((t - 1L) * p + 1L):(t * p)
			i2 <- (t * p + 1L):((t + 1L) * p)
			# d row block: [-i, i] => d'd contributes i to (t,t) and (t+1,t+1),
			# -i to (t,t+1) and (t+1,t)
			H[i1, i1] <- H[i1, i1] + lambda * Dblock
			H[i2, i2] <- H[i2, i2] + lambda * Dblock
			H[i1, i2] <- H[i1, i2] - lambda * Dblock
			H[i2, i1] <- H[i2, i1] - lambda * Dblock
		}
	}
	# tiny ridge for numerical safety when h is rank-deficient
	H <- H + diag(1e-10, nrow(H))

	beta_stack <- tryCatch(
		solve(H, g),
		error = function(e) MASS::ginv(H) %*% g)
	beta_mat <- matrix(beta_stack, nrow = p, ncol = T_per,
	                   dimnames = list(cov_names, paste0("t", seq_len(T_per))))

	# residual sum of squares
	rss <- sum(vapply(seq_len(T_per), function(t) {
		r <- y_list[[t]] - X_list[[t]] %*% beta_mat[, t]
		sum(r^2)
	}, numeric(1L)))

	res <- list(
		beta         = beta_mat,
		lambda       = lambda,
		residual_ss  = rss,
		intercept    = intercept,
		call         = match.call(),
		label        = "penalised ALS estimate"
	)
	class(res) <- c("als_dynamic_beta", "list")
	res
}

#' Print method for penalised ALS time-varying beta
#' @param x An \code{als_dynamic_beta} result.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for its side effect of printing a
#'   summary of the penalised ALS time-varying coefficient estimates.
#' @export
#' @method print als_dynamic_beta
print.als_dynamic_beta <- function(x, digits = 4, ...) {
	cli::cli_h2("Penalised ALS time-varying beta")
	cli::cli_ul(c(
		"lambda: {.val {x$lambda}}",
		"residual SS: {.val {signif(x$residual_ss, 5)}}",
		"label: {.emph {x$label}} (NOT a posterior summary)"
	))
	print(round(x$beta, digits))
	invisible(x)
}

#' Extract beta path from a penalised-ALS object
#' @param object An \code{als_dynamic_beta} object.
#' @param ... Ignored.
#' @return The \eqn{p \times T} matrix of estimates.
#' @method coef als_dynamic_beta
#' @export
coef.als_dynamic_beta <- function(object, ...) {
	object$beta
}
