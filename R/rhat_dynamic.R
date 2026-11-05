# multivariate Gelman-Rubin diagnostic for dynamic_beta paths: the length-T path
# at one iteration is one multivariate draw (Brooks and Gelman 1998)

#' Multivariate R-hat for dynamic_beta coefficient paths
#'
#' Given one or more multi-chain \code{lame} fits whose \code{$BETA}
#' is \code{[n_iter, p, T]}, computes the multivariate R-hat per
#' coefficient \eqn{k} treating the length-\eqn{T} path as one
#' multivariate observation per iteration.
#'
#' @param fit_list A list of fitted \code{lame} objects from
#'   \code{lame_parallel(..., n_chains = K, combine_method = "list")}, or any \code{ame_chain_list}
#'   produced by re-running \code{lame()} with different seeds.
#' @param coefs Optional character vector of coefficient names to
#'   subset (matches \code{dimnames(fit$BETA)[[2]]}).
#' @return Data frame with one row per coefficient: \code{coef},
#'   \code{rhat_mvt}, \code{rhat_max_univariate}, \code{n_chains},
#'   \code{n_iter_per_chain}, \code{n_periods}. \code{rhat_mvt} is the
#'   Brooks-Gelman multivariate statistic; \code{rhat_max_univariate}
#'   is the max over periods of the split rank-normalised R-hat of each
#'   \eqn{(k, t)} series, the statistic of \code{fit$rhat} and of the
#'   between-chain gate of \code{\link{lame}}.
#'
#' @details For chain \eqn{c}, let \eqn{\beta^{(c)}_t \in \mathbb{R}^T}
#'   be the path. With \eqn{m} chains and \eqn{n} iterations per
#'   chain, define
#'   \deqn{W_k = \tfrac{1}{m}\sum_c S_c^{(k)}, \quad
#'         B_k = \tfrac{n}{m-1}\sum_c (\bar\beta_c^{(k)} - \bar\beta^{(k)})(\bar\beta_c^{(k)} - \bar\beta^{(k)})'}
#'   where \eqn{S_c^{(k)}} is the within-chain sample covariance of
#'   path \eqn{k} in chain \eqn{c}. Then
#'   \eqn{V_k = ((n-1)/n) W_k + ((m+1)/(mn)) B_k} and
#'   \eqn{\hat R_k^{mvt} = \sqrt{\lambda_{\max}(W_k^{-1} V_k)}}.
#'   For nearly degenerate covariances we add a tiny ridge to \eqn{W_k}.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' # note: pass seed = to lame() -- an external set.seed() does not vary
#' # the sampler, so it would produce identical chains
#' # each fit is one chain (chains = 1); rhat_dynamic_beta() compares them
#' fit_list <- lapply(c(1L, 2L), function(s) {
#'   lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X, family = "binary", R = 0,
#'        dynamic_beta = "dyad", seed = s,
#'        nscan = 60, burn = 15, odens = 5, verbose = FALSE, chains = 1)
#' })
#' rhat_dynamic_beta(fit_list)
#' }
#'
#' @export
rhat_dynamic_beta <- function(fit_list, coefs = NULL) {
	# R-hat needs several chains, so a single fit is refused
	if (inherits(fit_list, "lame") || inherits(fit_list, "ame")) {
		cli::cli_abort(c(
			"{.arg fit_list} appears to be a single fit; R-hat needs at least two chains.",
			"i" = "Run {.code lame_parallel(..., n_chains = 4, combine_method = \"list\")} to produce a chain list, or",
			"i" = "Fit {.fn lame} under several seeds and wrap them as {.code list(fit1, fit2, fit3, fit4)}."
		))
	}
	if (!is.list(fit_list) || length(fit_list) < 2L) {
		cli::cli_abort(c(
			"{.arg fit_list} must be a list of at least two {.cls lame} fits (one per chain).",
			"i" = "Run {.code lame_parallel(..., n_chains = 4, combine_method = \"list\")} to produce a chain list, or",
			"i" = "Fit {.fn lame} under several seeds and wrap them as {.code list(fit1, fit2, fit3, fit4)}."
		))
	}
	if (!all(vapply(fit_list, function(f) inherits(f, "lame"), logical(1L)))) {
		cli::cli_abort(c(
			"All elements of {.arg fit_list} must be {.cls lame} fits.",
			"i" = "Found classes: {.cls {unique(unlist(lapply(fit_list, class)))}}."
		))
	}
	dims <- lapply(fit_list, function(f) dim(f$BETA))
	if (!all(vapply(dims, length, integer(1L)) == 3L)) {
		cli::cli_abort("Every fit must have a 3-D {.code $BETA} (dynamic_beta).")
	}
	if (length(unique(vapply(dims, function(d) d[2L], integer(1L)))) > 1L ||
			length(unique(vapply(dims, function(d) d[3L], integer(1L)))) > 1L) {
		cli::cli_abort("All fits must share the same p and T.")
	}
	n_iter <- dims[[1L]][1L]
	p <- dims[[1L]][2L]
	T_per <- dims[[1L]][3L]
	m <- length(fit_list)
	# identical chains (usually one seed for every fit) make every R-hat look converged
	if (all(vapply(
		fit_list[-1L],
		function(f) identical(f$BETA, fit_list[[1L]]$BETA),
		logical(1L)
	))) {
		cli::cli_warn(c(
			"All fits in {.arg fit_list} have identical draws; R-hat is meaningless.",
			"i" = "Pass a different {.arg seed} to each {.fn lame} call ({.fn set.seed} outside the call does not vary the sampler)."
		))
	}

	coef_nms <- dimnames(fit_list[[1L]]$BETA)[[2L]]
	if (is.null(coef_nms)) coef_nms <- paste0("v", seq_len(p))
	if (is.null(coefs)) coefs <- coef_nms

	out <- data.frame(
		coef = coefs,
		rhat_mvt = NA_real_,
		rhat_max_univariate = NA_real_,
		n_chains = m,
		n_iter_per_chain = n_iter,
		n_periods = T_per,
		stringsAsFactors = FALSE
	)

	for (i in seq_along(coefs)) {
		k <- match(coefs[i], coef_nms)
		if (is.na(k)) next

		# draws of this coefficient as [n_iter, T, chains]
		paths <- vapply(
			fit_list,
			function(f) f$BETA[, k, , drop = TRUE],
			matrix(0, n_iter, T_per)
		)
		if (length(dim(paths)) == 2L) {
			# a single period collapses to a matrix
			paths <- array(paths, c(n_iter, 1L, m))
		}

		# within-chain covariance
		W <- matrix(0, T_per, T_per)
		chain_means <- matrix(NA_real_, m, T_per)
		for (c in seq_len(m)) {
			P_c <- paths[, , c]
			if (is.null(dim(P_c))) P_c <- matrix(P_c, ncol = 1L)
			chain_means[c, ] <- colMeans(P_c)
			W <- W + stats::cov(P_c)
		}
		W <- W / m
		# small ridge so W stays invertible
		W <- W + diag(1e-10, T_per)

		grand <- colMeans(chain_means)
		dev <- chain_means - matrix(grand, m, T_per, byrow = TRUE)
		B <- (n_iter / (m - 1L)) * crossprod(dev)

		V <- ((n_iter - 1L) / n_iter) * W + ((m + 1L) / (m * n_iter)) * B
		eig <- max(Re(eigen(solve(W, V), only.values = TRUE)$values))
		out$rhat_mvt[i] <- sqrt(eig)

		# largest univariate split rank-normalised R-hat over periods, as held in fit$rhat
		rhats_t <- numeric(T_per)
		for (t in seq_len(T_per)) {
			rhats_t[t] <- .lame_rhat_chains(matrix(paths[, t, ], n_iter, m))
		}
		out$rhat_max_univariate[i] <- if (all(is.na(rhats_t))) NA_real_ else max(rhats_t, na.rm = TRUE)
	}
	out
}
