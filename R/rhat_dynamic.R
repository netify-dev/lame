# multivariate gelman-rubin diagnostic for dynamic_beta paths.
# treats the length-t path at one iteration as a multivariate
# observation and returns the brooks-gelman statistic per coefficient.

#' Multivariate split-R-hat for dynamic_beta coefficient paths
#'
#' Given one or more multi-chain \code{lame} fits whose \code{$BETA}
#' is \code{[n_iter, p, T]}, computes the multivariate R-hat per
#' coefficient \eqn{k} treating the length-\eqn{T} path as one
#' multivariate observation per iteration.
#'
#' @param fit_list A list of fitted \code{lame} objects from
#'   \code{lame_parallel(..., chains = K)}, or any \code{ame_chain_list}
#'   produced by re-running \code{lame()} with different seeds.
#' @param coefs Optional character vector of coefficient names to
#'   subset (matches \code{dimnames(fit$BETA)[[2]]}).
#' @return Data frame with one row per coefficient: \code{coef},
#'   \code{rhat_mvt}, \code{rhat_max_univariate}, \code{n_chains},
#'   \code{n_iter_per_chain}, \code{T}. \code{rhat_mvt} is the
#'   Brooks-Gelman multivariate statistic; \code{rhat_max_univariate}
#'   is the max over per-(k,t) split-R-hat values for comparison.
#'
#' @details For chain \eqn{c}, let \eqn{\beta^{(c)}_t \in \mathbb{R}^T}
#'   be the path. With \eqn{m} chains and \eqn{n} iterations per
#'   chain, define
#'   \deqn{W_k = \tfrac{1}{m}\sum_c S_c^{(k)}, \quad
#'         B_k = \tfrac{n}{m-1}\sum_c (\bar\beta_c^{(k)} - \bar\beta^{(k)})(\bar\beta_c^{(k)} - \bar\beta^{(k)})'}
#'   where \eqn{S_c^{(k)}} is the within-chain sample covariance of
#'   path \eqn{k} in chain \eqn{c}. Then
#'   \eqn{V_k = ((n-1)/n) W_k + ((m+1)/(mn)) B_k} and
#'   \eqn{\hat R_k^{mvt} = \sqrt{\lambda_\max(W_k^{-1} V_k)}}.
#'   For nearly degenerate covariances we add a tiny ridge to \eqn{W_k}.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit_list <- lapply(c(1L, 2L, 3L, 4L), function(s) {
#'   set.seed(s)
#'   lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'        dynamic_beta = "dyad",
#'        nscan = 200, burn = 50, odens = 5, verbose = FALSE)
#' })
#' rhat_dynamic_beta(fit_list)
#' }
#'
#' @export
rhat_dynamic_beta <- function(fit_list, coefs = NULL) {
	# common slip: passing a single fit object directly (either an `ame`
	# or `lame` fit). detect this and point the user at ame_parallel() or
	# running lame() under multiple seeds. the r-hat statistic is defined
	# only across chains so it cannot be computed from one fit.
	if (inherits(fit_list, "lame") || inherits(fit_list, "ame")) {
		cli::cli_abort(c(
			"{.arg fit_list} appears to be a single fit; R-hat needs at least two chains.",
			"i" = "Run {.code ame_parallel(..., n_chains = 4)} to produce a chain list, or",
			"i" = "Fit {.fn lame} under several seeds and wrap them as {.code list(fit1, fit2, fit3, fit4)}."))
	}
	if (!is.list(fit_list) || length(fit_list) < 2L) {
		cli::cli_abort(c(
			"{.arg fit_list} must be a list of at least two {.cls lame} fits (one per chain).",
			"i" = "Run {.code ame_parallel(..., n_chains = 4)} to produce a chain list, or",
			"i" = "Fit {.fn lame} under several seeds and wrap them as {.code list(fit1, fit2, fit3, fit4)}."))
	}
	if (!all(vapply(fit_list, function(f) inherits(f, "lame"), logical(1L)))) {
		cli::cli_abort(c(
			"All elements of {.arg fit_list} must be {.cls lame} fits.",
			"i" = "Found classes: {.cls {unique(unlist(lapply(fit_list, class)))}}."))
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
	p      <- dims[[1L]][2L]
	T_per  <- dims[[1L]][3L]
	m      <- length(fit_list)

	coef_nms <- dimnames(fit_list[[1L]]$BETA)[[2L]]
	if (is.null(coef_nms)) coef_nms <- paste0("v", seq_len(p))
	if (is.null(coefs)) coefs <- coef_nms

	out <- data.frame(
		coef = coefs,
		rhat_mvt = NA_real_,
		rhat_max_univariate = NA_real_,
		n_chains = m,
		n_iter_per_chain = n_iter,
		T = T_per,
		stringsAsFactors = FALSE
	)

	for (i in seq_along(coefs)) {
		k <- match(coefs[i], coef_nms)
		if (is.na(k)) next

		# extract [m chains] x [n iter] x [t] for this coef
		paths <- vapply(fit_list, function(f) f$BETA[, k, , drop = TRUE],
		                matrix(0, n_iter, T_per))
		# paths is [n_iter, t, m]
		if (length(dim(paths)) == 2L) {
			# t == 1 collapsed
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
		# small ridge to guarantee invertibility on degenerate problems
		W <- W + diag(1e-10, T_per)

		grand <- colMeans(chain_means)
		dev <- chain_means - matrix(grand, m, T_per, byrow = TRUE)
		B <- (n_iter / (m - 1L)) * crossprod(dev)

		V <- ((n_iter - 1L) / n_iter) * W + ((m + 1L) / (m * n_iter)) * B
		eig <- max(Re(eigen(solve(W, V), only.values = TRUE)$values))
		out$rhat_mvt[i] <- sqrt(eig)

		# univariate per-(k, t) for comparison: split-r-hat manually
		rhats_t <- numeric(T_per)
		for (t in seq_len(T_per)) {
			x_ct <- paths[, t, , drop = TRUE]
			rhats_t[t] <- .split_rhat_univariate(x_ct)
		}
		out$rhat_max_univariate[i] <- max(rhats_t, na.rm = TRUE)
	}
	out
}

#' @noRd
.split_rhat_univariate <- function(x_mat) {
	# x_mat: [n_iter, n_chains]
	if (is.null(dim(x_mat))) return(NA_real_)
	n <- nrow(x_mat); m <- ncol(x_mat)
	if (n < 4L || m < 2L) return(NA_real_)
	# split each chain in half (split-r-hat)
	half <- floor(n / 2L)
	splits <- vector("list", 2L * m)
	for (c in seq_len(m)) {
		splits[[2L * c - 1L]] <- x_mat[seq_len(half), c]
		splits[[2L * c]]      <- x_mat[(half + 1L):(2L * half), c]
	}
	M <- 2L * m
	N <- half
	chain_means <- vapply(splits, mean, numeric(1L))
	chain_vars  <- vapply(splits, stats::var,  numeric(1L))
	W <- mean(chain_vars)
	B <- N * stats::var(chain_means)
	V <- ((N - 1L) / N) * W + ((1L) / N) * B
	if (W <= 0 || !is.finite(W)) return(NA_real_)
	sqrt(V / W)
}
