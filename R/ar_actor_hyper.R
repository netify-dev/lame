# shared AR hyperparameters for actor deviations on their constrained subspace:
# each centred margin removes one dimension from every period's Gaussian normaliser
.lame_ab_hyper_update <- function(
	a_mat,
	b_mat,
	rho_ab,
	sigma_ab,
	prior,
	symmetric = FALSE,
	center = FALSE,
	rvar = TRUE,
	cvar = TRUE,
	stat_var_cap = Inf
) {
	paths <- list()
	if (isTRUE(rvar)) paths[[length(paths) + 1L]] <- a_mat
	if (!isTRUE(symmetric) && isTRUE(cvar)) paths[[length(paths) + 1L]] <- b_mat
	if (!length(paths)) return(list(rho_ab = rho_ab, sigma_ab = sigma_ab))
	N <- ncol(paths[[1L]])
	if (N < 2L) return(list(rho_ab = rho_ab, sigma_ab = sigma_ab))
	if (isTRUE(center) && any(vapply(paths, function(x) any(!is.finite(x)), logical(1))))
		stop("Centered actor hyperparameter updates require finite latent paths.", call. = FALSE)
	sigma_ab <- .lame_sample_sigma_ab(
		a_mat,
		b_mat,
		rho_ab,
		sigma_ab,
		symmetric,
		prior_shape = prior$sigma_ab_shape %||% 2,
		prior_scale = prior$sigma_ab_scale %||% 1,
		stat_var_cap = stat_var_cap,
		center = center,
		rvar = rvar,
		cvar = cvar
	)
	q1 <- p1 <- product <- lag_square <- 0
	for (x in paths) {
		first <- x[, 1L]
		first <- first[is.finite(first)]
		q1 <- q1 + sum(first^2)
		p1 <- p1 + length(first) - as.integer(isTRUE(center) && length(first) > 0L)
		lag <- x[, -N, drop = FALSE]
		current <- x[, -1L, drop = FALSE]
		ok <- is.finite(lag) & is.finite(current)
		product <- product + sum(lag[ok] * current[ok])
		lag_square <- lag_square + sum(lag[ok]^2)
	}
	prior_mean <- prior$rho_ab_mean %||% .8
	prior_precision <- 1 / (prior$rho_ab_sd %||% .15)^2
	variance <- 1 / (lag_square / sigma_ab^2 + prior_precision)
	mean <- variance * (product / sigma_ab^2 + prior_mean * prior_precision)
	rho_ab <- .lame_rho_uv_mh_stationary(
		rho_ab,
		mean,
		variance,
		q1,
		p1,
		sigma_ab,
		stat_var_cap = stat_var_cap
	)
	list(rho_ab = rho_ab, sigma_ab = sigma_ab)
}
