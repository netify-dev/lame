# jointly rescale both factor paths, their innovation scale and the interaction; the
# likelihood is unchanged and the ratio carries the full Jacobian
.lame_bipartite_scale_log_ratio <- function(delta, sigma, G, shape, scale, g_var = 1) {
	if (delta == 0) return(0)
	log_abs_expm1 <- function(x) {
		if (x > log(2)) x + log1p(-exp(-x)) else if (x > 0)
			log(expm1(x)) else log(-expm1(x))
	}
	g_max <- max(abs(G))
	log_g_square <- if (g_max == 0) -Inf else
		2 * log(g_max) + log(sum((G / g_max)^2)) - log(g_var)
	terms <- c(
		log(scale) - 2 * log(sigma) + log_abs_expm1(-2 * delta),
		log_g_square - log(2) + log_abs_expm1(-4 * delta)
	)
	largest <- max(terms)
	magnitude <- if (largest == Inf) Inf else if (largest == -Inf) 0 else
		exp(largest) * sum(exp(terms - largest))
	-2 * (shape + length(G)) * delta + sign(delta) * magnitude
}

.lame_bipartite_scale_step <- function(
	U,
	V,
	G,
	sigma,
	rho,
	shape,
	scale,
	max_abs = Inf,
	max_G = Inf,
	sigma_max = Inf,
	stat_var_cap = Inf,
	step_sd = 0.5,
	g_var = 1
) {
	stopifnot(
		length(dim(U)) == 3L,
		length(dim(V)) == 3L,
		dim(U)[3L] == dim(V)[3L],
		identical(dim(G), c(dim(U)[2L], dim(V)[2L])),
		all(is.finite(U)),
		all(is.finite(V)),
		all(is.finite(G)),
		length(sigma) == 1L,
		is.finite(sigma),
		sigma > 0,
		length(rho) == 1L,
		is.finite(rho),
		abs(rho) < 1,
		is.finite(shape),
		shape > 0,
		is.finite(scale),
		scale > 0,
		is.finite(step_sd),
		step_sd > 0,
		all(!is.na(c(max_abs, max_G, sigma_max, stat_var_cap))),
		all(c(max_abs, max_G, sigma_max, stat_var_cap) > 0)
	)
	result <- list(
		U = U,
		V = V,
		G = G,
		sigma = sigma,
		attempts = 1L,
		accepts = 0L,
		bound_rejections = 0L,
		numerical_rejections = 0L
	)
	delta <- stats::rnorm(1L, sd = step_sd)
	factor <- exp(delta)
	up <- U * factor
	vp <- V * factor
	gp <- (G / factor) / factor
	sp <- sigma * factor
	if (!is.finite(factor) || factor <= 0 || !is.finite(sp) || sp <= 0 ||
			any(!is.finite(up)) || any(!is.finite(vp)) || any(!is.finite(gp)) ||
			any(U != 0 & up == 0) || any(V != 0 & vp == 0) ||
			any(G != 0 & gp == 0)) {
		result$numerical_rejections <- 1L
		return(result)
	}
	if (max(abs(up), abs(vp)) > max_abs || max(abs(gp)) > max_G ||
			sp > sigma_max || 2 * log(sp) - log1p(-rho^2) > log(stat_var_cap)) {
		result$bound_rejections <- 1L
		return(result)
	}
	ratio <- .lame_bipartite_scale_log_ratio(delta, sigma, G, shape, scale, g_var)
	if (is.na(ratio)) {
		result$numerical_rejections <- 1L
		return(result)
	}
	if (log(stats::runif(1L)) < ratio) {
		result$U <- up
		result$V <- vp
		result$G <- gp
		result$sigma <- sp
		result$accepts <- 1L
	}
	result
}

# keep numerical failures distinct from rejections at prior bounds
.lame_report_bipartite_scale_health <- function(fit) {
	move <- fit[["bipartite_scale_move"]]
	if (!is.null(move) && isTRUE(unname(move["numerical_rejections"]) > 0)) {
		cli::cli_abort(
			c(
				"Bipartite scale proposals failed numerically on {unname(move['numerical_rejections'])} of {unname(move['attempts'])} attempts.",
				"i" = "A proposal that fails is rejected and the chain keeps its state, so it is not a draw from the posterior and no fit is returned. Prior-bound rejections are counted separately."
			),
			class = "lame_sampler_failure"
		)
	}
	invisible(move)
}
