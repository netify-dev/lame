# a coefficient-free normal model does not use g
.ame_check_normal_g <- function(g, family, p) {
	if (!identical(family, "normal") || p == 0L) return(invisible(NULL))
	if (!is.numeric(g) || length(g) != 1L || !is.finite(g) || g <= 0) {
		cli::cli_abort(
			c(
				"{.arg g} must be a finite positive numeric scalar for normal AME with regression coefficients.",
				"i" = "The data-derived default can be zero for constant normal responses; supply an explicit positive {.arg g}."
			),
			class = "lame_invalid_normal_g"
		)
	}
	invisible(NULL)
}

# refuse an improper posterior for rho when every observed reciprocal pair is
# exactly equal: the mean absorbs at most `bound` of the m pair contrasts, and
# the density of rho is not integrable at one once m - bound reaches
# .lame_rho_improper_dyads()
.ame_check_normal_reciprocity <- function(Y, X, rvar, cvar, dcor, symmetric, R, family) {
	if (!identical(family, "normal") || isTRUE(symmetric) || !isTRUE(dcor)) {
		return(invisible(NULL))
	}
	pairs <- which(upper.tri(Y) & !is.na(Y) & !is.na(t(Y)), arr.ind = TRUE)
	m <- nrow(pairs)
	if (m == 0L) return(invisible(NULL))
	all_equal <- all(Y[pairs] == t(Y)[pairs])
	# unequal reciprocal responses invalidate the higher-rank boundary proof
	if (R > 0 && !all_equal) return(invisible(NULL))
	n <- nrow(Y)
	# a dense asymmetric pair set fits neither the refusal nor the forest safeguard
	if (!all_equal && m > n - 1L) return(invisible(NULL))
	component_count <- .lame_pair_components(pairs, n)
	actor_bound <- if (isTRUE(rvar) || isTRUE(cvar)) n - component_count else 0L
	p <- dim(X)[3L]
	asymmetric_columns <- if (p == 0L) 0L else sum(vapply(seq_len(p), function(k) {
		x <- X[, , k]
		any(x[pairs] != t(x)[pairs])
	}, logical(1L)))
	rank_bound <- actor_bound + asymmetric_columns
	improper_from <- .lame_rho_improper_dyads()
	if (R > 0) {
		# shrinking all factor coordinates gives a conservative n*R volume power
		factor_bound <- n * R
		if (all_equal && m - (rank_bound + factor_bound) >= improper_from) {
			cli::cli_abort(
				c(
					"This directed normal model has an improper posterior under the reciprocity prior.",
					"i" = "All {m} observed reciprocal pairs are exactly equal; the actor/regression rank bound is {rank_bound}, and the factor bound is {factor_bound}. The posterior is improper once the pairs exceed the combined bound by {improper_from}.",
					"i" = "Use {.code symmetric = TRUE} for an undirected network, or {.code dcor = FALSE} for fixed zero reciprocity."
				),
				class = "lame_improper_normal_reciprocity"
			)
		}
		# a sufficient check only: it neither certifies propriety nor initializes factors
		return(invisible(NULL))
	}
	if (all_equal && m - rank_bound >= improper_from) {
		cli::cli_abort(
			c(
				"This directed normal model has an improper posterior under the reciprocity prior.",
				"i" = "All {m} observed reciprocal pairs are exactly equal, but their antisymmetric mean design has rank at most {rank_bound}; the posterior is improper once the pairs exceed that rank by {improper_from}.",
				"i" = "Use {.code symmetric = TRUE} for an undirected network, or {.code dcor = FALSE} for fixed zero reciprocity."
			),
			class = "lame_improper_normal_reciprocity"
		)
	}
	invisible(list(
		reciprocal_pairs = m,
		actor_rank_bound = actor_bound,
		antisymmetric_rank_bound = rank_bound
	))
}

# initialize a full-row-rank observed actor design with one exact Gaussian
# block draw; a design with cycles gets the same draw when the next rho
# conditional would be improper, otherwise the rho update reports the
# degenerate state
.ame_normal_reciprocity_initial_draw <- function(
	Z,
	EZ,
	Sab,
	rho,
	X,
	s2,
	UV,
	g,
	certificate,
	rvar,
	cvar
) {
	if (is.null(certificate) || certificate$reciprocal_pairs == 0L) return(NULL)
	forest <- certificate$actor_rank_bound == certificate$reciprocal_pairs
	residual <- Z - EZ
	ut <- upper.tri(residual)
	e1 <- residual[ut]
	e2 <- t(residual)[ut]
	# the actual complete latent pairs used by rrho_mh must be exactly degenerate
	if (!length(e1) || any(!is.finite(e1)) || any(!is.finite(e2)) ||
			!(all(e1 == e2) || all(e1 == -e2))) return(NULL)
	if (!forest && length(e1) < .lame_rho_improper_dyads()) return(NULL)
	active <- Sab
	if (!isTRUE(rvar)) active[1L, ] <- active[, 1L] <- 0
	if (!isTRUE(cvar)) active[2L, ] <- active[, 2L] <- 0
	direction <- c(1, -1)
	skew_variance <- drop(t(direction) %*% active %*% direction)
	if (!is.finite(skew_variance) || skew_variance <= 0) {
		if (!forest) return(NULL)
		cli::cli_abort(
			"The initial reciprocity refresh requires a finite positive actor contrast prior variance.",
			class = "lame_normal_reciprocity_initialization"
		)
	}
	if (dim(X)[3L] > 0L && (!is.numeric(g) || length(g) != 1L || !is.finite(g) || g <= 0)) {
		if (!forest) return(NULL)
		cli::cli_abort(
			"The initial reciprocity refresh requires a finite positive {.arg g}; supply {.arg g} explicitly for constant normal responses.",
			class = "lame_normal_reciprocity_initialization"
		)
	}
	# one exact draw only: any still-degenerate next rho conditional fails normally
	rbeta_ab_fc(Z, active, rho, X, s2, offset = UV, g = g)
}

# connected components of the actors joined by the pairs (two-column index matrix)
.lame_pair_components <- function(pairs, n) {
	parent <- seq_len(n)
	find_root <- function(actor) {
		while (parent[actor] != actor) {
			parent[actor] <<- parent[parent[actor]]
			actor <- parent[actor]
		}
		actor
	}
	comps <- n
	for (k in seq_len(nrow(pairs))) {
		a <- find_root(pairs[k, 1L]); b <- find_root(pairs[k, 2L])
		if (a != b) { parent[b] <- a; comps <- comps - 1L }
	}
	comps
}

# propriety bound of a directed normal panel whose reciprocal pairs are all
# exactly equal: near one the density of rho behaves as (1 - rho)^(-(M - r) / 2),
# M the pairs of the panel and r the rank of the pair contrasts over all mean
# parameters at once (a shared parameter counts once), so the posterior is
# improper from M - r = .lame_rho_improper_dyads(). r is an upper bound, so a
# proper panel is never refused; a pooled panel can be improper even when each
# period alone is proper
.lame_panel_reciprocity_bound <- function(
	Y,
	Xlist,
	rvar,
	cvar,
	R,
	dynamic_ab = FALSE,
	dynamic_uv = FALSE,
	dynamic_beta = FALSE
) {
	n <- dim(Y)[1L]; N <- dim(Y)[3L]
	pairs <- lapply(seq_len(N), function(t) {
		y <- Y[, , t]
		which(upper.tri(y) & !is.na(y) & !is.na(t(y)), arr.ind = TRUE)
	})
	M <- sum(vapply(pairs, nrow, 1L))
	actor <- if (!(isTRUE(rvar) || isTRUE(cvar))) 0L else if (isTRUE(dynamic_ab)) {
		sum(vapply(pairs, function(p) n - .lame_pair_components(p, n), 1L))
	} else n - .lame_pair_components(unique(do.call(rbind, pairs)), n)
	p <- if (length(Xlist)) dim(Xlist[[1L]])[3L] else 0L
	asym <- matrix(FALSE, N, p)
	for (t in seq_len(N)) for (k in seq_len(p)) {
		x <- Xlist[[t]][, , k]
		asym[t, k] <- any(x[pairs[[t]]] != t(x)[pairs[[t]]])
	}
	coef <- if (isTRUE(dynamic_beta)) sum(asym) else sum(colSums(asym) > 0L)
	factor <- if (R > 0) n * R * (if (isTRUE(dynamic_uv)) N else 1L) else 0L
	list(
		pairs = M,
		actor = actor,
		coef = as.integer(coef),
		factor = as.integer(factor),
		rank = actor + as.integer(coef) + as.integer(factor),
		improper = M - (actor + coef + factor) >= .lame_rho_improper_dyads()
	)
}
