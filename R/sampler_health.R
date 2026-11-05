# split each step's failure count into numerical failures and prior-bound rejections
.lame_sampler_health <- function(counts, uv_bounds = 0L, g_bounds = 0L) {
	total <- unlist(counts, use.names = TRUE)
	bounds <- setNames(numeric(length(total)), names(total))
	if ('UV' %in% names(bounds)) bounds['UV'] <- uv_bounds
	if ('G' %in% names(bounds)) bounds['G'] <- g_bounds
	if (any(
		!is.finite(total) | total < 0 | !is.finite(bounds) |
			bounds < 0 | bounds > total
	)) {
		stop('Invalid sampler failure and prior-bound counts.', call. = FALSE)
	}
	data.frame(
		step = names(total),
		numerical_failures = unname(total - bounds),
		prior_bound_rejections = unname(bounds),
		stringsAsFactors = FALSE
	)
}

# abort when a conditional failed numerically: it kept its state, so the chain is not a posterior draw.
# prior-bound rejections are Metropolis rejections, counted in fit$sampler_health (see .lame_bound_share_line)
.lame_report_sampler_health <- function(health, total_iters) {
	bad <- health$numerical_failures > 0
	if (any(bad)) {
		labels <- paste0(health$step[bad], ' (', health$numerical_failures[bad], ')')
		cli::cli_abort(
			c(
				'The sampler failed numerically in {.val {labels}} of {total_iters} sweeps.',
				'i' = 'A conditional update that fails keeps its previous state, so the chain is not a draw from the posterior and no fit is returned. The failure counts are in the {.code health} field of the condition.'
			),
			class = 'lame_sampler_failure',
			health = health
		)
	}
	invisible(health)
}

# share of sweeps in which a latent proposal fell outside the prior box (|u| <= uv_max_abs, |G| <= G_max),
# as a plain statement; NULL below one percent of the sweeps
.lame_bound_share_line <- function(x, min_share = 0.01) {
	h <- x[['sampler_health']]; it <- x[['sampler_health_iterations']]
	if (!is.data.frame(h) || !is.numeric(it) || length(it) != 1L || !is.finite(it) || it <= 0) return(NULL)
	share <- setNames(h$prior_bound_rejections / it, h$step)
	share <- share[intersect(names(share), c('UV', 'G'))]
	share <- share[is.finite(share) & share >= min_share]
	if (!length(share)) return(NULL)
	box <- c(
		UV = 'the positions fell outside the prior box |u| <= prior$uv_max_abs', G = 'the interaction matrix fell outside the prior box |G| <= prior$G_max'
	)
	paste0(
		'A conditional proposal for ',
		paste0(box[names(share)], ' on ', round(100 * share), '% of sweeps', collapse = ', and for '),
		'; that conditional kept its state, as a Metropolis rejection does.'
	)
}

.lame_pool_sampler_diagnostics <- function(combined, chains) {
	move_fields <- c(
		'level_move', 'dev_move', 'marginal_move', 'dynamic_marginal_move',
		'nc_ar_move', 'translation_move', 'rho_collapsed_move', 'uv_scale_move',
		'bipartite_scale_move', 'actor_factor_slice', 'actor_move', 'ab_hyper_move', 'uv_g_scale_move'
	)
	fields <- c('sampler_health', 'sampler_health_iterations', 'tryErrorChecks', move_fields)
	records <- lapply(chains, function(f) {
		setNames(lapply(fields, function(nm) f[[nm]]), fields)
	})
	combined$sampler_diagnostics_by_chain <- records
	# each chain adapts its own step sizes, so move diagnostics stay per chain
	for (nm in move_fields) combined[[nm]] <- NULL
	hs <- lapply(chains, `[[`, 'sampler_health')
	complete <- all(vapply(hs, is.data.frame, logical(1)))
	combined['sampler_health'] <- list(NULL) # exact name, no partial match
	if (complete) {
		steps <- hs[[1L]]$step
		complete <- all(vapply(hs, function(h) setequal(h$step, steps), logical(1)))
		if (complete) {
			total <- hs[[1L]]
			for (nm in c('numerical_failures', 'prior_bound_rejections'))
				total[[nm]] <- Reduce(`+`, lapply(hs, function(h) h[[nm]][match(steps, h$step)]))
			combined$sampler_health <- total
		}
	}
	counts <- lapply(chains, `[[`, 'tryErrorChecks')
	named <- !any(vapply(counts, is.null, logical(1))) &&
		all(vapply(counts, function(x) setequal(names(x), names(counts[[1L]])), logical(1)))
	combined$tryErrorChecks <- if (named) {
		keys <- names(counts[[1L]])
		setNames(lapply(keys, function(nm) sum(vapply(counts, function(x) x[[nm]], numeric(1)))), keys)
	} else NULL
	combined$mh_counters <- combined$tryErrorChecks
	iterations <- lapply(chains, `[[`, 'sampler_health_iterations')
	known <- all(vapply(iterations, function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, logical(1)))
	combined$sampler_health_iterations <- if (known) sum(unlist(iterations)) else NULL
	combined$refusal_share <- if (known && named && all(c('UV', 'G') %in% names(combined$tryErrorChecks))) {
		unlist(combined$tryErrorChecks[c('UV', 'G')]) / combined$sampler_health_iterations
	} else NULL
	combined
}

# abort on a numerical failure in either joint Hamiltonian kernel, which keep their own counters
.lame_report_joint_sampler_health <- function(fit) {
	fields <- c('marginal_move', 'dynamic_marginal_move')
	failures <- vapply(fields, function(nm) {
		move <- fit[[nm]]
		if (is.null(move) || !'numerical_rejections' %in% names(move)) return(0)
		unname(move['numerical_rejections'])
	}, numeric(1))
	bad <- failures > 0
	if (any(bad)) {
		labels <- paste0(fields[bad], ' (', failures[bad], ')')
		cli::cli_abort(
			c(
				'Joint Hamiltonian trajectories failed numerically: {.val {labels}}.',
				'i' = 'A trajectory that fails is rejected and the chain keeps its state, so it is not a draw from the posterior and no fit is returned. Valid prior-bound rejections have separate counts.'
			),
			class = 'lame_sampler_failure',
			failures = failures
		)
	}
	invisible(failures)
}

# quasi-separation check for binary probit fits: a fitted predictor beyond |8| (pnorm(-8) < 1e-15)
# is a probability the data cannot tell from 0 or 1, so its size comes from the prior.
# warn when such dyads are a sizeable share of the observed ones
.lame_separation_check <- function(fit, eta_max = 8, share_max = 0.05) {
	if (!identical(fit[['family', exact = TRUE]], 'binary')) return(invisible(NA_real_))
	EZ <- tryCatch(reconstruct_EZ(fit), error = function(e) NULL)
	if (is.null(EZ)) return(invisible(NA_real_))
	Y <- fit[['Y', exact = TRUE]]
	if (is.list(EZ)) {
		# lame stores one matrix per period, holding only that period's actors
		# when the actor set changes; Y is the cube over all actors
		eta <- unlist(lapply(seq_along(EZ), function(t) {
			y <- if (length(dim(Y)) == 3L && dim(Y)[3L] == length(EZ)) Y[, , t] else NULL
			.lame_observed_eta(EZ[[t]], y)
		}), use.names = FALSE)
	} else {
		if (length(dim(EZ)) == 3L) EZ <- matrix(EZ, nrow = dim(EZ)[1L])
		eta <- .lame_observed_eta(EZ, if (is.matrix(Y)) Y else NULL)
	}
	eta <- eta[is.finite(eta)]
	if (length(eta) == 0L) return(invisible(NA_real_))
	share <- mean(abs(eta) > eta_max)
	if (share >= share_max) {
		R_fit <- max(
			c(
				fit[['R', exact = TRUE]],
				fit[['R_row', exact = TRUE]],
				fit[['R_col', exact = TRUE]],
				0
			),
			na.rm = TRUE
		)
		pct <- sprintf('%.0f%%', 100 * share)
		hint <- if (R_fit > 0) {
			'The multiplicative term of rank {R_fit} reproduces those ties; a smaller {.arg R} (down to {.code R = 0} when few ties change over time) reproduces them without saturating.'
		} else {
			'A covariate, actor or block whose ties are all present or all absent separates the fit; check for one.'
		}
		cli::cli_warn(
			c(
				'The fit has separated: the fitted linear predictor exceeds {eta_max} in absolute value on {pct} of the observed dyads.',
				'i' = 'Those probabilities are numerically 0 or 1, so the data do not fix the size of the predictor on them: its posterior is its prior. The signs and the fitted probabilities are determined by the data; the size of the coefficients, effects and positions behind them is not.',
				'i' = hint
			),
			class = 'lame_separation'
		)
	}
	invisible(share)
}

# prior sd of the Xdyad coefficients `k` in ame() (Zellner g-prior when unipartite, N(0, g s20) when
# bipartite), for the separation message; NA when the design cannot be formed
.ame_separation_prior_sd <- function(Y, Xdyad, Xrow, Xcol, intercept, mode, g, prior, k) {
	tryCatch({
		bip <- identical(mode, 'bipartite')
		g_eff <- if (length(g) == 1L && is.na(g)) {
			if (bip) as.numeric(nrow(Y)) * ncol(Y) else sum(!is.na(Y))
		} else g[1L]
		if (bip) return(rep(sqrt(g_eff * (prior[['s20']] %||% 1)), length(k)))
		n <- nrow(Y)
		X <- precomputeX(design_array(Xrow, Xcol, Xdyad, intercept, n, warn = FALSE))
		V0 <- solve(.beta_prior_precision(attributes(X)$XX, attributes(X)$mX, g_eff))
		unname(sqrt(diag(V0))[dim(X)[3L] - dim(Xdyad)[3L] + k])
	}, error = function(e) rep(NA_real_, length(k)))
}

# fitted linear predictor on the observed cells of y, matching actors by
# name when the fitted matrix covers only some of y's actors
.lame_observed_eta <- function(e, y) {
	if (!is.null(y) && !is.null(rownames(e)) && !is.null(colnames(e)) &&
			all(rownames(e) %in% rownames(y)) && all(colnames(e) %in% colnames(y))) {
		y <- y[rownames(e), colnames(e), drop = FALSE]
	}
	if (!is.null(y) && identical(dim(y), dim(e))) return(e[!is.na(y)])
	if (nrow(e) == ncol(e)) e[row(e) != col(e)] else as.vector(e)
}
