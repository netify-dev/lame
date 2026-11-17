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
				.lame_failure_lines(health$step[bad]),
				'i' = 'A conditional update that fails keeps its previous state, so the chain is not a draw from the posterior and no fit is returned. The failure counts are in the {.code health} field of the condition.'
			),
			class = 'lame_sampler_failure',
			health = health,
			first_errors = .lame_fallback_env$messages[intersect(health$step[bad], names(.lame_fallback_env$messages))]
		)
	}
	invisible(health)
}

# exact fallbacks: a step whose numerical work failed and that kept its state where that leaves the
# posterior invariant (a block skipped given quantities it does not move, a proposal rejected), counted
# per step. they do not stop the fit; summary() states their share and a share above ten percent warns
.lame_fallback_env <- new.env(parent = emptyenv())

.lame_fallback_reset <- function() {
	.lame_fallback_env$counts <- numeric(0)
	.lame_fallback_env$messages <- list()
	invisible(NULL)
}

# the first error of a step that failed, with its sweep, for the abort that reports the failure
.lame_failure_note <- function(step, sweep, e) {
	msgs <- .lame_fallback_env$messages %||% list()
	if (is.null(msgs[[step]])) {
		msgs[[step]] <- list(sweep = as.integer(sweep), message = conditionMessage(e))
		.lame_fallback_env$messages <- msgs
	}
	invisible(NULL)
}

# one 'x' bullet per step with a recorded error; braces are escaped for cli
.lame_failure_lines <- function(steps) {
	msgs <- .lame_fallback_env$messages %||% list()
	steps <- intersect(steps, names(msgs))
	if (!length(steps)) return(character(0))
	lines <- vapply(steps, function(st) {
		m <- gsub('}', '}}', gsub('{', '{{', msgs[[st]]$message, fixed = TRUE), fixed = TRUE)
		paste0(st, ', first failure at sweep ', msgs[[st]]$sweep, ': ', m)
	}, character(1))
	setNames(lines, rep('x', length(lines)))
}

.lame_fallback_note <- function(step) {
	counts <- .lame_fallback_env$counts %||% numeric(0)
	counts[step] <- (if (step %in% names(counts)) counts[[step]] else 0) + 1
	.lame_fallback_env$counts <- counts
	invisible(NULL)
}

# the counts of the run, cleared for the next one; the first errors stay for the reports that follow
.lame_fallback_take <- function() {
	counts <- .lame_fallback_env$counts %||% numeric(0)
	.lame_fallback_env$counts <- numeric(0)
	counts
}

# warn when a step fell back more often than on `max_share` of the sweeps: the chain is valid, but that
# step did not move the state each time
.lame_report_fallbacks <- function(counts, total_iters, max_share = 0.1) {
	if (!length(counts) || !is.numeric(total_iters) || !is.finite(total_iters) || total_iters <= 0) return(invisible(counts))
	high <- counts / total_iters > max_share
	if (any(high)) {
		labels <- paste0(names(counts)[high], ' (', counts[high], ')')
		cli::cli_warn(
			c(
				'After a numerical failure, {.val {labels}} kept its state more often than on {round(100 * max_share)}% of the {total_iters} sweeps.',
				'i' = 'Keeping the state there is a valid move, so the draws are from the posterior, but that step did not move the chain each time; the counts are in {.code fit$exact_fallbacks}.'
			),
			class = 'lame_exact_fallbacks'
		)
	}
	invisible(counts)
}

# the exact fallbacks of a fit as a plain statement; NULL when there was none
.lame_fallback_share_line <- function(x) {
	counts <- x[['exact_fallbacks']]; it <- x[['sampler_health_iterations']]
	if (!is.numeric(counts) || is.null(names(counts)) || !any(counts > 0)) return(NULL)
	counts <- counts[counts > 0]
	sweeps <- if (is.numeric(it) && length(it) == 1L && is.finite(it) && it > 0) paste0(' in ', it, ' sweeps') else ''
	paste0(
		'After a numerical failure a step kept its state, a valid move, ',
		paste0(names(counts), ' ', counts, ' times', collapse = ', '),
		sweeps,
		' (fit$exact_fallbacks).'
	)
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
		'bipartite_scale_move', 'actor_factor_slice', 'actor_move', 'ab_hyper_move', 'uv_g_scale_move',
		'uv_split_move', 'uv_rotation_move', 'uv_rho_flip_move', 'ab_shift_move'
	)
	fields <- c('sampler_health', 'sampler_health_iterations', 'tryErrorChecks', 'exact_fallbacks', move_fields)
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
	fb <- lapply(chains, function(f) f[['exact_fallbacks']] %||% numeric(0))
	fb_steps <- unique(unlist(lapply(fb, names)))
	combined$exact_fallbacks <- setNames(
		vapply(fb_steps, function(nm) sum(vapply(fb, function(x) if (nm %in% names(x)) x[[nm]] else 0, numeric(1))), numeric(1)),
		fb_steps
	)
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
# is a probability the data cannot tell from 0 or 1, so its size comes from the prior, unless the
# covariates put it there through slopes the other dyads determine (.lame_pinned_slopes).
# warn when the dyads that remain are a sizeable share of the observed ones
.lame_separation_check <- function(fit, eta_max = 8, share_max = 0.05) {
	if (!identical(fit[['family', exact = TRUE]], 'binary')) return(invisible(NA_real_))
	EZ <- tryCatch(reconstruct_EZ(fit), error = function(e) {
		cli::cli_warn(
			c(
				'The separation check was not run: the fitted linear predictor could not be formed.',
				'x' = '{conditionMessage(e)}'
			),
			class = 'lame_separation_unchecked'
		)
		NULL
	})
	if (is.null(EZ)) return(invisible(NA_real_))
	Y <- fit[['Y', exact = TRUE]]
	if (is.list(EZ)) {
		# lame stores one matrix per period, holding only that period's actors
		# when the actor set changes; Y is the cube over all actors
		Xlist <- fit[['Xlist', exact = TRUE]]
		parts <- lapply(seq_along(EZ), function(t) {
			y <- if (length(dim(Y)) == 3L && dim(Y)[3L] == length(EZ)) Y[, , t] else NULL
			keep <- .lame_observed_cells(EZ[[t]], y)
			x <- if (is.list(Xlist) && length(Xlist) == length(EZ)) Xlist[[t]] else NULL
			list(eta = EZ[[t]][keep], x = .lame_design_rows(x, EZ[[t]], keep))
		})
	} else {
		if (length(dim(EZ)) == 3L) EZ <- matrix(EZ, nrow = dim(EZ)[1L])
		keep <- .lame_observed_cells(EZ, if (is.matrix(Y)) Y else NULL)
		parts <- list(list(eta = EZ[keep], x = .lame_design_rows(fit[['X', exact = TRUE]], EZ, keep)))
	}
	eta <- unlist(lapply(parts, `[[`, 'eta'), use.names = FALSE)
	# without a design or coefficient draws nothing is taken off the predictor
	pinned <- tryCatch(.lame_pinned_part(fit, parts), error = function(e) NULL)
	if (length(pinned) != length(eta)) pinned <- numeric(length(eta))
	ok <- is.finite(eta)
	if (!any(ok)) return(invisible(NA_real_))
	beyond <- abs(eta[ok]) > eta_max
	share <- mean(beyond & abs(eta[ok] - pinned[ok]) > eta_max)
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
				if (share < mean(beyond)) {
					'The fit has separated: on {pct} of the observed dyads the fitted linear predictor exceeds {eta_max} in absolute value beyond what the covariates account for.'
				} else {
					'The fit has separated: the fitted linear predictor exceeds {eta_max} in absolute value on {pct} of the observed dyads.'
				},
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
		g_def <- length(g) == 1L && is.na(g)
		g_eff <- if (g_def) {
			if (bip) as.numeric(nrow(Y)) * ncol(Y) else sum(!is.na(Y))
		} else g[1L]
		if (bip) {
			# default bipartite g of a slope per unit of its covariate, as in ame_bipartite()
			ms <- vapply(k, function(j) {
				x <- if (length(dim(Xdyad)) == 3L) Xdyad[, , j] else Xdyad
				m <- mean(x^2, na.rm = TRUE)
				if (g_def && is.finite(m) && m > 0) m else 1
			}, numeric(1))
			return(sqrt(g_eff / ms * (prior[['s20']] %||% 1)))
		}
		n <- nrow(Y)
		X <- precomputeX(design_array(Xrow, Xcol, Xdyad, intercept, n, warn = FALSE))
		V0 <- solve(.beta_prior_precision(attributes(X)$XX, attributes(X)$mX, g_eff))
		unname(sqrt(diag(V0))[dim(X)[3L] - dim(Xdyad)[3L] + k])
	}, error = function(e) rep(NA_real_, length(k)))
}

# observed cells of y on the fitted matrix e, matching actors by name when
# e covers only some of y's actors
.lame_observed_cells <- function(e, y) {
	if (!is.null(y) && !is.null(rownames(e)) && !is.null(colnames(e)) &&
			all(rownames(e) %in% rownames(y)) && all(colnames(e) %in% colnames(y))) {
		y <- y[rownames(e), colnames(e), drop = FALSE]
	}
	if (!is.null(y) && identical(dim(y), dim(e))) return(!is.na(y))
	if (nrow(e) == ncol(e)) row(e) != col(e) else matrix(TRUE, nrow(e), ncol(e))
}

# rows of the design array X on the cells `keep` of the fitted matrix e, one column per slice;
# NULL when X is not a design for e
.lame_design_rows <- function(X, e, keep) {
	if (length(dim(X)) != 3L || !identical(dim(X)[1:2], dim(e))) return(NULL)
	named <- !is.null(rownames(X)) && !is.null(rownames(e)) && !is.null(colnames(X)) && !is.null(colnames(e))
	if (named && !(identical(rownames(X), rownames(e)) && identical(colnames(X), colnames(e)))) return(NULL)
	p <- dim(X)[3L]
	matrix(X[rep(keep, p)], ncol = p, dimnames = list(NULL, dimnames(X)[[3L]]))
}

# covariate part of the predictor that the ties determine, on the dyads of `parts` (a list of
# list(eta, x), one per period, or one in all): static coefficients are judged on every period
# together, dynamic ones period by period. NULL when the fit carries no design or coefficient draws
.lame_pinned_part <- function(fit, parts) {
	B <- fit[['BETA', exact = TRUE]]
	if (is.null(B) || any(vapply(parts, function(p) is.null(p$x), logical(1)))) return(NULL)
	sym <- isTRUE(fit[['symmetric', exact = TRUE]])
	if (length(dim(B)) == 3L) {
		if (dim(B)[3L] != length(parts)) return(NULL)
		bm <- apply(B, c(2L, 3L), mean)
		out <- lapply(seq_along(parts), function(t) {
			b <- .beta_on_design(stats::setNames(bm[, t], dimnames(B)[[2L]]), colnames(parts[[t]]$x))
			if (length(b) != ncol(parts[[t]]$x)) return(NULL)
			.lame_pinned_slopes(parts[[t]]$eta, parts[[t]]$x, b, sym)
		})
		if (any(vapply(out, is.null, logical(1)))) return(NULL)
		return(unlist(out, use.names = FALSE))
	}
	x <- do.call(rbind, lapply(parts, `[[`, 'x'))
	b <- .beta_on_design(colMeans(B), colnames(x))
	if (length(b) != ncol(x)) return(NULL)
	.lame_pinned_slopes(unlist(lapply(parts, `[[`, 'eta'), use.names = FALSE), x, b, sym)
}

# slopes times covariates on each dyad, centred at the information-weighted mean of the covariates,
# where the ties determine that product to within `sd_max` on the probit scale; zero elsewhere.
# the information is that of the probit likelihood at the fitted predictor, the effects and positions
# held there: a dyad beyond the cut carries none, so a slope met only on such dyads is not determined,
# while a strong covariate with informative dyads along its range is
.lame_pinned_slopes <- function(eta, X, beta, symmetric = FALSE, sd_max = 1) {
	out <- numeric(length(eta))
	ok <- is.finite(eta) & rowSums(!is.finite(X)) == 0L
	if (!any(ok)) return(out)
	vary <- apply(X[ok, , drop = FALSE], 2L, function(x) max(x) > min(x)) & is.finite(beta)
	if (!any(vary)) return(out)
	X <- X[ok, vary, drop = FALSE]
	e <- eta[ok]
	# expected information of one probit observation at predictor e; a symmetric tie is stored twice
	w <- exp(2 * stats::dnorm(e, log = TRUE) - stats::pnorm(e, log.p = TRUE) - stats::pnorm(-e, log.p = TRUE))
	if (symmetric) w <- w / 2
	sw <- sum(w)
	if (!is.finite(sw) || sw <= 0) return(out)
	Xc <- sweep(X, 2L, colSums(w * X) / sw)
	ev <- eigen(crossprod(Xc * sqrt(w)), symmetric = TRUE)
	if (!is.finite(ev$values[1L]) || ev$values[1L] <= 0) return(out)
	# variance of the product on each dyad under that information
	lam <- pmax(ev$values, ev$values[1L] * 1e-12)
	v <- drop((Xc %*% ev$vectors)^2 %*% (1 / lam))
	out[ok] <- ifelse(v <= sd_max^2, drop(Xc %*% beta[vary]), 0)
	out
}
