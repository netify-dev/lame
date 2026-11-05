# split-sample jackknife bias correction for the static ALS estimator. a fixed-effects binary fit
# stretches the latent scale by about 1 + B / m (m cells per actor); refits on two interleaved halves
# carry twice the stretch, so 2 * full - mean(half_1, half_2) cancels the first-order term

# validate a user-supplied bias_correct value (NULL = the family default)
.ae_check_bias_correct <- function(bias_correct) {
	if (is.null(bias_correct)) return(NULL)
	allowed <- c("jackknife", "bootstrap", "none", "eb")
	if (!is.character(bias_correct) || length(bias_correct) != 1L ||
			is.na(bias_correct) || !bias_correct %in% allowed) {
		cli::cli_abort(c(
			"{.arg bias_correct} must be NULL or one of {.val {allowed}}.",
			"i" = "NULL uses the family and estimator defaults; none requests the uncorrected fit, jackknife the split-sample correction, bootstrap the parametric-bootstrap correction and eb the empirical-Bayes estimator of a static binary fit."
		))
	}
	bias_correct
}

# default correction of a fixed-effects fit: the jackknife for binary IRLS fits, none otherwise; where
# the jackknife cannot be computed, the size rule of .ae_size_state chooses between the uncorrected
# fit and the parametric bootstrap (ame_als_boot_bias.R)
.ae_resolve_bias_correct <- function(
	bias_correct,
	family,
	non_normal_method,
	n_time
) {
	if (!is.null(bias_correct)) return(bias_correct)
	if (identical(family, "binary") && identical(non_normal_method, "irls")) {
		"jackknife"
	} else {
		"none"
	}
}

# size of a binary panel for the fixed-effects stretch: mean observed cells per actor on the sparser
# side, times the share of the rarer outcome; `boundary` lists actors with no tie or no non-tie
.ae_binary_size <- function(
	Yarr,
	symmetric,
	bip = FALSE,
	row_names = NULL,
	col_names = NULL
) {
	obs <- !is.na(Yarr)
	pos <- obs & Yarr > 0
	pos[is.na(pos)] <- FALSE
	share <- if (any(obs)) sum(pos) / sum(obs) else 0
	side <- function(margin, nm, role) {
		n <- apply(obs, margin, sum)
		k <- apply(pos, margin, sum)
		act <- n > 0
		nm <- nm %||% paste0(if (margin == 1L) "row" else "col", seq_along(n))
		bnd <- act & (k == 0 | k == n)
		list(
			m = if (any(act)) mean(n[act]) else 0,
			boundary = if (any(bnd)) paste0(nm[bnd], role) else character(0)
		)
	}
	if (symmetric) {
		r <- side(1L, row_names, "")
		m <- r$m
		boundary <- r$boundary
	} else {
		r <- side(1L, row_names, if (bip) "" else " (sender)")
		cc <- side(2L, col_names, if (bip) "" else " (receiver)")
		m <- min(r$m, cc$m)
		boundary <- c(r$boundary, cc$boundary)
	}
	list(
		size = m * min(share, 1 - share),
		cells = m,
		share = share,
		boundary = boundary
	)
}

# relative bias allowed on the dyadic coefficients and on va and vb; the size floors follow from it
.ae_bias_tol <- c(slope = 0.02, va = 0.05, vb = 0.05)

# bias of the uncorrected binary ALS estimate, modelled as c / size: envelope constants `slope` and
# `vc` per rank band (rows: R = 0, 1-2, 3-4, 5-8) and additive-effect spread band (columns, upper
# edges `spread`), `kmax` the largest size measured, and `jk_floor0` the size from which the R = 0
# split-sample correction holds; a band whose floor exceeds `kmax` is corrected by the parametric bootstrap
.ae_size_basis <- list(
	jk_floor0 = 25,
	max_rank = c(0L, 2L, 4L, 8L),
	spread = c(0.7, 1.25, 2.25),
	slope = rbind(
		c(1.25, 1.25, 1.62),
		c(1.48, 1.87, 2.94),
		c(1.83, 2.48, 4.64),
		c(1.83, 2.90, 5.56)
	),
	vc = rbind(
		c(3.63, 3.63, 4.74),
		c(3.87, 4.40, 7.72),
		c(5.85, 5.85, 7.06),
		c(5.85, 5.85, 10.92)
	),
	kmax = rbind(
		c(191, 198, 174),
		c(193, 200, 174),
		c(189, 162, 175),
		c(194, 163, 176)
	)
)

# fitted additive-effect spread of a binary fit in latent noise sd units (unestimated effects left out)
.ae_fit_spread <- function(object) {
	sig <- if (identical(object$link, "logit")) pi / sqrt(3) else 1
	sdv <- function(x) {
		x <- x[is.finite(x)]
		if (length(x) > 1L) stats::sd(x) else NA_real_
	}
	s <- c(sdv(object$a), if (!isTRUE(object$symmetric)) sdv(object$b))
	if (!any(is.finite(s))) return(NA_real_)
	max(s, na.rm = TRUE) / sig
}

# constants and floor for rank R and fitted spread s; the floor is Inf beyond the measured ranks,
# spreads and sizes
.ae_size_band <- function(R, spread) {
	b <- .ae_size_basis
	tol <- .ae_bias_tol
	i <- which(R <= b$max_rank)[1L]
	j <- if (is.finite(spread)) which(spread <= b$spread)[1L] else 1L
	lo <- if (is.na(i) || i == 1L) 0L else b$max_rank[i - 1L] + 1L
	rank <- if (is.na(i)) {
		sprintf("R = %d", R)
	} else if (lo == b$max_rank[i]) {
		sprintf("R = %d", lo)
	} else sprintf("R = %d to %d", lo, b$max_rank[i])
	cs <- if (is.na(i) || is.na(j)) NA_real_ else b$slope[i, j]
	cv <- if (is.na(i) || is.na(j)) NA_real_ else b$vc[i, j]
	km <- if (is.na(i) || is.na(j)) NA_real_ else b$kmax[i, j]
	floor <- if (is.finite(cs) && is.finite(cv)) {
		ceiling(max(cs / tol[["slope"]], cv / min(tol[["va"]], tol[["vb"]])))
	} else Inf
	reach <- is.finite(floor) && is.finite(km) && floor <= km
	# largest spread band of this rank whose floor is measured
	ok_j <- if (is.na(i)) integer(0) else which(vapply(seq_along(b$spread), function(jj) {
		f <- max(b$slope[i, jj] / tol[["slope"]], b$vc[i, jj] / min(tol[["va"]], tol[["vb"]]))
		is.finite(f) && is.finite(b$kmax[i, jj]) && f <= b$kmax[i, jj]
	}, logical(1)))
	list(
		floor = if (reach) floor else Inf,
		slope = cs,
		vc = cv,
		rank = rank,
		spread_edge = if (is.na(j)) NA_real_ else b$spread[j],
		spread_max = if (length(ok_j)) b$spread[max(ok_j)] else NA_real_,
		rank_ok = !is.na(i),
		spread_ok = reach,
		kmax = km
	)
}

# size floor of an uncorrected binary ALS fit of rank R and additive-effect spread s (lowest band by default)
.ae_size_floor <- function(R, spread = 0) .ae_size_band(as.integer(R), spread)$floor

# size rule: the fit's size and spread, the floor for its band, the bias at its size and whether
# it is supported
.ae_size_state <- function(object) {
	spread <- .ae_fit_spread(object)
	band <- .ae_size_band(as.integer(object$R), spread)
	sz <- .ae_binary_size(
		object$Y,
		isTRUE(object$symmetric),
		identical(object$mode, "bipartite"),
		object$row_names,
		object$col_names
	)
	k <- sz$size
	pct <- function(cc) if (k > 0 && is.finite(cc)) 100 * cc / k else Inf
	list(
		size = k,
		floor = band$floor,
		spread = spread,
		band = band,
		boundary = sz$boundary,
		supported = is.finite(band$floor) && k >= band$floor,
		slope_pct = pct(band$slope),
		vc_pct = pct(band$vc),
		rank = band$rank
	)
}

# basis text of an uncorrected binary fit at or above the floor
.ae_size_basis_text <- function(st) {
	sprintf(
		"%.0f rarer-outcome cells per actor and a fitted additive-effect spread of %.2f, at or above the floor of %.0f for %s and spread up to %.2f; measured bias of the uncorrected estimate at this size up to about %.1f%% on the dyadic coefficients and %.1f%% on va and vb",
		st$size,
		st$spread,
		st$floor,
		st$rank,
		st$band$spread_edge,
		st$slope_pct,
		st$vc_pct
	)
}

# refusal of a binary fit the correction does not reach and the size rule does not support
.ae_size_refuse <- function(st, problem) {
	b <- st$band
	tol <- .ae_bias_tol
	tol_txt <- sprintf(
		"%.0f%% on the coefficients and %.0f%% on va / vb",
		100 * tol[["slope"]],
		100 * min(tol[["va"]], tol[["vb"]])
	)
	if (!b$rank_ok) {
		head <- sprintf("the size floor is not measured for %s", st$rank)
		why <- sprintf(
			"The finite-sample bias of the uncorrected estimate grows with the rank and was measured up to R = %d.",
			max(.ae_size_basis$max_rank)
		)
	} else if (!b$spread_ok) {
		head <- sprintf(
			"its fitted additive-effect spread is %.2f (sd of the estimated sender / receiver effects over the latent noise sd)",
			st$spread
		)
		why <- if (is.na(b$spread_edge) || !is.finite(b$slope)) {
			sprintf(
				"The bias of the uncorrected estimate grows with the spread of the effects; at %s it was measured to stay within %s from a size floor for spreads up to %.2f, and spreads above that were not measured.",
				st$rank,
				tol_txt,
				b$spread_max
			)
		} else {
			sprintf(
				"The bias of the uncorrected estimate grows with the spread of the effects; at %s it was measured to stay within %s from a size floor for spreads up to %.2f, and at spreads up to %.2f no measured size (up to %.0f rarer-outcome cells per actor) kept it there.",
				st$rank,
				tol_txt,
				b$spread_max,
				b$spread_edge,
				b$kmax
			)
		}
	} else {
		head <- sprintf("%.1f rarer-outcome cells per actor", st$size)
		why <- sprintf(
			"At this size and a fitted additive-effect spread of %.2f the measured bias curve puts the uncorrected estimate about %.0f%% high on the dyadic coefficients and %.0f%% on va / vb; binary ALS is supported from %.0f rarer-outcome cells per actor for %s and spreads up to %.2f, where the bias is within %s.",
			st$spread,
			st$slope_pct,
			st$vc_pct,
			st$floor,
			st$rank,
			b$spread_edge,
			tol_txt
		)
	}
	msg <- c(
		"Binary ALS is not supported for this network: {head}.",
		"x" = "{why}"
	)
	if (!is.null(problem)) {
		msg <- c(msg, "i" = "The split-sample correction cannot be computed: {problem}.")
	}
	cli::cli_abort(
		c(
			msg,
			"i" = "Use {.fn ame} / {.fn lame} (MCMC), which integrate the effects out and are calibrated at this size, or {.code bias_correct = \"none\"} for the uncorrected ALS estimate."
		),
		class = "lame_als_size"
	)
}

# fixed seed of the random split, so a fit is reproducible and independent of the caller's RNG stream
.ae_jk_seed <- 6886L

# the two half samples of an [n_row, n_col, T] panel, as the cells each half hides:
#   checkerboard: (i + j + t) mod 2, so every dyad stays observed in both halves (panels only)
#   random: each dyad-period by a fair coin from a fixed internal seed
#   periods: odd periods in one half, even in the other
# one-mode splits keep both directions of a dyad together, so an undirected panel stays undirected
.ae_jk_masks <- function(d, kind, bip) {
	Tt <- d[3L]
	if (identical(kind, "checkerboard")) {
		ij <- outer(seq_len(d[1L]), seq_len(d[2L]), "+")
		h <- array(FALSE, dim = d)
		for (t in seq_len(Tt)) h[, , t] <- (ij + t) %% 2L == 0L
		return(list(h, !h))
	}
	if (identical(kind, "periods")) {
		h <- array(FALSE, dim = d)
		for (t in seq_len(Tt)) h[, , t] <- t %% 2L == 0L
		return(list(h, !h))
	}
	old_kind <- RNGkind()
	had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
	if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
	on.exit({
		suppressWarnings(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]))
		if (had_seed) {
			assign(".Random.seed", old_seed, envir = globalenv())
		} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			rm(".Random.seed", envir = globalenv())
		}
	}, add = TRUE)
	suppressWarnings(set.seed(
		.ae_jk_seed,
		kind = "Mersenne-Twister",
		normal.kind = "Inversion",
		sample.kind = "Rejection"
	))
	h <- array(stats::runif(prod(d)) < 0.5, dim = d)
	if (!bip && d[1L] == d[2L]) {
		for (t in seq_len(Tt)) {
			m <- h[, , t]
			m[lower.tri(m)] <- t(m)[lower.tri(m)]
			h[, , t] <- m
		}
	}
	list(h, !h)
}

# per-actor observed and positive cell counts (rows senders, or all cells when symmetric; columns
# receivers); `no_tie` / `no_nontie` flag actors at one boundary value, whose effect hits the link cap
.ae_jk_actor_state <- function(Y, family, symmetric) {
	obs <- !is.na(Y)
	pos <- obs & Y > 0
	pos[is.na(pos)] <- FALSE
	side <- function(margin) {
		n <- apply(obs, margin, sum)
		p <- apply(pos, margin, sum)
		fam_bound <- family %in% c("binary", "poisson")
		no_tie <- fam_bound & n > 0 & p == 0
		no_nontie <- identical(family, "binary") & n > 0 & p == n
		list(
			n = n,
			no_tie = no_tie,
			no_nontie = no_nontie,
			bound = no_tie | no_nontie
		)
	}
	row <- side(1)
	col <- if (symmetric) row else side(2)
	list(row = row, col = col)
}

# connected components of the observed-dyad graph among actors with observed cells
.ae_jk_active_components <- function(Y, symmetric, bip) {
	obs <- apply(!is.na(Y), c(1, 2), sum)
	if (symmetric) obs <- obs + t(obs)
	comp <- .ae_components(obs, bip)
	if (bip) {
		lab <- c(comp$row[rowSums(obs) > 0], comp$col[colSums(obs) > 0])
	} else {
		lab <- comp$row[rowSums(obs) + colSums(obs) > 0]
	}
	length(unique(lab))
}

# half samples for the correction, or the reason there are none. the first admissible split
# (checkerboard, random, periods) keeps every actor estimable in both halves and neither half's graph
# disconnected; a panel with a multiplicative term needs five or more periods, and a single
# cross-section goes to .ae_jk_split_cs. unequal period halves are weighted as
# c0 full + c (half1 + half2), c = -1 / (1 / f1 + 1 / f2 - 2), c0 = 1 - 2 c, for observed shares f1, f2
.ae_jk_split <- function(Yarr, family, symmetric, bip, R = 0L) {
	d <- dim(Yarr)
	full <- .ae_jk_actor_state(Yarr, family, symmetric)
	names_row <- dimnames(Yarr)[[1L]] %||% paste0("row", seq_len(d[1L]))
	names_col <- dimnames(Yarr)[[2L]] %||% paste0("col", seq_len(d[2L]))
	sides <- if (symmetric) "row" else c("row", "col")
	nm_of <- function(side) if (side == "row") names_row else names_col
	role <- function(side, x) {
		if (symmetric || bip || !length(x)) return(x)
		paste0(x, if (side == "row") " (sender)" else " (receiver)")
	}
	few <- function(x) paste(utils::head(x, 3L), collapse = ", ")
	hint_actor <- "The correction needs every actor to keep observed cells with both outcome values in each half sample; sparse or small panels may not."
	active_row <- full$row$n > 0
	active_col <- full$col$n > 0
	base <- list(
		common_row = active_row,
		common_col = active_col,
		hint = hint_actor
	)

	if (d[3L] < 2L && R > 0L) {
		return(c(
			base,
			list(
				problem = sprintf(
					"a single cross-section with a multiplicative term (R = %d) cannot carry the correction: neither a checkerboard nor a random split of its cells corrects the fit without over-correcting it",
					R
				),
				hint = "The correction needs a panel with five or more periods when R > 0.",
				tried = NULL
			)
		))
	}
	if (R > 0L && d[3L] < 5L) {
		return(c(
			base,
			list(
				problem = sprintf(
					"with a multiplicative term (R = %d) and %d periods each dyad keeps one or two periods in a half sample, too few for its product cell, so the correction over-corrects",
					R,
					d[3L]
				),
				hint = "The correction needs five or more periods when R > 0, and two or more when R = 0.",
				tried = NULL
			)
		))
	}

	# at R = 0 the combination over-corrects a small panel, where the second-order term grows
	if (R == 0L && d[3L] >= 2L) {
		size <- .ae_binary_size(Yarr, symmetric, bip)$size
		jk_floor <- .ae_size_basis$jk_floor0
		if (size < jk_floor) {
			return(c(
				base,
				list(
					problem = sprintf(
						"a panel of %.1f rarer-outcome cells per actor is below the %.0f from which the correction holds at R = 0",
						size,
						jk_floor
					),
					hint = sprintf(
						"The correction needs %.0f or more rarer-outcome cells per actor at R = 0.",
						jk_floor
					),
					tried = NULL
				)
			))
		}
	}

	# an actor at a boundary value in the full panel: no split can help
	at_full <- unlist(lapply(sides, function(sd) role(sd, nm_of(sd)[full[[sd]]$bound])))
	if (length(at_full)) {
		n_nt <- sum(vapply(sides, function(sd) sum(full[[sd]]$no_tie), numeric(1)))
		n_nn <- sum(vapply(sides, function(sd) sum(full[[sd]]$no_nontie), numeric(1)))
		what <- if (identical(family, "poisson")) {
			sprintf("%d with only zero counts", n_nt)
		} else {
			sprintf("%d without a tie, %d without a non-tie", n_nt, n_nn)
		}
		return(c(
			base,
			list(
				problem = sprintf(
					"%d actor%s %s a single outcome value in the full panel (%s; %s), so %s effect sits at the link cap in every fit",
					length(at_full),
					if (length(at_full) > 1L) "s" else "",
					if (length(at_full) > 1L) "have" else "has",
					what,
					few(at_full),
					if (length(at_full) > 1L) "their" else "its"
				),
				tried = NULL
			)
		))
	}

	if (d[3L] < 2L) {
		return(.ae_jk_split_cs(
			Yarr,
			family,
			symmetric,
			bip,
			base,
			full,
			sides,
			nm_of,
			role,
			few
		))
	}

	n_comp <- .ae_jk_active_components(Yarr, symmetric, bip)
	kinds <- c("checkerboard", "random", "periods")
	tried <- vector("list", length(kinds))
	for (k in seq_along(kinds)) {
		hides <- .ae_jk_masks(d, kinds[k], bip)
		lost <- no_tie <- no_non <- character(0)
		graph <- FALSE
		for (h in hides) {
			Yh <- Yarr
			Yh[h] <- NA_real_
			st <- .ae_jk_actor_state(Yh, family, symmetric)
			for (sd in sides) {
				f <- full[[sd]]; hs <- st[[sd]]; nm <- nm_of(sd)
				ok <- f$n > 0
				lost <- c(lost, role(sd, nm[ok & hs$n == 0]))
				no_tie <- c(no_tie, role(sd, nm[ok & hs$no_tie]))
				no_non <- c(no_non, role(sd, nm[ok & hs$no_nontie]))
			}
			if (.ae_jk_active_components(Yh, symmetric, bip) > n_comp) graph <- TRUE
		}
		lost <- unique(lost); no_tie <- unique(no_tie); no_non <- unique(no_non)
		tried[[k]] <- data.frame(
			split = kinds[k],
			no_cell = length(lost),
			no_tie = length(no_tie),
			no_nontie = length(no_non),
			graph_split = graph,
			examples = few(unique(c(lost, no_tie, no_non))),
			stringsAsFactors = FALSE
		)
		if (!length(lost) && !length(no_tie) && !length(no_non) && !graph) {
			n_obs <- vapply(hides, function(h) sum(!h & !is.na(Yarr)), numeric(1))
			f <- n_obs / sum(!is.na(Yarr))
			c1 <- if (identical(kinds[k], "periods") && all(f > 0)) {
				-1 / (1 / f[1L] + 1 / f[2L] - 2)
			} else -0.5
			return(c(
				base,
				list(
					problem = NULL,
					split = kinds[k],
					hide = hides,
					weights = c(full = 1 - 2 * c1, half = c1),
					tried = do.call(rbind, tried[seq_len(k)])
				)
			))
		}
	}
	tried <- do.call(rbind, tried)
	desc <- vapply(seq_len(nrow(tried)), function(k) {
		r <- tried[k, ]
		parts <- c(
			if (r$no_cell > 0) sprintf("%d without a cell", r$no_cell),
			if (r$no_tie > 0) sprintf("%d without a tie", r$no_tie),
			if (r$no_nontie > 0) sprintf("%d without a non-tie", r$no_nontie)
		)
		txt <- if (length(parts)) {
			sprintf("%s split: %s in a half (%s)", r$split, paste(parts, collapse = ", "), r$examples)
		} else sprintf("%s split", r$split)
		if (r$graph_split) txt <- paste0(txt, if (length(parts)) "; " else ": ", "a half's observed graph falls apart")
		txt
	}, character(1))
	c(
		base,
		list(
			problem = paste0(
				"no admissible split keeps every actor estimable in both half samples (",
				paste(desc, collapse = "; "),
				")"
			),
			tried = tried
		)
	)
}

# checkerboard over a fixed, outcome-independent balanced permutation of the actors (attempt-th draw
# from an internal seed; the caller's RNG state is restored)
.ae_jk_permuted_checkerboard <- function(d, bip, attempt) {
	old_kind <- RNGkind()
	had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
	if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
	on.exit({
		suppressWarnings(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]))
		if (had_seed) {
			assign(".Random.seed", old_seed, envir = globalenv())
		} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			rm(".Random.seed", envir = globalenv())
		}
	}, add = TRUE)
	suppressWarnings(set.seed(
		629401L, kind = "Mersenne-Twister",
		normal.kind = "Inversion", sample.kind = "Rejection"
	))
	for (k in seq_len(attempt)) {
		rows <- sample.int(d[1L])
		cols <- if (bip) sample.int(d[2L]) else rows
	}
	h <- array(FALSE, d)
	for (t in seq_len(d[3L])) h[, , t] <- (outer(rows, cols, "+") + t) %% 2L == 0L
	list(h, !h)
}

# half samples of a single cross-section at R = 0, or the reason there are none: a checkerboard of
# (i + j) parity, whose halves split into two components with unidentified relative levels (handled by
# .ae_jk_anchor_cs); admissible when every actor keeps a cell, a tie and a non-tie in each half, and
# used from .ae_size_basis$jk_floor0 rarer-outcome cells per actor (see ?ame_als)
.ae_jk_split_cs <- function(
	Yarr,
	family,
	symmetric,
	bip,
	base,
	full,
	sides,
	nm_of,
	role,
	few
) {
	d <- dim(Yarr)
	size <- .ae_binary_size(Yarr, symmetric, bip)$size
	floor_cs <- .ae_size_basis$jk_floor0
	if (size < floor_cs) {
		return(c(
			base,
			list(
				problem = sprintf(
					"a single cross-section of %.1f rarer-outcome cells per actor is below the %.0f from which its checkerboard correction holds",
					size,
					floor_cs
				),
				hint = sprintf(
					"The correction of a single cross-section needs %.0f or more rarer-outcome cells per actor and R = 0.",
					floor_cs
				),
				tried = NULL
			)
		))
	}
	ncomp <- function(Y) {
		obs <- apply(!is.na(Y), c(1, 2), sum)
		if (symmetric) {
			obs <- obs + t(obs)
			comp <- .ae_components(obs, FALSE)
			return(length(unique(comp$row[rowSums(obs) > 0])))
		}
		comp <- .ae_components(obs, TRUE)
		length(unique(c(comp$row[rowSums(obs) > 0], comp$col[colSums(obs) > 0])))
	}
	c_full <- ncomp(Yarr)
	trials <- vector("list", 32L)
	for (attempt in seq_len(32L)) {
		hides <- if (attempt == 1L) .ae_jk_masks(d, "checkerboard", bip) else
			.ae_jk_permuted_checkerboard(d, bip, attempt - 1L)
		lost <- no_tie <- no_non <- character(0)
		graph <- FALSE
		for (h in hides) {
			Yh <- Yarr
			Yh[h] <- NA_real_
			st <- .ae_jk_actor_state(Yh, family, symmetric)
			for (sd in sides) {
				f <- full[[sd]]; hs <- st[[sd]]; nm <- nm_of(sd)
				ok <- f$n > 0
				lost <- c(lost, role(sd, nm[ok & hs$n == 0]))
				no_tie <- c(no_tie, role(sd, nm[ok & hs$no_tie]))
				no_non <- c(no_non, role(sd, nm[ok & hs$no_nontie]))
			}
			if (ncomp(Yh) > 2L * c_full) graph <- TRUE
		}
		lost <- unique(lost); no_tie <- unique(no_tie); no_non <- unique(no_non)
		trials[[attempt]] <- data.frame(
			split = "checkerboard",
			attempt = attempt,
			no_cell = length(lost),
			no_tie = length(no_tie),
			no_nontie = length(no_non),
			graph_split = graph,
			examples = few(unique(c(lost, no_tie, no_non))),
			stringsAsFactors = FALSE
		)
		if (!length(lost) && !length(no_tie) && !length(no_non) && !graph) {
			return(c(
				base,
				list(
					problem = NULL,
					split = "checkerboard",
					cross_section = TRUE,
					hide = hides,
					checkerboard_attempt = attempt,
					weights = c(full = 2, half = -0.5),
					tried = do.call(rbind, trials[seq_len(attempt)])
				)
			))
		}
	}
	tried <- do.call(rbind, trials)
	best <- which.min(tried$no_cell + tried$no_tie + tried$no_nontie + tried$graph_split)
	z <- tried[best, ]
	parts <- c(
		if (z$no_cell) sprintf("%d without a cell", z$no_cell),
		if (z$no_tie) sprintf("%d without a tie", z$no_tie),
		if (z$no_nontie) sprintf("%d without a non-tie", z$no_nontie)
	)
	txt <- if (length(parts)) {
		sprintf(
			"all 32 balanced checkerboard partitions fail; the best leaves %s in a half (%s)",
			paste(parts, collapse = ", "),
			z$examples
		)
	} else "all 32 balanced checkerboard partitions fail"
	if (z$graph_split) txt <- paste0(txt, "; a half's observed graph falls apart")
	c(base, list(problem = txt, tried = tried))
}

# raw additive effects of an assembled fit with the node-covariate parts folded back in, so that
# mu + x'beta + row_i + col_j + o_ij is its linear predictor
.ae_jk_raw_effects <- function(object, p_dyad) {
	beta <- as.numeric(object$beta)
	a <- as.numeric(object$a); b <- as.numeric(object$b)
	pr <- if (is.null(object$W_row)) 0L else ncol(object$W_row)
	pc <- if (is.null(object$W_col)) 0L else ncol(object$W_col)
	if (isTRUE(object$symmetric)) {
		W_node <- cbind(object$W_row, object$W_col)
		if (!is.null(W_node) && ncol(W_node) > 0L) {
			a <- a + as.numeric(W_node %*% beta[p_dyad + seq_len(pr + pc)])
		}
		b <- a
	} else {
		if (pr > 0L) {
			a <- a + as.numeric(object$W_row %*% beta[p_dyad + seq_len(pr)])
		}
		if (pc > 0L) {
			b <- b + as.numeric(object$W_col %*% beta[p_dyad + pr + seq_len(pc)])
		}
	}
	list(row = a, col = b)
}

# warm start of a half-sample refit: the full fit's state with node-covariate parts folded into the
# effects and the Poisson latent shift undone, so it reproduces the fitted linear predictor
.ae_jk_init <- function(object, p_dyad) {
	# a withheld product is restored to its variational means
	object <- .ae_product_restore(object)
	eff <- .ae_jk_raw_effects(object, p_dyad)
	list(
		mu = object$mu + .ae_default(object$pois_latent_var, 0) / 2,
		beta = as.numeric(object$beta)[seq_len(p_dyad)],
		a = eff$row,
		b = eff$col,
		U = object$U,
		V = object[["V"]],
		L = object$L
	)
}

# one half-sample refit of `object`'s estimator on `Yh` (and working response `Zh`), warm-started
# from `object`
.ae_jk_half <- function(object, Yh, Zh, max_iter, tol, init = NULL) {
	X <- object[["X"]]
	prep <- .ae_prep_from_canonical(
		Yh,
		X,
		object$W_row,
		object$W_col,
		object$mode,
		object$symmetric,
		object$row_names,
		object$col_names,
		object$time_names
	)
	if (is.null(init)) init <- .ae_jk_init(object, prep$p_dyad)
	lrm <- .ae_default(object$lowrank_method, "mm")
	nnm <- .ae_default(object$non_normal_method, "transform")
	lsv <- .ae_default(object$linear_solver, "eigen")
	sym <- isTRUE(object$symmetric)
	irls_iters <- NA_integer_
	if (!is.null(Zh)) {
		Z <- Zh
		if (sym) for (t in seq_len(prep$Tt)) Z[, , t] <- .ae_symmetrize(Z[, , t])
		fit <- .ae_bcd_fit(
			Z,
			X,
			object$R,
			sym,
			max_iter,
			tol,
			init = init,
			lowrank_method = lrm,
			linear_solver = lsv,
			certify = TRUE,
			unit_rules = identical(object$family, "normal")
		)
		out <- .ae_assemble(
			fit,
			prep,
			object$family,
			object$mode,
			sym,
			object$R,
			object$longitudinal,
			call = NULL,
			Zwork = Z
		)
	} else if (identical(nnm, "irls")) {
		ir <- .ae_irls_run(
			prep,
			object$family,
			object$link,
			object$R,
			sym,
			max_iter,
			tol,
			lowrank_method = lrm,
			init = init,
			linear_solver = lsv,
			warm = TRUE
		)
		fit <- ir$fit
		irls_iters <- ir$irls_iters
		out <- .ae_assemble(
			fit,
			prep,
			object$family,
			object$mode,
			sym,
			object$R,
			object$longitudinal,
			call = NULL,
			Zwork = ir$Z,
			link = object$link,
			Wobs = ir$W
		)
		out$objective <- ir$best_obj
		out$irls_iterations <- ir$irls_iters
		out$irls_diagnostics <- ir$stationarity
	} else {
		Z <- .ae_working_response(Yh, object$family)
		if (sym) for (t in seq_len(prep$Tt)) Z[, , t] <- .ae_symmetrize(Z[, , t])
		fit <- .ae_bcd_fit(
			Z,
			X,
			object$R,
			sym,
			max_iter,
			tol,
			init = init,
			lowrank_method = lrm,
			linear_solver = lsv,
			certify = TRUE
		)
		out <- .ae_assemble(
			fit,
			prep,
			object$family,
			object$mode,
			sym,
			object$R,
			object$longitudinal,
			call = NULL
		)
	}
	out$irls_iters <- irls_iters
	out
}

# a converged half-sample refit, or an error of class "lame_jk_nonconvergence"; an unconverged half
# is continued once with five times the budget so the correction never combines a last iterate
.ae_jk_half_converged <- function(object, Yh, Zh, max_iter, tol, half) {
	budget <- as.integer(max_iter)
	init <- NULL
	attempts <- vector("list", 2L)
	for (attempt in 1:2) {
		h <- suppressWarnings(.ae_jk_half(object, Yh, Zh, budget, tol, init = init))
		attempts[[attempt]] <- list(
			max_iter = budget,
			converged = isTRUE(h$converged),
			iterations = h$iterations,
			irls_iterations = h$irls_iters
		)
		if (isTRUE(h$converged)) {
			h$attempts <- attempts[seq_len(attempt)]
			return(h)
		}
		init <- .ae_jk_init(h, dim(object[["X"]])[3L])
		budget <- 5L * as.integer(max_iter)
	}
	cli::cli_abort(
		c(
			"Half {half} of the split-sample bias correction did not converge within {.arg max_iter} = {max_iter} iterations, nor within {budget} from its last iterate.",
			"x" = "The correction combines the full fit with both halves, so a half that has not converged cannot be used.",
			"i" = "Increase {.arg max_iter}, or use {.code bias_correct = \"none\"} for the uncorrected fit."
		),
		class = "lame_jk_nonconvergence",
		half = half,
		attempts = attempts
	)
}

# double-centred multiplicative surface the fit's linear predictor carries (NULL at R = 0)
.ae_jk_product <- function(object) {
	# a withheld product is restored to its variational means
	object <- .ae_product_restore(object)
	if (object$R < 1L) return(NULL)
	O <- object$UVPM
	if (is.null(O)) {
		if (is.null(object$U)) return(NULL)
		U <- unname(object$U)
		return(if (isTRUE(object$symmetric)) {
			U %*% (object$L * t(U))
		} else {
			tcrossprod(U, unname(object[["V"]]))
		})
	}
	O <- unname(O)
	nr <- nrow(O); nc <- ncol(O)
	O - outer(rowMeans(O), rep(1, nc)) - outer(rep(1, nr), colMeans(O)) + mean(O)
}

# one fit's surface in the common gauge, double-centred over the actors with observed cells (Sr, Sc):
# intercept plus centring level `m`, dyadic coefficients, centred row and column effects and product
.ae_jk_anchor <- function(f, p_dyad, Sr, Sc) {
	eff <- .ae_jk_raw_effects(f, p_dyad)
	Fs <- outer(eff$row, eff$col, "+")
	O <- .ae_jk_product(f)
	if (!is.null(O)) Fs <- Fs + O
	g <- mean(Fs[Sr, Sc])
	rr <- rowMeans(Fs[, Sc, drop = FALSE]) - g
	cc <- colMeans(Fs[Sr, , drop = FALSE]) - g
	Q <- if (is.null(O)) NULL else Fs - outer(rr, cc, "+") - g
	list(
		m = unname(f$coefficients[1L]) + g,
		beta = as.numeric(f$beta)[seq_len(p_dyad)],
		rr = rr,
		cc = cc,
		Q = Q
	)
}

# the three anchored fits of a single cross-section: the full fit and twice the split surface, whose
# value at each cell comes from the half observing that cell's component block (levels across
# blocks are not identified and not used); each copy keeps its own half's dyadic coefficients
.ae_jk_anchor_cs <- function(object, halves, split, p_dyad) {
	Sr <- split$common_row
	Sc <- split$common_col
	surface <- function(f) {
		eff <- .ae_jk_raw_effects(f, p_dyad)
		Fs <- unname(f$coefficients[1L]) + outer(eff$row, eff$col, "+")
		O <- .ae_jk_product(f)
		if (!is.null(O)) Fs <- Fs + O
		Fs
	}
	own1 <- !split$hide[[1L]][, , 1L]
	Fh <- surface(halves[[1L]])
	Fh[!own1] <- surface(halves[[2L]])[!own1]
	g <- mean(Fh[Sr, Sc])
	rr <- rowMeans(Fh[, Sc, drop = FALSE]) - g
	cc <- colMeans(Fh[Sr, , drop = FALSE]) - g
	Q <- if (object$R > 0L) Fh - outer(rr, cc, "+") - g else NULL
	beta_of <- function(f) as.numeric(f$beta)[seq_len(p_dyad)]
	list(
		.ae_jk_anchor(object, p_dyad, Sr, Sc),
		list(m = g, beta = beta_of(halves[[1L]]), rr = rr, cc = cc, Q = Q),
		list(m = g, beta = beta_of(halves[[2L]]), rr = rr, cc = cc, Q = Q)
	)
}

# noise-corrected variances and covariance of a cross-section's split surface: the spread of its
# covariate-orthogonal effects less the full fit's sampling-noise term
.ae_jk_vc_cs <- function(object, anchored_split) {
	sym <- isTRUE(object$symmetric)
	bip <- identical(object$mode, "bipartite")
	obs_count <- apply(!is.na(object$Y), c(1, 2), sum)
	if (sym) obs_count <- obs_count + t(obs_count)
	active_row <- as.numeric(rowSums(obs_count) > 0)
	active_col <- as.numeric(colSums(obs_count) > 0)
	if (sym) {
		W_node <- object$W_row
		if (!is.null(object$W_col)) {
			W_node <- if (is.null(W_node)) object$W_col else cbind(W_node, object$W_col)
		}
		ea <- .ae_decompose(anchored_split$rr, W_node, active_row)$resid
		eb <- ea
	} else {
		ea <- .ae_decompose(anchored_split$rr, object$W_row, active_row)$resid
		eb <- .ae_decompose(anchored_split$cc, object$W_col, active_col)$resid
	}
	noise <- object$vc_noise
	net <- function(v, nm) {
		if (is.finite(v) && is.finite(noise[[nm]])) v - noise[[nm]] else v
	}
	va <- if (sum(active_row) > 1) stats::var(ea[active_row == 1]) else NA_real_
	vb <- if (sym) {
		NA_real_
	} else if (sum(active_col) > 1) {
		stats::var(eb[active_col == 1])
	} else NA_real_
	both <- active_row == 1 & active_col == 1
	cab <- if (sym || bip || sum(both) <= 2) NA_real_ else stats::cov(ea[both], eb[both])
	list(vc = c(va = net(va, "va"), vb = if (sym) vb else net(vb, "vb")), cab = cab)
}

# noise-corrected additive variances before the floor at zero
.ae_jk_unfloored_vc <- function(object) {
	one <- function(nm) {
		raw <- object$vc_raw[[nm]]
		noise <- object$vc_noise[[nm]]
		if (is.finite(raw) && is.finite(noise)) raw - noise else raw
	}
	c(va = one("va"), vb = one("vb"))
}

# split-sample jackknife of an assembled static ALS fit: the full and half fits are put in a common
# gauge (surfaces double-centred over the actors with observed cells) and combined as
# 2 * full - mean(halves). `Z` is the working response of a working-scale refit, NULL otherwise.
# `explicit` means the caller asked for the correction, so a split that cannot carry it is refused
# ("lame_jk_split"); by default the size rule (.ae_size_state) applies instead
.ae_jackknife <- function(object, max_iter, tol, Z = NULL, explicit = TRUE) {
	# a withheld product is restored to its variational means
	object <- .ae_product_restore(object)
	Yarr <- object$Y
	sym <- isTRUE(object$symmetric)
	bip <- identical(object$mode, "bipartite")
	# by default refuse a binary fit whose additive-effect spread is beyond the measured band
	if (!explicit && identical(object$family, "binary")) {
		spread <- .ae_fit_spread(object)
		if (is.finite(spread) && spread > max(.ae_size_basis$spread)) {
			st <- .ae_size_state(object)
			.ae_size_refuse(st, NULL)
		}
	}
	split <- .ae_jk_split(Yarr, object$family, sym, bip, object$R)
	if (!is.null(split$problem)) {
		problem <- split$problem
		hint <- split$hint
		if (explicit) {
			cli::cli_abort(
				c(
					"The split-sample bias correction cannot be computed for this panel.",
					"x" = "{problem}.",
					"i" = "{hint}",
					"i" = "Use {.code bias_correct = \"none\"} for the uncorrected fit, or {.fn ame} / {.fn lame} (MCMC)."
				),
				class = "lame_jk_split"
			)
		}
		# by default a binary network below the size floor for its rank and spread is corrected by the
		# parametric bootstrap (an unmeasured rank is refused); at or above the floor it is returned
		# uncorrected
		st <- .ae_size_state(object)
		if (!isTRUE(st$supported)) {
			if (!isTRUE(st$band$rank_ok)) .ae_size_refuse(st, problem)
			out <- .ae_boot_bias(object, max_iter, tol)
			out$bias_correct$size <- st$size
			out$bias_correct$floor <- st$floor
			out$bias_correct$spread <- st$spread
			out$bias_correct$correction_problem <- problem
			return(out)
		}
		object$bias_correct <- list(
			method = "none",
			applied = FALSE,
			basis = .ae_size_basis_text(st),
			size = st$size,
			floor = st$floor,
			spread = st$spread,
			correction_problem = problem,
			tried = split$tried
		)
		return(object)
	}

	halves <- lapply(1:2, function(k) {
		hide <- split$hide[[k]]
		Yh <- Yarr
		Yh[hide] <- NA_real_
		Zh <- NULL
		if (!is.null(Z)) {
			Zh <- Z
			Zh[hide] <- NA_real_
		}
		.ae_jk_half_converged(object, Yh, Zh, max_iter, tol, k)
	})
	fits <- list(object, halves[[1L]], halves[[2L]])
	w <- split$weights
	jk <- function(x) w[["full"]] * x[[1L]] + w[["half"]] * (x[[2L]] + x[[3L]])
	R <- object$R
	p_dyad <- dim(object[["X"]])[3L]
	Sr <- split$common_row
	Sc <- split$common_col

	# common gauge: each surface double-centred over the actors with observed cells (for a single
	# cross-section, each cell's half value from the half that observes it)
	cross_section <- isTRUE(split$cross_section)
	anchored <- if (cross_section) {
		.ae_jk_anchor_cs(object, halves, split, p_dyad)
	} else {
		lapply(fits, .ae_jk_anchor, p_dyad = p_dyad, Sr = Sr, Sc = Sc)
	}
	pick <- function(nm) lapply(anchored, `[[`, nm)
	# combine on the actors with observed cells; the others keep their full-fit value times the shrink
	# measured on the rest
	combine <- function(nm, S) {
		full <- anchored[[1L]][[nm]]
		comb <- jk(pick(nm))
		den <- sum(full[S]^2)
		kappa <- if (den > 0) sum(comb[S] * full[S]) / den else 1
		out <- full * kappa
		out[S] <- comb[S]
		list(value = out, kappa = kappa)
	}
	m <- jk(pick("m"))
	beta_dyad <- jk(pick("beta"))
	rr <- combine("rr", Sr)
	cc <- if (sym) rr else combine("cc", Sc)
	Fs <- outer(rr$value, cc$value, "+")
	kappa_Q <- NA_real_
	if (R > 0L && !is.null(anchored[[1L]]$Q)) {
		Qc <- combine("Q", outer(Sr, Sc, "&"))
		Fs <- Fs + Qc$value
		kappa_Q <- Qc$kappa
	}
	if (sym) Fs <- (Fs + t(Fs)) / 2

	# back to the package gauge: effects summing to zero over all actors and a fully double-centred product
	g <- mean(Fs)
	row_raw <- rowMeans(Fs) - g
	col_raw <- colMeans(Fs) - g
	mu <- m + g
	U <- object$U; V <- object[["V"]]; L <- object$L
	O_new <- NULL
	if (R > 0L && !is.null(object$U)) {
		O_lin <- Fs - outer(row_raw, col_raw, "+") - g
		if (sym) {
			eg <- eigen((O_lin + t(O_lin)) / 2, symmetric = TRUE)
			idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
			U <- eg$vectors[, idx, drop = FALSE]
			L <- eg$values[idx]
			V <- U
			O_new <- U %*% (L * t(U))
		} else {
			sv <- svd(O_lin)
			dd <- sqrt(pmax(sv$d[seq_len(R)], 0))
			U <- sweep(sv$u[, seq_len(R), drop = FALSE], 2, dd, "*")
			V <- sweep(sv$v[, seq_len(R), drop = FALSE], 2, dd, "*")
			O_new <- tcrossprod(U, V)
		}
		dimnames(U) <- list(object$row_names, paste0("dim", seq_len(R)))
		dimnames(V) <- list(object$col_names, paste0("dim", seq_len(R)))
	}

	# node-covariate coefficients and covariate-orthogonal effects, by the estimator's own decomposition
	obs_count <- apply(!is.na(Yarr), c(1, 2), sum)
	if (sym) obs_count <- obs_count + t(obs_count)
	active_row <- as.numeric(rowSums(obs_count) > 0)
	active_col <- as.numeric(colSums(obs_count) > 0)
	if (sym) {
		W_node <- object$W_row
		if (!is.null(object$W_col)) {
			W_node <- if (is.null(W_node)) object$W_col else cbind(W_node, object$W_col)
		}
		dec <- .ae_decompose(row_raw, W_node, active_row)
		a <- dec$resid; b <- a
		mu <- mu + 2 * dec$intercept
		beta <- c(beta_dyad, dec$beta)
		row_part <- col_part <- row_raw - dec$intercept
	} else {
		dec_r <- .ae_decompose(row_raw, object$W_row, active_row)
		dec_c <- .ae_decompose(col_raw, object$W_col, active_col)
		a <- dec_r$resid; b <- dec_c$resid
		mu <- mu + dec_r$intercept + dec_c$intercept
		beta <- c(beta_dyad, dec_r$beta, dec_c$beta)
		row_part <- row_raw - dec_r$intercept
		col_part <- col_raw - dec_c$intercept
	}
	names(beta) <- names(object$beta)
	names(a) <- names(object$a); names(b) <- names(object$b)
	coefs <- c(intercept = mu, beta)

	# linear predictor, fitted values and residuals from the corrected parts
	Tt <- length(object$EZ)
	lat <- .ae_default(object$pois_latent_var, 0)
	add <- outer(row_part, col_part, "+")
	if (!is.null(O_new)) add <- add + O_new
	EZ <- fitted_list <- resid_list <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		ez <- mu + add
		if (p_dyad > 0L) for (k in seq_len(p_dyad)) {
			ez <- ez + beta_dyad[k] * object[["X"]][, , k, t]
		}
		ez <- unname(ez)
		if (!bip) diag(ez) <- NA_real_
		# the full fit's working response is its residual plus its linear predictor (marginal scale for
		# Poisson)
		resid_list[[t]] <- object$residuals[[t]] + (unname(object$EZ[[t]]) - ez)
		fit_t <- .ae_link_inverse(ez, object$family, object$link, latent_var = lat)
		if (!bip) diag(fit_t) <- NA_real_
		dimnames(ez) <- dimnames(object$EZ[[t]])
		dimnames(fit_t) <- dimnames(object$fitted[[t]])
		EZ[[t]] <- ez
		fitted_list[[t]] <- fit_t
	}
	names(EZ) <- names(object$EZ)
	names(fitted_list) <- names(object$fitted)
	names(resid_list) <- names(object$residuals)

	# variance components: the noise-corrected variances are combined before the floor at zero, then
	# floored once
	if (cross_section) {
		vs <- .ae_jk_vc_cs(object, anchored[[2L]])
		vcu <- jk(list(.ae_jk_unfloored_vc(object), vs$vc, vs$vc))
		cab <- jk(list(object$VC[["cab"]], vs$cab, vs$cab))
	} else {
		vcu <- jk(lapply(fits, .ae_jk_unfloored_vc))
		cab <- jk(lapply(fits, function(f) f$VC[["cab"]]))
	}
	VC <- object$VC
	VC[["va"]] <- if (is.finite(vcu[["va"]])) max(vcu[["va"]], 0) else NA_real_
	VC[["vb"]] <- if (is.finite(vcu[["vb"]])) max(vcu[["vb"]], 0) else NA_real_
	VC[["cab"]] <- cab
	vc_noise <- if (cross_section) object$vc_noise else jk(lapply(fits, `[[`, "vc_noise"))

	raw <- list(
		coefficients = object$coefficients,
		a = object$a,
		b = object$b,
		VC = object$VC,
		O = NULL
	)
	if (R > 0L) {
		raw$O <- .ae_jk_product(object)
		if (!is.null(raw$O)) dimnames(raw$O) <- list(object$row_names, object$col_names)
	}

	out <- object
	out$coefficients <- coefs
	out$mu <- unname(mu)
	out$beta <- beta
	out$a <- a
	out$b <- b
	out$U <- U; out[["V"]] <- V; out$L <- L
	if (!is.null(O_new)) {
		Om <- O_new
		dimnames(Om) <- list(object$row_names, object$col_names)
		out$UVPM <- Om
		if (sym) out$ULUPM <- Om
	}
	out$VC <- VC
	out$vc_noise <- vc_noise
	out$vc_raw <- c(
		va = unname(vcu[["va"]] + vc_noise[["va"]]),
		vb = unname(vcu[["vb"]] + vc_noise[["vb"]])
	)
	out$EZ <- EZ
	out$fitted <- fitted_list
	out$YPM <- fitted_list
	out$residuals <- resid_list
	# deviance of the corrected fit from its rebuilt fitted values; the full fit's stays in raw
	raw$deviance <- object$deviance
	out$deviance <- .ae_fit_deviance(
		object$Y,
		fitted_list,
		resid_list,
		object$family,
		object$link,
		eta = EZ,
		latent_var = lat
	)
	out$raw <- raw
	h1 <- halves[[1L]]; h2 <- halves[[2L]]
	half_conv <- c(half1 = isTRUE(h1$converged), half2 = isTRUE(h2$converged))
	# the reported estimate is converged when all three fits are
	out$converged <- isTRUE(object$converged) && all(half_conv)
	out$bias_correct <- list(
		method = "jackknife",
		applied = TRUE,
		split = split$split,
		weights = w,
		tried = split$tried,
		converged = half_conv,
		irls_diagnostics = lapply(fits[-1L], `[[`, "irls_diagnostics"),
		iterations = c(half1 = h1$iterations, half2 = h2$iterations),
		irls_iterations = c(half1 = h1$irls_iters, half2 = h2$irls_iters),
		shrink = c(a = rr$kappa, b = cc$kappa, UV = kappa_Q)
	)
	.ae_update_response_rho(out, out$link)
}
