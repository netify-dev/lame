# parametric-bootstrap bias correction of the static binary ALS estimator
# simulates the incidental-parameter bias twice and takes the secant step toward the fixed point
# theta + bias(theta) = estimate: corrected = estimate - b0 / (1 + s), s = 1 - b1 / b0 in [0, 1]

# replicates per round and their seed (the caller's RNG state is left alone)
.ae_bb_reps <- 30L
.ae_bb_seed <- 7741L

# share of a round's replicates that must refit and converge
.ae_bb_min_valid <- 0.8

# linear predictor of a parameter set, one matrix per period
.ae_bb_eta <- function(object, coef, a, b, O) {
	X <- object[["X"]]
	p_dyad <- dim(X)[3L]
	beta <- as.numeric(coef[-1L])
	pr <- if (is.null(object$W_row)) 0L else ncol(object$W_row)
	pc <- if (is.null(object$W_col)) 0L else ncol(object$W_col)
	a <- .ae_na_zero(as.numeric(a)); b <- .ae_na_zero(as.numeric(b))
	if (isTRUE(object$symmetric)) {
		W_node <- cbind(object$W_row, object$W_col)
		if (!is.null(W_node) && ncol(W_node) > 0L) {
			a <- a + as.numeric(W_node %*% beta[p_dyad + seq_len(pr + pc)])
		}
		b <- a
	} else {
		if (pr > 0L) a <- a + as.numeric(object$W_row %*% beta[p_dyad + seq_len(pr)])
		if (pc > 0L) b <- b + as.numeric(object$W_col %*% beta[p_dyad + pr + seq_len(pc)])
	}
	base <- unname(coef[1L]) + outer(a, b, "+")
	if (!is.null(O)) base <- base + O
	lapply(seq_len(dim(object$Y)[3L]), function(t) {
		ez <- base
		if (p_dyad > 0L) for (k in seq_len(p_dyad)) ez <- ez + beta[k] * X[, , k, t]
		ez
	})
}

# centred additive effects drawn from (va, vb, cab), as in .ae_boot_draw_ab
.ae_bb_draw_effects <- function(object, va, vb, cab) {
	nr <- object$dims$n_row; nc <- object$dims$n_col
	va <- max(.ae_default(va, 0), 0, na.rm = TRUE)
	vb <- max(.ae_default(vb, 0), 0, na.rm = TRUE)
	cab <- if (is.finite(cab)) cab else 0
	ctr <- function(x) if (length(x) > 1L) x - mean(x) else x
	if (isTRUE(object$symmetric)) {
		a <- ctr(stats::rnorm(nr, 0, sqrt(va)))
		return(list(a = a, b = a))
	}
	if (identical(object$mode, "bipartite")) {
		return(list(
			a = ctr(stats::rnorm(nr, 0, sqrt(va))),
			b = ctr(stats::rnorm(nc, 0, sqrt(vb)))
		))
	}
	S <- matrix(c(va, cab, cab, vb), 2, 2)
	eg <- eigen(S, symmetric = TRUE)
	root <- eg$vectors %*% (sqrt(pmax(eg$values, 0)) * t(eg$vectors))
	ab <- matrix(stats::rnorm(2 * nr), nr, 2) %*% root
	list(a = ctr(ab[, 1L]), b = ctr(ab[, 2L]))
}

# one simulated outcome panel with the fit's observed-cell pattern
.ae_bb_draw_y <- function(object, EZ, rho) {
	Y <- array(NA_real_, dim(object$Y))
	for (t in seq_along(EZ)) {
		ez <- EZ[[t]]
		ez[!is.finite(ez)] <- 0
		yt <- as.matrix(.ae_boot_draw_response(object, ez, rho, 1))
		storage.mode(yt) <- "double"
		if (isTRUE(object$symmetric)) {
			yt[lower.tri(yt)] <- 0
			yt <- yt + t(yt)
		}
		yt[is.na(object$Y[, , t])] <- NA_real_
		Y[, , t] <- yt
	}
	Y
}

# targets of the correction: coefficients, noise-corrected va / vb before their zero floor,
# cab and the double-centred product surface
.ae_bb_targets <- function(f) {
	vc <- .ae_jk_unfloored_vc(f)
	list(
		coef = f$coefficients,
		va = unname(vc[["va"]]),
		vb = unname(vc[["vb"]]),
		cab = unname(f$VC[["cab"]]),
		O = .ae_jk_product(f)
	)
}

# one round of replicates at `psi`: the mean bias of every target and the slope of the
# estimated effects on the drawn ones; rounds sharing `seed` differ only by their parameters
.ae_bb_round <- function(object, psi, n_rep, seed, max_iter, tol) {
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
		seed,
		kind = "Mersenne-Twister",
		normal.kind = "Inversion",
		sample.kind = "Rejection"
	))
	sym <- isTRUE(object$symmetric)
	rho <- unname(object$VC[["rho"]])
	if (!is.finite(rho) || sym || identical(object$mode, "bipartite")) rho <- 0
	active <- .ae_bb_active(object)
	EZ0 <- .ae_bb_eta(
		object,
		psi$coef,
		rep(0, object$dims$n_row),
		rep(0, object$dims$n_col),
		psi$O
	)
	recs <- vector("list", n_rep)
	for (r in seq_len(n_rep)) {
		eff <- .ae_bb_draw_effects(object, psi$va, psi$vb, psi$cab)
		ab <- outer(eff$a, eff$b, "+")
		Yb <- .ae_bb_draw_y(object, lapply(EZ0, `+`, ab), rho)
		fb <- tryCatch(
			suppressWarnings(suppressMessages(
				.ae_product_restore(ame_als_refit(
					object,
					Y_new = Yb,
					bias_correct = "none",
					max_iter = max_iter,
					tol = tol
				))
			)),
			error = function(e) NULL
		)
		if (is.null(fb) || !isTRUE(fb$converged) ||
				!all(is.finite(fb$coefficients))) next
		tg <- .ae_bb_targets(fb)
		ga <- eff$a[active$row]; gb <- eff$b[active$col]
		no_cab <- sym || identical(object$mode, "bipartite")
		both <- if (no_cab) NULL else active$row & active$col
		slope <- function(est, drawn) {
			ok <- is.finite(est) & is.finite(drawn)
			den <- sum(drawn[ok]^2)
			if (sum(ok) > 1L && den > 0) sum(est[ok] * drawn[ok]) / den else NA_real_
		}
		recs[[r]] <- list(
			coef = tg$coef,
			va = tg$va - stats::var(ga),
			vb = if (sym) NA_real_ else tg$vb - stats::var(gb),
			cab = if (no_cab) NA_real_ else
				tg$cab - stats::cov(eff$a[both], eff$b[both]),
			ka = slope(as.numeric(fb$a)[active$row], ga),
			kb = if (sym) NA_real_ else slope(as.numeric(fb$b)[active$col], gb),
			O = tg$O
		)
	}
	recs <- recs[!vapply(recs, is.null, logical(1))]
	if (!length(recs)) {
		return(list(n_valid = 0L, n_rep = n_rep))
	}
	mean_of <- function(nm) mean(vapply(recs, `[[`, numeric(1), nm), na.rm = TRUE)
	C <- do.call(rbind, lapply(recs, `[[`, "coef"))
	bias_O <- if (is.null(psi$O)) NULL else
		Reduce(`+`, lapply(recs, `[[`, "O")) / length(recs) - psi$O
	list(
		n_valid = length(recs),
		n_rep = n_rep,
		coef = colMeans(C) - psi$coef,
		coef_sd = apply(C, 2L, stats::sd),
		va = mean_of("va"),
		vb = mean_of("vb"),
		cab = mean_of("cab"),
		ka = mean_of("ka"),
		kb = mean_of("kb"),
		O = bias_O
	)
}

# actors with observed cells on each side (others have no estimated effect)
.ae_bb_active <- function(object) {
	obs <- apply(!is.na(object$Y), c(1, 2), sum)
	if (isTRUE(object$symmetric)) obs <- obs + t(obs)
	list(row = rowSums(obs) > 0, col = colSums(obs) > 0)
}

# secant factor s in [0, 1] (0 where b0 is zero or either bias is missing)
.ae_bb_secant <- function(b0, b1) {
	s <- 1 - b1 / b0
	s[!is.finite(s)] <- 0
	pmin(pmax(s, 0), 1)
}

# overall stretch of a product-surface bias against its generating surface
.ae_bb_stretch <- function(bias, O) {
	if (is.null(O) || is.null(bias)) return(NA_real_)
	den <- sum(O^2)
	if (den > 0) sum(bias * O) / den else NA_real_
}

# bias correction of an assembled static binary IRLS fit, before boundary actors are marked
# errors with class "lame_als_size" when the replicates do not refit
.ae_boot_bias <- function(object, max_iter, tol, n_rep = .ae_bb_reps) {
	est <- .ae_bb_targets(object)
	sym <- isTRUE(object$symmetric)
	# undirected fits with a multiplicative term are refused: the replicates do not
	# reproduce how far the fitted product is flattened
	if (sym && object$R > 0L) {
		cli::cli_abort(
			c(
				"Binary ALS cannot correct this undirected network with {.code R = {object$R}}.",
				"x" = "The split-sample correction does not reach it and the bootstrap bias correction does not hold for undirected fits with a multiplicative term: the replicates do not reproduce how far the fitted product is flattened, so the correction moves the dyadic slope away from the truth.",
				"i" = "Use the default empirical-Bayes estimator, {.fn ame} / {.fn lame} (MCMC), {.code R = 0} (corrected by the bootstrap), a panel of five or more periods (split-sample correction), or {.code bias_correct = \"none\"} for the uncorrected ALS estimate."
			),
			class = "lame_als_size"
		)
	}
	psi0 <- list(
		coef = est$coef,
		va = est$va,
		vb = if (sym) est$va else est$vb,
		cab = est$cab,
		O = est$O
	)
	refuse <- function(round, rd) {
		cli::cli_abort(
			c(
				"Binary ALS cannot correct this network: {rd$n_valid} of {rd$n_rep} simulated replicates of round {round} of the bootstrap bias correction refit and converged.",
				"i" = "The correction needs at least {round(100 * .ae_bb_min_valid)}% of them; networks simulated from this fit are too sparse for the fixed-effects estimator.",
				"i" = "Use {.fn ame} / {.fn lame} (MCMC), which integrate the effects out, or {.code bias_correct = \"none\"} for the uncorrected ALS estimate."
			),
			class = "lame_als_size"
		)
	}
	r0 <- .ae_bb_round(object, psi0, n_rep, .ae_bb_seed, max_iter, tol)
	if (r0$n_valid < .ae_bb_min_valid * n_rep) refuse(1L, r0)
	step <- function(x, b) if (is.null(x)) NULL else x - b
	psi1 <- list(
		coef = psi0$coef - r0$coef,
		va = psi0$va - r0$va,
		vb = psi0$vb - (if (sym) r0$va else r0$vb),
		cab = psi0$cab - .ae_na_zero(r0$cab),
		O = step(psi0$O, r0$O)
	)
	r1 <- .ae_bb_round(object, psi1, n_rep, .ae_bb_seed, max_iter, tol)
	if (r1$n_valid < .ae_bb_min_valid * n_rep) refuse(2L, r1)

	s_coef <- .ae_bb_secant(r0$coef, r1$coef)
	s_va <- .ae_bb_secant(r0$va, r1$va)
	s_vb <- if (sym) s_va else .ae_bb_secant(r0$vb, r1$vb)
	s_cab <- .ae_bb_secant(r0$cab, r1$cab)
	s_O <- .ae_bb_secant(.ae_bb_stretch(r0$O, psi0$O), .ae_bb_stretch(r1$O, psi1$O))
	coef <- est$coef - r0$coef / (1 + s_coef)
	va <- est$va - r0$va / (1 + s_va)
	vb <- if (sym) va else est$vb - r0$vb / (1 + s_vb)
	cab <- if (is.finite(r0$cab)) est$cab - r0$cab / (1 + s_cab) else est$cab
	O <- if (is.null(est$O)) NULL else est$O - r0$O / (1 + s_O)
	ka <- if (is.finite(r0$ka) && r0$ka > 0) r0$ka else 1
	kb <- if (sym) ka else if (is.finite(r0$kb) && r0$kb > 0) r0$kb else 1
	a <- as.numeric(object$a) / ka
	b <- if (sym) a else as.numeric(object$b) / kb

	out <- .ae_bb_rebuild(object, coef, a, b, O, va, vb, cab)
	out$bias_correct <- list(
		method = "bootstrap",
		applied = TRUE,
		replicates = c(round1 = r0$n_rep, round2 = r1$n_rep),
		valid = c(round1 = r0$n_valid, round2 = r1$n_valid),
		bias = list(
			coef = r0$coef,
			va = r0$va,
			vb = r0$vb,
			cab = r0$cab,
			effect_scale = c(a = ka, b = kb),
			O_stretch = .ae_bb_stretch(r0$O, psi0$O)
		),
		secant = list(coef = s_coef, va = s_va, vb = s_vb, cab = s_cab, O = s_O),
		correction = list(
			coef = est$coef - coef,
			va = est$va - va,
			vb = est$vb - vb,
			cab = est$cab - cab,
			effect_scale = c(a = ka, b = kb),
			O = if (is.null(O)) NULL else est$O - O
		),
		replicate_sd = r0$coef_sd
	)
	out
}

# corrected fit from its corrected parts; predictor, fitted values, residuals and
# deviance are rebuilt, and the uncorrected values kept in `raw`
.ae_bb_rebuild <- function(object, coef, a, b, O, va, vb, cab) {
	# restore a withheld product (.ae_withhold_product) to its variational means
	object <- .ae_product_restore(object)
	R <- object$R
	sym <- isTRUE(object$symmetric)
	U <- object$U; V <- object[["V"]]; L <- object$L
	O_new <- NULL
	if (R > 0L && !is.null(O)) {
		if (sym) {
			eg <- eigen((O + t(O)) / 2, symmetric = TRUE)
			idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
			U <- eg$vectors[, idx, drop = FALSE]
			L <- eg$values[idx]
			V <- U
			O_new <- U %*% (L * t(U))
		} else {
			sv <- svd(O)
			dd <- sqrt(pmax(sv$d[seq_len(R)], 0))
			U <- sweep(sv$u[, seq_len(R), drop = FALSE], 2, dd, "*")
			V <- sweep(sv$v[, seq_len(R), drop = FALSE], 2, dd, "*")
			O_new <- tcrossprod(U, V)
		}
		dimnames(U) <- list(object$row_names, paste0("dim", seq_len(R)))
		dimnames(V) <- list(object$col_names, paste0("dim", seq_len(R)))
	}
	EZ <- .ae_bb_eta(object, coef, a, b, O_new)
	bip <- identical(object$mode, "bipartite")
	lat <- .ae_default(object$pois_latent_var, 0)
	fitted_list <- resid_list <- vector("list", length(EZ))
	for (t in seq_along(EZ)) {
		ez <- unname(EZ[[t]])
		if (!bip) diag(ez) <- NA_real_
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

	raw <- list(
		coefficients = object$coefficients,
		a = object$a,
		b = object$b,
		VC = object$VC,
		O = NULL,
		deviance = object$deviance
	)
	if (R > 0L) {
		raw$O <- .ae_jk_product(object)
		if (!is.null(raw$O)) dimnames(raw$O) <- list(object$row_names, object$col_names)
	}
	out <- object
	names(coef) <- names(object$coefficients)
	out$coefficients <- coef
	out$mu <- unname(coef[1L])
	out$beta <- stats::setNames(as.numeric(coef[-1L]), names(object$beta))
	out$a <- stats::setNames(a, names(object$a))
	out$b <- stats::setNames(b, names(object$b))
	out$U <- U; out[["V"]] <- V; out$L <- L
	if (!is.null(O_new)) {
		Om <- O_new
		dimnames(Om) <- list(object$row_names, object$col_names)
		out$UVPM <- Om
		if (sym) out$ULUPM <- Om
	}
	VC <- object$VC
	VC[["va"]] <- if (is.finite(va)) max(va, 0) else NA_real_
	VC[["vb"]] <- if (sym) NA_real_ else if (is.finite(vb)) max(vb, 0) else NA_real_
	VC[["cab"]] <- cab
	out$VC <- VC
	out$vc_raw <- c(
		va = unname(va + object$vc_noise[["va"]]),
		vb = unname(vb + object$vc_noise[["vb"]])
	)
	out$EZ <- EZ
	out$fitted <- fitted_list
	out$YPM <- fitted_list
	out$residuals <- resid_list
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
	.ae_update_response_rho(out, out$link)
}

# apply the stored correction to a warm-start refit or bootstrap replicate, so a
# replicate does not repeat the replicate rounds
.ae_bb_apply_stored <- function(out, record) {
	cr <- record$correction
	if (is.null(cr)) return(out)
	sym <- isTRUE(out$symmetric)
	est <- .ae_bb_targets(out)
	coef <- est$coef - cr$coef
	va <- est$va - cr$va
	vb <- if (sym) va else est$vb - cr$vb
	cab <- if (is.finite(cr$cab)) est$cab - cr$cab else est$cab
	O <- if (is.null(est$O) || is.null(cr$O)) est$O else est$O - cr$O
	a <- as.numeric(out$a) / cr$effect_scale[["a"]]
	b <- if (sym) a else as.numeric(out$b) / cr$effect_scale[["b"]]
	res <- .ae_bb_rebuild(out, coef, a, b, O, va, vb, cab)
	res$bias_correct <- c(
		record[c("method", "applied")],
		list(stored = TRUE, correction = cr)
	)
	res
}
