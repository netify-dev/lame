# estimand of fitted(), predict(type = "response") and residuals() for a dynamic ALS fit:
# normal gives the posterior mean of the linear predictor, counts the plug-in mean
# exp(eta + ve / 2), and an empirical-Bayes binary fit the predictive probability
# Phi(m / sqrt(1 + v)) (the penalised fit uses the plug-in Phi(m)); fit$fitted and
# fit$YPM hold the plug-in values. v is the stored cavity variance of a cell plus the
# coefficients' part x' V x formed at the requested design

# regressors' contribution to the period-t linear predictor, `coef` times the
# intercept, dyadic and nodal covariates, on the design `Xn` [n_row, n_col, p_dyad, T]
.lda_xbeta_cells <- function(object, Xn, t, coef) {
	p_dyad <- dim(Xn)[3L]
	n_row <- object$dims$n_row
	n_col <- object$dims$n_col
	dyn_node_names <- character(0)
	if (!is.null(object$beta_dynamic_mask) && !is.null(object$beta_dynamic_groups)) {
		bm <- object$beta_dynamic_mask
		bg <- object$beta_dynamic_groups
		dyn_node_names <- names(bm)[bm & bg %in% c("row", "col")]
	}
	if (isTRUE(object$symmetric) && p_dyad > 0L) {
		for (k in seq_len(p_dyad)) Xn[, , k, t] <- .ae_symmetrize(Xn[, , k, t])
	}
	xb <- matrix(coef[1L], n_row, n_col)
	for (k in seq_len(p_dyad)) xb <- xb + coef[k + 1L] * Xn[, , k, t]
	cn <- rownames(object$coef_path)
	names(coef) <- cn
	if (isTRUE(object$symmetric)) {
		W_node <- .lda_bind_node_matrices(object$W_row, object$W_col)
		W_node_dyn <- .lda_bind_node_matrices(
			.lda_node_matrix(object, "row", t, dynamic = TRUE),
			.lda_node_matrix(object, "col", t, dynamic = TRUE)
		)
		if (!is.null(W_node)) {
			for (nm in intersect(colnames(W_node), cn)) {
				W_use <- if (nm %in% dyn_node_names && !is.null(W_node_dyn) && nm %in% colnames(W_node_dyn)) W_node_dyn else W_node
				g <- as.numeric(coef[[nm]] * W_use[, nm])
				xb <- xb + outer(g, g, "+")
			}
		}
	} else {
		if (!is.null(object$W_row)) {
			for (nm in intersect(colnames(object$W_row), cn)) {
				W_row_t <- .lda_node_matrix(object, "row", t, dynamic = nm %in% dyn_node_names)
				g <- as.numeric(coef[[nm]] * W_row_t[, nm])
				xb <- xb + matrix(g, n_row, n_col)
			}
		}
		if (!is.null(object$W_col)) {
			for (nm in intersect(colnames(object$W_col), cn)) {
				W_col_t <- .lda_node_matrix(object, "col", t, dynamic = nm %in% dyn_node_names)
				g <- as.numeric(coef[[nm]] * W_col_t[, nm])
				xb <- xb + matrix(g, n_row, n_col, byrow = TRUE)
			}
		}
	}
	xb
}

# empirical-Bayes binary fit, whose fitted values are predictive probabilities
.lde_is_binary_eb <- function(object) {
	identical(object$family, "binary") && identical(object$bias_correct$method, "eb") &&
		isTRUE(object$bias_correct$dynamic)
}

# record stored with the fit (R/lame_dynamic_als_eb.R): each cell's predictor variance
# apart from the coefficients [n_row, n_col, T], and each period's coefficient covariance [p, p, T]
.lde_predictive_record <- function(ctx, st, ev, Vc, nr, nc, Tt, p, hyp, rho, sigma) {
	va <- ev$lap$vcell
	vd <- ev$lap$vdiag
	na <- length(ctx$ia)
	nb <- length(ctx$ib)
	if (is.null(va) || is.null(vd) || !identical(as.integer(dim(va)), c(na, nb, Tt))) {
		cli::cli_abort(
			"The posterior variance of the additive predictor is not available for the predictive probability.",
			class = "lame_predictive_missing"
		)
	}
	# an actor with no observed cell keeps its prior: level variance plus the stationary
	# variance of its deviation path
	sab <- as.matrix(hyp$sab)
	stat <- if (length(rho) == 1L && length(sigma) == 1L && is.finite(rho) && abs(rho) < 1 && is.finite(sigma)) sigma^2 / (1 - rho^2) else 0
	prior_a <- sab[1L, 1L] + stat
	prior_b <- (if (nrow(sab) >= 2L) sab[2L, 2L] else sab[1L, 1L]) + stat
	var_effects <- array(0, c(nr, nc, Tt))
	for (t in seq_len(Tt)) {
		dg <- vd[[t]]
		wa <- rep(prior_a, nr)
		wb <- rep(prior_b, nc)
		wa[ctx$ia] <- dg[seq_len(na)]
		wb[ctx$ib] <- if (ctx$sym) dg[seq_len(na)] else dg[na + seq_len(nb)]
		v <- outer(wa, wb, "+")
		v[ctx$ia, ctx$ib] <- va[, , t]
		var_effects[, , t] <- v + if (is.null(st$Vpos)) 0 else {
			vp <- matrix(0, nr, nc)
			vp[ctx$ia, ctx$ib] <- st$Vpos[, , t]
			vp
		}
	}
	m <- ncol(Vc)
	unit <- vector("list", m)
	for (k in seq_len(m)) {
		e <- numeric(m)
		e[k] <- 1
		unit[[k]] <- .lde_coef_path(ctx, e)
	}
	coef_cov <- array(0, c(p, p, Tt), dimnames = list(ctx$names, ctx$names, ctx$time_names))
	for (t in seq_len(Tt)) {
		M <- vapply(unit, function(u) u[, t], numeric(p))
		coef_cov[, , t] <- M %*% Vc %*% t(M)
	}
	list(var_effects = var_effects, coef_cov = coef_cov)
}

# predictor variance of every cell of a design `Xn`, as a list over periods
.lde_predictive_var <- function(object, Xn) {
	pr <- object$eb$predictive
	if (is.null(pr)) {
		cli::cli_abort(
			c(
				"This fit does not carry the posterior variance of its linear predictor, which the predictive probability needs.",
				"i" = "Refit it with this version of the package, or use {.code fitted(fit, plug_in = TRUE)}."
			),
			class = "lame_predictive_missing"
		)
	}
	cn <- rownames(object$coef_path)
	p <- length(cn)
	n_row <- object$dims$n_row
	n_col <- object$dims$n_col
	lapply(seq_len(object$n_time), function(t) {
		R <- vapply(seq_len(p), function(k) {
			e <- numeric(p)
			e[k] <- 1
			as.numeric(.lda_xbeta_cells(object, Xn, t, e))
		}, numeric(n_row * n_col))
		pr$var_effects[, , t] + matrix(rowSums((R %*% pr$coef_cov[, , t]) * R), n_row, n_col)
	})
}

# predictive probabilities of the cells with predictor `EZ` (a list over periods) at design `Xn`
.lde_predictive_response <- function(object, EZ, Xn = object[["X"]]) {
	v <- .lde_predictive_var(object, Xn)
	out <- lapply(seq_along(EZ), function(t) .ae_eb_integrate(EZ[[t]], v[[t]], .ae_default(object$link, "probit")))
	names(out) <- names(EZ)
	out
}

# fitted values; `plug_in = TRUE` gives the plug-in Phi(m) of an empirical-Bayes binary fit
#' @export
fitted.lame_dynamic_als <- function(object, plug_in = FALSE, ...) {
	if (!isTRUE(plug_in) && .lde_is_binary_eb(object)) return(.lde_predictive_response(object, object$EZ))
	object$fitted
}

# response less the fitted values above (the working-scale residuals are the engine's)
#' @export
residuals.lame_dynamic_als <- function(object, type = c("response", "working"), plug_in = FALSE, ...) {
	type <- match.arg(type)
	if (identical(type, "working")) return(object$residuals)
	fv <- fitted.lame_dynamic_als(object, plug_in = plug_in)
	lapply(seq_len(object$n_time), function(t) {
		r <- object$Y[, , t] - fv[[t]]
		dimnames(r) <- dimnames(fv[[t]])
		r
	})
}
