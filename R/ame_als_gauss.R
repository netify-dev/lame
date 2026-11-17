# cell of the EB normal estimator (ame_als() and lame_als() with family = "normal"): y = eta + e,
# e ~ N(0, s2). score (y - eta) / s2 and information 1 / s2 do not depend on eta, so the Laplace
# approximation is exact; the dyadic correlation gives each directed pair covariance s2 [[1, rho], [rho, 1]]

# prior of s2 in ame() and lame() for the normal family: scaled inverse chi-square with one degree
# of freedom on s20, by default the variance of the observed responses (R/residual_variance_prior.R)
.ae_gauss_ve_prior <- function(Yarr, s20 = NULL) {
	c(nu0 = 1, s20 = .lame_residual_prior_scale(Yarr, "normal", s20))
}

# additive-effects prior of ame() and lame() for normal data: scale .lame_ab_scale (R/normal_scale.R)
# and the weight of one actor (.ae_eb_prior_weight; the sampler's is round(4 + 3 n / 100), n the
# number of row actors)
.ae_eb_hyperprior_norm <- function(Yarr) {
	nact <- dim(Yarr)[1L]
	list(Sab0 = diag(2) * .lame_ab_scale(Yarr), eta0 = .ae_eb_prior_weight(round(4 + 3 * nact / 100), min(dim(Yarr)[1:2])))
}

# floors for the noise variance and effect covariance: relative to the response scale for a normal
# fit, else 1e-8
.ae_ve_floor <- function(ab) if (!is.null(ab) && identical(ab$family, "normal")) 1e-8 * ab$ve_prior[["s20"]] else 1e-8
.ae_sab_floor <- function(ab) if (!is.null(ab) && identical(ab$family, "normal")) 1e-8 * ab$Sab0[1L, 1L] else 1e-8


# cell terms at `eta` and noise variance `ve`: log likelihood l, score s, information w, its eta
# derivative hp (zero) and noise second moment ee, all zero off the observed cells
.ae_gauss_terms <- function(y, eta, ve) {
	ok <- is.finite(y) & is.finite(eta)
	r <- y - eta
	r[!ok] <- 0
	r2 <- r * r
	z0 <- array(0, .ae_dims(y))
	list(l = (-0.5 * log(2 * pi * ve) - 0.5 / ve * r2) * ok, s = r / ve, w = ok / ve, hp = z0, ee = r2, ve = ve)
}

# Gaussian working response z = y and weight 1 / ve (no increment cap)
.ae_irls_gauss <- function(y, terms) {
	z <- y
	z[!is.finite(z)] <- 0
	w <- terms$w
	list(z = z, w = w)
}

# expected working problem for cells with predictor N(m, v): z = y, weight 1 / ve, expected log
# likelihood and ee, the expected noise second moment for the EM image of ve (zero off the observed cells)
.ae_gauss_cells <- function(y, m, v, ve) {
	obs <- is.finite(y) & is.finite(m)
	ix <- which(obs)
	vv <- pmax(rep_len(v, length(y))[ix], 0)
	r2 <- (y[ix] - m[ix])^2
	z <- w <- ell <- ee <- array(0, .ae_dims(y))
	z[ix] <- y[ix]
	w[ix] <- 1 / ve
	ell[ix] <- -0.5 * log(2 * pi * ve) - 0.5 * (r2 + vv) / ve
	ee[ix] <- r2 + vv
	list(z = z, w = w, ell = ell, ee = ee)
}

# cavity information weight for the factor-row covariance (see .ae_mf_info_w): 1 / (ve + v)
.ae_gauss_info_w <- function(y, m, v, ve) {
	out <- array(0, .ae_dims(y))
	ix <- which(is.finite(y) & is.finite(m))
	out[ix] <- 1 / (ve + pmax(rep_len(v, length(y))[ix], 0))
	out
}


# residual variance of a state on the observed cells, corrected for the mean model's degrees of
# freedom (.ae_assemble)
.ae_gauss_resid_var <- function(prep, st, R, symmetric) {
	y <- prep$Yarr
	nr <- prep$nr; nc <- prep$nc; Tt <- prep$Tt
	Om <- matrix(0, nr, nc)
	if (R > 0L && !is.null(st$U)) {
		Om <- if (symmetric && !is.null(st$L)) st$U %*% (st$L * t(st$U)) else if (!is.null(st[["V"]])) tcrossprod(st$U, st[["V"]]) else Om
	}
	eta <- array(st$mu + outer(st$a, st$b, "+") + Om, dim(y))
	if (prep$p_dyad > 0L) {
		xmat <- matrix(aperm(prep$X, c(1L, 2L, 4L, 3L)), nr * nc * Tt, prep$p_dyad)
		eta <- eta + as.vector(xmat %*% st$beta)
	}
	ok <- is.finite(y) & is.finite(eta)
	cnt <- rowSums(ok, dims = 2L)
	if (symmetric) cnt <- cnt + t(cnt)
	df <- .ae_df_mean(prep$p_dyad, sum(rowSums(cnt) > 0), sum(colSums(cnt) > 0), R, symmetric)
	if (symmetric) df <- 2L * df
	ve <- sum((y[ok] - eta[ok])^2) / max(sum(ok) - df, 1)
	max(ve, 1e-8 * stats::var(y[ok]), 1e-12)
}

# start of the penalised normal fit: a short fixed-effects fit and its residual variance (the first
# weight 1 / ve); a supplied start keeps its state and variance, and `ve_hold` is a variance held fixed
.ae_gauss_seed <- function(prep, Z0, R, symmetric, max_iter, tol, init, lowrank_method, linear_solver, ve_hold = NULL) {
	st <- init
	if (is.null(st) || is.null(st$mu)) {
		fe <- .ae_bcd_fit(
			Z0,
			prep$X,
			R,
			symmetric,
			min(max_iter, 50L),
			tol,
			lowrank_method = lowrank_method,
			linear_solver = linear_solver,
			unit_rules = TRUE
		)
		st <- list(mu = fe$mu, beta = fe$beta, a = fe$a, b = fe$b, U = fe$U, V = fe[["V"]], L = fe$L)
	}
	ve0 <- .ae_gauss_resid_var(prep, st, R, symmetric)
	ve <- if (!is.null(ve_hold)) ve_hold else if (!is.null(init$ve)) init$ve else ve0
	st$ve <- ve
	list(init = st, ve = ve)
}


# reweighting steps of the first mode of a normal fit, which only starts the variational stage
# (the working problem is linear, so the hyperparameters settle quickly)
.ae_gauss_start_iters <- 4L
# mode tolerance for the evidence (the certificate threshold is its square root)
.ae_gauss_evidence_tol <- 1e-4

# certified posterior mode at the fitted noise variance, solved from the variational state with the
# variance held fixed; its Laplace marginal likelihood is the rank evidence (fit$eb$log_marginal) and
# a mode that does not settle is refused by the caller
.ae_gauss_evidence <- function(
	prep,
	ir,
	ab,
	R,
	symmetric,
	link,
	max_iter,
	tol,
	lowrank_method = "mm",
	linear_solver = "eigen"
) {
	f <- ir$fit
	# the evidence needs the mode only to 1e-2 of its curvature scale, since the log marginal is quadratic
	# in the distance from the mode
	tol <- max(tol, .ae_gauss_evidence_tol)
	run <- .ae_irls_run(
		prep,
		"normal",
		link,
		R,
		symmetric,
		max_iter,
		tol,
		lowrank_method = lowrank_method,
		init = .ae_pln_mode_init(f, f$ve),
		linear_solver = linear_solver,
		warm = TRUE,
		ab_prior = ab,
		ve_fixed = f$ve,
		warm_retry = TRUE,
		vcov_mode = "evidence"
	)
	ev <- run$fit$eb_vcov$evidence
	if (!is.null(ev) && is.finite(ev$log_marginal)) {
		ir$fit$eb_vcov$evidence <- ev
		if (!is.null(ir$fit$eb_mf)) ir$fit$eb_mf$record$mode_log_marginal <- ev$log_marginal
	}
	ir$stationarity <- run$stationarity
	ir$fit$converged <- isTRUE(run$stationarity$pass)
	ir$irls_iters <- ir$irls_iters + run$irls_iters
	ir$selection_converged <- isTRUE(ir$selection_converged) && isTRUE(run$stationarity$pass)
	ir
}


# signal of the fitted product against the noise per loading: the squared norm of the double-centred
# mean product, summed over periods, over ve and the number of loadings (see .ae_withhold_product)
.ae_gauss_information <- function(cx, Pm, ve, loadings) {
	if (is.null(Pm)) return(NULL)
	P <- Pm - rowMeans(Pm) - rep(colMeans(Pm), each = nrow(Pm)) + mean(Pm)
	signal <- cx$Tt * sum(P^2)
	c(signal = signal, loadings = loadings, per_loading = signal / (ve * loadings), threshold = .ae_product_threshold)
}

# the same statistic for a binary or count fit, weighting each observed cell by its Fisher weight
.ae_glm_information <- function(cx, Pm, W, loadings, events) {
	if (is.null(Pm)) return(NULL)
	P <- Pm - rowMeans(Pm) - rep(colMeans(Pm), each = nrow(Pm)) + mean(Pm)
	Wa <- array(W, c(nrow(P), ncol(P), length(W) %/% length(P)))
	signal <- sum(Wa * as.vector(P^2), na.rm = TRUE)
	c(
		events = events,
		loadings = loadings,
		per_loading = events / loadings,
		signal = signal,
		signal_per_loading = signal / loadings,
		threshold = .ae_product_threshold
	)
}

# information-weighted signal per loading below which the data barely determine an EB latent product;
# a dynamic fit withholds such a product (R/lame_dynamic_als_eb.R), a static fit reports it with its
# standard errors (see ?ame_als)
.ae_product_threshold <- 2.1

# statistic the threshold reads: noise-variance based for normal, information weighted for binary
# and counts (NA for a record without one)
.ae_product_statistic <- function(family, info) {
	if (is.null(info) || !"signal" %in% names(info)) return(NA_real_)
	if (identical(family, "normal")) info[["per_loading"]] else if ("signal_per_loading" %in% names(info)) info[["signal_per_loading"]] else NA_real_
}

.ae_product_is_withheld <- function(object) !is.null(object$eb$product_withheld)

# the fit with its variational means back in the withheld slots, for internal computations that start
# from the fitted state; the withholding record stays. `state = FALSE` keeps the reported
# (identified) values of the intercept, effects and positions (.ae_report_fold)
.ae_product_restore <- function(object, state = TRUE) {
	object <- .lde_intercept_restore(object)
	pw <- object$eb$product_withheld
	if (!is.null(pw) && !is.null(pw$means)) {
		m <- pw$means
		# the means were withheld as reported; a fit already in its fitted state takes that state's
		rf <- object$eb$report_fold
		if (!is.null(rf) && identical(rf$view, "state")) {
			k <- intersect(names(m), names(rf$state))
			m[k] <- rf$state[k]
		}
		for (k in names(m)) object[k] <- list(m[[k]])
	}
	if (isTRUE(state)) .ae_report_set(object, "state") else object
}

# slots the identified report rewrites
.ae_report_slots <- c("mu", "coefficients", "a", "b", "U", "V", "L", "UVPM", "ULUPM")

# identified report of an empirical-bayes fit, as the sampler stores its draws: the positions are
# centred over the actors and the row, column and grand means the product loses go to a, b and the
# intercept, so the linear predictor is unchanged. fit$eb$report_fold keeps the fitted state the
# internal computations start from beside the reported values (.ae_report_set)
.ae_report_fold <- function(out) {
	R <- out$R
	U <- out$U; V <- out[["V"]]; Om <- out$UVPM
	if (is.null(out$eb) || is.null(R) || R < 1L || !is.matrix(U) || !is.matrix(V) || !is.matrix(Om) ||
			!all(is.finite(U)) || !all(is.finite(V))) {
		return(out)
	}
	sym <- isTRUE(out$symmetric)
	cen <- if (sym) .lame_center_uv_period(U, U, diag(out$L, R)) else .lame_center_uv_period(U, V, diag(1, R))
	if (is.null(cen)) return(out)
	slots <- intersect(.ae_report_slots, names(out))
	state <- out[slots]
	a_row <- cen$a_row; b_col <- cen$b_col
	# one effect per actor: the row and column halves agree, and averaging keeps that under rounding
	if (sym) a_row <- b_col <- (a_row + b_col) / 2
	out$mu <- out$mu + cen$const
	out$coefficients[1L] <- out$coefficients[1L] + cen$const
	out$a <- out$a + a_row
	out$b <- out$b + b_col
	Oc <- Om - outer(a_row, rep(1, ncol(Om))) - outer(rep(1, nrow(Om)), b_col) - cen$const
	out$UVPM <- Oc
	if (!is.null(out$ULUPM)) out$ULUPM <- Oc
	if (sym) {
		# the centred product keeps its rank: orthonormal positions and eigenvalues again, each
		# column signed as the centred column it loads on most
		eg <- eigen(cen$U %*% (out$L * t(cen$U)), symmetric = TRUE)
		idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
		Un <- eg$vectors[, idx, drop = FALSE]
		ld <- crossprod(cen$U, Un)
		sgn <- vapply(seq_len(R), function(r) if (ld[which.max(abs(ld[, r])), r] < 0) -1 else 1, 0)
		Un <- sweep(Un, 2L, sgn, "*")
		dimnames(Un) <- dimnames(U)
		out$U <- Un
		out[["V"]] <- Un
		out$L <- eg$values[idx]
	} else {
		out$U <- cen$U
		out[["V"]] <- cen[["V"]]
	}
	out$eb$report_fold <- list(
		view = "report",
		state = state,
		report = out[slots],
		const = cen$const,
		a_row = a_row,
		b_col = b_col
	)
	out
}

# the fit with its fitted state ("state") or its reported values ("report") in the slots of the
# identified report; a fit with no record, or already in that view, is returned as it is
.ae_report_set <- function(object, view) {
	rf <- object$eb$report_fold
	if (is.null(rf) || identical(rf$view, view)) return(object)
	m <- rf[[view]]
	for (k in names(m)) object[k] <- list(m[[k]])
	object$eb$report_fold$view <- view
	object
}
.ae_report_view <- function(object) .ae_report_set(object, "report")

# sentence print() and summary() give for a withheld product
.ae_product_withheld_text <- function(pw) {
	k <- pw$rank_chosen
	unit <- if (isTRUE(pw$moving)) " noise variances of signal per loading of a period" else
	if (identical(.ae_default(pw$family, "normal"), "normal")) " noise variances of signal per loading" else
		" units of signal per loading, the product weighted by the information of each cell"
	paste0(
		"The latent product is not identified by the data at this rank (",
		signif(pw$per_loading, 2),
		unit,
		", below ",
		pw$threshold,
		"); it is withheld (the positions {.code U}, {.code V} and the product {.code UVPM} are NA). ",
		if (isTRUE(is.finite(k))) paste0(
			"The rank evidence ({.code fit$eb$product_withheld$log_marginal}) chooses rank ",
			k,
			". "
		) else "",
		"Fitted values, predictions and residuals carry the variational posterior of the predictor."
	)
}
