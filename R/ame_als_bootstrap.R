# bootstrap inference for the fast AME estimator (ame_als() / lame_als()), adapted from the social
# influence regression work of Hoff and Minhas (sir::boot_sir(); Minhas and Hoff 2025, Political Analysis)

# internal helpers

# a fresh additive contribution for the full parametric bootstrap: redrawing a and b from Sab per
# replicate gives the intercept and node-covariate standard errors their actor-level component
.ae_boot_draw_ab <- function(object) {
	nr <- object$dims$n_row; nc <- object$dims$n_col
	va <- as.numeric(object$VC["va"])
	vb <- as.numeric(object$VC["vb"])
	cab <- as.numeric(object$VC["cab"])
	if (!is.finite(va)) va <- 0
	if (!is.finite(vb)) vb <- 0
	if (!is.finite(cab)) cab <- 0
	# draws are centred to the fit's sum-to-zero gauge, so the refit intercept does not absorb
	# identification noise; empirical-Bayes fits keep uncentred effects (mu is the population intercept)
	eb_fit <- identical(object$bias_correct$method, "eb")
	ctr <- function(x) if (length(x) > 1L && !eb_fit) x - mean(x) else x
	# intercept the drawn effects induce through the node covariates ("intercept_shift"); the replicate
	# intercept of a fixed-effects fit is measured against it
	shift <- function(a, b) {
		if (eb_fit) return(0)
		if (object$symmetric) {
			W <- cbind(object$W_row, object$W_col)
			return(2 * .ae_decompose(a, W)$intercept)
		}
		.ae_decompose(a, object$W_row)$intercept + .ae_decompose(b, object$W_col)$intercept
	}
	finish <- function(a, b) {
		out <- outer(a, b, "+")
		attr(out, "intercept_shift") <- as.numeric(shift(a, b))
		out
	}
	if (object$symmetric) {
		a <- ctr(stats::rnorm(nr, 0, sqrt(max(va, 0))))
		return(finish(a, a))
	}
	if (identical(object$mode, "bipartite")) {
		a <- ctr(stats::rnorm(nr, 0, sqrt(max(va, 0))))
		b <- ctr(stats::rnorm(nc, 0, sqrt(max(vb, 0))))
		return(finish(a, b))
	}
	# unipartite directed: joint (a_i, b_i) draw preserving cab
	Sab <- matrix(c(va, cab, cab, vb), 2, 2)
	eg <- eigen((Sab + t(Sab)) / 2, symmetric = TRUE)
	eg$values[eg$values < 0] <- 0
	Sroot <- eg$vectors %*% (sqrt(eg$values) * t(eg$vectors))
	ab <- matrix(stats::rnorm(2 * nr), nr, 2) %*% Sroot
	finish(ctr(ab[, 1]), ctr(ab[, 2]))
}

# additive effects of a replicate: an actor a binary refit leaves out at the boundary records +/-Inf,
# so its interval is open on that side; actors the original fit left out stay NA, and `ok` is FALSE
# when any other effect is not finite
.ae_boot_effects <- function(fit_b, object) {
	a <- as.numeric(fit_b$a)
	b <- as.numeric(fit_b$b)
	sym <- isTRUE(fit_b$symmetric)
	orig_a <- is.na(object$a)
	orig_b <- is.na(object$b)
	ne <- fit_b$not_estimated
	listed_a <- rep(FALSE, length(a))
	listed_b <- rep(FALSE, length(b))
	if (!is.null(ne) && nrow(ne)) {
		row_side <- ne$side %in% c("sender", "row", "actor")
		bound <- ifelse(
			ne$outcome == "only ties",
			Inf,
			ifelse(ne$outcome == "no tie", -Inf, NA_real_)
		)
		ia <- match(ne$actor[row_side], fit_b$row_names)
		keep <- !is.na(ia)
		a[ia[keep]] <- bound[row_side][keep]
		listed_a[ia[keep]] <- TRUE
		if (sym) {
			b <- a
			listed_b <- listed_a
		} else {
			ib <- match(ne$actor[!row_side], fit_b$col_names)
			keep <- !is.na(ib)
			b[ib[keep]] <- bound[!row_side][keep]
			listed_b[ib[keep]] <- TRUE
		}
	}
	a[orig_a] <- NA_real_
	b[orig_b] <- NA_real_
	# every effect the original fit estimates is a number in the replicate, and
	# an infinite one only where the refit reports the actor as not estimated
	ok <- !anyNA(a[!orig_a]) && !anyNA(b[!orig_b]) &&
		!any(is.infinite(a) & !listed_a) && !any(is.infinite(b) & !listed_b)
	list(a = a, b = b, ok = ok)
}

# replicate sd of each column over valid replicates: NA with no value, Inf when a replicate is at the boundary
.ae_boot_col_sd <- function(M, v) {
	apply(M[v, , drop = FALSE], 2, function(x) {
		x <- x[!is.na(x)]
		if (!length(x)) return(NA_real_)
		if (any(is.infinite(x))) return(Inf)
		if (length(x) < 2L) return(NA_real_)
		stats::sd(x)
	})
}

# orthogonal Procrustes rotation R minimizing || M_src R - M_tgt ||
.ae_procrustes <- function(M_src, M_tgt) {
	sv <- tryCatch(svd(crossprod(M_src, M_tgt)), error = function(e) NULL)
	if (is.null(sv)) return(diag(ncol(M_src)))
	sv$u %*% t(sv$v)
}

# minimal "prep" list (as .ae_assemble() consumes) from canonical arrays, for ame_als_refit()
.ae_prep_from_canonical <- function(
	Yarr,
	X,
	W_row,
	W_col,
	mode,
	symmetric,
	row_names,
	col_names,
	time_names
) {
	d <- dim(Yarr)
	xn <- if (dim(X)[3] > 0 && !is.null(dimnames(X))) {
		dimnames(X)[[3]]
	} else {
		character(0)
	}
	list(
		Yarr = Yarr,
		X = X,
		p_dyad = dim(X)[3],
		x_names_dyad = xn,
		W_row = W_row,
		W_col = W_col,
		nr = d[1],
		nc = d[2],
		Tt = d[3],
		bip = identical(mode, "bipartite"),
		row_names = row_names,
		col_names = col_names,
		time_names = time_names
	)
}

# dynamic refits

# dynamic ALS fits carry time-indexed parameters, so they are refit by the fitter that produced them
.ae_is_dynamic_als <- function(object) {
	inherits(object, "lame_dynamic_als")
}

# static seed for a dynamic refit: each trajectory collapsed to its time mean (the shape the dynamic
# initialiser expands back into paths); symmetric fits keep their eigenvalue sign pattern
.ae_static_seed <- function(object) {
	# restore a withheld product (.ae_withhold_product) to its variational means
	object <- .ae_product_restore(object)
	flatten <- function(x) {
		if (is.null(x) || length(dim(x)) != 3L) return(x)
		out <- apply(x, c(1, 2), mean, na.rm = TRUE)
		out[!is.finite(out)] <- 0
		out
	}
	U <- flatten(object$U)
	V <- flatten(object[["V"]])
	L <- object$L
	if (isTRUE(object$symmetric) && is.null(L) && is.matrix(U) &&
			is.matrix(V) && ncol(U) > 0L) {
		L <- vapply(seq_len(ncol(U)), function(r) {
			cs <- sum(U[, r] * V[, r], na.rm = TRUE)
			if (!is.finite(cs) || cs >= 0) 1 else -1
		}, numeric(1))
	}
	cube <- function(x) if (length(dim(x)) == 3L) x else NULL
	path <- function(x) if (is.matrix(x)) x else NULL
	mu <- object$mu
	coefficients <- object$coefficients
	coef_path <- path(object$coef_path)
	if (identical(object$family, "poisson") && identical(object$link, "log") &&
			!identical(object$bias_correct$method, "eb")) {
		# IRLS works on the count log mean while reported Poisson predictors subtract half the latent
		# variance, so the warm start adds it back (empirical-Bayes count fits report the latent predictor)
		shift <- .ae_default(object$pois_latent_var, 0) / 2
		mu <- mu + shift
		if ("intercept" %in% names(coefficients)) {
			coefficients["intercept"] <- coefficients["intercept"] + shift
		}
		if (!is.null(coef_path) && "intercept" %in% rownames(coef_path)) {
			coef_path["intercept", ] <- coef_path["intercept", ] + shift
		}
	}
	list(
		mu = mu,
		beta = object$beta,
		coefficients = coefficients,
		a = as.numeric(object$a),
		b = as.numeric(object$b),
		U = U,
		V = V,
		L = L,
		G = object[["G"]],
		s2 = object$s2,
		# the trajectories themselves, used by a fitter that can seed paths
		coef_path = coef_path,
		a_path = path(object$a_dynamic),
		b_path = path(object$b_dynamic),
		U_path = cube(object$U),
		V_path = cube(object[["V"]]),
		# what an empirical-bayes fit adds: the levels of the additive paths and
		# the variance components, which a refit starts its updates from
		a_level = object$a_level,
		b_level = object$b_level,
		eb = object$eb
	)
}

# per-slice inputs for a dynamic fitter, rebuilt from the stored canonical arrays (with the
# per-slice actor sets of a changing-composition panel)
.ae_dynamic_refit_data <- function(object, Yarr, Xarr) {
	Tt <- dim(Yarr)[3L]
	rn <- dimnames(Yarr)[[1L]] %||% object$row_names
	cn <- dimnames(Yarr)[[2L]] %||% object$col_names
	tn <- dimnames(Yarr)[[3L]] %||% object$time_names
	changing <- isTRUE(object$changing_composition)
	row_pres <- object$row_presence %||% object$actor_presence
	col_pres <- object$col_presence %||% object$actor_presence
	p <- if (!is.null(Xarr) && length(dim(Xarr)) == 4L) dim(Xarr)[3L] else 0L
	xn <- if (p > 0L) dimnames(Xarr)[[3L]] else NULL
	Ylist <- vector("list", Tt)
	Xlist <- if (p > 0L) vector("list", Tt) else NULL
	for (t in seq_len(Tt)) {
		ri <- if (changing && !is.null(row_pres)) {
			which(row_pres[, t])
		} else seq_len(dim(Yarr)[1L])
		ci <- if (changing && !is.null(col_pres)) {
			which(col_pres[, t])
		} else seq_len(dim(Yarr)[2L])
		Ylist[[t]] <- matrix(
			Yarr[ri, ci, t],
			length(ri),
			length(ci),
			dimnames = list(rn[ri], cn[ci])
		)
		if (p > 0L) {
			Xlist[[t]] <- array(
				Xarr[ri, ci, , t],
				dim = c(length(ri), length(ci), p),
				dimnames = list(rn[ri], cn[ci], xn)
			)
		}
	}
	names(Ylist) <- tn
	if (!is.null(Xlist)) names(Xlist) <- tn
	list(Y = Ylist, Xdyad = Xlist)
}

# warm restart of a dynamic als fit: the fitter that produced it, run on the
# same configuration and seeded from the stored trajectories
.ae_dynamic_refit <- function(
	object,
	Y_new = NULL,
	X_new = NULL,
	max_iter = NULL,
	tol = NULL,
	verbose = FALSE,
	bias_correct = NULL
) {
	# restore a withheld product (.ae_withhold_product) to its variational means
	object <- .ae_product_restore(object)
	cfg <- object$als_refit
	if (is.null(cfg) || is.null(cfg$args) || is.null(cfg$fitter)) {
		cli::cli_abort(c(
			"This dynamic ALS fit does not carry the arguments a refit needs.",
			"i" = "Refit the model with the current version of {.fn lame} before refitting or bootstrapping it."
		))
	}
	fitter <- switch(
		cfg$fitter,
		lame_dynamic_als = lame_dynamic_als,
		cli::cli_abort("Unknown dynamic ALS fitter {.val {cfg$fitter}}.")
	)
	Yarr <- if (is.null(Y_new)) object$Y else Y_new
	Xarr <- if (is.null(X_new)) object[["X"]] else X_new
	dat <- .ae_dynamic_refit_data(object, Yarr, Xarr)
	args <- cfg$args
	if (!is.null(bias_correct)) args$bias_correct <- bias_correct
	args$Y <- dat$Y
	args$Xdyad <- dat$Xdyad
	args$verbose <- isTRUE(verbose)
	args$stability <- "none"
	args$seed <- cfg$seed
	# an empirical-Bayes count fit refit with another estimator starts cold: its posterior-mean paths
	# and latent predictor are not a valid start for the penalised estimators
	if (!(identical(object$family, "poisson") && identical(object$bias_correct$method, "eb") &&
			!is.null(bias_correct) && !identical(bias_correct, "eb"))) args$start <- .ae_static_seed(object)
	else args$start <- NULL
	if ("bootstrap" %in% names(formals(fitter))) args$bootstrap <- 0L
	if (!is.null(max_iter)) args$max_iter <- max_iter
	if (!is.null(tol)) args$tol <- tol
	out <- do.call(fitter, args)
	out$call <- object$call
	out
}

# warm-start refit

#' Refit a fast AME model with a warm start
#'
#' @description
#' Refits an AME model by iterative block coordinate descent, initialised
#' (\dQuote{warm-started}) from an existing \code{ame_als} fit. This is the
#' workhorse of \code{\link{ame_als_bootstrap}}: every bootstrap replicate
#' is refit from the original point estimate rather than from a cold random
#' start, which prevents replicates from converging to different local optima
#' or rotations and is essential for meaningful bootstrap standard errors.
#'
#' @details
#' The estimation algorithm is the iterative block coordinate descent estimator
#' of the Social Influence Regression model of Hoff & Minhas
#' (\code{sir::sir_alsfit()}), adapted to the AME model, with each bootstrap
#' replicate warm-started from the original point estimate.
#'
#' A fit from one of the dynamic ALS routes (\code{lame(..., method = "als")}
#' with \code{dynamic_uv}, \code{dynamic_ab} or \code{dynamic_beta}) stores
#' trajectories rather than a single
#' value per parameter. Such a fit is refit through the dynamic fitter that
#' produced it, seeded from the stored trajectories, and the result is again a
#' dynamic fit of the same class. \code{Z_new} is not available on that route,
#' which runs its own working-response loop.
#'
#' A static fit carrying the split-sample jackknife (see \code{\link{ame_als}})
#' is refit with it: the refit on the new data is followed by its own two
#' half-sample refits, warm-started from it with the same \code{max_iter} and
#' \code{tol}, and the corrected values are returned, with the refit's
#' uncorrected values in \code{$raw}.
#'
#' @param object an \code{ame_als} object supplying the warm-start values
#'   (\code{mu}, \code{beta}, \code{a}, \code{b}, \code{U}, \code{V}) and the
#'   model configuration (\code{family}, \code{mode}, \code{symmetric},
#'   \code{R}).
#' @param Y_new optional canonical outcome array \code{[n_row, n_col, T]} to
#'   refit on. Defaults to \code{object$Y}.
#' @param X_new optional canonical design array \code{[n_row, n_col, p, T]}.
#'   Defaults to \code{object$X}.
#' @param Z_new optional canonical \emph{working-response} array
#'   \code{[n_row, n_col, T]}. When supplied, the model is refit by a Gaussian
#'   block coordinate descent directly on \code{Z_new}, bypassing the
#'   family transform / IRLS reweighting. This is used by the parametric
#'   bootstrap of a non-normal \code{transform} fit, whose estimator is a
#'   Gaussian fit to a fixed transformed response. \code{Y_new} is still used
#'   for the observed-cell pattern.
#' @param max_iter maximum block coordinate descent iterations (default 30;
#'   fewer are needed than for a cold start). On the dynamic route, left
#'   unset it keeps the original fit's iteration budget.
#' @param tol convergence tolerance (default \code{1e-5}). On the dynamic
#'   route, left unset it keeps the original fit's tolerance.
#' @param verbose logical; print progress (default \code{FALSE}).
#' @param bias_correct \code{NULL} (default) applies the bias correction
#'   \code{object} was fit with (\code{object$bias_correct$method}), so a
#'   refit of a jackknife-corrected fit is corrected the same way, a refit
#'   of a fit corrected by the parametric bootstrap has the bias that fit
#'   estimated removed (\code{fit$bias_correct$stored} is \code{TRUE}) and a
#'   refit of an empirical-Bayes fit (\code{"eb"}) is an empirical-Bayes fit,
#'   which re-estimates the covariance of the actor effects starting from the
#'   fit's; \code{"jackknife"}, \code{"bootstrap"} (re-estimating the bias),
#'   \code{"none"} or \code{"eb"} override it. On a refit the
#'   jackknife is computed or refused: a panel the correction cannot carry
#'   (a single cross-section at \code{R > 0} or below 25 rarer-outcome cells
#'   per actor, fewer than five periods at \code{R > 0}, fewer than 25
#'   rarer-outcome cells per actor at \code{R = 0}, or an actor without both
#'   outcome values in the panel or, under every split, in a half sample;
#'   see \code{\link{ame_als}}) raises an error of class
#'   \code{"lame_jk_split"}. A refit of an uncorrected fit keeps its
#'   \code{bias_correct} record. Dynamic Poisson fits preserve their estimator
#'   (the empirical-Bayes estimator, or the penalised fit with the
#'   variance/covariance correction of directed panels) and override it in the
#'   same way (see \code{bias_correct} in \code{\link{lame}}); an override of
#'   the empirical-Bayes estimator by a penalised one starts cold.
#'
#' @return An object of class \code{"ame_als"}; see \code{\link{ame_als}}.
#'
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}. The iterative block
#' coordinate descent estimator refit here originates with that work
#' (implemented in \code{sir::sir_alsfit()}).
#'
#' @seealso \code{\link{ame_als}}, \code{\link{ame_als_bootstrap}}.
#'
#' @examples
#' Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
#' fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' refit <- ame_als_refit(fit, verbose = FALSE)
#' coef(refit)
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
ame_als_refit <- function(
	object,
	Y_new = NULL,
	X_new = NULL,
	Z_new = NULL,
	max_iter = 30,
	tol = 1e-5,
	verbose = FALSE,
	bias_correct = NULL
) {
	if (!inherits(object, "ame_als")) {
		cli::cli_abort("{.arg object} must be an {.cls ame_als} object.")
	}
	# restore a withheld product to its variational means, where the refit starts
	object <- .ae_product_restore(object)
	bias_correct <- .ae_check_bias_correct(bias_correct)
	# a dynamic fit is refit by the fitter that produced it, which decides what
	# it supports; a static fit needs the binary irls estimator
	if (identical(bias_correct, "eb") && !.ae_is_dynamic_als(object) &&
			!((object$family %in% c("binary", "poisson") &&
					identical(.ae_default(object$non_normal_method, "transform"), "irls")) ||
					identical(object$family, "normal"))) {
		cli::cli_abort(
			c(
				"{.code bias_correct = \"eb\"} is available for binary, Poisson and normal fits only.",
				"i" = "This fit's family is {.val {object$family}}."
			),
			class = "lame_als_eb"
		)
	}
	if (.ae_is_dynamic_als(object)) {
		# a dynamic fit is refit dynamically: the static block coordinate
		# descent has no time dimension to put the trajectories in
		if (!is.null(Z_new)) {
			cli::cli_abort(c(
				"{.arg Z_new} is not available for a dynamic ALS fit.",
				"i" = "The dynamic estimator runs its own working-response loop; refit on {.arg Y_new}."
			))
		}
		out <- .ae_dynamic_refit(
			object,
			Y_new = Y_new,
			X_new = X_new,
			max_iter = if (missing(max_iter)) NULL else max_iter,
			tol = if (missing(tol)) NULL else tol,
			verbose = verbose,
			bias_correct = bias_correct
		)
		if (verbose) {
			cli::cli_text("Warm-start refit (dynamic): converged = {.val {out$converged}} ({out$iterations} iter)")
		}
		return(out)
	}
	Yarr <- if (is.null(Y_new)) object$Y else Y_new
	X <- if (is.null(X_new)) object[["X"]] else X_new
	# accept a plain matrix for a cross-sectional refit (canonical form is 3-D)
	if (is.matrix(Yarr)) {
		Yarr <- array(Yarr, dim = c(nrow(Yarr), ncol(Yarr), 1L))
	}
	if (is.matrix(X)) {
		X <- array(X, dim = c(nrow(X), ncol(X), 1L, 1L))
	}
	if (!is.null(Z_new) && is.matrix(Z_new)) {
		Z_new <- array(Z_new, dim = c(nrow(Z_new), ncol(Z_new), 1L))
	}

	prep <- .ae_prep_from_canonical(
		Yarr,
		X,
		object$W_row,
		object$W_col,
		object$mode,
		object$symmetric,
		object$row_names,
		object$col_names,
		object$time_names
	)

	# init$beta holds the dyadic coefficients only (the head of object$beta); node-covariate
	# coefficients are recovered in .ae_assemble()
	init <- list(
		mu = object$mu,
		beta = object$beta[seq_len(prep$p_dyad)],
		a = .ae_na_zero(object$a),
		b = .ae_na_zero(object$b),
		U = object$U,
		V = object[["V"]],
		L = object$L
	)
	# an empirical-Bayes count fit starts its lognormal variance at the fitted value
	if (identical(object$family, "poisson") && !is.null(object$eb$ve)) init$ve <- object$eb$ve

	lrm <- .ae_default(object$lowrank_method, "mm")
	nnm <- .ae_default(object$non_normal_method, "transform")
	lsv <- .ae_default(object$linear_solver, "eigen")
	# a normal fit carrying the empirical-Bayes record (or asked for it) is refit with the estimator
	# of ame_als_gauss.R
	norm_eb <- identical(object$family, "normal") && is.null(Z_new) &&
		(identical(bias_correct, "eb") || (is.null(bias_correct) && identical(object$bias_correct$method, "eb")))
	if (norm_eb) {
		nnm <- "irls"
		if (is.null(object$link)) object$link <- "identity"
	}

	if (!is.null(Z_new)) {
		# working-scale refit: Z_new is the Gaussian working response, so one Gaussian
		# BCD pass stands in for the family transform
		Z <- Z_new
		if (object$symmetric) for (t in seq_len(prep$Tt)) {
			Z[, , t] <- .ae_symmetrize(Z[, , t])
		}
		fit <- .ae_bcd_fit(
			Z,
			X,
			object$R,
			object$symmetric,
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
			object$symmetric,
			object$R,
			object$longitudinal,
			call = match.call(),
			Zwork = Z
		)
		out$gaussian_diagnostics <- fit$gaussian_diagnostics
		out$lowrank_method <- lrm
		out$meta$lowrank_method <- lrm
		out$non_normal_method <- nnm
		out$link <- object$link
		out$linear_solver <- lsv
		out$obs_weights <- NULL
		out <- .ae_refit_correct(
			out,
			object,
			bias_correct,
			max_iter,
			tol,
			Z = Z_new
		)
		if (verbose) {
			cli::cli_text("Warm-start refit (working scale): converged = {.val {out$converged}} ({fit$iterations} iter)")
		}
		# withhold a product the data do not determine, as the fit does (.ae_withhold_product)
		return(.ae_withhold_product(out, function(r) NA_real_))
	}

	.ae_check_als_response(prep$Yarr, object$family)
	# an empirical-Bayes fit is refit with its own estimator: every actor stays in (the penalty keeps
	# effects finite) and the actor-effect covariance starts from the fit's
	eb_refit <- identical(nnm, "irls") && (object$family %in% c("binary", "poisson") || norm_eb) &&
		(identical(bias_correct, "eb") ||
				(is.null(bias_correct) && identical(object$bias_correct$method, "eb")))
	ab_prior <- NULL
	if (eb_refit) {
		ab_prior <- .ae_eb_prior(prep, isTRUE(object$symmetric), object$family, s20 = object$eb$ve_prior[["s20"]])
		if (is.matrix(object$eb$sab) && all(dim(object$eb$sab) == dim(ab_prior$sab))) {
			ab_prior$sab <- object$eb$sab
		}
	}
	db <- NULL
	if (identical(nnm, "irls") && identical(object$family, "binary") && !eb_refit) {
		# an actor at the boundary of this data is left out on that side, as
		# in the original fit
		db <- .ae_drop_boundary(prep$Yarr, isTRUE(object$symmetric))
		if (length(db$row) || length(db$col)) prep$Yarr <- db$Y else db <- NULL
	}
	if (identical(nnm, "irls")) {
		if (eb_refit && identical(object$eb$estimator, "variational")) {
			# a variational fit is refit with its own estimator from its own variational state; the Laplace
			# evidence is formed unless lame.refit_evidence is FALSE (bootstrap replicates), and a converged
			# fit refit on its own data keeps its noise variance
			same_data <- isTRUE(object$converged) && identical(dim(Yarr), dim(object$Y)) &&
				identical(as.vector(Yarr), as.vector(object$Y)) &&
				identical(as.vector(X), as.vector(object[["X"]]))
			vb <- .ae_mf_refit(
				prep,
				object,
				ab_prior,
				object$link,
				max_iter,
				tol,
				mode_args = if (isTRUE(getOption("lame.refit_evidence", TRUE)))
					list(lowrank_method = lrm, linear_solver = lsv),
				ve_hold = same_data
			)
			ir <- list(
				fit = vb$fit,
				W = vb$W,
				Z = vb$Z,
				eta = vb$eta,
				irls_iters = vb$irls_iters,
				best_obj = vb$objective,
				stationarity = list(
					pass = vb$converged,
					stage = "variational",
					termination = vb$record$termination
				)
			)
			if (!isTRUE(vb$converged)) .ae_mf_refuse(vb$record, object$R)
		} else {
			# refit an irls fit by re-running the reweighting loop, warm-started
			ir <- .ae_irls_run(
				prep,
				object$family,
				object$link,
				object$R,
				object$symmetric,
				max_iter,
				tol,
				lowrank_method = lrm,
				init = init,
				linear_solver = lsv,
				ab_prior = ab_prior
			)
			if (norm_eb && eb_refit) {
				# a normal fit without the empirical-Bayes estimator gets all of it: the variational stage (with
				# the dyadic correlation of a directed network) and the evidence
				gp <- if (!isTRUE(object$symmetric) && identical(object$mode, "unipartite")) list(gauss = TRUE)
				ir <- .ae_mf_stage(
					prep,
					ir,
					ab_prior,
					object$R,
					isTRUE(object$symmetric),
					object$link,
					max_iter,
					tol,
					pair = gp
				)
				if (!isTRUE(ir$selection_converged)) .ae_mf_refuse(ir$fit$eb_mf$record, object$R)
				ir <- .ae_gauss_evidence(
					prep,
					ir,
					ab_prior,
					object$R,
					isTRUE(object$symmetric),
					object$link,
					max_iter,
					tol,
					lrm,
					lsv
				)
			}
		}
		fit <- ir$fit
		out <- .ae_assemble(
			fit,
			prep,
			object$family,
			object$mode,
			object$symmetric,
			object$R,
			object$longitudinal,
			call = match.call(),
			Zwork = ir$Z,
			link = object$link,
			Wobs = ir$W
		)
		out$objective <- ir$best_obj
		out$irls_iterations <- ir$irls_iters
		out$irls_diagnostics <- ir$stationarity
		if (eb_refit) out$iterations <- ir$irls_iters
	} else {
		Z <- .ae_working_response(Yarr, object$family)
		if (object$symmetric) for (t in seq_len(prep$Tt)) {
			Z[, , t] <- .ae_symmetrize(Z[, , t])
		}
		fit <- .ae_bcd_fit(
			Z,
			X,
			object$R,
			object$symmetric,
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
			object$symmetric,
			object$R,
			object$longitudinal,
			call = match.call()
		)
	}

	out$gaussian_diagnostics <- fit$gaussian_diagnostics
	out$lowrank_method <- lrm
	out$meta$lowrank_method <- lrm
	out$non_normal_method <- nnm
	out$link <- object$link
	out$linear_solver <- lsv
	out$obs_weights <- if (identical(nnm, "irls")) ir$W else NULL
	out <- .ae_refit_correct(out, object, bias_correct, max_iter, tol)
	db <- .ae_refit_boundary(
		db,
		object$not_estimated,
		out$row_names,
		out$col_names,
		isTRUE(object$symmetric)
	)
	out <- .ae_mark_not_estimated(out, db)

	if (verbose) {
		cli::cli_text("Warm-start refit: converged = {.val {out$converged}} ({fit$iterations} iter)")
	}
	# withhold a product the data do not determine, as the fit does (.ae_withhold_product); the
	# evidence of the smaller ranks is not refit here
	.ae_withhold_product(out, function(r) NA_real_)
}

# bias correction of a static refit: the one `object` carries unless the caller names one; a
# requested jackknife that cannot be refit errors, so a bootstrap replicate is dropped and counted
.ae_refit_correct <- function(
	out,
	object,
	bias_correct,
	max_iter,
	tol,
	Z = NULL
) {
	correction <- if (is.null(bias_correct)) {
		.ae_default(object$bias_correct$method, "none")
	} else {
		bias_correct
	}
	if (identical(correction, "eb")) {
		# the empirical-bayes refit is an estimator, not a correction
		out$bias_correct <- list(
			method = "eb",
			applied = FALSE,
			requested = if (identical(bias_correct, "eb")) "eb",
			default = is.null(bias_correct) &&
				isTRUE(object$bias_correct$default)
		)
		out
	} else if (identical(correction, "jackknife")) {
		.ae_jackknife(out, max_iter, tol, Z = Z, explicit = TRUE)
	} else if (identical(correction, "bootstrap")) {
		# a refit asked for the correction runs it; refits of a corrected fit (warm-start refits,
		# bootstrap replicates) drop its stored correction
		if (identical(bias_correct, "bootstrap") ||
				is.null(object$bias_correct$correction)) {
			.ae_boot_bias(out, max_iter, tol)
		} else {
			.ae_bb_apply_stored(out, object$bias_correct)
		}
	} else {
		# a refit of an uncorrected fit keeps its record (size basis or opt-out); one asked to drop the
		# correction records the request
		out$bias_correct <- if (is.null(bias_correct) &&
				identical(object$bias_correct$method, "none")) {
			object$bias_correct
		} else if (identical(bias_correct, "none")) {
			list(method = "none", applied = FALSE, requested = "none")
		} else {
			list(method = "none", applied = FALSE)
		}
		out
	}
}

# bootstrap

# draw on the stated response scale, with independent cells in bipartite networks even when the
# two actor sets have equal size
.ae_boot_draw_response <- function(object, ez, rho, s2, working = FALSE) {
	nr <- nrow(ez)
	nc <- ncol(ez)
	draw_normal <- function(eta, variance) {
		if (identical(object$mode, "bipartite")) {
			eta + matrix(stats::rnorm(nr * nc, sd = sqrt(variance)), nr, nc)
		} else simZ(eta, rho, variance)
	}
	if (working || identical(object$family, "normal")) {
		return(draw_normal(ez, s2))
	}
	if (identical(object$family, "poisson")) {
		z_lat <- if (is.finite(s2) && s2 > 0) draw_normal(ez, s2) else ez
		return(matrix(stats::rpois(nr * nc, .als_pois_rate(z_lat)), nr, nc))
	}
	eta <- if (identical(object$link, "logit")) {
		stats::qnorm(stats::plogis(ez, log.p = TRUE), log.p = TRUE)
	} else ez
	1 * (draw_normal(eta, 1) > 0)
}

# replicate accounting, single-replicate refits and checks on the replicate distribution

# fewest replicates: `bootstrap = N` needs .ae_boot_min_requested so a 95% percentile interval has
# tails, and ame_als_bootstrap() needs .ae_boot_min_valid valid replicates
.ae_boot_min_requested <- 30L
.ae_boot_min_valid <- 10L

# share of the requested replicates that may be dropped before the bootstrap is refused; 5 percent
# keeps a 95% percentile interval at 90% or more in the worst case
.ae_boot_max_fail <- 0.05

# separation of adjacent singular values (eigenvalues for symmetric fits), in bootstrap sds of the gap,
# below which the columns rotate between replicates and get no column-wise standard error
.ae_boot_sv_z <- 5

# a bootstrap of fewer than .ae_boot_min_requested replicates is refused
.ae_check_bootstrap_n <- function(n) {
	min_req <- .ae_boot_min_requested
	if (n > 0L && n < min_req) {
		cli::cli_abort(
			c(
				"{.arg bootstrap} = {n} is below the minimum of {.val {min_req}} replicates.",
				"x" = "A 95% percentile interval from fewer replicates is not returned: its end points are the smallest and largest replicate.",
				"i" = "Use {.code bootstrap = 200}, any value of at least {.val {min_req}}, or {.code bootstrap = 0} for the point fit."
			),
			class = "lame_als_bootstrap_replicates"
		)
	}
	invisible(n)
}

# objective that compares two refits of one replicate: the penalised IRLS deviance, or the Gaussian
# residual sum of squares (uncorrected when a split-sample correction rebuilt the fit)
.ae_boot_objective <- function(fit) {
	obj <- fit$objective
	if (is.null(obj) || !is.finite(obj)) obj <- fit$raw$deviance
	if (is.null(obj) || !is.finite(obj)) obj <- fit$deviance
	if (is.null(obj)) NA_real_ else as.numeric(obj)
}

# distance between two refits of one replicate in deviance units: squared difference of fitted
# predictors weighted by the mean IRLS weight, summed over observed cells (as .ae_irls_distance)
.ae_boot_distance <- function(f, ref) {
	ez_f <- tryCatch(array(unlist(f$EZ), dim(f$Y)), error = function(e) NULL)
	ez_r <- tryCatch(array(unlist(ref$EZ), dim(ref$Y)), error = function(e) NULL)
	if (is.null(ez_f) || is.null(ez_r)) return(Inf)
	w <- if (!is.null(f$obs_weights) && !is.null(ref$obs_weights)) {
		(f$obs_weights + ref$obs_weights) / 2
	} else 1
	sum(w * (ez_f - ez_r)^2, na.rm = TRUE)
}

# what the refit of a replicate checks the point estimate's optimum against
.ae_boot_optimum_check <- function(object, dyn) {
	if (isTRUE(dyn)) return("none: a dynamic fit is refit by its own fitter")
	if (identical(object$bias_correct$method, "eb")) {
		return("none: an empirical-Bayes refit starts from the variational state of the fit")
	}
	if (object$R < 1L) return("none: without a multiplicative term the optimum is unique")
	k <- switch(
		as.character(.ae_default(object$multistart, "none")),
		cheap = 4L,
		full = 8L,
		1L
	)
	if (k > 1L) {
		sprintf("warm start against the cold start and %d random starts of the original multistart", k - 1L)
	} else "warm start against the cold start of the original fit"
}

# refit one replicate: warm from the point estimate (the budget, then four times it), checked against
# the original fit's cold and multistart starts; a check start wins only with a lower objective and a
# fit more than a deviance unit away. returns list(fit, attempts, iterations, start, error, last_fit)
.ae_boot_refit <- function(object, Y_b, X_b, Z_b, dyn, max_iter, tol) {
	family <- object$family
	Rlat <- object$R
	gaussian <- identical(family, "normal") || !is.null(Z_b)
	eb <- identical(object$bias_correct$method, "eb")
	check <- !dyn && Rlat > 0L && !eb
	# replicates are summarised, not reported, so a withheld product is restored to its variational means
	refit_with <- function(seed_object, budget) .ae_product_restore(refit_raw(seed_object, budget))
	refit_raw <- function(seed_object, budget) {
		suppressWarnings(if (dyn) {
			.ae_dynamic_refit(
				seed_object,
				Y_new = Y_b,
				X_new = X_b,
				max_iter = budget,
				tol = tol,
				verbose = FALSE
			)
		} else if (!is.null(Z_b)) {
			ame_als_refit(
				seed_object,
				Y_new = object$Y,
				Z_new = Z_b,
				max_iter = budget,
				tol = tol,
				verbose = FALSE
			)
		} else {
			ame_als_refit(
				seed_object,
				Y_new = Y_b,
				X_new = X_b,
				max_iter = budget,
				tol = tol,
				verbose = FALSE
			)
		})
	}
	attempts <- 0L
	last_budget <- as.integer(max_iter)
	first_error <- NULL
	last_fit <- NULL
	# check starts run at the budget; the four-times budget is for the start a replicate depends on
	run_start <- function(seed_object, escalate = TRUE) {
		for (budget in if (escalate) c(max_iter, 4L * max_iter) else max_iter) {
			attempts <<- attempts + 1L
			last_budget <<- as.integer(budget)
			fit <- tryCatch(refit_with(seed_object, budget), error = function(e) e)
			if (inherits(fit, "error")) {
				if (is.null(first_error)) first_error <<- fit
				# an unconverged split-sample half may finish with more budget
				if (!inherits(fit, "lame_jk_nonconvergence")) return(NULL)
				next
			}
			last_fit <<- fit
			if (isTRUE(fit$converged)) {
				attr(fit, "boot_budget") <- as.integer(budget)
				return(fit)
			}
		}
		NULL
	}
	found <- list()
	warm <- run_start(object)
	if (!is.null(warm)) found$warm <- warm
	if (!dyn && (check || (is.null(warm) && gaussian))) {
		cold <- object
		cold$mu <- 0
		cold$beta[] <- 0
		cold$a[] <- 0
		cold$b[] <- 0
		cold$U <- cold$V <- cold$L <- NULL
		f <- run_start(cold, escalate = is.null(warm))
		if (!is.null(f)) found$cold <- f
		# random starts of the original multistart, or seven for retrying an unfinished Gaussian refit
		k_orig <- switch(
			as.character(.ae_default(object$multistart, "none")),
			cheap = 4L,
			full = 8L,
			1L
		)
		n_random <- if (check && k_orig > 1L) {
			k_orig - 1L
		} else if (!length(found) && Rlat > 0L && gaussian) 7L else 0L
		if (n_random > 0L) {
			Zsize <- if (!is.null(Z_b)) Z_b else .ae_working_response(Y_b, family)
			starts <- .ae_multistart_starts(
				Zsize,
				Rlat,
				isTRUE(object$symmetric),
				n_random + 1L,
				6886L
			)
			for (j in 2:(n_random + 1L)) {
				alternate <- cold
				alternate$U <- starts[[j]]$U
				alternate$V <- starts[[j]]$V
				alternate$L <- starts[[j]]$L
				f <- run_start(alternate, escalate = !length(found))
				if (!is.null(f)) found[[paste0("random", j - 1L)]] <- f
			}
		}
	}
	if (!length(found)) {
		return(list(
			fit = NULL,
			attempts = attempts,
			iterations = last_budget,
			start = NA_character_,
			error = first_error,
			last_fit = last_fit
		))
	}
	ref_label <- if ("warm" %in% names(found)) "warm" else names(found)[1L]
	best <- found[[ref_label]]
	best_label <- ref_label
	best_obj <- .ae_boot_objective(best)
	for (nm in setdiff(names(found), ref_label)) {
		f <- found[[nm]]
		o <- .ae_boot_objective(f)
		lower <- is.finite(o) &&
			(!is.finite(best_obj) || o < best_obj - 1e-6 * max(1, abs(best_obj)))
		other <- gaussian || .ae_boot_distance(f, best) > 1
		if (lower && other) {
			best <- f
			best_label <- nm
			best_obj <- o
		}
	}
	list(
		fit = best,
		attempts = attempts,
		iterations = .ae_default(attr(best, "boot_budget"), last_budget),
		start = best_label,
		error = first_error,
		last_fit = last_fit
	)
}

# refusal of a bootstrap whose dropped replicates exceed the cap
.ae_boot_abort_failed <- function(
	fail_reason,
	R,
	max_fail,
	first_failed,
	max_iter,
	attempts,
	b
) {
	n_fail <- sum(nzchar(fail_reason))
	n_nc <- sum(fail_reason == "nonconverged")
	n_sp <- sum(fail_reason == "split")
	n_er <- sum(fail_reason %in% c("error", "nonfinite"))
	share <- round(100 * .ae_boot_max_fail)
	budget <- 4L * as.integer(max_iter)
	cli::cli_abort(
		c(
			"More than {share}% of the bootstrap replicates could not be refit: {n_fail} of the first {b} (of {R}) were dropped, and at most {max_fail} may be.",
			"x" = "{n_nc} did not converge within {.arg max_iter} = {max_iter} iterations nor within {budget}; {n_sp} could not carry the split-sample correction; {n_er} failed to refit or returned non-finite values.",
			"i" = "Increase {.arg max_iter}, or simplify the model. No interval is returned from the remaining replicates: dropping many replicates keeps the ones that are easiest to fit."
		),
		class = "lame_als_bootstrap_convergence",
		replicate = first_failed$replicate,
		fit = first_failed$fit,
		n_failed = n_fail,
		n_total = R,
		max_failed = max_fail,
		reasons = c(nonconverged = n_nc, split = n_sp, error = n_er)
	)
}

# singular values of the multiplicative term, largest first (signed eigenvalues for symmetric
# fits), from the r x r product (U'U)(V'V)
.ae_boot_sv <- function(fit) {
	if (isTRUE(fit$symmetric)) return(sort(as.numeric(fit$L), decreasing = TRUE))
	m <- crossprod(fit$U) %*% crossprod(fit[["V"]])
	ev <- Re(eigen(m, only.values = TRUE)$values)
	sort(sqrt(pmax(ev, 0)), decreasing = TRUE)
}

# columns whose singular value is separated from both neighbours by .ae_boot_sv_z bootstrap sds of the gap
.ae_boot_separation <- function(object, sv_rep) {
	Rlat <- ncol(sv_rep)
	sv0 <- .ae_boot_sv(object)
	# column of each sorted singular value
	idx <- if (isTRUE(object$symmetric)) order(as.numeric(object$L), decreasing = TRUE) else seq_len(Rlat)
	resolved <- rep(TRUE, Rlat)
	gap <- z <- sd_gap <- rep(NA_real_, max(Rlat - 1L, 0L))
	for (k in seq_len(max(Rlat - 1L, 0L))) {
		gap[k] <- sv0[k] - sv0[k + 1L]
		sd_gap[k] <- stats::sd(sv_rep[, k] - sv_rep[, k + 1L], na.rm = TRUE)
		z[k] <- if (is.finite(sd_gap[k]) && sd_gap[k] > 0) gap[k] / sd_gap[k] else
		if (is.finite(gap[k]) && gap[k] > 0) Inf else 0
		if (!(z[k] >= .ae_boot_sv_z)) resolved[idx[c(k, k + 1L)]] <- FALSE
	}
	list(
		threshold = .ae_boot_sv_z,
		values = sv0,
		gap = gap,
		sd_gap = sd_gap,
		z = z,
		resolved = resolved
	)
}

# principal angles, in degrees, between the column spaces of `a` and `b` (both n x r): all r
# of them, 90 for a direction one of the two lacks
.ae_boot_angles <- function(a, b) {
	r <- ncol(b)
	qa <- qr(a); qb <- qr(b)
	ra <- min(qa$rank, r); rb <- min(qb$rank, r)
	out <- rep(90, r)
	if (ra < 1L || rb < 1L) return(out)
	sv <- svd(
		crossprod(
			qr.Q(qa)[, seq_len(ra), drop = FALSE],
			qr.Q(qb)[, seq_len(rb), drop = FALSE]
		),
		nu = 0,
		nv = 0
	)$d
	ang <- acos(pmin(pmax(sv, 0), 1)) * 180 / pi
	out[seq_along(ang)] <- sort(ang)
	out
}

# subspace uncertainty, invariant to column rotation: principal angles (degrees) of each
# replicate's U and V spaces to the point estimate's
.ae_boot_subspace <- function(object, U_raw, V_raw) {
	Rlat <- ncol(object$U)
	angles <- function(raw, point) {
		if (!nrow(raw)) return(matrix(NA_real_, 0L, Rlat))
		t(apply(raw, 1L, function(r) .ae_boot_angles(matrix(r, nrow(point), Rlat), point)))
	}
	au <- angles(U_raw, object$U)
	av <- if (isTRUE(object$symmetric)) au else angles(V_raw, object[["V"]])
	largest <- function(a) c(
		median = stats::median(a[, Rlat]),
		q95 = unname(stats::quantile(a[, Rlat], 0.95))
	)
	list(
		unit = "degrees",
		U = au,
		V = av,
		largest = rbind(U = largest(au), V = largest(av))
	)
}

# the point estimate must lie inside the central 95% of its replicates: a coefficient outside is an
# error carrying the bootstrap, a variance component outside has its interval withheld (vc_uncentred)
.ae_boot_check_centred <- function(out) {
	slack <- function(x) 1e-8 * pmax(1, abs(x))
	outside <- function(est, lo, hi) {
		which(is.finite(est) & is.finite(lo) & is.finite(hi) &
				(est < lo - slack(est) | est > hi + slack(est)))
	}
	ib <- outside(out$point_est, out$ci_lo, out$ci_hi)
	if (length(ib)) {
		items <- sprintf(
			"%s (estimate %s, replicate 2.5%% to 97.5%%: %s to %s)",
			out$param_names[ib],
			format(signif(out$point_est[ib], 4)),
			format(signif(out$ci_lo[ib], 4)),
			format(signif(out$ci_hi[ib], 4))
		)
		cli::cli_abort(
			c(
				"The point estimate lies outside the 95% interval of its own bootstrap replicates for {length(items)} coefficient{?s}.",
				stats::setNames(items, rep("x", length(items))),
				"i" = "The replicates are not distributed around the fitted optimum: a replicate refit landed at another local optimum, or the estimator is biased here by more than about two standard errors.",
				"i" = "Refit the point estimate with {.code multistart = \"full\"} so that it is the best of several starts. The bootstrap is attached to this error as {.code boot}."
			),
			class = "lame_als_bootstrap_optimum",
			boot = out
		)
	}
	iv <- outside(out$vc_point_est, out$vc_ci_lo, out$vc_ci_hi)
	if (length(iv)) {
		out$vc_uncentred <- data.frame(
			component = names(out$vc_point_est)[iv],
			estimate = unname(out$vc_point_est[iv]),
			lo = unname(out$vc_ci_lo[iv]),
			hi = unname(out$vc_ci_hi[iv]),
			stringsAsFactors = FALSE
		)
		out$vc_se[iv] <- NA_real_
		out$vc_ci_lo[iv] <- NA_real_
		out$vc_ci_hi[iv] <- NA_real_
	}
	out
}

# quantile levels of a bootstrap interval

# probabilities of the end points of a `level` interval: nominal for a parametric bootstrap, and the
# expanded levels Phi(sqrt(B / (B - 1)) * qt(alpha / 2, B - 1)) for B blocks (Hesterberg 2015)
.ae_level_probs <- function(level, type, n_blocks) {
	a <- (1 - level) / 2
	if (identical(type, "block") && is.numeric(n_blocks) && length(n_blocks) == 1L &&
			is.finite(n_blocks) && n_blocks > 2) {
		a <- stats::pnorm(sqrt(n_blocks / (n_blocks - 1)) * stats::qt(a, n_blocks - 1))
	}
	c(a, 1 - a)
}

# the same for a bootstrap object
.ae_boot_probs <- function(x, level = 0.95) {
	.ae_level_probs(level, x$type, x$n_blocks)
}

# parameters with no standard error or interval

# additive paths of a fit with time-varying effects and a sampler of replicate paths (actor level
# plus a stationary AR(1) deviation), so replicates carry the variability of the effects
.ae_boot_ab_paths <- function(object) {
	Tt <- object$n_time
	nr <- object$dims$n_row
	nc <- object$dims$n_col
	orig <- lapply(seq_len(Tt), function(t) {
		outer(as.numeric(object$a_dynamic[, t]), as.numeric(object$b_dynamic[, t]), "+")
	})
	prof <- list(a = .ae_path_lag_profile(object$a_dynamic))
	if (!object$symmetric) prof$b <- .ae_path_lag_profile(object$b_dynamic)
	rho_w <- .ae_path_decay_rho(prof)
	if (!is.finite(rho_w)) rho_w <- .ae_default(object$rho_ab, 0.8)
	lev <- lapply(prof, function(pf) .ae_path_level_var(pf, rho_w))
	dev_pars <- list(
		rho = rho_w,
		tau2_a = max(.ae_default(lev$a$tau2, 0), 0, na.rm = TRUE),
		tau2_b = if (object$symmetric) 0 else max(.ae_default(lev$b$tau2, 0), 0, na.rm = TRUE)
	)
	# one replicate's additive paths: level plus stationary AR(1) deviation
	draw <- function() {
		ar_path <- function(n, tau2) {
			out <- matrix(0, n, Tt)
			if (tau2 <= 0) return(out)
			sd_stat <- sqrt(tau2)
			out[, 1L] <- stats::rnorm(n, 0, sd_stat)
			if (Tt > 1L) for (t in 2:Tt) {
				out[, t] <- dev_pars$rho * out[, t - 1L] +
					stats::rnorm(n, 0, sd_stat * sqrt(1 - dev_pars$rho^2))
			}
			out
		}
		va_l <- max(.ae_default(as.numeric(object$VC["va"]), 0), 0, na.rm = TRUE)
		vb_l <- max(.ae_default(as.numeric(object$VC["vb"]), 0), 0, na.rm = TRUE)
		a_new <- stats::rnorm(nr, 0, sqrt(va_l)) + ar_path(nr, dev_pars$tau2_a)
		b_new <- if (object$symmetric) a_new else {
			stats::rnorm(nc, 0, sqrt(vb_l)) + ar_path(nc, dev_pars$tau2_b)
		}
		lapply(seq_len(Tt), function(t) outer(a_new[, t], b_new[, t], "+"))
	}
	list(orig = orig, draw = draw)
}

# interval of a variance (va, vb) from replicate quantiles q: the basic interval on the log scale
# (est^2 / q) for fixed effects, 2 est - q for empirical Bayes; percentile ends when a quantile is <= 0
.ae_vc_interval <- function(q, est, eb = FALSE) {
	if (!is.finite(est) || est <= 0 || length(q) != 2L || any(!is.finite(q))) return(q)
	if (eb) return(pmax(c(2 * est - q[2L], 2 * est - q[1L]), 0))
	if (any(q <= 0)) return(q)
	c(est^2 / q[2L], est^2 / q[1L])
}

# kind of each entry of fit$coefficients: intercept, dyadic, then node-covariate (sender, receiver)
.ae_coef_kinds <- function(object) {
	p <- length(object$coefficients)
	pd <- dim(object[["X"]])[3L]
	pd <- if (length(pd) == 1L && is.finite(pd)) min(as.integer(pd), p - 1L) else 0L
	c("intercept", rep("dyad", pd), rep("node", p - 1L - pd))
}

# parameters the bootstrap of scheme `type` gives no standard error or interval for, one row each
# with the reason; NULL when it covers all (?ame_als_bootstrap, "Block bootstrap")
.ae_boot_withheld <- function(object, type) {
	nm <- names(object$coefficients)
	kinds <- .ae_coef_kinds(object)
	eb <- identical(object$bias_correct$method, "eb")
	rows <- list()
	add <- function(parameter, kind, reason) {
		if (length(parameter)) {
			rows[[length(rows) + 1L]] <<- data.frame(
				parameter = parameter,
				kind = kind,
				reason = reason,
				stringsAsFactors = FALSE
			)
		}
	}
	if (identical(type, "block")) {
		# a block resample redraws time slices and holds actors fixed, so it misses the variability over actors
		add(
			nm[kinds == "node"],
			"coefficient",
			"slice resampling holds the actors fixed, and the sampling variability of a node-covariate coefficient is that over actors"
		)
		vc <- c("va", "cab", "vb")
		add(
			vc[vc %in% names(object$VC) & is.finite(object$VC[vc])],
			"vc",
			"slice resampling holds the actors fixed, and the sampling variability of the variances of the actor effects is that over actors"
		)
		# a resample repeats some slices and omits others, so the refit absorbs more noise and the
		# replicate residual variances sit below the estimate
		if (!identical(object$family, "binary") && is.finite(object$VC["ve"])) {
			add(
				"ve", "vc",
				"the refit to resampled slices, which repeat some and omit others, absorbs more of the noise than the original fit, so the replicate residual variances sit below the estimate"
			)
		}
		if (eb) {
			add(
				"intercept", "coefficient",
				"the intercept of this estimator is the population intercept, whose sampling variability includes that over actors, which slice resampling holds fixed"
			)
			# empirical-Bayes effects (and factors) are shrunk toward their prior mean, and the replicates
			# vary around the shrunk estimates
			shrunk <- c(
				"a",
				if (!isTRUE(object$symmetric)) "b",
				if (object$R > 0L && is.matrix(object$U)) c("U", if (!isTRUE(object$symmetric)) "V")
			)
			add(
				shrunk,
				"effects",
				"the effects of this estimator are shrunk toward their prior mean, and the replicates of a resample vary around the shrunk estimates, not around the effects"
			)
		}
	}
	if (.ae_product_is_withheld(object)) {
		add(
			c("U", if (!isTRUE(object$symmetric)) "V"),
			"effects",
			"the latent product is not identified by the data at this rank and is withheld (fit$eb$product_withheld)"
		)
	}
	if (identical(object$family, "poisson") && identical(object$non_normal_method, "irls") &&
			!eb) {
		add(
			"intercept", "coefficient",
			"the fixed-effects Poisson intercept carries a finite-sample bias that every replicate reproduces, so the interval is displaced from the truth by about twice the bias"
		)
	}
	if (!length(rows)) return(NULL)
	out <- do.call(rbind, rows)
	out[!duplicated(out[c("parameter", "kind")]), , drop = FALSE]
}

# names of the parameters of `kind` ("coefficient" or "vc") with no bootstrap interval
.ae_boot_withheld_names <- function(x, kind) {
	w <- x$withheld
	if (is.null(w) || !nrow(w)) return(character(0))
	w$parameter[w$kind == kind]
}

# set the replicate summaries of withheld parameters to NA (the replicates stay in the object)
.ae_boot_apply_withheld <- function(out, withheld) {
	if (is.null(withheld) || !nrow(withheld)) return(out)
	out$withheld <- withheld
	w <- .ae_boot_withheld_names(out, "coefficient")
	out$se[w] <- NA_real_
	out$ci_lo[w] <- NA_real_
	out$ci_hi[w] <- NA_real_
	v <- .ae_boot_withheld_names(out, "vc")
	out$vc_se[v] <- NA_real_
	out$vc_ci_lo[v] <- NA_real_
	out$vc_ci_hi[v] <- NA_real_
	e <- .ae_boot_withheld_names(out, "effects")
	if ("a" %in% e) out$se_a[] <- NA_real_
	if ("b" %in% e) out$se_b[] <- NA_real_
	if ("U" %in% e && !is.null(out$se_U)) out$se_U[] <- out$se_U_raw[] <- NA_real_
	if ("V" %in% e && !is.null(out$se_V)) out$se_V[] <- out$se_V_raw[] <- NA_real_
	if (any(c("U", "V") %in% e) && !is.null(out$se_L)) {
		out$se_L[] <- out$ci_L_lo[] <- out$ci_L_hi[] <- NA_real_
	}
	out
}

# one sentence per reason on the parameters without an interval, for print() and summary()
.ae_boot_withheld_lines <- function(x) {
	w <- x$withheld
	if (is.null(w) || !nrow(w)) return(character(0))
	by <- split(w$parameter, w$reason)
	vapply(names(by), function(r) {
		paste0(
			"No standard error or interval for ",
			paste(by[[r]], collapse = ", "),
			": ",
			r,
			"."
		)
	}, character(1L), USE.NAMES = FALSE)
}

# refuse the bootstrap of an uncorrected binary fixed-effects fit: its replicates reproduce the
# incidental-parameter bias of the dyadic coefficient
.ae_boot_refuse_uncorrected <- function() {
	cli::cli_abort(
		c(
			"A bootstrap is not available for an uncorrected binary fixed-effects fit.",
			"x" = "The fit is biased by the incidental-parameters problem and every replicate reproduces the bias, so the interval is centred away from the truth.",
			"i" = "Bootstrap the default empirical-Bayes fit (omit {.arg bias_correct}), or a fit corrected with {.code bias_correct = \"jackknife\"} or {.code \"bootstrap\"}, or use MCMC with {.fn ame} / {.fn lame}."
		),
		class = c("lame_als_bootstrap_bias", "lame_als_bootstrap_unsupported")
	)
}

# the same refusal from the ame_als() / lame_als() arguments, before anything is fitted
.ae_check_bootstrap_bias <- function(bootstrap, family, non_normal_method, bias_correct) {
	if (bootstrap > 0L && is.character(family) && length(family) == 1L &&
			family %in% c("binary", "bin") && identical(non_normal_method, "irls") &&
			identical(bias_correct, "none")) {
		.ae_boot_refuse_uncorrected()
	}
	invisible(bootstrap)
}

#' Bootstrap uncertainty for the fast AME estimator
#'
#' @description
#' Computes bootstrap standard errors and percentile confidence intervals for
#' a fast AME fit produced by \code{\link{ame_als}} or
#' \code{\link{lame_als}}. Two strategies are available:
#' \describe{
#'   \item{\code{parametric}}{(default) Simulates fresh outcomes from the
#'     fitted model and refits. For \code{normal} and IRLS (\code{binary},
#'     \code{poisson}) fits the simulation is on the calibrated response scale;
#'     for a non-normal \code{transform} fit -- whose estimator is a Gaussian
#'     fit to a fixed transformed response -- it is on the Gaussian working
#'     scale. The fitted residual variance is df-corrected (mean-model degrees
#'     of freedom) at fit time, so simulating with it keeps the replicates
#'     centred on the point estimate.
#'     Available for both cross-sectional and longitudinal fits.}
#'   \item{\code{block}}{Resamples time slices with replacement. With
#'     \code{block_length = 1} slices are resampled independently;
#'     \code{block_length > 1} draws contiguous blocks (a moving-block
#'     bootstrap), which preserves short-run temporal dependence. Longitudinal
#'     fits only, and only when the panel supplies at least
#'     \code{ceiling(T / block_length) >= 10} resampling blocks: below that the
#'     ends of the interval are the extreme replicates (a plain percentile
#'     interval covered the dyadic coefficient in 0.83 of the panels at five
#'     blocks), so the call warns and
#'     uses \code{parametric} instead. The returned object records the
#'     scheme actually used in \code{type} and the one asked for in
#'     \code{requested_type}. The interval end points are read at expanded
#'     percentile levels, and the scheme returns no standard error or interval
#'     for the parameters that vary over actors and not over time; see
#'     \emph{Block bootstrap} below.}
#' }
#'
#' Every replicate is refit by \code{\link{ame_als_refit}}, warm-started
#' from the original point estimate so that replicates do not drift to
#' different local optima. An empirical-Bayes fit
#' (\code{fit$bias_correct$method == "eb"}, the default for normal, binary and Poisson
#' fits; see \code{\link{ame_als}}) has every replicate refit with the same estimator,
#' re-estimating the covariance of the actor effects from the fitted one, and
#' draws the replicate effects uncentred, since the fitted intercept is the
#' population intercept; the intervals then describe the empirical-Bayes
#' estimator, and agree with the Laplace standard errors of
#' \code{\link{vcov.ame_als}} where those are calibrated. A fit carrying the
#' split-sample jackknife
#' (\code{fit$bias_correct$method == "jackknife"}) has every replicate
#' corrected the same
#' way, so the intervals describe the corrected estimator; a replicate whose
#' split cannot be refit is dropped and counted like any failed replicate. A
#' fit carrying the parametric-bootstrap correction
#' (\code{fit$bias_correct$method == "bootstrap"}) has the bias it estimated
#' removed from every replicate rather than re-estimated, so a replicate costs
#' one refit; the bias moves slowly with the parameters, and the replicate
#' spread is that of the corrected estimator to first order. A binary
#' replicate in which an actor has a single outcome value on a side records
#' that actor's effect there as \code{-Inf} (no tie) or \code{Inf} (only
#' ties), its maximum-likelihood value, so the actor's interval is open on
#' that side and its standard error is \code{Inf}; an actor the original fit
#' leaves out is \code{NA} throughout. A dynamic ALS fit is refit dynamically, on the
#' parametric route only: a block resample reorders the time axis the
#' temporal penalties are defined along. Only a penalised normal-family dynamic
#' fit (\code{bias_correct = "none"}) has a bootstrap, and its time-varying
#' factors are not summarised. Every replicate redraws the additive paths (an
#' actor level with the fitted level variance plus an AR(1) deviation), so the
#' standard errors carry the variability of the actors: the intercept's
#' interval is for the population intercept (it covered 0.90 and 0.82 of 100
#' panels of 20 actors and 8 periods and 40 of 30 and 12, against 0.08 and
#' 0.07 when the fitted paths were held fixed, which covers the realized
#' intercept only) and the dyadic coefficient's covered 0.91 and 0.95. The
#' bootstrap of \code{\link{lame}} (\code{bootstrap = N}) draws its replicates
#' the same way. The
#' default dynamic empirical-Bayes fit of every family refuses a bootstrap, with
#' an error of class \code{"lame_dynamic_eb_unsupported"}, and carries standard
#' errors in \code{vcov()} and \code{confint()} instead (Laplace ones for a
#' binary or Poisson fit, the exact Gaussian marginal ones for a normal fit). A
#' penalised binary or Poisson dynamic fit refuses it with an error of class
#' \code{"lame_als_bootstrap_unsupported"}. These refusals come before the
#' replicate count is checked.
#' Standard errors are replicate column standard deviations; confidence
#' intervals use the percentile method.
#'
#' \strong{Replicates.} At least 10 valid replicates are needed
#' (\code{R} below 10 is an error of class
#' \code{"lame_als_bootstrap_replicates"}; \code{bootstrap = N} in
#' \code{\link{ame_als}} and \code{\link{lame_als}} needs 30, and a 95 percent
#' percentile interval is best read from 200 or more). A replicate whose refit has
#' not converged is retried on the same data with four times the budget and,
#' failing that, dropped: it never enters an interval. A replicate that errors,
#' returns non-finite values or cannot carry the split-sample correction is
#' dropped too, and every dropped replicate is counted
#' (\code{n_failed}, with the reason in \code{fail_reason}). A dropped replicate
#' is usually an extreme one, so the nominal level of a percentile interval can
#' move by at most the share dropped; more than 5 percent of the \code{R}
#' replicates dropped is an error of class
#' \code{"lame_als_bootstrap_convergence"}, raised at the replicate that crosses
#' the cap, with no interval from the rest.
#'
#' \strong{One optimum per replicate.} The replicates of a fit with a
#' multiplicative term (\code{R > 0}) are refit from the point estimate and
#' also from the start the original fit used (the cold spectral start and, for a
#' fit with \code{multistart}, its random starts); the converged refit with the
#' lowest objective is kept, so a replicate cannot sit at a worse local optimum
#' than the original procedure finds on its data. The start each replicate kept
#' is in \code{boot_refit_start}, and \code{optimum_check} says what was
#' compared. Without a multiplicative term the optimum is unique and a replicate
#' is refit once. (An empirical-Bayes refit starts from the variational state of
#' the fit and is not compared with other starts.) The point estimate of every
#' coefficient must lie inside the central 95 percent of its own replicates: if
#' it does not, the replicates are distributed around another optimum, or the
#' estimator is biased by more than two standard errors, and the call is an
#' error of class \code{"lame_als_bootstrap_optimum"} that carries the
#' bootstrap as its \code{boot} field. A variance component whose estimate lies
#' outside its replicate interval is the bias of its estimator in a regime a
#' percentile interval cannot describe (the variance of the lognormal layer of a
#' count network without extra-Poisson noise is the usual case): its standard
#' error and interval are \code{NA}, it is named in \code{vc_uncentred} with the
#' estimate and the replicate range, and the rest of the bootstrap stands.
#'
#' \strong{Block bootstrap.} A statistic that varies over \eqn{B} resampling
#' blocks has a percentile interval that falls short of its level: in the
#' simulations below a 95 percent interval covered the dyadic coefficient
#' in 0.90 of the panels over ten to thirty blocks. The end points of a
#' block interval are therefore read at the expanded percentile levels
#' \eqn{\Phi(\sqrt{B / (B - 1)}\, t_{B - 1, \alpha / 2})}{Phi(sqrt(B / (B - 1)) * qt(alpha / 2, B - 1))}
#' (Hesterberg, 2015): 0.9 and 99.1 percent for a 95 percent interval at ten
#' blocks, 1.6 and 98.4 percent at twenty. \code{n_blocks} records \eqn{B}, and
#' \code{confint()}, \code{tidy()} and the interval columns of \code{print()}
#' and \code{summary()} use these levels. The number of replicates matters
#' more here than for a parametric bootstrap, since the end points sit in the
#' tails of the replicates: use 200 or more.
#'
#' A block resample redraws the time slices and holds the actors fixed, so it
#' carries the variability of the data over time and none of the variability
#' over actors. The parameters whose sampling variability is that over actors
#' have no standard error and no interval under \code{type = "block"}: they
#' are \code{NA}, named in \code{withheld} with the reason, and \code{print()}
#' and \code{summary()} say so. \code{type = "parametric"} gives them. They are
#' the node-covariate coefficients, \code{va}, \code{vb} and \code{cab}, and,
#' for an empirical-Bayes fit, the intercept (the population intercept) and the
#' additive effects and factors, which are shrunk toward their prior mean: the
#' replicates vary around the shrunk estimates, not around the effects. The
#' residual variance \code{ve} has none either: a refit to a resample, which
#' repeats some slices and omits others, absorbs more of the noise than the
#' original fit, so its replicates sit below the estimate. The intercept of a
#' fixed-effects fit is the mean of its linear predictor over the observed
#' actors (the additive effects sum to zero), and its interval is for that
#' quantity and not for a population intercept.
#'
#' Coverage of the 95 percent interval was measured on directed normal-family
#' panels of 15 to 40 actors and 10 to 30 periods with a time-varying dyadic
#' covariate, static node covariates, random actor effects (standard
#' deviation 0.6) and reciprocal errors (correlation 0.3), at block lengths 1
#' to 3, with 100 to 200 replicates a bootstrap and 80 to 150 panels a cell
#' (1720 bootstraps of simulated panels in all, 530 of them with 200 replicates).
#' Coverage of the 95 percent interval, over all panels with ten or more blocks
#' (plain percentile), at the levels above with 200 replicates (expanded), and
#' for the parametric scheme on the same panels:
#' \tabular{lccc}{
#'   \strong{parameter} \tab \strong{plain} \tab \strong{expanded} \tab \strong{parametric} \cr
#'   dyadic coefficients \tab 0.90 \tab 0.95 \tab 0.93 \cr
#'   intercept, fixed-effects fit \tab 0.89 \tab 0.93 \tab 0.92 \cr
#'   \code{rho} \tab 0.90 \tab 0.94 \tab 0.94 \cr
#'   additive effects \tab 0.90 \tab 0.94 \tab 0.93 \cr
#'   multiplicative term (\code{R = 1}) \tab 0.91 to 0.92 \tab 0.94 \tab 0.93 \cr
#'   \code{ve} \tab 0.87 \tab 0.89 \tab 0.94 \cr
#'   node-covariate coefficients \tab 0.07 to 0.11 \tab 0.11 to 0.13 \tab 0.90 \cr
#'   \code{va}, \code{vb}, \code{cab} \tab 0.12 to 0.14 \tab 0.15 to 0.17 \tab 0.92 to 0.98
#' }
#' The first five rows keep their intervals (for a fit by the fixed-effects
#' estimator, \code{bias_correct = "none"}). The parametric column is the
#' interval of the parametric scheme with its own corrections. The intercept of a
#' fixed-effects fit is the mean of its linear predictor, so a replicate's
#' intercept carries the projection of the redrawn actor effects on the node
#' covariates, which the realized intercept does not: the replicate intercept
#' is taken relative to it (the interval was 2.6 times the sampling standard
#' deviation, coverage 1.00, before; now 0.92 of 500 panels, standard error 0.94
#' of the sampling standard deviation), and the node-covariate coefficients keep
#' that variation (0.90). The variances \code{va} and \code{vb} are biased and
#' their replicates reproduce the bias (a percentile interval covered 0.82 to
#' 0.89 of the fixed-effects panels and 0.65 to 0.75 of the empirical-Bayes
#' ones), so their parametric intervals are bias-corrected: the basic interval on
#' the log scale for a fixed-effects fit (0.92 to 0.98 at 15, 20 and 40 actors)
#' and the basic interval for an empirical-Bayes fit, whose estimate is biased
#' upward by an amount that does not scale with it (0.91 to 0.97 at 20 and 30
#' actors).
#' On binary empirical-Bayes panels of 20
#' actors and 12 periods (80 panels) the block intervals covered the dyadic
#' coefficient in 0.91 of the panels and \code{rho} in 0.94, and the intercept in
#' 0.19, the node-covariate coefficients in 0.25 and 0.26, \code{va},
#' \code{vb} and \code{cab} in 0.43, 0.28 and 0.33 and the additive effects in
#' 0.67 and 0.63 (0.93 for the parametric scheme). The same holds for the
#' empirical-Bayes normal fit, the default for the normal family (20 actors, 12
#' to 20 periods, 160 panels, 200 replicates): the block intervals covered the
#' dyadic coefficient in 0.96 to 1.00 of the panels and \code{rho} in 0.93 to
#' 0.98, and the intercept in 0.03 to 0.13, the node-covariate coefficients in
#' 0.05 to 0.10, \code{va}, \code{vb} and \code{cab} in 0.02 to 0.18 and the
#' additive effects in 0.27 to 0.39; \code{ve} covered 0.91 to 0.95 and has no
#' interval. The fixed-effects rows above are the fits made with
#' \code{bias_correct = "none"}.
#'
#' \strong{Fixed-effects fits.} A binary fit by the fixed-effects estimator
#' (\code{bias_correct = "none"}) has no bootstrap: an error of class
#' \code{"lame_als_bootstrap_bias"}, raised before anything is fitted by
#' \code{\link{ame_als}} and \code{\link{lame_als}}. The dyadic coefficient of
#' such a fit carries the incidental-parameters bias, every replicate is
#' refit at the biased estimate and reproduces it, and the percentile
#' interval sits about twice the bias from the truth. A 95 percent
#' interval covered the dyadic coefficient in 0.54, 0.72, 0.76 and 0.84 of 100
#' simulated panels of 30 actors and one period, 60 and two, 100 and two, and 30
#' and four; the bias relative to the standard error falls slowly with the number
#' of actors (0.85, 0.71, 0.57, 0.40) and the fit passes the size rule of
#' \code{\link{ame_als}} at 100 actors and two periods. The default
#' empirical-Bayes estimator, or the split-sample or parametric-bootstrap
#' correction, is the fit to bootstrap. The 95 percent interval of a corrected fit
#' covered the dyadic coefficient in 0.940 (jackknife) and 0.940 (bootstrap
#' correction) of 150 simulated panels of 30 actors and four periods, and in
#' 0.945 and 0.953 of 145 and 150 panels of 60 actors and two (200 replicates; Monte
#' Carlo standard error 0.019). Over those four cells the intervals for the dyadic
#' coefficient, \code{va} and \code{vb} covered 0.944, 0.940 and 0.943 of 595
#' panels (standard error 0.010), and the standard error of the dyadic
#' coefficient was 0.99 to 1.04 times its sampling standard deviation; 5 of
#' the 600 panels were refused by the refit rule or the size rule. A replicate of a
#' jackknife-corrected fit repeats the split-sample correction, and a replicate of
#' a bootstrap-corrected fit has the stored bias removed. Scaling that
#' stored bias by each replicate's own estimate changed the coverage of 60 panels
#' in each of two cells by at most 0.02 and is not used; repeating the
#' bootstrap correction in every replicate costs about 60 refits a replicate. A fixed-effects Poisson fit
#' has no incidental-parameters bias in its dyadic coefficient (0.93 to 0.94)
#' or its variance components, but its intercept does: it has no interval
#' under either scheme (0.61, 0.85 and 0.73 of panels of 30 actors and one
#' period, 60 and two, and 100 and two).
#'
#' @details
#' \strong{Choice of inference.} The Social Influence Regression paper of
#' Hoff & Minhas (2025) derives its primary standard errors from the observed
#' Hessian (classical \eqn{-H^{-1}} and the sandwich/robust estimator
#' \eqn{H^{-1} S H^{-1}}). Two features of the AME model make a Hessian-based
#' variance awkward here. First, without an explicit gauge fix the rank-\code{R}
#' multiplicative term is identified only up to a full \eqn{R\times R}
#' rotation/reflection (\eqn{U \to U R}, \eqn{V \to V R^{-\top}}), leaving the
#' joint Hessian rank-deficient. Second, even with a gauge fixed, the non-normal
#' families are fit on a Gaussian working response, so a Hessian computed from
#' that working objective is not calibrated to the family likelihood. The
#' bootstrap side-steps both issues and is the recommended uncertainty tool
#' here, mirroring \code{sir::boot_sir()}.
#'
#' \strong{Alignment-sensitive quantities.} The regression coefficients
#' \code{beta}, the additive effects \code{a}, \code{b} and the variance
#' components are rotation-invariant and are aggregated directly. The
#' multiplicative factors \code{U}, \code{V} are \emph{not}: each replicate is
#' Procrustes-aligned to the original fit before its standard errors are
#' computed. The object stores both the raw and aligned replicate factors so
#' the effect of alignment can be inspected. The alignment identifies a column
#' only when its singular value is separated from its neighbours': two adjacent
#' singular values (eigenvalues, for a symmetric model) count as separated when
#' their gap in the point estimate is at least 5 bootstrap standard deviations of
#' that gap (\code{sv_separation$threshold}; the rotation the alignment applies to
#' the pair is then under 4 degrees at the median and under 10 at the 90th
#' percentile, against 9 and 25 when the two are equal). The \code{U} and \code{V}
#' columns that are not separated from a neighbour have \code{NA} for the
#' column-wise standard errors (\code{se_U}, \code{se_V}, and their unaligned
#' \code{_raw} versions) and \code{NA} intervals in \code{confint()}. What is
#' reported for every fit with a multiplicative term, separated or not, is the
#' uncertainty of the subspaces the factors span, which a rotation of the columns
#' does not change: the principal angles, in degrees, between each replicate's
#' row-factor (and column-factor) space and the point estimate's
#' (\code{subspace}). For a symmetric model the eigenvalues \code{L} of the
#' multiplicative term are aggregated and reported as well.
#'
#' \strong{The parametric bootstrap} is a full parametric bootstrap of the AME
#' model: each replicate draws fresh additive random effects \code{a}, \code{b}
#' from their fitted dispersion (\code{va}, \code{vb}, \code{cab}) and fresh
#' residuals, holding \code{mu}, \code{beta} and the multiplicative term fixed.
#' Regenerating the additive effects (rather than holding the fitted \code{a},
#' \code{b} fixed) is what gives the intercept and node-covariate standard
#' errors their actor-level sampling-variability component.
#'
#' \strong{Assumptions.} The block bootstrap treats the time slices as
#' exchangeable replicates of the static-effects model --- appropriate for that
#' model, but not for strongly trended or serially dependent series (use the
#' dynamic \code{\link{lame}} there). The block bootstrap needs at least ten
#' resampling blocks, so \code{parametric} (the default) is
#' substituted, with a warning, on a panel that supplies fewer than ten
#' blocks. The finite-sample
#' (incidental-parameters) stretch of a binary IRLS fit is removed from the
#' point estimate and from every replicate by the empirical-Bayes estimator
#' (\code{bias_correct}, on by default for binary, Poisson and normal fits) or, when
#' requested, the split-sample jackknife or the parametric-bootstrap
#' correction; a binary fit that carries none of them has no bootstrap
#' (\emph{Fixed-effects fits} above). The MCMC
#' \code{\link{ame}} / \code{\link{lame}} path gives posterior summaries when
#' that is the target.
#'
#' @param object an \code{ame_als} object from \code{\link{ame_als}}
#'   or \code{\link{lame_als}}.
#' @param R integer number of bootstrap replicates (default 200); at least 10.
#'   A block bootstrap reads the ends of its interval from the tails of the
#'   replicates, so use 200 or more.
#' @param type \code{"parametric"} (default) or \code{"block"}; see Description.
#' @param block_length block length for the block bootstrap: \code{1} (default)
#'   resamples time slices independently; an integer \code{> 1} draws contiguous
#'   blocks of that many slices (a moving-block bootstrap), appropriate when the
#'   series has short-run temporal dependence. Capped at \code{T}; ignored for
#'   the parametric bootstrap.
#' @param seed optional integer random seed.
#' @param verbose logical; print progress (default \code{TRUE}).
#' @param max_iter positive integer iteration budget for each refit (default
#'   200). A replicate that does not converge is retried on the same data
#'   with four times this budget. A static Gaussian or working-scale refit
#'   that still fails also tries a fresh spectral start and retains it only
#'   if it converges with a lower objective. An unfinished low-rank Gaussian
#'   refit also tries seven starts from the fixed internal seed 6886, retaining
#'   only converged fits with no worse objective. A replicate that still has not
#'   converged is dropped, and more than 5 percent dropped is an error.
#' @param tol positive convergence tolerance for each refit (default
#'   \code{1e-6}).

#'
#' @return An object of class \code{"boot_ame"} with components including
#'   \code{coefs} (replicate intercept + regression coefficients), \code{se},
#'   \code{ci_lo}, \code{ci_hi}, \code{point_est}, \code{param_names};
#'   \code{vc_*} for the variance components; \code{a_coefs}, \code{b_coefs},
#'   \code{se_a}, \code{se_b}; \code{U_aligned}, \code{V_aligned},
#'   \code{U_raw}, \code{V_raw}, \code{se_U}, \code{se_V} (when \code{R > 0}
#'   and the factors are static; \code{NA} for the columns whose singular values
#'   are not separated), \code{sv_separation} (the singular values of the point
#'   estimate, the gaps and their z-scores in bootstrap standard deviations, the
#'   threshold, and which columns are \code{resolved}) and \code{subspace} (the
#'   principal angles, in degrees, of the replicate factor spaces to the point
#'   estimate's, and their \code{largest} angle by side);
#'   \code{n_valid}, \code{n_total}, \code{n_failed}, \code{fail_reason}
#'   (per replicate: \code{""}, \code{"nonconverged"}, \code{"split"},
#'   \code{"error"} or \code{"nonfinite"}), \code{max_fail_share},
#'   \code{boot_refit_start} (the start each replicate kept),
#'   \code{optimum_check}, \code{n_other_start} (replicates that kept a start
#'   other than the warm one), \code{vc_uncentred} (the variance components with
#'   no interval because their replicates are not centred on the estimate, when
#'   there are any), \code{withheld} (a data frame of the parameters with no
#'   standard error or interval, their \code{kind} --- \code{"coefficient"} or
#'   \code{"vc"} or \code{"effects"} --- and the \code{reason}; \code{NULL}
#'   when every parameter has one), \code{n_blocks} (the number of resampling blocks of a block
#'   bootstrap, \code{NA} otherwise), \code{estimator_eb} (whether the fit is
#'   empirical Bayes, which sets the interval of \code{va} and \code{vb}),
#'   \code{type}, \code{family}.
#'
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}. The block and parametric
#' bootstrap follow the inference scheme of \code{sir::boot_sir()} developed
#' for the SIR estimator.
#'
#' Hesterberg, T. C. (2015). What teachers should know about the bootstrap:
#' resampling in the undergraduate statistics curriculum. \emph{The American
#' Statistician} 69(4), 371--386.
#'
#' @seealso \code{\link{ame_als}} and \code{\link{lame_als}}, both of
#'   which accept \code{bootstrap = N}, \code{bootstrap_type},
#'   \code{bootstrap_block_length}, \code{bootstrap_seed} to do this
#'   work in a single call when you are about to fit the model
#'   anyway. \code{ame_als_bootstrap()} is the post-hoc path
#'   (bootstrap an \emph{already-fit} object without refitting from
#'   scratch). \code{\link{ame_als_refit}} for warm-start refits;
#'   \code{\link{confint.boot_ame}} for the bootstrap CI extractor.
#'
#' @examples
#' \donttest{
#' # twelve periods supply the ten resampling blocks the block scheme needs
#' Y <- replicate(12, { m <- matrix(rnorm(225), 15, 15); diag(m) <- NA; m },
#'                simplify = FALSE)
#' fit <- lame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' # post-hoc bootstrap of an existing fit:
#' bt <- ame_als_bootstrap(fit, R = 50, type = "block",
#'                             seed = 1, verbose = FALSE)
#' # equivalent one-shot call:
#' # fit_b <- lame_als(Y, R = 1, family = "normal", verbose = FALSE,
#' #                   bootstrap = 50, bootstrap_type = "block",
#' #                   bootstrap_seed = 1)
#' print(bt)
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
ame_als_bootstrap <- function(
	object,
	R = 200,
	type = c("parametric", "block"),
	block_length = 1,
	seed = NULL,
	verbose = TRUE,
	max_iter = 200L,
	tol = 1e-6
) {
	# a withheld product: replicates are drawn from its variational means and positions get no interval
	product_withheld <- .ae_product_is_withheld(object)
	object <- .ae_product_restore(object)

	if (!inherits(object, "ame_als")) {
		cli::cli_abort("{.arg object} must be an {.cls ame_als} object from {.fn ame_als} or {.fn lame_als}.")
	}
	if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 ||
			max_iter > .Machine$integer.max %/% 4L || max_iter != as.integer(max_iter) || length(tol) != 1L ||
			!is.finite(tol) || tol <= 0) {
		cli::cli_abort("Supply a positive integer {.arg max_iter} and positive finite {.arg tol}.")
	}
	if (length(R) != 1L || !is.numeric(R) || !is.finite(R) || R != round(R) || R < 1) {
		cli::cli_abort("{.arg R} must be a single positive whole number of replicates.")
	}
	R <- as.integer(R)
	# estimators without a bootstrap are refused before the replicate count is checked; a dynamic
	# empirical-Bayes fit carries its own standard errors
	dyn <- .ae_is_dynamic_als(object)
	if (dyn && identical(object$bias_correct$method, "eb")) {
		se_kind <- if (identical(object$family, "normal")) "the exact Gaussian standard errors given the hyperparameters" else "Laplace standard errors"
		cli::cli_abort(
			c(
				"A bootstrap is not offered for the dynamic empirical-Bayes estimator.",
				"i" = "It carries {se_kind}: {.code vcov(fit)}, {.code confint(fit)}, {.code summary(fit)} and {.code fit$eb$se}.",
				"i" = "The parametric bootstrap of the penalised normal fit is {.code bias_correct = \"none\"}; for posterior summaries use {.fn lame} with {.code method = \"mcmc\"}, optionally warm-started with {.code start_vals = als_start_vals(fit)}."
			),
			class = c("lame_dynamic_eb_unsupported", "lame_als_eb")
		)
	}
	# dynamic non-normal fits are refused: their replicates are not calibrated
	if (dyn && !identical(object$family, "normal")) {
		fam <- object$family
		cli::cli_abort(
			c(
				"A bootstrap is not available for a {.val {fam}} dynamic ALS fit.",
				"i" = "Each replicate refits through the IRLS reweighting loop and lands at a different local optimum: on a 60-actor binary {.arg dynamic_uv} panel the replicate standard deviations came out 8 to 11 times the MCMC posterior standard deviation, with 14% of replicates not converging.",
				"i" = "For uncertainty on this fit use MCMC warm-started from it: {.code lame(..., start_vals = als_start_vals(fit))}.",
				"i" = "The parametric bootstrap is available for normal-family dynamic fits and for every static ALS fit."
			),
			class = "lame_als_bootstrap_unsupported"
		)
	}
	# as is the uncorrected fixed-effects fit of a binary network
	if (identical(object$family, "binary") && identical(object$non_normal_method, "irls") &&
			!identical(object$bias_correct$method, "eb") && !isTRUE(object$bias_correct$applied)) {
		.ae_boot_refuse_uncorrected()
	}
	min_valid <- .ae_boot_min_valid
	min_req <- .ae_boot_min_requested
	if (R < min_valid) {
		cli::cli_abort(
			c(
				"{.arg R} = {R} cannot give the {.val {min_valid}} valid replicates a bootstrap needs.",
				"i" = "Use {.code R = 200}; any value of at least {.val {min_valid}} is accepted here, and {.code bootstrap = N} in {.fn ame_als} and {.fn lame_als} needs {.val {min_req}}."
			),
			class = "lame_als_bootstrap_replicates"
		)
	}
	if (!isTRUE(object$converged)) {
		cli::cli_abort(
			c(
				"Bootstrap intervals require a converged original fit.",
				"i" = "Refit with a larger iteration budget and inspect the convergence diagnostics."
			),
			class = "lame_als_bootstrap_convergence"
		)
	}
	type <- match.arg(type)
	# the requested scheme, kept even when a short panel or a dynamic fit uses another
	requested_type <- type
	if (length(block_length) != 1L || !is.finite(block_length) ||
			block_length < 1) {
		cli::cli_abort("{.arg block_length} must be a single positive integer.")
	}
	block_length <- as.integer(block_length)
	# seed without disturbing the global RNG stream, as in ame_als()
	if (!is.null(seed)) {
		if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			old_seed <- get(".Random.seed", envir = globalenv())
			on.exit(
				assign(".Random.seed", old_seed, envir = globalenv()),
				add = TRUE
			)
		} else {
			on.exit(
				if (exists(
					".Random.seed",
					envir = globalenv(),
					inherits = FALSE
				))
					rm(".Random.seed", envir = globalenv()),
				add = TRUE
			)
		}
		set.seed(seed)
	}
	if (dyn && identical(type, "block")) {
		# a block resample reorders the time axis that the smoothness penalties and per-slice actor sets follow
		cli::cli_warn(c(
			"{.arg type} = {.val block} is not available for a dynamic ALS fit.",
			"i" = "Resampling time slices breaks the temporal structure the dynamic estimator smooths along; using {.val parametric}."
		))
		type <- "parametric"
	}
	# block_length capped at T, the length the moving-block resampler uses
	blk <- min(block_length, object$n_time)

	Tt <- object$n_time
	Rlat <- object$R
	nr <- object$dims$n_row
	nc <- object$dims$n_col
	family <- object$family
	rho_hat <- as.numeric(object$VC["rho"])
	if (!is.finite(rho_hat)) rho_hat <- 0
	s2_hat <- object$s2
	nnm <- .ae_default(object$non_normal_method, "transform")
	# a non-normal `transform` fit is a Gaussian fit to a fixed transformed response, so it is
	# resampled on that working scale (z = ez + noise); IRLS and normal fits use the response scale
	working_boot <- identical(type, "parametric") &&
		family %in% c("binary", "poisson", "ordinal") &&
		!identical(nnm, "irls")

	# object$s2 is already df-corrected in .ae_assemble(), so simulating with it centres the refits
	# on the point estimate
	s2_sim <- s2_hat

	# fitted additive contribution, swapped per replicate for a fresh draw
	swap_ab <- !isTRUE(object$dynamic_ab)
	ab_orig <- if (swap_ab) {
		outer(as.numeric(object$a), as.numeric(object$b), "+")
	} else {
		matrix(0, nr, nc)
	}
	# time-varying additive effects are redrawn as whole paths (.ae_boot_ab_paths)
	dyn_ab <- isTRUE(object$dynamic_ab) && !is.null(object$a_dynamic) &&
		is.matrix(object$a_dynamic) && ncol(object$a_dynamic) == Tt
	ab_paths <- if (dyn_ab) .ae_boot_ab_paths(object)
	ab_path_orig <- ab_paths$orig
	draw_ab_path <- ab_paths$draw

	if (type == "block" && Tt < 2L) {
		cli::cli_abort(c(
			"Block bootstrap resamples time slices and needs {.arg T} > 1.",
			"i" = "This fit has a single time slice; use {.code type = \"parametric\"}.",
			"i" = "(A nodal/dyadic resample is not offered: it would duplicate actors and not yield a valid network.)"
		))
	}
	# the block scheme needs at least ten blocks; with fewer, the parametric one is used
	min_blocks <- .ae_min_boot_blocks
	n_blocks <- if (identical(type, "block")) as.integer(ceiling(Tt / blk)) else NA_integer_
	if (identical(type, "block") && n_blocks < min_blocks) {
		cli::cli_warn(c(
			"A block bootstrap needs at least {.val {min_blocks}} resampling blocks; this fit gives {n_blocks}.",
			"i" = "{Tt} time slice{?s} at {.arg block_length} = {blk} {?is/are} too few: the ends of the interval would be the extreme replicates.",
			"i" = "Using {.code type = \"parametric\"}; {.code fit$requested_type} records what was asked for."
		))
		type <- "parametric"
		working_boot <- identical(type, "parametric") &&
			family %in% c("binary", "poisson", "ordinal") &&
			!identical(nnm, "irls")
	}

	if (identical(type, "parametric") && .ae_rho_required(object)) {
		cli::cli_abort(
			c(
				"Dyadic dependence is not estimated from the available reciprocal pairs.",
				"i" = "Use {.fn ame} or {.fn lame} (MCMC) to integrate uncertainty in dyadic dependence."
			),
			class = "lame_als_prediction"
		)
	}

	# headline parameters: intercept and regression coefficients
	point_coef <- object$coefficients
	param_names <- names(point_coef)
	n_par <- length(point_coef)

	boot_coefs <- matrix(NA_real_, R, n_par)
	colnames(boot_coefs) <- param_names
	vc_coefs <- matrix(NA_real_, R, 5)
	colnames(vc_coefs) <- names(object$VC)
	a_coefs <- matrix(NA_real_, R, nr)
	b_coefs <- matrix(NA_real_, R, nc)
	boot_converged <- logical(R) # did the warm-start refit converge
	# the replicates need no marginal likelihood (see .ae_mf_refit)
	op_evidence <- options(lame.refit_evidence = FALSE)
	on.exit(options(op_evidence), add = TRUE)
	boot_refit_iterations <- integer(R)
	boot_refit_attempts <- integer(R)
	boot_refit_start <- rep("warm", R)
	fail_reason <- character(R) # why a replicate was dropped ("" = kept)
	n_split_dropped <- 0L # replicates the jackknife refused
	max_fail <- floor(.ae_boot_max_fail * R + 1e-9) # replicates that may be dropped
	first_failed <- NULL
	# a dynamic fit's factors are trajectories, with no single point estimate to align replicates to
	uv_track <- Rlat > 0 && is.matrix(object$U) && is.matrix(object[["V"]]) && !product_withheld
	U_raw <- V_raw <- U_aln <- V_aln <- L_coefs <- uv_ok <- sv_boot <- NULL
	if (uv_track) {
		U_raw <- matrix(NA_real_, R, nr * Rlat)
		U_aln <- matrix(NA_real_, R, nr * Rlat)
		V_raw <- matrix(NA_real_, R, nc * Rlat)
		V_aln <- matrix(NA_real_, R, nc * Rlat)
		sv_boot <- matrix(NA_real_, R, Rlat) # singular values of each replicate's product
		uv_ok <- logical(R) # finite U and V
		if (object$symmetric) L_coefs <- matrix(NA_real_, R, Rlat)
	}

	if (verbose) {
		cli::cli_h3("Bootstrap for fast AME fit")
		cli::cli_text("Type: {.field {type}} | Replicates: {.val {R}} | Family: {.field {family}}")
		cli::cli_text("Network: {nr} x {nc}, {Tt} time slice{?s}, R = {Rlat}")
	}

	for (b in seq_len(R)) {
		if (verbose && b %% 25 == 0) cli::cli_inform("Bootstrap {.val {b}}/{.val {R}}")

		Y_b <- NULL; X_b <- object[["X"]]; Z_b <- NULL
		intercept_shift <- 0

		if (type == "block") {
			if (blk <= 1L) {
				t_idx <- sample.int(Tt, Tt, replace = TRUE)
			} else {
				# moving-block resample: contiguous blocks of `blk` slices, wrapped circularly so every
				# slice has equal inclusion probability
				n_blk <- ceiling(Tt / blk)
				starts <- sample.int(Tt, n_blk, replace = TRUE)
				t_idx <- unlist(lapply(
					starts,
					function(s) ((s - 1L + seq_len(blk) - 1L) %% Tt) + 1L
				))[seq_len(Tt)]
			}
			Y_b <- object$Y[, , t_idx, drop = FALSE]
			X_b <- object[["X"]][, , , t_idx, drop = FALSE]
		} else if (working_boot) {
			# non-normal transform: resample on the Gaussian working scale, swapping in fresh additive
			# effects and adding residual noise with the fit's variance and reciprocity
			Z_b <- array(NA_real_, dim = c(nr, nc, Tt))
			ab_new <- if (swap_ab) .ae_boot_draw_ab(object) else ab_orig
			intercept_shift <- .ae_default(attr(ab_new, "intercept_shift"), 0)
			ab_path_new <- if (dyn_ab) draw_ab_path() else NULL
			for (t in seq_len(Tt)) {
				ez <- if (dyn_ab) {
					object$EZ[[t]] - ab_path_orig[[t]] + ab_path_new[[t]]
				} else object$EZ[[t]] - ab_orig + ab_new
				ez[!is.finite(ez)] <- 0
				zt <- .ae_boot_draw_response(object, ez, rho_hat, s2_sim, working = TRUE)
				storage.mode(zt) <- "double"
				if (object$symmetric) {
					zt[lower.tri(zt)] <- 0
					zt <- zt + t(zt)
				}
				zt[is.na(object$Y[, , t])] <- NA_real_
				Z_b[, , t] <- zt
			}
		} else {
			# parametric, response scale: simulate outcomes from the fitted model with fresh additive
			# effects; mu, beta and the multiplicative term stay fixed
			Y_b <- array(NA_real_, dim = c(nr, nc, Tt))
			ab_new <- if (swap_ab) .ae_boot_draw_ab(object) else ab_orig
			intercept_shift <- .ae_default(attr(ab_new, "intercept_shift"), 0)
			ab_path_new <- if (dyn_ab) draw_ab_path() else NULL
			for (t in seq_len(Tt)) {
				ez <- if (dyn_ab) {
					object$EZ[[t]] - ab_path_orig[[t]] + ab_path_new[[t]]
				} else object$EZ[[t]] - ab_orig + ab_new
				ez[!is.finite(ez)] <- 0 # the unipartite diagonal of EZ is NA
				yt <- .ae_boot_draw_response(object, ez, rho_hat, s2_sim)
				yt <- as.matrix(yt)
				storage.mode(yt) <- "double"
				# symmetric model: mirror the upper triangle so each tie has one draw and the right variance
				if (object$symmetric) {
					yt[lower.tri(yt)] <- 0
					yt <- yt + t(yt)
				}
				# keep the observed-cell pattern of the data, so a replicate carries no more information than
				# the original fit had
				yt[is.na(object$Y[, , t])] <- NA_real_
				Y_b[, , t] <- yt
			}
		}

		# refit from the warm start (and the original fit's start, see .ae_boot_refit); per-replicate
		# warnings are suppressed since the original fit already raised them
		res <- tryCatch(
			.ae_boot_refit(object, Y_b, X_b, Z_b, dyn, max_iter, tol),
			error = function(e) list(
				fit = NULL,
				attempts = 0L,
				iterations = max_iter,
				start = NA_character_,
				error = e,
				last_fit = NULL
			)
		)
		boot_refit_attempts[b] <- res$attempts
		boot_refit_iterations[b] <- res$iterations
		fit_b <- res$fit
		reason <- ""
		if (is.null(fit_b)) {
			reason <- if (inherits(res$error, "lame_jk_split")) {
				n_split_dropped <- n_split_dropped + 1L
				"split"
			} else if (!is.null(res$error) && !inherits(res$error, "lame_jk_nonconvergence")) {
				"error"
			} else "nonconverged"
		} else {
			boot_refit_start[b] <- res$start
			cf <- fit_b$coefficients
			eff_b <- .ae_boot_effects(fit_b, object)
			if (all(is.finite(cf)) && eff_b$ok) {
				boot_coefs[b, ] <- cf
				boot_coefs[b, 1L] <- cf[[1L]] - intercept_shift
				vc_coefs[b, ] <- fit_b$VC
				a_coefs[b, ] <- eff_b$a
				b_coefs[b, ] <- eff_b$b
				boot_converged[b] <- isTRUE(fit_b$converged)
				if (uv_track && is.matrix(fit_b$U) && is.matrix(fit_b[["V"]]) &&
						all(is.finite(fit_b$U)) && all(is.finite(fit_b[["V"]]))) {
					U_raw[b, ] <- c(fit_b$U)
					V_raw[b, ] <- c(fit_b[["V"]])
					# Procrustes-align the replicate U, V to the original fit
					src <- rbind(fit_b$U, fit_b[["V"]])
					tgt <- rbind(object$U, object[["V"]])
					Rrot <- .ae_procrustes(src, tgt)
					U_aln[b, ] <- c(fit_b$U %*% Rrot)
					V_aln[b, ] <- c(fit_b[["V"]] %*% Rrot)
					sv_boot[b, ] <- .ae_boot_sv(fit_b)
					if (object$symmetric && !is.null(fit_b$L)) {
						L_coefs[b, ] <- fit_b$L
					}
					uv_ok[b] <- TRUE
				}
			} else {
				reason <- "nonfinite"
			}
		}
		if (nzchar(reason)) {
			fail_reason[b] <- reason
			if (is.null(first_failed)) first_failed <- list(
				replicate = b,
				fit = res$last_fit,
				error = res$error
			)
			# stop once the dropped share exceeds the cap, since later replicates cannot bring it back
			if (sum(nzchar(fail_reason)) > max_fail) {
				.ae_boot_abort_failed(
					fail_reason,
					R,
					max_fail,
					first_failed,
					max_iter,
					boot_refit_attempts[b],
					b
				)
			}
		}
	}

	valid <- stats::complete.cases(boot_coefs)
	n_valid <- sum(valid)
	n_nonconverged <- sum(valid & !boot_converged)
	n_failed <- R - n_valid
	if (n_valid < min_valid) {
		cli::cli_abort(
			c(
				"Only {.val {n_valid}} of {.val {R}} bootstrap replicates are valid; a bootstrap needs at least {.val {min_valid}}.",
				"i" = "Increase {.arg R}, or simplify the model."
			),
			class = "lame_als_bootstrap_replicates"
		)
	}
	if (n_failed > 0L && verbose) {
		cli::cli_alert_info("{n_failed} of {R} replicate{?s} could not be refit and {?was/were} dropped ({.val {sum(fail_reason == 'nonconverged')}} not converged, {.val {sum(fail_reason == 'split')}} split-sample correction not carried, {.val {sum(fail_reason %in% c('error', 'nonfinite'))}} refit error or non-finite).")
	}

	col_sd <- function(M, v) apply(M[v, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
	col_q <- function(M, v, pr) apply(
		M[v, , drop = FALSE],
		2,
		stats::quantile,
		probs = pr,
		na.rm = TRUE
	)

	# end points of the 95% interval: expanded percentile levels for a block resample
	n_blocks_used <- if (identical(type, "block")) n_blocks else NA_integer_
	pr95 <- .ae_level_probs(0.95, type, n_blocks_used)
	ci_lo_vec <- col_q(boot_coefs, valid, pr95[1])
	ci_hi_vec <- col_q(boot_coefs, valid, pr95[2])

	# variance intervals: bias-corrected for va and vb of a parametric bootstrap (.ae_vc_interval)
	vc_ci <- col_q(vc_coefs, valid, pr95)
	if (identical(type, "parametric")) {
		for (v in intersect(c("va", "vb"), colnames(vc_coefs))) {
			vc_ci[, v] <- .ae_vc_interval(
				vc_ci[, v],
				as.numeric(object$VC[v]),
				identical(object$bias_correct$method, "eb")
			)
		}
	}
	out <- list(
		coefs = boot_coefs,
		se = col_sd(boot_coefs, valid),
		ci_lo = ci_lo_vec,
		ci_hi = ci_hi_vec,
		point_est = point_coef,
		param_names = param_names,
		vc_coefs = vc_coefs,
		vc_se = col_sd(vc_coefs, valid),
		vc_ci_lo = vc_ci[1L, ],
		vc_ci_hi = vc_ci[2L, ],
		vc_point_est = object$VC,
		a_coefs = a_coefs,
		b_coefs = b_coefs,
		se_a = .ae_boot_col_sd(a_coefs, valid),
		se_b = .ae_boot_col_sd(b_coefs, valid),
		point_a = object$a,
		point_b = object$b,
		n_valid = n_valid,
		n_total = R,
		n_nonconverged = n_nonconverged,
		n_failed = n_failed,
		fail_reason = fail_reason,
		max_fail_share = .ae_boot_max_fail,
		boot_converged = boot_converged,
		boot_refit_attempts = boot_refit_attempts,
		boot_refit_start = boot_refit_start,
		boot_refit_iterations = boot_refit_iterations,
		refit_tol = tol,
		n_split_dropped = n_split_dropped,
		optimum_check = .ae_boot_optimum_check(object, dyn),
		n_other_start = sum(valid & boot_refit_start != "warm"),
		uv_summary = uv_track,
		type = type,
		requested_type = requested_type,
		block_length = if (type == "block") blk else NA_integer_,
		n_blocks = n_blocks_used,
		estimator_eb = identical(object$bias_correct$method, "eb"),
		family = family,
		mode = object$mode,
		symmetric = object$symmetric,
		R = Rlat,
		# actor names in dims, so confint() labels intervals with them
		dims = c(
			object$dims,
			list(
				row_names = names(object$a) %||% object$row_names,
				col_names = names(object$b) %||% object$col_names
			)
		),
		n_time = Tt,
		call = match.call()
	)

	if (uv_track) {
		vuv <- valid & uv_ok
		out$n_valid_uv <- sum(vuv)
		if (out$n_valid_uv > 0) {
			out$U_raw <- U_raw; out$V_raw <- V_raw
			out$U_aligned <- U_aln; out$V_aligned <- V_aln
			out$point_U <- object$U; out$point_V <- object[["V"]]
			# columns whose singular values are not separated from a neighbour's have no column-wise standard error
			sep <- .ae_boot_separation(object, sv_boot[vuv, , drop = FALSE])
			out$sv_separation <- sep
			out$subspace <- .ae_boot_subspace(
				object,
				U_raw[vuv, , drop = FALSE],
				V_raw[vuv, , drop = FALSE]
			)
			na_cols <- function(se, n) {
				for (k in which(!sep$resolved)) se[(k - 1L) * n + seq_len(n)] <- NA_real_
				se
			}
			out$se_U <- na_cols(col_sd(U_aln, vuv), nr)
			out$se_V <- na_cols(col_sd(V_aln, vuv), nc)
			out$se_U_raw <- na_cols(col_sd(U_raw, vuv), nr)
			out$se_V_raw <- na_cols(col_sd(V_raw, vuv), nc)
			if (object$symmetric && !is.null(L_coefs)) {
				out$L_coefs <- L_coefs
				out$se_L <- col_sd(L_coefs, vuv)
				out$ci_L_lo <- col_q(L_coefs, vuv, pr95[1])
				out$ci_L_hi <- col_q(L_coefs, vuv, pr95[2])
				out$point_L <- object$L
			}
		}
	}

	out <- .ae_boot_apply_withheld(out, .ae_boot_withheld(object, type))
	class(out) <- "boot_ame"

	# a coefficient outside the 95% interval of its replicates is an error; a variance component
	# outside its interval has none (NA, named in vc_uncentred)
	.ae_boot_check_centred(out)
}

#' @rdname ame_als_bootstrap
#' @export
boot_ame <- function(
	object,
	R = 200,
	type = c("parametric", "block"),
	block_length = 1,
	seed = NULL,
	verbose = TRUE,
	max_iter = 200L,
	tol = 1e-6
) {
	ame_als_bootstrap(
		object,
		R = R,
		type = type,
		block_length = block_length,
		seed = seed,
		verbose = verbose,
		max_iter = max_iter,
		tol = tol
	)
}

# S3 methods for boot_ame objects

#' Print bootstrap results for a fast AME fit
#'
#' @param x a \code{boot_ame} object.
#' @param digits number of digits to display.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @method print boot_ame
#' @export
print.boot_ame <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	.lame_cli_stdout({
		drop_note <- if (isTRUE(x$n_failed > 0)) {
			paste0(" (", x$n_failed, " dropped: not converged, refit error or split correction not carried)")
		} else {
			""
		}
		cli::cli_text("")
		cli::cli_text("{.strong Bootstrap results: fast AME estimator}")
		cli::cli_text("Type: {.field {x$type}} | Family: {.field {x$family}} | Replicates: {.val {x$n_valid}}/{.val {x$n_total}} valid{drop_note}")
		tab <- cbind(
			Estimate = x$point_est,
			`Boot SE` = x$se,
			`2.5%` = x$ci_lo,
			`97.5%` = x$ci_hi
		)
		rownames(tab) <- x$param_names
		cli::cli_text("")
		cli::cli_verbatim(utils::capture.output(print(round(tab, digits))))
		cli::cli_text("")
		cli::cli_text(
			cli::col_grey("beta, a, b and variance components aggregate directly; U/V are Procrustes-aligned.")
		)
		.ae_print_block_levels(x)
		.ae_print_withheld(x)
	})
	invisible(x)
}

# levels the end points of a block interval are read at
.ae_print_block_levels <- function(x) {
	if (identical(x$type, "block") && is.finite(x$n_blocks %||% NA_real_)) {
		pr <- 100 * .ae_boot_probs(x, 0.95)
		cli::cli_text(
			cli::col_grey("Block resample of {x$n_blocks} blocks: the interval ends are the {signif(pr[1], 2)}% and {signif(pr[2], 3)}% quantiles of the replicates.")
		)
	}
	invisible(NULL)
}

# one line per reason on the parameters without a standard error or interval
.ae_print_withheld <- function(x) {
	lines <- .ae_boot_withheld_lines(x)
	for (l in lines) cli::cli_text("{l}")
	if (length(lines) && identical(x$type, "block")) {
		cli::cli_text("{.code type = \"parametric\"} gives intervals for them.")
	}
	invisible(NULL)
}

#' Summarize bootstrap results for a fast AME fit
#'
#' @param object a \code{boot_ame} object.
#' @param ... ignored.
#' @return An object of class \code{"summary.boot_ame"}, printed as a set of
#'   tables.
#' @method summary boot_ame
#' @export
summary.boot_ame <- function(object, ...) {
	structure(list(boot = object), class = "summary.boot_ame")
}

#' @rdname summary.boot_ame
#' @param x a \code{summary.boot_ame} object.
#' @method print summary.boot_ame
#' @export
print.summary.boot_ame <- function(x, ...) {
	.lame_cli_stdout({
		object <- x$boot
		cli::cli_h1("Bootstrap results: fast AME estimator")
		cli::cli_text("Type: {.field {object$type}} | Family: {.field {object$family}} | Mode: {.field {object$mode}}")
		if (identical(object$type, "block") && isTRUE(object$block_length > 1L)) {
			cli::cli_text("Moving-block resample: block length {.val {object$block_length}}")
		}
		.ae_print_block_levels(object)
		cli::cli_text("Replicates: {.val {object$n_valid}} valid of {.val {object$n_total}} ({round(100 * object$n_valid / object$n_total, 1)}%); every one converged")
		if (isTRUE(object$n_failed > 0)) {
			cli::cli_text("Dropped (not in the intervals): {.val {object$n_failed}}, at most {.val {floor(object$max_fail_share * object$n_total)}} allowed")
		}
		if (!is.null(object$optimum_check)) {
			cli::cli_text("Optimum check of each replicate: {object$optimum_check}")
		}
		cli::cli_text("Network: {object$dims$n_row} x {object$dims$n_col}, {object$n_time} time slice{?s}, R = {object$R}")

		cli::cli_rule()
		cli::cli_text("{.strong Regression coefficients}")
		tab <- cbind(
			Estimate = object$point_est,
			`Boot SE` = object$se,
			`2.5%` = object$ci_lo,
			`97.5%` = object$ci_hi
		)
		rownames(tab) <- object$param_names
		covers0 <- object$ci_lo <= 0 & object$ci_hi >= 0
		tab <- cbind(round(tab, 4), sig = ifelse(is.na(covers0) | covers0, "", "*"))
		cli::cli_verbatim(utils::capture.output(print(noquote(tab))))
		cli::cli_text(cli::col_grey("* = 95% percentile CI excludes zero"))


		cli::cli_rule()
		cli::cli_text("{.strong Variance components}")
		vtab <- cbind(
			Estimate = object$vc_point_est,
			`Boot SE` = object$vc_se,
			`2.5%` = object$vc_ci_lo,
			`97.5%` = object$vc_ci_hi
		)
		cli::cli_verbatim(utils::capture.output(print(round(vtab, 4))))
		if (!is.null(object$vc_uncentred) && nrow(object$vc_uncentred)) {
			un <- object$vc_uncentred
			items <- paste0(
				un$component,
				" (estimate ",
				signif(un$estimate, 4),
				", replicates ",
				signif(un$lo, 4),
				" to ",
				signif(un$hi, 4),
				")"
			)
			cli::cli_text("No interval for {paste(items, collapse = '; ')}: the replicates are not centred on the estimate, which a percentile interval cannot describe.")
		}

		.ae_print_withheld(object)

		cli::cli_rule()
		cli::cli_text("{.strong Additive effects} (rotation-invariant, aggregated directly)")
		# an actor at the boundary in some replicates has an infinite SE; the mean is over the others,
		# with the count shown
		se_line <- function(se, lbl) {
			if (!any(is.finite(se))) {
				cli::cli_text("{lbl}: no standard error")
				return(invisible(NULL))
			}
			n_inf <- sum(is.infinite(se))
			txt <- sprintf(
				"%s: mean |SE| = %s",
				lbl,
				round(mean(se[is.finite(se)]), 4)
			)
			if (n_inf > 0L) {
				txt <- sprintf(
					"%s (%d actor%s at the boundary in some replicates: SE Inf, interval open on that side)",
					txt,
					n_inf,
					if (n_inf > 1L) "s" else ""
				)
			}
			cli::cli_text("{txt}")
		}
		se_line(object$se_a, "Sender a")
		se_line(object$se_b, "Receiver b")

		if (object$R > 0) {
			if (isFALSE(object$uv_summary)) {
				cli::cli_text("{.strong Multiplicative factors}: time-varying, so replicates carry no common point estimate to align to.")
			} else if (is.null(object$se_U)) {
				cli::cli_text("{.strong Multiplicative factors}: no replicate produced finite U/V.")
			} else if (object$symmetric) {
				cli::cli_text("{.strong Multiplicative factors} (symmetric eigenmodel; {object$n_valid_uv} replicate{?s})")
				cli::cli_text("Eigenvalues L  point: {paste(sprintf('%.3f', object$point_L), collapse = ', ')}")
				cli::cli_text("Eigenvalues L  boot SE: {paste(sprintf('%.3f', object$se_L), collapse = ', ')}")
				cli::cli_text(
					cli::col_grey("L (eigenvalues) is the primary multiplicative-uncertainty summary; U columns are sign/eigenspace-identified only.")
				)
			} else {
				cli::cli_text("{.strong Multiplicative factors} (Procrustes-aligned; {object$n_valid_uv} replicate{?s})")
				if (any(object$sv_separation$resolved) && any(is.finite(object$se_U))) {
					cli::cli_text("U: mean SE = {round(mean(object$se_U, na.rm = TRUE), 4)} (raw, unaligned: {round(mean(object$se_U_raw, na.rm = TRUE), 4)})")
					cli::cli_text("V: mean SE = {round(mean(object$se_V, na.rm = TRUE), 4)} (raw, unaligned: {round(mean(object$se_V_raw, na.rm = TRUE), 4)})")
				}
			}
		}
		if (object$R > 0 && !isFALSE(object$uv_summary) && !is.null(object$sv_separation)) {
			sp <- object$sv_separation
			if (!all(sp$resolved)) {
				cols <- paste(which(!sp$resolved), collapse = ", ")
				gaps <- paste(sprintf("%.1f", sp$z), collapse = ", ")
				cli::cli_text("Factor column(s) {cols} not separated, so without a column-wise standard error (NA). Adjacent singular values are separated when their gap is at least {sp$threshold} bootstrap standard deviations of the gap; the gaps are {gaps}.")
			}
			la <- object$subspace$largest
			u_txt <- paste(sprintf("%.1f", la["U", c("median", "q95")]), collapse = ", ")
			v_txt <- paste(sprintf("%.1f", la["V", c("median", "q95")]), collapse = ", ")
			cli::cli_text("Subspace, largest principal angle to the point estimate (degrees; median, 95th percentile): U {u_txt}; V {v_txt}")
		}
	})
	invisible(x)
}

#' Confidence intervals from a fast AME bootstrap
#'
#' Bootstrap confidence intervals for the intercept and regression
#' coefficients of a fast AME fit.
#'
#' @param object a \code{boot_ame} object.
#' @param parm character vector of parameter names, or integer indices. If
#'   \code{NULL} (default), all coefficients are returned.
#' @param level confidence level (default 0.95).
#' @param ci_type interval type: \code{"percentile"} (default) uses the
#'   replicate quantiles directly (for a block bootstrap, those at the
#'   expanded percentile levels of \code{\link{ame_als_bootstrap}}) and
#'   matches the intervals reported by
#'   \code{summary()} and \code{tidy()} on the same fit; \code{"basic"}
#'   reflects the replicate quantiles about the point estimate
#'   (\eqn{2\hat\theta - q}), which corrects first-order bias and can be
#'   preferable for the IRLS binary/Poisson point estimators on small
#'   samples.
#' @param which which uncertainty channel to return: \code{"all"} (default,
#'   intercept + regression coefficients + variance components + sender effects
#'   + receiver effects + multiplicative U/V), or a single channel:
#'   \code{"beta"}, \code{"vc"}, \code{"a"}, \code{"b"}, \code{"U"},
#'   \code{"V"}.
#' @param ... ignored.
#' @return A matrix with one row per parameter and lower/upper bound columns.
#'   The rows of the \code{U} and \code{V} columns whose singular values are
#'   not separated from a neighbour's, of a variance component with no
#'   interval (\code{vc_uncentred}), and of a parameter the bootstrap returns
#'   no interval for (\code{withheld}), are \code{NA}; see
#'   \code{\link{ame_als_bootstrap}}.
#' @method confint boot_ame
#' @export
confint.boot_ame <- function(
	object,
	parm = NULL,
	level = 0.95,
	ci_type = c("percentile", "basic"),
	which = c("all", "beta", "vc", "a", "b", "U", "V"),
	...
) {
	ci_type <- match.arg(ci_type)
	which <- match.arg(which)
	if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
		cli::cli_abort("{.arg level} must be a single number in (0, 1).")
	}
	a <- (1 - level) / 2
	pct <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")
	# end-point probabilities: nominal, or the expanded ones of a block resample
	pr <- .ae_boot_probs(object, level)

	# interval matrix from a draws matrix and point estimates; columns that are NA in every
	# replicate (cab on a bipartite fit) give NA quantiles
	.boot_ci <- function(draws, point_est, names_, prefix = "") {
		if (is.null(draws) || nrow(draws) == 0L || ncol(draws) == 0L) {
			return(matrix(numeric(0), 0, 2, dimnames = list(NULL, pct)))
		}
		q <- t(apply(draws, 2, function(col) {
			# infinite replicate values are boundary estimates and stay in
			col <- col[!is.na(col)]
			if (length(col) == 0L) return(c(NA_real_, NA_real_))
			stats::quantile(col, probs = pr)
		}))
		out <- if (ci_type == "percentile") {
			q
		} else {
			# basic bootstrap interval: 2 * theta_hat minus the upper/lower percentile
			cbind(2 * point_est - q[, 2], 2 * point_est - q[, 1])
		}
		rownames(out) <- if (nzchar(prefix)) paste0(prefix, names_) else names_
		colnames(out) <- pct
		out
	}

	parts <- list()
	if (which %in% c("all", "beta")) {
		parts$beta <- .boot_ci(
			object$coefs,
			object$point_est,
			object$param_names
		)
		# parameters without a bootstrap interval (see withheld)
		parts$beta[intersect(
			.ae_boot_withheld_names(object, "coefficient"),
			rownames(parts$beta)
		), ] <- NA_real_
	}
	if (which %in% c("all", "vc") && !is.null(object$vc_coefs)) {
		vc <- .boot_ci(
			object$vc_coefs,
			object$vc_point_est,
			colnames(object$vc_coefs)
		)
		# va and vb of a parametric bootstrap get a bias-corrected interval (.ae_vc_interval)
		if (ci_type == "percentile" && identical(object$type, "parametric")) {
			for (v in intersect(c("va", "vb"), rownames(vc))) {
				vc[v, ] <- .ae_vc_interval(vc[v, ], as.numeric(object$vc_point_est[v]), isTRUE(object$estimator_eb))
			}
		}
		# variance components are non-negative, but the basic interval can cross zero near the
		# boundary, so clip va, vb and ve at zero
		non_neg <- intersect(c("va", "vb", "ve"), rownames(vc))
		if (length(non_neg) > 0L) {
			vc[non_neg, 1] <- pmax(vc[non_neg, 1], 0, na.rm = TRUE)
			vc[non_neg, 2] <- pmax(vc[non_neg, 2], 0, na.rm = TRUE)
		}
		# no interval for components whose replicates are not centred on the estimate, or that the
		# bootstrap withholds
		vc[intersect(
			c(object$vc_uncentred$component, .ae_boot_withheld_names(object, "vc")),
			rownames(vc)
		), ] <- NA_real_
		parts$vc <- vc
	}
	if (which %in% c("all", "a") && !is.null(object$a_coefs)) {
		nm_a <- if (!is.null(object$dims$row_names)) object$dims$row_names
		else colnames(object$a_coefs) %||% paste0("a", seq_len(ncol(object$a_coefs)))
		parts$a <- .boot_ci(
			object$a_coefs,
			object$point_a,
			nm_a,
			prefix = "a_"
		)
		if ("a" %in% .ae_boot_withheld_names(object, "effects")) parts$a[] <- NA_real_
	}
	if (which %in% c("all", "b") && !is.null(object$b_coefs) &&
			nrow(object$b_coefs) > 0L && ncol(object$b_coefs) > 0L) {
		nm_b <- if (!is.null(object$dims$col_names)) object$dims$col_names
		else colnames(object$b_coefs) %||% paste0("b", seq_len(ncol(object$b_coefs)))
		parts$b <- .boot_ci(
			object$b_coefs,
			object$point_b,
			nm_b,
			prefix = "b_"
		)
		if ("b" %in% .ae_boot_withheld_names(object, "effects")) parts$b[] <- NA_real_
	}
	# U and V are stored flat: u_aligned is [r_boot, n_row * r_lat], each row c(point_u) of one
	# replicate (column-major)
	.uv_names <- function(actors, Rlat, label) {
		# column-major order matches c(matrix(...)): position i + (r - 1) * n
		as.vector(outer(
			actors,
			paste0("d", seq_len(Rlat)),
			FUN = function(a, d) paste0(label, "_", d, "_", a)
		))
	}
	# columns not separated from a neighbour get no column-wise interval (see sv_separation)
	.na_unresolved <- function(ci, n, Rlat) {
		res <- object$sv_separation$resolved
		if (is.null(res)) return(ci)
		for (k in which(!res)) ci[(k - 1L) * n + seq_len(n), ] <- NA_real_
		ci
	}
	if (which %in% c("all", "U") && !is.null(object$U_aligned) &&
			is.matrix(object$U_aligned) && nrow(object$U_aligned) > 0L) {
		nr <- nrow(object$point_U)
		Rlat <- ncol(object$point_U)
		actors <- if (!is.null(rownames(object$point_U))) rownames(object$point_U)
		else paste0("a", seq_len(nr))
		parts$U <- .na_unresolved(
			.boot_ci(
				object$U_aligned,
				as.vector(object$point_U),
				.uv_names(actors, Rlat, "U")
			),
			nr,
			Rlat
		)
		if ("U" %in% .ae_boot_withheld_names(object, "effects")) parts$U[] <- NA_real_
	}
	if (which %in% c("all", "V") && !is.null(object$V_aligned) &&
			is.matrix(object$V_aligned) && nrow(object$V_aligned) > 0L &&
			ncol(object$V_aligned) > 0L) {
		nc <- nrow(object$point_V)
		Rlat <- ncol(object$point_V)
		actors <- if (!is.null(rownames(object$point_V))) rownames(object$point_V)
		else paste0("a", seq_len(nc))
		parts[["V"]] <- .na_unresolved(
			.boot_ci(
				object$V_aligned,
				as.vector(object$point_V),
				.uv_names(actors, Rlat, "V")
			),
			nc,
			Rlat
		)
		if ("V" %in% .ae_boot_withheld_names(object, "effects")) parts[["V"]][] <- NA_real_
	}

	ci <- do.call(rbind, parts)
	if (is.null(ci) || nrow(ci) == 0L) {
		return(matrix(numeric(0), 0, 2, dimnames = list(NULL, pct)))
	}

	if (!is.null(parm)) {
		if (is.character(parm)) {
			bad <- setdiff(parm, rownames(ci))
			if (length(bad)) {
				cli::cli_abort(c(
					"Unknown parameter name{?s} in {.arg parm}: {.val {bad}}.",
					"i" = "Available: {.val {head(rownames(ci), 20)}}{cli::qty(max(0, nrow(ci) - 20))}{?/, ...}."
				))
			}
		} else if (is.numeric(parm) &&
				(any(parm < 1) || any(parm > nrow(ci)))) {
			cli::cli_abort("{.arg parm} index is out of range (1:{nrow(ci)}).")
		}
		ci <- ci[parm, , drop = FALSE]
	}
	ci
}

#' Bootstrap covariance of the regression coefficients
#'
#' Returns the sample covariance matrix of the bootstrap replicate
#' intercept + regression coefficients of a \code{boot_ame} object -- the
#' covariance underlying the reported standard errors and confidence intervals.
#'
#' @param object a \code{boot_ame} object.
#' @param ... ignored.
#' @return A covariance matrix over the intercept and regression coefficients;
#'   the rows and columns of a parameter the bootstrap returns no standard
#'   error for (\code{withheld}) are \code{NA}.
#' @method vcov boot_ame
#' @export
vcov.boot_ame <- function(object, ...) {
	valid <- stats::complete.cases(object$coefs)
	V <- stats::cov(object$coefs[valid, , drop = FALSE])
	dimnames(V) <- list(object$param_names, object$param_names)
	# parameters without a bootstrap standard error get no covariance
	w <- intersect(.ae_boot_withheld_names(object, "coefficient"), object$param_names)
	V[w, ] <- NA_real_
	V[, w] <- NA_real_
	V
}

#' fitted/residuals are not defined for a bootstrap object
#'
#' A \code{boot_ame} stores replicate-level distributional summaries, not
#' per-cell predictions or residuals. These methods exist only to refuse the
#' call clearly (instead of falling through to \code{stats::fitted.default} /
#' \code{stats::residuals.default}, which would silently return \code{NULL}).
#' Apply \code{\link{fitted}} / \code{\link{residuals}} to the underlying
#' \code{ame_als} fit instead.
#'
#' @param object a \code{boot_ame} object.
#' @param ... ignored.
#' @return Never returns; raises an error.
#' @name boot_ame-no-fitted
NULL

#' @rdname boot_ame-no-fitted
#' @method fitted boot_ame
#' @export
fitted.boot_ame <- function(object, ...) {
	cli::cli_abort(c(
		"{.fn fitted} is not defined for a {.cls boot_ame} object.",
		"i" = "A bootstrap aggregates per-replicate fits; call {.fn fitted} on the underlying {.cls ame_als} fit."
	))
}

#' @rdname boot_ame-no-fitted
#' @method residuals boot_ame
#' @export
residuals.boot_ame <- function(object, ...) {
	cli::cli_abort(c(
		"{.fn residuals} is not defined for a {.cls boot_ame} object.",
		"i" = "Call {.fn residuals} on the underlying {.cls ame_als} fit."
	))
}

#' Point estimates from a fast AME bootstrap
#'
#' Returns the intercept and regression coefficient point estimates carried by
#' a \code{boot_ame} object (the estimates the bootstrap quantifies). For their
#' bootstrap standard errors and intervals see \code{\link{confint.boot_ame}}
#' and \code{\link{vcov.boot_ame}}.
#'
#' @param object a \code{boot_ame} object.
#' @param ... ignored.
#' @return A named numeric vector of coefficient point estimates.
#' @method coef boot_ame
#' @export
coef.boot_ame <- function(object, ...) {
	cf <- object$point_est
	names(cf) <- object$param_names
	cf
}

# sampler_describe generic

#' Describe the estimator behind a fitted object
#'
#' @description
#' Reports which estimation routine produced a fitted object and what
#' uncertainty information, if any, is available. Methods exist for
#' \code{ame_als} fits (the fast block coordinate descent estimator) and
#' \code{boot_ame} bootstrap objects.
#'
#' @param object a fitted object (\code{ame_als} or \code{boot_ame}).
#' @param verbose logical; if \code{TRUE} (default) print the description.
#' @return The input \code{object}, invisibly.
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' sampler_describe(fit)
#' @export
sampler_describe <- function(object, verbose = TRUE) {
	UseMethod("sampler_describe", object)
}

#' @rdname sampler_describe
#' @method sampler_describe ame_als
#' @export
sampler_describe.ame_als <- function(object, verbose = TRUE) {
	if (verbose) {
		.lame_cli_stdout({
			cli::cli_h3("Estimator")
			cli::cli_text("Method: {.field iterative block coordinate descent} (fast, MCMC-free)")
			cli::cli_text("Origin: Hoff & Minhas Social Influence Regression estimator (sir::sir_alsfit())")
			cli::cli_text("Uncertainty available: no (point estimate only)")
			cli::cli_text("For standard errors, run {.fn ame_als_bootstrap}.")
			cli::cli_text(cli::col_grey(object$meta$approximation_note))
		})
	}
	invisible(object)
}

#' @rdname sampler_describe
#' @method sampler_describe boot_ame
#' @export
sampler_describe.boot_ame <- function(object, verbose = TRUE) {
	if (verbose) {
		.lame_cli_stdout({
			cli::cli_h3("Estimator")
			cli::cli_text("Method: {.field {object$type}} bootstrap of the iterative block coordinate descent estimator")
			cli::cli_text("Origin: inference scheme of sir::boot_sir() (Hoff & Minhas SIR)")
			cli::cli_text("Valid replicates: {.val {object$n_valid}}/{.val {object$n_total}}")
			cli::cli_text("Standard errors available for: beta, a, b, variance components{ifelse(object$R > 0, ', U, V (Procrustes-aligned; per column where the singular values are separated) and the principal angles of the factor subspaces', '')}")
			.ae_print_withheld(object)
		})
	}
	invisible(object)
}

# sampler_describe() for the Bayesian MCMC fits (ame, lame): sampler, family, dimensions and
# dynamic options from the fit
#' @rdname sampler_describe
#' @method sampler_describe ame
#' @export
sampler_describe.ame <- function(object, verbose = TRUE) {
	if (verbose) {
		.lame_cli_stdout({
			cli::cli_h3("Estimator")
			cli::cli_text("Method: {.field Bayesian MCMC} (Gibbs / Metropolis-within-Gibbs)")
			cli::cli_text("Origin: Hoff (2005, 2009) additive and multiplicative effects sampler")
			# an amen fit (class "ame") carries no family or R: report them as unknown rather than as
			# normal and 0
			fam <- object[["family", exact = TRUE]]
			rnk <- object[["R", exact = TRUE]]
			if (is.null(fam) && is.null(rnk)) {
				cli::cli_text("Family: {.val unknown} (no {.pkg lame} metadata; looks like an {.pkg amen} fit)")
				cli::cli_text("Mode: {.val {object$mode %||% \"unipartite\"}}")
				cli::cli_text("Latent rank R: {.val {if (is.matrix(object$U)) ncol(object$U) else 0L}} (recovered from {.code fit$U})")
			} else {
				cli::cli_text("Family: {.val {fam %||% \"normal\"}}")
				cli::cli_text("Mode: {.val {object$mode %||% \"unipartite\"}}")
				cli::cli_text("Latent rank R: {.val {rnk %||% 0}}")
			}
			n_iter <- if (!is.null(object$BETA)) dim(object$BETA)[1L] else NA_integer_
			cli::cli_text("Stored draws: {.val {n_iter}}")
			cli::cli_text("Uncertainty available: yes (posterior draws)")
			has_ll <- !is.null(object$log_lik) || !is.null(object$log_lik_chunks)
			cli::cli_text("log_lik stored: {.val {has_ll}} (needed for {.fn loo::loo})")
		})
	}
	invisible(object)
}

#' @rdname sampler_describe
#' @method sampler_describe lame
#' @export
sampler_describe.lame <- function(object, verbose = TRUE) {
	if (verbose) {
		.lame_cli_stdout({
			cli::cli_h3("Estimator")
			cli::cli_text("Method: {.field Bayesian MCMC} (Gibbs / Metropolis-within-Gibbs)")
			cli::cli_text("Origin: Hoff (2005, 2009) AME sampler; Sewell & Chen (2015) and Durante & Dunson (2014) dynamic extensions")
			cli::cli_text("Family: {.val {object$family %||% \"normal\"}}")
			cli::cli_text("Mode: {.val {object$mode %||% \"unipartite\"}}")
			cli::cli_text("Periods T: {.val {length(object$EZ %||% object$YPM %||% list(NULL))}}")
			cli::cli_text("Latent rank R: {.val {object$R %||% 0}}")
			dyn_msg <- character()
			if (isTRUE(object$dynamic_beta)) {
				kind <- object$dynamic_beta_kind %||% "ar1"
				dyn_msg <- c(dyn_msg, paste0("dynamic_beta (", kind, ")"))
			}
			if (isTRUE(object$dynamic_ab)) dyn_msg <- c(dyn_msg, "dynamic_ab")
			if (isTRUE(object$dynamic_uv)) dyn_msg <- c(dyn_msg, "dynamic_uv")
			if (length(dyn_msg) == 0L) dyn_msg <- "none"
			cli::cli_text("Dynamic components: {.val {paste(dyn_msg, collapse = \", \")}}")
			n_iter <- if (!is.null(object$BETA)) dim(object$BETA)[1L] else NA_integer_
			cli::cli_text("Stored draws: {.val {n_iter}}")
			has_ll <- !is.null(object$log_lik) || !is.null(object$log_lik_chunks)
			cli::cli_text("log_lik stored: {.val {has_ll}} (needed for {.fn loo::loo})")
		})
	}
	invisible(object)
}

#' Update an `ame_als` / `lame_als` fit
#'
#' S3 method for \code{\link[stats]{update}} that re-fits an ALS fit
#' with modified arguments. Two routes:
#' \itemize{
#'   \item if \code{...} contains only the warm-start refit knobs
#'     (\code{Y_new}, \code{X_new}, \code{Z_new}, \code{max_iter},
#'     \code{tol}, \code{verbose}, \code{bias_correct}), it delegates to
#'     \code{\link{ame_als_refit}} for a fast warm-started refit (a
#'     dynamic fit is refit through its own fitter and comes back dynamic);
#'   \item otherwise (e.g. changing \code{R}, \code{family},
#'     \code{lambda}, \code{bootstrap}) it re-evaluates the original
#'     \code{object$call} with the overrides, matching the
#'     \code{\link{update.ame}} semantics on the MCMC side.
#' }
#' The split exists because \code{ame_als_refit()} is a genuinely
#' warm-started single-fit refit (faster, narrower argument set),
#' while changing \code{R} / \code{family} requires a cold-start
#' refit through \code{\link{ame_als}}.
#'
#' @param object A fitted \code{ame_als} / \code{lame_als} object.
#' @param ... Named arguments. See \code{\link{ame_als_refit}} for
#'   the warm-start argument list and \code{\link{ame_als}} for the
#'   cold-start argument list.
#' @param evaluate Logical: if \code{TRUE} (default), evaluate the
#'   updated call and return the new fit; if \code{FALSE}, return the
#'   unevaluated call (cold-start path only).
#'
#' @return A new fitted \code{ame_als} object, or (when
#'   \code{evaluate = FALSE} on the cold-start path) the modified
#'   call.
#'
#' @examples
#' \donttest{
#' Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
#' fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' # warm-start refit with the same arguments (fast)
#' fit_w <- update(fit, max_iter = 50)
#' # cold-start refit with a different R (full restart)
#' fit_R2 <- update(fit, R = 2)
#' }
#'
#' @method update ame_als
#' @export
update.ame_als <- function(object, ..., evaluate = TRUE) {
	if (!inherits(object, "ame_als")) {
		cli::cli_abort("{.arg object} must be an {.cls ame_als} object.")
	}
	dots <- list(...)
	# warm-start arguments of ame_als_refit() (excluding `object`)
	warm_args <- c(
		"Y_new", "X_new", "Z_new", "max_iter", "tol", "verbose",
		"bias_correct"
	)
	is_warm <- length(dots) > 0L && all(names(dots) %in% warm_args)
	if (is_warm || length(dots) == 0L) {
		# warm-start path: delegate to ame_als_refit()
		if (!isTRUE(evaluate)) {
			cli::cli_abort(c(
				"{.code evaluate = FALSE} is not supported on the warm-start refit path.",
				"i" = "Pass a non-warm argument (e.g. {.code R = 2}) to take the cold-start branch."
			))
		}
		return(do.call(ame_als_refit, c(list(object = object), dots)))
	}
	# cold-start path: re-evaluate the original call with the overrides
	if (is.null(object$call)) {
		cli::cli_abort(c(
			"No {.code call} stored on the fit object.",
			"i" = "Cold-start {.fn update} needs the original call. Refit explicitly with the new arguments."
		))
	}
	cl <- object$call
	nms <- names(dots)
	if (is.null(nms) || any(nms == "")) {
		cli::cli_abort("All arguments to {.fn update.ame_als} must be named.")
	}
	for (nm in nms) cl[[nm]] <- dots[[nm]]
	if (!isTRUE(evaluate)) return(cl)
	eval(cl, envir = parent.frame())
}

#' @rdname update.ame_als
#' @method update lame_als
#' @export
update.lame_als <- function(object, ..., evaluate = TRUE) {
	update.ame_als(object, ..., evaluate = evaluate)
}
