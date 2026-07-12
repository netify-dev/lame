# ame_als_bootstrap.r
#
# bootstrap uncertainty quantification for the fast ame estimator
# (ame_als() / lame_als()).
#
# the point estimator and this bootstrap-based inference are adapted from the
# social influence regression (sir) work of hoff & minhas --- the iterative
# block coordinate descent estimator and the block/parametric bootstrap of
# sir::boot_sir(). see minhas & hoff (2025), "decomposing network dynamics:
# social influence regression", political analysis. this is a port and
# adaptation, not original lame methodology.

# --------------------------------------------------------------------------
# internal helpers
# --------------------------------------------------------------------------

# draw a fresh additive contribution outer(a*, b*) for the full parametric
# bootstrap. the additive effects are random effects with dispersion
# sab = [[va, cab], [cab, vb]]; regenerating them per replicate (rather than
# holding the fitted a, b fixed) is what gives the intercept and node-covariate
# standard errors their actor-level sampling-variability component.
.ae_boot_draw_ab <- function(object) {
	nr <- object$dims$n_row; nc <- object$dims$n_col
	va  <- as.numeric(object$VC["va"])
	vb  <- as.numeric(object$VC["vb"])
	cab <- as.numeric(object$VC["cab"])
	if (!is.finite(va))  va  <- 0
	if (!is.finite(vb))  vb  <- 0
	if (!is.finite(cab)) cab <- 0
	if (object$symmetric) {
		a <- stats::rnorm(nr, 0, sqrt(max(va, 0)))
		return(outer(a, a, "+"))
	}
	if (identical(object$mode, "bipartite")) {
		a <- stats::rnorm(nr, 0, sqrt(max(va, 0)))
		b <- stats::rnorm(nc, 0, sqrt(max(vb, 0)))
		return(outer(a, b, "+"))
	}
	# unipartite directed: joint (a_i, b_i) draw preserving the cab covariance
	Sab <- matrix(c(va, cab, cab, vb), 2, 2)
	eg  <- eigen((Sab + t(Sab)) / 2, symmetric = TRUE)
	eg$values[eg$values < 0] <- 0
	Sroot <- eg$vectors %*% (sqrt(eg$values) * t(eg$vectors))
	ab <- matrix(stats::rnorm(2 * nr), nr, 2) %*% Sroot
	outer(ab[, 1], ab[, 2], "+")
}

# orthogonal procrustes rotation: returns rrot minimizing || m_src rrot - m_tgt ||
.ae_procrustes <- function(M_src, M_tgt) {
	sv <- tryCatch(svd(crossprod(M_src, M_tgt)), error = function(e) NULL)
	if (is.null(sv)) return(diag(ncol(M_src)))
	sv$u %*% t(sv$v)
}

# build a minimal "prep" list (the structure .ae_assemble() consumes) from
# canonical arrays --- used by ame_als_refit()
.ae_prep_from_canonical <- function(Yarr, X, W_row, W_col, mode, symmetric,
                                    row_names, col_names, time_names) {
	d <- dim(Yarr)
	xn <- if (dim(X)[3] > 0 && !is.null(dimnames(X))) {
		dimnames(X)[[3]]
	} else {
		character(0)
	}
	list(Yarr = Yarr, X = X, p_dyad = dim(X)[3], x_names_dyad = xn,
	     W_row = W_row, W_col = W_col,
	     nr = d[1], nc = d[2], Tt = d[3],
	     bip = identical(mode, "bipartite"),
	     row_names = row_names, col_names = col_names, time_names = time_names)
}

# --------------------------------------------------------------------------
# warm-start refit
# --------------------------------------------------------------------------

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
#'   fewer are needed than for a cold start).
#' @param tol convergence tolerance (default \code{1e-5}).
#' @param verbose logical; print progress (default \code{FALSE}).
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
ame_als_refit <- function(object, Y_new = NULL, X_new = NULL, Z_new = NULL,
                              max_iter = 30, tol = 1e-5, verbose = FALSE) {
	if (!inherits(object, "ame_als")) {
		cli::cli_abort("{.arg object} must be an {.cls ame_als} object.")
	}
	Yarr <- if (is.null(Y_new)) object$Y else Y_new
	X    <- if (is.null(X_new)) object$X else X_new
	# accept a plain matrix for a cross-sectional refit (canonical form is 3d)
	if (is.matrix(Yarr)) {
		Yarr <- array(Yarr, dim = c(nrow(Yarr), ncol(Yarr), 1L))
	}
	if (is.matrix(X)) {
		X <- array(X, dim = c(nrow(X), ncol(X), 1L, 1L))
	}
	if (!is.null(Z_new) && is.matrix(Z_new)) {
		Z_new <- array(Z_new, dim = c(nrow(Z_new), ncol(Z_new), 1L))
	}

	prep <- .ae_prep_from_canonical(Yarr, X, object$W_row, object$W_col,
	                                object$mode, object$symmetric,
	                                object$row_names, object$col_names,
	                                object$time_names)

	# the bcd's init$beta is the dyadic-covariate coefficient vector only
	# (node-covariate coefficients are recovered separately in .ae_assemble);
	# object$beta concatenates dyad + row + col, so pass just the dyadic head
	init <- list(mu = object$mu,
	             beta = object$beta[seq_len(prep$p_dyad)],
	             a = object$a, b = object$b,
	             U = object$U, V = object$V, L = object$L)

	lrm <- .ae_default(object$lowrank_method, "mm")
	nnm <- .ae_default(object$non_normal_method, "transform")
	lsv <- .ae_default(object$linear_solver, "eigen")

	if (!is.null(Z_new)) {
		# working-scale refit: z_new is the (gaussian) working response itself,
		# so a single gaussian bcd pass replaces the family transform / irls
		Z <- Z_new
		if (object$symmetric) for (t in seq_len(prep$Tt)) {
			Z[, , t] <- .ae_symmetrize(Z[, , t])
		}
		fit <- .ae_bcd_fit(Z, X, object$R, object$symmetric, max_iter, tol,
		                   init = init, lowrank_method = lrm,
		                   linear_solver = lsv)
		out <- .ae_assemble(fit, prep, object$family, object$mode,
		                    object$symmetric, object$R, object$longitudinal,
		                    call = match.call(), Zwork = Z)
		out$lowrank_method <- lrm
		out$meta$lowrank_method <- lrm
		out$non_normal_method <- nnm
		out$link <- object$link
		out$linear_solver <- lsv
		out$obs_weights <- NULL
		if (verbose) {
			cli::cli_text("Warm-start refit (working scale): converged = {.val {fit$converged}} ({fit$iterations} iter)")
		}
		return(out)
	}

	if (identical(nnm, "irls")) {
		# refit an irls fit by re-running the reweighting loop, warm-started
		ir <- .ae_irls_run(prep, object$family, object$link, object$R,
		                   object$symmetric, max_iter, tol,
		                   lowrank_method = lrm, init = init,
		                   linear_solver = lsv)
		fit <- ir$fit
		out <- .ae_assemble(fit, prep, object$family, object$mode,
		                    object$symmetric, object$R, object$longitudinal,
		                    call = match.call(),
		                    Zwork = ir$Z, link = object$link, Wobs = ir$W)
	} else {
		Z <- .ae_working_response(Yarr, object$family)
		if (object$symmetric) for (t in seq_len(prep$Tt)) {
			Z[, , t] <- .ae_symmetrize(Z[, , t])
		}
		fit <- .ae_bcd_fit(Z, X, object$R, object$symmetric, max_iter, tol,
		                   init = init, lowrank_method = lrm,
		                   linear_solver = lsv)
		out <- .ae_assemble(fit, prep, object$family, object$mode,
		                    object$symmetric, object$R, object$longitudinal,
		                    call = match.call())
	}

	if (verbose) {
		cli::cli_text("Warm-start refit: converged = {.val {fit$converged}} ({fit$iterations} iter)")
	}

	out$lowrank_method <- lrm
	out$meta$lowrank_method <- lrm
	out$non_normal_method <- nnm
	out$link <- object$link
	out$linear_solver <- lsv
	out$obs_weights <- if (identical(nnm, "irls")) ir$W else NULL
	out
}

# --------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------

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
#'     fits only (requires \code{T > 1}); with few slices its intervals are
#'     necessarily coarse.}
#' }
#'
#' Every replicate is refit by \code{\link{ame_als_refit}}, warm-started
#' from the original point estimate so that replicates do not drift to
#' different local optima. Replicates that error or return non-finite values
#' are dropped and counted (\code{n_valid} / \code{n_total}). Standard errors
#' are replicate column standard deviations; confidence intervals use the
#' percentile method.
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
#' the effect of alignment can be inspected. The alignment is well-determined
#' only when the multiplicative singular values are well separated; with
#' near-equal singular values the per-column \code{U}/\code{V} standard errors
#' reflect an unstable rotation, and the subspace --- or, for symmetric models,
#' the eigenvalues \code{L} --- is the summary of record. For a symmetric model
#' the eigenvalues \code{L} of the multiplicative term are aggregated and
#' reported as the primary multiplicative-uncertainty summary.
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
#' dynamic \code{\link{lame}} there); with only a few time slices it is
#' necessarily coarse, its intervals correspondingly imprecise and somewhat
#' anti-conservative, so \code{parametric} is the default. A variance component
#' fit with \code{R > 0} can still carry a small residual parametric-bootstrap
#' bias (the low-rank refit re-absorbs simulated noise); \code{summary()} flags
#' any point estimate that falls outside its interval. A binary IRLS fit can
#' carry a finite-sample (incidental-parameters) bias in the point estimator
#' itself, which the bootstrap reproduces rather than removes. The MCMC
#' \code{\link{ame}} / \code{\link{lame}} path gives posterior summaries when
#' that is the target.
#'
#' @param object an \code{ame_als} object from \code{\link{ame_als}}
#'   or \code{\link{lame_als}}.
#' @param R integer number of bootstrap replicates (default 200).
#' @param type \code{"parametric"} (default) or \code{"block"}; see Description.
#' @param block_length block length for the block bootstrap: \code{1} (default)
#'   resamples time slices independently; an integer \code{> 1} draws contiguous
#'   blocks of that many slices (a moving-block bootstrap), appropriate when the
#'   series has short-run temporal dependence. Capped at \code{T}; ignored for
#'   the parametric bootstrap.
#' @param seed optional integer random seed.
#' @param verbose logical; print progress (default \code{TRUE}).
#'
#' @return An object of class \code{"boot_ame"} with components including
#'   \code{coefs} (replicate intercept + regression coefficients), \code{se},
#'   \code{ci_lo}, \code{ci_hi}, \code{point_est}, \code{param_names};
#'   \code{vc_*} for the variance components; \code{a_coefs}, \code{b_coefs},
#'   \code{se_a}, \code{se_b}; \code{U_aligned}, \code{V_aligned},
#'   \code{U_raw}, \code{V_raw}, \code{se_U}, \code{se_V} (when \code{R > 0});
#'   and \code{n_valid}, \code{n_total}, \code{type}, \code{family}.
#'
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}. The block and parametric
#' bootstrap follow the inference scheme of \code{sir::boot_sir()} developed
#' for the SIR estimator.
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
#' Y <- replicate(6, { m <- matrix(rnorm(225), 15, 15); diag(m) <- NA; m },
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
ame_als_bootstrap <- function(object, R = 200,
                                  type = c("parametric", "block"),
                                  block_length = 1, seed = NULL,
                                  verbose = TRUE) {

	if (!inherits(object, "ame_als")) {
		cli::cli_abort("{.arg object} must be an {.cls ame_als} object from {.fn ame_als} or {.fn lame_als}.")
	}
	type <- match.arg(type)
	if (length(block_length) != 1L || !is.finite(block_length) ||
	    block_length < 1) {
		cli::cli_abort("{.arg block_length} must be a single positive integer.")
	}
	block_length <- as.integer(block_length)
	# seed the bootstrap without leaving the global rng stream mutated, so a
	# scripted analysis is not silently perturbed (parity with ame_als())
	if (!is.null(seed)) {
		if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			old_seed <- get(".Random.seed", envir = globalenv())
			on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
			        add = TRUE)
		} else {
			on.exit(if (exists(".Random.seed", envir = globalenv(),
			                   inherits = FALSE))
				rm(".Random.seed", envir = globalenv()), add = TRUE)
		}
		set.seed(seed)
	}
	# block_length capped at t; effective length used by the moving-block resampler
	blk <- min(block_length, object$n_time)

	Tt   <- object$n_time
	Rlat <- object$R
	nr   <- object$dims$n_row
	nc   <- object$dims$n_col
	family <- object$family
	bip <- identical(object$mode, "bipartite")
	rho_hat <- as.numeric(object$VC["rho"])
	if (!is.finite(rho_hat)) rho_hat <- 0
	s2_hat  <- object$s2
	nnm <- .ae_default(object$non_normal_method, "transform")
	# a non-normal `transform` fit has no calibrated generative model: its
	# estimator is a gaussian fit to a fixed transformed response, so the
	# parametric bootstrap resamples on that working scale (z = ez + noise) and
	# refits directly. irls and normal fits have a genuine response model and
	# are simulated on the response scale.
	working_boot <- identical(type, "parametric") &&
		family %in% c("binary", "poisson", "ordinal") &&
		!identical(nnm, "irls")

	# simulated residual variance: object$s2 (ve) is already df-corrected in
	# .ae_assemble (sse / (n_obs - df_mean), via .ae_df_mean), so simulating
	# with it directly centres the refits on the point estimate -- each
	# refit's flexible mean re-absorbs its nominal share of the fresh noise,
	# which the df correction has already added back.
	s2_sim <- s2_hat

	# the fitted additive contribution, swapped out per replicate for a fresh
	# draw in the full parametric bootstrap
	ab_orig <- outer(as.numeric(object$a), as.numeric(object$b), "+")

	if (type == "block" && Tt < 2L) {
		cli::cli_abort(c(
			"Block bootstrap resamples time slices and needs {.arg T} > 1.",
			"i" = "This fit has a single time slice; use {.code type = \"parametric\"}.",
			"i" = "(A nodal/dyadic resample is not offered: it would duplicate actors and not yield a valid network.)"))
	}
	# moving-block bootstrap with too few blocks per replicate is degenerate
	if (type == "block" && blk > 1L && ceiling(Tt / blk) < 3L) {
		cli::cli_warn(c(
			"Moving-block bootstrap yields fewer than 3 blocks per replicate.",
			"i" = "{.arg block_length} = {block_length} is large relative to {Tt} time slices; standard errors will be unstable."))
	}

	# headline parameters: intercept + regression coefficients
	point_coef <- object$coefficients
	param_names <- names(point_coef)
	n_par <- length(point_coef)

	boot_coefs <- matrix(NA_real_, R, n_par)
	colnames(boot_coefs) <- param_names
	vc_coefs <- matrix(NA_real_, R, 5)
	colnames(vc_coefs) <- names(object$VC)
	a_coefs <- matrix(NA_real_, R, nr)
	b_coefs <- matrix(NA_real_, R, nc)
	boot_converged <- logical(R)               # did the warm-start refit converge?
	U_raw <- V_raw <- U_aln <- V_aln <- L_coefs <- uv_ok <- NULL
	if (Rlat > 0) {
		U_raw <- matrix(NA_real_, R, nr * Rlat)
		U_aln <- matrix(NA_real_, R, nr * Rlat)
		V_raw <- matrix(NA_real_, R, nc * Rlat)
		V_aln <- matrix(NA_real_, R, nc * Rlat)
		uv_ok <- logical(R)                    # did u,v come back finite?
		if (object$symmetric) L_coefs <- matrix(NA_real_, R, Rlat)
	}

	if (verbose) {
		cli::cli_h3("Bootstrap for fast AME fit")
		cli::cli_text("Type: {.field {type}} | Replicates: {.val {R}} | Family: {.field {family}}")
		cli::cli_text("Network: {nr} x {nc}, {Tt} time slice{?s}, R = {Rlat}")
	}

	for (b in seq_len(R)) {
		if (verbose && b %% 25 == 0) cli::cli_inform("Bootstrap {.val {b}}/{.val {R}}")

		Y_b <- NULL; X_b <- object$X; Z_b <- NULL

		if (type == "block") {
			if (blk <= 1L) {
				t_idx <- sample.int(Tt, Tt, replace = TRUE)
			} else {
				# moving-block resample: contiguous blocks of `blk` slices,
				# wrapped circularly so every slice has equal inclusion
				# probability (a plain non-circular resample under-samples the
				# boundary slices)
				n_blk  <- ceiling(Tt / blk)
				starts <- sample.int(Tt, n_blk, replace = TRUE)
				t_idx  <- unlist(lapply(starts,
					function(s) ((s - 1L + seq_len(blk) - 1L) %% Tt) + 1L))[seq_len(Tt)]
			}
			Y_b <- object$Y[, , t_idx, drop = FALSE]
			X_b <- object$X[, , , t_idx, drop = FALSE]
		} else if (working_boot) {
			# non-normal transform: resample on the gaussian working scale.
			# object$ez is the fitted working-response predictor; swap the
			# additive contribution for a fresh draw, then add gaussian
			# residual noise with the fit's variance and reciprocity.
			Z_b <- array(NA_real_, dim = c(nr, nc, Tt))
			ab_new <- .ae_boot_draw_ab(object)
			for (t in seq_len(Tt)) {
				ez <- object$EZ[[t]] - ab_orig + ab_new
				ez[!is.finite(ez)] <- 0
				zt <- as.matrix(simY_nrm(ez, rho_hat, s2_sim))
				storage.mode(zt) <- "double"
				if (object$symmetric) {
					zt[lower.tri(zt)] <- 0
					zt <- zt + t(zt)
				}
				zt[is.na(object$Y[, , t])] <- NA_real_
				Z_b[, , t] <- zt
			}
		} else {
			# parametric, response scale: simulate fresh outcomes from the
			# calibrated fitted model (normal, or an irls glm fit). the additive
			# random effects are regenerated per replicate; mu, beta and the
			# multiplicative term stay fixed.
			Y_b <- array(NA_real_, dim = c(nr, nc, Tt))
			ab_new <- .ae_boot_draw_ab(object)
			for (t in seq_len(Tt)) {
				ez <- object$EZ[[t]] - ab_orig + ab_new
				ez[!is.finite(ez)] <- 0        # ez diagonal is na for unipartite
				yt <- if (family == "normal") {
					simY_nrm(ez, rho_hat, s2_sim)
				} else if (family == "poisson") {
					# irls poisson: ez is the log-mean linear predictor
					simY_pois(ez)
				} else {
					# irls binary: ez is the link-scale predictor; map to its
					# probit-equivalent so the simulator (probit threshold)
					# reproduces the fitted rate for either link
					p   <- .ae_link_inverse(ez, family, object$link)
					p[!is.finite(p)] <- 0.5
					eta <- stats::qnorm(pmin(pmax(p, 1e-6), 1 - 1e-6))
					simY_bin(eta, rho_hat)
				}
				yt <- as.matrix(yt)
				storage.mode(yt) <- "double"
				# symmetric model: mirror the upper triangle so the simulated
				# network is symmetric with the correct per-tie variance
				# (otherwise the refit averages two independent draws and the
				# bootstrap ses are deflated)
				if (object$symmetric) {
					yt[lower.tri(yt)] <- 0
					yt <- yt + t(yt)
				}
				# re-impose the original observed-cell pattern: cells that were
				# missing in the data (the unipartite diagonal, and any genuinely
				# unobserved dyads) must stay missing, so each replicate carries
				# the same information content as the data. simulating those cells
				# would give the refit more data than the original fit and bias
				# the bootstrap standard errors downward.
				yt[is.na(object$Y[, , t])] <- NA_real_
				Y_b[, , t] <- yt
			}
		}

		rep_ok <- tryCatch({
			# suppress per-replicate warnings (e.g. a rank-deficient design):
			# they would repeat once per replicate; the original fit already
			# surfaced them, and the bootstrap reports its own diagnostics
			fit_b <- suppressWarnings(if (!is.null(Z_b)) {
				ame_als_refit(object, Y_new = object$Y, Z_new = Z_b,
				                  max_iter = 50, tol = 1e-5, verbose = FALSE)
			} else {
				ame_als_refit(object, Y_new = Y_b, X_new = X_b,
				                  max_iter = 50, tol = 1e-5, verbose = FALSE)
			})
			cf <- fit_b$coefficients
			vc <- fit_b$VC
			ab_ok <- all(is.finite(fit_b$a)) && all(is.finite(fit_b$b))
			if (all(is.finite(cf)) && ab_ok) {
				boot_coefs[b, ] <- cf
				vc_coefs[b, ]   <- vc
				a_coefs[b, ]    <- fit_b$a
				b_coefs[b, ]    <- fit_b$b
				boot_converged[b] <- isTRUE(fit_b$converged)
				if (Rlat > 0 && !is.null(fit_b$U) && !is.null(fit_b$V) &&
				    all(is.finite(fit_b$U)) && all(is.finite(fit_b$V))) {
					U_raw[b, ] <- c(fit_b$U)
					V_raw[b, ] <- c(fit_b$V)
					# procrustes-align replicate (u, v) to the original fit
					src <- rbind(fit_b$U, fit_b$V)
					tgt <- rbind(object$U, object$V)
					Rrot <- .ae_procrustes(src, tgt)
					U_aln[b, ] <- c(fit_b$U %*% Rrot)
					V_aln[b, ] <- c(fit_b$V %*% Rrot)
					if (object$symmetric && !is.null(fit_b$L)) {
						L_coefs[b, ] <- fit_b$L
					}
					uv_ok[b] <- TRUE
				}
				TRUE
			} else {
				FALSE
			}
		}, error = function(e) FALSE)
		if (!isTRUE(rep_ok)) next
	}

	valid <- stats::complete.cases(boot_coefs)
	n_valid <- sum(valid)
	n_nonconverged <- sum(valid & !boot_converged)

	if (n_valid == 0) {
		cli::cli_abort(c(
			"All {R} bootstrap replicates failed to refit.",
			"i" = "The model may be misspecified or the resampled data too sparse.",
			"i" = "Inspect the fit, refit with a smaller latent dimension, or use {.fn ame}/{.fn lame}."))
	}
	if (n_valid < 10) {
		# distinguish "few replicates requested" from "many requested but
		# most failed" so the warning matches the actual situation.
		if (n_valid == R) {
			cli::cli_warn(c(
				"Only {.val {n_valid}} bootstrap replicate{?s} total -- standard errors are unreliable.",
				"i" = "Increase {.arg R} to >= 50 (>= 200 recommended)."))
		} else {
			cli::cli_warn(c(
				"Only {.val {n_valid}} of {.val {R}} bootstrap replicates succeeded.",
				"i" = "Standard errors are unreliable; increase the replicate count or simplify the model."))
		}
	} else if (n_valid < R && verbose) {
		cli::cli_alert_info("{R - n_valid} of {R} replicate{?s} failed and were dropped.")
	}
	if (n_nonconverged > 0) {
		msg <- "{n_nonconverged} of {n_valid} refit{?s} did not fully converge (kept; warm-started and near-optimal)."
		if (n_nonconverged > 0.10 * n_valid) {
			# > 10% non-convergence: surface a real warning + actionable hint
			frac <- round(100 * n_nonconverged / n_valid)
			cli::cli_warn(c(
				"!" = msg,
				"i" = "{frac}% non-convergence: bumping {.arg max_iter} (default 200) often helps.",
				"i" = "Per-replicate flags are at {.code fit$bootstrap$boot_converged}; intervals on non-converged replicates may be biased.",
				"i" = "If the rate stays high, simplify the model ({.arg R} = 0) or use {.fn ame} (MCMC) for calibrated intervals."))
		} else if (verbose) {
			cli::cli_alert_info(msg)
		}
	}

	col_sd <- function(M, v) apply(M[v, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
	col_q  <- function(M, v, pr) apply(M[v, , drop = FALSE], 2, stats::quantile,
	                                   probs = pr, na.rm = TRUE)

	ci_lo_vec <- col_q(boot_coefs, valid, 0.025)
	ci_hi_vec <- col_q(boot_coefs, valid, 0.975)
	# flag coefficients whose point estimate falls outside their own ci
	outside <- which(is.finite(point_coef) &
	                 is.finite(ci_lo_vec) &
	                 is.finite(ci_hi_vec) &
	                 (point_coef < ci_lo_vec | point_coef > ci_hi_vec))
	if (length(outside) > 0L) {
		bad <- param_names[outside]
		cli::cli_warn(c(
			"!" = "Bootstrap CIs do not contain the point estimate for: {.code {bad}}.",
			"i" = "This is a sign the warm-start refits drifted to a different local optimum",
			"i" = "(common with non-normal IRLS at {.code R > 0} on small/sparse data).",
			"i" = "These intervals are unreliable. Refit with {.code R = 0}, or use {.fn ame}/{.fn lame} (MCMC) for calibrated uncertainty."))
	}

	out <- list(
		coefs = boot_coefs, se = col_sd(boot_coefs, valid),
		ci_lo = ci_lo_vec, ci_hi = ci_hi_vec,
		point_est = point_coef, param_names = param_names,
		vc_coefs = vc_coefs, vc_se = col_sd(vc_coefs, valid),
		vc_ci_lo = col_q(vc_coefs, valid, 0.025),
		vc_ci_hi = col_q(vc_coefs, valid, 0.975),
		vc_point_est = object$VC,
		a_coefs = a_coefs, b_coefs = b_coefs,
		se_a = col_sd(a_coefs, valid), se_b = col_sd(b_coefs, valid),
		point_a = object$a, point_b = object$b,
		n_valid = n_valid, n_total = R, n_nonconverged = n_nonconverged,
		type = type, block_length = if (type == "block") blk else NA_integer_,
		family = family, mode = object$mode, symmetric = object$symmetric,
		R = Rlat,
		# enrich dims with the underlying actor names so confint() labels
		# intervals with real names instead of generic "a_a1...a_an".
		dims = c(object$dims,
		         list(row_names = names(object$a) %||% object$row_names,
		              col_names = names(object$b) %||% object$col_names)),
		n_time = Tt, call = match.call())

	if (Rlat > 0) {
		vuv <- valid & uv_ok
		out$n_valid_uv <- sum(vuv)
		if (out$n_valid_uv > 0) {
			out$U_raw <- U_raw; out$V_raw <- V_raw
			out$U_aligned <- U_aln; out$V_aligned <- V_aln
			out$se_U <- col_sd(U_aln, vuv); out$se_V <- col_sd(V_aln, vuv)
			out$se_U_raw <- col_sd(U_raw, vuv); out$se_V_raw <- col_sd(V_raw, vuv)
			out$point_U <- object$U; out$point_V <- object$V
			if (object$symmetric && !is.null(L_coefs)) {
				out$L_coefs <- L_coefs
				out$se_L <- col_sd(L_coefs, vuv)
				out$ci_L_lo <- col_q(L_coefs, vuv, 0.025)
				out$ci_L_hi <- col_q(L_coefs, vuv, 0.975)
				out$point_L <- object$L
			}
		}
	}

	class(out) <- "boot_ame"
	out
}

#' @rdname ame_als_bootstrap
#' @export
boot_ame <- function(object, R = 200, type = c("parametric", "block"),
                     block_length = 1, seed = NULL, verbose = TRUE) {
	ame_als_bootstrap(object, R = R, type = type,
	                      block_length = block_length, seed = seed,
	                      verbose = verbose)
}

# --------------------------------------------------------------------------
# s3 methods for boot_ame objects
# --------------------------------------------------------------------------

#' Print bootstrap results for a fast AME fit
#'
#' @param x a \code{boot_ame} object.
#' @param digits number of digits to display.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @method print boot_ame
#' @export
print.boot_ame <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
	nc_note <- if (isTRUE(x$n_nonconverged > 0)) {
		paste0(" (", x$n_nonconverged, " not fully converged)")
	} else {
		""
	}
	cli::cli_text("")
	cli::cli_text("{.strong Bootstrap results: fast AME estimator}")
	cli::cli_text("Type: {.field {x$type}} | Family: {.field {x$family}} | Replicates: {.val {x$n_valid}}/{.val {x$n_total}} valid{nc_note}")
	tab <- cbind(Estimate = x$point_est, `Boot SE` = x$se,
	             `2.5%` = x$ci_lo, `97.5%` = x$ci_hi)
	rownames(tab) <- x$param_names
	cli::cli_text("")
	print(round(tab, digits))
	cli::cli_text("")
	cli::cli_text(cli::col_grey("beta, a, b and variance components aggregate directly; U/V are Procrustes-aligned."))
	if (identical(x$type, "block")) {
		cli::cli_text(cli::col_yellow(
			"! Block bootstrap on a static-effects fit: SEs for the intercept, node-covariate coefficients and additive variance components (va, vb, cab) can severely under-cover because slice resampling does not perturb actor-level structure. Use {.code type = \"parametric\"} for those parameters."))
	}
	invisible(x)
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
	object <- x$boot
	cli::cli_h1("Bootstrap results: fast AME estimator")
	cli::cli_text("Type: {.field {object$type}} | Family: {.field {object$family}} | Mode: {.field {object$mode}}")
	if (identical(object$type, "block")) {
		cli::cli_text(cli::col_yellow(
			"! Block bootstrap on a static-effects fit: intervals for the intercept, node-covariate coefficients and additive variance components (va, vb, cab) can severely under-cover (slice resampling does not perturb actor-level structure). Prefer {.code type = \"parametric\"} for those parameters."))
	}
	if (identical(object$type, "block") && isTRUE(object$block_length > 1L)) {
		cli::cli_text("Moving-block resample: block length {.val {object$block_length}}")
	}
	cli::cli_text("Replicates: {.val {object$n_valid}} valid of {.val {object$n_total}} ({round(100 * object$n_valid / object$n_total, 1)}%)")
	if (isTRUE(object$n_nonconverged > 0)) {
		cli::cli_text("Refits not fully converged (kept, warm-started): {.val {object$n_nonconverged}}")
	}
	cli::cli_text("Network: {object$dims$n_row} x {object$dims$n_col}, {object$n_time} time slice{?s}, R = {object$R}")

	cli::cli_rule()
	cli::cli_text("{.strong Regression coefficients}")
	tab <- cbind(Estimate = object$point_est, `Boot SE` = object$se,
	             `2.5%` = object$ci_lo, `97.5%` = object$ci_hi)
	rownames(tab) <- object$param_names
	covers0 <- object$ci_lo <= 0 & object$ci_hi >= 0
	tab <- cbind(round(tab, 4), sig = ifelse(covers0, "", "*"))
	print(noquote(tab))
	cli::cli_text(cli::col_grey("* = 95% percentile CI excludes zero"))

	out_coef <- object$point_est < object$ci_lo | object$point_est > object$ci_hi
	if (any(out_coef, na.rm = TRUE)) {
		cli::cli_text(cli::col_yellow(
			"! point estimate outside its bootstrap CI for: {paste(object$param_names[which(out_coef)], collapse = ', ')} -- the replicate distribution is biased; prefer {.code confint(., ci_type = \"basic\")}."))
	}

	cli::cli_rule()
	cli::cli_text("{.strong Variance components}")
	vtab <- cbind(Estimate = object$vc_point_est, `Boot SE` = object$vc_se,
	              `2.5%` = object$vc_ci_lo, `97.5%` = object$vc_ci_hi)
	print(round(vtab, 4))
	out_vc <- object$vc_point_est < object$vc_ci_lo |
		object$vc_point_est > object$vc_ci_hi
	if (any(out_vc, na.rm = TRUE)) {
		cli::cli_text(cli::col_yellow(
			"! point estimate outside its bootstrap CI for: {paste(names(object$vc_point_est)[which(out_vc)], collapse = ', ')} -- a residual parametric-bootstrap bias; treat that interval as approximate."))
	}

	cli::cli_rule()
	cli::cli_text("{.strong Additive effects} (rotation-invariant, aggregated directly)")
	cli::cli_text("Sender a: mean |SE| = {round(mean(object$se_a, na.rm = TRUE), 4)}")
	cli::cli_text("Receiver b: mean |SE| = {round(mean(object$se_b, na.rm = TRUE), 4)}")

	if (object$R > 0) {
		if (is.null(object$se_U)) {
			cli::cli_text("{.strong Multiplicative factors}: no replicate produced finite U/V.")
		} else if (object$symmetric) {
			cli::cli_text("{.strong Multiplicative factors} (symmetric eigenmodel; {object$n_valid_uv} replicate{?s})")
			cli::cli_text("Eigenvalues L  point: {paste(sprintf('%.3f', object$point_L), collapse = ', ')}")
			cli::cli_text("Eigenvalues L  boot SE: {paste(sprintf('%.3f', object$se_L), collapse = ', ')}")
			cli::cli_text(cli::col_grey("L (eigenvalues) is the primary multiplicative-uncertainty summary; U columns are sign/eigenspace-identified only."))
		} else {
			cli::cli_text("{.strong Multiplicative factors} (Procrustes-aligned; {object$n_valid_uv} replicate{?s})")
			cli::cli_text("U: mean SE = {round(mean(object$se_U, na.rm = TRUE), 4)} (raw, unaligned: {round(mean(object$se_U_raw, na.rm = TRUE), 4)})")
			cli::cli_text("V: mean SE = {round(mean(object$se_V, na.rm = TRUE), 4)} (raw, unaligned: {round(mean(object$se_V_raw, na.rm = TRUE), 4)})")
		}
	}
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
#'   replicate quantiles directly and matches the intervals reported by
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
#' @method confint boot_ame
#' @export
confint.boot_ame <- function(object, parm = NULL, level = 0.95,
                             ci_type = c("percentile", "basic"),
                             which = c("all", "beta", "vc", "a", "b", "U", "V"),
                             ...) {
	ci_type <- match.arg(ci_type)
	which   <- match.arg(which)
	if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
		cli::cli_abort("{.arg level} must be a single number in (0, 1).")
	}
	a  <- (1 - level) / 2
	pr <- c(a, 1 - a)
	pct <- paste0(format(100 * pr, trim = TRUE), "%")

	# helper: build a ci matrix from a draws-matrix + point-estimate vector.
	# quantile() is computed column-wise with na.rm so a column that's na in
	# every replicate (e.g. cab on a bipartite fit) just returns na quantiles
	# instead of removing every row from the matrix.
	.boot_ci <- function(draws, point_est, names_, prefix = "") {
		if (is.null(draws) || nrow(draws) == 0L || ncol(draws) == 0L) {
			return(matrix(numeric(0), 0, 2, dimnames = list(NULL, pct)))
		}
		q <- t(apply(draws, 2, function(col) {
			col <- col[is.finite(col)]
			if (length(col) == 0L) return(c(NA_real_, NA_real_))
			stats::quantile(col, probs = pr)
		}))
		out <- if (ci_type == "percentile") {
			q
		} else {
			# basic bootstrap interval: 2 * theta_hat - upper/lower percentile
			cbind(2 * point_est - q[, 2], 2 * point_est - q[, 1])
		}
		rownames(out) <- if (nzchar(prefix)) paste0(prefix, names_) else names_
		colnames(out) <- pct
		out
	}

	parts <- list()
	if (which %in% c("all", "beta")) {
		parts$beta <- .boot_ci(object$coefs, object$point_est,
		                       object$param_names)
	}
	if (which %in% c("all", "vc") && !is.null(object$vc_coefs)) {
		vc <- .boot_ci(object$vc_coefs, object$vc_point_est,
		               colnames(object$vc_coefs))
		# variance components are non-negative by construction; the basic
		# bootstrap (2 theta_hat - q) can push the lower endpoint below 0
		# near a boundary. clip to [0, inf) for va/vb/ve and document.
		non_neg <- intersect(c("va", "vb", "ve"), rownames(vc))
		if (length(non_neg) > 0L) {
			vc[non_neg, 1] <- pmax(vc[non_neg, 1], 0, na.rm = TRUE)
			vc[non_neg, 2] <- pmax(vc[non_neg, 2], 0, na.rm = TRUE)
		}
		parts$vc <- vc
	}
	if (which %in% c("all", "a") && !is.null(object$a_coefs)) {
		nm_a <- if (!is.null(object$dims$row_names)) object$dims$row_names
		        else colnames(object$a_coefs) %||% paste0("a", seq_len(ncol(object$a_coefs)))
		parts$a <- .boot_ci(object$a_coefs, object$point_a,
		                    nm_a, prefix = "a_")
	}
	if (which %in% c("all", "b") && !is.null(object$b_coefs) &&
	    nrow(object$b_coefs) > 0L && ncol(object$b_coefs) > 0L) {
		nm_b <- if (!is.null(object$dims$col_names)) object$dims$col_names
		        else colnames(object$b_coefs) %||% paste0("b", seq_len(ncol(object$b_coefs)))
		parts$b <- .boot_ci(object$b_coefs, object$point_b,
		                    nm_b, prefix = "b_")
	}
	# u/v are stored flat: u_aligned is an [r_boot, n_row * r_lat] matrix where
	# each row is c(point_u) for that replicate (column-major flatten).
	.uv_names <- function(actors, Rlat, label) {
		# column-major flatten matches c(matrix(...)): position i + (r-1)*n
		as.vector(outer(actors, paste0("d", seq_len(Rlat)),
		                FUN = function(a, d) paste0(label, "_", d, "_", a)))
	}
	if (which %in% c("all", "U") && !is.null(object$U_aligned) &&
	    is.matrix(object$U_aligned) && nrow(object$U_aligned) > 0L) {
		nr   <- nrow(object$point_U)
		Rlat <- ncol(object$point_U)
		actors <- if (!is.null(rownames(object$point_U))) rownames(object$point_U)
		          else paste0("a", seq_len(nr))
		parts$U <- .boot_ci(object$U_aligned,
		                    as.vector(object$point_U),
		                    .uv_names(actors, Rlat, "U"))
	}
	if (which %in% c("all", "V") && !is.null(object$V_aligned) &&
	    is.matrix(object$V_aligned) && nrow(object$V_aligned) > 0L &&
	    ncol(object$V_aligned) > 0L) {
		nc   <- nrow(object$point_V)
		Rlat <- ncol(object$point_V)
		actors <- if (!is.null(rownames(object$point_V))) rownames(object$point_V)
		          else paste0("a", seq_len(nc))
		parts$V <- .boot_ci(object$V_aligned,
		                    as.vector(object$point_V),
		                    .uv_names(actors, Rlat, "V"))
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
					"i" = "Available: {.val {head(rownames(ci), 20)}}{cli::qty(max(0, nrow(ci) - 20))}{?/, ...}."))
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
#' @return A covariance matrix over the intercept and regression coefficients.
#' @method vcov boot_ame
#' @export
vcov.boot_ame <- function(object, ...) {
	valid <- stats::complete.cases(object$coefs)
	V <- stats::cov(object$coefs[valid, , drop = FALSE])
	dimnames(V) <- list(object$param_names, object$param_names)
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
		"i" = "A bootstrap aggregates per-replicate fits; call {.fn fitted} on the underlying {.cls ame_als} fit."))
}

#' @rdname boot_ame-no-fitted
#' @method residuals boot_ame
#' @export
residuals.boot_ame <- function(object, ...) {
	cli::cli_abort(c(
		"{.fn residuals} is not defined for a {.cls boot_ame} object.",
		"i" = "Call {.fn residuals} on the underlying {.cls ame_als} fit."))
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

# --------------------------------------------------------------------------
# sampler_describe generic
# --------------------------------------------------------------------------

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
#' @export
sampler_describe <- function(object, verbose = TRUE) {
	UseMethod("sampler_describe", object)
}

#' @rdname sampler_describe
#' @method sampler_describe ame_als
#' @export
sampler_describe.ame_als <- function(object, verbose = TRUE) {
	if (verbose) {
		cli::cli_h3("Estimator")
		cli::cli_text("Method: {.field iterative block coordinate descent} (fast, MCMC-free)")
		cli::cli_text("Origin: Hoff & Minhas Social Influence Regression estimator (sir::sir_alsfit())")
		cli::cli_text("Uncertainty available: no (point estimate only)")
		cli::cli_text("For standard errors, run {.fn ame_als_bootstrap}.")
		cli::cli_text(cli::col_grey(object$meta$approximation_note))
	}
	invisible(object)
}

#' @rdname sampler_describe
#' @method sampler_describe boot_ame
#' @export
sampler_describe.boot_ame <- function(object, verbose = TRUE) {
	if (verbose) {
		cli::cli_h3("Estimator")
		cli::cli_text("Method: {.field {object$type}} bootstrap of the iterative block coordinate descent estimator")
		cli::cli_text("Origin: inference scheme of sir::boot_sir() (Hoff & Minhas SIR)")
		cli::cli_text("Valid replicates: {.val {object$n_valid}}/{.val {object$n_total}}")
		cli::cli_text("Standard errors available for: beta, a, b, variance components{ifelse(object$R > 0, ', U, V (Procrustes-aligned)', '')}")
	}
	invisible(object)
}

# describe() for the bayesian mcmc path (ame and lame fits). prints the
# sampler / family / dimensions / dynamic options surfaced on the fit
# object. shares the same generic so `sampler_describe(fit)` works
# regardless of which constructor produced the fit.
#' @rdname sampler_describe
#' @method sampler_describe ame
#' @export
sampler_describe.ame <- function(object, verbose = TRUE) {
	if (verbose) {
		cli::cli_h3("Estimator")
		cli::cli_text("Method: {.field Bayesian MCMC} (Gibbs / Metropolis-within-Gibbs)")
		cli::cli_text("Origin: Hoff (2005, 2009) additive and multiplicative effects sampler")
		cli::cli_text("Family: {.val {object$family %||% \"normal\"}}")
		cli::cli_text("Mode: {.val {object$mode %||% \"unipartite\"}}")
		cli::cli_text("Latent rank R: {.val {object$R %||% 0}}")
		n_iter <- if (!is.null(object$BETA)) dim(object$BETA)[1L] else NA_integer_
		cli::cli_text("Stored draws: {.val {n_iter}}")
		cli::cli_text("Uncertainty available: yes (posterior draws)")
		has_ll <- !is.null(object$log_lik) || !is.null(object$log_lik_chunks)
		cli::cli_text("log_lik stored: {.val {has_ll}} (needed for {.fn loo::loo})")
	}
	invisible(object)
}

#' @rdname sampler_describe
#' @method sampler_describe lame
#' @export
sampler_describe.lame <- function(object, verbose = TRUE) {
	if (verbose) {
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
#'     \code{tol}, \code{verbose}), it delegates to
#'     \code{\link{ame_als_refit}} for a fast warm-started refit;
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
	# warm-start knobs exposed by ame_als_refit() (excluding `object`)
	warm_args <- c("Y_new", "X_new", "Z_new", "max_iter", "tol", "verbose")
	is_warm <- length(dots) > 0L && all(names(dots) %in% warm_args)
	if (is_warm || length(dots) == 0L) {
		# warm-start path: delegate to ame_als_refit
		if (!isTRUE(evaluate)) {
			cli::cli_abort(c(
				"{.code evaluate = FALSE} is not supported on the warm-start refit path.",
				"i" = "Pass a non-warm argument (e.g. {.code R = 2}) to take the cold-start branch."))
		}
		return(do.call(ame_als_refit, c(list(object = object), dots)))
	}
	# cold-start path: re-evaluate the original call with overrides
	if (is.null(object$call)) {
		cli::cli_abort(c(
			"No {.code call} stored on the fit object.",
			"i" = "Cold-start {.fn update} needs the original call. Refit explicitly with the new arguments."))
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
