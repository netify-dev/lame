#' Simulate networks from a fitted AME model
#'
#' @description
#' Generates networks from aligned joint posterior draws. A fit that did not
#' save them can use the conditional approximation by asking for it with
#' \code{type = "approximate"}; see Details.
#'
#' @param object fitted model object of class "ame"
#' @param nsim number of networks to simulate (default: 100)
#' @param seed Optional random seed. An explicit seed is local to this call:
#'   the caller's RNG state is restored on return or error. With NULL, the
#'   caller's current RNG stream advances normally.
#' @param newdata optional list containing new covariate data:
#'   \describe{
#'     \item{Xdyad}{dyadic covariates (n x n x p array or nA x nB x p for bipartite)}
#'     \item{Xrow}{row/sender covariates (n x p matrix or nA x p for bipartite)}
#'     \item{Xcol}{column/receiver covariates (n x p matrix or nB x p for bipartite)}
#'   }
#'   If NULL, uses covariates from original model fit
#' @param burn_in number of initial MCMC samples to discard (default: 0, assumes
#'   burn-in already removed)
#' @param thin thinning interval for MCMC samples (default: 1, use every sample)
#' @param return_latent logical: return latent Z matrices in addition to Y? (default: FALSE)
#' @param type simulation route: \code{"posterior"} (default) replicates the
#'   network from aligned saved posterior draws; \code{"approximate"} is the
#'   conditional approximation. \code{"auto"} is accepted as an alias of
#'   \code{"posterior"}. See Details.
#' @param ... additional arguments (not currently used)
#'
#' @return A list with components:
#'   \describe{
#'     \item{Y}{list of nsim simulated networks in the same format as the original data}
#'     \item{method}{joint posterior or conditional approximation}
#'     \item{draw_index}{sampled coefficient-row indices for joint replication}
#'     \item{Z}{if return_latent=TRUE, list of nsim latent Z matrices}
#'     \item{family}{the family of the model (binary, normal, etc.)}
#'     \item{mode}{network mode (unipartite or bipartite)}
#'   }
#'
#' @details
#' With \code{type = "posterior"}, every replicate uses the same saved MCMC
#' iteration for coefficients, actor effects, positions, coupling and observation
#' variance. Saved blocks may be thinned differently when their iteration indices
#' are available; only shared iterations are sampled. Missing joint storage is an
#' error (class \code{lame_simulate_draws}) that names \code{type = "approximate"}.
#' Store effects and positions with \code{posterior_options(save_ab = TRUE,
#' save_UV = TRUE, thin_ab = 1, thin_UV = 1)} when fitting.
#'
#' \code{type = "approximate"} is the conditional approximation: it substitutes
#' posterior summaries for the latent trajectories and does not propagate their
#' joint posterior uncertainty. It runs only when requested. The former default
#' \code{type = "auto"} chose between the two routes; it is now an alias of
#' \code{"posterior"}. The returned \code{method} identifies the route.
#'
#' Ordinal fits with explicit cutpoints use the corresponding saved cutpoint
#' draw. Data-induced ordinal fits match the observed category frequencies, so
#' these replications assess dependence conditional on those frequencies.
#'
#' @examples
#' \donttest{
#' # Fit a model
#' data(YX_bin)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_bin$Y, Xdyad = YX_bin$X, burn = 50, nscan = 250, odens = 5,
#'            family = "binary", verbose = FALSE,
#'            posterior_opts = posterior_options(save_ab = TRUE, save_UV = TRUE), chains = 1)
#'
#' # Simulate 10 networks from posterior
#' sims <- simulate(fit, nsim = 10)
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @importFrom stats simulate
#' @method simulate ame
#' @export
simulate.ame <- function(
	object,
	nsim = 100,
	seed = NULL,
	newdata = NULL,
	burn_in = 0,
	thin = 1,
	return_latent = FALSE,
	type = c("posterior", "approximate", "auto"),
	...
) {
	.lame_require_draws(object, "simulate")

	if (!is.null(seed)) {
		had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
		old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
		on.exit({
			if (had_seed) {
				assign(".Random.seed", old_seed, envir = globalenv())
			} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
				rm(".Random.seed", envir = globalenv())
			}
		}, add = TRUE)
		set.seed(seed)
	}

	type <- match.arg(type)
	# "auto" is an alias for the posterior route
	if (identical(type, "auto")) type <- "posterior"
	if (identical(type, "posterior")) {
		joint <- .simulate_joint_state(object, burn_in, thin)
		designs <- if (is.null(newdata)) NULL else list(.joint_new_design(object, newdata))
		return(.simulate_joint(object, joint, nsim, 1L, return_latent, designs))
	}

	fit <- object
	family <- if (!is.null(fit$family)) fit$family else "normal"
	symmetric <- if (!is.null(fit$symmetric)) {
		if (is.na(fit$symmetric)) FALSE else fit$symmetric
	} else FALSE
	mode <- if (!is.null(fit$mode)) fit$mode else "unipartite"
	bip <- (mode == "bipartite")

	if (bip) {
		nA <- fit$nA
		nB <- fit$nB
		n_row <- nA
		n_col <- nB
	} else {
		n <- length(fit$APM)
		n_row <- n_col <- n
	}

	BETA <- fit$BETA
	VC <- fit$VC
	stored_index <- seq_len(nrow(BETA))

	# burn-in and thinning
	if (burn_in > 0 && burn_in < nrow(BETA)) {
		BETA <- BETA[-(1:burn_in), , drop = FALSE]
		VC <- VC[-(1:burn_in), , drop = FALSE]
		stored_index <- stored_index[-(1:burn_in)]
	}
	if (thin > 1) {
		keep <- seq(1, nrow(BETA), by = thin)
		BETA <- BETA[keep, , drop = FALSE]
		VC <- VC[keep, , drop = FALSE]
		stored_index <- stored_index[keep]
	}

	n_mcmc <- nrow(BETA)

	# covariates
	if (!is.null(newdata)) {
		# design array the fit used, so X aligns with BETA
		has_int <- "intercept" %in% colnames(BETA)
		if (bip) {
			Xdyad <- newdata$Xdyad
			if (!is.null(Xdyad) && length(dim(Xdyad)) == 2) {
				Xdyad <- array(Xdyad, dim = c(dim(Xdyad), 1))
			}
			# intercept slice first, then the dyadic covariates, each named after its coefficient
			p_dyad <- if (is.null(Xdyad)) 0L else dim(Xdyad)[3]
			X <- array(1, dim = c(n_row, n_col, has_int + p_dyad))
			if (p_dyad > 0) X[, , has_int + seq_len(p_dyad)] <- Xdyad
			if (dim(X)[3] != ncol(BETA)) {
				cli::cli_abort(paste0(
					"newdata covariates produce ",
					dim(X)[3],
					" design columns but the fit has ",
					ncol(BETA),
					" coefficients. Supply the dyadic covariates the model was fit with."
				))
			}
			dimnames(X) <- list(NULL, NULL, colnames(BETA))
		} else {
			X <- design_array(
				newdata$Xrow,
				newdata$Xcol,
				newdata$Xdyad,
				intercept = has_int,
				n = n_row,
				warn = FALSE
			)
			if (dim(X)[3] != ncol(BETA)) {
				cli::cli_abort(paste0(
					"newdata covariates produce ",
					dim(X)[3],
					" design columns but the fit has ",
					ncol(BETA),
					" coefficients. Supply the covariates the model was fit with."
				))
			}
		}
	} else {
		if (!is.null(fit[["X"]])) {
			X <- fit[["X"]]
		} else {
			if (ncol(BETA) > 0L) cli::cli_warn(paste0(
				"Original covariates not stored in model. Using zero covariates. ",
				"Consider providing covariates via newdata argument."
			))
			X <- array(0, dim = c(n_row, n_col, 1))
		}
	}

	Y_sims <- vector("list", nsim)
	if (return_latent) {
		Z_sims <- vector("list", nsim)
	}

	# simulate networks
	for (i in 1:nsim) {
		s <- sample(1:n_mcmc, 1)
		beta <- BETA[s, ]

		# variance components of this draw
		vc_names <- colnames(VC)
		if (!is.null(vc_names)) {
			s2 <- if ("ve" %in% vc_names) VC[s, "ve"] else VC[s, ncol(VC)]
			rho <- if ("rho" %in% vc_names) VC[s, "rho"] else 0
			a_var <- if ("va" %in% vc_names) VC[s, "va"] else VC[s, 1]
			b_var <- if ("vb" %in% vc_names) VC[s, "vb"] else a_var

			a_var <- max(a_var, 0.01)
			b_var <- max(b_var, 0.01)
			s2 <- max(s2, 0.01)
		} else {
			if (symmetric) {
				s2 <- VC[s, ncol(VC)]
				rho <- 0
				a_var <- VC[s, 1]
				b_var <- a_var
			} else {
				if (bip) {
					s2 <- VC[s, 3]
					rho <- 0
					a_var <- VC[s, 1]
					b_var <- VC[s, 2]
				} else {
					s2 <- VC[s, ncol(VC)]
					rho <- VC[s, ncol(VC) - 1]
					a_var <- VC[s, 1]
					b_var <- VC[s, 3]
				}
			}
			a_var <- max(a_var, 0.01)
			b_var <- max(b_var, 0.01)
			s2 <- max(s2, 0.01)
		}

		# additive effects from saved draws when available, else from the posterior means
		if (!is.null(fit$a_samples)) {
			a <- fit$a_samples[, sample.int(ncol(fit$a_samples), 1)]
		} else if (!is.null(fit$APM) && length(fit$APM) == n_row) {
			# add the posterior spread to the shrunken means, so a convex link does not
			# under-predict the mean
			a <- as.numeric(fit$APM)
			if (!is.null(fit$APM_sd) && length(fit$APM_sd) == length(a)) {
				a <- a + rnorm(length(a), 0, pmax(as.numeric(fit$APM_sd), 0))
			}
		} else if (bip) {
			a <- rnorm(nA, mean = 0, sd = sqrt(a_var))
		} else {
			a <- rnorm(n, mean = 0, sd = sqrt(a_var))
		}
		if (isTRUE(symmetric) && !bip) {
			b <- a
		} else if (!is.null(fit$b_samples)) {
			b <- fit$b_samples[, sample.int(ncol(fit$b_samples), 1)]
		} else if (!is.null(fit$BPM) && length(fit$BPM) == n_col) {
			b <- as.numeric(fit$BPM)
			if (!is.null(fit$BPM_sd) && length(fit$BPM_sd) == length(b)) {
				b <- b + rnorm(length(b), 0, pmax(as.numeric(fit$BPM_sd), 0))
			}
		} else if (bip) {
			b <- rnorm(nB, mean = 0, sd = sqrt(b_var))
		} else {
			b <- rnorm(n, mean = 0, sd = sqrt(b_var))
		}

		# linear predictor
		if (length(beta) > 0 && any(beta != 0)) {
			if (length(dim(X)) == 3 && dim(X)[3] > 0) {
				# symmetric fits store one merged nodal coefficient per covariate, so slices map to coefficients by name
				beta_slice <- .beta_on_design(beta, dimnames(X)[[3]])
				EZ <- array(0, dim = c(n_row, n_col))
				for (k in 1:min(dim(X)[3], length(beta_slice))) {
					EZ <- EZ + beta_slice[k] * X[, , k]
				}
			} else {
				EZ <- matrix(beta[1], n_row, n_col) # intercept only
			}
		} else {
			EZ <- matrix(0, n_row, n_col)
		}

		EZ <- EZ + outer(a, b, "+")

		# multiplicative effects from a saved draw, else the posterior-mean U and V
		U_s <- fit$U; V_s <- fit[["V"]]
		# symmetric fits store U and L but no V: the term U L U' uses V = U
		if (is.null(V_s) && !is.null(U_s) && isTRUE(symmetric) && !bip) V_s <- U_s
		if (!is.null(fit$U_samples)) {
			uidx <- sample.int(dim(fit$U_samples)[3], 1)
			U_s <- fit$U_samples[, , uidx]
			V_s <- if (!is.null(fit$V_samples)) fit$V_samples[, , uidx] else U_s
		}
		# symmetric fits store no V: the term is U L U'
		if (is.null(V_s) && isTRUE(symmetric) && !bip) V_s <- U_s
		if (!is.null(U_s) && !is.null(V_s)) {
			if (bip && !is.null(fit[["G"]])) {
				EZ <- EZ + U_s %*% fit[["G"]] %*% t(V_s)
			} else if (!bip && isTRUE(symmetric) && !is.null(fit$L)) {
				L_mat <- if (is.matrix(fit$L)) fit$L else
					diag(fit$L, nrow = length(fit$L))
				EZ <- EZ + U_s %*% L_mat %*% t(V_s)
			} else if (!bip) {
				EZ <- EZ + U_s %*% t(V_s)
			}
		}

		if (symmetric && !bip) {
			EZ <- (EZ + t(EZ)) / 2
		}

		# simulate Y
		if (family == "binary") {
			Y_sim <- simY_bin(EZ, rho, bipartite = bip)
		} else if (family == "normal") {
			Y_sim <- simY_nrm(EZ, rho, s2, bipartite = bip)
		} else if (family == "ordinal") {
			# simY_ord takes Y (or YPM) as the category template
			y_for_sim <- fit$Y %||% fit$YPM
			if (is.null(y_for_sim)) cli::cli_abort("simulate.ame: ordinal needs `fit$Y` to be present.")
			Y_sim <- if (identical(fit$ordinal_cutpoints, "explicit"))
				.simY_ord_explicit(
					EZ,
					rho,
					.ordinal_saved_cutpoints(fit, stored_index[s]),
					fit$ordinal_levels,
					symmetric
				)
			else simY_ord(EZ, rho, y_for_sim, bipartite = bip)
		} else if (family == "poisson") {
			# latent overdispersion layer: Y ~ Poisson(exp(Z)), Z ~ N(EZ, s2 * Sigma_rho)
			Y_sim <- simY_pois(simZ(EZ, rho, s2, bipartite = bip), bipartite = bip)
		} else if (family == "cbin") {
			# censored binary: an frn draw thresholded at zero, as in the fitter's GOF draws
			odmax <- if (!is.null(fit$odmax)) fit$odmax else rep(Inf, n_row)
			Y_sim <- 1 * (simY_frn(EZ, rho, odmax, bipartite = bip) > 0)
		} else if (family == "frn") {
			odmax <- if (!is.null(fit$odmax)) fit$odmax else rep(Inf, n_row)
			Y_sim <- simY_frn(EZ, rho, odmax, bipartite = bip)
		} else {
			cli::cli_abort(paste0("Family ", family, " not yet implemented for simulation"))
		}

		if (!bip) {
			diag(Y_sim) <- NA
		}
		# simY_* draws independent noise per cell, so mirror the upper triangle for symmetric fits
		if (isTRUE(symmetric) && !bip && nrow(Y_sim) == ncol(Y_sim)) {
			Y_sim[lower.tri(Y_sim)] <- t(Y_sim)[lower.tri(Y_sim)]
		}

		Y_sims[[i]] <- Y_sim
		if (return_latent) {
			Z_sims[[i]] <- simZ(EZ, rho, s2, bipartite = bip)
		}
	}

	result <- list(
		Y = Y_sims,
		family = family,
		method = "conditional approximation",
		mode = mode
	)

	if (return_latent) {
		result$Z <- Z_sims
	}

	class(result) <- "ame.sim"
	return(result)
}
