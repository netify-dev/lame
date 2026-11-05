#' Simulate longitudinal networks from a fitted LAME model
#'
#' Generates network replicates at fitted periods from aligned posterior
#' trajectories. A fit that did not save them, and the legacy trajectory
#' scenarios (\code{start_from}, extended \code{n_time}), use the conditional
#' approximation by asking for it with \code{type = "approximate"}; see Details.
#'
#' @param object fitted model object of class "lame"
#' @param nsim number of network trajectories to simulate (default: 100)
#' @param seed Optional random seed. An explicit seed is local to this call:
#'   the caller's RNG state is restored on return or error. With NULL, the
#'   caller's current RNG stream advances normally.
#' @param newdata optional list containing new covariate data:
#'   \describe{
#'     \item{Xdyad}{list of T dyadic covariate arrays (n x n x p or nA x nB x p)}
#'     \item{Xrow}{list of T row/sender covariate matrices (n x p or nA x p)}
#'     \item{Xcol}{list of T column/receiver covariate matrices (n x p or nB x p)}
#'   }
#'   If NULL, uses covariates from original model fit
#' @param n_time number of time periods to simulate. If NULL, uses same as original data
#' @param burn_in number of initial MCMC samples to discard (default: 0)
#' @param thin thinning interval for MCMC samples (default: 1, use every sample)
#' @param return_latent logical: return latent Z matrices in addition to Y? (default: FALSE)
#' @param start_from character: how to initialize the simulation
#'   \describe{
#'     \item{"posterior"}{use saved trajectories for joint replication; posterior summaries for the approximation (default)}
#'     \item{"random"}{random initialization}
#'     \item{"data"}{use first time point from original data}
#'   }
#' @param type simulation route: \code{"posterior"} (default) replicates the
#'   network from aligned saved posterior trajectories; \code{"approximate"} is
#'   the conditional approximation. \code{"auto"} is accepted as an alias of
#'   \code{"posterior"}. See Details.
#' @param ... additional arguments (not currently used)
#'
#' @return A list with components:
#'   \describe{
#'     \item{Y}{list of nsim simulated longitudinal network trajectories,
#'            each element is a list of T networks}
#'     \item{method}{joint posterior or conditional approximation}
#'     \item{draw_index}{sampled coefficient-row indices for joint replication}
#'     \item{Z}{if return_latent=TRUE, list of nsim latent Z trajectories}
#'     \item{family}{the family of the model (binary, normal, etc.)}
#'     \item{mode}{network mode (unipartite or bipartite)}
#'     \item{n_time}{number of time periods}
#'   }
#'
#' @details
#' With \code{type = "posterior"}, every replicate uses the same saved MCMC
#' iteration for coefficients, actor effects, positions, coupling and observation
#' variance and period-specific reciprocity. Saved blocks may be thinned differently when their iteration indices
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
#' Joint longitudinal replication uses saved trajectories at fitted periods;
#' it does not evolve the final state through the historical covariates. Use
#' \code{predict(fit, h = ...)} for future periods. The legacy \code{start_from}
#' and extended \code{n_time} options apply to the conditional approximation
#' only; with \code{type = "posterior"} they are an error (class
#' \code{lame_simulate_periods}). The approximation extends \code{n_time} past
#' the fitted periods only for static coefficients; a \code{dynamic_beta} fit
#' has coefficients for its fitted periods alone, so a longer \code{n_time} is
#' the same error. In the approximation, \code{newdata$Xdyad}
#' must be a list with one array per period; any other shape is an error (class
#' \code{lame_simulate_xlist}).
#' Ordinal fits with explicit cutpoints use the corresponding saved cutpoint
#' draw. Data-induced ordinal fits match the observed category frequencies, so
#' these replications assess dependence conditional on those frequencies.
#'
#' @examples
#' \donttest{
#' # Create simple longitudinal network data
#' set.seed(1)
#' n <- 10
#' nms <- paste0("n", 1:n)
#' Y_list <- list(
#'   matrix(rnorm(n * n), n, n, dimnames = list(nms, nms)),
#'   matrix(rnorm(n * n), n, n, dimnames = list(nms, nms))
#' )
#' diag(Y_list[[1]]) <- diag(Y_list[[2]]) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(Y_list, family = "normal",
#'             nscan = 50, burn = 10, odens = 1, verbose = FALSE, plot = FALSE,
#'             posterior_opts = posterior_options(save_ab = TRUE, save_UV = TRUE), chains = 1)
#'
#' # Simulate 10 network trajectories from posterior
#' sims <- simulate(fit, nsim = 10)
#' }
#'
#' @author Shahryar Minhas
#' @method simulate lame
#' @export
simulate.lame <- function(
	object,
	nsim = 100,
	seed = NULL,
	newdata = NULL,
	n_time = NULL,
	burn_in = 0,
	thin = 1,
	return_latent = FALSE,
	start_from = "posterior",
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
	# "auto" is the posterior route
	if (identical(type, "auto")) type <- "posterior"
	if (identical(type, "posterior")) {
		T_fit <- object$n_time %||% length(object$Xlist)
		nt <- n_time %||% T_fit
		compatible <- identical(start_from, "posterior") && length(nt) == 1L && is.finite(nt) && nt >= 1L && nt == as.integer(nt) && nt <= T_fit
		if (!compatible) cli::cli_abort(
			c(
				"Posterior replication uses the fitted periods and {.code start_from = \"posterior\"}.",
				"i" = "Use {.code predict(h = ...)} for future trajectories, or {.code type = \"approximate\"} for the legacy {.arg start_from} and {.arg n_time} options."
			),
			class = "lame_simulate_periods"
		)
		joint <- .simulate_joint_state(object, burn_in, thin)
		designs <- if (is.null(newdata)) NULL else lapply(
			seq_len(nt),
			function(t) .joint_new_design(object, newdata, t)
		)
		return(.simulate_joint(object, joint, nsim, nt, return_latent, designs))
	}

	fit <- object
	# dynamic positions factored from the posterior mean products
	if (isTRUE(fit$dynamic_uv) && !identical(fit$uv_positions, "posterior_mean_product")) {
		fit <- .lame_dynamic_positions(fit)
	}
	family <- if (!is.null(fit$family)) fit$family else "normal"
	symmetric <- if (!is.null(fit$symmetric)) fit$symmetric else FALSE
	mode <- if (!is.null(fit$mode)) fit$mode else "unipartite"
	bip <- (mode == "bipartite")
	dynamic_uv <- if (!is.null(fit$dynamic_uv)) fit$dynamic_uv else FALSE
	dynamic_ab <- if (!is.null(fit$dynamic_ab)) fit$dynamic_ab else FALSE

	if (bip) {
		nA <- fit$nA
		nB <- fit$nB
		n_row <- nA
		n_col <- nB
	} else {
		if (is.matrix(fit$APM)) {
			n <- nrow(fit$APM)
		} else {
			n <- length(fit$APM)
		}
		n_row <- n_col <- n
	}

	# number of periods
	if (is.null(n_time)) {
		if (!is.null(fit$n_time)) {
			n_time <- fit$n_time
		} else if (!is.null(fit$Xlist)) {
			n_time <- length(fit$Xlist)
		} else {
			n_time <- 10 # default
			cli::cli_warn("Could not determine number of time periods from model. Using n_time=10")
		}
	}

	BETA <- fit$BETA
	VC <- fit$VC
	stored_index <- seq_len(dim(BETA)[1L])
	# 3-d BETA under dynamic_beta; the first dimension is the iteration
	beta_is_dyn <- length(dim(BETA)) == 3L

	# burn-in and thinning
	.dim1 <- function(x) dim(x)[1]
	if (burn_in > 0 && burn_in < .dim1(BETA)) {
		if (beta_is_dyn) BETA <- BETA[-(1:burn_in), , , drop = FALSE]
		else BETA <- BETA[-(1:burn_in), , drop = FALSE]
		VC <- VC[-(1:burn_in), , drop = FALSE]
		stored_index <- stored_index[-(1:burn_in)]
	}
	if (thin > 1) {
		keep <- seq(1, .dim1(BETA), by = thin)
		if (beta_is_dyn) BETA <- BETA[keep, , , drop = FALSE]
		else BETA <- BETA[keep, , drop = FALSE]
		VC <- VC[keep, , drop = FALSE]
		stored_index <- stored_index[keep]
	}

	n_mcmc <- .dim1(BETA)
	rho_period <- NULL
	if (isTRUE(fit$dynamic_rho) && !bip && !symmetric) {
		if (is.matrix(fit$RHO) && nrow(fit$RHO) == dim(fit$BETA)[1L]) {
			rho_period <- fit$RHO[stored_index, , drop = FALSE]
		} else if (is.numeric(fit$rho_path) && length(fit$rho_path)) {
			rho_period <- matrix(rep(fit$rho_path, each = n_mcmc), nrow = n_mcmc)
		} else cli::cli_abort("Period-specific reciprocity requires saved RHO draws or a rho_path summary.")
		if (ncol(rho_period) < n_time)
			cli::cli_abort("Simulation beyond the fitted reciprocity periods requires a future reciprocity model.")
		if (any(!is.finite(rho_period) | abs(rho_period) > 1))
			cli::cli_abort("Period-specific reciprocity must be finite and between -1 and 1.")
	}
	# per-period coefficients exist for the fitted periods only
	if (beta_is_dyn && n_time > dim(BETA)[3L]) {
		cli::cli_abort(
			c(
				"Simulation beyond the {dim(BETA)[3L]} fitted period{?s} of a {.arg dynamic_beta} fit needs future coefficient paths.",
				"i" = "Use {.code predict(h = ...)} for future periods, or set {.arg n_time} to at most {dim(BETA)[3L]}."
			),
			class = "lame_simulate_periods"
		)
	}

	# covariates
	if (!is.null(newdata)) {
		Xdyad_new <- newdata$Xdyad
		# one array per period: a bare array's length is its cell count, which the
		# cycling below would read as a number of periods
		if (!is.null(Xdyad_new) && (!is.list(Xdyad_new) || is.data.frame(Xdyad_new) || length(Xdyad_new) < 1L)) {
			cli::cli_abort(
				c(
					"{.code newdata$Xdyad} must be a list with one dyadic covariate array per period.",
					"i" = "It is {.cls {class(Xdyad_new)[1L]}} of length {length(Xdyad_new)}; the simulation has {n_time} period{?s}."
				),
				class = "lame_simulate_xlist"
			)
		}
		# periods supplied in newdata, recycled over the simulated periods
		nd_periods <- max(
			0L,
			vapply(
				newdata[intersect(c("Xdyad", "Xrow", "Xcol"), names(newdata))],
				function(x) if (is.list(x) && !is.data.frame(x)) length(x) else 0L,
				integer(1)
			)
		)
		# the design of each period, built as the fit's: intercept, nodal and dyadic layers in coefficient order
		Xlist <- lapply(seq_len(n_time), function(t) {
			tt <- if (nd_periods > 0L) ((t - 1L) %% nd_periods) + 1L else NULL
			.joint_new_design(fit, newdata, tt)
		})
	} else if (!is.null(fit$Xlist)) {
		Xlist <- fit$Xlist
		# extend or truncate to n_time
		orig_length <- length(Xlist)
		if (orig_length < n_time) {
			for (t in (orig_length + 1):n_time) {
				cycle_index <- ((t - 1) %% orig_length) + 1
				if (cycle_index > 0 && cycle_index <= orig_length) {
					Xlist[[t]] <- Xlist[[cycle_index]]
				} else {
					cli::cli_abort(paste0("Invalid cycle index: ", cycle_index, " for orig_length: ", orig_length))
				}
			}
		} else if (orig_length > n_time) {
			Xlist <- Xlist[1:n_time]
		}
	} else {
		cli::cli_warn("No covariates found. Using zero covariates.")
		Xlist <- lapply(1:n_time, function(t) {
			array(0, dim = c(n_row, n_col, 1))
		})
	}

	# the loop below reads Xlist[[t]] for every period
	if (length(Xlist) != n_time) {
		cli::cli_abort(
			c(
				"The covariate list has {length(Xlist)} element{?s} but the simulation has {n_time} period{?s}.",
				"i" = "Supply {.code newdata$Xdyad} with one array per period."
			),
			class = "lame_simulate_xlist"
		)
	}

	# temporal correlation parameters, drawn per iteration to carry their uncertainty
	rho_ab_samples <- if (!is.null(fit$rho_ab)) fit$rho_ab else rep(0, n_mcmc)
	rho_uv_samples <- if (!is.null(fit$rho_uv)) fit$rho_uv else rep(0, n_mcmc)
	# dynamic additive effects are a per-actor level plus an AR(1) deviation with
	# innovation sd sigma_ab; the level is the fit's posterior mean (or the time
	# average of the stored paths)
	sigma_ab_samples <- if (!is.null(fit$sigma_ab)) fit$sigma_ab else rep(NA_real_, n_mcmc)
	a_level_fit <- if (!is.null(fit$a_level)) as.numeric(fit$a_level) else if (is.matrix(fit$a_dynamic)) rowMeans(fit$a_dynamic, na.rm = TRUE) else NULL
	b_level_fit <- if (!is.null(fit$b_level)) as.numeric(fit$b_level) else if (is.matrix(fit$b_dynamic)) rowMeans(fit$b_dynamic, na.rm = TRUE) else NULL

	# the dynamic_uv AR(1) keeps the stationary spread of the reported positions
	u_marginal_sd <- v_marginal_sd <- 0
	if (dynamic_uv && !is.null(fit$U) && length(dim(fit$U)) == 3) {
		u_marginal_sd <- stats::sd(as.vector(fit$U), na.rm = TRUE)
		v_marginal_sd <- if (!is.null(fit[["V"]]) && length(dim(fit[["V"]])) == 3) {
			stats::sd(as.vector(fit[["V"]]), na.rm = TRUE)
		} else {
			u_marginal_sd
		}
	}

	Y_sims <- vector("list", nsim)
	if (return_latent) {
		Z_sims <- vector("list", nsim)
	}

	# simulate trajectories
	for (i in 1:nsim) {
		s <- sample(1:n_mcmc, 1)
		# beta_per_t is p x T (coefficient by period): the dynamic_beta draw, or the
		# static draw repeated across periods
		if (beta_is_dyn) {
			beta_per_t <- BETA[s, , , drop = FALSE]
			beta_per_t <- matrix(beta_per_t, nrow = dim(BETA)[2], ncol = dim(BETA)[3])
			beta <- as.numeric(apply(beta_per_t, 1, mean)) # period-averaged coefficients
		} else {
			beta <- BETA[s, ]
			beta_per_t <- matrix(beta, nrow = length(beta), ncol = n_time)
		}

		# per-iteration temporal parameters
		rho_ab <- rho_ab_samples[min(s, length(rho_ab_samples))]
		sigma_ab_s <- sigma_ab_samples[min(s, length(sigma_ab_samples))]
		rho_uv <- rho_uv_samples[min(s, length(rho_uv_samples))]

		# extract variance components
		vc_names <- colnames(VC)
		if (!is.null(vc_names)) {
			s2 <- if ("ve" %in% vc_names) VC[s, "ve"] else VC[s, ncol(VC)]
			rho <- if ("rho" %in% vc_names) VC[s, "rho"] else 0
			sigma_a <- sqrt(if ("va" %in% vc_names) VC[s, "va"] else VC[s, 1])
			sigma_b <- sqrt(if ("vb" %in% vc_names) VC[s, "vb"] else sigma_a^2)
		} else {
			if (symmetric) {
				s2 <- VC[s, ncol(VC)]
				rho <- 0
				sigma_a <- sqrt(VC[s, 1])
				sigma_b <- sigma_a
			} else {
				if (bip) {
					s2 <- VC[s, 3]
					rho <- 0
					sigma_a <- sqrt(VC[s, 1])
					sigma_b <- sqrt(VC[s, 2])
				} else {
					s2 <- VC[s, ncol(VC)]
					rho <- VC[s, ncol(VC) - 1]
					sigma_a <- sqrt(VC[s, 1])
					sigma_b <- sqrt(VC[s, 3])
				}
			}
		}

		Y_t_list <- vector("list", n_time)
		if (return_latent) {
			Z_t_list <- vector("list", n_time)
		}

		# initialize effects
		U_t <- V_t <- NULL
		if (start_from == "posterior") {
			if (bip) {
				a_t <- fit$APM
				b_t <- fit$BPM
			} else {
				a_t <- fit$APM
				b_t <- if (symmetric) fit$APM else fit$BPM
			}
			# static effects: a posterior draw restores the spread the shrunk posterior
			# means lack (dynamic effects get theirs below)
			if (!dynamic_ab && !is.null(fit$a_samples)) {
				a_t <- fit$a_samples[, sample.int(ncol(fit$a_samples), 1)]
				if (!(symmetric && !bip)) b_t <- if (!is.null(fit$b_samples)) fit$b_samples[, sample.int(ncol(fit$b_samples), 1)] else b_t
				if (symmetric && !bip) b_t <- a_t
			} else if (!dynamic_ab && !is.null(a_t)) {
				if (!is.null(fit$APM_sd) && length(fit$APM_sd) == length(a_t)) {
					a_t <- a_t + rnorm(length(a_t), 0, pmax(as.numeric(fit$APM_sd), 0))
				}
				if (symmetric && !bip) {
					b_t <- a_t
				} else if (!is.null(b_t) && !is.null(fit$BPM_sd) && length(fit$BPM_sd) == length(b_t)) {
					b_t <- b_t + rnorm(length(b_t), 0, pmax(as.numeric(fit$BPM_sd), 0))
				}
			}
			# under dynamic effects the posterior start is the actor level
			a_lev <- if (dynamic_ab && !is.null(a_level_fit)) a_level_fit else a_t
			b_lev <- if (dynamic_ab && !is.null(b_level_fit)) { if (symmetric) a_lev else b_level_fit } else b_t
			if (dynamic_ab) {
				# jitter the level by its posterior sd: a plug-in mean has less spread than
				# the effects, and a convex link would then under-predict the mean of the
				# simulated networks
				a_jit <- if (is.matrix(fit$a_dynamic_sd) && nrow(fit$a_dynamic_sd) == length(a_lev))
					rnorm(length(a_lev), 0, rowMeans(fit$a_dynamic_sd, na.rm = TRUE)) else rep(0, length(a_lev))
				b_jit <- if (is.matrix(fit$b_dynamic_sd) && nrow(fit$b_dynamic_sd) == length(b_lev))
					rnorm(length(b_lev), 0, rowMeans(fit$b_dynamic_sd, na.rm = TRUE)) else rep(0, length(b_lev))
				a_jit[!is.finite(a_jit)] <- 0
				b_jit[!is.finite(b_jit)] <- 0
				a_lev <- a_lev + a_jit
				b_lev <- if (!bip && symmetric) a_lev else b_lev + b_jit
				# the last fitted period's deviation from the level is the state
				# the trajectory continues from
				a_dev <- if (is.matrix(fit$a_dynamic) && !is.null(a_level_fit) && nrow(fit$a_dynamic) == length(a_lev))
					fit$a_dynamic[, ncol(fit$a_dynamic)] - a_level_fit else 0
				b_dev <- if (is.matrix(fit$b_dynamic) && !is.null(b_level_fit) && nrow(fit$b_dynamic) == length(b_lev))
					fit$b_dynamic[, ncol(fit$b_dynamic)] - b_level_fit else 0
				a_dev[!is.finite(a_dev)] <- 0
				b_dev[!is.finite(b_dev)] <- 0
				a_t <- a_lev + a_dev
				b_t <- if (!bip && symmetric) a_t else b_lev + b_dev
			}

			if (!is.null(fit$U)) {
				# dynamic U/V are n x r x T; take the last slice and keep the matrix shape
				# when R = 1
				uv_sd_mean <- 0
				if (length(dim(fit$U)) == 3) {
					T_fit <- dim(fit$U)[3]
					U_t <- matrix(
						fit$U[, , T_fit],
						nrow = dim(fit$U)[1],
						ncol = dim(fit$U)[2]
					)
					# absent actors start at the origin; the last period's posterior sd
					# spreads the start state
					U_t[!is.finite(U_t)] <- 0
					if (!is.null(fit$U_sd) && length(dim(fit$U_sd)) == 3 && all(dim(fit$U_sd) == dim(fit$U))) {
						u_sd <- matrix(fit$U_sd[, , T_fit], nrow(U_t), ncol(U_t))
						u_sd[!is.finite(u_sd)] <- 0
						U_t <- U_t + matrix(rnorm(length(U_t), 0, as.numeric(u_sd)), nrow(U_t), ncol(U_t))
						uv_sd_mean <- mean(u_sd)
					}
				} else {
					U_t <- fit$U
				}
				if (!is.null(fit[["V"]])) {
					if (length(dim(fit[["V"]])) == 3) {
						V_t <- matrix(
							fit[["V"]][, , dim(fit[["V"]])[3]],
							nrow = dim(fit[["V"]])[1],
							ncol = dim(fit[["V"]])[2]
						)
						V_t[!is.finite(V_t)] <- 0
						# receiver positions use their own posterior sd when stored, else the sender average
						if (!is.null(fit$V_sd) && length(dim(fit$V_sd)) == 3 && all(dim(fit$V_sd) == dim(fit[["V"]]))) {
							v_sd <- matrix(fit$V_sd[, , dim(fit[["V"]])[3]], nrow(V_t), ncol(V_t))
							v_sd[!is.finite(v_sd)] <- 0
							V_t <- V_t + matrix(rnorm(length(V_t), 0, as.numeric(v_sd)), nrow(V_t), ncol(V_t))
						} else if (uv_sd_mean > 0) {
							V_t <- V_t + matrix(rnorm(length(V_t), 0, uv_sd_mean), nrow(V_t), ncol(V_t))
						}
					} else {
						V_t <- fit[["V"]]
					}
				} else {
					V_t <- U_t # symmetric case
				}
			}
		} else if (start_from == "random") {
			if (bip) {
				a_t <- rnorm(nA, 0, sigma_a)
				b_t <- rnorm(nB, 0, sigma_b)
			} else {
				a_t <- rnorm(n, 0, sigma_a)
				b_t <- if (symmetric) a_t else rnorm(n, 0, sigma_b)
			}
			# a fresh actor draws its level from the level variance; the
			# deviation then starts from its stationary distribution
			a_lev <- a_t
			b_lev <- b_t
			if (dynamic_ab && is.finite(sigma_ab_s) && abs(rho_ab) < 1) {
				sd_stat <- sigma_ab_s / sqrt(1 - rho_ab^2)
				a_t <- a_lev + rnorm(length(a_lev), 0, sd_stat)
				b_t <- if (!bip && symmetric) a_t else b_lev + rnorm(length(b_lev), 0, sd_stat)
			}

			if (!is.null(fit$U)) {
				R <- ncol(fit$U)
				# scale random init to match the fitted latent factor scale
				u_scale <- sd(as.vector(fit$U), na.rm = TRUE)
				if (is.na(u_scale) || u_scale < 0.01) u_scale <- 0.1
				if (bip) {
					U_t <- matrix(rnorm(nA * R, 0, u_scale), nA, R)
					V_t <- matrix(rnorm(nB * R, 0, u_scale), nB, R)
				} else {
					U_t <- matrix(rnorm(n * R, 0, u_scale), n, R)
					V_t <- if (symmetric) U_t else matrix(rnorm(n * R, 0, u_scale), n, R)
				}
			}
		}

		for (t in 1:n_time) {
			X_t <- Xlist[[t]]
			if (!is.null(rho_period)) rho <- rho_period[s, t]

			# dynamic additive effects: AR(1) on the deviation from the level; without a
			# stored innovation sd, the effect variance sets it through the stationary form
			if (dynamic_ab && t > 1) {
				innov_a <- if (is.finite(sigma_ab_s)) sigma_ab_s else sigma_a * sqrt(1 - rho_ab^2)
				innov_b <- if (is.finite(sigma_ab_s)) sigma_ab_s else sigma_b * sqrt(1 - rho_ab^2)
				if (bip) {
					a_t <- a_lev + rho_ab * (a_t - a_lev) + rnorm(nA, 0, innov_a)
					b_t <- b_lev + rho_ab * (b_t - b_lev) + rnorm(nB, 0, innov_b)
				} else {
					a_t <- a_lev + rho_ab * (a_t - a_lev) + rnorm(n, 0, innov_a)
					if (!symmetric) {
						b_t <- b_lev + rho_ab * (b_t - b_lev) + rnorm(n, 0, innov_b)
					} else {
						b_t <- a_t
					}
				}
			}

			# dynamic multiplicative effects AR(1)
			if (dynamic_uv && t > 1 && !is.null(U_t)) {
				innov_scale <- sqrt(max(1 - rho_uv^2, 0))
				U_t <- rho_uv * U_t + matrix(
					rnorm(
						prod(dim(U_t)),
						0,
						u_marginal_sd * innov_scale
					),
					nrow(U_t),
					ncol(U_t)
				)
				if (!symmetric) {
					V_t <- rho_uv * V_t + matrix(
						rnorm(
							prod(dim(V_t)),
							0,
							v_marginal_sd * innov_scale
						),
						nrow(V_t),
						ncol(V_t)
					)
				} else {
					V_t <- U_t
				}
			}

			# linear predictor (per-period beta from beta_per_t)
			beta_t_here <- beta_per_t[, t]
			if (length(beta_t_here) > 0 && any(beta_t_here != 0)) {
				if (length(dim(X_t)) == 3 && dim(X_t)[3] > 0) {
					EZ <- array(0, dim = c(n_row, n_col))
					for (k in 1:min(dim(X_t)[3], length(beta_t_here))) {
						EZ <- EZ + beta_t_here[k] * X_t[, , k]
					}
				} else {
					EZ <- matrix(beta_t_here[1], n_row, n_col)
				}
			} else {
				EZ <- matrix(0, n_row, n_col)
			}

			EZ <- EZ + outer(a_t, b_t, "+")

			# multiplicative effects; a symmetric fit's term is U L U', so V mirrors U and
			# L has to be applied
			if (!is.null(U_t) && !is.null(V_t)) {
				if (bip && !is.null(fit[["G"]])) {
					EZ <- EZ + U_t %*% fit[["G"]] %*% t(V_t)
				} else if (!bip && symmetric && !is.null(fit$L)) {
					L_mat <- if (is.matrix(fit$L)) fit$L else
						diag(fit$L, nrow = length(fit$L))
					EZ <- EZ + U_t %*% L_mat %*% t(V_t)
				} else if (!bip) {
					EZ <- EZ + U_t %*% t(V_t)
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
				# simY_ord needs Y for the category lookup
				y_for_sim <- if (is.list(fit$Y)) fit$Y[[t]] else
				if (length(dim(fit$Y)) == 3L) fit$Y[, , t] else fit$Y
				if (is.null(y_for_sim))
					cli::cli_abort("simulate.lame: ordinal needs `fit$Y` for the category template.")
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
				# Poisson draws carry the latent overdispersion layer (Y ~ Pois(exp(Z)),
				# Z ~ N(EZ, s2)) and the fit's period_exposure when stored
				Z_od <- simZ(EZ, rho, s2, bipartite = bip)
				exposure_t <- if (!is.null(fit$period_exposure))
					fit$period_exposure[min(t, length(fit$period_exposure))] else 1
				if (isTRUE(exposure_t != 1)) {
					Y_sim <- simY_pois_exposure(Z_od, exposure_t, bipartite = bip)
				} else {
					Y_sim <- simY_pois(Z_od, bipartite = bip)
				}
			} else if (family == "cbin") {
				# censored binary: an frn draw thresholded at zero, as in the fitter's GOF
				# draws; an actor absent from the period cannot be nominated
				odmax <- if (!is.null(fit$odmax)) fit$odmax else rep(Inf, n_row)
				Y_sim <- 1 * (simY_frn(EZ, rho, odmax, YO = .sim_absent_yo(fit, t, n_row, bip), bipartite = bip) > 0)
			} else if (family == "frn") {
				odmax <- if (!is.null(fit$odmax)) fit$odmax else rep(Inf, n_row)
				Y_sim <- simY_frn(EZ, rho, odmax, YO = .sim_absent_yo(fit, t, n_row, bip), bipartite = bip)
			} else {
				cli::cli_abort(paste0("Family ", family, " not yet implemented for simulation"))
			}

			if (!bip) {
				diag(Y_sim) <- NA
			}
			# a symmetric fit's draws are symmetric: mirror the upper triangle, since the
			# simY_* draws add independent noise per cell
			if (isTRUE(symmetric) && !bip && nrow(Y_sim) == ncol(Y_sim)) {
				Y_sim[lower.tri(Y_sim)] <- t(Y_sim)[lower.tri(Y_sim)]
			}

			# a unipartite fit under changing composition reports each period
			# on its own actors; the draw is returned on the same grid
			Y_t_list[[t]] <- .sim_period_grid(Y_sim, fit, t, bip)
			if (return_latent) {
				Z_t_list[[t]] <- .sim_period_grid(simZ(EZ, rho, s2, bipartite = bip), fit, t, bip)
			}
		}

		Y_sims[[i]] <- Y_t_list
		if (return_latent) {
			Z_sims[[i]] <- Z_t_list
		}
	}

	result <- list(
		Y = Y_sims,
		family = family,
		method = "conditional approximation",
		mode = mode,
		n_time = n_time
	)

	if (return_latent) {
		result$Z <- Z_sims
	}

	class(result) <- "lame.sim"
	return(result)
}

# template for the nomination draw under changing composition: actors absent
# from period t are left out of contention, as in the sampler (see
# .lame_nom_predictive_yo); NULL when every actor is present
.sim_absent_yo <- function(fit, t, n, bip) {
	if (isTRUE(bip) || !is.list(fit$YPM) || length(fit$YPM) < t) return(NULL)
	yp <- fit$YPM[[t]]
	if (is.null(yp) || nrow(yp) == n || is.null(rownames(yp))) return(NULL)
	an <- names(fit$APM) %||% rownames(fit$APM) %||% rownames(yp)
	if (length(an) != n || !all(rownames(yp) %in% an)) return(NULL)
	out <- !(an %in% rownames(yp))
	YO <- matrix(0, n, n)
	YO[out, ] <- NA
	YO[, out] <- NA
	YO
}
# cut a simulated period matrix, drawn over every actor of the fit, to the
# actors the fit reports for that period (fit$YPM[[t]]) when the two differ
.sim_period_grid <- function(M, fit, t, bip) {
	if (isTRUE(bip) || !is.list(fit$YPM) || length(fit$YPM) < t) return(M)
	yp <- fit$YPM[[t]]
	if (is.null(yp) || identical(dim(yp), dim(M)) || is.null(rownames(yp)) || is.null(colnames(yp))) return(M)
	an <- names(fit$APM) %||% rownames(fit$APM) %||% rownames(yp)
	if (length(an) != nrow(M) || !all(rownames(yp) %in% an) || !all(colnames(yp) %in% an)) return(M)
	dimnames(M) <- list(an, an)
	M[rownames(yp), colnames(yp), drop = FALSE]
}
