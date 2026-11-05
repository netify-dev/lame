#' AME model fitting for bipartite (rectangular) networks
#'
#' Internal function that handles AME model fitting for bipartite networks.
#' This function contains the core logic for rectangular adjacency matrices
#' with distinct row and column node sets.
#'
#' @inheritParams ame
#' @param R_row integer: dimension of row node multiplicative effects
#' @param R_col integer: dimension of column node multiplicative effects
#' @return An ame object with posterior samples and estimates
#' @keywords internal
#' @noRd
ame_bipartite <- function(
	Y,
	Xdyad = NULL,
	Xrow = NULL,
	Xcol = NULL,
	rvar = TRUE,
	cvar = TRUE,
	R_row = 0,
	R_col = 0,
	family = "normal",
	intercept = !(family == "ordinal"),
	odmax = NULL,
	prior = list(),
	g = NA,
	seed = 6886,
	nscan = 10000,
	burn = 500,
	odens = 25,
	verbose = TRUE,
	gof = TRUE,
	custom_gof = NULL,
	start_vals = NULL,
	periodic_save = FALSE,
	out_file = NULL,
	save_interval = 0.25,
	model.name = NULL,
	posterior_opts = NULL,
	use_sparse_matrices = FALSE
) {

	# parameter setup
	# seed locally and restore the global RNG stream on exit
	.had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
	.old_seed <- if (.had_seed) get(".Random.seed", envir = globalenv()) else NULL
	on.exit({
		if (.had_seed) {
			assign(".Random.seed", .old_seed, envir = globalenv())
		} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			rm(".Random.seed", envir = globalenv())
		}
	}, add = TRUE)
	set.seed(seed)
	nA <- nrow(Y)
	nB <- ncol(Y)
	if (is.null(rownames(Y))) {
		rownames(Y) <- paste0("rowNode", 1:nA)
	}
	if (is.null(colnames(Y))) {
		colnames(Y) <- paste0("colNode", 1:nB)
	}

	n_eff <- floor(nscan / odens)

	Sab0 <- prior$Sab0
	if (is.null(Sab0)) Sab0 <- diag(2)
	eta0 <- prior$eta0
	# degrees of freedom as in lame(): round(4 + 3 n / 100), n the row actors
	if (is.null(eta0)) eta0 <- round(4 + 3 * nA / 100)
	etaab <- prior$etaab
	if (is.null(etaab)) etaab <- eta0
	s20 <- prior$s20
	if (is.null(s20)) s20 <- 1
	s2u0 <- prior$s2u0
	if (is.null(s2u0)) s2u0 <- 1
	# the bilinear term U G V' needs both ranks positive
	if (xor(R_row > 0L, R_col > 0L)) {
		cli::cli_warn(c(
			"!" = "Only one of {.arg R_row} ({R_row}) / {.arg R_col} ({R_col}) is positive.",
			"i" = "The multiplicative term U G V' needs both ranks > 0; it is dropped from this fit.",
			"i" = "Set both {.arg R_row} and {.arg R_col} (or pass {.arg R}) to include latent factors."
		))
		R_row <- 0L
		R_col <- 0L
	}
	# bipartite U and V rows have a unit N(0, I) prior with the scale carried by G, so Suv0 is refused
	.lame_refuse_bipartite_suv0(prior)

	# covariate setup
	if (!is.null(Xdyad)) {
		Xdyad <- as.array(Xdyad)
		if (length(dim(Xdyad)) == 2) {
			Xdyad <- array(Xdyad, c(nA, nB, 1))
		}
		p_dyad <- dim(Xdyad)[3]
	} else {
		p_dyad <- 0
	}

	if (!is.null(Xrow)) {
		Xrow <- as.matrix(Xrow)
		p_row <- ncol(Xrow)
	} else {
		p_row <- 0
	}

	if (!is.null(Xcol)) {
		Xcol <- as.matrix(Xcol)
		p_col <- ncol(Xcol)
	} else {
		p_col <- 0
	}

	p <- p_dyad + p_row + p_col + as.numeric(intercept)

	# a cell with a missing covariate becomes an unobserved tie (data-augmented), as in
	# ame_unipartite() and lame(); the covariate is then filled with 0
	if ((p_dyad + p_row + p_col) > 0) {
		n_na_before <- sum(is.na(Y))
		Y <- .ame_propagate_cov_na(
			Y,
			if (p_row > 0) Xrow,
			if (p_col > 0) Xcol,
			if (p_dyad > 0) Xdyad
		)
		.lame_inform_cov_na(sum(is.na(Y)) - n_na_before)
		if (p_dyad > 0 && anyNA(Xdyad)) Xdyad[is.na(Xdyad)] <- 0
		if (p_row > 0 && anyNA(Xrow)) Xrow[is.na(Xrow)] <- 0
		if (p_col > 0 && anyNA(Xcol)) Xcol[is.na(Xcol)] <- 0
	}

	# covariate names for the beta columns (.lame_apply_suffix() does not double-suffix)
	beta_names <- character(0)
	if (intercept) beta_names <- c(beta_names, "intercept")
	if (p_dyad > 0) {
		dyad_nms <- dimnames(Xdyad)[[3]]
		if (is.null(dyad_nms)) dyad_nms <- paste0("Xdyad", 1:p_dyad)
		beta_names <- c(beta_names, .lame_apply_suffix(dyad_nms, "dyad"))
	}
	if (p_row > 0) {
		row_nms <- colnames(Xrow)
		if (is.null(row_nms)) row_nms <- paste0("Xrow", 1:p_row)
		beta_names <- c(beta_names, .lame_apply_suffix(row_nms, "row"))
	}
	if (p_col > 0) {
		col_nms <- colnames(Xcol)
		if (is.null(col_nms)) col_nms <- paste0("Xcol", 1:p_col)
		beta_names <- c(beta_names, .lame_apply_suffix(col_nms, "col"))
	}

	# MCMC initialization
	# slots of start_vals replace the cold-start defaults below; a slot that does not fit
	# the panel is ignored and reported
	sv_in <- if (is.list(start_vals)) .start_vals_plain(start_vals) else list()
	sv_skipped <- character(0)
	sv_fits <- function(nm, ok) {
		x <- sv_in[[nm]]
		if (is.null(x)) return(NULL)
		if (isTRUE(ok(x))) return(x)
		sv_skipped <<- c(sv_skipped, nm)
		NULL
	}
	sv_dims <- function(x, d) is.numeric(x) && length(dim(x)) == length(d) && all(dim(x) == d)
	sv_effect <- function(nm, n) {
		x <- sv_in[[nm]]
		# a per-period start (from a dynamic fit) enters as the actor means
		if (is.numeric(x) && length(dim(x)) == 2L && nrow(x) == n && ncol(x) >= 1L) {
			x <- rowMeans(x)
			sv_in[[nm]] <<- x
		}
		sv_fits(nm, function(x) is.numeric(x) && is.null(dim(x)) && length(x) == n)
	}
	Z <- sv_fits("Z", function(x) sv_dims(x, c(nA, nB)))
	if (is.null(Z)) {
		if (family == "normal") {
			Z <- Y
		} else if (family == "binary") {
			Z <- matrix(rnorm(nA * nB), nA, nB)
			# start Z on the side of zero its observed tie implies
			Z[!is.na(Y) & Y == 1] <- abs(Z[!is.na(Y) & Y == 1])
			Z[!is.na(Y) & Y == 0] <- -abs(Z[!is.na(Y) & Y == 0])
		} else if (family == "ordinal") {
			Z <- matrix(qnorm((Y + 0.5) / (max(Y, na.rm = TRUE) + 1)), nA, nB)
		} else if (family == "cbin") {
			Z <- matrix(rnorm(nA * nB), nA, nB)
			Z[!is.na(Y) & Y == 1] <- abs(Z[!is.na(Y) & Y == 1])
			Z[!is.na(Y) & Y == 0] <- -abs(Z[!is.na(Y) & Y == 0])
		} else if (family == "frn") {
			Z <- matrix(rnorm(nA * nB), nA, nB)
			for (i in 1:nA) {
				if (any(Y[i, ] > 0, na.rm = TRUE)) {
					ranked <- which(Y[i, ] > 0)
					Z[i, ranked] <- sort(abs(rnorm(length(ranked))), decreasing = TRUE)
				}
			}
		} else if (family == "poisson") {
			# zero cells start at the mean rate (see .ame_pois_init_z); the square-bipartite diagonal is a real cell
			Z <- .ame_pois_init_z(Y, unipartite = FALSE)
		} else {
			Z <- Y
		}
	}

	beta <- sv_fits("beta", function(x) is.numeric(x) && is.null(dim(x)) && length(x) == p)
	if (is.null(beta)) beta <- rep(0, p)
	a <- sv_effect("a", nA)
	if (is.null(a)) a <- rep(0, nA)
	b <- sv_effect("b", nB)
	if (is.null(b)) b <- rep(0, nB)
	Sab <- sv_fits("Sab", function(x) sv_dims(x, c(2L, 2L)) && all(is.finite(diag(x))) && all(diag(x) > 0))
	sab_start_given <- !is.null(Sab)
	if (is.null(Sab)) Sab <- Sab0
	s2 <- sv_fits("s2", function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0)
	# unit of the response (R/normal_scale.R): 1 where the link fixes the latent scale, the response
	# variance for the normal family; the starting s2, G and the prior variance of G are in this unit
	unit_var <- .lame_unit_var(Y, family)
	unit_sd <- sqrt(unit_var)
	if (is.null(s2)) s2 <- unit_var

	# a supplied slot is kept as is; a drawn one carries the actor names
	U <- if (R_row > 0) sv_fits("U", function(x) sv_dims(x, c(nA, R_row))) else NULL
	if (R_row > 0 && is.null(U)) {
		U <- matrix(rnorm(nA * R_row, 0, 0.1), nA, R_row)
		rownames(U) <- rownames(Y)
	}

	V <- if (R_col > 0) sv_fits("V", function(x) sv_dims(x, c(nB, R_col))) else NULL
	if (R_col > 0 && is.null(V)) {
		V <- matrix(rnorm(nB * R_col, 0, 0.1), nB, R_col)
		rownames(V) <- colnames(Y)
	}

	if (R_row > 0 && R_col > 0) {
		G <- sv_fits("G", function(x) sv_dims(x, c(R_row, R_col)) && all(is.finite(x)))
		if (is.null(G)) G <- diag(min(R_row, R_col), R_row, R_col) * unit_sd
	} else {
		G <- NULL
	}
	if (length(sv_skipped) > 0L) {
		cli::cli_warn(c(
			"{.arg start_vals} slot{?s} {.val {sv_skipped}} do{?es/} not match this bipartite panel and {?is/are} ignored.",
			"i" = "Expected {.code Z} {nA} x {nB}, {.code beta} of length {p}, {.code a} of length {nA}, {.code b} of length {nB}, {.code U} {nA} x {R_row}, {.code V} {nB} x {R_col}, {.code G} {R_row} x {R_col}, {.code Sab} 2 x 2 and {.code s2} a positive scalar."
		))
	}

	# scale the additive-effects prior to the data, as the unipartite sampler does
	if (is.null(prior$Sab0)) {
		vscale_ab <- switch(
			family,
			normal = {
				Ec <- Y - mean(Y, na.rm = TRUE)
				a_tmp <- rowMeans(Ec, na.rm = TRUE); a_tmp[!is.finite(a_tmp)] <- 0
				b_tmp <- colMeans(Ec, na.rm = TRUE); b_tmp[!is.finite(b_tmp)] <- 0
				mean(c(stats::var(a_tmp), stats::var(b_tmp)))
			},
			poisson = .ame_pois_effect_scale(Y),
			binary = {
				YB <- 1 * (Y > 0)
				ybar <- mean(YB, na.rm = TRUE)
				if (is.finite(ybar) && ybar > 0 && ybar < 1) {
					E <- (YB - ybar) / dnorm(qnorm(ybar))
					a_tmp <- rowMeans(E, na.rm = TRUE); a_tmp[!is.finite(a_tmp)] <- 0
					b_tmp <- colMeans(E, na.rm = TRUE); b_tmp[!is.finite(b_tmp)] <- 0
					mean(c(stats::var(a_tmp), stats::var(b_tmp)))
				} else 1
			},
			# cbin and frn read the scale off the rows below their cap (a missing odmax is the largest outdegree)
			{
				sc_ab <- .lame_probit_ab_scale(
					Y,
					family,
					if (is.null(odmax)) max(rowSums(Y > 0, na.rm = TRUE)) else odmax
				)
				sc_ab$vscale
			}
		)
		if (!is.finite(vscale_ab) || vscale_ab <= 0) vscale_ab <- unit_var
		Sab0 <- diag(2) * vscale_ab
		# cbin and frn take one scale per side (see .lame_probit_ab_scale)
		if (family %in% c("cbin", "frn")) Sab0 <- diag(c(sc_ab$vscale_row, sc_ab$vscale_col))
		# a normal response starts its effect covariance at the prior scale of its own unit
		if (identical(family, "normal") && !sab_start_given) Sab <- Sab0
	}

	if (family == "binary" && is.null(prior$eta0)) {
		ydist <- table(Y)
		ymode <- as.numeric(names(ydist)[ydist == max(ydist)])[1]
		YB <- 1 * (Y != ymode)
		ybar <- mean(YB, na.rm = TRUE)
		mu_bin <- qnorm(ybar)
		E <- (YB - ybar) / dnorm(qnorm(ybar))
		a_tmp <- rowMeans(E, na.rm = TRUE)
		b_tmp <- colMeans(E, na.rm = TRUE)
		a_tmp[is.na(a_tmp)] <- 0
		b_tmp[is.na(b_tmp)] <- 0
		PHAT <- pnorm(mu_bin + outer(a_tmp, b_tmp, "+"))
		vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
		eta0 <- round(4 * vdfmlt)
		if (is.null(prior$etaab)) etaab <- eta0
	}

	# the fit reports the prior it used
	prior$Sab0 <- Sab0
	prior$eta0 <- eta0
	prior$etaab <- etaab

	BETA <- matrix(NA_real_, n_eff, p)
	if (length(beta_names) == p) colnames(BETA) <- beta_names
	VC <- matrix(NA_real_, n_eff, 5)
	# same column order as the unipartite and lame samplers; accessors read ve from the last column
	colnames(VC) <- c("va", "cab", "vb", "rho", "ve")

	U_samples <- NULL
	V_samples <- NULL
	a_samples <- NULL
	b_samples <- NULL

	if (!is.null(posterior_opts)) {
		# normalise partially specified lists so a missing field is FALSE, not NULL
		posterior_opts$save_UV <- isTRUE(posterior_opts$save_UV)
		posterior_opts$save_ab <- isTRUE(posterior_opts$save_ab)
		if (is.null(posterior_opts$thin_UV)) posterior_opts$thin_UV <- 1
		if (is.null(posterior_opts$thin_ab)) posterior_opts$thin_ab <- 1
		if (posterior_opts$save_UV && R_row > 0 && R_col > 0) {
			n_UV_samples <- ceiling(n_eff / posterior_opts$thin_UV)
			U_samples <- array(NA_real_, dim = c(nA, R_row, n_UV_samples))
			V_samples <- array(NA_real_, dim = c(nB, R_col, n_UV_samples))
			# the multiplicative mean is U G V' with a full G, so G is saved per draw too
			G_samples <- array(NA_real_, dim = c(R_row, R_col, n_UV_samples))
		}
		if (posterior_opts$save_ab && (rvar || cvar)) {
			n_ab_samples <- ceiling(n_eff / posterior_opts$thin_ab)
			if (rvar) a_samples <- matrix(NA_real_, nA, n_ab_samples)
			if (cvar) b_samples <- matrix(NA_real_, nB, n_ab_samples)
		}
	}

	UV_sample_idx <- 0
	ab_sample_idx <- 0
	APS <- rep(0, nA)
	APS2 <- rep(0, nA)
	BPS2 <- rep(0, nB)
	BPS <- rep(0, nB)
	YPS <- array(0, dim = c(nA, nB))
	UVPS <- array(0, dim = c(nA, nB))

	# custom GOF errors are reported once after fitting
	.cg_tracker <- .new_custom_gof_err_tracker()
	if (gof) {
		n_base_stats <- 3
		custom_gof_names <- NULL
		n_custom_stats <- 0
		if (!is.null(custom_gof)) {
			if (is.function(custom_gof)) {
				tryCatch({
					test_output <- custom_gof(Y)
					n_custom_stats <- length(test_output)
					custom_gof_names <- names(test_output)
					if (is.null(custom_gof_names)) {
						custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
					}
				}, error = function(e) {
					cli::cli_warn(c(
						"!" = "{.arg custom_gof} failed on the initial test call; the function will be skipped for this fit.",
						"i" = "Error: {conditionMessage(e)}",
						"i" = "Check that {.arg custom_gof} accepts a single matrix {.arg Y} and returns a numeric scalar / vector."
					))
					custom_gof <<- NULL
				})
			} else if (is.list(custom_gof)) {
				n_custom_stats <- length(custom_gof)
				custom_gof_names <- names(custom_gof)
				if (is.null(custom_gof_names)) {
					custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
				}
			}
		}

		n_total_stats <- n_base_stats + n_custom_stats
		GOF <- matrix(NA, n_eff + 1, n_total_stats)

		base_names <- .lame_bipartite_gof_names(family)
		all_names <- c(base_names, custom_gof_names)
		colnames(GOF) <- all_names

		# first row holds the observed GOF statistics
		GOF[1, 1:n_base_stats] <- gof_stats_bipartite(Y, warn_square = FALSE, family = family)
		if (!is.null(custom_gof)) {
			if (is.function(custom_gof)) {
				tryCatch({
					custom_vals <- custom_gof(Y)
					GOF[1, (n_base_stats + 1):(n_base_stats + length(custom_vals))] <- custom_vals
				}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
			} else if (is.list(custom_gof)) {
				for (i in seq_along(custom_gof)) {
					if (is.function(custom_gof[[i]])) {
						tryCatch({
							GOF[1, n_base_stats + i] <- custom_gof[[i]](Y)
						}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
					}
				}
			}
		}
	} else {
		GOF <- NULL
		custom_gof <- NULL
	}

	if (p > 0) {
		X <- array(0, c(nA, nB, p))
		idx <- 1

		if (intercept) {
			X[, , idx] <- 1
			idx <- idx + 1
		}

		if (p_dyad > 0) {
			for (k in 1:p_dyad) {
				X[, , idx] <- Xdyad[, , k]
				idx <- idx + 1
			}
		}

		if (p_row > 0) {
			for (k in 1:p_row) {
				X[, , idx] <- Xrow[, k]
				idx <- idx + 1
			}
		}

		if (p_col > 0) {
			for (k in 1:p_col) {
				for (j in 1:nB) {
					X[, j, idx] <- Xcol[j, k]
				}
				idx <- idx + 1
			}
		}
		# slice names let predict() hold the intercept and nodal slices at their fitted values
		if (length(beta_names) == p) dimnames(X) <- list(rownames(Y), colnames(Y), beta_names)
		.lame_check_covariate_scale(X)
		.lame_check_design_rank(X, Y)
	} else {
		X <- NULL
	}

	if (verbose) {
		cli::cli_h2("Running MCMC for bipartite network")
		cli::cli_text("Dimensions: {.val {nA}} x {.val {nB}}")
		cli::cli_text("R_row = {.val {R_row}}, R_col = {.val {R_col}}")
		cli::cli_text("Settings: Burn-in = {.val {burn}}, MCMC = {.val {nscan}}, Thin = {.val {odens}}")

		if (gof) {
			n_gof_stats <- 3
			if (!is.null(custom_gof)) {
				if (is.function(custom_gof)) {
					n_gof_stats <- n_gof_stats + length(custom_gof(Y))
				} else if (is.list(custom_gof)) {
					n_gof_stats <- n_gof_stats + length(custom_gof)
				}
			}
		}

		start_time <- Sys.time()

		if (burn > 0) {
			cli::cli_progress_bar("Burn-in", total = burn, clear = TRUE)
		}
	}

	# rank-family bookkeeping for the row constraints
	if (family %in% c("cbin", "frn")) {
		odobs_bip <- apply(Y > 0, 1, sum, na.rm = TRUE)
		# a hidden cell holds a nomination when its current latent is positive; cells of an absent
		# actor never do (.lame_nom_room)
		nom_hid_bip <- .lame_nom_hidden_cells(Y, TRUE)
		.odobs_room_bip <- function(Zc) .lame_nom_room(odobs_bip, nom_hid_bip, Zc)
		if (!is.null(odmax) && !(length(odmax) %in% c(1L, nA))) {
			cli::cli_abort(c(
				"{.arg odmax} must be a single value or one value per row.",
				"i" = "Got length {length(odmax)} for {nA} rows."
			))
		}
		odmax_vec <- if (is.null(odmax)) {
			rep(max(odobs_bip, na.rm = TRUE), nA)
		} else if (length(odmax) == 1L) {
			rep(odmax, nA)
		} else odmax
		odmax_vec <- pmin(odmax_vec, nB)
		odmax <- odmax_vec
		# posterior-predictive template (.lame_nom_predictive_yo)
		YO_pred <- .lame_nom_predictive_yo(Y, TRUE)
		# ranked receiver positions by row
		if (family == "frn") {
			ymax_global <- max(Y, na.rm = TRUE)
			if (is.finite(ymax_global) && ymax_global >= 1) {
				YL_bip <- matrix(NA_integer_, nA, ymax_global)
				for (i in seq_len(nA)) {
					yi <- Y[i, ]
					if (any(yi > 0, na.rm = TRUE)) {
						YL_bip[i, ] <- match(seq_len(ymax_global), yi)
					}
				}
			} else {
				YL_bip <- matrix(NA_integer_, nA, 1L)
			}
		} else {
			YL_bip <- NULL
		}
	} else {
		odobs_bip <- NULL; odmax_vec <- NULL; YL_bip <- NULL
	}

	# latent-array update by family; ordinal, cbin and frn respect their row and nomination constraints
	update_Z_fn <- switch(
		family,
		# missing cells are imputed from their conditional so every sufficient statistic sees a full matrix
		normal = function(Z, EZ, s2, Y) {
			Zn <- Y
			miss <- is.na(Y)
			if (any(miss)) Zn[miss] <- EZ[miss] + rnorm(sum(miss), 0, sqrt(s2))
			Zn
		},
		binary = function(Z, EZ, s2, Y) {
			# the probit scale is 1; missing cells draw from N(EZ, 1)
			Z_new <- EZ + matrix(rnorm(nA * nB), nA, nB)
			obs <- !is.na(Y)
			idx0 <- obs & Y == 0
			idx1 <- obs & Y == 1
			# exact one-sided truncated-normal draws on the log.p scale (clipping at 0 would leave an atom
			# at 0 and attenuate the coefficients)
			if (any(idx1)) {
				u <- pmax(runif(sum(idx1)), .Machine$double.xmin)
				Z_new[idx1] <- EZ[idx1] + qnorm(
					log(u) +
						pnorm(-EZ[idx1], lower.tail = FALSE, log.p = TRUE),
					lower.tail = FALSE,
					log.p = TRUE
				)
			}
			if (any(idx0)) {
				u <- pmax(runif(sum(idx0)), .Machine$double.xmin)
				Z_new[idx0] <- EZ[idx0] + qnorm(
					log(u) +
						pnorm(-EZ[idx0], log.p = TRUE),
					log.p = TRUE
				)
			}
			Z_new
		},
		poisson = function(Z, EZ, s2, Y) {
			rZ_pois_bip_fc(Z, EZ, s2, Y)
		},
		ordinal = function(Z, EZ, s2, Y) {
			rZ_ord_bip_fc(Z, EZ, Y)
		},
		cbin = function(Z, EZ, s2, Y) {
			rZ_cbin_bip_fc(Z, EZ, Y, odmax_vec, .odobs_room_bip(Z))
		},
		frn = function(Z, EZ, s2, Y) {
			rZ_frn_bip_fc(Z, EZ, Y, YL_bip, odmax_vec, .odobs_room_bip(Z))
		},
		cli::cli_abort(paste("Unknown family:", family))
	)

	needs_s2 <- family %in% c("normal", "poisson")
	s2_prior_scale <- .lame_residual_prior_scale(Y, family, prior$s20)
	# coefficient prior N(0, g s20 I): s20 is the response variance for a normal response, 1 on the
	# latent scale of the other families, or the user's prior$s20
	s20_beta <- if (identical(family, "normal") && is.null(prior$s20)) s2_prior_scale else s20

	UV_cache <- if (!is.null(U) && !is.null(V) && !is.null(G)) {
		U %*% G %*% t(V)
	} else {
		matrix(0, nA, nB)
	}

	# EZ_cache from the current parameter values
	EZ_cache <- matrix(0, nA, nB)
	if (p > 0) for (k in 1:p) EZ_cache <- EZ_cache + beta[k] * X[, , k]
	if (rvar) EZ_cache <- EZ_cache + a
	if (cvar) for (j in 1:nB) EZ_cache[, j] <- EZ_cache[, j] + b[j]
	EZ_cache <- EZ_cache + UV_cache

	if (p > 0) {
		XtX <- matrix(0, p, p)
		for (k1 in 1:p) {
			for (k2 in k1:p) {
				XtX[k1, k2] <- sum(X[, , k1] * X[, , k2], na.rm = TRUE)
				if (k1 != k2) XtX[k2, k1] <- XtX[k1, k2]
			}
		}
	}

	iter_save <- 0

	# cbin / frn: each sender's level moves with its latent row once per sweep (.lame_nom_row_shift)
	nom_shift <- family %in% c("cbin", "frn") && isTRUE(rvar)
	if (nom_shift) nom_shift_masks <- .lame_nom_shift_masks(array(Y, c(nA, nB, 1L)), TRUE)
	# and each receiver's level with its latent column (.lame_nom_col_shift); prior$actor_move = FALSE
	# switches it off, prior$actor_move_every sets the sweeps between moves (default 1)
	nom_col_move <- family %in% c("cbin", "frn") && isTRUE(cvar) && min(nA, nB) >= 3L
	if (!is.null(prior[["actor_move"]])) nom_col_move <- nom_col_move && isTRUE(prior[["actor_move"]])
	nom_col_stats <- c(visits = 0, col_calls = 0, col_skips = 0)
	if (nom_col_move) nom_col_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))
	# the levels are not centred, so the intercept and the level means share a flat ridge
	# exact draw along it (.lame_nom_ridge_bip) for every family but binary, which draws it in its own update
	ridge_move <- !identical(family, "binary") && isTRUE(intercept) && p > 0 && (isTRUE(rvar) || isTRUE(cvar))
	if (!is.null(prior[["actor_move"]])) ridge_move <- ridge_move && isTRUE(prior[["actor_move"]])
	ridge_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))
	ridge_calls <- 0
	# cbin / frn level move (.lame_nom_level_shift): rows free above move with the intercept against the rows with room
	nom_lvl_move <- family %in% c("cbin", "frn") && ridge_move && isTRUE(rvar) && nA >= 3L
	nom_lvl_stats <- c(visits = 0, calls = 0, free_rows = 0)
	if (nom_lvl_move) nom_lvl_masks <- .lame_nom_shift_masks(array(Y, c(nA, nB, 1L)), TRUE)

	# count failed conditional proposals over all sweeps
	tryErrorChecks <- list(U = 0L, V = 0L, G = 0L)

	# collapsed per-actor update (R/binary_actor_collapsed_move_bip.R): an exact elliptical slice draw of each
	# actor's level and position with the latent array integrated out, for actors with few events
	# prior$actor_move = FALSE switches it off, prior$actor_move_every sets the sweeps between visits
	if (!is.null(prior[["actor_move"]]) &&
			(!is.logical(prior[["actor_move"]]) || length(prior[["actor_move"]]) != 1L ||
					is.na(prior[["actor_move"]]))) cli::cli_abort("prior$actor_move must be TRUE or FALSE.")
	if (!is.null(prior[["actor_move_every"]]) &&
			(!is.numeric(prior[["actor_move_every"]]) || length(prior[["actor_move_every"]]) != 1L ||
					!is.finite(prior[["actor_move_every"]]) || prior[["actor_move_every"]] < 1 ||
					prior[["actor_move_every"]] != floor(prior[["actor_move_every"]])))
		cli::cli_abort("prior$actor_move_every must be a positive whole number.")
	actor_move <- identical(family, "binary") &&
		(isTRUE(rvar) || isTRUE(cvar) || (R_row > 0 && R_col > 0)) &&
		nA >= 2L && nB >= 2L && isTRUE(s2 == 1)
	if (!is.null(prior[["actor_move"]])) actor_move <- actor_move && isTRUE(prior[["actor_move"]])
	actor_move_uv <- actor_move && R_row > 0 && R_col > 0 && !is.null(U) && !is.null(V) && !is.null(G)
	actor_move_stats <- c(visits = 0, level_calls = 0, level_evals = 0, uv_calls = 0, uv_evals = 0, stuck = 0)
	# exact draw along the scale ridge of the positions and G (R/binary_bipartite_scale_moves.R)
	# prior$uv_g_scale_move = FALSE switches it off
	uv_g_scale_move <- R_row > 0 && R_col > 0 && !identical(prior[["uv_g_scale_move"]], FALSE)
	# exact draw along the ridge between row / column covariates and the additive effects
	# (R/bipartite_cov_ridge.R); prior$cov_ridge_move = FALSE switches it off
	cov_ridge <- p > 0 && (isTRUE(rvar) || isTRUE(cvar)) && !identical(prior[["cov_ridge_move"]], FALSE)
	uv_g_scale_stats <- c(visits = 0, drawn = 0)
	if (actor_move) {
		am_cx <- .lame_actor_move_context_bip(array(Y, c(nA, nB, 1L)))
		am_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 10L))
		# with an intercept a level translates together with it, otherwise alone
		am_comp <- isTRUE(intercept) && p > 0
		# the uncentred levels and the intercept also share a flat ridge
		if (am_comp) actor_move_stats <- c(actor_move_stats, ridge_calls = 0)
	}
	# translations of the additive levels with their latent rows and columns for Poisson and ordinal
	# fits (R/latent_translation_move.R), plus a location draw for the ordinal family
	# prior$actor_move = FALSE switches them off, prior$actor_move_every sets the sweeps between visits
	ltm_on <- family %in% c("poisson", "ordinal") && nA >= 2L && nB >= 2L &&
		(identical(family, "poisson") || isTRUE(rvar) || isTRUE(cvar)) && !identical(prior[["actor_move"]], FALSE)
	ltm_stats <- c(visits = 0, level_calls = 0, level_evals = 0, stuck = 0, skipped = 0, location_calls = 0)
	if (ltm_on) {
		ltm_cx <- .lame_ltm_context(array(Y, c(nA, nB, 1L)), family, "bip")
		# Poisson visits every 2nd sweep by default (as in lame()), other families every 10th
		ltm_every <- max(
			1L,
			as.integer(prior[["actor_move_every"]] %||% (if (identical(family, "poisson")) 2L else 10L))
		)
		ltm_comp <- isTRUE(intercept) && p > 0
		ltm_loc <- identical(family, "ordinal")
		# Poisson: the coefficients move with the latent array too, and the residual variance with the latent residuals
		ltm_coef <- identical(family, "poisson") && p > 0L && length(dim(X)) == 3L && dim(X)[3L] == p &&
			all(is.finite(X))
		if (ltm_coef) ltm_Xk <- lapply(seq_len(p), function(k) array(X[, , k], c(nA, nB, 1L)))
		ltm_scale <- identical(family, "poisson")
		ltm_uv <- identical(family, "poisson") && R_row > 0 && R_col > 0 && !is.null(U) && !is.null(V) && !is.null(G)
		ltm_stats <- c(
			ltm_stats,
			coef_calls = 0,
			coef_evals = 0,
			coef_stuck = 0,
			scale_attempts = 0,
			scale_accepts = 0,
			rho_attempts = 0,
			rho_accepts = 0,
			uv_calls = 0,
			uv_evals = 0,
			uv_stuck = 0
		)
	}

	# exact scale draws of the whole state and of the rank-likelihood latent array (R/ordinal_scale_move.R)
	# prior$ordinal_scale_move = FALSE switches them off
	osm_on <- identical(family, "ordinal") && !identical(prior[["ordinal_scale_move"]], FALSE)
	osm_stats <- c(visits = 0, orbit_moved = 0)
	# exact draw of the scale of a binary probit predictor with Z integrated out (R/binary_scale_orbit.R)
	# prior$binary_scale_move = FALSE switches it off, prior$binary_scale_every sets the sweeps between visits
	bso_on <- identical(family, "binary") && !identical(prior[["binary_scale_move"]], FALSE) && isTRUE(s2 == 1)
	bso_stats <- c(visits = 0, evals = 0, moved = 0)
	if (bso_on) {
		bso_every <- max(1L, as.integer(prior[["binary_scale_every"]] %||% 5L))
		bso_cx <- .lame_bso_context(array(Y, c(nA, nB, 1L)), "bip")
	}
	# Poisson residual variance moves with the latent residuals along their conditional law
	# every sweep (R/poisson_period_move.R)
	ppm_scale_on <- identical(family, "poisson")
	ppm_stats <- c(visits = 0, scale_attempts = 0, scale_accepts = 0, scale_rho_attempts = 0, scale_rho_accepts = 0)
	# the context also serves the coefficient translations
	if (ppm_scale_on || (ltm_on && ltm_coef)) {
		ppm_cx <- .lame_ppm_context(array(Y, c(nA, nB, 1L)), "bip", ltm = if (ltm_on) ltm_cx else NULL)
	}
	if (ppm_scale_on && ltm_on) ltm_scale <- FALSE

	# MCMC loop
	.chain_acc <- NULL
	for (iter in 1:(nscan + burn)) {
		EZ_curr <- EZ_cache

		Z <- update_Z_fn(Z, EZ_curr, s2, Y)
		if (nom_shift) {
			nom_sh <- .lame_nom_row_shift(
				array(Z, c(nA, nB, 1L)),
				a,
				b,
				Sab,
				odmax_vec,
				nom_shift_masks,
				TRUE,
				rep(TRUE, nA)
			)
			Z[] <- nom_sh$Z
			a <- nom_sh$a
		}
		if (nom_col_move && (iter %% nom_col_every == 0L)) {
			nom_cs <- .lame_nom_col_shift(
				array(Z, c(nA, nB, 1L)),
				array(Y, c(nA, nB, 1L)),
				family,
				if (family == "frn") list(YL_bip) else NULL,
				odmax_vec,
				a,
				b,
				Sab,
				TRUE,
				rep(TRUE, nB)
			)
			Z[] <- nom_cs$Z
			b <- nom_cs$b
			nom_col_stats <- nom_col_stats + c(1, nom_cs$calls, nom_cs$skips)
		}

		if (p > 0) {
			resid <- Z
			if (rvar) resid <- resid - a
			if (cvar) {
				for (j in 1:nB) {
					resid[, j] <- resid[, j] - b[j]
				}
			}
			if (!is.null(U) && !is.null(V) && !is.null(G)) {
				resid <- resid - UV_cache
			}

			XtX <- matrix(0, p, p)
			Xty <- rep(0, p)

			for (k1 in 1:p) {
				for (k2 in 1:p) {
					XtX[k1, k2] <- sum(X[, , k1] * X[, , k2], na.rm = TRUE)
				}
				Xty[k1] <- sum(X[, , k1] * resid, na.rm = TRUE)
			}

			if (length(g) == 1L && is.na(g)) g <- nA * nB

			V_beta <- solve(XtX / s2 + diag(p) / (g * s20_beta))
			m_beta <- V_beta %*% (Xty / s2)

			beta <- m_beta + t(chol(V_beta)) %*% rnorm(p)
		}

		if (rvar) {
			resid <- Z
			if (p > 0) {
				for (k in 1:p) {
					resid <- resid - beta[k] * X[, , k]
				}
			}
			if (cvar) {
				for (j in 1:nB) {
					resid[, j] <- resid[, j] - b[j]
				}
			}
			if (!is.null(U) && !is.null(V) && !is.null(G)) {
				resid <- resid - UV_cache
			}

			for (i in 1:nA) {
				prec_a <- nB / s2 + 1 / Sab[1, 1]
				mean_a <- sum(resid[i, ], na.rm = TRUE) / s2 / prec_a
				a[i] <- rnorm(1, mean_a, sqrt(1 / prec_a))
			}
		}

		if (cvar) {
			resid <- Z
			if (p > 0) {
				for (k in 1:p) {
					resid <- resid - beta[k] * X[, , k]
				}
			}
			if (rvar) resid <- resid - a
			if (!is.null(U) && !is.null(V) && !is.null(G)) {
				resid <- resid - UV_cache
			}

			for (j in 1:nB) {
				prec_b <- nA / s2 + 1 / Sab[2, 2]
				mean_b <- sum(resid[, j], na.rm = TRUE) / s2 / prec_b
				b[j] <- rnorm(1, mean_b, sqrt(1 / prec_b))
			}
		}

		if (rvar || cvar) {
			# va and vb are drawn independently (no row-column covariance)
			if (rvar) {
				Sab[1, 1] <- 1 / rgamma(1, (etaab + nA) / 2, (etaab * Sab0[1, 1] + sum(a^2)) / 2)
			}
			if (cvar) {
				Sab[2, 2] <- 1 / rgamma(1, (etaab + nB) / 2, (etaab * Sab0[2, 2] + sum(b^2)) / 2)
			}
			Sab[1, 2] <- Sab[2, 1] <- 0
		}

		if (cov_ridge) {
			cr <- .lame_bip_cov_ridge(X, beta, a, b, g * s20_beta, Sab[1, 1], Sab[2, 2], isTRUE(rvar), isTRUE(cvar))
			beta[] <- cr$beta
			a <- cr$a
			b <- cr$b
			rm(cr)
		}

		if (R_row > 0 && R_col > 0) {
			resid <- Z
			if (p > 0) {
				for (k in 1:p) {
					resid <- resid - beta[k] * X[, , k]
				}
			}
			if (rvar) resid <- resid - a
			if (cvar) {
				for (j in 1:nB) {
					resid[, j] <- resid[, j] - b[j]
				}
			}

			# exact row draws given the other factor and G (R/bipartite_uv_gibbs.R)
			U[] <- .lame_bip_factor_draw(array(resid, c(nA, nB, 1L)), V %*% t(G), s2)
			V[] <- .lame_bip_factor_draw(array(t(resid), c(nB, nA, 1L)), U %*% G, s2)

			# conjugate draw of G, which carries the scale of the multiplicative term under the unit
			# N(0, I) prior on U and V; sample_G_bip_cpp builds the precision and mean
			if (R_row > 0 && R_col > 0) {
				G_new <- tryCatch(
					sample_G_bip_cpp(
						array(resid, c(nA, nB, 1)),
						array(U, c(nA, R_row, 1)),
						array(V, c(nB, R_col, 1)),
						lambdaG = 1 / .lame_weight_prior_var(Y, family),
						s2 = s2
					),
					error = function(e) {
						tryErrorChecks$G <<- tryErrorChecks$G + 1L
						G
					}
				)
				if (is.matrix(G_new) && identical(dim(G_new), as.integer(c(R_row, R_col))) &&
						all(is.finite(G_new))) {
					G <- G_new
				} else {
					tryErrorChecks$G <- tryErrorChecks$G + 1L
				}
			}
			if (uv_g_scale_move && is.matrix(G)) {
				rg <- .lame_bip_uv_g_scale_step(U, V, G, g_var = unit_var)
				U <- rg$U
				V <- rg$V
				G <- rg$G
				uv_g_scale_stats <- uv_g_scale_stats + c(1, rg$drawn)
				rm(rg)
			}
			UV_cache <- U %*% G %*% t(V)
		}

		if (needs_s2) {
			# residual variance on the Gaussian latent scale (normal and Poisson)
			resid <- Z
			if (p > 0) {
				for (k in 1:p) {
					resid <- resid - beta[k] * X[, , k]
				}
			}
			if (rvar) resid <- resid - a
			if (cvar) {
				for (j in 1:nB) {
					resid[, j] <- resid[, j] - b[j]
				}
			}
			if (!is.null(U) && !is.null(V) && !is.null(G)) {
				resid <- resid - UV_cache
			}

			shape <- (sum(is.finite(resid)) + 1) / 2
			rate <- (sum(resid^2, na.rm = TRUE) + s2_prior_scale) / 2
			s2 <- 1 / rgamma(1, shape, rate)
		}

		# recompute EZ_cache from the current parameter values
		EZ_cache <- matrix(0, nA, nB)
		if (p > 0) for (k in 1:p) EZ_cache <- EZ_cache + beta[k] * X[, , k]
		if (rvar) EZ_cache <- EZ_cache + a
		if (cvar) for (j in 1:nB) EZ_cache[, j] <- EZ_cache[, j] + b[j]
		EZ_cache <- EZ_cache + UV_cache

		if (nom_lvl_move && (iter %% ridge_every == 0L)) {
			nom_lv <- .lame_nom_level_shift(
				array(Z, c(nA, nB, 1L)),
				a,
				b,
				Sab,
				odmax_vec,
				nom_lvl_masks,
				TRUE,
				rep(TRUE, nA),
				as.numeric(beta),
				diag(p) / (g * s20_beta),
				center = FALSE
			)
			Z[] <- nom_lv$Z
			a <- nom_lv$a
			beta[] <- nom_lv$beta
			EZ_cache <- matrix(0, nA, nB)
			for (k in 1:p) EZ_cache <- EZ_cache + beta[k] * X[, , k]
			if (rvar) EZ_cache <- EZ_cache + a
			if (cvar) for (j in 1:nB) EZ_cache[, j] <- EZ_cache[, j] + b[j]
			EZ_cache <- EZ_cache + UV_cache
			nom_lvl_stats <- nom_lvl_stats + c(1, nom_lv$calls, nom_lv$free)
		}
		if (ridge_move && (iter %% ridge_every == 0L)) {
			rg <- .lame_nom_ridge_bip(
				as.numeric(beta),
				a,
				b,
				diag(p) / (g * s20_beta),
				Sab[1L, 1L],
				Sab[2L, 2L],
				rvar,
				cvar
			)
			if (rg$calls > 0) {
				beta[] <- rg$beta
				a <- rg$a
				b <- rg$b
				ridge_calls <- ridge_calls + 1
			}
		}

		# collapsed per-actor update; the next sweep redraws the latent array given the moved state
		if (actor_move && (iter %% am_every == 0L)) {
			am_out <- .lame_actor_move_visit_bip(
				EZ_cache,
				as.numeric(beta),
				a,
				b,
				NULL,
				NULL,
				if (actor_move_uv) U else NULL,
				if (actor_move_uv) V else NULL,
				if (actor_move_uv) G else NULL,
				Sab[1L, 1L],
				Sab[2L, 2L],
				diag(p) / (g * s20_beta),
				am_comp,
				am_cx,
				isTRUE(rvar),
				isTRUE(cvar),
				free_means = TRUE
			)
			if (p > 0) beta[] <- am_out$beta
			a <- am_out$la
			b <- am_out$lb
			if (actor_move_uv) {
				U <- am_out$U
				V <- am_out$V
				UV_cache <- U %*% G %*% t(V)
			}
			EZ_cache <- matrix(am_out$EZ, nA, nB)
			actor_move_stats[c("visits", names(am_out$stats))] <-
				actor_move_stats[c("visits", names(am_out$stats))] + c(1, am_out$stats)
			rm(am_out)
		}

		# translations of the levels with their latent rows and columns (Poisson, ordinal)
		if (ltm_on && (iter %% ltm_every == 0L)) {
			ltm_st <- .lame_ltm_visit(
				ltm_cx,
				array(Z, c(nA, nB, 1L)),
				array(EZ_cache, c(nA, nB, 1L)),
				a,
				b,
				NULL,
				NULL,
				as.numeric(beta),
				Sab,
				diag(p) / (g * s20_beta),
				centred = ltm_comp,
				lev = c(isTRUE(rvar), isTRUE(cvar))
			)
			ltm_stats[c("visits", names(ltm_st$stats))] <- ltm_stats[c("visits", names(ltm_st$stats))] + c(1, ltm_st$stats)
			if (ltm_loc) {
				ltm_st <- .lame_ltm_location(
					ltm_cx,
					ltm_st$Z,
					ltm_st$EZ,
					ltm_st$la,
					ltm_st$lb,
					NULL,
					NULL,
					ltm_st$beta,
					Sab,
					diag(p) / (g * s20_beta),
					centred = ltm_comp,
					lev = c(isTRUE(rvar), isTRUE(cvar))
				)
				ltm_stats["location_calls"] <- ltm_stats["location_calls"] + 1
			}
			ltm_Z <- ltm_st$Z
			ltm_EZ <- ltm_st$EZ
			ltm_beta <- ltm_st$beta
			if (ltm_coef) {
				ltm_cf <- .lame_ppm_static(ppm_cx, ltm_Z, ltm_EZ, ltm_beta, ltm_Xk, diag(p) / (g * s20_beta))
				ltm_Z <- ltm_cf$Z; ltm_EZ <- ltm_cf$EZ; ltm_beta <- ltm_cf$beta
				ltm_stats[names(ltm_cf$stats)] <- ltm_stats[names(ltm_cf$stats)] + ltm_cf$stats
				rm(ltm_cf)
			}
			if (ltm_uv) {
				ltm_uo <- .lame_ltm_uv_visit(ltm_cx, ltm_Z, ltm_EZ, U, V, G)
				ltm_Z <- ltm_uo$Z; ltm_EZ <- ltm_uo$EZ
				U <- ltm_uo$U; V <- ltm_uo$V
				UV_cache <- U %*% G %*% t(V)
				ltm_stats[names(ltm_uo$stats)] <- ltm_stats[names(ltm_uo$stats)] + ltm_uo$stats
				rm(ltm_uo)
			}
			if (ltm_scale) {
				ltm_sc <- .lame_ltm_scale(
					ltm_cx,
					ltm_Z,
					ltm_EZ,
					s2,
					0,
					c(shape0 = 1 / 2, rate0 = s2_prior_scale / 2),
					two = FALSE
				)
				ltm_Z <- ltm_sc$Z; s2 <- ltm_sc$s2
				ltm_stats[names(ltm_sc$stats)] <- ltm_stats[names(ltm_sc$stats)] + ltm_sc$stats
				rm(ltm_sc)
			}
			Z[] <- ltm_Z
			EZ_cache <- matrix(ltm_EZ, nA, nB)
			if (p > 0) beta[] <- ltm_beta
			a <- ltm_st$la
			b <- ltm_st$lb
			rm(ltm_st, ltm_Z, ltm_EZ, ltm_beta)
		}

		# Poisson residual variance with the latent residuals carried along their conditional law
		if (ppm_scale_on) {
			ppm_sc <- .lame_ppm_scale(
				ppm_cx,
				array(Z, c(nA, nB, 1L)),
				array(EZ_cache, c(nA, nB, 1L)),
				s2,
				0,
				c(shape0 = 1 / 2, rate0 = s2_prior_scale / 2),
				FALSE,
				.lame_rho_prior_log
			)
			Z[] <- ppm_sc$Z
			s2 <- ppm_sc$s2
			ppm_stats <- ppm_stats + c(1, ppm_sc$stats)
			# the record of the translations is the work of the scale move that ran
			if (ltm_on) {
				ltm_stats[c("scale_attempts", "scale_accepts")] <- ltm_stats[c("scale_attempts", "scale_accepts")] +
					ppm_sc$stats[c("scale_attempts", "scale_accepts")]
			}
			rm(ppm_sc)
		}

		# scale of the whole state and the shift and scale of the rank-likelihood latent array
		if (osm_on) {
			osm_ab <- .lame_osm_prior_ab(a, b, Sab, rvar, cvar, bip = TRUE)
			osm_o <- .lame_osm_visit(
				array(Z, c(nA, nB, 1L)),
				array(EZ_cache, c(nA, nB, 1L)),
				array(UV_cache, c(nA, nB, 1L)),
				"bip",
				0,
				s2,
				A_pri = if (is.null(osm_ab)) 0 else (if (p > 0) sum(beta^2) / (g * s20_beta) else 0) + osm_ab[1L],
				B_pri = 0,
				d_th = if (is.null(osm_ab)) 0 else p + osm_ab[2L],
				orbit = !is.null(osm_ab)
			)
			Z[] <- osm_o$Z
			osm_stats <- osm_stats + c(1, osm_o$c != 1)
			if (osm_o$c != 1) {
				if (p > 0) beta[] <- osm_o$c * beta
				a <- osm_o$c * a
				b <- osm_o$c * b
				EZ_cache <- osm_o$c * (EZ_cache - UV_cache) + UV_cache
			}
			rm(osm_o, osm_ab)
		}

		# scale of a binary predictor, with the latent array integrated out
		if (bso_on && (iter %% bso_every == 0L)) {
			bso_ab <- .lame_osm_prior_ab(a, b, Sab, rvar, cvar, bip = TRUE)
			if (!is.null(bso_ab)) {
				bso_o <- .lame_bso_visit(
					array(Z, c(nA, nB, 1L)),
					array(EZ_cache, c(nA, nB, 1L)),
					array(UV_cache, c(nA, nB, 1L)),
					bso_cx,
					0,
					A_pri = (if (p > 0) sum(beta^2) / (g * s20_beta) else 0) + bso_ab[1L],
					B_pri = 0,
					d_th = p + bso_ab[2L]
				)
				bso_stats <- bso_stats + c(1, bso_o$evals, bso_o$c != 1)
				if (bso_o$c != 1) {
					Z[] <- bso_o$Z
					if (p > 0) beta[] <- bso_o$c * beta
					a <- bso_o$c * a
					b <- bso_o$c * b
					EZ_cache <- bso_o$c * (EZ_cache - UV_cache) + UV_cache
				}
				rm(bso_o)
			}
			rm(bso_ab)
		}

		if (iter > burn && (iter - burn) %% odens == 0) {
			iter_save <- iter_save + 1
			.chain_acc <- .lame_chain_acc_update(
				.chain_acc,
				iter_save,
				n_eff,
				if (rvar) a else NULL,
				if (cvar) b else NULL,
				U,
				V,
				G,
				actors = list(rownames(Y), colnames(Y))
			)

			if (iter_save <= n_eff) {
				BETA[iter_save, ] <- beta
				VC[iter_save, ] <- c(Sab[1, 1], Sab[1, 2], Sab[2, 2], 0, s2) # no rho for bipartite

				if (rvar) { APS <- APS + a; APS2 <- APS2 + a^2 }
				if (cvar) { BPS <- BPS + b; BPS2 <- BPS2 + b^2 }
				if (!is.null(U) && !is.null(V) && !is.null(G)) {
					UVPS <- UVPS + U %*% G %*% t(V)
				}

				# save posterior samples if requested
				if (!is.null(posterior_opts)) {
					# U / V samples (thinned)
					if (posterior_opts$save_UV && !is.null(U) && !is.null(V)) {
						if (iter_save %% posterior_opts$thin_UV == 0) {
							UV_sample_idx <- UV_sample_idx + 1
							if (UV_sample_idx <= dim(U_samples)[3]) {
								U_samples[, , UV_sample_idx] <- U
								V_samples[, , UV_sample_idx] <- V
								if (!is.null(G)) G_samples[, , UV_sample_idx] <- G
							}
						}
					}

					# a / b samples (thinned)
					if (posterior_opts$save_ab) {
						if (iter_save %% posterior_opts$thin_ab == 0) {
							ab_sample_idx <- ab_sample_idx + 1
							if (rvar && !is.null(a_samples) && ab_sample_idx <= ncol(a_samples)) {
								a_samples[, ab_sample_idx] <- a
							}
							if (cvar && !is.null(b_samples) && ab_sample_idx <= ncol(b_samples)) {
								b_samples[, ab_sample_idx] <- b
							}
						}
					}
				}
			}

			EZ_curr <- matrix(0, nA, nB)
			if (p > 0) {
				for (k in 1:p) {
					EZ_curr <- EZ_curr + beta[k] * X[, , k]
				}
			}
			if (rvar) EZ_curr <- EZ_curr + a
			if (cvar) {
				for (j in 1:nB) {
					EZ_curr[, j] <- EZ_curr[, j] + b[j]
				}
			}
			if (!is.null(U) && !is.null(V) && !is.null(G)) {
				EZ_curr <- EZ_curr + U %*% G %*% t(V)
			}

			if (family == "normal") {
				Y_sim <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
			} else if (family == "binary") {
				Y_sim <- 1 * (EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB) > 0)
			} else if (family == "ordinal") {
				# categories from the observed frequencies, as in the unipartite sampler
				Y_sim <- simY_ord(EZ_curr, 0, Y, bipartite = TRUE)
			} else if (family == "cbin") {
				# nominations from a latent draw around the predictor; hidden cells may hold a nomination,
				# cells of an absent actor stay out (.lame_nom_predictive_yo)
				Y_sim <- 1 * (simY_frn(EZ_curr, 0, odmax_vec, YO = YO_pred, bipartite = TRUE) > 0)
			} else if (family == "frn") {
				Y_sim <- simY_frn(EZ_curr, 0, odmax_vec, YO = YO_pred, bipartite = TRUE)
			} else if (family == "poisson") {
				# include the latent overdispersion: y ~ Pois(exp(z)), z ~ N(EZ, s2)
				z_sim <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
				Y_sim <- simY_pois(z_sim, bipartite = TRUE)
			} else {
				Y_sim <- EZ_curr
			}

			if (iter_save <= n_eff) {
				YPS <- YPS + Y_sim
			}

			if (gof && iter_save <= n_eff) {
				# compare the replicate with the observed statistics on the observed cells only
				Y_sim[is.na(Y)] <- NA
				gof_vals <- gof_stats_bipartite(Y_sim, warn_square = FALSE, family = family)
				GOF[iter_save + 1, 1:3] <- gof_vals

				if (!is.null(custom_gof)) {
					if (is.function(custom_gof)) {
						tryCatch({
							custom_vals <- custom_gof(Y_sim)
							GOF[iter_save + 1, (4):(3 + length(custom_vals))] <- custom_vals
						}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
					} else if (is.list(custom_gof)) {
						for (i in seq_along(custom_gof)) {
							if (is.function(custom_gof[[i]])) {
								tryCatch({
									GOF[iter_save + 1, 3 + i] <- custom_gof[[i]](Y_sim)
								}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
							}
						}
					}
				}
			}
		}

		if (verbose) {
			if (iter <= burn) {
				cli::cli_progress_update(id = NULL)
				if (iter == burn) {
					cli::cli_progress_done()
					cli::cli_progress_bar("Sampling", total = nscan, clear = TRUE)
				}
			} else if (iter > burn) {
				cli::cli_progress_update(id = NULL)
			}
		}
	}

	# output assembly
	if (verbose) {
		cli::cli_progress_done()

		end_time <- Sys.time()
		elapsed <- difftime(end_time, start_time, units = "secs")
		cli::cli_alert_success("MCMC complete in {.val {round(elapsed, 1)}} seconds")
	}

	chain_summary <- .lame_chain_acc_finish(.chain_acc)
	APM <- APS / n_eff
	APM_sd <- sqrt(pmax(APS2 / n_eff - APM^2, 0))
	BPM_sd <- sqrt(pmax(BPS2 / n_eff - (BPS / n_eff)^2, 0))
	names(APM) <- rownames(Y)
	BPM <- BPS / n_eff
	names(BPM) <- colnames(Y)
	YPM <- YPS / n_eff
	rownames(YPM) <- rownames(Y)
	colnames(YPM) <- colnames(Y)
	UVPM <- UVPS / n_eff
	rownames(UVPM) <- rownames(Y)
	colnames(UVPM) <- colnames(Y)

	fit <- list(
		BETA = BETA,
		VC = VC,
		APM = APM,
		BPM = BPM,
		U = U,
		V = V,
		G = G,
		UVPM = UVPM,
		YPM = YPM,
		GOF = GOF,
		X = X,
		Y = Y,
		family = family,
		# per-row nomination cap of cbin / frn, for simulate()
		odmax = odmax,
		R = max(R_row, R_col),
		R_row = R_row,
		R_col = R_col,
		rvar = rvar,
		cvar = cvar,
		symmetric = FALSE,
		mode = "bipartite",
		nA = nA,
		nB = nB,
		# numerical failures, distinct from ordinary Metropolis rejections
		tryErrorChecks = tryErrorChecks,
		mh_counters = tryErrorChecks,
		sampler_health = .lame_sampler_health(tryErrorChecks),
		sampler_health_iterations = nscan + burn,
		# resolved g and prior, for prior_summary()
		g = g,
		prior = prior,
		# the chain's final state, to continue it (start_vals of ame())
		start_vals = list(
			Z = Z,
			beta = as.vector(beta),
			a = as.vector(a),
			b = as.vector(b),
			U = U,
			V = V,
			G = G,
			s2 = as.numeric(s2),
			Sab = Sab
		),
		# posterior sd of the additive effects, for simulate()
		APM_sd = APM_sd,
		BPM_sd = BPM_sd
	)

	# collapsed per-actor update record: visits and slice likelihood evaluations
	if (actor_move && actor_move_stats[["visits"]] > 0) fit$actor_move <- actor_move_stats
	if (ltm_on && ltm_stats[["visits"]] > 0) fit$actor_move <- ltm_stats
	# rank-likelihood scale record: visits and visits at which the whole-state scale moved
	if (osm_on && osm_stats[["visits"]] > 0) fit$ordinal_scale_move <- osm_stats
	if (bso_on && bso_stats[["visits"]] > 0) fit$binary_scale_move <- bso_stats
	if (ppm_scale_on && ppm_stats[["visits"]] > 0) fit$poisson_period_move <- ppm_stats
	if (nom_col_move && nom_col_stats[["visits"]] > 0) fit$actor_move <- nom_col_stats
	if (nom_lvl_move && nom_lvl_stats[["visits"]] > 0)
		fit$actor_move <- c(
			if (is.null(fit$actor_move)) c(visits = nom_lvl_stats[["visits"]]) else fit$actor_move,
			level_calls = nom_lvl_stats[["calls"]],
			level_rows = nom_lvl_stats[["free_rows"]]
		)
	if (ridge_move && ridge_calls > 0)
		fit$actor_move <- c(
			if (is.null(fit$actor_move)) c(visits = ridge_calls) else fit$actor_move,
			ridge_calls = ridge_calls
		)
	if (uv_g_scale_stats[["visits"]] > 0) fit$uv_g_scale_move <- uv_g_scale_stats
	if (!is.null(chain_summary)) fit$chain_summary <- chain_summary

	if (!is.null(posterior_opts)) {
		if (!is.null(U_samples)) {
			actual_UV_samples <- min(UV_sample_idx, dim(U_samples)[3])
			if (actual_UV_samples > 0) {
				fit$U_samples <- U_samples[, , 1:actual_UV_samples, drop = FALSE]
				fit$V_samples <- V_samples[, , 1:actual_UV_samples, drop = FALSE]
				if (exists("G_samples") && !is.null(G_samples))
					fit$G_samples <- G_samples[, , 1:actual_UV_samples, drop = FALSE]
			}
		}
		if (!is.null(a_samples)) {
			actual_ab_samples <- min(ab_sample_idx, ncol(a_samples))
			if (actual_ab_samples > 0 && rvar) {
				fit$a_samples <- a_samples[, 1:actual_ab_samples, drop = FALSE]
			}
		}
		if (!is.null(b_samples)) {
			actual_ab_samples <- min(ab_sample_idx, ncol(b_samples))
			if (actual_ab_samples > 0 && cvar) {
				fit$b_samples <- b_samples[, 1:actual_ab_samples, drop = FALSE]
			}
		}
	}

	# one post-MCMC warning for errors raised by a user custom_gof
	if (exists(".cg_tracker", inherits = FALSE))
		.maybe_warn_custom_gof(.cg_tracker, where = "GOF")

	.lame_report_sampler_health(fit$sampler_health, fit$sampler_health_iterations)
	if (!is.null(s2_prior_scale)) fit$s2_prior_scale <- s2_prior_scale
	class(fit) <- "ame"
	# Poisson fitted count at the posterior means of the predictor and s2 (see .poisson_fitted_mean)
	if (identical(family, "poisson")) fit$YPM <- .poisson_fitted_from_fit(fit)
	.lame_separation_check(fit)

	if (use_sparse_matrices) {
		fit <- compact_ame(fit, use_sparse_matrices = use_sparse_matrices)
	}

	return(fit)
}
