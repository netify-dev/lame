#' AME model fitting for unipartite (square) networks
#'
#' Internal function that handles AME model fitting for unipartite networks.
#' This function contains the core logic extracted from the main ame() function
#' specifically for square adjacency matrices.
#'
#' @inheritParams ame
#' @return An ame object with posterior samples and estimates
#' @keywords internal
#' @noRd
ame_unipartite <- function(
	Y,
	Xdyad = NULL,
	Xrow = NULL,
	Xcol = NULL,
	rvar = TRUE,
	cvar = TRUE,
	dcor = !symmetric,
	nvar = TRUE,
	R = 0,
	family = "normal",
	intercept = !(family == "ordinal"),
	symmetric = FALSE,
	odmax = rep(max(apply(Y > 0, 1, sum, na.rm = TRUE)), nrow(Y)),
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
	use_sparse_matrices = FALSE,
	ordinal_cutpoints = c("data_induced", "explicit")
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
	diag(Y) <- NA
	if (is.element(family, c("binary", "cbin"))) { Y <- 1 * (Y > 0) }

	if (is.element(family, c("cbin", "frn"))) {
		odobs <- apply(Y > 0, 1, sum, na.rm = TRUE)
		# a hidden cell holds a nomination when its latent value is positive; rows hold their
		# non-nominations below zero only while short of odmax (.lame_nom_hidden_cells, .lame_nom_room)
		nom_hid_cells <- .lame_nom_hidden_cells(Y, FALSE)
		.odobs_room <- function(Zc) .lame_nom_room(odobs, nom_hid_cells, Zc)
		# reject an NA or non-positive odmax when nominations are observed
		if (!is.null(odmax)) {
			if (anyNA(odmax)) {
				cli::cli_abort(c(
					"{.arg odmax} contains {.val NA}.",
					"i" = "Supply a positive integer (scalar or per-row vector), or pass {.code odmax = NULL} to infer from {.arg Y}."
				))
			}
			if (any(odmax <= 0L) && any(Y > 0, na.rm = TRUE)) {
				cli::cli_abort(c(
					"{.arg odmax} <= 0 but {.arg Y} contains nominations.",
					"i" = "{.code odmax = 0} prohibits any nominations; inconsistent with the observed data."
				))
			}
		}
		# a NULL odmax defaults to the largest observed outdegree, as in lame() and the bipartite path
		if (is.null(odmax)) {
			odmax <- max(odobs, na.rm = TRUE)
			cli::cli_inform(c(
				"i" = "{.arg odmax} not supplied for family {.val {family}}; using the maximum observed outdegree ({odmax})."
			))
		}
		if (!is.null(odmax) && !(length(odmax) %in% c(1L, nrow(Y)))) {
			cli::cli_abort(c(
				"{.arg odmax} must be a single value or one value per row.",
				"i" = "Got length {length(odmax)} for {nrow(Y)} rows."
			))
		}
		if (length(odmax) == 1) { odmax <- rep(odmax, nrow(Y)) }
	}

	n <- nrow(Y)
	if (is.null(rownames(Y))) {
		rownames(Y) <- paste0("node", 1:n)
	}
	if (is.null(colnames(Y))) {
		colnames(Y) <- paste0("node", 1:n)
	}

	if (symmetric) { Xcol <- Xrow ; rvar <- cvar <- nvar }

	# `g` is the scalar Zellner g-prior scale that rbeta_ab_fc() expects
	if (length(g) > 1L) {
		cli::cli_abort(c(
			"{.arg g} must be a single scalar (length {length(g)} supplied).",
			"i" = "Per-coefficient g-priors are not currently supported in the unipartite path."
		))
	}
	if (length(g) == 1L && is.na(g)) {
		if (family == "normal") {
			g <- sum(!is.na(Y)) * var(c(Y), na.rm = TRUE)
		} else {
			g <- sum(!is.na(Y))
		}
	}

	if (!is.list(prior)) {
		if (!is.null(prior)) {
			cli::cli_warn(c(
				"{.arg prior} must be a named list; got {.cls {class(prior)[1]}}.",
				"i" = "Ignoring it and fitting with the default priors."
			))
		}
		prior <- list()
	}
	vscale <- 1
	vdfmlt <- 1
	if (!(family %in% c("binary", "ordinal", "cbin", "frn"))) {
		# scale the additive and multiplicative effects priors to the data, as the binary path does
		if (is.null(prior$Sab0)) {
			if (family == "poisson") {
				# moments on the log-rate scale (the scale of a, b and uv' under the log link)
				vscale_ab <- .ame_pois_effect_scale(Y)
			} else {
				Ec <- Y - mean(Y, na.rm = TRUE)
				a_tmp <- rowMeans(Ec, na.rm = TRUE); a_tmp[is.na(a_tmp)] <- 0
				b_tmp <- colMeans(Ec, na.rm = TRUE); b_tmp[is.na(b_tmp)] <- 0
				vscale_ab <- mean(c(stats::var(a_tmp), stats::var(b_tmp)), na.rm = TRUE)
				if (!is.finite(vscale_ab) || vscale_ab <= 0) {
					vscale_ab <- max(stats::var(c(Y), na.rm = TRUE), 1e-12)
				}
			}
		}
		if (is.null(prior$Sab0)) {
			prior$Sab0 <- diag(2) * vscale_ab
		}
		if (is.null(prior$eta0)) {
			prior$eta0 <- round(4 + 3 * n / 100)
		}
	}

	# design matrix
	pr <- length(Xrow) / n
	pc <- length(Xcol) / n
	pd <- length(Xdyad) / n^2

	# a dyad with a missing covariate is an unobserved tie (data-augmented), as in lame() and
	# the bipartite sampler; a symmetric dyad missing in one triangle is filled from the other
	if ((pr + pc + pd) > 0) {
		if (isTRUE(symmetric) && pd > 0) Xdyad <- .lame_fill_symmetric_dyadic(Xdyad)
		n_na_before <- sum(is.na(Y))
		Y <- .ame_propagate_cov_na(Y, Xrow, Xcol, Xdyad, symmetric = symmetric)
		.lame_inform_cov_na(sum(is.na(Y)) - n_na_before)
	}
	# multiplicative-effects prior shared with lame() (R/uv_prior.R), read off the observed cells
	if (R > 0) {
		uv_default <- .lame_uv_prior(Y, family, R, odmax, directed = !symmetric)
		if (is.null(prior$Suv0)) prior$Suv0 <- uv_default$Suv0
		if (is.null(prior$kappa0)) prior$kappa0 <- uv_default$kappa0
	}

	s2_prior_scale <- .lame_residual_prior_scale(Y, family, prior$s20)
	# unit of the response (R/normal_scale.R): 1 where the link fixes the scale, the response variance
	# for the normal family; the starting values and the eigen-weight prior variance are in this unit
	unit_var <- .lame_unit_var(Y, family)
	unit_sd <- sqrt(unit_var)
	X <- design_array(Xrow, Xcol, Xdyad, intercept, n, warn = FALSE)
	.lame_check_covariate_scale(X)
	.lame_check_design_rank(X, Y)
	.ame_check_normal_g(g, family, dim(X)[3L])
	normal_rho_certificate <- .ame_check_normal_reciprocity(Y, X, rvar, cvar, dcor, symmetric, R, family)

	Xnames <- dimnames(X)[[3]]
	if (is.null(Xnames) || length(Xnames) == 0) {
		if (dim(X)[3] > 0) {
			Xnames <- paste0("X", 1:dim(X)[3])
		} else {
			Xnames <- character(0)
		}
	}

	if (family == "ordinal" && identical(ordinal_cutpoints, "data_induced") &&
			isTRUE(any(apply(X, 3, function(x) {var(c(x))}) == 0))) {
		cli::cli_warn("An intercept is not estimable using this procedure")
	}

	# Z initialization
	if (family == "frn") {
		# the sampler works on ranks renumbered to consecutive integers within each row
		# the fit stores the ranks the user supplied
		Y_user <- Y
		ymx <- max(apply(1 * (Y > 0), 1, sum, na.rm = TRUE))
		YL <- NULL
		warn <- FALSE
		for (i in 1:nrow(Y)) {
			yi <- Y[i, ] ; rnkd <- which(!is.na(yi) & yi > 0)
			if (length(yi[rnkd]) > length(unique(yi[rnkd]))) {warn <- TRUE}
			yi[rnkd] <- rank(yi[rnkd], ties.method = "random")
			Y[i, ] <- yi
			YL <- rbind(YL, match(1:ymx, yi))
		}
		if (warn) {cli::cli_warn("Random reordering used to break ties in ranks")}
	}

	if (family == "normal") { Z <- Y }
	if (family == "ordinal") { Z <- matrix(zscores(Y), nrow(Y), ncol(Y)) }
	if (family == "binary") {
		Z <- matrix(zscores(Y), nrow(Y), ncol(Y))
		if (sum(Y == 0, na.rm = TRUE) > 0 && sum(Y == 1, na.rm = TRUE) > 0) {
			z01 <- .5 * (max(Z[Y == 0], na.rm = TRUE) + min(Z[Y == 1], na.rm = TRUE))
			Z <- Z - z01
		}
	}
	if (family == "poisson") {
		# start the latent log-rates near the observed rates, with zero cells at the mean rate (see .ame_pois_init_z)
		Z <- .ame_pois_init_z(Y)
	}

	if (is.element(family, c("cbin", "frn"))) {
		Z <- Y
		for (i in 1:nrow(Y)) {
			yi <- Y[i, ]
			zi <- zscores(yi)
			rnkd <- which(!is.na(yi) & yi > 0)
			if (length(rnkd) > 0 && min(zi[rnkd]) < 0) {
				zi[rnkd] <- zi[rnkd] - min(zi[rnkd]) + 1e-3
			}

			if (length(rnkd) < odmax[i]) {
				urnkd <- which(!is.na(yi) & yi == 0)
				if (length(urnkd) > 0 && max(zi[urnkd]) > 0) {
					zi[urnkd] <- zi[urnkd] - max(zi[urnkd]) - 1e-3
				}
			}

			Z[i, ] <- zi
		}
	}

	# MCMC initialization
	mu <- mean(Z, na.rm = TRUE)
	a <- if (rvar) rowMeans(Z, na.rm = TRUE) else rep(0, nrow(Z))
	b <- if (cvar) colMeans(Z, na.rm = TRUE) else rep(0, ncol(Z))
	a[is.na(a)] <- 0 ; b[is.na(b)] <- 0
	ZA <- mu + outer(a, b, "+")
	Z[is.na(Z)] <- ZA[is.na(Z)]

	# count in-loop sampler failures, reported as fit$tryErrorChecks
	tryErrorChecks <- list(beta = 0L, Sab = 0L, s2 = 0L, rho = 0L, UV = 0L)
	if (!is.null(start_vals)) {
		# the state is plain numeric: an ALS start carries actor labels (see .start_vals_plain)
		start_vals <- .start_vals_plain(start_vals)
		# exact `[[` lookups: `$b` would partial-match `beta` and `$rho` would match `rho_uv` / `rho_ab`
		if (!is.null(start_vals[["beta"]])) beta <- start_vals[["beta"]] else beta <- rep(0, dim(X)[3])
		if (!is.null(start_vals[["s2"]])) s2 <- start_vals[["s2"]] else s2 <- unit_var
		if (!is.null(start_vals[["a"]])) a <- start_vals[["a"]]
		if (!is.null(start_vals[["b"]])) b <- start_vals[["b"]]
		if (!is.null(start_vals[["rho"]])) rho <- start_vals[["rho"]] else rho <- 0
		# a missing Sab is seeded from the current a / b spread, as in a cold start
		if (!is.null(start_vals[["Sab"]])) {
			Sab <- start_vals[["Sab"]]
		} else {
			Sab <- matrix(c(max(0.1 * unit_var, var(a)), 0, 0, max(0.1 * unit_var, var(b))), 2, 2)
		}
		if (!is.null(start_vals[["U"]])) U <- start_vals[["U"]] else U <- matrix(0, nrow(Y), R)
		if (!is.null(start_vals[["V"]])) V <- start_vals[["V"]] else V <- matrix(0, nrow(Y), R)
		if (!is.null(start_vals[["Z"]])) Z <- start_vals[["Z"]]
	} else {
		beta <- rep(0, dim(X)[3])
		if (dim(X)[3] > 0 && !all(is.na(Z))) {
			if (is.null(attributes(X)$XX)) {
				X <- precomputeX(X)
			}
			tryCatch({
				beta <- solve(attributes(X)$XX + diag(dim(X)[3])) %*%
					crossprod(attributes(X)$mX, c(Z))
				beta <- as.vector(beta)
			}, error = function(e) {
				beta <- rep(0, dim(X)[3])
			})
		}

		s2 <- unit_var

		# ordinal, cbin and frn have no intercept: a fit without covariates takes the same
		# data-scaled additive-effects prior as one with covariates (and as lame())
		if (dim(X)[3] == 0 && family %in% c("ordinal", "cbin", "frn")) {
			sc <- .lame_probit_ab_scale(Y, family, odmax)
			if (is.null(prior$Sab0)) {
				prior$Sab0 <- if (family %in% c("cbin", "frn")) {
					diag(c(sc$vscale_row, sc$vscale_col))
				} else diag(2) * sc$vscale
			}
			if (is.null(prior$eta0)) prior$eta0 <- round(4 * sc$vdfmlt)
		}
		if (dim(X)[3] > 0) {
			E_tmp <- Z - Xbeta(X, beta)
			a <- rowMeans(E_tmp, na.rm = TRUE) * rvar
			b <- colMeans(E_tmp, na.rm = TRUE) * cvar
			a[is.na(a)] <- 0
			b[is.na(b)] <- 0

			if (family %in% c("binary", "ordinal", "cbin", "frn")) {
				if (family == "binary") {
					# the probit moments the factor prior reads too (R/uv_prior.R)
					sc <- .lame_binary_probit_scale(Y)
					vscale <- sc$vscale
					vdfmlt <- sc$vdfmlt
				} else {
					# ordinal, cbin and frn: the probit-moment scale, since the rank-normal start values attenuate
					# the actor effects; cbin and frn read it off the rows below their cap
					sc <- .lame_probit_ab_scale(Y, family, odmax)
					vscale <- sc$vscale
					vdfmlt <- sc$vdfmlt
				}

				if (is.null(prior$Sab0)) {
					# cbin and frn take one scale per side (see .lame_probit_ab_scale)
					prior$Sab0 <- if (family %in% c("cbin", "frn")) {
						diag(c(sc$vscale_row, sc$vscale_col))
					} else diag(2) * vscale
				}
				if (is.null(prior$eta0)) {
					prior$eta0 <- round(4 * vdfmlt)
				}
			}
		}

		E_init <- Z - Xbeta(X, beta) - mu - outer(a, b, "+")
		if (dcor && !all(is.na(E_init))) {
			rho <- cor(
				c(E_init[upper.tri(E_init)]),
				c(t(E_init)[upper.tri(E_init)]),
				use = "complete.obs"
			)
			if (is.na(rho)) rho <- 0
			# keep the start strictly inside (-1, 1): an exactly symmetric Y gives rho = 1 and a singular dyadic covariance
			rho <- min(max(rho, -0.9), 0.9)
		} else {
			rho <- 0
		}
		U <- V <- matrix(0, nrow(Y), R)
		if (R > 0) {
			E_init <- Z - mu - outer(a, b, "+")
			if (!all(is.na(E_init)) && any(E_init != 0, na.rm = TRUE)) {
				E_init[is.na(E_init)] <- 0
				tryCatch({
					svd_E <- svd(E_init)
					# a normal response takes the singular-vector sign that makes the largest coordinate of u positive,
					# so the start is the same for y and c y
					if (identical(family, "normal") && length(svd_E$d) >= R) {
						for (j in seq_len(R)) {
							if (svd_E$u[which.max(abs(svd_E$u[, j])), j] < 0) {
								svd_E$u[, j] <- -svd_E$u[, j]
								svd_E$v[, j] <- -svd_E$v[, j]
							}
						}
					}
					if (length(svd_E$d) >= R) {
						scaling_factor <- 0.5
						# a symmetric fit has unit-free coordinates and weights in the unit of the response (V = U L),
						# so the start splits the singular values that way
						split_u <- if (symmetric) 1 / sqrt(unit_sd) else 1
						split_v <- if (symmetric) sqrt(unit_sd) else 1
						U <- svd_E$u[, 1:R, drop = FALSE] %*% diag(
							sqrt(svd_E$d[1:R] * scaling_factor) * split_u,
							nrow = R
						)
						V <- svd_E$v[, 1:R, drop = FALSE] %*% diag(
							sqrt(svd_E$d[1:R] * scaling_factor) * split_v,
							nrow = R
						)
					}
				}, error = function(e) {
					U <- matrix(0, nrow(Y), R)
					V <- matrix(0, nrow(Y), R)
				})
			}
		}
	}

	if (!exists("Sab") || is.null(Sab)) {
		if (rvar && cvar) {
			Sab <- tryCatch({
				S <- cov(cbind(a, b)) * tcrossprod(c(rvar, cvar))
				if (any(is.na(S)) || det(S) <= 0) {
					diag(c(max(0.1 * unit_var, var(a, na.rm = TRUE)), max(0.1 * unit_var, var(b, na.rm = TRUE))))
				} else {
					S
				}
			}, error = function(e) {
				diag(c(0.1, 0.1) * unit_var)
			})
		} else if (rvar && !cvar) {
			Sab <- matrix(c(max(0.1 * unit_var, var(a)), 0, 0, 0.1 * unit_var), 2, 2)
		} else if (!rvar && cvar) {
			Sab <- matrix(c(0.1 * unit_var, 0, 0, max(0.1 * unit_var, var(b))), 2, 2)
		} else {
			Sab <- diag(c(0.1, 0.1) * unit_var)
		}
	}

	# storage allocation
	# normalise a partially specified posterior_opts so its flags are never NULL
	if (!is.null(posterior_opts)) {
		posterior_opts$save_UV <- isTRUE(posterior_opts$save_UV)
		posterior_opts$save_ab <- isTRUE(posterior_opts$save_ab)
		if (is.null(posterior_opts$thin_UV)) posterior_opts$thin_UV <- 1
		if (is.null(posterior_opts$thin_ab)) posterior_opts$thin_ab <- 1
	}
	n_samples <- floor((nscan + burn * (burn > 0)) / odens) - floor(burn / odens)
	BETA <- matrix(NA_real_, nrow = n_samples, ncol = dim(X)[3] - pr * symmetric)
	VC <- matrix(NA_real_, nrow = n_samples, ncol = 5 - 3 * symmetric)
	sample_idx <- 0
	U_samples <- NULL ; V_samples <- NULL
	a_samples <- NULL ; b_samples <- NULL

	if (!is.null(posterior_opts)) {
		if (posterior_opts$save_UV && R > 0) {
			n_UV_samples <- ceiling(n_samples / posterior_opts$thin_UV)
			U_samples <- array(NA_real_, dim = c(n, R, n_UV_samples))
			# symmetric fits store V = U L per draw too, so U L U' can be rebuilt draw by draw
			V_samples <- array(NA_real_, dim = c(n, R, n_UV_samples))
		}
		if (posterior_opts$save_ab && (rvar || cvar)) {
			n_ab_samples <- ceiling(n_samples / posterior_opts$thin_ab)
			if (rvar) a_samples <- matrix(NA_real_, n, n_ab_samples)
			# a symmetric fit has one effect per actor, so b_samples repeats a_samples, as in lame()
			if (cvar) b_samples <- matrix(NA_real_, n, n_ab_samples)
		}
	}

	UV_sample_idx <- 0
	ab_sample_idx <- 0

	UVPS <- U %*% t(V) * 0
	APS <- BPS <- rep(0, nrow(Y))
	# running sums of squares give the per-actor posterior sd of the effects
	APS2 <- BPS2 <- rep(0, nrow(Y))
	YPS <- matrix(0, nrow(Y), ncol(Y)) ; dimnames(YPS) <- dimnames(Y)

	if (gof) {
		n_base_stats <- 5
		custom_gof_names <- NULL
		n_custom_stats <- 0
		# custom_gof failures are counted and reported as one warning at the end of the run
		.cg_tracker <- .new_custom_gof_err_tracker()
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
		GOF <- matrix(NA_real_, nrow = n_samples + 1, ncol = n_total_stats)

		base_gof <- gof_stats(Y, mode = "unipartite")
		GOF[1, 1:n_base_stats] <- base_gof

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

		gof_idx <- 1
		rownames(GOF) <- c("obs", rep("", n_samples))
		base_names <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
		all_names <- c(base_names, custom_gof_names)
		colnames(GOF) <- all_names
	} else {
		GOF <- matrix(gof_stats(Y, mode = "unipartite"), 1, 5)
		rownames(GOF) <- "obs"
		colnames(GOF) <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
		custom_gof <- NULL
	}

	# pre-MCMC setup
	names(APS) <- names(BPS) <- rownames(U) <- rownames(V) <- rownames(Y)

	if (!symmetric) {
		colnames(VC) <- c("va", "cab", "vb", "rho", "ve")
		colnames(BETA) <- dimnames(X)[[3]]
	}

	if (symmetric) {
		colnames(VC) <- c("va", "ve")
		rb <- intercept + seq(1, pr, length = pr) ; cb <- intercept + pr + seq(1, pr, length = pr)
		bnames <- dimnames(X)[[3]]
		# build the exclude index explicitly: `-c(0, integer(0))` would silently empty the vector
		excl_idx <- c(if (intercept) 1L else integer(0), rb, cb)
		bni <- if (intercept) bnames[1L] else character(0)
		bnn <- if (length(rb) > 0L) gsub("row", bnames[rb], replacement = "node") else character(0)
		bnd <- if (length(excl_idx) > 0L) bnames[-excl_idx] else bnames
		colnames(BETA) <- c(bni, bnn, bnd)
	}

	Xr <- apply(X, c(1, 3), sum)
	Xc <- apply(X, c(2, 3), sum)
	mX <- apply(X, 3, c)
	mXt <- apply(aperm(X, c(2, 1, 3)), 3, c)
	# apply() drops the cell extent of a coefficient-free design; keep it as an n^2 x 0 matrix
	if (dim(X)[3L] == 0L) {
		mX <- mXt <- matrix(numeric(0), nrow = n * n, ncol = 0L)
	}
	XX <- t(mX) %*% mX
	XXt <- t(mX) %*% mXt

	X_precomp <- X
	attributes(X_precomp) <- c(
		attributes(X),
		list(Xr = Xr, Xc = Xc, mX = mX, mXt = mXt, XX = XX, XXt = XXt)
	)


	save_points <- (burn:(nscan + burn))[(burn:(nscan + burn)) %% odens == 0]
	save_points <- save_points[
		round(quantile(
			1:length(save_points),
			probs = seq(save_interval, 1, save_interval)
		))]

	beta_running_mean <- rep(0, dim(X)[3] - pr * symmetric)
	vc_running_mean <- rep(0, 5 - 3 * symmetric)

	# explicit-cutpoint state (ordinal_cutpoints = "explicit"): the explicit-alpha Z sampler
	# plus a Cowles (1996) MH update of alpha after each Z sweep
	ordinal_cutpoints <- match.arg(ordinal_cutpoints)
	use_explicit_cutpoints <- (family == "ordinal" &&
			identical(ordinal_cutpoints, "explicit"))
	.lame_explicit_ordinal_support(family, ordinal_cutpoints, FALSE, symmetric, dcor)
	if (use_explicit_cutpoints) {
		ord_lvls <- sort(unique(c(Y[is.finite(Y)])))
		K_ord <- length(ord_lvls)
		if (K_ord < 2L) {
			cli::cli_abort(c(
				"{.arg ordinal_cutpoints} = {.val explicit} requires at least 2 distinct observed categories.",
				"i" = "Found {K_ord}; check {.arg Y} or use a different family."
			))
		}
		Y_int <- matrix(
			match(Y, ord_lvls),
			nrow(Y),
			ncol(Y),
			dimnames = dimnames(Y)
		)
		alpha <- .init_alpha_from_data(Y)
		tau_prop <- 0.5
		# a continued chain keeps the cutpoints and proposal scale it stopped with
		if (is.list(start_vals) && is.numeric(start_vals[["alpha"]]) && length(start_vals[["alpha"]]) == length(alpha)) alpha <- start_vals[["alpha"]]
		if (is.list(start_vals) && is.numeric(start_vals[["tau_prop"]]) && length(start_vals[["tau_prop"]]) == 1L) tau_prop <- start_vals[["tau_prop"]]
		tau_target <- 0.30
		alpha_accept_n <- 0L
		alpha_accept_total <- 0L
		ALPHA <- if (K_ord >= 3L) {
			matrix(
				NA_real_,
				nrow = nscan / odens,
				ncol = K_ord - 2L,
				dimnames = list(NULL, paste0("alpha_", seq.int(2L, K_ord - 1L)))
			)
		} else NULL
	} else {
		alpha <- numeric(0)
		Y_int <- NULL
		tau_prop <- NA_real_
		K_ord <- 0L
		ALPHA <- NULL
	}

	# symmetric ordinal draws Z on the upper triangle and mirrors it (R/rZ_ord_sym_fc.R)
	# explicit cutpoints are handled inline in the MCMC loop
	update_Z_fn <- switch(
		family,
		normal = rZ_nrm_fc,
		# symmetric fits carry one latent per unordered pair (the *_sym_fc samplers)
		binary = if (symmetric) function(Z, EZ, rho, Y) rZ_bin_sym_fc(Z, EZ, Y, rho) else function(Z, EZ, rho, Y) {
			# compiled kernel (same full conditional and sweep order as rZ_bin_fc), with the R sampler as fallback
			z_new <- tryCatch(
				rZ_bin_fused_cpp(Z, EZ, rho, Y),
				error = function(e) NULL
			)
			if (is.null(z_new)) return(rZ_bin_fc(Z, EZ, rho, Y))
			dimnames(z_new) <- dimnames(Z)
			z_new
		},
		ordinal = if (symmetric) function(Z, EZ, rho, Y) rZ_ord_sym_fc(Z, EZ, Y, rho) else rZ_ord_fc,
		cbin = function(Z, EZ, rho, Y) rZ_cbin_fc(Z, EZ, rho, Y, odmax, .odobs_room(Z)),
		frn = function(Z, EZ, rho, Y) rZ_frn_fc(Z, EZ, rho, Y, YL, odmax, .odobs_room(Z)),
		# a symmetric fit carries one latent per unordered pair
		poisson = if (symmetric) function(Z, EZ, rho, s2, Y) rZ_pois_sym_fc(Z, EZ, s2, Y, rho = rho) else rZ_pois_fc,
		cli::cli_abort(paste("Unknown family:", family))
	)

	# cbin / frn: a hidden cell may hold a nomination in the replicate, as in the sampler
	# only the diagonal and absent actors stay out (.lame_nom_predictive_yo)
	YO_pred <- if (family %in% c("cbin", "frn")) .lame_nom_predictive_yo(Y, FALSE) else NULL
	simulate_Y_fn <- switch(
		family,
		normal = function(EZ, rho, s2) simY_nrm(EZ, rho, s2),
		binary = function(EZ, rho, s2) simY_bin(EZ, rho),
		ordinal = function(EZ, rho, s2) {
			if (use_explicit_cutpoints) .simY_ord_explicit(EZ, rho, alpha, ord_lvls, symmetric)
			else simY_ord(EZ, rho, Y)
		},
		cbin = function(EZ, rho, s2) 1 * (simY_frn(EZ, rho, odmax, YO = YO_pred) > 0),
		frn = function(EZ, rho, s2) simY_frn(EZ, rho, odmax, YO = YO_pred),
		# posterior-predictive draws include the latent layer: y ~ Pois(exp(z)), z ~ N(EZ, s2 Sigma_rho)
		poisson = function(EZ, rho, s2) simY_pois(simZ(EZ, rho, s2)),
		cli::cli_abort(paste("Unknown family:", family))
	)

	needs_s2 <- family %in% c("normal", "poisson")

	Xbeta_cache <- if (dim(X)[3] > 0) Xbeta(X, beta) else matrix(0, n, n)
	ab_outer_cache <- outer(a, b, "+")
	UV_cache <- U %*% t(V)
	EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache

	if (verbose) {
		cli::cli_h2("Running MCMC for unipartite network")
		cli::cli_text("Dimensions: {.val {n}} x {.val {n}} nodes")
		cli::cli_text("R = {.val {R}}, symmetric = {.val {symmetric}}")
		cli::cli_text("Settings: Burn-in = {.val {burn}}, MCMC = {.val {nscan}}, Thin = {.val {odens}}")

		if (gof) {
			n_gof_stats <- 5
			if (!is.null(custom_gof)) {
				if (is.function(custom_gof)) {
					tryCatch({
						n_gof_stats <- n_gof_stats + length(custom_gof(Y))
					}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
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

	# collapse the binary auxiliary variables before proposing reciprocity (prior in R/rho_prior.R)
	rho_collapsed_enabled <- identical(family, "binary") && isTRUE(dcor) && !symmetric
	rho_collapsed_attempts <- rho_collapsed_accepts <- 0L
	rho_collapsed_context <- if (rho_collapsed_enabled)
		.lame_binary_dynamic_collapsed_rho_context(array(Y, c(n, n, 1L))) else NULL
	if (rho_collapsed_enabled) tryErrorChecks$rho_collapsed <- 0L

	# symmetric binary and ordinal fits count each dyad once, with unit latent variance
	sym_gaussian_context <- if (family %in% c("binary", "ordinal") && isTRUE(symmetric)) {
		q <- if (dim(X)[3L]) .beta_prior_precision(
			attributes(X_precomp)$XX,
			attributes(X_precomp)$mX,
			g
		) else matrix(0, 0L, 0L)
		.lame_symmetric_gaussian_context(list(X), q, actor = isTRUE(nvar), center = FALSE)
	} else NULL

	# collapsed per-actor update (R/binary_actor_collapsed_move.R): exact elliptical slice draws of
	# each actor's level and static position with the latent array integrated out, for actors with few events
	# prior$actor_move = FALSE switches it off, prior$actor_move_every sets the sweeps between visits (default 10)
	am_lev <- if (symmetric) rep(isTRUE(nvar), 2L) else c(isTRUE(rvar), isTRUE(cvar))
	actor_move <- identical(family, "binary") && n >= 3L && (any(am_lev) || R > 0L)
	if (!is.null(prior[["actor_move"]])) actor_move <- actor_move && isTRUE(prior[["actor_move"]])
	actor_move_stats <- c(
		visits = 0, level_calls = 0, level_evals = 0,
		uv_calls = 0, uv_evals = 0, stuck = 0
	)
	if (actor_move) {
		am_cx <- if (symmetric) .lame_actor_move_context_sym(array(Y, c(n, n, 1L)))
		else .lame_actor_move_context(array(Y, c(n, n, 1L)))
		am_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 10L))
		# factor prior of rSuv_fc: inverse Wishart with kappa0 + n df and scale kappa0 Suv0 + W'W
		am_uv_prior <- if (!symmetric && R > 0L) {
			list(df = prior$kappa0, scale = prior$kappa0 * prior$Suv0)
		}
	}

	# exact scale draws of the whole state and of the rank-likelihood latent array (R/ordinal_scale_move.R)
	# prior$ordinal_scale_move = FALSE switches them off
	osm_on <- identical(family, "ordinal") && !identical(prior[["ordinal_scale_move"]], FALSE)
	if (osm_on) {
		osm_mode <- if (symmetric) "sym" else "dir"
		osm_iV0 <- if (dim(X)[3L]) .beta_prior_precision(
			attributes(X_precomp)$XX,
			attributes(X_precomp)$mX,
			g
		) else matrix(0, 0L, 0L)
		osm_stats <- c(visits = 0, orbit_moved = 0)
	}

	# exact draw of the scale of a binary probit predictor with Z integrated out (R/binary_scale_orbit.R)
	# prior$binary_scale_move = FALSE switches it off, prior$binary_scale_every sets the sweeps between visits
	bso_on <- identical(family, "binary") && !identical(prior[["binary_scale_move"]], FALSE)
	if (bso_on) {
		bso_every <- max(1L, as.integer(prior[["binary_scale_every"]] %||% 5L))
		bso_cx <- .lame_bso_context(array(Y, c(n, n, 1L)), if (symmetric) "sym" else "dir")
		bso_iV0 <- if (dim(X)[3L]) .beta_prior_precision(
			attributes(X_precomp)$XX,
			attributes(X_precomp)$mX,
			g
		) else matrix(0, 0L, 0L)
		bso_stats <- c(visits = 0, evals = 0, moved = 0)
	}

	# translations of the additive levels with their latent rows and columns for Poisson and ordinal
	# fits (R/latent_translation_move.R); Poisson also moves the coefficients, dyadic variance and
	# correlation, ordinal draws the location of the latent array
	# prior$actor_move = FALSE switches these off, prior$actor_move_every sets the sweeps between visits (default 10)
	ltm_on <- family %in% c("poisson", "ordinal") && n >= 3L && !identical(prior[["actor_move"]], FALSE) &&
		(identical(family, "poisson") || any(am_lev))
	uv_split_stats <- c(visits = 0, moved = 0)
	# positions of directed cbin, frn and ordinal fits move with the latent rows and columns they bound
	# (R/sparse_position_move.R); prior$sparse_position_move = FALSE omits it
	spm_on <- family %in% c("cbin", "frn", "ordinal") && !symmetric && R > 0L && n >= 3L &&
		!identical(prior[["sparse_position_move"]], FALSE)
	spm_stats <- c(visits = 0, spm_calls = 0, spm_evals = 0, spm_stuck = 0)
	spm_every <- max(1L, as.integer(prior[["sparse_position_every"]] %||% 1L))
	# ridge of a row or column covariate with the additive effects of its margin
	# (R/covariate_ridge_move.R); prior$cov_ridge_move = FALSE omits it
	cr_on <- !symmetric && isTRUE(rvar) && isTRUE(cvar) && dim(X)[3L] > 0L && !identical(
		prior[["cov_ridge_move"]],
		FALSE
	)
	cr_stats <- c(visits = 0, ridge_calls = 0, ridge_moved = 0)
	if (cr_on) cr_P0 <- .beta_prior_precision(XX, mX, g)
	if (spm_on) {
		spm_ctx <- .lame_spm_context(
			array(Y, c(n, n, 1L)),
			family,
			if (identical(family, "frn")) list(YL) else NULL,
			if (family %in% c("cbin", "frn")) odmax else NULL,
			Y_int = if (use_explicit_cutpoints) array(Y_int, c(n, n, 1L)) else NULL
		)
		spm_on <- isTRUE(spm_ctx$ok)
		spm_uv_prior <- list(df = prior$kappa0, scale = prior$kappa0 * prior$Suv0)
	}
	ltm_stats <- c(
		visits = 0, level_calls = 0, level_evals = 0, stuck = 0, skipped = 0, location_calls = 0,
		coef_calls = 0, coef_evals = 0, coef_stuck = 0, scale_attempts = 0, scale_accepts = 0, rho_attempts = 0, rho_accepts = 0
	)
	if (ltm_on) {
		ltm_cx <- .lame_ltm_context(
			array(Y, c(n, n, 1L)),
			family,
			if (symmetric) "sym" else "dir",
			Y_int = if (use_explicit_cutpoints) Y_int else NULL
		)
		ltm_on <- isTRUE(ltm_cx$ok)
	}
	if (ltm_on) {
		# Poisson visits every 2nd sweep by default (as in lame()), other families every 10th
		ltm_every <- max(
			1L,
			as.integer(prior[["actor_move_every"]] %||% (if (identical(family, "poisson")) 2L else 10L))
		)
		ltm_loc <- identical(family, "ordinal")
		ltm_coef <- identical(family, "poisson") && dim(X)[3L] > 0L && dim(X)[3L] == length(beta) && all(is.finite(X))
		if (ltm_coef) {
			ltm_Xk <- lapply(seq_len(dim(X)[3L]), function(k) array(X[, , k], c(n, n, 1L)))
			ltm_Xz <- if (symmetric) lapply(ltm_Xk, .lame_ppm_sym_array) else ltm_Xk
			ltm_iV0 <- .beta_prior_precision(XX, mX, g)
		}
		ltm_scale <- identical(family, "poisson")
		# Poisson: static positions move with the latent array too, under the rSuv_fc factor prior
		ltm_uv <- identical(family, "poisson") && R > 0L
		if (ltm_uv && !symmetric) {
			ltm_uv_prior <- list(df = prior$kappa0, scale = prior$kappa0 * prior$Suv0)
		}
		ltm_stats <- c(ltm_stats, uv_calls = 0, uv_evals = 0, uv_stuck = 0)
	}
	# Poisson dyadic variance and correlation move with the latent residuals along their
	# conditional law every sweep (R/poisson_period_move.R)
	ppm_scale_on <- identical(family, "poisson")
	ppm_stats <- c(visits = 0, scale_attempts = 0, scale_accepts = 0, scale_rho_attempts = 0, scale_rho_accepts = 0)
	# the context also serves the coefficient translations
	if (ppm_scale_on || (ltm_on && ltm_coef)) {
		ppm_cx <- .lame_ppm_context(
			array(Y, c(n, n, 1L)),
			if (symmetric) "sym" else "dir",
			ltm = if (ltm_on) ltm_cx else NULL
		)
	}
	if (ppm_scale_on && ltm_on) ltm_scale <- FALSE

	# cbin / frn receiver-level move (R/nomination_actor_move.R), the column counterpart of
	# raSab_cbin_fc / raSab_frn_fc; prior$actor_move = FALSE switches it off,
	# prior$actor_move_every sets the sweeps between moves (default 1)
	nom_col_move <- family %in% c("cbin", "frn") && !symmetric && n >= 3L && isTRUE(cvar)
	if (!is.null(prior[["actor_move"]])) nom_col_move <- nom_col_move && isTRUE(prior[["actor_move"]])
	nom_col_stats <- c(visits = 0, col_calls = 0, col_skips = 0)
	if (nom_col_move) nom_col_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))
	# level move (.lame_nom_level_shift): rows free above move with the intercept against the rows with room
	nom_lvl_move <- family %in% c("cbin", "frn") && !symmetric && n >= 3L && isTRUE(rvar) && intercept &&
		dim(X)[3L] > 0L && all(mX[, 1L] == 1)
	if (!is.null(prior[["actor_move"]])) nom_lvl_move <- nom_lvl_move && isTRUE(prior[["actor_move"]])
	nom_lvl_stats <- c(visits = 0, calls = 0, free_rows = 0)
	if (nom_lvl_move) {
		nom_lvl_every <- max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))
		nom_lvl_iV0 <- .beta_prior_precision(XX, mX, g)
		nom_lvl_masks <- .lame_nom_shift_masks(array(Y, c(n, n, 1L)), FALSE)
	}
	# collapsed reciprocity update (.lame_nom_rho_collapsed), with the non-nominations integrated out
	nom_rho_move <- family %in% c("cbin", "frn") && !symmetric && n >= 3L && isTRUE(dcor)
	if (!is.null(prior[["actor_move"]])) nom_rho_move <- nom_rho_move && isTRUE(prior[["actor_move"]])
	if (identical(prior[["actor_move_rho"]], FALSE)) nom_rho_move <- FALSE
	nom_rho_stats <- c(calls = 0, attempts = 0, accepts = 0)
	if (nom_rho_move) nom_rho_every <- 10L * max(1L, as.integer(prior[["actor_move_every"]] %||% 1L))

	.lame_explicit_ordinal_support(family, ordinal_cutpoints, FALSE, symmetric, dcor, rho)
	# MCMC loop
	.chain_acc <- NULL
	for (s in 1:(nscan + burn)) {

		EZ <- EZ_cache
		rho_eff <- rho
		Z <- if (use_explicit_cutpoints) {
			# explicit-cutpoint ordinal: rZ_ord_(sym_)explicit_fc ignores rho (unit-variance latent)
			if (symmetric) {
				rZ_ord_sym_explicit_fc(Z, EZ, Y_int, alpha)
			} else {
				rZ_ord_explicit_fc(Z, EZ, Y_int, alpha, rho = if (dcor) rho_eff else 0)
			}
		} else if (family %in% c("normal", "poisson")) {
			update_Z_fn(Z, EZ, rho_eff, s2, Y)
		} else {
			update_Z_fn(Z, EZ, rho_eff, Y)
		}

		# directed correlated models draw alpha given Z under the flat prior; independent and symmetric
		# models use the marginalised Cowles move, with Robbins-Monro tuning of tau_prop in burn-in only
		if (use_explicit_cutpoints) {
			alpha_step <- if (!symmetric && dcor) .sample_alpha_given_z(alpha, Y_int, Z)
			else sample_alpha_cowles(
				alpha = alpha,
				Y_int = Y_int,
				EZ = EZ,
				tau_prop = tau_prop,
				symmetric = symmetric
			)
			alpha <- alpha_step$alpha
			alpha_accept_total <- alpha_accept_total + isTRUE(alpha_step$accept)
			alpha_accept_n <- alpha_accept_n + 1L
			if ((symmetric || !dcor) && s <= burn && alpha_accept_n >= 25L) {
				acc_rate <- alpha_accept_total / alpha_accept_n
				log_tau <- log(tau_prop) + (acc_rate - tau_target) / sqrt(s)
				tau_prop <- min(max(exp(log_tau), 0.01), 2.0)
				alpha_accept_n <- 0L
				alpha_accept_total <- 0L
			}
			# refresh the latent array at the updated cutpoints
			Z <- if (symmetric) rZ_ord_sym_explicit_fc(Z, EZ, Y_int, alpha)
			else rZ_ord_explicit_fc(Z, EZ, Y_int, alpha, rho = if (dcor) rho_eff else 0)
		}

		if (needs_s2) {
			if (symmetric) {
				# symmetric path: e_ij = e_ji, so rs2_fc whitening would double-count the degrees of freedom
				# draw s2 from the upper-triangle residuals instead
				resid <- Z - EZ
				ut <- upper.tri(resid)
				e <- resid[ut]; e <- e[is.finite(e)]
				n_e <- length(e)
				if (n_e > 1) {
					ssr <- sum(e^2)
					nu0 <- 1
					s2 <- 1 / rgamma(
						1,
						shape = (n_e + nu0) / 2,
						rate = (ssr + nu0 * s2_prior_scale) / 2
					)
				}
			} else {
				s2 <- rs2_fc(Z, rho_eff, offset = EZ, s20 = s2_prior_scale)
			}
		}

		if (s == 1L && !is.null(normal_rho_certificate)) {
			initial_draw <- .ame_normal_reciprocity_initial_draw(
				Z,
				EZ,
				Sab,
				rho,
				X_precomp,
				s2,
				UV_cache,
				g,
				normal_rho_certificate,
				rvar,
				cvar
			)
			if (!is.null(initial_draw)) {
				beta <- initial_draw$beta
				a <- initial_draw$a * rvar
				b <- initial_draw$b * cvar
				Xbeta_cache <- if (dim(X)[3L] > 0L) Xbeta(X, beta) else matrix(0, n, n)
				ab_outer_cache <- outer(a, b, "+")
				EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
				EZ <- EZ_cache
			}
		}
		if (nom_rho_move && s %% nom_rho_every == 0L && isTRUE(s2 == 1)) {
			nom_r <- .lame_nom_rho_collapsed(
				array(Z, c(n, n, 1L)),
				array(EZ, c(n, n, 1L)),
				rho,
				array(Y, c(n, n, 1L)),
				family,
				if (family == "frn") list(YL) else NULL,
				odmax
			)
			nom_rho_stats <- nom_rho_stats + c(1, 4, nom_r$accepts)
			if (nom_r$accepts > 0L) {
				Z[] <- nom_r$Z
				rho <- nom_r$rho
			}
		}
		if (dcor) {
			# rho mixes slowest (random-walk MH), so its kernel runs several times per scan
			for (k_rho in 1:4) {
				rho <- rrho_mh(Z, rho, s2, offset = EZ)
			}
		}

		if (rho_collapsed_enabled && s %% 25L == 0L) {
			step <- tryCatch({
				z <- .ame_binary_collapsed_rho_mh(rho, EZ, rho_collapsed_context)
				z$Z <- if (z$rho != rho)
					.lame_binary_dynamic_collapsed_rho_z(EZ, Y, z$rho, rho_collapsed_context) else Z
				z
			}, error = function(e) {
				tryErrorChecks$rho_collapsed <<- tryErrorChecks$rho_collapsed + 1L
				NULL
			})
			if (!is.null(step)) {
				rho <- step$rho; Z <- step$Z
				rho_collapsed_attempts <- rho_collapsed_attempts + step$attempts
				rho_collapsed_accepts <- rho_collapsed_accepts + step$accepts
			}
		}

		if (symmetric) { rho <- min(.9999, 1 - 1 / sqrt(s)) }
		if (dcor || symmetric) {
			Z <- .lame_refresh_diagonal(Z, diag(EZ), rho, s2)
		}

		# a symmetric fit has one effect per actor: a rank-one covariance ties the row and column
		# effects so the shared effect keeps the reported prior variance (see .lame_sab_gls)
		tmp <- if (!is.null(sym_gaussian_context))
			.lame_symmetric_gaussian_draw(
				array(Z - UV_cache, c(n, n, 1L)),
				Sab[1L, 1L],
				sym_gaussian_context
			)
		else {
			Sab_active <- .lame_sab_gls(Sab, symmetric)
			if (!isTRUE(rvar)) Sab_active[1L, ] <- Sab_active[, 1L] <- 0
			if (!isTRUE(cvar)) Sab_active[2L, ] <- Sab_active[, 2L] <- 0
			rbeta_ab_fc(Z, Sab_active, rho, X_precomp, s2, offset = UV_cache, g = g)
		}
		beta_new <- tmp$beta
		a_new <- tmp$a * rvar
		b_new <- tmp$b * cvar

		if (symmetric) { a_new <- b_new <- (a_new + b_new) / 2 }

		if (dim(X)[3] > 0 && !isTRUE(all.equal(beta, beta_new))) {
			Xbeta_cache <- Xbeta(X, beta_new)
		}
		if (!isTRUE(all.equal(a, a_new)) || !isTRUE(all.equal(b, b_new))) {
			ab_outer_cache <- outer(a_new, b_new, "+")
		}

		beta <- beta_new
		a <- a_new
		b <- b_new

		EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache

		if (R > 0) {
			shrink <- (s > .5 * burn)
			offset <- Xbeta_cache + ab_outer_cache

			if (symmetric) {
				E <- Z - offset ; E <- .5 * (E + t(E))
				UV <- rUV_sym_fc(E, U, V, s2, shrink, l_var = unit_var)
			}
			if (!symmetric) {
				E <- Z - offset

				# rSuv_fc scales Suv0 by kappa0 internally
				Suv <- rSuv_fc(U, V, Suv0 = prior$Suv0, kappa0 = prior$kappa0)

				UV <- rUV_fc(E, U, V, Suv, rho, s2, offset = 0)
			}

			U <- UV$U ; V <- UV[["V"]]
			# scale split of the factors, with the covariance integrated out (R/uv_scale_split.R)
			if (!symmetric && !identical(prior[["uv_split_move"]], FALSE)) {
				uv_sp <- .lame_uv_split_draw(U, V, prior$kappa0, prior$kappa0 * prior$Suv0)
				U[] <- uv_sp$U; V[] <- uv_sp$V
				uv_split_stats <- uv_split_stats + c(1, uv_sp$c != 1)
				rm(uv_sp)
			}
			UV_cache <- U %*% t(V)
			EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
		}

		# normal / ordinal / poisson draw Sab by rSab_fc every scan; binary / cbin / frn use the joint
		# raSab_* updates when rvar, cvar and directed, and rSab_fc alone otherwise
		if (is.element(family, c("normal", "ordinal", "poisson"))) {
			Sab <- rSab_fc(
				a,
				b,
				Sab0 = prior$Sab0,
				eta0 = prior$eta0,
				tiny = 1e-8 * unit_var,
				rvar = rvar,
				cvar = cvar,
				symmetric = symmetric
			)
		}

		if (family == "binary") {
			if (rvar & cvar & !symmetric) {
				Sab <- rSab_fc(a, b, Sab0 = prior$Sab0, eta0 = prior$eta0, tiny = 1e-8 * unit_var)
				tmp <- raSab_bin_fc(Z, Y, a, b, Sab, Sab0 = prior$Sab0, eta0 = prior$eta0)
				Z <- tmp$Z;Sab <- tmp$Sab;a <- tmp$a
				ab_outer_cache <- outer(a, b, "+")
				EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			} else {
				Sab <- rSab_fc(
					a,
					b,
					Sab0 = prior$Sab0,
					eta0 = prior$eta0,
					tiny = 1e-8 * unit_var,
					rvar = rvar,
					cvar = cvar,
					symmetric = symmetric
				)
			}
		}

		if (family == "cbin") {
			if (rvar & cvar & !symmetric) {
				tmp <- raSab_cbin_fc(
					Z,
					Y,
					a,
					b,
					Sab,
					odmax,
					.odobs_room(Z),
					Sab0 = prior$Sab0,
					eta0 = prior$eta0
				)
				Z <- tmp$Z;Sab <- tmp$Sab;a <- tmp$a
				ab_outer_cache <- outer(a, b, "+")
				EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			} else {
				Sab <- rSab_fc(
					a,
					b,
					Sab0 = prior$Sab0,
					eta0 = prior$eta0,
					tiny = 1e-8 * unit_var,
					rvar = rvar,
					cvar = cvar,
					symmetric = symmetric
				)
			}
		}

		if (family == "frn") {
			if (rvar & cvar & !symmetric) {
				tmp <- raSab_frn_fc(
					Z,
					Y,
					YL,
					a,
					b,
					Sab,
					odmax,
					.odobs_room(Z),
					Sab0 = prior$Sab0,
					eta0 = prior$eta0
				)
				Z <- tmp$Z;Sab <- tmp$Sab;a <- tmp$a
				ab_outer_cache <- outer(a, b, "+")
				EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			} else {
				Sab <- rSab_fc(
					a,
					b,
					Sab0 = prior$Sab0,
					eta0 = prior$eta0,
					tiny = 1e-8 * unit_var,
					rvar = rvar,
					cvar = cvar,
					symmetric = symmetric
				)
			}
		}

		# cbin / frn: each receiver's level moves with its latent column (.lame_nom_col_shift)
		if (nom_col_move && (s %% nom_col_every == 0L)) {
			nom_cs <- .lame_nom_col_shift(
				array(Z, c(n, n, 1L)),
				array(Y, c(n, n, 1L)),
				family,
				if (family == "frn") list(YL) else NULL,
				odmax,
				a,
				b,
				Sab,
				FALSE,
				rep(TRUE, n)
			)
			Z[] <- nom_cs$Z
			b <- nom_cs$b
			ab_outer_cache <- outer(a, b, "+")
			EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			nom_col_stats <- nom_col_stats + c(1, nom_cs$calls, nom_cs$skips)
		}
		if (nom_lvl_move && (s %% nom_lvl_every == 0L)) {
			nom_lv <- .lame_nom_level_shift(
				array(Z, c(n, n, 1L)),
				a,
				b,
				Sab,
				odmax,
				nom_lvl_masks,
				FALSE,
				rep(TRUE, n),
				as.numeric(beta),
				nom_lvl_iV0,
				center = FALSE
			)
			Z[] <- nom_lv$Z
			a <- nom_lv$a
			beta[] <- nom_lv$beta
			Xbeta_cache <- Xbeta(X, beta)
			ab_outer_cache <- outer(a, b, "+")
			EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			nom_lvl_stats <- nom_lvl_stats + c(1, nom_lv$calls, nom_lv$free)
		}

		# ridge of a row or column covariate with the additive effects of its margin (R/covariate_ridge_move.R)
		if (cr_on) {
			cr_o <- .lame_cov_ridge_visit(
				list(X),
				as.numeric(beta),
				cr_P0,
				a,
				b,
				NULL,
				NULL,
				Sab,
				NA_real_,
				NA_real_,
				TRUE,
				TRUE,
				center = FALSE
			)
			beta[] <- cr_o$beta
			a <- cr_o$a_lev; b <- cr_o$b_lev
			Xbeta_cache <- Xbeta(X, beta)
			ab_outer_cache <- outer(a, b, "+")
			EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			cr_stats[c("visits", names(cr_o$stats))] <- cr_stats[c("visits", names(cr_o$stats))] + c(1, cr_o$stats)
			rm(cr_o)
		}
		# positions with the latent rows and columns of the actors they bound (R/sparse_position_move.R)
		if (spm_on && (s %% spm_every == 0L)) {
			spm_o <- .lame_spm_visit(
				spm_ctx,
				array(Z, c(n, n, 1L)),
				U,
				V,
				shrink = TRUE,
				uv_prior = spm_uv_prior,
				alpha = if (use_explicit_cutpoints) alpha else NULL
			)
			Z[] <- spm_o$Z
			U <- spm_o$U; V <- spm_o$V
			UV_cache <- U %*% t(V)
			EZ_cache <- EZ_cache + matrix(spm_o$shift, n, n)
			spm_stats[c("visits", names(spm_o$stats))] <- spm_stats[c("visits", names(spm_o$stats))] + c(1, spm_o$stats)
			rm(spm_o)
		}

		# collapsed per-actor update with the latent array integrated out
		# (a symmetric fit draws its positions only after the first half of burn-in)
		if (actor_move && (s %% am_every == 0L) && (any(am_lev) || !symmetric || s > .5 * burn)) {
			Z3 <- array(Z, c(n, n, 1L))
			EZ3 <- array(EZ_cache, c(n, n, 1L))
			if (symmetric) {
				# symmetric positions are V = U diag(L); a state not of that form leaves the positions to the sampler
				am_Ld <- NULL
				if (R > 0L) {
					am_Ld <- vapply(seq_len(R), function(r) {
						pv <- which.max(abs(U[, r]))
						if (U[pv, r] == 0) 1 else V[pv, r] / U[pv, r]
					}, numeric(1))
					if (!all(is.finite(am_Ld)) ||
							max(abs(V - U %*% diag(am_Ld, R))) > 1e-8 * max(1, max(abs(V))))
						am_Ld <- NULL
				}
				am_out <- .lame_actor_move_visit_sym(
					Z3,
					EZ3,
					array(Y, c(n, n, 1L)),
					rho,
					beta,
					a,
					NULL,
					if (is.null(am_Ld)) NULL else U,
					am_Ld,
					Sab[1L, 1L],
					NULL,
					(s > .5 * burn),
					am_cx,
					centred = FALSE,
					lev = am_lev[1L]
				)
				a <- am_out$la; b <- a
				if (!is.null(am_Ld)) { U <- am_out$U; V <- U %*% diag(am_Ld, R) }
			} else {
				am_out <- .lame_actor_move_visit(
					Z3,
					EZ3,
					rho,
					beta,
					a,
					b,
					NULL,
					NULL,
					if (R > 0L) U else NULL,
					if (R > 0L) V else NULL,
					Sab,
					NULL,
					TRUE,
					am_cx,
					centred = FALSE,
					uv_prior = am_uv_prior,
					lev = am_lev
				)
				a <- am_out$la; b <- am_out$lb
				if (R > 0L) { U <- am_out$U; V <- am_out$V }
			}
			Z[] <- am_out$Z
			actor_move_stats[c("visits", names(am_out$stats))] <-
				actor_move_stats[c("visits", names(am_out$stats))] + c(1, am_out$stats)
			ab_outer_cache <- outer(a, b, "+")
			UV_cache <- U %*% t(V)
			EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			rm(am_out, Z3, EZ3)
		}

		# translations of the levels with their latent rows and columns (Poisson, ordinal)
		if (ltm_on && (s %% ltm_every == 0L)) {
			ltm_st <- .lame_ltm_visit(
				ltm_cx,
				array(Z, c(n, n, 1L)),
				array(EZ_cache, c(n, n, 1L)),
				a,
				b,
				NULL,
				NULL,
				as.numeric(beta),
				Sab,
				NULL,
				centred = FALSE,
				lev = am_lev,
				alpha = if (use_explicit_cutpoints) alpha else NULL
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
					NULL,
					centred = FALSE,
					lev = am_lev,
					alpha = if (use_explicit_cutpoints) alpha else NULL
				)
				ltm_stats["location_calls"] <- ltm_stats["location_calls"] + 1
			}
			ltm_Z <- ltm_st$Z
			ltm_EZ <- ltm_st$EZ
			ltm_beta <- ltm_st$beta
			if (ltm_coef) {
				ltm_cf <- .lame_ppm_static(ppm_cx, ltm_Z, ltm_EZ, ltm_beta, ltm_Xk, ltm_iV0, Xz = ltm_Xz)
				ltm_Z <- ltm_cf$Z; ltm_EZ <- ltm_cf$EZ; ltm_beta <- ltm_cf$beta
				ltm_stats[names(ltm_cf$stats)] <- ltm_stats[names(ltm_cf$stats)] + ltm_cf$stats
				rm(ltm_cf)
			}
			if (ltm_uv) {
				ltm_uo <- NULL
				if (symmetric) {
					# symmetric positions are V = U diag(L); a state not of that form leaves the positions to the sampler
					ltm_Ld <- vapply(seq_len(R), function(r) {
						pv <- which.max(abs(U[, r]))
						if (U[pv, r] == 0) 1 else V[pv, r] / U[pv, r]
					}, numeric(1))
					if (all(is.finite(ltm_Ld)) && max(
						abs(V - U %*% diag(ltm_Ld, R))
					) <= 1e-8 * max(1, max(abs(V))) && s > .5 * burn) {
						ltm_uo <- .lame_ltm_uv_visit(
							ltm_cx,
							array(ltm_Z, c(n, n, 1L)),
							array(ltm_EZ, c(n, n, 1L)),
							U,
							NULL,
							NULL,
							ltm_Ld,
							shrink = TRUE
						)
						U <- ltm_uo$U; V <- U %*% diag(ltm_Ld, R)
					}
				} else {
					ltm_uo <- .lame_ltm_uv_visit(
						ltm_cx,
						array(ltm_Z, c(n, n, 1L)),
						array(ltm_EZ, c(n, n, 1L)),
						U,
						V,
						uv_prior = ltm_uv_prior
					)
					U <- ltm_uo$U; V <- ltm_uo$V
				}
				if (!is.null(ltm_uo)) {
					ltm_Z <- ltm_uo$Z; ltm_EZ <- ltm_uo$EZ
					ltm_stats[names(ltm_uo$stats)] <- ltm_stats[names(ltm_uo$stats)] + ltm_uo$stats
				}
				rm(ltm_uo)
			}
			if (ltm_scale) {
				# priors of the Poisson scale move: rrho_mh()'s default rho prior and the one-observation
				# inverse gamma of s2, counting the diagonal unless symmetric
				ltm_two <- !symmetric && isTRUE(dcor)
				ltm_sc <- .lame_ltm_scale(
					ltm_cx,
					ltm_Z,
					ltm_EZ,
					s2,
					if (ltm_two) rho else 0,
					c(shape0 = 1 / 2, rate0 = s2_prior_scale / 2),
					two = ltm_two,
					rho_logprior = .lame_rho_prior_log,
					diag_in_target = !symmetric
				)
				ltm_Z <- ltm_sc$Z; s2 <- ltm_sc$s2
				if (ltm_two) rho <- ltm_sc$rho
				ltm_stats[names(ltm_sc$stats)] <- ltm_stats[names(ltm_sc$stats)] + ltm_sc$stats
				rm(ltm_sc)
			}
			Z[] <- ltm_Z
			a <- ltm_st$la
			b <- if (symmetric) a else ltm_st$lb
			if (ltm_coef) {
				beta[] <- ltm_beta
				Xbeta_cache <- Xbeta(X, beta)
			}
			ab_outer_cache <- outer(a, b, "+")
			if (ltm_uv) UV_cache <- U %*% t(V)
			EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			rm(ltm_st, ltm_Z, ltm_EZ, ltm_beta)
		}

		# Poisson dyadic variance and correlation with the latent residuals along their conditional law
		if (ppm_scale_on) {
			ppm_two <- !symmetric && isTRUE(dcor)
			ppm_sc <- .lame_ppm_scale(
				ppm_cx,
				array(Z, c(n, n, 1L)),
				array(EZ_cache, c(n, n, 1L)),
				s2,
				if (ppm_two) rho else 0,
				c(shape0 = 1 / 2, rate0 = s2_prior_scale / 2),
				ppm_two,
				.lame_rho_prior_log,
				diag_in_target = !symmetric
			)
			Z[] <- ppm_sc$Z
			s2 <- ppm_sc$s2
			if (ppm_two) rho <- ppm_sc$rho
			ppm_stats <- ppm_stats + c(1, ppm_sc$stats)
			# the record of the translations is the work of the scale move that ran
			if (ltm_on) {
				ltm_stats[c("scale_attempts", "scale_accepts", "rho_attempts", "rho_accepts")] <-
					ltm_stats[c("scale_attempts", "scale_accepts", "rho_attempts", "rho_accepts")] +
					ppm_sc$stats[c("scale_attempts", "scale_accepts", "scale_rho_attempts", "scale_rho_accepts")]
			}
			rm(ppm_sc)
		}

		# scale of the whole state and the shift and scale of the rank-likelihood latent array
		if (osm_on) {
			osm_ab <- .lame_osm_prior_ab(
				a,
				b,
				Sab,
				if (symmetric) nvar else rvar,
				if (symmetric) nvar else cvar,
				symmetric = symmetric
			)
			osm_o <- .lame_osm_visit(
				array(Z, c(n, n, 1L)),
				array(EZ_cache, c(n, n, 1L)),
				if (R > 0L) array(UV_cache, c(n, n, 1L)) else NULL,
				osm_mode,
				rho,
				s2,
				A_pri = if (is.null(osm_ab)) 0 else drop(crossprod(beta, osm_iV0 %*% beta)) + osm_ab[1L],
				B_pri = 0,
				d_th = if (is.null(osm_ab)) 0 else length(beta) + osm_ab[2L] + (if (use_explicit_cutpoints) length(alpha) - 1L else 0L),
				affine = !use_explicit_cutpoints,
				orbit = !is.null(osm_ab)
			)
			# explicit cutpoints scale with the state (the first stays at zero)
			if (use_explicit_cutpoints && osm_o$c != 1) {
				osm_alpha <- c(alpha[1L], osm_o$c * alpha[-1L])
				if (.lame_osm_ex_ok(osm_o$Z[, , 1L], Y_int, osm_alpha)) alpha <- osm_alpha else osm_o$c <- 1
			}
			if (osm_o$c != 1 || !use_explicit_cutpoints) Z[] <- osm_o$Z
			osm_stats <- osm_stats + c(1, osm_o$c != 1)
			if (osm_o$c != 1) {
				beta[] <- osm_o$c * beta
				a <- osm_o$c * a
				b <- if (symmetric) a else osm_o$c * b
				Xbeta_cache <- osm_o$c * Xbeta_cache
				ab_outer_cache <- osm_o$c * ab_outer_cache
				EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
			}
			rm(osm_o, osm_ab)
		}

		# scale of a binary predictor, with the latent array integrated out
		if (bso_on && (s %% bso_every == 0L)) {
			bso_ab <- .lame_osm_prior_ab(
				a,
				b,
				Sab,
				if (symmetric) nvar else rvar,
				if (symmetric) nvar else cvar,
				symmetric = symmetric
			)
			if (!is.null(bso_ab)) {
				bso_o <- .lame_bso_visit(
					array(Z, c(n, n, 1L)),
					array(EZ_cache, c(n, n, 1L)),
					if (R > 0L) array(UV_cache, c(n, n, 1L)) else NULL,
					bso_cx,
					rho,
					A_pri = drop(crossprod(beta, bso_iV0 %*% beta)) + bso_ab[1L],
					B_pri = 0,
					d_th = length(beta) + bso_ab[2L]
				)
				bso_stats <- bso_stats + c(1, bso_o$evals, bso_o$c != 1)
				if (bso_o$c != 1) {
					Z[] <- bso_o$Z
					beta[] <- bso_o$c * beta
					a <- bso_o$c * a
					b <- if (symmetric) a else bso_o$c * b
					Xbeta_cache <- bso_o$c * Xbeta_cache
					ab_outer_cache <- bso_o$c * ab_outer_cache
					EZ_cache <- Xbeta_cache + ab_outer_cache + UV_cache
				}
				rm(bso_o)
			}
			rm(bso_ab)
		}


		if (verbose) {
			if (s <= burn) {
				cli::cli_progress_update(id = NULL)
				if (s == burn) {
					cli::cli_progress_done()
					cli::cli_progress_bar("Sampling", total = nscan, clear = TRUE)
				}
			} else if (s > burn) {
				cli::cli_progress_update(id = NULL)
			}
		}

		if (s %% odens == 0 & s > burn) {
			sample_idx <- sample_idx + 1

			# store the draw in its identified form: positions centred over actors, with the removed row,
			# column and grand means added to a, b and the intercept (as in lame()); the sampler state is unchanged
			# disabled margins keep raw factor coordinates since they cannot absorb the offsets
			beta_rep <- beta ; a_rep <- a ; b_rep <- b
			U_rep <- U ; V_rep <- V ; UV_rep <- UV_cache
			if (isTRUE(rvar) && isTRUE(cvar) && R > 0 &&
					!is.null(U) && !is.null(V) && NCOL(U) > 0L &&
					all(is.finite(U)) && all(is.finite(V))) {
				cen_rep <- .lame_center_uv_period(U, V, diag(1, ncol(U)))
				if (!is.null(cen_rep)) {
					a_row <- cen_rep$a_row ; b_col <- cen_rep$b_col
					cst <- cen_rep$const
					if (symmetric) {
						# one effect per actor: under V = U L the row and column halves agree, and averaging keeps that under rounding
						a_row <- (a_row + b_col) / 2 ; b_col <- a_row
					}
					if (isTRUE(intercept) && length(beta_rep) > 0L &&
							is.finite(beta_rep[1L])) {
						beta_rep[1L] <- beta_rep[1L] + cst
					} else if (symmetric) {
						# a shared actor effect appears twice in each predictor
						a_row <- a_row + cst / 2
						b_col <- b_col + cst / 2
					} else {
						a_row <- a_row + cst
					}
					a_rep <- a_rep + a_row ; b_rep <- b_rep + b_col
					U_rep <- cen_rep$U ; V_rep <- cen_rep[["V"]]
					UV_rep <- U_rep %*% t(V_rep)
				}
			}

			# summaries the chain policy compares across the run
			.chain_acc <- .lame_chain_acc_update(
				.chain_acc,
				sample_idx,
				floor(nscan / odens),
				a_rep,
				b_rep,
				if (R > 0) U_rep else NULL,
				if (R > 0) V_rep else NULL,
				NULL,
				symmetric = symmetric,
				actors = list(rownames(Y), colnames(Y))
			)

			if (symmetric) {
				br <- beta_rep[rb] ; bc <- beta_rep[cb] ; bn <- (br + bc) / 2
				# same `-c(0, integer(0))` guard as in the colnames(beta) build
				excl_idx <- c(if (intercept) 1L else integer(0), rb, cb)
				bint <- if (intercept) beta_rep[1L] else numeric(0)
				bdyad <- if (length(excl_idx) > 0L) beta_rep[-excl_idx] else beta_rep
				sbeta <- c(bint, bn, bdyad)
				BETA[sample_idx, ] <- sbeta

				VC[sample_idx, ] <- c(Sab[1, 1], s2)

				beta_running_mean <- beta_running_mean + (sbeta - beta_running_mean) / sample_idx
				vc_running_mean <- vc_running_mean + (c(Sab[1, 1], s2) - vc_running_mean) / sample_idx
			}

			if (!symmetric) {
				BETA[sample_idx, ] <- beta_rep
				VC[sample_idx, ] <- c(Sab[upper.tri(Sab, diag = TRUE)], rho, s2)

				beta_running_mean <- beta_running_mean + (beta_rep - beta_running_mean) / sample_idx
				vc_running_mean <- vc_running_mean + (c(
					Sab[upper.tri(Sab, diag = TRUE)],
					rho,
					s2
				) - vc_running_mean) / sample_idx
			}

			# explicit cutpoints alpha_2..alpha_{k-1} (alpha_1 = 0 is the anchor); none to store when k = 2
			if (use_explicit_cutpoints && !is.null(ALPHA)) {
				ALPHA[sample_idx, ] <- alpha[-1L]
			}

			if (!is.null(posterior_opts)) {
				if (posterior_opts$save_UV && R > 0) {
					if (sample_idx %% posterior_opts$thin_UV == 0) {
						UV_sample_idx <- UV_sample_idx + 1
						if (UV_sample_idx <= dim(U_samples)[3]) {
							U_samples[, , UV_sample_idx] <- U_rep
							if (!is.null(V_samples)) {
								V_samples[, , UV_sample_idx] <- V_rep
							}
						}
					}
				}

				if (posterior_opts$save_ab) {
					if (sample_idx %% posterior_opts$thin_ab == 0) {
						ab_sample_idx <- ab_sample_idx + 1
						if (rvar && !is.null(a_samples) && ab_sample_idx <= ncol(a_samples)) {
							a_samples[, ab_sample_idx] <- a_rep
						}
						if (cvar && !is.null(b_samples) && ab_sample_idx <= ncol(b_samples)) {
							b_samples[, ab_sample_idx] <- b_rep
						}
					}
				}
			}

			running_sample_count <- sample_idx

			UVPS <- UVPS + UV_rep
			APS <- APS + a_rep
			BPS <- BPS + b_rep
			APS2 <- APS2 + a_rep^2
			BPS2 <- BPS2 + b_rep^2

			EZ <- EZ_cache
			if (symmetric) { EZ <- (EZ + t(EZ)) / 2 }

			rho_eff <- rho
			Ys <- simulate_Y_fn(EZ, rho_eff, s2)

			if (symmetric) { Ys[lower.tri(Ys)] <- 0 ; Ys <- Ys + t(Ys) }

			YPS <- YPS + Ys

			if (gof) {
				Ys[is.na(Y)] <- NA
				gof_idx <- gof_idx + 1
				GOF[gof_idx, 1:5] <- gof_stats(Ys, mode = "unipartite")

				if (!is.null(custom_gof)) {
					if (is.function(custom_gof)) {
						tryCatch({
							custom_vals <- custom_gof(Ys)
							GOF[gof_idx, (6):(5 + length(custom_vals))] <- custom_vals
						}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
					} else if (is.list(custom_gof)) {
						for (i in seq_along(custom_gof)) {
							if (is.function(custom_gof[[i]])) {
								tryCatch({
									GOF[gof_idx, 5 + i] <- custom_gof[[i]](Ys)
								}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
							}
						}
					}
				}
			}

		}

		if (periodic_save & s %in% save_points & !is.null(out_file)) {
			start_vals <- list(Z = Z, beta = beta, a = a, b = b, U = U, V = V, rho = rho, s2 = s2, Sab = Sab)

			fit <- list(
				APM = APS / (sum(s > burn & s %% odens == 0)),
				BPM = BPS / (sum(s > burn & s %% odens == 0)),
				BETA = BETA,
				VC = VC,
				GOF = GOF,
				X = X,
				Y = if (exists("Y_user", inherits = FALSE)) Y_user else Y,
				start_vals = start_vals,
				symmetric = symmetric,
				model.name = model.name
			)
			save(fit, file = out_file)
			rm(list = c('fit', 'start_vals'))
		}

	} # end mcmc

	# output assembly
	if (verbose) {
		cli::cli_progress_done()

		if (exists("start_time")) {
			end_time <- Sys.time()
			elapsed <- difftime(end_time, start_time, units = "secs")
			cli::cli_alert_success("MCMC complete in {.val {round(elapsed, 1)}} seconds")
		}
	}

	# trim matrices to the actual sample size
	if (exists("sample_idx") && sample_idx > 0) {
		BETA <- BETA[1:sample_idx, , drop = FALSE]
		VC <- VC[1:sample_idx, , drop = FALSE]
		if (use_explicit_cutpoints && !is.null(ALPHA)) {
			ALPHA <- ALPHA[1:sample_idx, , drop = FALSE]
		}
		if (gof && exists("gof_idx")) {
			GOF <- GOF[1:gof_idx, , drop = FALSE]
		}
	}

	APM <- APS / nrow(VC)
	BPM <- BPS / nrow(VC)
	APM_sd <- sqrt(pmax(APS2 / nrow(VC) - APM^2, 0))
	BPM_sd <- sqrt(pmax(BPS2 / nrow(VC) - BPM^2, 0))
	if (!is.null(rownames(Y))) {
		names(APM) <- rownames(Y)
		names(BPM) <- rownames(Y)
	} else {
		names(APM) <- paste0("Actor", 1:length(APM))
		names(BPM) <- paste0("Actor", 1:length(BPM))
	}
	UVPM <- UVPS / nrow(VC)
	YPM <- YPS / nrow(VC)
	names(APM) <- names(BPM) <- rownames(UVPM) <- colnames(UVPM) <- dimnames(Y)[[1]]

	dimnames(YPM) <- dimnames(Y)
	rownames(BETA) <- NULL

	if (!symmetric) {
		if (R > 0) {
			UDV <- svd(UVPM)
			U <- UDV$u[, seq(1, R, length = R)] %*% diag(sqrt(UDV$d)[seq(1, R, length = R)], nrow = R)
			V <- UDV$v[, seq(1, R, length = R)] %*% diag(sqrt(UDV$d)[seq(1, R, length = R)], nrow = R)
			rownames(U) <- rownames(V) <- rownames(Y)
		} else {
			U <- NULL
			V <- NULL
		}
		start_vals_final <- list(Z = Z, beta = beta, a = a, b = b, U = U, V = V, rho = rho, s2 = s2, Sab = Sab)
		if (use_explicit_cutpoints) start_vals_final <- c(start_vals_final, list(alpha = alpha, tau_prop = tau_prop))
		# resolved g and prior (with sampler defaults filled in), for prior_summary()
		.resolved_prior <- prior
		fit <- list(
			BETA = BETA,
			VC = VC,
			APM = APM,
			BPM = BPM,
			U = U,
			V = V,
			L = NULL,
			UVPM = UVPM,
			YPM = YPM,
			GOF = GOF,
			X = X,
			Y = if (exists("Y_user", inherits = FALSE)) Y_user else Y,
			start_vals = start_vals_final,
			model.name = model.name,
			mode = "unipartite",
			family = family,
			symmetric = symmetric,
			odmax = odmax,
			R = R,
			n = n,
			p = dim(X)[3],
			rvar = rvar,
			cvar = cvar,
			dcor = dcor,
			g = g,
			prior = .resolved_prior,
			tryErrorChecks = tryErrorChecks,
			mh_counters = tryErrorChecks,
			X_names = colnames(BETA),
			actor_names = rownames(Y)
		)

		if (!is.null(posterior_opts)) {
			if (posterior_opts$save_UV && !is.null(U_samples) && UV_sample_idx > 0) {
				actual_UV_samples <- min(UV_sample_idx, dim(U_samples)[3])
				if (actual_UV_samples > 0) {
					fit$U_samples <- U_samples[, , 1:actual_UV_samples, drop = FALSE]
					if (!is.null(V_samples)) {
						fit$V_samples <- V_samples[, , 1:actual_UV_samples, drop = FALSE]
					}
				}
			}
			if (posterior_opts$save_ab && ab_sample_idx > 0) {
				if (rvar && !is.null(a_samples)) {
					actual_ab_samples <- min(ab_sample_idx, ncol(a_samples))
					if (actual_ab_samples > 0) {
						fit$a_samples <- a_samples[, 1:actual_ab_samples, drop = FALSE]
					}
				}
				if (cvar && !is.null(b_samples)) {
					actual_ab_samples <- min(ab_sample_idx, ncol(b_samples))
					if (actual_ab_samples > 0) {
						fit$b_samples <- b_samples[, 1:actual_ab_samples, drop = FALSE]
					}
				}
			}
		}
	}

	if (symmetric) {
		ULUPM <- UVPM
		if (R > 0) {
			eULU <- eigen(ULUPM)
			# keep the R eigenvalues largest in magnitude, so U L U' is the rank-R truncation of ULUPM
			eR <- which(rank(-abs(eULU$val), ties.method = "first") <= R)
			U <- eULU$vec[, eR, drop = FALSE]
			L <- diag(eULU$val[eR], nrow = R)
			rownames(U) <- rownames(ULUPM)
		} else {
			U <- NULL
			L <- NULL
		}
		YPM <- .5 * (YPM + t(YPM))
		start_vals_final <- list(Z = Z, beta = beta, a = a, b = b, U = U, V = U, rho = rho, s2 = s2, Sab = Sab)
		if (use_explicit_cutpoints) start_vals_final <- c(start_vals_final, list(alpha = alpha, tau_prop = tau_prop))
		fit <- list(
			BETA = BETA,
			VC = VC,
			APM = APM,
			BPM = NULL,
			U = U,
			V = NULL,
			L = L,
			ULUPM = ULUPM,
			YPM = YPM,
			GOF = GOF,
			X = X,
			Y = if (exists("Y_user", inherits = FALSE)) Y_user else Y,
			start_vals = start_vals_final,
			model.name = model.name,
			mode = "unipartite",
			family = family,
			symmetric = symmetric,
			odmax = odmax,
			R = R,
			call = match.call(),
			n = n,
			p = dim(X)[3],
			rvar = rvar,
			cvar = cvar,
			dcor = dcor,
			g = g,
			prior = prior, # parity with the asymmetric path
			tryErrorChecks = tryErrorChecks,
			mh_counters = tryErrorChecks,
			X_names = colnames(BETA),
			actor_names = rownames(Y)
		)

		if (!is.null(posterior_opts)) {
			if (posterior_opts$save_UV && !is.null(U_samples) && UV_sample_idx > 0) {
				actual_UV_samples <- min(UV_sample_idx, dim(U_samples)[3])
				if (actual_UV_samples > 0) {
					fit$U_samples <- U_samples[, , 1:actual_UV_samples, drop = FALSE]
					if (!is.null(V_samples)) {
						fit$V_samples <- V_samples[, , 1:actual_UV_samples, drop = FALSE]
					}
				}
			}
			if (posterior_opts$save_ab && ab_sample_idx > 0) {
				if (rvar && !is.null(a_samples)) {
					actual_ab_samples <- min(ab_sample_idx, ncol(a_samples))
					if (actual_ab_samples > 0) {
						fit$a_samples <- a_samples[, 1:actual_ab_samples, drop = FALSE]
					}
				}
				if (cvar && !is.null(b_samples)) {
					actual_ab_samples <- min(ab_sample_idx, ncol(b_samples))
					if (actual_ab_samples > 0) {
						fit$b_samples <- b_samples[, 1:actual_ab_samples, drop = FALSE]
					}
				}
			}
		}
	}

	if (!is.null(fit$RHO)) fit$RHO <- as.numeric(fit$RHO)
	if (!is.null(fit$s2)) fit$s2 <- as.numeric(fit$s2)

	# explicit-cutpoint posterior and bookkeeping; alpha_post_mean prepends the anchor 0
	fit$ordinal_cutpoints <- ordinal_cutpoints
	if (use_explicit_cutpoints) fit$ordinal_cutpoint_update <-
		if (!symmetric && dcor) "conditional_gibbs" else "marginal_cowles"
	if (use_explicit_cutpoints) {
		fit$ordinal_levels <- ord_lvls
		if (!is.null(ALPHA)) {
			fit$ALPHA <- ALPHA
			fit$alpha_post_mean <- c(0, colMeans(ALPHA, na.rm = TRUE))
		} else {
			fit$alpha_post_mean <- 0
		}
	}

	# one post-MCMC warning for errors raised by a user custom_gof
	if (exists(".cg_tracker", inherits = FALSE))
		.maybe_warn_custom_gof(.cg_tracker, where = "GOF")

	if (!is.null(s2_prior_scale)) fit$s2_prior_scale <- s2_prior_scale
	class(fit) <- "ame"
	# Poisson fitted count at the posterior means of the predictor and s2 (see .poisson_fitted_mean)
	if (identical(family, "poisson")) fit$YPM <- .poisson_fitted_from_fit(fit)
	# chain summary, and posterior sd of the additive effects for simulate()
	chain_summary <- .lame_chain_acc_finish(.chain_acc)
	if (!is.null(chain_summary)) fit$chain_summary <- chain_summary
	fit$APM_sd <- stats::setNames(APM_sd, names(APM))
	fit$BPM_sd <- stats::setNames(BPM_sd, names(BPM))

	if (use_sparse_matrices) {
		fit <- compact_ame(fit, use_sparse_matrices = use_sparse_matrices)
	}

	if (rho_collapsed_enabled) fit$rho_collapsed_move <- c(
		attempts = rho_collapsed_attempts,
		accepts = rho_collapsed_accepts,
		rate = if (rho_collapsed_attempts) rho_collapsed_accepts / rho_collapsed_attempts else NA_real_
	)
	# collapsed per-actor update record: visits and slice likelihood evaluations
	if (actor_move && actor_move_stats["visits"] > 0) fit$actor_move <- actor_move_stats
	if (ltm_on && ltm_stats["visits"] > 0) fit$actor_move <- ltm_stats
	if (uv_split_stats[["visits"]] > 0) fit$uv_split_move <- uv_split_stats
	if (spm_on && spm_stats[["visits"]] > 0) fit$sparse_position_move <- spm_stats
	if (cr_on && cr_stats[["ridge_moved"]] > 0) fit$cov_ridge_move <- cr_stats
	# rank-likelihood scale record: visits and visits at which the whole-state scale moved
	if (osm_on && osm_stats[["visits"]] > 0) fit$ordinal_scale_move <- osm_stats
	if (bso_on && bso_stats[["visits"]] > 0) fit$binary_scale_move <- bso_stats
	if (ppm_scale_on && ppm_stats[["visits"]] > 0) fit$poisson_period_move <- ppm_stats
	# receiver move record of the nomination families: moves made, columns drawn and columns left to the sampler
	if ((nom_col_move && nom_col_stats["visits"] > 0) || (nom_lvl_move && nom_lvl_stats["visits"] > 0) || (nom_rho_move && nom_rho_stats[["calls"]] > 0)) {
		nom_rec <- if (nom_col_move) nom_col_stats else c(visits = nom_lvl_stats[["visits"]])
		if (nom_lvl_move) nom_rec <- c(
			nom_rec,
			level_calls = nom_lvl_stats[["calls"]],
			level_rows = nom_lvl_stats[["free_rows"]]
		)
		if (nom_rho_move && nom_rho_stats[["calls"]] > 0) nom_rec <- c(
			nom_rec,
			rho_attempts = nom_rho_stats[["attempts"]],
			rho_accepts = nom_rho_stats[["accepts"]]
		)
		fit$actor_move <- nom_rec
	}
	fit$sampler_health <- .lame_sampler_health(tryErrorChecks)
	fit$sampler_health_iterations <- s
	.lame_report_sampler_health(fit$sampler_health, fit$sampler_health_iterations)
	.lame_separation_check(fit)
	fit
}


# prior scale of the additive and multiplicative effects of a Poisson fit: the variance across actors
# of the centred log row / column rates, with a half-count pseudo-observation (1 when undefined)
.ame_pois_effect_scale <- function(Y) {
	# Y may be a matrix or an [n, m, t] array; totals run over every dimension but the actor's own
	Y <- pmax(Y, 0)
	obs <- !is.na(Y)
	a_tmp <- log((apply(Y, 1, sum, na.rm = TRUE) + 0.5) / (apply(obs, 1, sum) + 1))
	b_tmp <- log((apply(Y, 2, sum, na.rm = TRUE) + 0.5) / (apply(obs, 2, sum) + 1))
	a_tmp <- a_tmp - mean(a_tmp)
	b_tmp <- b_tmp - mean(b_tmp)
	vscale <- mean(c(stats::var(a_tmp), stats::var(b_tmp)), na.rm = TRUE)
	if (!is.finite(vscale) || vscale <= 0) vscale <- 1
	vscale
}

# starting latent log-rates log(y + c), c the mean count capped at 1; missing cells and a
# unipartite diagonal start at the mean log-rate
.ame_pois_init_z <- function(Y, unipartite = TRUE) {
	ybar <- mean(Y, na.rm = TRUE)
	if (!is.finite(ybar) || ybar <= 0) ybar <- 1
	cc <- min(ybar, 1)
	Z <- log(pmax(Y, 0) + cc)
	Z[is.na(Z)] <- log(cc)
	if (unipartite && nrow(Z) == ncol(Z)) diag(Z) <- log(cc)
	Z
}
