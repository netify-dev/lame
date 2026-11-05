# recovery grid: simulate from the model's own generative process with known parameters, fit, score
# recovery and check the post-processing invariants. test_recovery_grid.R runs a fast subset; the full
# grid runs from the source tree with
#   source("tests/testthat/helper_recovery_grid.R")
#   df <- grid_run(grid_cells("G1_static_uni"), seeds = 1:2, out_csv = "g1.csv")
#   grid_systematic(df)
# for the groups G1_static_uni, G2_bipartite, G3_dynamic_uni, G4_als, G5_kinds (the snap and t
# transition kinds) and G6_dynrho; see inst/recovery/README.md
`%||%` = function(a, b) if (is.null(a)) b else a

# cell definition

grid_cell = function(
	id,
	family,
	estimator = "lame",
	mode = "unipartite",
	symmetric = FALSE,
	T = 4L,
	R = 0L,
	R_row = 0L,
	R_col = 0L,
	dynamic_ab = FALSE,
	dynamic_uv = FALSE,
	dynamic_beta = FALSE,
	dynamic_beta_kind = "ar1",
	dynamic_G = FALSE,
	dcor = FALSE,
	intercept = TRUE,
	missing = 0,
	density = "dense",
	n = 30L,
	nA = 24L,
	nB = 18L,
	odmax = 4L,
	nscan = 2000L,
	burn = 500L,
	odens = 10L,
	expect_error = NULL,
	group = "",
	notes = "",
	dynamic_uv_kind = "ar1",
	dynamic_rho = FALSE,
	b0 = NULL,
	rate = NULL,
	sd_a = NULL
) {
	if (estimator %in% c("ame", "ame_als")) T = 1L
	# the ordinal latent location is not identified: the package default is
	# no intercept and the level scores are not meaningful
	if (identical(family, "ordinal")) intercept = FALSE
	# a nomination cap that saturates every row leaves the sender effects
	# weakly identified; a third of the possible targets keeps most rows below it
	if (family %in% c("cbin", "frn") && identical(odmax, 4L)) odmax = 10L
	list(
		id = id,
		family = family,
		estimator = estimator,
		mode = mode,
		symmetric = symmetric,
		T = as.integer(T),
		R = as.integer(R),
		R_row = as.integer(R_row),
		R_col = as.integer(R_col),
		dynamic_ab = dynamic_ab,
		dynamic_uv = dynamic_uv,
		dynamic_beta = dynamic_beta,
		dynamic_beta_kind = dynamic_beta_kind,
		dynamic_G = dynamic_G,
		dcor = dcor,
		intercept = intercept,
		missing = missing,
		density = density,
		n = as.integer(n),
		nA = as.integer(nA),
		nB = as.integer(nB),
		odmax = as.integer(odmax),
		nscan = as.integer(nscan),
		burn = as.integer(burn),
		odens = as.integer(odens),
		expect_error = expect_error,
		group = group,
		notes = notes,
		dynamic_uv_kind = dynamic_uv_kind,
		dynamic_rho = dynamic_rho,
		b0 = b0,
		rate = rate,
		sd_a = sd_a
	)
}

# simulation

.ar1_path = function(n, TT, rho, sd_stat) {
	# stationary AR(1) paths, one row per unit
	out = matrix(0, n, TT)
	out[, 1] = rnorm(n, 0, sd_stat)
	if (TT > 1) for (t in 2:TT) {
		out[, t] = rho * out[, t - 1] + rnorm(n, 0, sd_stat * sqrt(1 - rho^2))
	}
	out
}

# the other two transition kinds: t has innovation variance sigma^2 / lambda with lambda ~ Gamma(nu/2, nu/2),
# and the (nu - 2) / nu factor keeps the stationary variance at sd_stat^2 so the product scale matches the ar1
# cells; snap resets the position to a fresh stationary draw with probability pi_snap
.ar1_path_t = function(n, TT, rho, sd_stat, nu = 4) {
	out = matrix(0, n, TT)
	out[, 1] = rnorm(n, 0, sd_stat)
	sd_innov = sd_stat * sqrt(1 - rho^2) * sqrt((nu - 2) / nu)
	if (TT > 1) for (t in 2:TT) {
		lam = rgamma(n, nu / 2, rate = nu / 2)
		out[, t] = rho * out[, t - 1] + rnorm(n, 0, sd_innov / sqrt(lam))
	}
	out
}
# snap paths of all R coordinates at once: one reset indicator per actor-period replaces the whole position
# with a fresh draw of sd reset_sd (the stationary sd unless the cell sets snap_reset_sd)
.uv_paths_snap = function(n, R, TT, rho, sd_stat, pi_snap = 0.15, reset_sd = sd_stat) {
	out = array(0, c(n, R, TT))
	out[, , 1] = rnorm(n * R, 0, sd_stat)
	if (TT > 1) for (t in 2:TT) {
		drift = rho * matrix(out[, , t - 1], n, R) + rnorm(n * R, 0, sd_stat * sqrt(1 - rho^2))
		snap = runif(n) < pi_snap
		fresh = matrix(rnorm(n * R, 0, reset_sd), n, R)
		drift[snap, ] = fresh[snap, ]
		out[, , t] = drift
	}
	out
}
.uv_path = function(kind, n, TT, rho, sd_stat) {
	switch(kind, t = .ar1_path_t(n, TT, rho, sd_stat), .ar1_path(n, TT, rho, sd_stat))
}

grid_sim = function(cell, seed) {
	set.seed(seed)
	uni = identical(cell$mode, "unipartite")
	sym = isTRUE(cell$symmetric)
	fam = cell$family
	n_a = if (uni) cell$n else cell$nA
	n_b = if (uni) cell$n else cell$nB
	# zero-padded names keep the alphabetical order the fitters sort by equal
	# to the generating order at any size
	wd = max(2L, nchar(max(n_a, n_b)))
	rn = if (uni) sprintf("a%0*d", wd, seq_len(n_a)) else sprintf("r%0*d", wd, seq_len(n_a))
	cn = if (uni) rn else sprintf("c%0*d", wd, seq_len(n_b))
	is_lame = cell$estimator %in% c("lame", "lame_als", "lame_dynamic_als")
	TT = cell$T
	Tsim = TT + as.integer(is_lame) # one extra period for forecast scoring

	# regression truth
	b1 = 0.7
	# a cell may set its own generating intercept (the density of a binary cell)
	b0 = cell$b0 %||% switch(
		fam,
		normal = 0.2,
		binary = if (cell$density == "sparse") -2.1 else -0.5,
		poisson = 0.3,
		ordinal = 0,
		cbin = 0.2,
		frn = 0.2
	)
	if (!isTRUE(cell$intercept)) b0 = 0
	# a target event rate per actor-period (directed binary or cbin, no product): the intercept that puts the
	# expected events per row at the rate given the latent variance; a cbin cap above the rate rarely binds
	if (!is.null(cell$rate) && fam %in% c("binary", "cbin") && uni && !sym) {
		b0 = stats::qnorm(cell$rate / (n_b - 1)) * sqrt(b1^2 + 2 * 0.5^2 + 1)
	}
	beta_path = rep(b1, Tsim)
	if (!isFALSE(cell$dynamic_beta)) {
		if (identical(cell$dynamic_beta_kind, "rw1")) {
			beta_path = b1 + cumsum(c(0, rnorm(Tsim - 1, 0, 0.15)))
		} else {
			beta_path = b1 + drop(.ar1_path(1, Tsim, 0.7, 0.3))
		}
	}

	# additive effects: level plus an optional AR(1) deviation
	sd_a = cell$sd_a %||% 0.5
	a_lev = rnorm(n_a, 0, sd_a)
	b_lev = if (sym) a_lev else rnorm(n_b, 0, sd_a)
	rho_ab = 0.7
	sd_dev = 0.3
	a_path = matrix(a_lev, n_a, Tsim)
	b_path = matrix(b_lev, n_b, Tsim)
	if (isTRUE(cell$dynamic_ab)) {
		a_path = a_path + .ar1_path(n_a, Tsim, rho_ab, sd_dev)
		b_path = if (sym) a_path else b_path + .ar1_path(n_b, Tsim, rho_ab, sd_dev)
	}

	# multiplicative term
	R_a = if (uni) cell$R else cell$R_row
	R_b = if (uni) cell$R else cell$R_col
	has_uv = R_a > 0L && R_b > 0L
	rho_uv = 0.8
	sd_uv = 0.6
	product = vector("list", Tsim)
	U_path = V_path = G_path = NULL
	L = NULL
	if (has_uv) {
		U0 = matrix(rnorm(n_a * R_a, 0, sd_uv), n_a, R_a)
		V0 = if (sym) U0 else matrix(rnorm(n_b * R_b, 0, sd_uv), n_b, R_b)
		if (sym) L = rep(c(1, -1), length.out = R_a)
		G0 = if (uni) (if (sym) diag(L, R_a) else diag(1, R_a)) else diag(1, R_a, R_b)
		U_path = array(U0, c(n_a, R_a, Tsim))
		V_path = array(V0, c(n_b, R_b, Tsim))
		G_path = array(G0, c(R_a, R_b, Tsim))
		if (isTRUE(cell$dynamic_uv)) {
			kind = cell$sim_kind %||% cell$dynamic_uv_kind %||% "ar1"
			if (identical(kind, "snap")) {
				reset_sd = cell$snap_reset_sd %||% sd_uv
				U_path = .uv_paths_snap(n_a, R_a, Tsim, rho_uv, sd_uv, reset_sd = reset_sd)
				V_path = if (sym) U_path else .uv_paths_snap(n_b, R_b, Tsim, rho_uv, sd_uv, reset_sd = reset_sd)
			} else {
				for (r in seq_len(R_a)) U_path[, r, ] = .uv_path(kind, n_a, Tsim, rho_uv, sd_uv)
				if (sym) V_path = U_path else for (r in seq_len(R_b)) V_path[, r, ] = .uv_path(
					kind,
					n_b,
					Tsim,
					rho_uv,
					sd_uv
				)
			}
		}
		for (t in seq_len(Tsim)) {
			product[[t]] = U_path[, , t, drop = TRUE] %*% matrix(G_path[, , t], R_a, R_b) %*% t(
				matrix(V_path[, , t], n_b, R_b)
			)
			product[[t]] = matrix(product[[t]], n_a, n_b)
		}
	} else {
		for (t in seq_len(Tsim)) product[[t]] = matrix(0, n_a, n_b)
	}

	# noise
	s2 = if (fam == "poisson") 0.3 else 1
	rho_e = if (isTRUE(cell$dcor) && uni && !sym) 0.5 else 0
	# a period-varying dyadic correlation for the dynamic_rho cells
	rho_e_path = if (isTRUE(cell$dynamic_rho) && rho_e > 0 && !isTRUE(cell$flat_rho)) {
		0.5 + 0.3 * sin(2 * pi * (seq_len(Tsim) - 1) / TT)
	} else rep(rho_e, Tsim)

	X = Y = eta = E = vector("list", Tsim)
	for (t in seq_len(Tsim)) {
		x = matrix(rnorm(n_a * n_b), n_a, n_b)
		if (sym) x = (x + t(x)) / 2
		if (uni) diag(x) = NA
		et = b0 + beta_path[t] * x + outer(a_path[, t], b_path[, t], "+") + product[[t]]
		if (uni && sym) et = (et + t(et)) / 2
		# dyadic noise
		if (uni && sym) {
			e = matrix(0, n_a, n_b)
			e[upper.tri(e)] = rnorm(n_a * (n_a - 1) / 2, 0, sqrt(s2))
			e = e + t(e)
		} else if (uni) {
			# directed pairs (e_ij, e_ji) with correlation rho_e
			ut = upper.tri(matrix(0, n_a, n_b))
			m = sum(ut)
			z1 = rnorm(m); z2 = rnorm(m)
			e = matrix(0, n_a, n_b); e[ut] = sqrt(s2) * z1
			e_low = matrix(0, n_a, n_b)
			e_low[ut] = sqrt(s2) * (rho_e_path[t] * z1 + sqrt(1 - rho_e_path[t]^2) * z2)
			e = e + t(e_low)
		} else {
			e = matrix(rnorm(n_a * n_b, 0, sqrt(s2)), n_a, n_b)
		}
		z = et + e
		y = switch(
			fam,
			normal = z,
			binary = 1 * (z > 0),
			poisson = {
				yy = matrix(rpois(n_a * n_b, exp(pmin(ifelse(is.na(z), 0, z), 6))), n_a, n_b)
				# symmetric data carry one count per unordered pair
				if (uni && sym) yy[lower.tri(yy)] = t(yy)[lower.tri(yy)]
				yy
			},
			ordinal = z, # cut below with pooled thresholds
			cbin = ,
			frn = {
				yy = matrix(0, n_a, n_b)
				for (i in seq_len(n_a)) {
					zi = z[i, ]
					if (uni) zi[i] = -Inf
					zi[is.na(zi)] = -Inf
					top = order(zi, decreasing = TRUE)[seq_len(min(cell$odmax, n_b))]
					top = top[zi[top] > 0]
					if (length(top)) {
						if (fam == "cbin") yy[i, top] = 1 else yy[i, top] = rank(zi[top])
					}
				}
				yy
			}
		)
		if (uni) { diag(y) = NA; diag(et) = NA; diag(e) = NA }
		dimnames(y) = list(rn, cn); dimnames(et) = list(rn, cn)
		X[[t]] = array(x, c(n_a, n_b, 1), dimnames = list(rn, cn, "x1"))
		Y[[t]] = y; eta[[t]] = et; E[[t]] = e
	}
	if (fam == "ordinal") {
		zall = unlist(Y)
		cuts = quantile(zall, c(0.25, 0.5, 0.75), na.rm = TRUE)
		Y = lapply(Y, function(z) {
			y = matrix(as.numeric(cut(z, c(-Inf, cuts, Inf), labels = FALSE)), n_a, n_b)
			if (uni) diag(y) = NA
			dimnames(y) = list(rn, cn); y
		})
	}
	# symmetric outcomes are symmetric by construction
	if (uni && sym) for (t in seq_len(Tsim)) {
		yy = Y[[t]]; d = yy - t(yy); stopifnot(all(is.na(d) | abs(d) < 1e-12))
	}
	# missing cells (MCAR, off-diagonal), on the fitted periods only
	if (cell$missing > 0) for (t in seq_len(TT)) {
		yy = Y[[t]]
		k = which(!is.na(yy))
		drop = sample(k, round(cell$missing * length(k)))
		yy[drop] = NA
		if (uni && sym) { m = is.na(yy); m = m | t(m); yy[m] = NA }
		Y[[t]] = yy
	}

	# identified effects: level plus the row / column means of the product
	a_ident = sapply(seq_len(Tsim), function(t) a_path[, t] + rowMeans(product[[t]]))
	b_ident = sapply(seq_len(Tsim), function(t) b_path[, t] + colMeans(product[[t]]))
	b0_ident = b0 + mean(a_path[, seq_len(TT)]) + mean(b_path[, seq_len(TT)]) +
		mean(sapply(seq_len(TT), function(t) mean(product[[t]])))

	list(
		Y = Y[seq_len(TT)],
		X = X[seq_len(TT)],
		eta = eta[seq_len(TT)],
		E = E[seq_len(TT)],
		Y_fc = if (is_lame) Y[[Tsim]] else NULL,
		X_fc = if (is_lame) X[[Tsim]] else NULL,
		eta_fc = if (is_lame) eta[[Tsim]] else NULL,
		truth = list(
			b0 = b0,
			b0_ident = b0_ident,
			b1 = b1,
			beta_path = beta_path[seq_len(TT)],
			a_lev = a_lev,
			b_lev = b_lev,
			a_path = a_path[, seq_len(TT), drop = FALSE],
			b_path = b_path[, seq_len(TT), drop = FALSE],
			a_ident = matrix(a_ident, n_a)[, seq_len(TT), drop = FALSE],
			b_ident = matrix(b_ident, n_b)[, seq_len(TT), drop = FALSE],
			product = product[seq_len(TT)],
			has_uv = has_uv,
			L = L,
			s2 = s2,
			rho_e = rho_e,
			rho_ab = rho_ab,
			sd_dev = sd_dev,
			sigma_ab_innov = sd_dev * sqrt(1 - rho_ab^2),
			rho_uv = rho_uv,
			sd_uv = sd_uv,
			sigma_uv_innov = sd_uv * sqrt(1 - rho_uv^2) *
				(if (identical(cell$sim_kind %||% cell$dynamic_uv_kind %||% "ar1", "t")) sqrt((4 - 2) / 4) else 1),
			rho_path = rho_e_path[seq_len(TT)],
			uv_kind = cell$sim_kind %||% cell$dynamic_uv_kind %||% "ar1",
			pi_snap = 0.15,
			t_nu = 4
		),
		uni = uni,
		sym = sym,
		n_a = n_a,
		n_b = n_b,
		TT = TT
	)
}

# fitting

grid_fit = function(cell, sim, seed) {
	n_warn = 0L; n_msg = 0L; warn_txt = character(0); msg_txt = character(0)
	fam = cell$family
	uni = sim$uni
	odmax = if (fam %in% c("cbin", "frn")) cell$odmax else NULL
	dyn_beta = if (isFALSE(cell$dynamic_beta)) FALSE else "dyad"
	t0 = Sys.time()
	fit = withCallingHandlers(
		tryCatch({
			if (cell$estimator == "ame") {
				if (uni) {
					ame(
						sim$Y[[1]],
						sim$X[[1]],
						family = fam,
						R = cell$R,
						symmetric = cell$symmetric,
						dcor = cell$dcor,
						intercept = cell$intercept,
						odmax = odmax,
						nscan = cell$nscan,
						burn = cell$burn,
						odens = cell$odens,
						verbose = FALSE,
						plot = FALSE,
						gof = FALSE,
						seed = seed + 11L
					)
				} else {
					ame(
						sim$Y[[1]],
						sim$X[[1]],
						family = fam,
						mode = "bipartite",
						R_row = cell$R_row,
						R_col = cell$R_col,
						intercept = cell$intercept,
						odmax = odmax,
						nscan = cell$nscan,
						burn = cell$burn,
						odens = cell$odens,
						verbose = FALSE,
						plot = FALSE,
						gof = FALSE,
						seed = seed + 11L
					)
				}
			} else if (cell$estimator == "ame_als") {
				ame_als(
					sim$Y[[1]],
					sim$X[[1]],
					family = fam,
					R = if (uni) cell$R else max(cell$R_row, cell$R_col),
					mode = cell$mode,
					symmetric = cell$symmetric,
					verbose = FALSE
				)
			} else if (cell$estimator == "lame_als") {
				# the pooled-static ALS path; MCMC-only arguments are not passed
				args = list(
					sim$Y,
					sim$X,
					family = fam,
					symmetric = cell$symmetric,
					intercept = cell$intercept,
					method = "als",
					verbose = FALSE,
					plot = FALSE,
					seed = seed + 11L
				)
				if (uni) args$R = cell$R else { args$mode = "bipartite"; args$R_row = cell$R_row; args$R_col = cell$R_col }
				do.call(lame, args)
			} else if (cell$estimator == "lame_dynamic_als") {
				args = list(
					sim$Y,
					sim$X,
					family = fam,
					symmetric = cell$symmetric,
					dynamic_ab = cell$dynamic_ab,
					dynamic_uv = cell$dynamic_uv,
					dynamic_beta = !isFALSE(cell$dynamic_beta),
					dynamic_beta_kind = cell$dynamic_beta_kind,
					verbose = FALSE,
					seed = seed + 11L
				)
				if (uni) args$R = cell$R else { args$mode = "bipartite"; args$R_row = cell$R_row; args$R_col = cell$R_col }
				do.call(lame:::lame_dynamic_als, args)
			} else {
				args = list(
					sim$Y,
					sim$X,
					family = fam,
					symmetric = cell$symmetric,
					dcor = cell$dcor,
					intercept = cell$intercept,
					odmax = odmax,
					dynamic_ab = cell$dynamic_ab,
					dynamic_uv = cell$dynamic_uv,
					dynamic_beta = dyn_beta,
					dynamic_beta_kind = cell$dynamic_beta_kind,
					nscan = cell$nscan,
					burn = cell$burn,
					odens = cell$odens,
					verbose = FALSE,
					plot = FALSE,
					gof = FALSE,
					seed = seed + 11L
				)
				if (!identical(cell$dynamic_uv_kind %||% "ar1", "ar1")) args$dynamic_uv_kind = cell$dynamic_uv_kind
				if (isTRUE(cell$dynamic_rho)) args$dynamic_rho = TRUE
				if (uni) args$R = cell$R else {
					args$mode = "bipartite"; args$R_row = cell$R_row; args$R_col = cell$R_col
				}
				do.call(lame, args)
			}
		}, error = function(e) structure(list(message = conditionMessage(e)), class = "grid_error")),
		warning = function(w) { n_warn <<- n_warn + 1L; warn_txt <<- c(warn_txt, conditionMessage(w)); invokeRestart("muffleWarning") },
		message = function(m) { n_msg <<- n_msg + 1L; msg_txt <<- c(msg_txt, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	list(
		fit = fit,
		runtime = as.numeric(difftime(Sys.time(), t0, units = "secs")),
		n_warn = n_warn,
		n_msg = n_msg,
		warn_txt = unique(warn_txt),
		msg_txt = unique(msg_txt)
	)
}

# accessors that hide the ame / lame / ALS layout differences

.as_list = function(x, TT) {
	if (is.null(x)) return(NULL)
	if (is.list(x)) return(x)
	if (length(dim(x)) == 3L) return(lapply(seq_len(dim(x)[3]), function(t) x[, , t]))
	rep(list(x), TT)
}

.get_EZ = function(fit, cell, sim) {
	if (inherits(fit, "ame_als")) return(.as_list(fit$EZ, sim$TT))
	if (cell$estimator == "ame") {
		ez = tryCatch(predict(fit, type = "link"), error = function(e) NULL)
		return(.as_list(ez, 1L))
	}
	.as_list(fit$EZ, sim$TT)
}

.get_YPM = function(fit, sim) .as_list(fit$YPM, sim$TT)

.get_beta_draws = function(fit) {
	B = fit$BETA
	if (is.null(B)) return(NULL)
	B
}

.obs_mask = function(sim) lapply(sim$Y, function(y) !is.na(y))

.cor_safe = function(x, y) {
	ok = is.finite(x) & is.finite(y)
	if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
	cor(x[ok], y[ok])
}

# scoring

grid_score = function(cell, sim, fitted, seed) {
	fit = fitted$fit
	tr = sim$truth
	TT = sim$TT
	row = list(
		id = cell$id,
		group = cell$group,
		seed = seed,
		family = cell$family,
		estimator = cell$estimator,
		mode = cell$mode,
		symmetric = cell$symmetric,
		T = TT,
		R = if (sim$uni) cell$R else paste0(cell$R_row, "x", cell$R_col),
		dynamic_ab = cell$dynamic_ab,
		dynamic_uv = cell$dynamic_uv,
		dynamic_beta = if (isFALSE(cell$dynamic_beta)) "FALSE" else cell$dynamic_beta_kind,
		dynamic_G = cell$dynamic_G,
		dcor = cell$dcor,
		intercept = cell$intercept,
		missing = cell$missing,
		density = cell$density,
		n = sim$n_a,
		runtime_s = round(fitted$runtime, 1),
		n_warn = fitted$n_warn,
		n_msg = fitted$n_msg,
		status = "ok",
		error = "",
		notes = cell$notes
	)
	na_scores = list(
		beta1_est = NA, beta1_sd = NA, beta1_z = NA, beta_path_cor = NA,
		icpt_est = NA, icpt_z = NA, level_gap = NA, ez_cor = NA, ypm_ratio = NA, ypm_gap = NA,
		va_ratio = NA, vb_ratio = NA, ve_ratio = NA, rho_est = NA, rho_gap = NA,
		cor_a = NA, cor_b = NA, cor_a_path = NA, cor_a_dev = NA, uv_cor = NA,
		rho_ab_est = NA, sigma_ab_ratio_innov = NA, sigma_ab_ratio_stat = NA,
		rho_uv_est = NA, sigma_uv_ratio_innov = NA,
		inv_pred = NA, inv_fit = NA, inv_res = NA, inv_coef = NA, inv_ci = NA, inv_sim = NA,
		fc_cor = NA, fc_ratio = NA, inv_summary = NA, inv_latent = NA,
		tryerr = NA, tryerr_frac = NA, terminated_early = NA,
		rho_path_cor = NA, pi_snap_est = NA, vab_period_ratio = NA
	)
	row = c(row, na_scores)
	if (inherits(fit, "grid_error")) {
		row$status = "error"; row$error = substr(fit$message, 1, 200)
		if (!is.null(cell$expect_error) && grepl(cell$expect_error, fit$message)) row$status = "expected_error"
		return(as.data.frame(row, stringsAsFactors = FALSE))
	}
	if (!is.null(cell$expect_error)) { row$status <- "ran_but_expected_error"; }
	is_als = inherits(fit, "ame_als")
	obs = .obs_mask(sim)
	fam = cell$family

	# coefficients
	if (is_als) {
		bb = fit$beta
		if (is.matrix(bb)) bb = rowMeans(bb)
		b1 = unname(if (!is.null(names(bb)) && "x1_dyad" %in% names(bb)) bb["x1_dyad"] else bb[length(bb)])
		# ALS has no posterior sd: the gap is scaled by the MCMC posterior sd
		# of this coefficient on the same design (T = 1, n = 30), shrunk by sqrt(T)
		sd_ref = switch(fam, normal = 0.037, binary = 0.065, poisson = 0.038, 0.05) / sqrt(TT)
		row$beta1_est = b1; row$beta1_sd = sd_ref; row$beta1_z = (b1 - tr$b1) / sd_ref
		row$icpt_est = unname(fit$mu %||% NA)
	} else {
		B = .get_beta_draws(fit)
		if (!is.null(B) && length(dim(B)) == 3L) {
			k = which(dimnames(B)[[2]] == "x1_dyad")
			bp = apply(B[, k, , drop = FALSE], 3, mean); bs = apply(B[, k, , drop = FALSE], 3, sd)
			row$beta1_est = mean(bp)
			row$beta1_z = mean(abs((bp - tr$beta_path) / pmax(bs, 1e-8)))
			row$beta_path_cor = .cor_safe(bp, tr$beta_path)
			if ("intercept" %in% dimnames(B)[[2]]) {
				ki = which(dimnames(B)[[2]] == "intercept")
				row$icpt_est = mean(B[, ki, ])
			}
		} else if (!is.null(B) && "x1_dyad" %in% colnames(B)) {
			row$beta1_est = mean(B[, "x1_dyad"]); row$beta1_sd = sd(B[, "x1_dyad"])
			row$beta1_z = (row$beta1_est - tr$b1) / max(row$beta1_sd, 1e-8)
			if ("intercept" %in% colnames(B)) {
				row$icpt_est = mean(B[, "intercept"])
				# identified only when the reporting is identified: lame, or ame without a product
				if ((cell$estimator == "lame" || !tr$has_uv) && fam != "ordinal") {
					row$icpt_z = (row$icpt_est - tr$b0_ident) / max(sd(B[, "intercept"]), 1e-8)
				}
			}
		}
	}

	# linear predictor and first moment
	EZ = .get_EZ(fit, cell, sim)
	if (!is.null(EZ)) {
		ez = unlist(lapply(seq_len(TT), function(t) EZ[[t]][obs[[t]]]))
		et = unlist(lapply(seq_len(TT), function(t) sim$eta[[t]][obs[[t]]]))
		row$ez_cor = .cor_safe(ez, et)
		row$level_gap = if (fam == "ordinal") NA_real_ else mean(ez, na.rm = TRUE) - mean(et, na.rm = TRUE)
	}
	YPM = .get_YPM(fit, sim)
	if (!is.null(YPM)) {
		yp = unlist(lapply(seq_len(TT), function(t) YPM[[t]][obs[[t]]]))
		yo = unlist(lapply(seq_len(TT), function(t) sim$Y[[t]][obs[[t]]]))
		if (fam == "normal") row$ypm_gap = (mean(yp) - mean(yo)) / sd(yo) else row$ypm_ratio = mean(yp) / mean(yo)
	}

	# effects and variance components
	a_ident_lev = rowMeans(tr$a_ident); b_ident_lev = rowMeans(tr$b_ident)
	APM = fit$APM; BPM = fit$BPM %||% (if (sim$sym) APM else NULL)
	if (is_als) { APM = fit$a; BPM = fit$b %||% (if (sim$sym) APM else NULL) }
	if (!is.null(APM) && length(APM) == sim$n_a) row$cor_a = .cor_safe(APM, a_ident_lev)
	if (!is.null(BPM) && length(BPM) == sim$n_b && !sim$sym) row$cor_b = .cor_safe(BPM, b_ident_lev)
	VC = fit$VC
	if (!is.null(VC) && is.matrix(VC) && !is.null(colnames(VC))) {
		# under dynamic_ab the reported va / vb are the variances of the actor
		# levels; the AR(1) deviation is carried by sigma_ab / rho_ab, so the
		# level plus the stationary deviation variance is the per-period target
		if (isTRUE(cell$dynamic_ab)) {
			pm = sapply(tr$product, rowMeans); pmc = sapply(tr$product, colMeans)
			a_lev_id = tr$a_lev + rowMeans(matrix(pm, nrow = sim$n_a))
			b_lev_id = tr$b_lev + rowMeans(matrix(pmc, nrow = sim$n_b))
			if ("va" %in% colnames(VC)) row$va_ratio = mean(VC[, "va"]) / var(a_lev_id)
			if ("vb" %in% colnames(VC) && !sim$sym) row$vb_ratio = mean(VC[, "vb"]) / var(b_lev_id)
			if ("va" %in% colnames(VC) && !is.null(fit$sigma_ab) && !is.null(fit$rho_ab)) {
				dev_var = mean(fit$sigma_ab, na.rm = TRUE)^2 / (1 - mean(fit$rho_ab, na.rm = TRUE)^2)
				row$vab_period_ratio = (mean(VC[, "va"]) + dev_var) / mean(apply(tr$a_ident, 2, var))
			}
		} else {
			if ("va" %in% colnames(VC)) row$va_ratio = mean(VC[, "va"]) / var(a_ident_lev)
			if ("vb" %in% colnames(VC) && !sim$sym) row$vb_ratio = mean(VC[, "vb"]) / var(b_ident_lev)
		}
		if ("ve" %in% colnames(VC) && fam %in% c("normal", "poisson")) {
			e_real = unlist(lapply(seq_len(TT), function(t) sim$E[[t]][obs[[t]]]))
			row$ve_ratio = mean(VC[, "ve"]) / (if (fam == "normal") var(e_real) else tr$s2)
		}
		if ("rho" %in% colnames(VC) && isTRUE(cell$dcor)) {
			row$rho_est = mean(VC[, "rho"]); row$rho_gap = row$rho_est - tr$rho_e
		}
	}
	# a period-varying dyadic correlation: the posterior mean path against the
	# generating path, scored by the mean absolute gap rather than the pooled one
	if (isTRUE(cell$dynamic_rho) && !is.null(fit$rho_path) && !is.null(tr$rho_path) &&
			length(fit$rho_path) == length(tr$rho_path)) {
		row$rho_path_cor = .cor_safe(fit$rho_path, tr$rho_path)
		row$rho_gap = mean(abs(fit$rho_path - tr$rho_path))
		row$rho_est = mean(fit$rho_path)
	}
	if (identical(cell$dynamic_uv_kind, "snap") && !is.null(fit$pi_snap)) {
		row$pi_snap_est = mean(fit$pi_snap, na.rm = TRUE)
	}
	# multiplicative term
	if (tr$has_uv) {
		P = fit$UVPM %||% fit$ULUPM
		if (!is.null(P)) {
			P = .as_list(P, TT)
			cc = sapply(seq_len(min(TT, length(P))), function(t) .cor_safe(P[[t]][obs[[t]]], tr$product[[t]][obs[[t]]]))
			row$uv_cor = mean(cc, na.rm = TRUE)
		}
	}
	# dynamic pieces
	if (isTRUE(cell$dynamic_ab) && !is.null(fit$a_dynamic)) {
		ad = fit$a_dynamic; at = tr$a_ident
		if (all(dim(ad) == dim(at))) {
			row$cor_a_path = .cor_safe(as.vector(sweep(ad, 2, colMeans(ad))), as.vector(sweep(at, 2, colMeans(at))))
			# the per-period deviations alone (actor level removed), which is
			# the part the dynamic kernel adds
			row$cor_a_dev = .cor_safe(as.vector(ad - rowMeans(ad)), as.vector(at - rowMeans(at)))
		}
		# the ALS estimators echo the penalty persistence, not an estimate
		if (!is.null(fit$rho_ab) && !is_als) row$rho_ab_est = mean(fit$rho_ab, na.rm = TRUE)
		if (!is.null(fit$sigma_ab)) {
			row$sigma_ab_ratio_innov = mean(fit$sigma_ab, na.rm = TRUE) / tr$sigma_ab_innov
			row$sigma_ab_ratio_stat = mean(fit$sigma_ab, na.rm = TRUE) / tr$sd_dev
		}
	}
	if (isTRUE(cell$dynamic_uv)) {
		if (!is.null(fit$rho_uv) && !is_als) row$rho_uv_est = mean(fit$rho_uv, na.rm = TRUE)
		# symmetric fits fix sigma_uv at sqrt(1 - rho_uv^2) by the gauge, so the score does not apply there
		if (!is.null(fit$sigma_uv) && !isTRUE(cell$symmetric)) row$sigma_uv_ratio_innov = mean(fit$sigma_uv, na.rm = TRUE) / tr$sigma_uv_innov
	}

	# invariants (MCMC fits)
	if (!is_als) {
		note = character(0)
		safe = function(expr, label) tryCatch(
			suppressWarnings(suppressMessages(expr)),
			error = function(e) { note <<- c(note, paste0(label, ": ", substr(conditionMessage(e), 1, 80))); NA_real_ }
		)
		if (cell$estimator == "lame") {
			row$inv_pred = safe({
				p = predict(fit, newdata = sim$X, type = "link")
				max(sapply(seq_len(TT), function(t) max(abs(p[[t]][obs[[t]]] - EZ[[t]][obs[[t]]]))))
			}, "predict_newdata")
			row$inv_fit = safe(
				{ f = fitted(fit); max(
					sapply(seq_len(TT), function(t) max(abs(f[[t]][obs[[t]]] - YPM[[t]][obs[[t]]])))
				) }, "fitted")
			row$inv_res = safe(
				{ r = residuals(fit); max(
					sapply(
						seq_len(TT),
						function(t) max(abs(r[[t]][obs[[t]]] - (sim$Y[[t]][obs[[t]]] - YPM[[t]][obs[[t]]])))
					)
				) }, "residuals")
			row$fc_cor = safe({
				p = predict(fit, h = 1L, newdata = list(sim$X_fc), type = "link")
				o = !is.na(sim$Y_fc); .cor_safe(p[[1]][o], sim$eta_fc[o])
			}, "forecast_link")
			if (fam %in% c("normal", "binary", "poisson")) row$fc_ratio = safe({
				p = predict(fit, h = 1L, newdata = list(sim$X_fc), type = "response")
				o = !is.na(sim$Y_fc)
				if (fam == "normal") (mean(p[[1]][o]) - mean(sim$Y_fc[o])) / sd(sim$Y_fc[o]) else mean(p[[1]][o]) / mean(sim$Y_fc[o])
			}, "forecast_response")
			row$inv_sim = safe({
				s = simulate(fit, nsim = 100L, seed = 6886L)
				ys = mean(
					sapply(
						s$Y,
						function(traj) mean(unlist(lapply(seq_len(TT), function(t) traj[[t]][obs[[t]]])), na.rm = TRUE)
					)
				)
				yp = mean(unlist(lapply(seq_len(TT), function(t) YPM[[t]][obs[[t]]])))
				if (fam == "normal") ys - yp else ys / yp
			}, "simulate")
		} else {
			row$inv_pred = safe({
				p1 = predict(fit, type = "link"); p2 = predict(fit, type = "link", newdata = sim$X[[1]])
				max(abs(p1[obs[[1]]] - p2[obs[[1]]]))
			}, "predict_newdata")
			row$inv_fit = safe(max(abs(fitted(fit)[obs[[1]]] - YPM[[1]][obs[[1]]])), "fitted")
			row$inv_res = safe(
				max(abs(residuals(fit)[obs[[1]]] - (sim$Y[[1]][obs[[1]]] - YPM[[1]][obs[[1]]]))),
				"residuals"
			)
			row$inv_sim = safe({
				s = simulate(fit, nsim = 100L, seed = 6886L, newdata = list(Xdyad = sim$X[[1]]))
				ys = mean(sapply(s$Y, function(m) mean(m[obs[[1]]], na.rm = TRUE)))
				yp = mean(YPM[[1]][obs[[1]]])
				if (fam == "normal") ys - yp else ys / yp
			}, "simulate")
		}
		B = fit$BETA
		if (!is.null(B)) {
			row$inv_coef = safe({
				cf = coef(fit)
				if (length(dim(B)) == 3L) max(
					abs(as.matrix(cf) - apply(B, c(2, 3), mean))
				) else max(abs(cf - colMeans(B)))
			}, "coef")
			if (length(dim(B)) != 3L) {
				row$inv_ci = safe(
					{ ci = confint(fit); cf = coef(fit); as.numeric(
						all(cf >= ci[names(cf), 1] & cf <= ci[names(cf), 2])
					) }, "confint")
			}
		}
		row$inv_summary = safe({ invisible(capture.output(print(summary(fit)))); 1 }, "summary")
		if (tr$has_uv) row$inv_latent = safe(as.numeric(nrow(latent_positions(fit)) > 0), "latent_positions")
		te = fit$tryErrorChecks
		if (!is.null(te)) { row$tryerr = sum(unlist(te), na.rm = TRUE); row$tryerr_frac = row$tryerr / (cell$nscan + cell$burn) }
		row$terminated_early = as.numeric(isTRUE(fit$terminated_early))
		if (length(note)) row$notes = paste(c(row$notes, note), collapse = " | ")
	}
	if (is_als) {
		note = character(0)
		safe = function(expr, label) tryCatch(
			suppressWarnings(suppressMessages(expr)),
			error = function(e) { note <<- c(note, paste0(label, ": ", substr(conditionMessage(e), 1, 80))); NA_real_ }
		)
		row$inv_fit = safe(
			{ f = .as_list(fitted(fit, plug_in = TRUE), TT); max(
				sapply(seq_len(TT), function(t) max(abs(f[[t]][obs[[t]]] - YPM[[t]][obs[[t]]])))
			) }, "fitted")
		row$inv_res = safe({ r = .as_list(residuals(fit), TT); f = .as_list(fitted(fit), TT)
			max(
				sapply(seq_len(TT), function(t) max(abs(r[[t]][obs[[t]]] - (sim$Y[[t]][obs[[t]]] - f[[t]][obs[[t]]]))))
			) }, "residuals")
		row$inv_coef = safe({ cf = coef(fit); as.numeric(!all(is.finite(cf))) }, "coef")
		row$inv_summary = safe({ invisible(capture.output(print(summary(fit)))); 1 }, "summary")
		if (length(note)) row$notes = paste(c(row$notes, note), collapse = " | ")
	}
	if (length(fitted$warn_txt)) row$notes = paste(
		c(row$notes, paste0("warn: ", substr(paste(fitted$warn_txt, collapse = " // "), 1, 300))),
		collapse = " | "
	)
	if (length(fitted$msg_txt)) row$notes = paste(
		c(row$notes, paste0("msg: ", substr(paste(gsub("\\s+", " ", fitted$msg_txt), collapse = " // "), 1, 400))),
		collapse = " | "
	)
	as.data.frame(lapply(row, function(v) if (is.null(v)) NA else v), stringsAsFactors = FALSE)
}

# flags

grid_flags = function(df) {
	f = function(cond) ifelse(is.na(cond), FALSE, cond)
	flags = list(
		ERR = df$status == "error",
		BETA = f(abs(df$beta1_z) > 3),
		ICPT = f(abs(df$icpt_z) > 3),
		LEVEL = f(abs(df$level_gap) > 0.3),
		EZ = f(df$ez_cor < 0.6),
		YPM = f(df$ypm_ratio < 0.9 | df$ypm_ratio > 1.1) | f(abs(df$ypm_gap) > 0.1),
		VA = f(df$va_ratio < 0.5 | df$va_ratio > 2),
		VB = f(df$vb_ratio < 0.5 | df$vb_ratio > 2),
		VE = f(df$ve_ratio < 0.6 | df$ve_ratio > 1.6),
		RHO = f(abs(df$rho_gap) > 0.25),
		# an invariant that errored (its label appears in the notes) is a
		# finding in its own right, not a missing value
		INVERR = grepl(
			"(predict_newdata|fitted|residuals|simulate|forecast_link|forecast_response|summary|latent)\\s*:",
			df$notes %||% ""
		),
		CORA = f(df$cor_a < 0.6),
		CORB = f(df$cor_b < 0.6),
		UV = f(df$uv_cor < 0.5),
		DYNAB = f(abs(df$rho_ab_est - 0.7) > 0.3) | f(df$cor_a_path < 0.5),
		DYNUV = f(abs(df$rho_uv_est - 0.8) > 0.3),
		INV = f(df$inv_pred > 1e-3) | f(df$inv_fit > 1e-8) | f(df$inv_res > 1e-8) | f(df$inv_coef > 1e-8) |
			f(df$inv_ci == 0) | f(df$inv_summary != 1) | f(df$inv_latent == 0) |
			f(df$inv_sim < 0.85 | df$inv_sim > 1.15) & df$family != "normal" | f(
				abs(df$inv_sim) > 0.1 & df$family == "normal"
			),
		FC = f(df$fc_cor < 0.5) | f(df$fc_ratio < 0.7 | df$fc_ratio > 1.4) & df$family != "normal" | f(
			abs(df$fc_ratio) > 0.25 & df$family == "normal"
		),
		TRY = f(df$tryerr_frac > 0.05) | f(df$terminated_early == 1)
	)
	fm = as.data.frame(flags)
	df$flags = apply(fm, 1, function(r) paste(names(flags)[r], collapse = ","))
	df$n_flags = rowSums(fm)
	df
}

# flags that fire in every seed of a cell are systematic
grid_systematic = function(df) {
	df = grid_flags(df)
	sp = split(df, df$id)
	out = lapply(sp, function(d) {
		fl = lapply(strsplit(d$flags, ","), function(x) x[nzchar(x)])
		common = Reduce(intersect, fl)
		data.frame(
			id = d$id[1],
			group = d$group[1],
			family = d$family[1],
			estimator = d$estimator[1],
			seeds = nrow(d),
			systematic = paste(common, collapse = ","),
			any = paste(sort(unique(unlist(fl))), collapse = ","),
			stringsAsFactors = FALSE
		)
	})
	out = do.call(rbind, out)
	out[order(out$systematic == "", out$id), ]
}

# runner

grid_run = function(cells, seeds = 1:2, out_csv, verbose = TRUE) {
	done = if (file.exists(out_csv)) read.csv(out_csv, stringsAsFactors = FALSE) else NULL
	for (cell in cells) for (seed in seeds) {
		if (!is.null(done) && any(done$id == cell$id & done$seed == seed)) next
		if (verbose) cat(sprintf("[%s] %s seed %d ... ", format(Sys.time(), "%H:%M:%S"), cell$id, seed))
		sim = tryCatch(
			grid_sim(cell, seed),
			error = function(e) structure(list(message = conditionMessage(e)), class = "grid_error")
		)
		if (inherits(sim, "grid_error")) {
			row = data.frame(
				id = cell$id,
				group = cell$group,
				seed = seed,
				family = cell$family,
				estimator = cell$estimator,
				status = "sim_error",
				error = substr(sim$message, 1, 200),
				stringsAsFactors = FALSE
			)
			fitted = list(runtime = 0)
		} else {
			fitted = grid_fit(cell, sim, seed)
			row = tryCatch(grid_score(cell, sim, fitted, seed), error = function(e) {
				data.frame(
					id = cell$id,
					group = cell$group,
					seed = seed,
					family = cell$family,
					estimator = cell$estimator,
					status = "score_error",
					error = substr(conditionMessage(e), 1, 200),
					stringsAsFactors = FALSE
				)
			})
		}
		if (verbose) cat(
			sprintf(
				"%s (%.0fs) %s\n",
				row$status,
				fitted$runtime,
				if (!is.null(row$ez_cor)) sprintf(
					"ez_cor %.2f ypm %.2f va %.2f",
					row$ez_cor,
					row$ypm_ratio %||% NA,
					row$va_ratio %||% NA
				) else ""
			)
		)
		if (file.exists(out_csv)) {
			old = read.csv(out_csv, stringsAsFactors = FALSE)
			for (nm in setdiff(names(old), names(row))) row[[nm]] = NA
			for (nm in setdiff(names(row), names(old))) old[[nm]] = NA
			write.csv(rbind(old, row[names(old)]), out_csv, row.names = FALSE)
		} else {
			write.csv(row, out_csv, row.names = FALSE)
		}
	}
	invisible(read.csv(out_csv, stringsAsFactors = FALSE))
}

# refit one cell for diagnosis: returns sim, fitted, and the scored row
grid_repro = function(cell, seed) {
	sim = grid_sim(cell, seed)
	fitted = grid_fit(cell, sim, seed)
	list(sim = sim, fitted = fitted, row = grid_score(cell, sim, fitted, seed))
}

# cell lists per group

grid_cells = function(group) {
	cells = list()
	add = function(...) cells[[length(cells) + 1]] <<- grid_cell(..., group = group)
	fams = c("normal", "binary", "ordinal", "poisson", "cbin", "frn")
	if (group == "G1_static_uni") {
		for (est in c("ame", "lame")) for (fam in fams) {
			p = paste0(est, "_", fam)
			add(paste0(p, "_asym_R0"), fam, est)
			add(paste0(p, "_asym_R2"), fam, est, R = 2)
			add(paste0(p, "_asym_R0_dcor"), fam, est, dcor = TRUE)
			add(paste0(p, "_asym_R0_miss20"), fam, est, missing = 0.2)
			if (!fam %in% c("cbin", "frn")) {
				add(paste0(p, "_sym_R0"), fam, est, symmetric = TRUE)
				add(paste0(p, "_sym_R2"), fam, est, symmetric = TRUE, R = 2)
			}
		}
		for (est in c("ame", "lame")) {
			add(paste0(est, "_normal_asym_R0_noint"), "normal", est, intercept = FALSE)
			add(paste0(est, "_binary_asym_R0_noint"), "binary", est, intercept = FALSE)
			add(paste0(est, "_binary_asym_R0_sparse"), "binary", est, density = "sparse")
			add(paste0(est, "_binary_asym_R2_sparse"), "binary", est, R = 2, density = "sparse")
			add(paste0(est, "_binary_sym_R0_miss20"), "binary", est, symmetric = TRUE, missing = 0.2)
		}
	}
	if (group == "G2_bipartite") {
		for (est in c("ame", "lame")) for (fam in fams) {
			p = paste0(est, "_bip_", fam)
			add(paste0(p, "_R00"), fam, est, mode = "bipartite")
			add(paste0(p, "_R22"), fam, est, mode = "bipartite", R_row = 2, R_col = 2)
			add(paste0(p, "_R00_miss20"), fam, est, mode = "bipartite", missing = 0.2)
		}
		for (fam in c("normal", "binary", "poisson", "ordinal")) {
			p = paste0("lame_bip_", fam)
			add(paste0(p, "_dynab_R00"), fam, "lame", mode = "bipartite", dynamic_ab = TRUE, T = 6, nscan = 2500)
			add(
				paste0(p, "_dynuv_R22"),
				fam,
				"lame",
				mode = "bipartite",
				R_row = 2,
				R_col = 2,
				dynamic_uv = TRUE,
				T = 6,
				nscan = 2500
			)
			add(
				paste0(p, "_dynab_dynuv_R22"),
				fam,
				"lame",
				mode = "bipartite",
				R_row = 2,
				R_col = 2,
				dynamic_ab = TRUE,
				dynamic_uv = TRUE,
				T = 6,
				nscan = 2500
			)
		}
		add(
			"lame_bip_binary_dynab_R00_sparse", "binary", "lame", mode = "bipartite", dynamic_ab = TRUE, T = 6, density = "sparse", nscan = 2500
		)
	}
	if (group == "G3_dynamic_uni") {
		for (fam in c("normal", "binary", "poisson", "ordinal")) for (sym in c(FALSE, TRUE)) {
			p = paste0("lame_", fam, if (sym) "_sym" else "_asym")
			add(paste0(p, "_dynab_R0"), fam, "lame", symmetric = sym, dynamic_ab = TRUE, T = 6, nscan = 2500)
			add(paste0(p, "_dynuv_R2"), fam, "lame", symmetric = sym, R = 2, dynamic_uv = TRUE, T = 6, nscan = 2500)
			add(
				paste0(p, "_dynab_dynuv_R2"),
				fam,
				"lame",
				symmetric = sym,
				R = 2,
				dynamic_ab = TRUE,
				dynamic_uv = TRUE,
				T = 6,
				nscan = 2500
			)
			add(
				paste0(p, "_dynbeta_ar1_R0"),
				fam,
				"lame",
				symmetric = sym,
				dynamic_beta = TRUE,
				dynamic_beta_kind = "ar1",
				T = 6,
				nscan = 2500
			)
			add(
				paste0(p, "_dynbeta_rw1_R0"),
				fam,
				"lame",
				symmetric = sym,
				dynamic_beta = TRUE,
				dynamic_beta_kind = "rw1",
				T = 6,
				nscan = 2500
			)
			add(
				paste0(p, "_dynab_dynbeta_R0"),
				fam,
				"lame",
				symmetric = sym,
				dynamic_ab = TRUE,
				dynamic_beta = TRUE,
				T = 6,
				nscan = 2500
			)
		}
		for (fam in c("cbin", "frn")) {
			add(paste0("lame_", fam, "_asym_dynab_R0"), fam, "lame", dynamic_ab = TRUE, T = 6, nscan = 2500)
			add(paste0("lame_", fam, "_asym_dynuv_R2"), fam, "lame", R = 2, dynamic_uv = TRUE, T = 6, nscan = 2500)
		}
		add(
			"lame_binary_asym_dynab_R0_sparse", "binary", "lame", dynamic_ab = TRUE, T = 6, density = "sparse", nscan = 3000
		)
		add(
			"lame_binary_asym_dynab_dynuv_R2_sparse", "binary", "lame", R = 2, dynamic_ab = TRUE, dynamic_uv = TRUE, T = 6, density = "sparse", nscan = 3000
		)
		add("lame_binary_asym_dynab_R0_dcor", "binary", "lame", dynamic_ab = TRUE, dcor = TRUE, T = 6, nscan = 2500)
		add("lame_normal_asym_dynab_R0_miss20", "normal", "lame", dynamic_ab = TRUE, T = 6, missing = 0.2, nscan = 2500)
		add(
			"lame_binary_asym_dynuv_R2_miss20", "binary", "lame", R = 2, dynamic_uv = TRUE, T = 6, missing = 0.2, nscan = 2500
		)
	}
	if (group == "G4_als") {
		for (fam in c("normal", "binary", "poisson")) {
			for (est in c("ame_als", "lame_als")) {
				p = paste0(est, "_", fam)
				# binary ALS refuses below its size floor: every 30-actor
				# cross-section, and a four-period panel with a multiplicative term
				cs_ref = if (fam == "binary" && est == "ame_als") "Binary ALS is not supported" else NULL
				pn_ref = if (fam == "binary") "Binary ALS is not supported" else NULL
				add(paste0(p, "_asym_R0"), fam, est, expect_error = cs_ref)
				add(paste0(p, "_asym_R2"), fam, est, R = 2, expect_error = pn_ref)
				add(paste0(p, "_sym_R0"), fam, est, symmetric = TRUE, expect_error = cs_ref)
				add(paste0(p, "_sym_R2"), fam, est, symmetric = TRUE, R = 2, expect_error = pn_ref)
				add(paste0(p, "_bip_R00"), fam, est, mode = "bipartite", expect_error = cs_ref)
				add(paste0(p, "_bip_R22"), fam, est, mode = "bipartite", R_row = 2, R_col = 2, expect_error = pn_ref)
				add(paste0(p, "_asym_R0_miss20"), fam, est, missing = 0.2, expect_error = cs_ref)
			}
			add(paste0("lame_dynamic_als_", fam, "_asym_dynab_R0"), fam, "lame_dynamic_als", dynamic_ab = TRUE, T = 6)
			add(
				paste0("lame_dynamic_als_", fam, "_asym_dynuv_R2"),
				fam,
				"lame_dynamic_als",
				R = 2,
				dynamic_uv = TRUE,
				T = 6
			)
			add(
				paste0("lame_dynamic_als_", fam, "_bip_dynab_R00"),
				fam,
				"lame_dynamic_als",
				mode = "bipartite",
				dynamic_ab = TRUE,
				T = 6
			)
			# the same data through MCMC, for the cross-estimator comparison
			add(paste0("lame_", fam, "_asym_R2_alsref"), fam, "lame", R = 2, notes = "mcmc reference for the als cells")
			add(paste0("ame_", fam, "_asym_R2_alsref"), fam, "ame", R = 2, notes = "mcmc reference for the als cells")
		}
		add("ame_als_ordinal_asym_R0", "ordinal", "ame_als", expect_error = "ordinal|not supported|ALS")
		add("ame_als_cbin_asym_R0", "cbin", "ame_als", expect_error = "cbin|not supported|ALS")
	}
	if (group == "G5_kinds") {
		# the two other transition kinds of dynamic_uv, generated as the kind
		# describes them (heavy-tailed innovations; resets), directed and symmetric
		for (kind in c("snap", "t")) for (fam in c("normal", "binary", "poisson")) for (sym in c(FALSE, TRUE)) {
			p = paste0("lame_", fam, if (sym) "_sym" else "_asym")
			# symmetric fits take ar1 only: the kinds abort there
			ee = if (sym) "not available for symmetric" else NULL
			add(
				paste0(p, "_dynuv_", kind, "_R2"),
				fam,
				"lame",
				symmetric = sym,
				R = 2,
				dynamic_uv = TRUE,
				dynamic_uv_kind = kind,
				T = 6,
				nscan = 2500,
				expect_error = ee
			)
			if (fam != "poisson") {
				add(
					paste0(p, "_dynab_dynuv_", kind, "_R2"),
					fam,
					"lame",
					symmetric = sym,
					R = 2,
					dynamic_ab = TRUE,
					dynamic_uv = TRUE,
					dynamic_uv_kind = kind,
					T = 6,
					nscan = 2500,
					expect_error = ee
				)
			}
		}
		# the ar1 kind fitted to snap- and t-generated paths: the default's robustness
		for (kind in c("snap", "t")) {
			add(
				paste0("lame_binary_asym_dynuv_", kind, "data_ar1fit_R2"),
				"binary",
				"lame",
				R = 2,
				dynamic_uv = TRUE,
				T = 6,
				nscan = 2500,
				notes = paste0("data from the ", kind, " kind, fitted with ar1")
			)
			cells[[length(cells)]]$dynamic_uv_kind = "ar1"
			cells[[length(cells)]]$sim_kind = kind
		}
	}
	if (group == "G6_dynrho") {
		# a period-varying dyadic correlation (directed normal only)
		add("lame_normal_asym_dynrho_R0", "normal", "lame", dcor = TRUE, dynamic_rho = TRUE, T = 6, nscan = 2500)
		add("lame_normal_asym_dynrho_R2", "normal", "lame", R = 2, dcor = TRUE, dynamic_rho = TRUE, T = 6, nscan = 2500)
		add(
			"lame_normal_asym_dynab_dynrho_R0", "normal", "lame", dynamic_ab = TRUE, dcor = TRUE, dynamic_rho = TRUE, T = 6, nscan = 2500
		)
		add("lame_normal_asym_dynrho_R0_T10", "normal", "lame", dcor = TRUE, dynamic_rho = TRUE, T = 10, nscan = 2500)
		# a constant true rho fitted with dynamic_rho: the path stays flat
		add(
			"lame_normal_asym_dynrho_flat_R0", "normal", "lame", dcor = TRUE, dynamic_rho = TRUE, T = 6, nscan = 2500,
			notes = "constant generating rho 0.5 under a dynamic_rho fit"
		)
		cells[[length(cells)]]$flat_rho = TRUE
	}
	cells
}
