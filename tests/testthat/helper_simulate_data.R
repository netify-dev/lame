# helper function to generate simulated network data with known parameters
# for testing parameter recovery across families, modes, and temporal settings.
#
# this file is auto-loaded by testthat before any test files run.

simulate_test_network = function(
		n = 20,
		n_time = 1,
		family = "normal",
		R = 0,
		beta_intercept = 0,
		beta_dyad = 0.5,
		rvar = TRUE,
		cvar = TRUE,
		rho = 0,
		s2 = 1,
		sd_a = 0.5,
		sd_b = 0.5,
		rho_uv = NULL,
		rho_ab = NULL,
		sigma_uv = 0.3,
		sigma_ab = 0.3,
		mode = "unipartite",
		nA = NULL,
		nB = NULL,
		odmax = 3,
		symmetric = FALSE,
		seed = 42
) {
	set.seed(seed)

	bip = identical(mode, "bipartite")

	if (bip) {
		if (is.null(nA)) nA = n
		if (is.null(nB)) nB = max(n - 5, 8)
		nr = nA
		nc = nB
	} else {
		nr = n
		nc = n
	}

	# generate dyadic covariate (one continuous predictor)
	gen_Xdyad = function() {
		X = matrix(rnorm(nr * nc), nr, nc)
		if (!bip) diag(X) = NA
		X
	}

	# generate row and column covariates
	gen_Xrow = function() matrix(rnorm(nr), nr, 1)
	gen_Xcol = function() matrix(rnorm(nc), nc, 1)

	# generate additive effects
	gen_ab = function(a_prev = NULL, b_prev = NULL) {
		if (!is.null(rho_ab) && !is.null(a_prev)) {
			a = rho_ab * a_prev + rnorm(nr, 0, sigma_ab)
			b = rho_ab * b_prev + rnorm(nc, 0, sigma_ab)
		} else {
			a = if (rvar) rnorm(nr, 0, sd_a) else rep(0, nr)
			b = if (cvar) rnorm(nc, 0, sd_b) else rep(0, nc)
		}
		if (symmetric && !bip) b = a
		list(a = a, b = b)
	}

	# generate latent factors
	gen_UV = function(U_prev = NULL, V_prev = NULL) {
		if (R == 0) return(list(U = matrix(0, nr, 1), V = matrix(0, nc, 1)))
		if (!is.null(rho_uv) && !is.null(U_prev)) {
			U = rho_uv * U_prev + matrix(rnorm(nr * R, 0, sigma_uv), nr, R)
			V = rho_uv * V_prev + matrix(rnorm(nc * R, 0, sigma_uv), nc, R)
		} else {
			U = matrix(rnorm(nr * R, 0, 1), nr, R)
			V = matrix(rnorm(nc * R, 0, 1), nc, R)
		}
		if (symmetric && !bip) V = U
		list(U = U, V = V)
	}

	# compute Z = intercept + beta*X + outer(a,b) + UV' + noise
	compute_Z = function(Xd, a, b, U, V) {
		eta = beta_intercept + beta_dyad * Xd
		eta[is.na(Xd)] = beta_intercept
		eta = eta + outer(a, rep(1, nc)) + outer(rep(1, nr), b)
		if (R > 0) eta = eta + tcrossprod(U, V)
		if (symmetric && !bip) eta = (eta + t(eta)) / 2
		Z = eta + matrix(rnorm(nr * nc, 0, sqrt(s2)), nr, nc)
		if (!bip) diag(Z) = NA
		Z
	}

	# convert Z to Y via family-specific link
	Z_to_Y = function(Z, Xd) {
		Y = switch(family,
			"normal" = {
				Y = Z
				if (!bip) diag(Y) = NA
				Y
			},
			"binary" = {
				Y = ifelse(Z > 0, 1, 0)
				if (!bip) diag(Y) = NA
				Y
			},
			"poisson" = {
				# z is on log scale for poisson
				lambda = exp(Z)
				lambda[is.na(lambda)] = 1
				lambda[lambda > 1000] = 1000
				Y = matrix(rpois(length(lambda), lambda), nr, nc)
				if (!bip) diag(Y) = NA
				Y
			},
			"tobit" = {
				Y = Z
				Y[Y < 0] = 0
				if (!bip) diag(Y) = NA
				Y
			},
			"ordinal" = {
				# map continuous Z to 4 ordinal categories (0-3)
				breaks = quantile(Z, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
				Y = matrix(0L, nr, nc)
				Y[Z > breaks[1]] = 1L
				Y[Z > breaks[2]] = 2L
				Y[Z > breaks[3]] = 3L
				if (!bip) diag(Y) = NA
				Y
			},
			"cbin" = {
				# censored binary: binary with outdegree constraint
				Y = ifelse(Z > 0, 1, 0)
				if (!bip) diag(Y) = NA
				# enforce outdegree constraint
				for (i in 1:nr) {
					row_i = Y[i, ]
					row_i[is.na(row_i)] = 0
					if (sum(row_i) > odmax) {
						ones = which(row_i == 1)
						drop = sample(ones, sum(row_i) - odmax)
						Y[i, drop] = 0
					}
				}
				Y
			},
			"frn" = {
				# fixed rank nomination
				if (!bip) diag(Z) = -Inf
				Y = matrix(0, nr, nc)
				for (i in 1:nr) {
					zi = Z[i, ]
					zi[is.na(zi)] = -Inf
					top = order(zi, decreasing = TRUE)[1:min(odmax, sum(is.finite(zi)))]
					top = top[zi[top] > 0]
					if (length(top) > 0) Y[i, top] = seq_along(top)
				}
				if (!bip) diag(Y) = NA
				Y
			},
			"rrl" = {
				# row rank likelihood
				if (!bip) diag(Z) = -Inf
				Y = matrix(0, nr, nc)
				for (i in 1:nr) {
					zi = Z[i, ]
					valid = which(!is.na(zi) & is.finite(zi))
					n_nom = min(odmax, length(valid))
					top = valid[order(zi[valid], decreasing = TRUE)[1:n_nom]]
					Y[i, top] = seq(n_nom, 1)
				}
				if (!bip) diag(Y) = NA
				Y
			}
		)
		Y
	}

	# for poisson, adjust Z computation to use log link
	compute_Z_pois = function(Xd, a, b, U, V) {
		# for poisson, EZ is on the log scale
		eta = beta_intercept + beta_dyad * Xd
		eta[is.na(Xd)] = beta_intercept
		eta = eta + outer(a, rep(1, nc)) + outer(rep(1, nr), b)
		if (R > 0) eta = eta + tcrossprod(U, V)
		if (symmetric && !bip) eta = (eta + t(eta)) / 2
		if (!bip) diag(eta) = NA
		eta
	}

	# generate data for all time periods
	Y_list = list()
	Xdyad_list = list()
	Xrow_list = list()
	Xcol_list = list()
	true_a = list()
	true_b = list()
	true_U = list()
	true_V = list()

	a_prev = NULL
	b_prev = NULL
	U_prev = NULL
	V_prev = NULL

	for (t in seq_len(n_time)) {
		Xd = gen_Xdyad()
		Xr = gen_Xrow()
		Xc = gen_Xcol()

		ab = gen_ab(a_prev, b_prev)
		uv = gen_UV(U_prev, V_prev)

		if (family == "poisson") {
			Z = compute_Z_pois(Xd, ab$a, ab$b, uv$U, uv$V)
			Y = Z_to_Y(Z, Xd)
		} else {
			Z = compute_Z(Xd, ab$a, ab$b, uv$U, uv$V)
			Y = Z_to_Y(Z, Xd)
		}

		Y_list[[t]] = Y

		# format covariates as 3D arrays (n x n x 1) for Xdyad
		Xdyad_list[[t]] = array(Xd, dim = c(nr, nc, 1))
		Xrow_list[[t]] = Xr
		Xcol_list[[t]] = Xc

		true_a[[t]] = ab$a
		true_b[[t]] = ab$b
		true_U[[t]] = uv$U
		true_V[[t]] = uv$V

		a_prev = ab$a
		b_prev = ab$b
		U_prev = uv$U
		V_prev = uv$V
	}

	# for cross-sectional (n_time=1), unwrap from lists
	if (n_time == 1) {
		Y_out = Y_list[[1]]
		Xdyad_out = Xdyad_list[[1]]
		Xrow_out = Xrow_list[[1]]
		Xcol_out = Xcol_list[[1]]
	} else {
		# for longitudinal, lame() expects lists as input
		# (it will internally convert to arrays via list_to_array)
		Y_out = Y_list
		names(Y_out) = paste0("t", seq_len(n_time))

		# Xdyad: list of T 3D arrays (nr x nc x 1)
		Xdyad_out = Xdyad_list
		names(Xdyad_out) = paste0("t", seq_len(n_time))

		# Xrow: list of T matrices (nr x 1)
		Xrow_out = Xrow_list
		names(Xrow_out) = paste0("t", seq_len(n_time))

		# Xcol: list of T matrices (nc x 1)
		Xcol_out = Xcol_list
		names(Xcol_out) = paste0("t", seq_len(n_time))
	}

	# name actors consistently
	row_actors = paste0("actor", 1:nr)
	col_actors = if (bip) paste0("col", 1:nc) else row_actors

	if (n_time == 1) {
		rownames(Y_out) = row_actors
		colnames(Y_out) = col_actors
		rownames(Xdyad_out) = row_actors
		colnames(Xdyad_out) = col_actors
		rownames(Xrow_out) = row_actors
		rownames(Xcol_out) = col_actors
	} else {
		for (t in seq_len(n_time)) {
			rownames(Y_out[[t]]) = row_actors
			colnames(Y_out[[t]]) = col_actors
			rownames(Xdyad_out[[t]]) = row_actors
			colnames(Xdyad_out[[t]]) = col_actors
			rownames(Xrow_out[[t]]) = row_actors
			rownames(Xcol_out[[t]]) = col_actors
		}
	}

	list(
		Y = Y_out,
		Xdyad = Xdyad_out,
		Xrow = Xrow_out,
		Xcol = Xcol_out,
		true_params = list(
			beta_intercept = beta_intercept,
			beta_dyad = beta_dyad,
			a = true_a,
			b = true_b,
			U = true_U,
			V = true_V,
			rho = rho,
			s2 = s2,
			rho_uv = rho_uv,
			rho_ab = rho_ab,
			sigma_uv = sigma_uv,
			sigma_ab = sigma_ab
		),
		family = family,
		mode = mode,
		n_time = n_time
	)
}
