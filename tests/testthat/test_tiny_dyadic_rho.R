tiny_uv_context = function(rank = 1L, reciprocal = 1e13) {
	U = array(seq_len(3L * rank * 2L) / 40, c(3L, rank, 2L))
	V = array(rev(seq_len(3L * rank * 2L)) / 30, c(3L, rank, 2L))
	E = array(0.1, c(3L, 3L, 2L))
	E[2:3, 1L, 1L] = reciprocal * c(1, -0.3)
	list(U = U, V = V, E = E, ar = 0.35, sigma = 0.8, s2 = 1.3)
}

tiny_uv_oracle = function(cx, rho, route, snap = 0) {
	rank = dim(cx$U)[2L]
	H = diag((1 - cx$ar^2) / cx$sigma^2, rank)
	h = numeric(rank)
	forward = if (route == "snap") as.numeric(snap < 0.5) else if (route == "t") 1.7 else 1
	H = H + diag(forward * cx$ar^2 / cx$sigma^2, rank)
	h = h + forward * cx$ar * cx$U[1L, , 2L] / cx$sigma^2
	pair_cov = cx$s2 * matrix(c(1, rho, rho, 1), 2L)
	for (j in 2:3) {
		design = rbind(cx$V[j, , 1L], numeric(rank))
		response = c(cx$E[1L, j, 1L], cx$E[j, 1L, 1L] - sum(cx$U[j, , 1L] * cx$V[1L, , 1L]))
		H = H + crossprod(design, solve(pair_cov, design))
		h = h + drop(crossprod(design, solve(pair_cov, response)))
	}
	v = cx$V[1L, , 1L]
	H = H + tcrossprod(v) / (cx$s2 * (1 + rho))
	h = h + v * cx$E[1L, 1L, 1L] / (cx$s2 * (1 + rho))
	covariance = solve(H)
	list(mean = drop(solve(H, h)), covariance = covariance)
}

tiny_uv_call = function(cx, rho, route, prefix = "", snap = 0) {
	suffix = switch(route, gaussian = "", snap = "_snap", t = "_t")
	fn = get(paste0(prefix, "rUV_dynamic", suffix, "_fc_cpp"), mode = "function")
	common = list(
		U_current = cx$U,
		V_current = cx$V,
		ET = cx$E,
		rho_uv = cx$ar,
		sigma_uv = cx$sigma,
		s2 = cx$s2
	)
	if (route == "snap") {
		common = c(
			common,
			list(
				kappa = 2,
				pi_snap = 0.2,
				delta_u_current = matrix(snap, 3L, 2L),
				delta_v_current = matrix(snap, 3L, 2L)
			)
		)
	}
	if (route == "t") {
		common = c(
			common,
			list(
				nu = 5,
				lambda_u = matrix(1.7, 3L, 2L),
				lambda_v = matrix(1.7, 3L, 2L)
			)
		)
	}
	do.call(fn, c(common, list(shrink = FALSE, symmetric = FALSE, rho = rho)))
}

test_that("nonzero tiny rho enters the exact missing-cell conditional", {
	Y = array(0, c(2L, 2L, 1L))
	Y[1L, 2L, 1L] = NA_real_
	Y[2L, 1L, 1L] = 1e11
	zero = array(0, dim(Y))
	for (rho in c(-1e-11, 0, 1e-11, -0.4, 0.4)) {
		set.seed(47332971)
		expected = rho * 1e11 + sqrt((1 - rho) * (1 + rho)) * sqrt(1.3) * rnorm(1L)
		rng = .Random.seed
		set.seed(47332971)
		actual = rZ_nrm_batch_cpp(zero, zero, rho, 1.3, Y)
		expect_equal(actual$Z[1L, 2L, 1L], expected, tolerance = 1e-12)
		expect_identical(.Random.seed, rng)
	}
})

test_that("all directed dynamic factor routes retain tiny rho in the first conditional", {
	cx = tiny_uv_context()
	for (route in c("gaussian", "snap", "t")) {
		for (rho in c(-1e-13, 0, 1e-13)) {
			law = tiny_uv_oracle(cx, rho, route)
			reference = cx
			reference$E[1L, , 1L] = 0
			reference$E[2:3, 1L, 1L] = drop(reference$U[2:3, , 1L] * reference$V[1L, , 1L])
			reference$U[1L, , 2L] = 0
			expect_equal(tiny_uv_oracle(reference, rho, route)$mean, 0, tolerance = 1e-14)
			set.seed(47332972)
			noise = abs(tiny_uv_call(reference, rho, route)$U[1L, 1L, 1L])
			set.seed(47332972)
			actual = abs(tiny_uv_call(cx, rho, route)$U[1L, 1L, 1L])
			# a one-dimensional reporting rotation changes only the sign
			# the paired zero-mean law has identical covariance and native rng
			expected = abs(law$mean + c(-noise, noise))
			expect_lt(min(abs(actual - expected)), 1e-10)
		}
	}
})

test_that("invalid directed rho fails before a random draw", {
	cx = tiny_uv_context(reciprocal = 1)
	zero = array(0, c(2L, 2L, 1L))
	for (rho in c(-1, 1, -1.1, 1.1, Inf, NaN)) {
		for (route in c("gaussian", "snap", "t")) {
			set.seed(47332973)
			rng = .Random.seed
			expect_error(tiny_uv_call(cx, rho, route), "strictly between")
			expect_identical(.Random.seed, rng)
		}
		set.seed(47332973)
		rng = .Random.seed
		expect_error(rZ_nrm_batch_cpp(zero, zero, rho, 1, zero), "strictly between")
		expect_identical(.Random.seed, rng)
	}
})
