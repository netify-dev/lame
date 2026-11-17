skip_on_cran()

library(lame)

# symmetric panel u_i' L u_j whose paths reset with probability pi_true
snap_sym_panel = function(family = "binary", n = 9L, N = 4L, seed = 1L, pi_true = 0.15) {
	set.seed(seed)
	X = array(0, c(n, 1L, N)); D = matrix(0L, n, N)
	X[, 1, 1] = rnorm(n)
	for (t in 2:N) {
		D[, t] = as.integer(runif(n) < pi_true)
		X[, 1, t] = ifelse(D[, t] == 1L, rnorm(n), 0.9 * X[, 1, t - 1] + rnorm(n, 0, sqrt(0.19)))
	}
	nm = sprintf("a%02d", 1:n)
	Y = lapply(seq_len(N), function(t) {
		eta = -0.3 + 0.8 * tcrossprod(X[, 1, t])
		e = matrix(rnorm(n * n), n, n)
		e[lower.tri(e)] = t(e)[lower.tri(e)]
		y = switch(
			family,
			binary = 1 * (eta + e > 0),
			normal = eta + 0.5 * e,
			poisson = matrix(rpois(n * n, exp(eta)), n, n)
		)
		y[lower.tri(y)] = t(y)[lower.tri(y)]
		diag(y) = NA
		dimnames(y) = list(nm, nm)
		y
	})
	names(Y) = paste0("t", seq_len(N))
	list(Y = Y, D = D)
}

# rank AUC of reset probabilities against the planted resets
snap_auc = function(score, truth) {
	r = rank(score)
	n1 = sum(truth == 1); n0 = sum(truth == 0)
	(sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# paths (n x R x nt) and resets (n x nt) drawn from the snap prior; rho, sg, kap, pis scalars
snap_prior_paths = function(n, R, nt, rho, sg, kap, pis) {
	X = array(0, c(n, R, nt)); D = matrix(0, n, nt)
	X[, , 1] = rnorm(n * R, 0, sg / sqrt(1 - rho^2))
	for (t in 2:nt) {
		D[, t] = 1 * (runif(n) < pis)
		X[, , t] = ifelse(matrix(D[, t] == 1, n, R), rnorm(n * R, 0, kap), rho * X[, , t - 1] + rnorm(n * R, 0, sg))
	}
	list(X = X, D = D)
}

# z of the chains' means against the prior means
snap_rep_z = function(out, target) (colMeans(out) - target) / (apply(out, 2, sd) / sqrt(nrow(out)))

test_that("the symmetric snap path visit keeps the snap prior of the paths when the data are drawn from them", {
	# successive conditionals from exact prior draws of the unit-variance paths of U_t L U_t': probit, Poisson and
	# Gaussian pairs given the paths, then one visit, which moves both triangles and the diagonal
	n = 4L; R = 2L; nt = 4L; rho = .8; sg = sqrt(1 - rho^2); kap = 1.5; pis = .2; s2 = .5
	Ld = c(1.3, -.8)
	vt = numeric(nt); vt[1] = 1
	for (t in 2:nt) vt[t] = pis * kap^2 + (1 - pis) * (rho^2 * vt[t - 1] + sg^2)
	ez = function(U) {
		E = array(0, c(n, n, nt))
		for (t in seq_len(nt)) E[, , t] = -.3 + U[, , t] %*% diag(Ld) %*% t(U[, , t])
		E
	}
	pair_sym = function(A) {
		for (t in seq_len(nt)) { a = A[, , t]; a[lower.tri(a)] = t(a)[lower.tri(a)]; A[, , t] = a }
		A
	}
	for (kind in 0:2) {
		set.seed(171 + kind)
		K = 500
		out = matrix(0, K, 1 + nt)
		ez_err = 0
		for (k in seq_len(K)) {
			pu = snap_prior_paths(n, R, nt, rho, sg, kap, pis)
			U = pu$X; du = pu$D
			for (it in 1:10) {
				EZ = ez(U)
				cnt = numeric(0)
				mask = array(1, dim(EZ))
				Z = EZ
				if (kind == 0L) {
					mask = pair_sym(2 * array(rbinom(length(EZ), 1, pnorm(EZ)), dim(EZ)) - 1)
				} else if (kind == 1L) {
					cnt = pair_sym(array(rpois(length(EZ), exp(EZ)), dim(EZ)))
				} else {
					Z = EZ + pair_sym(array(rnorm(length(EZ), 0, sqrt(s2)), dim(EZ)))
				}
				for (t in seq_len(nt)) diag(mask[, , t]) = 0
				o = snap_path_visit_mode_cpp(EZ, Z, numeric(0), mask, cnt, numeric(nt), U, U, matrix(0, 0, 0), Ld,
					rho, sg, kap, pis, du, du, 1, Inf, s2, numeric(nt), kind, 2L)
				ez_err = max(ez_err, abs(o$EZ - ez(o$U)))
				if (kind == 1L) ez_err = max(ez_err, abs(o$Z - ez(o$U)))
				U = o$U; du = o$delta_u
			}
			out[k, ] = c(mean(du[, -1]), apply(U^2, 3, mean))
		}
		expect_lt(ez_err, 1e-10)
		expect_equal(du[, 1], rep(0, n))
		expect_lt(max(abs(snap_rep_z(out, c(pis, vt)))), 4, label = paste("kind", kind))
	}
	# the R visit returns V_t = U_t diag(L) and the resets on both sides
	EZ = ez(U)
	o = .lame_dpm_snap_visit(NULL, EZ, EZ, 0, U, U, rho, sg, list(pi = pis, kappa = kap, delta_u = du, delta_v = du),
		family = "gaussian", s2 = s2, mode = "sym", Ld = Ld)
	expect_equal(o$V[, , 2], o$U[, , 2] %*% diag(Ld), tolerance = 1e-12)
	expect_identical(o$delta_v, o$delta_u)
})

test_that("the symmetric snap sweep, its gauge moves and its hyperparameters keep the joint prior", {
	# successive conditionals from exact prior draws of (paths, resets, L, rho, pi) in the unit-variance gauge
	# (sigma^2 = 1 - rho^2, linked reset scale 1): Gaussian pairs given the state, then the Gibbs sweep, the path
	# move, the draw of L, the split of each column against its weight, the period signs, the rotations and boosts
	# of the coordinate pairs, the persistence over the drifting transitions and its non-centred draw with the resets
	# held, and the reset rate
	n = 4L; R = 2L; nt = 4L; s2 = .6
	a = (-.99 - .9) / .1; b = (.99 - .9) / .1
	rho_mean = .9 + .1 * (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
	set.seed(181)
	K = 300
	out = matrix(0, K, 5)
	mask = array(1, c(n, n, nt)); for (t in seq_len(nt)) diag(mask[, , t]) = 0
	for (k in seq_len(K)) {
		rho = qnorm(runif(1, pnorm(a), pnorm(b))) * .1 + .9
		pis = rbeta(1, 1, 9); Ld = rnorm(R)
		pu = snap_prior_paths(n, R, nt, rho, .lame_sym_uv_innov_sd(rho), 1, pis)
		U = pu$X; du = pu$D
		for (it in 1:8) {
			sg = .lame_sym_uv_innov_sd(rho)
			kap = .lame_snap_kappa(NULL, rho, sg)
			EZ = array(0, c(n, n, nt))
			for (t in seq_len(nt)) EZ[, , t] = U[, , t] %*% diag(Ld) %*% t(U[, , t])
			E = EZ + array(rnorm(length(EZ), 0, sqrt(s2)), dim(EZ))
			for (t in seq_len(nt)) { e = E[, , t]; e[lower.tri(e)] = t(e)[lower.tri(e)]; diag(e) = 0; E[, , t] = e }
			o = rUV_dynamic_snap_sym_fc_cpp(U, Ld, E, rho, sg, s2, kap, pis, du)
			U = o$U; du = o$delta_u
			for (t in seq_len(nt)) EZ[, , t] = U[, , t] %*% diag(Ld) %*% t(U[, , t])
			o = snap_path_visit_mode_cpp(EZ, E, numeric(0), mask, numeric(0), numeric(0), U, U, matrix(0, 0, 0), Ld,
				rho, sg, kap, pis, du, du, 1, Inf, s2, numeric(nt), 2L, 2L)
			U = o$U; du = o$delta_u
			Ld = .lame_sym_dyn_L_draw(E, U, Ld, s2, 1)
			kw = .lame_uv_dyn_kind_weights("snap", NULL, NULL, du, du, sg^2 / kap^2)
			o = .lame_uv_sym_split_dynamic_draw(U, Ld, rho, sg, 1, w = kw$w_u, k = kw$k_u)
			U = o$U; Ld = o$Ld
			U = .lame_uv_sym_sign_draw(U, rho, keep = 1 - du)$U
			U = .lame_snap_sym_orient_draw(U, Ld, rho, sg, du, sg^2 / kap^2)$U
			rho = .lame_dpm_rho_unit_mh(U, rho, .9, .1, keep = 1 - du)
			o = .lame_uv_sym_nc_rho(U, Ld, E, s2, rho, .9, .1, keep = 1 - du)
			U = o$U; rho = o$rho
			pis = .lame_snap_pi_draw(du, matrix(0, 0, nt), 1, 9)
		}
		out[k, ] = c(mean(U^2), mean(Ld^2), rho, pis, mean(du[, -1]))
	}
	expect_lt(max(abs(snap_rep_z(out, c(1, 1, rho_mean, .1, .1)))), 4)
	# without resets the persistence draws are the ar1 ones
	set.seed(5); r1 = .lame_dpm_rho_unit_mh(U, .8, .9, .1, with_u = TRUE)
	set.seed(5); r2 = .lame_dpm_rho_unit_mh(U, .8, .9, .1, with_u = TRUE, keep = matrix(1, n, nt))
	expect_equal(r1, r2)
	set.seed(6); n1 = .lame_uv_sym_nc_rho(U, Ld, E, s2, .8, .9, .1)
	set.seed(6); n2 = .lame_uv_sym_nc_rho(U, Ld, E, s2, .8, .9, .1, keep = matrix(1, n, nt))
	expect_equal(n1, n2)
})

test_that("the persistence draw of symmetric snap paths that carries the latent array of a count panel keeps the joint law", {
	# draws of (rho, resets, unit-variance paths, latent array, counts) from the joint law, then one non-centred
	# persistence draw with the reset positions held and the latent array carried (R/uv_nc_step.R); every statistic
	# keeps its mean
	n = 6L; R = 2L; N = 4L; s2 = .4; L = c(.8, -.5); off = c(0, .3, -.2, .1); pis = .25
	rtn = function() { repeat { r = rnorm(1, .9, .1); if (abs(r) < .99) return(r) } }
	prod_of = function(U) {
		P = array(0, c(n, n, N))
		for (t in 1:N) { Ut = matrix(U[, , t], n, R); P[, , t] = Ut %*% diag(L) %*% t(Ut) }
		P
	}
	run = function(drift = FALSE, ndraw = 4000L, seed = 43L) {
		set.seed(seed)
		S = matrix(0, ndraw, 11)
		for (r in seq_len(ndraw)) {
			rho = rtn()
			keep = matrix(1, n, N); keep[, -1] = 1 * (runif(n * (N - 1)) > pis)
			U = .lame_uv_unit_paths(array(rnorm(n * R * N), c(n, R, N)), rho, keep)
			P = prod_of(U); Z = Y = array(0, c(n, n, N))
			for (t in 1:N) {
				z = P[, , t] + matrix(rnorm(n * n, 0, sqrt(s2)), n)
				z[lower.tri(z)] = t(z)[lower.tri(z)]
				Z[, , t] = z
				Y[, , t] = matrix(rpois(n * n, exp(z + off[t])), n)
			}
			stat = function(rho, U, Z) {
				e = Z - prod_of(U)
				c(rho, rho^2, sum(U^2), sum(U[, , 1] * U[, , 2]), sum(U[, , 3] * U[, , 4]), e[1, 2, 1], e[3, 4, 3]^2,
					sum(Z[, , 4]), sum(Y[, , 1] * Z[, , 1]), sum(Y[, , 3] * Z[, , 3]), e[2, 2, 2]^2)
			}
			yy = Y
			for (t in 1:N) yy[, , t][lower.tri(yy[, , t], TRUE)] = NA
			cx = .lame_ltm_context(yy, "poisson", "sym", offset = off)
			o = .lame_uv_sym_nc_rho_carry(U, L, Z, cx, rho, .9, .1, keep = if (!drift) keep)
			S[r, ] = stat(o$rho, o$U, o$Z) - stat(rho, U, Z)
		}
		sdv = apply(S, 2L, sd)
		max(abs(colMeans(S)[sdv > 0] / (sdv[sdv > 0] / sqrt(ndraw))))
	}
	expect_lt(run(), 4)
	# the check fails with the reset positions drawn as drifts
	expect_gt(run(TRUE, 4000L, 44L), 5)
})

test_that("the symmetric snap path visit draws a Poisson pair from its posterior", {
	# two actors, rank one, three periods, counts on their pair: the resets and the product of the two paths
	# against importance sampling from the snap prior
	nt = 3L; rho = .8; sg = sqrt(1 - rho^2); kap = 1; pis = .2; l1 = 1.4
	y = c(2, 0, 3)
	prior_draw = function(M) {
		x = matrix(0, M, nt); d = matrix(0, M, nt)
		x[, 1] = rnorm(M)
		for (t in 2:nt) {
			d[, t] = runif(M) < pis
			x[, t] = ifelse(d[, t] == 1, rnorm(M, 0, kap), rho * x[, t - 1] + rnorm(M, 0, sg))
		}
		list(x = x, d = d)
	}
	set.seed(191)
	M = 1e6
	p1 = prior_draw(M); p2 = prior_draw(M)
	lw = rowSums(dpois(matrix(y, M, nt, byrow = TRUE), exp(l1 * p1$x * p2$x), log = TRUE))
	w = exp(lw - max(lw)); w = w / sum(w)
	ex = c(colSums(w * p1$d)[-1], colSums(w * p2$d)[-1], colSums(w * l1 * p1$x * p2$x))
	rm(p1, p2)
	U = array(0, c(2, 1, nt)); du = matrix(0, 2, nt)
	mask = array(0, c(2, 2, nt)); mask[1, 2, ] = mask[2, 1, ] = 1
	cnt = array(0, c(2, 2, nt)); cnt[1, 2, ] = cnt[2, 1, ] = y
	n_it = 100000
	acc = numeric(length(ex))
	for (it in seq_len(n_it)) {
		EZ = array(0, c(2, 2, nt))
		for (t in seq_len(nt)) EZ[, , t] = l1 * tcrossprod(U[, , t])
		o = snap_path_visit_mode_cpp(EZ, EZ, numeric(0), mask, cnt, numeric(nt), U, U, matrix(0, 0, 0), l1,
			rho, sg, kap, pis, du, du, 1, Inf, 1, numeric(nt), 1L, 2L)
		U = o$U; du = o$delta_u
		acc = acc + c(du[1, -1], du[2, -1], l1 * U[1, 1, ] * U[2, 1, ])
	}
	est = acc / n_it
	expect_lt(max(abs(est[1:4] - ex[1:4])), .012)
	expect_lt(max(abs(est[5:7] - ex[5:7])), .02)
})

test_that("symmetric snap fits run in every family the kind serves, with one path move per actor", {
	for (fam in c("binary", "normal", "poisson")) {
		d = snap_sym_panel(fam, seed = 3)
		f = suppressWarnings(suppressMessages(lame(
			d$Y,
			family = fam,
			symmetric = TRUE,
			R = 1,
			dynamic_uv = TRUE,
			dynamic_uv_kind = "snap",
			burn = 10,
			nscan = 30,
			odens = 5,
			seed = 3,
			chains = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)))
		expect_true(all(is.finite(f$pi_snap)), info = fam)
		expect_equal(dim(f$snap_prob), c(9L, 4L), info = fam)
		expect_null(f$snap_prob_v)
		expect_true(all(is.na(f$snap_prob[, 1])), info = fam)
		expect_equal(unname(f$actor_move["dpm_calls"]), 40 * 9 / if (fam == "poisson") 2 else 1, info = fam)
		expect_true(all(abs(f$sigma_uv - sqrt(1 - f$rho_uv^2)) < 1e-12), info = fam)
	}
})

test_that("two chains of a symmetric snap fit return a gated fit that recovers the planted resets", {
	d = snap_sym_panel("normal", n = 16L, N = 5L, seed = 11L, pi_true = .15)
	f = suppressWarnings(suppressMessages(lame(
		d$Y,
		family = "normal",
		symmetric = TRUE,
		R = 1,
		dynamic_uv = TRUE,
		dynamic_uv_kind = "snap",
		burn = 300,
		nscan = 1200,
		odens = 10,
		seed = 6886,
		chains = 2,
		cores = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_s3_class(f, "lame")
	expect_equal(f$chains, 2L)
	expect_true(all(c("rho_uv", "pi_snap") %in% f$rhat$parameter))
	score = f$snap_prob[, -1]
	truth = d$D[, -1]
	expect_gt(snap_auc(score, truth), .65)
})

test_that("reporting, forecasts, simulation and resumption work on a symmetric snap fit", {
	d = snap_sym_panel("poisson", seed = 5)
	ck = tempfile(fileext = ".rds")
	f = suppressWarnings(suppressMessages(lame(
		d$Y,
		family = "poisson",
		symmetric = TRUE,
		R = 1,
		dynamic_uv = TRUE,
		dynamic_uv_kind = "snap",
		keep_snap_draws = "draws",
		burn = 20,
		nscan = 60,
		odens = 4,
		seed = 2,
		chains = 1,
		checkpoint_path = ck,
		checkpoint_every = 20L,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(summary(f)$pi_snap, mean(f$pi_snap))
	expect_match(paste(capture.output(print(summary(f))), collapse = " "), "pi_snap =")
	expect_no_error(capture.output(print(f)))
	expect_no_error(latent_positions(f))
	expect_no_error(trace_plot(f, params = "dynamic"))
	fc = predict(f, h = 2)
	expect_length(fc, 2L)
	expect_true(all(is.finite(unlist(fc))))
	expect_true(isSymmetric(unname(fc[[1]])))
	expect_no_error(simulate(f, nsim = 2, seed = 1))
	expect_equal(dim(f$snap_draws), c(15L, 9L, 4L))
	expect_null(f$snap_draws_v)
	expect_equal(nrow(snap_index_draws(f, side = "u")), 15L * 4L)
	st = f$start_vals$checkpoint_state$values
	expect_equal(dim(st$delta_uv), c(9L, 4L))
	g = suppressWarnings(suppressMessages(lame_resume(ck, nscan_more = 20L)))
	expect_equal(length(g$pi_snap), 5L)
	unlink(ck)
})
