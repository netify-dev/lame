skip_on_cran()

library(lame)

# bipartite panel whose row and column paths reset with probability pi_true
snap_bip_panel = function(family = "binary", nA = 8L, nB = 10L, N = 4L, seed = 1L, pi_true = 0.15) {
	set.seed(seed)
	paths = function(n) {
		X = array(0, c(n, 1L, N)); D = matrix(0L, n, N)
		X[, 1, 1] = rnorm(n, 0, 0.8)
		for (t in 2:N) {
			D[, t] = as.integer(runif(n) < pi_true)
			X[, 1, t] = ifelse(D[, t] == 1L, rnorm(n, 0, 0.8), 0.9 * X[, 1, t - 1] + rnorm(n, 0, 0.8 * sqrt(0.19)))
		}
		list(X = X, D = D)
	}
	pu = paths(nA); pv = paths(nB)
	Y = lapply(seq_len(N), function(t) {
		eta = -0.3 + tcrossprod(pu$X[, 1, t], pv$X[, 1, t])
		y = switch(
			family,
			binary = 1 * (eta + matrix(rnorm(nA * nB), nA, nB) > 0),
			normal = eta + matrix(rnorm(nA * nB, 0, 0.5), nA, nB),
			poisson = matrix(rpois(nA * nB, exp(eta)), nA, nB)
		)
		dimnames(y) = list(sprintf("r%02d", 1:nA), sprintf("c%02d", 1:nB))
		y
	})
	names(Y) = paste0("t", seq_len(N))
	list(Y = Y, D_u = pu$D, D_v = pv$D)
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

test_that("the bipartite snap path visit keeps the snap prior of the paths when the data are drawn from them", {
	# successive conditionals from exact prior draws: independent short chains started from the snap prior, each step
	# data given the paths and then one visit, keep the prior at every step. probit, Poisson and Gaussian cells,
	# rectangular ranks through G
	nA = 4L; nB = 3L; RA = 2L; RB = 1L; nt = 4L; rho = .8; sg = .4; kap = 2; pis = .2; s2 = .5
	G = matrix(c(1.2, -.7), RA, RB)
	vt = numeric(nt); vt[1] = sg^2 / (1 - rho^2)
	for (t in 2:nt) vt[t] = pis * kap^2 + (1 - pis) * (rho^2 * vt[t - 1] + sg^2)
	ez = function(U, V) {
		E = array(0, c(nA, nB, nt))
		for (t in seq_len(nt)) E[, , t] = -.3 + matrix(U[, , t], nA) %*% G %*% t(matrix(V[, , t], nB))
		E
	}
	for (kind in 0:2) {
		set.seed(71 + kind)
		K = 500
		out = matrix(0, K, 1 + 2 * nt)
		ez_err = 0
		for (k in seq_len(K)) {
			pu = snap_prior_paths(nA, RA, nt, rho, sg, kap, pis); pv = snap_prior_paths(nB, RB, nt, rho, sg, kap, pis)
			U = pu$X; V = pv$X; du = pu$D; dv = pv$D
			for (it in 1:10) {
				EZ = ez(U, V)
				cnt = numeric(0)
				mask = array(1, dim(EZ))
				Z = EZ
				if (kind == 0L) {
					mask = 2 * array(rbinom(length(EZ), 1, pnorm(EZ)), dim(EZ)) - 1
				} else if (kind == 1L) {
					cnt = array(rpois(length(EZ), exp(EZ)), dim(EZ))
				} else {
					Z = EZ + array(rnorm(length(EZ), 0, sqrt(s2)), dim(EZ))
				}
				o = snap_path_visit_mode_cpp(EZ, Z, numeric(0), mask, cnt, numeric(nt), U, V, G, numeric(0), rho, sg, kap,
					pis, du, dv, 1, Inf, s2, numeric(nt), kind, 1L)
				ez_err = max(ez_err, abs(o$EZ - ez(o$U, o$V)))
				U = o$U; V = o$V; du = o$delta_u; dv = o$delta_v
			}
			out[k, ] = c(mean(c(du[, -1], dv[, -1])), apply(U^2, 3, mean), apply(V^2, 3, mean))
		}
		expect_lt(ez_err, 1e-10)
		expect_equal(c(du[, 1], dv[, 1]), rep(0, nA + nB))
		expect_lt(max(abs(snap_rep_z(out, c(pis, vt, vt)))), 4, label = paste("kind", kind))
	}
})

test_that("the bipartite snap sweep, its gauge and its hyperparameters keep the joint prior", {
	# successive conditionals from exact prior draws of (paths, resets, G, rho, sigma, pi): Gaussian data given the
	# state, then the Gibbs sweep with its sign moves, the path move, the draw of G, the joint scale of the paths,
	# sigma and G, the scale split and common gauge with the reset weights, the persistence and innovation scale
	# given the resets and the reset rate
	nA = 3L; nB = 4L; R = 2L; nt = 4L; s2 = .6
	pr = list(rho_uv_mean = .9, rho_uv_sd = .1, sigma_uv_shape = 3, sigma_uv_scale = 2)
	a = (-.99 - .9) / .1; b = (.99 - .9) / .1
	rho_mean = .9 + .1 * (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
	set.seed(81)
	K = 300
	out = matrix(0, K, 7)
	for (k in seq_len(K)) {
		rho = qnorm(runif(1, pnorm(a), pnorm(b))) * .1 + .9
		sg = sqrt(1 / rgamma(1, 3, 2)); pis = rbeta(1, 1, 9)
		G = matrix(rnorm(R * R), R, R)
		kap = .lame_snap_kappa(NULL, rho, sg)
		pu = snap_prior_paths(nA, R, nt, rho, sg, kap, pis); pv = snap_prior_paths(nB, R, nt, rho, sg, kap, pis)
		U = pu$X; V = pv$X; du = pu$D; dv = pv$D
		for (it in 1:8) {
			EZ = array(0, c(nA, nB, nt))
			for (t in seq_len(nt)) EZ[, , t] = U[, , t] %*% G %*% t(V[, , t])
			E = EZ + array(rnorm(length(EZ), 0, sqrt(s2)), dim(EZ))
			kap = .lame_snap_kappa(NULL, rho, sg)
			o = .lame_snap_bip_kernel(U, V, E, G, rho, sg, s2, kap, pis, du, dv, Inf)
			U = o$U; V = o$V; du = o$delta_u; dv = o$delta_v
			for (t in seq_len(nt)) EZ[, , t] = U[, , t] %*% G %*% t(V[, , t])
			o = snap_path_visit_mode_cpp(EZ, E, numeric(0), array(1, dim(E)), numeric(0), numeric(0), U, V, G, numeric(0),
				rho, sg, kap, pis, du, dv, 1, Inf, s2, numeric(nt), 2L, 1L)
			U = o$U; V = o$V; du = o$delta_u; dv = o$delta_v
			G = sample_G_bip_cpp(E, U, V, lambdaG = 1, s2 = rep(s2, nt))
			o = .lame_uv_bip_scale_step(U, V, G, rho, sg, pr$sigma_uv_shape, pr$sigma_uv_scale, g_var = 1)
			U = o$U; V = o$V; G = o$G; sg = o$sigma
			kw = .lame_uv_dyn_kind_weights("snap", NULL, NULL, du, dv, sg^2 / .lame_snap_kappa(NULL, rho, sg)^2)
			o = .lame_uv_split_dynamic_draw(U, V, rho, sg, kw$w_u, kw$w_v, kw$k_u, kw$k_v)
			U = o$U; V = o$V
			o = .lame_uv_gauge_bip_dynamic_draw(U, V, G, rho, sg, w_u = kw$w_u, w_v = kw$w_v, k_u = kw$k_u, k_v = kw$k_v)
			U = o$U; V = o$V; G = o$G
			h = .lame_snap_hyper_step(U, V, rho, sg, du, dv, pr, Inf)
			rho = h$rho; sg = h$sigma
			pis = .lame_snap_pi_draw(du, dv, 1, 9)
		}
		# coordinates in units of their stationary sd, N(0, 1) under the linked reset scale
		std = (1 - rho^2) / sg^2
		out[k, ] = c(mean(U^2) * std, mean(V^2) * std, mean(G^2), rho, pis, sg^2, mean(c(du[, -1], dv[, -1])))
	}
	expect_lt(max(abs(snap_rep_z(out, c(1, 1, 1, rho_mean, .1, 1, .1)))), 4)
})

test_that("the persistence of bipartite snap paths with the paths following keeps the joint prior under both variance priors", {
	# draws of (rho, sigma, resets, paths) from the snap prior, then Gaussian data given the state and the
	# non-centred persistence draw with the resets held (R/uv_product_scale.R) and the persistence and innovation
	# scale given the resets, eight times; every statistic keeps its mean under the inverse gamma on sigma_uv^2 and
	# under the one on the stationary variance, at the linked and at a fixed reset scale
	nA = 3L; nB = 4L; R = 2L; nt = 4L; s2 = .6
	G = matrix(c(1.1, -.4, .3, .9), R, R)
	a = (-.99 - .9) / .1; b = (.99 - .9) / .1
	ez = function(U, V) {
		E = array(0, c(nA, nB, nt))
		for (t in seq_len(nt)) E[, , t] = U[, , t] %*% G %*% t(V[, , t])
		E
	}
	stat = function(U, V, rho, sg) {
		x = c(rho, log(sg), mean(U^2), mean(V^2), mean(ez(U, V)^2))
		c(x, x^2)
	}
	run = function(stationary, kfix, seed, nc_draw = .lame_uv_nc_rho_draw, K = 400) {
		set.seed(seed)
		pr = list(rho_uv_mean = .9, rho_uv_sd = .1, snap_kappa = kfix)
		if (stationary) { pr$uv_var_shape = 3; pr$uv_var_scale = 2 } else { pr$sigma_uv_shape = 3; pr$sigma_uv_scale = 2 }
		out = matrix(0, K, 10)
		moved = 0
		for (k in seq_len(K)) {
			rho = qnorm(runif(1, pnorm(a), pnorm(b))) * .1 + .9
			v = 1 / rgamma(1, 3, 2)
			sg = sqrt(if (stationary) v * (1 - rho^2) else v)
			kap = .lame_snap_kappa(kfix, rho, sg)
			pu = snap_prior_paths(nA, R, nt, rho, sg, kap, .25); pv = snap_prior_paths(nB, R, nt, rho, sg, kap, .25)
			U = pu$X; V = pv$X; du = pu$D; dv = pv$D
			s0 = stat(U, V, rho, sg)
			for (it in 1:8) {
				EZ = ez(U, V)
				E = EZ + array(rnorm(length(EZ), 0, sqrt(s2)), dim(EZ))
				cx = .lame_uvps_context(E, "normal", "bip")
				o = nc_draw(cx, E, EZ, E, U, V, rho, sg, .9, .1, .lame_uv_ig_shape(pr), .lame_uv_var_scale(pr), s2 = s2,
					G = G, kind = "snap", kappa = kfix, delta_u = du, delta_v = dv, stationary = stationary)
				U = o$U; V = o$V; rho = o$rho; sg = o$sigma
				moved = moved + o$moved
				h = .lame_snap_hyper_step(U, V, rho, sg, du, dv, pr, Inf)
				rho = h$rho; sg = h$sigma
			}
			out[k, ] = stat(U, V, rho, sg) - s0
		}
		list(z = max(abs(colMeans(out) / (apply(out, 2, sd) / sqrt(K)))), moved = moved)
	}
	for (stationary in c(FALSE, TRUE)) for (kfix in list(NULL, 1.5)) {
		r = run(stationary, kfix, 301 + 2 * stationary + is.null(kfix))
		lab = paste("stationary", stationary, "kappa", if (is.null(kfix)) "linked" else kfix)
		expect_gt(r$moved, 1000, label = lab)
		expect_lt(r$z, 4, label = lab)
	}
	# the check fails with the resets drawn as drifts
	as_drift = function(...) {
		a = list(...)
		a$kind = "ar1"
		do.call(.lame_uv_nc_rho_draw, a)
	}
	expect_gt(run(TRUE, 1.5, 309, as_drift, K = 1200)$z, 5)
})

test_that("the bipartite snap path visit draws a Poisson line from its posterior", {
	# one row and one column actor, rank one, three periods: the resets and positions of both paths against
	# importance sampling from the snap prior
	nt = 3L; rho = .8; sg = .4; kap = 2; pis = .2; g = 1.3
	y = c(2, 0, 3)
	v1 = sg^2 / (1 - rho^2)
	prior_draw = function(M) {
		x = matrix(0, M, nt); d = matrix(0, M, nt)
		x[, 1] = rnorm(M, 0, sqrt(v1))
		for (t in 2:nt) {
			d[, t] = runif(M) < pis
			x[, t] = ifelse(d[, t] == 1, rnorm(M, 0, kap), rho * x[, t - 1] + rnorm(M, 0, sg))
		}
		list(x = x, d = d)
	}
	set.seed(91)
	M = 1e6
	pu = prior_draw(M); pv = prior_draw(M)
	lw = rowSums(dpois(matrix(y, M, nt, byrow = TRUE), exp(g * pu$x * pv$x), log = TRUE))
	w = exp(lw - max(lw)); w = w / sum(w)
	ex = c(colSums(w * pu$d)[-1], colSums(w * pv$d)[-1], colSums(w * g * pu$x * pv$x))
	rm(pu, pv)
	U = V = array(0, c(1, 1, nt)); du = dv = matrix(0, 1, nt)
	mask = array(1, c(1, 1, nt)); cnt = array(y, c(1, 1, nt))
	n_it = 100000
	acc = numeric(length(ex))
	for (it in seq_len(n_it)) {
		EZ = array(g * U * V, c(1, 1, nt))
		o = snap_path_visit_mode_cpp(EZ, EZ, numeric(0), mask, cnt, numeric(nt), U, V, matrix(g, 1, 1), numeric(0),
			rho, sg, kap, pis, du, dv, 1, Inf, 1, numeric(nt), 1L, 1L)
		U = o$U; V = o$V; du = o$delta_u; dv = o$delta_v
		acc = acc + c(du[1, -1], dv[1, -1], g * U[1, 1, ] * V[1, 1, ])
	}
	est = acc / n_it
	expect_lt(max(abs(est[1:4] - ex[1:4])), .012)
	expect_lt(max(abs(est[5:7] - ex[5:7])), .02)
})

test_that("bipartite snap fits run in every family the kind serves, with a path move on each side", {
	for (fam in c("binary", "normal", "poisson")) {
		d = snap_bip_panel(fam, seed = 3)
		f = suppressWarnings(suppressMessages(lame(
			d$Y,
			family = fam,
			mode = "bipartite",
			R_row = 1,
			R_col = 1,
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
		moves = f$actor_move
		expect_true(all(is.finite(f$pi_snap)), info = fam)
		expect_equal(dim(f$snap_prob), c(8L, 4L), info = fam)
		expect_equal(dim(f$snap_prob_v), c(10L, 4L), info = fam)
		expect_true(all(is.na(f$snap_prob[, 1])) && all(is.na(f$snap_prob_v[, 1])), info = fam)
		expect_equal(unname(moves["dpm_calls"]), 40 * (8 + 10) / if (fam == "poisson") 2 else 1, info = fam)
	}
})

test_that("two chains of a bipartite snap fit return a gated fit that recovers the planted resets", {
	d = snap_bip_panel("normal", nA = 10L, nB = 12L, N = 5L, seed = 11L, pi_true = .15)
	f = suppressWarnings(suppressMessages(lame(
		d$Y,
		family = "normal",
		mode = "bipartite",
		R_row = 1,
		R_col = 1,
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
	expect_true(all(c("rho_uv", "sigma_uv", "pi_snap") %in% f$rhat$parameter))
	score = c(f$snap_prob[, -1], f$snap_prob_v[, -1])
	truth = c(d$D_u[, -1], d$D_v[, -1])
	expect_gt(snap_auc(score, truth), .7)
})

test_that("reporting, forecasts, simulation and resumption work on a bipartite snap fit", {
	d = snap_bip_panel("binary", seed = 5)
	ck = tempfile(fileext = ".rds")
	f = suppressWarnings(suppressMessages(lame(
		d$Y,
		family = "binary",
		mode = "bipartite",
		R_row = 1,
		R_col = 1,
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
	expect_true(all(is.finite(unlist(fc))))
	expect_no_error(simulate(f, nsim = 2, seed = 1))
	expect_equal(dim(f$snap_draws), c(15L, 8L, 4L))
	expect_equal(dim(f$snap_draws_v), c(15L, 10L, 4L))
	expect_equal(nrow(snap_index_draws(f, side = "v")), 15L * 4L)
	st = f$start_vals$checkpoint_state$values
	expect_equal(dim(st$delta_uv), c(8L, 4L))
	expect_equal(dim(st$delta_uv_v), c(10L, 4L))
	g = suppressWarnings(suppressMessages(lame_resume(ck, nscan_more = 20L)))
	expect_equal(length(g$pi_snap), 5L)
	expect_equal(dim(g$snap_prob_v), c(10L, 4L))
	unlink(ck)
})
