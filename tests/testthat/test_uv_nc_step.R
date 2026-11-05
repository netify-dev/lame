# the non-centred draws of the dynamic position hyperparameters (R/uv_nc_step.R)

# stationarity: draw a state from the prior and the latent array around the position product,
# apply one step, and check that every statistic keeps its mean
ncs_exact = function(mode, mutant = 'none', ndraw = 6000L, seed = 31L) {
	set.seed(seed)
	n = 6L; R = 2L; N = 4L
	if (mode == 'sym') {
		pm = .9; psd = .1; s2 = 1.2; L = c(.8, -.5)
		rtn = function() { repeat { r = rnorm(1, pm, psd); if (abs(r) < .99) return(r) } }
		draw = function() {
			rho = rtn(); U = lame:::.lame_uv_unit_paths(array(rnorm(n * R * N), c(n, R, N)), rho)
			E = array(0, c(n, n, N))
			for (t in 1:N) { Ut = matrix(U[, , t], n, R); Z = Ut %*% diag(L) %*% t(Ut) + matrix(rnorm(n * n, 0, sqrt(s2)), n); Z[lower.tri(Z)] = NA; E[, , t] = Z }
			list(rho = rho, U = U, E = E)
		}
		step = function(s) {
			o = lame:::.lame_uv_sym_nc_rho(
				s$U,
				if (mutant == 'wrongL') L * 1.3 else L,
				s$E,
				if (mutant == 'wrongs2') 2 * s2 else s2,
				s$rho,
				pm,
				psd
			)
			s$rho = o$rho; s$U = o$U; s
		}
		stat = function(s) c(
			s$rho,
			s$rho^2,
			sum(s$U^2),
			sum(s$U[, , 1] * s$U[, , 2]),
			s$U[1, 1, 3] * s$U[2, 1, 4],
			sum(s$U[, , 4]^2),
			s$rho * sum(s$U[, , 1] * s$U[, , 2]),
			sum(s$U[, 1, 2] * s$U[, 1, 4])
		)
	} else {
		sh = 3; sc = 1.5; smax = 3; cap = 4; gv = 1.3; gcap = 5; rho = .85
		rsig = function() { repeat { s = sqrt(1 / rgamma(1, sh, rate = sc)); if (s < smax && s^2 / (1 - rho^2) <= cap) return(s) } }
		draw = function() {
			sg = rsig(); w = sg^2 / (1 - rho^2)
			U = array(rnorm(n * R * N), c(n, R, N)); V = array(rnorm(n * R * N), c(n, R, N))
			U[, , 1] = sqrt(w) * U[, , 1]; V[, , 1] = sqrt(w) * V[, , 1]
			for (t in 2:N) { U[, , t] = rho * U[, , t - 1] + sg * U[, , t]; V[, , t] = rho * V[, , t - 1] + sg * V[, , t] }
			repeat { G = matrix(rnorm(R * R, 0, sqrt(gv)), R); if (max(abs(G)) <= gcap) break }
			list(U = U, V = V, G = G, sigma = sg)
		}
		step = function(s) {
			o = lame:::.lame_uv_bip_scale_step(
				s$U,
				s$V,
				s$G,
				rho,
				s$sigma,
				if (mutant == 'wrongshape') sh + 1 else sh,
				sc,
				smax,
				cap,
				if (mutant == 'wrongG') 4 * gv else gv,
				gcap
			)
			for (t in 1:N) stopifnot(
				max(abs(o$U[, , t] %*% o$G %*% t(o$V[, , t]) - s$U[, , t] %*% s$G %*% t(s$V[, , t]))) < 1e-8
			)
			s$U = o$U; s$V = o$V; s$G = o$G; s$sigma = o$sigma; s
		}
		stat = function(s) c(
			s$sigma,
			s$sigma^2,
			log(s$sigma),
			s$G[1, 1],
			s$G[1, 2]^2,
			sum(s$G^2),
			sum(s$U^2),
			sum(s$V[, , 2]^2),
			s$U[1, 1, 3] * s$V[2, 2, 3],
			sum(s$U[, , 1]^2) * s$sigma,
			sum(s$G^2) * s$sigma^2,
			max(abs(s$G))
		)
	}
	S0 = S1 = NULL
	for (r in seq_len(ndraw)) { s = draw(); S0 = rbind(S0, stat(s)); S1 = rbind(S1, stat(step(s))) }
	d = S1 - S0
	max(abs(colMeans(d) / (apply(d, 2L, sd) / sqrt(ndraw))))
}

test_that('the non-centred draw of the symmetric persistence and the bipartite scale step leave the prior of the state unchanged', {
	skip_on_cran()
	expect_lt(ncs_exact('sym'), 4.2)
	expect_lt(ncs_exact('bip'), 4.2)
})

test_that('wrong variance, wrong weights, wrong inverse-gamma shape and a wrong prior of G are caught', {
	skip_on_cran()
	expect_gt(ncs_exact('sym', 'wrongL', 5000L), 5)
	expect_gt(ncs_exact('sym', 'wrongs2', 8000L), 4.5)
	expect_gt(ncs_exact('bip', 'wrongshape', 5000L), 5)
	expect_gt(ncs_exact('bip', 'wrongG', 5000L), 5)
})

test_that('the paths of the unit gauge reproduce their innovations', {
	set.seed(2)
	eps = array(rnorm(5 * 2 * 4), c(5, 2, 4))
	U = lame:::.lame_uv_unit_paths(eps, .7)
	expect_equal((U[, , 3] - .7 * U[, , 2]) / sqrt(1 - .49), eps[, , 3])
	expect_equal(U[, , 1], eps[, , 1])
})

test_that('lame() runs the steps for symmetric and bipartite dynamic positions and the switch omits them', {
	skip_on_cran()
	set.seed(3)
	n = 10L; N = 4L
	nm = sprintf('a%02d', 1:n)
	sym = setNames(
		lapply(
			1:N,
			function(t) { Y = matrix(rbinom(n * n, 1, .3), n); Y[lower.tri(Y)] = t(Y)[lower.tri(Y)]; diag(Y) = NA; dimnames(Y) = list(nm, nm); Y }
		),
		paste0('p', 1:N)
	)
	f = suppressWarnings(
		suppressMessages(
			lame(
				sym,
				family = 'binary',
				symmetric = TRUE,
				R = 1L,
				dynamic_uv = TRUE,
				burn = 30,
				nscan = 50,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1
			)
		)
	)
	expect_equal(unname(f$uv_nc_move['sym_calls']), 80)
	f0 = suppressWarnings(
		suppressMessages(
			lame(
				sym,
				family = 'binary',
				symmetric = TRUE,
				R = 1L,
				dynamic_uv = TRUE,
				burn = 30,
				nscan = 50,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1,
				prior = list(uv_nc_move = FALSE)
			)
		)
	)
	expect_null(f0$uv_nc_move)
	rn = sprintf('r%02d', 1:10); cn = sprintf('c%02d', 1:8)
	bip = setNames(
		lapply(1:N, function(t) { Y = matrix(rbinom(80, 1, .3), 10, 8); dimnames(Y) = list(rn, cn); Y }),
		paste0('p', 1:N)
	)
	fb = suppressWarnings(
		suppressMessages(
			lame(
				bip,
				family = 'binary',
				mode = 'bipartite',
				R_row = 1L,
				R_col = 1L,
				dynamic_uv = TRUE,
				burn = 30,
				nscan = 50,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE,
				seed = 1
			)
		)
	)
	expect_equal(unname(fb$uv_nc_move['bip_calls']), 80)
	expect_error(
		lame(
			bip,
			family = 'binary',
			mode = 'bipartite',
			R_row = 1L,
			R_col = 1L,
			dynamic_uv = TRUE,
			burn = 5,
			nscan = 5,
			odens = 1,
			verbose = FALSE,
			gof = FALSE,
			plot = FALSE,
			prior = list(uv_nc_move = 'x')
		),
		'TRUE or FALSE'
	)
})

# count-panel version: the Poisson factor of the observed cells is the data term,
# and the joint law of (persistence, paths, latent array, counts) stays stationary
ncc_exact = function(mutant = 'none', ndraw = 5000L, seed = 41L) {
	set.seed(seed)
	n = 6L; R = 2L; N = 4L; pm = .9; psd = .1; s2 = .4; L = c(.8, -.5); off = c(0, .3, -.2, .1)
	rtn = function() { repeat { r = rnorm(1, pm, psd); if (abs(r) < .99) return(r) } }
	prod_of = function(U) { P = array(0, c(n, n, N)); for (t in 1:N) { Ut = matrix(U[, , t], n, R); P[, , t] = Ut %*% diag(L) %*% t(Ut) }; P }
	draw = function() {
		rho = rtn(); U = lame:::.lame_uv_unit_paths(array(rnorm(n * R * N), c(n, R, N)), rho)
		P = prod_of(U); Z = Y = array(0, c(n, n, N))
		for (t in 1:N) { z = P[, , t] + matrix(rnorm(n * n, 0, sqrt(s2)), n); z[lower.tri(z)] = t(z)[lower.tri(z)]; Z[, , t] = z; Y[, , t] = matrix(
			rpois(n * n, exp(z + off[t])),
			n
		) }
		list(rho = rho, U = U, Z = Z, Y = Y)
	}
	step = function(s) {
		yy = s$Y; for (t in 1:N) yy[, , t][lower.tri(yy[, , t], TRUE)] = NA
		cx = lame:::.lame_ltm_context(yy, 'poisson', 'sym', offset = off)
		if (mutant == 'wrongoff') cx$off_ok = cx$off_ok + .5
		o = lame:::.lame_uv_sym_nc_rho_carry(s$U, L, s$Z, cx, s$rho, pm, psd)
		s$rho = o$rho; s$U = o$U; s$Z = if (mutant == 'nocarry') s$Z else o$Z
		s
	}
	stat = function(s) { e = s$Z - prod_of(s$U); c(
		s$rho,
		s$rho^2,
		sum(s$U^2),
		sum(s$U[, , 1] * s$U[, , 2]),
		e[1, 2, 1],
		e[3, 4, 3]^2,
		s$Z[2, 2, 2],
		sum(s$Z[, , 4]),
		sum(s$Y[, , 1] * s$Z[, , 1]),
		e[2, 2, 2]^2
	) }
	S0 = S1 = NULL
	for (r in seq_len(ndraw)) { s = draw(); S0 = rbind(S0, stat(s)); S1 = rbind(S1, stat(step(s))) }
	d = S1 - S0
	sdv = apply(d, 2L, sd)
	max(abs(colMeans(d)[sdv > 0] / (sdv[sdv > 0] / sqrt(ndraw))))
}

test_that('the persistence draw that carries the latent array of a count panel leaves the joint law of the state and the counts unchanged', {
	skip_on_cran()
	expect_lt(ncc_exact(), 4.2)
	expect_gt(ncc_exact('nocarry', 6000L), 4.5)
	expect_gt(ncc_exact('wrongoff', 3000L), 6)
})
