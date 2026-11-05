# shared states and panels for the position-path move tests

dpmt_ez = function(mode, base, U, V, G, Ld) {
	nr = dim(base)[1]; nc = dim(base)[2]
	out = base
	for (t in seq_len(dim(base)[3])) {
		u = matrix(U[, , t], nr, dim(U)[2])
		out[, , t] = base[, , t] + if (mode == 'sym') u %*% diag(Ld, length(Ld)) %*% t(u)
		else if (mode == 'bip') u %*% G %*% t(matrix(V[, , t], nc, dim(V)[2]))
		else u %*% t(matrix(V[, , t], nc, dim(V)[2]))
	}
	out
}
# a state of every mode: counts with missing cells and absent actor-periods, positions, predictor, latent array
dpmt_state = function(mode, rk = 1L, seed = 3L, nt = 3L) {
	set.seed(seed)
	nr = if (mode == 'bip') 5L else 6L; nc = if (mode == 'bip') 6L else nr
	Y = array(rpois(nr * nc * nt, .7), c(nr, nc, nt)); Y[1:2, , ] = 0
	if (mode == 'sym') for (t in seq_len(nt)) { m = Y[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Y[, , t] = m }
	if (mode != 'bip') for (t in seq_len(nt)) diag(Y[, , t]) = NA
	Y[3, 4, 2] = NA; Y[4, 3, 2] = if (mode == 'sym') NA else Y[4, 3, 2]
	if (mode != 'bip') Y[5, , 3] = NA
	if (mode == 'sym') for (t in seq_len(nt)) { m = Y[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Y[, , t] = m }
	off = c(0, .2, -.3)[seq_len(nt)]
	cx = .lame_ltm_context(Y, 'poisson', mode, offset = off)
	U = array(rnorm(nr * rk * nt, 0, .6), c(nr, rk, nt))
	V = if (mode == 'sym') NULL else array(rnorm(nc * rk * nt, 0, .6), c(nc, rk, nt))
	G = if (mode == 'bip') matrix(rnorm(rk * rk, 0, .8), rk, rk) else NULL
	Ld = if (mode == 'sym') c(0.9, -0.6)[seq_len(rk)] else NULL
	base = array(rnorm(nr * nc * nt, -.5, .4), c(nr, nc, nt))
	if (mode == 'sym') base = array(
		vapply(seq_len(nt), function(t) { m = base[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; m }, numeric(nr * nc)),
		c(nr, nc, nt)
	)
	EZ = dpmt_ez(mode, base, U, V, G, Ld)
	e = array(rnorm(nr * nc * nt, 0, .4), c(nr, nc, nt))
	if (mode == 'sym') e = array(
		vapply(seq_len(nt), function(t) { m = e[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; m }, numeric(nr * nc)),
		c(nr, nc, nt)
	)
	list(
		mode = mode,
		cx = cx,
		base = base,
		U = U,
		V = V,
		G = G,
		Ld = Ld,
		EZ = EZ,
		Z = EZ + e,
		rho_uv = .8,
		sigma_uv = if (mode == 'sym') sqrt(1 - .8^2) else .5
	)
}
dpmt_panel = function(mode, nt = 3L, binary = FALSE) {
	set.seed(5)
	nr = if (mode == 'bip') 7L else 8L; nc = if (mode == 'bip') 9L else 8L
	rn = sprintf('r%02d', 1:nr); cn = if (mode == 'bip') sprintf('c%02d', 1:nc) else rn
	Y = lapply(1:nt, function(t) {
		y = matrix(rpois(nr * nc, if (binary) .5 else .4), nr, nc, dimnames = list(rn, cn)); y[1:2, ] = 0; y[, 1:2] = 0
		if (binary) y[y > 0] = 1
		if (mode == 'sym') y[lower.tri(y)] = t(y)[lower.tri(y)]
		if (mode != 'bip') diag(y) = NA
		y
	})
	X = lapply(
		1:nt,
		function(t) { x = matrix(
			rnorm(nr * nc),
			nr,
			nc,
			dimnames = list(rn, cn)
		); if (mode == 'sym') x[lower.tri(x)] = t(x)[lower.tri(x)]; array(
			x,
			c(nr, nc, 1),
			dimnames = list(rn, cn, 'x')
		) }
	)
	list(Y = Y, X = X)
}
dpmt_run = function(p, mode, prior = list(), family = 'poisson', rk = 1L, ...) suppressMessages(suppressWarnings(lame(
	p$Y,
	Xdyad = p$X,
	mode = if (mode == 'bip') 'bipartite' else 'unipartite',
	symmetric = mode == 'sym',
	family = family,
	dynamic_uv = TRUE,
	R = if (mode == 'bip') 0L else rk,
	R_row = if (mode == 'bip') rk else NULL,
	R_col = if (mode == 'bip') rk else NULL,
	burn = 10,
	nscan = 30,
	odens = 5,
	seed = 3,
	verbose = FALSE,
	plot = FALSE,
	gof = FALSE,
	prior = modifyList(list(actor_move_every = 2, als_start = FALSE), prior),
	...
)))
# binary state of every mode: missing cells, absent actor-periods, inactive actors, positions, predictor and a
# latent array consistent with the signs of the data
dpmt_state_bin = function(mode, rk = 1L, seed = 3L, nt = 3L, rho = .4) {
	set.seed(seed)
	nr = if (mode == 'bip') 5L else 6L; nc = if (mode == 'bip') 6L else nr
	Y = array(rbinom(nr * nc * nt, 1, .45), c(nr, nc, nt)); Y[1:2, , ] = 0
	if (mode == 'sym') for (t in seq_len(nt)) { m = Y[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Y[, , t] = m }
	if (mode != 'bip') for (t in seq_len(nt)) diag(Y[, , t]) = NA
	Y[3, 4, 2] = NA; Y[4, 3, 2] = if (mode == 'sym') NA else Y[4, 3, 2]
	if (mode != 'bip') Y[5, , 3] = NA
	if (mode == 'sym') for (t in seq_len(nt)) { m = Y[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; Y[, , t] = m }
	cx = .lame_dpm_binary_context(Y, mode)
	U = array(rnorm(nr * rk * nt, 0, .6), c(nr, rk, nt))
	V = if (mode == 'sym') NULL else array(rnorm(nc * rk * nt, 0, .6), c(nc, rk, nt))
	G = if (mode == 'bip') matrix(rnorm(rk * rk, 0, .8), rk, rk) else NULL
	Ld = if (mode == 'sym') c(0.9, -0.6)[seq_len(rk)] else NULL
	base = array(rnorm(nr * nc * nt, -.3, .4), c(nr, nc, nt))
	if (mode == 'sym') base = array(
		vapply(seq_len(nt), function(t) { m = base[, , t]; m[lower.tri(m)] = t(m)[lower.tri(m)]; m }, numeric(nr * nc)),
		c(nr, nc, nt)
	)
	EZ = dpmt_ez(mode, base, U, V, G, Ld)
	Z0 = array(0, dim(EZ))
	Z = if (mode == 'dir') .lame_binary_z_given_factor(Z0, EZ, Z0, rho, cx$i1, cx$i0, cx$ina, cx$idg)
	else if (mode == 'sym') { for (t in seq_len(nt)) Z0[, , t] = rZ_bin_sym_fc(
		Z0[, , t],
		EZ[, , t],
		Y[, , t],
		rho
	); Z0 }
	else rZ_bin_bip_batch_cpp(Z0, EZ, Y)
	list(
		mode = mode,
		cx = cx,
		Y = Y,
		base = base,
		U = U,
		V = V,
		G = G,
		Ld = Ld,
		EZ = EZ,
		Z = Z,
		rho = if (mode == 'bip') 0 else rho,
		rho_uv = .8,
		sigma_uv = if (mode == 'sym') sqrt(1 - .8^2) else .5
	)
}
