# fixed references use 90-digit inversion of the original dyad likelihood
.joint_root_reference = function() list(
	list(
		n = 2L,
		rho = 0.99999999999999645,
		prior = matrix(c(1.1000000000000001), 1L),
		Sab = matrix(c(1, 0.20000000000000001, 0.20000000000000001, 0.69999999999999996), 2L),
		z = c(0.56648771063468495, 0.68640770679600682),
		mean = c(0.29830843272159352, -0.052182303068021678, 0.032613939417513546),
		covariance = matrix(
			c(
				0.47619047619047578, 0, 0, 0, 0.507692307692309, 0.50769230769230689, 0, 0.50769230769230689, 0.50769230769230822
			),
			3L
		)
	),
	list(
		n = 3L,
		rho = -0.99999999999999645,
		prior = matrix(c(1.1000000000000001), 1L),
		Sab = matrix(c(1, 0.20000000000000001, 0.20000000000000001, 0.69999999999999996), 2L),
		z = c(
			0.56648771063468495, 0.84697345837817473, 0.81687462335142302, 0.9072636290606505, 0.65337565025596978, 0.48320730939498541
		),
		mean = c(
			0.71236373017931431, -0.083525231834090691, 0.088378397088270652, 0.0058296325964979755, -0.037716670709060075
		),
		covariance = matrix(
			c(
				5.9211894646674973e-16, -8.5865865581780984e-48, -2.1171634588652062e-31, -6.8957204197963707e-48, -1.8271410672398359e-31, -8.5865865581780984e-48, 0.16176470588235486, 3.0814422713139391e-18, -0.16176470588235103, -3.0814422713137843e-18, -2.1171634588652062e-31, 3.0814422713139391e-18, 0.16176470588235484, -3.0814422713137843e-18, -0.16176470588235101, -6.8957204197963707e-48, -0.16176470588235103, -3.0814422713137843e-18, 0.16176470588235434, 3.0814422713139087e-18, -1.8271410672398359e-31, -3.0814422713137843e-18, -0.16176470588235101, 3.0814422713139087e-18, 0.16176470588235431
			),
			5L
		)
	),
	list(
		n = 3L,
		rho = 0.99999999999999645,
		prior = matrix(c(1.1000000000000001), 1L),
		Sab = matrix(c(1, 0.20000000000000001, 0.20000000000000001, 0.69999999999999996), 2L),
		z = c(
			0.56648771063468495, 0.84697345837817473, 0.81687462335142302, 0.9072636290606505, 0.65337565025596978, 0.48320730939498541
		),
		mean = c(
			0.52124175378974269, -0.11244824865261913, 0.15390557290930632, 0.059904866047139994, -0.098250668858039542
		),
		covariance = matrix(
			c(
				0.24390243902438993, -3.6156683904836278e-34, -2.7354922707289823e-17, -3.6156683904836201e-34, -2.7354922707289807e-17, -3.6156683904836278e-34, 0.33673469387755162, 4.4508295697278342e-18, 0.33673469387755023, 4.450829569727778e-18, -2.7354922707289823e-17, 4.4508295697278342e-18, 0.33673469387755156, 4.450829569727778e-18, 0.33673469387755023, -3.6156683904836201e-34, 0.33673469387755023, 4.450829569727778e-18, 0.33673469387755128, 4.4508295697278149e-18, -2.7354922707289807e-17, 4.450829569727778e-18, 0.33673469387755023, 4.4508295697278149e-18, 0.33673469387755123
			),
			5L
		)
	),
	list(
		n = 3L,
		rho = 0.69999999999999996,
		prior = matrix(c(1), 1L),
		Sab = matrix(c(1, 0, 0, 9.9999999999999995e-21), 2L),
		z = c(
			0.90929742682568171, 0.14112000805986721, -0.7568024953079282, -0.27941549819892586, 0.65698659871878906, 0.98935824662338179
		),
		mean = c(
			0.21565510217154102, 0.76198275454168063, -0.55582482438777514, -1.2468792349178605e-20, 1.1716905922651417e-20
		),
		covariance = matrix(
			c(
				0.22077922077922077, -2.2689758377220654e-34, -6.8723600262025712e-18, 2.0168674113085024e-54, -7.5595960288228271e-37, -2.2689758377220654e-34, 0.15887850467289724, 5.2455268184692462e-18, 7.4766355140186903e-21, -4.6626905053059959e-38, -6.8723600262025712e-18, 5.2455268184692462e-18, 0.15887850467289719, -4.6626905053059959e-38, 7.4766355140186903e-21, 2.0168674113085024e-54, 7.4766355140186903e-21, -4.6626905053059959e-38, 9.9999999999999995e-21, 4.7758293190376851e-57, -7.5595960288228271e-37, -4.6626905053059959e-38, 7.4766355140186903e-21, 4.7758293190376851e-57, 9.9999999999999995e-21
			),
			5L
		)
	)
)

test_that("joint coefficient and actor roots match independent boundary targets", {
	for (f in .joint_root_reference()) {
		n = f$n
		X = list(array(1, c(n, n, 1L)))
		Z = array(0, c(n, n, 1L))
		Z[, , 1L][row(Z[, , 1L]) != col(Z[, , 1L])] = f$z
		zero = matrix(0, n, 1L)
		cx = lame:::.lame_joint_levels_context(X, f$prior)
		withr::local_seed(47339331L)
		before = .Random.seed
		value = lame:::.lame_joint_levels_parameters(Z, zero, zero, f$Sab, f$rho, cx)
		expect_identical(.Random.seed, before)
		sd = sqrt(diag(f$covariance))
		expect_lt(max(abs(as.numeric(value$mean) - f$mean) / sd), 2e-7)
		expect_lt(max(abs(tcrossprod(value$noise_root) - f$covariance) / outer(sd, sd)), 2e-7)
		changed = Z
		diag(changed[, , 1L]) = 1e9
		expect_identical(
			lame:::.lame_joint_levels_parameters(
				changed,
				zero,
				zero,
				f$Sab,
				f$rho,
				cx
			),
			value
		)
	}
})

test_that("joint level draws retain centered support and completed-call replay", {
	f = .joint_root_reference()[[3L]]
	n = f$n
	X = list(array(1, c(n, n, 1L)))
	Z = array(0, c(n, n, 1L))
	Z[, , 1L][row(Z[, , 1L]) != col(Z[, , 1L])] = f$z
	zero = matrix(0, n, 1L)
	cx = lame:::.lame_joint_levels_context(X, f$prior)
	draw = function() lame:::.lame_joint_levels_draw(Z, zero, zero, f$Sab, f$rho, cx)
	withr::local_seed(47339332L)
	initial = .Random.seed
	whole = lapply(1:12, function(i) draw())
	ending = .Random.seed
	assign(".Random.seed", initial, envir = .GlobalEnv)
	prefix = lapply(1:5, function(i) draw())
	checkpoint = unserialize(serialize(list(prefix = prefix, RNG = .Random.seed), NULL))
	assign(".Random.seed", checkpoint$RNG, envir = .GlobalEnv)
	suffix = lapply(1:7, function(i) draw())
	expect_identical(whole, c(checkpoint$prefix, suffix))
	expect_identical(.Random.seed, ending)
	expect_lt(max(vapply(whole, function(x) abs(sum(x$a)), 0)), 1e-12)
	expect_lt(max(vapply(whole, function(x) abs(sum(x$b)), 0)), 1e-12)
	expect_true(all(vapply(whole, function(x) all(is.finite(x$q)), TRUE)))
	before = .Random.seed
	expect_error(lame:::.lame_joint_levels_draw(Z, zero, zero, f$Sab, 1, cx), "Invalid")
	expect_identical(.Random.seed, before)
})


test_that("deterministic joint moments and validation preserve an absent RNG state", {
	had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	old_seed = if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
	withr::defer({
		if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
		else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
			rm(".Random.seed", envir = .GlobalEnv)
	})
	if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
	X = list(array(1, c(3L, 3L, 1L)))
	Z = array(sin(seq_len(9L)), c(3L, 3L, 1L))
	zero = matrix(0, 3L, 1L)
	cx = lame:::.lame_joint_levels_context(X, matrix(1))
	expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
	value = lame:::.lame_joint_levels_parameters(Z, zero, zero, diag(2), .7, cx)
	expect_true(all(is.finite(value$mean)))
	expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
	expect_error(lame:::.lame_joint_levels_draw(Z, zero, zero, diag(2), 1, cx), "Invalid")
	expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

# static columns under their prior, coefficient paths as an offset, against the dense gaussian conditional
test_that("joint levels with the coefficient paths as an offset match the dense Gaussian conditional", {
	set.seed(58317)
	n = 4L; nt = 2L; k = n - 1L; rho = .35
	Sab = matrix(c(.6, .12, .12, .45), 2)
	H = qr.Q(qr(contr.helmert(n)))
	X = lapply(1:nt, function(t) {
		x = array(1, c(n, n, 3L))
		x[, , 2] = rnorm(n * n)
		x[, , 3] = rnorm(n * n)
		x
	})
	path = cbind(cumsum(c(rnorm(1, 0, 1.2), rnorm(1, 0, .5))), c(.4, -.2))
	Z = array(rnorm(n * n * nt), c(n, n, nt))
	for (case in c("rw1", "mixed", "empty")) {
		dyn = switch(case, rw1 = 1L, mixed = c(1L, 3L), empty = c(1L, 2L, 3L))
		st = setdiff(1:3, dyn)
		vs = c(2.2, 1.4)[seq_along(st)]
		offset = array(0, c(n, n, nt))
		for (t in 1:nt) for (j in seq_along(dyn)) {
			coef_t = if (j == 1) path[t, 1] else path[t, 2]
			offset[, , t] = offset[, , t] + coef_t * X[[t]][, , dyn[j]]
		}
		cx = lame:::.lame_joint_levels_static_context(X, st, vs)
		m = lame:::.lame_joint_levels_parameters(Z, matrix(0, n, nt), matrix(0, n, nt), Sab, rho, cx, offset)
		# dense reference over (static coefficients, sender contrasts, receiver contrasts)
		ps = length(st)
		q = ps + 2L * k
		cells = which(row(diag(n)) != col(diag(n)))
		D = NULL; y = NULL
		for (t in 1:nt) {
			ij = arrayInd(cells, c(n, n))
			Xt = if (ps) matrix(X[[t]][, , st], n * n, ps)[cells, , drop = FALSE] else NULL
			Dt = cbind(Xt, H[ij[, 1], , drop = FALSE], H[ij[, 2], , drop = FALSE])
			D = rbind(D, Dt)
			y = c(y, (Z[, , t] - offset[, , t])[cells])
		}
		cellid = unlist(lapply(1:nt, function(t) cells + (t - 1) * n * n))
		Sig = diag(length(y))
		for (a in seq_along(y)) {
			ij = arrayInd(cellid[a] - n * n * ((cellid[a] - 1) %/% (n * n)), c(n, n))
			t = (cellid[a] - 1) %/% (n * n)
			b = which(cellid == ij[2] + (ij[1] - 1) * n + t * n * n)
			Sig[a, b] = rho
		}
		P0 = matrix(0, q, q)
		if (ps) P0[1:ps, 1:ps] = diag(1 / vs, ps)
		iS = solve(Sab)
		for (h in 1:k) {
			ix = ps + c(h, k + h)
			P0[ix, ix] = iS
		}
		Pq = crossprod(D, solve(Sig, D)) + P0
		Vq = solve(Pq)
		mq = Vq %*% crossprod(D, solve(Sig, y))
		got_m = as.numeric(m$mean)
		got_V = tcrossprod(m$noise_root)
		if (case == "empty") {
			# the stand-in coefficient keeps its prior, independent of the levels
			expect_equal(got_m[1], 0, tolerance = 1e-10)
			expect_equal(got_V[1, 1], 1, tolerance = 1e-10)
			expect_lt(max(abs(got_V[1, -1])), 1e-10)
			got_m = got_m[-1]
			got_V = got_V[-1, -1]
		}
		expect_equal(got_m, as.numeric(mq), tolerance = 1e-8)
		expect_equal(got_V, Vq, tolerance = 1e-8)
		zero = matrix(0, n, nt)
		dr = lame:::.lame_joint_levels_draw(Z, zero, zero, Sab, rho, cx, offset)
		expect_length(dr$beta, ps)
		expect_lt(abs(sum(dr$a)), 1e-10)
	}
})
