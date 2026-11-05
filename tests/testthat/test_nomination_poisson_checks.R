skip_on_cran()

# nomination (cbin / frn) and Poisson checks: absent actors' cells are not nominations, square bipartite slices stay
# bipartite, odmax below an observed out-degree, non-integer counts, the prediction summary, bipartite rank bounds and the symmetric reciprocity ramp on resume

# a cap-2 slice with actor 5 absent: rows 2 to 4 are below the cap, so their
# observed non-nominations sit below zero exactly
.absent_cbin_slice = function() {
	matrix(
		c(
			NA, 1, 1, 0, NA,
			0, NA, 1, 0, NA,
			0, 0, NA, 0, NA,
			1, 0, 0, NA, NA,
			NA, NA, NA, NA, NA
		),
		5,
		5,
		byrow = TRUE
	)
}

test_that("an absent actor's cells are not counted as nominations (sampler, exact)", {
	Y = .absent_cbin_slice()
	n = nrow(Y); cap = 2; mu = 0.3
	EZ = matrix(mu, n, n)
	hid = lame:::.lame_nom_hidden_cells(Y, FALSE)
	expect_false(any(hid))
	odobs = rowSums(Y > 0, na.rm = TRUE)
	set.seed(11)
	Z = matrix(rnorm(n * n), n, n)
	room = odobs < cap
	zero = !is.na(Y) & Y == 0 & matrix(room, n, n)
	one = !is.na(Y) & Y == 1 & matrix(room, n, n)
	S = 1500; acc = matrix(0, n, n); n_pos = 0
	for (s in seq_len(S)) {
		Z = lame:::rZ_cbin_fc(Z, EZ, 0, Y, cap, lame:::.lame_nom_room(odobs, hid, Z))
		n_pos = n_pos + sum(Z[zero] > 0)
		acc = acc + Z
	}
	# rows below the cap: every observed non-nomination stays below zero
	expect_identical(n_pos, 0)
	# and with rho = 0 the cells are independent truncated normals
	gm = acc / S
	expect_lt(max(abs(gm[zero] - (mu - dnorm(mu) / pnorm(-mu)))), 0.06)
	expect_lt(max(abs(gm[one] - (mu + dnorm(mu) / pnorm(mu)))), 0.06)
})

test_that("lame() cbin with an absent actor keeps observed non-nominations below zero", {
	skip_on_cran()
	set.seed(5)
	n = 10; nm = sprintf("a%02d", 1:n)
	mk = function() {
		z = matrix(rnorm(n * n, 0.6), n, n); y = matrix(0, n, n)
		for (i in 1:n) {
			o = order(z[i, ], decreasing = TRUE); o = o[o != i][1:2]
			y[i, o[z[i, o] > 0]] = 1
		}
		diag(y) = NA; dimnames(y) = list(nm, nm); y
	}
	Y = list(mk(), mk(), mk())
	Y[[2]]["a10", ] = NA; Y[[2]][, "a10"] = NA
	Y[[2]][1:4, ] = 0; diag(Y[[2]]) = NA; Y[[2]][1:4, "a10"] = NA; Y[[2]][1:4, 5] = 1
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		family = "cbin",
		odmax = 2,
		R = 0,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 1
	)))
	Z = fit$start_vals$Z[, , 2]; y = Y[[2]]
	room = rowSums(y > 0, na.rm = TRUE) < 2
	cells = !is.na(y) & y == 0 & matrix(room, n, n)
	expect_gt(sum(cells), 20)
	expect_true(all(Z[cells] < 0))
})

test_that("square bipartite nomination slices keep the diagonal as ordinary cells", {
	Y = matrix(
		c(
			1, 0, NA,
			0, 1, 0,
			0, 0, 1
		),
		3,
		3,
		byrow = TRUE
	)
	Y2 = Y; Y2[2, 2] = NA
	# a hidden diagonal cell is a hidden cell in a bipartite slice only
	expect_true(lame:::.lame_nom_hidden_cells(Y2, TRUE)[2, 2])
	expect_false(lame:::.lame_nom_hidden_cells(Y2, FALSE)[2, 2])
	expect_identical(lame:::.lame_nom_predictive_yo(Y2, TRUE)[2, 2], 0)
	# row 1 is at its cap of 1 with its nomination on the diagonal: the
	# hidden cell (1, 3) can never outrank it while positive (exact share 0)
	EZ = matrix(0.5, 3, 3)
	for (fam in c("cbin", "frn")) {
		set.seed(2)
		Z = matrix(rnorm(9), 3, 3)
		bad = 0
		for (s in 1:1500) {
			odobs = rowSums(Y > 0, na.rm = TRUE)
			Z = if (fam == "cbin") lame:::rZ_cbin_bip_fc(Z, EZ, Y, 1, odobs)
			else lame:::rZ_frn_bip_fc(Z, EZ, Y, matrix(c(1, 2, 3), 3, 1), 1, odobs)
			bad = bad + (Z[1, 3] > 0 && Z[1, 3] > Z[1, 1])
			# the observed diagonal nomination stays positive
			bad = bad + (Z[1, 1] <= 0)
		}
		expect_identical(bad, 0, label = fam)
	}
})

test_that("odmax below an observed out-degree is refused in ame() and lame()", {
	set.seed(5); n = 8
	Y = replicate(2, {
		y = matrix(rbinom(n * n, 1, 0.4), n, n); diag(y) = NA
		dimnames(y) = list(letters[1:n], letters[1:n]); y
	}, simplify = FALSE)
	expect_error(
		suppressWarnings(suppressMessages(lame(
			Y,
			family = "cbin",
			odmax = 2,
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		))),
		class = "lame_odmax_below_outdegree"
	)
	expect_error(
		suppressWarnings(suppressMessages(ame(
			Y[[1]],
			family = "frn",
			odmax = 2,
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		))),
		class = "lame_odmax_below_outdegree"
	)
	expect_error(
		suppressWarnings(suppressMessages(ame(
			Y[[1]],
			family = "cbin",
			odmax = 2,
			mode = "bipartite",
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		))),
		class = "lame_odmax_below_outdegree"
	)
	# the message names the rows and the largest out-degree
	expect_error(
		suppressWarnings(suppressMessages(ame(
			Y[[1]],
			family = "cbin",
			odmax = 2,
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		))),
		"out-degree"
	)
})

test_that("non-integer and negative poisson counts are refused at every entry", {
	set.seed(6); n = 8
	Yp = matrix(rpois(n * n, 3) + 0.5, n, n); diag(Yp) = NA
	dimnames(Yp) = list(letters[1:n], letters[1:n])
	quiet = function(expr) suppressWarnings(suppressMessages(expr))
	expect_error(
		quiet(ame(
			Yp,
			family = "poisson",
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		)),
		class = "lame_poisson_counts"
	)
	expect_error(
		quiet(ame(
			Yp,
			family = "poisson",
			mode = "bipartite",
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		)),
		class = "lame_poisson_counts"
	)
	Ys = (Yp + t(Yp)) / 2
	expect_error(
		quiet(ame(
			Ys,
			family = "poisson",
			symmetric = TRUE,
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		)),
		class = "lame_poisson_counts"
	)
	expect_error(
		quiet(lame(
			list(Yp, Yp),
			family = "poisson",
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		class = "lame_poisson_counts"
	)
	expect_error(
		quiet(lame(list(Yp, Yp), family = "poisson", method = "als", R = 0, verbose = FALSE)),
		class = "lame_poisson_counts"
	)
	Yn = round(Yp); Yn[1, 2] = -1
	expect_error(
		quiet(lame(
			list(Yn, Yn),
			family = "poisson",
			R = 0,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		"non-negative counts"
	)
})

test_that("the poisson log mass is finite and exact at non-integer counts", {
	y = c(0, 3, 2.5, 7.25)
	eta = c(0.2, 0.4, 0.7, 1.9)
	ref = y * eta - exp(eta) - lgamma(y + 1)
	got = lame:::.poisson_log_mass(y, eta)
	expect_equal(got, ref, tolerance = 1e-10)
	# integer counts keep the dpois() value exactly
	expect_identical(got[2], dpois(3, 3, log = TRUE) - 3 * (expm1(0.4 - log(3)) - (0.4 - log(3))))
	expect_identical(lame:::.poisson_log_mass(-1, 0), -Inf)
	# compiled kernel: with rho = 0, s2 = 1 the (1, 2) entry is the two cells'
	# latent and count terms
	Z = matrix(c(0, 0.7, 1.9, 0), 2, 2); Y = matrix(c(NA, 2.5, 7.25, NA), 2, 2)
	ld = lame:::ldZgbme_cpp(Z, Y, matrix(0, 2, 2), 0, 1)
	expect_equal(ld[1, 2], -0.5 * (0.7^2 + 1.9^2) + sum(ref[3:4]), tolerance = 1e-10)
	# the latent draw moves on non-integer counts
	set.seed(6); n = 6
	Yp = matrix(rpois(n * n, 3) + 0.5, n, n); diag(Yp) = NA
	Z0 = log(Yp + 1); Z0[is.na(Z0)] = 0; Z = Z0
	for (s in 1:30) Z = lame:::rZ_pois_fc(Z, matrix(1, n, n), 0, 1, Yp)
	expect_gt(mean((Z != Z0)[!is.na(Yp)]), 0.9)
})

test_that("poisson fitted() and predict() return the stored mean counts silently", {
	S = 6L; n = 3L; rn = letters[seq_len(n)]
	x = array(1, c(n, n, 1), dimnames = list(rn, rn, "intercept"))
	y = matrix(2, n, n, dimnames = list(rn, rn)); diag(y) = NA_real_
	f = list(
		BETA = cbind(intercept = seq(-.2, .2, length.out = S)),
		VC = cbind(va = rep(0, S), cab = 0, vb = 0, rho = 0, ve = .6),
		family = "poisson",
		mode = "unipartite",
		R = 0L,
		rvar = FALSE,
		cvar = FALSE,
		symmetric = FALSE,
		intercept = TRUE,
		nA = n,
		nB = n,
		APM = setNames(numeric(n), rn),
		BPM = setNames(numeric(n), rn),
		X = x,
		YPM = y
	)
	class(f) = "ame"
	# no message, warning or note, on the first call or any later one
	expect_silent({ a = fitted(f) })
	expect_silent({ b = predict(f) })
	expect_identical(a, y); expect_identical(b, y)
})

test_that("bipartite ranks are bounded by their own mode", {
	set.seed(7)
	Yb = matrix(rnorm(12 * 4), 12, 4, dimnames = list(paste0("r", 1:12), paste0("c", 1:4)))
	fit = suppressWarnings(suppressMessages(ame(
		Yb,
		mode = "bipartite",
		family = "normal",
		R_row = 5,
		R_col = 2,
		nscan = 10,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE
	)))
	expect_identical(dim(fit$U), c(12L, 5L))
	expect_identical(dim(fit$V), c(4L, 2L))
	expect_error(
		suppressWarnings(suppressMessages(ame(
			Yb,
			mode = "bipartite",
			R_row = 2,
			R_col = 4,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		))),
		"R_col"
	)
	expect_error(
		suppressWarnings(suppressMessages(ame(
			Yb,
			mode = "bipartite",
			R_row = 12,
			R_col = 1,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			gof = FALSE
		))),
		"R_row"
	)
	fit_l = suppressWarnings(suppressMessages(lame(
		list(Yb, Yb),
		mode = "bipartite",
		family = "normal",
		R_row = 5,
		R_col = 2,
		nscan = 10,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_true(all(is.finite(fit_l$BETA)))
	expect_error(
		suppressWarnings(suppressMessages(lame(
			list(Yb, Yb),
			mode = "bipartite",
			R_row = 2,
			R_col = 4,
			nscan = 10,
			burn = 2,
			odens = 1,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		))),
		"R_col"
	)
})

test_that("a resumed symmetric chain continues the reciprocity ramp", {
	skip_on_cran()
	set.seed(8); n = 8
	Ys = replicate(2, {
		y = matrix(rnorm(n * n), n, n); y = (y + t(y)) / 2; diag(y) = NA
		dimnames(y) = list(letters[1:n], letters[1:n]); y
	}, simplify = FALSE)
	ck = withr::local_tempfile(fileext = ".rds")
	f1 = suppressWarnings(suppressMessages(lame(
		Ys,
		family = "normal",
		symmetric = TRUE,
		R = 0,
		nscan = 40,
		burn = 20,
		odens = 1,
		checkpoint_path = ck,
		checkpoint_every = 50L,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(f1$start_vals$rho, 1 - 1 / sqrt(60))
	expect_identical(f1$start_vals$sym_ramp_sweeps, 60)
	# checkpoint after 49 sweeps plus four more makes 53; a ramp restarted from 1 would read 1 - 1 / sqrt(4) = 0.5
	f2 = suppressWarnings(suppressMessages(lame_resume(ck, nscan_more = 4)))
	expect_equal(f2$start_vals$rho, 1 - 1 / sqrt(53))
	expect_identical(f2$start_vals$sym_ramp_sweeps, 53)
})
