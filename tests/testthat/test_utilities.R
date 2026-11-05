skip_on_cran()

# utility functions
library(lame)

test_that("get_start_vals returns appropriate starting values", {
	set.seed(6886)
	n = 20
	Tn = 1 # single time point

	# binary network as a 3-d array
	Y_bin_2d = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y_bin_2d) = NA
	Y_bin = array(Y_bin_2d, dim = c(n, n, Tn))

	starts_bin = get_start_vals(
		start_vals = NULL,
		Y = Y_bin,
		family = "binary",
		xP = 0, # no exogenous covariates
		rvar = TRUE,
		cvar = TRUE,
		R = 2 # number of latent dimensions
	)
	expect_true(is.list(starts_bin))
	expect_true("Z" %in% names(starts_bin))
	expect_equal(dim(starts_bin$Z), dim(Y_bin))

	# normal network as a 3-d array
	Y_norm_2d = matrix(rnorm(n * n), n, n)
	diag(Y_norm_2d) = NA
	Y_norm = array(Y_norm_2d, dim = c(n, n, Tn))

	starts_norm = get_start_vals(
		start_vals = NULL,
		Y = Y_norm,
		family = "normal",
		xP = 0,
		rvar = TRUE,
		cvar = TRUE,
		R = 2
	)
	# for the normal family Z starts at Y with the diagonal filled, so compare the off-diagonal cells
	expect_equal(dim(starts_norm$Z), dim(Y_norm))
	for (t in 1:Tn) {
		non_diag = !diag(n)
		expect_equal(starts_norm$Z[, , t][non_diag], Y_norm[, , t][non_diag])
	}

	# the other components
	expect_true(is.numeric(starts_bin$beta))
	expect_true(is.numeric(starts_bin$s2))
	expect_true(is.numeric(starts_bin$rho))
	expect_true(is.matrix(starts_bin$U))
	expect_true(is.matrix(starts_bin$V))
	expect_true(is.matrix(starts_bin$Sab))
})

test_that("mhalf function works correctly", {
	# a positive definite matrix
	A = matrix(c(4, 2, 2, 3), 2, 2)
	A_half = mhalf(A)

	# A_half %*% A_half equals A
	A_reconstructed = A_half %*% A_half
	expect_equal(A_reconstructed, A, tolerance = 1e-10)

	# the identity matrix
	I = diag(3)
	I_half = mhalf(I)
	expect_equal(I_half, I)

	# a singular matrix is handled
	B = matrix(c(1, 1, 1, 1), 2, 2)
	expect_no_error(mhalf(B))
})

test_that("rSab_fc generates valid covariance samples", {
	set.seed(6886)

	n = 20
	a = rnorm(n, 0, 1)
	b = rnorm(n, 0, 1)
	Sab0 = diag(2)
	eta0 = 4 # at least 4 for a proper prior

	# a draw with the default parameters
	Sab = rSab_fc(a, b, Sab0, eta0)

	# a 2 x 2 matrix
	expect_equal(dim(Sab), c(2, 2))

	# symmetric
	expect_equal(Sab[1, 2], Sab[2, 1])

	# positive definite
	eigenvals = eigen(Sab)$values
	expect_true(all(eigenvals > 0))
})

test_that("check_format validates input correctly", {
	n = 20
	Tn = 3

	# valid list of matrices
	Y_list = list()
	for (t in 1:Tn) {
		Y_t = matrix(rnorm(n * n), n, n)
		diag(Y_t) = NA
		rownames(Y_t) = colnames(Y_t) = paste0("Node", 1:n)
		Y_list[[t]] = Y_t
	}

	# check_format() returns nothing; it errors only on a bad format
	expect_no_error(check_format(Y_list))

	# check_format() expects a list and errors on a single matrix
	Y_single = matrix(rnorm(n * n), n, n)
	diag(Y_single) = NA
	expect_error(check_format(Y_single))

	# an array is refused as well
	Y_array = array(rnorm(n * n * Tn), dim = c(n, n, Tn))
	expect_error(check_format(Y_array))
})

test_that("zscores function computes correct scores", {
	# a simple vector
	y = c(3, 1, 4, 1, 5, 9, 2, 6)
	z = zscores(y)

	# rank order is preserved
	expect_equal(order(y), order(z))

	# ties
	y_ties = c(1, 2, 2, 3)
	z_ties = zscores(y_ties)
	expect_equal(z_ties[2], z_ties[3]) # tied values share a z-score

	# NA values
	y_na = c(1, NA, 3, 4)
	z_na = zscores(y_na)
	expect_true(is.na(z_na[2]))
	expect_false(any(is.na(z_na[-2])))
})

test_that("plotting functions run without error", {
	skip_on_cran()

	set.seed(6886)
	n = 15

	# a simple network and its fit
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA
	# row names so APM and BPM are named
	rownames(Y) = colnames(Y) = paste0("Actor", 1:n)

	fit = ame(
		Y,
		R = 2,
		family = "binary",
		burn = 50,
		nscan = 200,
		verbose = FALSE,
		gof = TRUE
	)

	# the plotting functions run without error
	expect_no_error({
		pdf(NULL) # null device, nothing is drawn

		# ab_plot, if it accepts ame objects
		if (exists("ab_plot")) {
			tryCatch(
				ab_plot(fit, effect = "sender"),
				error = function(e) {
					if (grepl("class 'ame' or 'lame'", e$message)) {
						# ab_plot may accept lame objects only
						NULL
					} else {
						stop(e)
					}
				}
			)
		}

		# uv_plot
		if (exists("uv_plot")) {
			tryCatch(
				uv_plot(fit),
				error = function(e) {
					if (grepl("class 'ame' or 'lame'", e$message)) {
						NULL
					} else {
						stop(e)
					}
				}
			)
		}

		# gof_plot
		if (exists("gof_plot")) {
			tryCatch(
				gof_plot(fit),
				error = function(e) {
					if (grepl("class 'ame' or 'lame'", e$message)) {
						NULL
					} else {
						stop(e)
					}
				}
			)
		}

		# trace_plot
		if (exists("trace_plot")) {
			tryCatch(
				trace_plot(list(BETA = fit$BETA)),
				error = function(e) {
					if (grepl("class 'ame' or 'lame'", e$message)) {
						NULL
					} else {
						stop(e)
					}
				}
			)
		}

		dev.off()
	})
})

test_that("summary and print methods work", {
	skip_on_cran()

	set.seed(6886)
	n = 15

	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA

	fit = ame(
		Y,
		R = 1,
		family = "binary",
		burn = 50,
		nscan = 200,
		verbose = FALSE
	)

	expect_no_error({
		s = summary(fit)
	})

	expect_no_error({
		capture.output(print(fit))
	})
})

test_that("effective sample size calculation works", {
	skip_on_cran()
	skip_if_not(
		requireNamespace("coda", quietly = TRUE),
		message = "coda package not available"
	)

	set.seed(6886)

	# an AR(1) sequence of moderate correlation
	n = 1000
	rho = 0.5
	x = numeric(n)
	x[1] = rnorm(1)
	for (i in 2:n) {
		x[i] = rho * x[i - 1] + rnorm(1, 0, sqrt(1 - rho^2))
	}

	ess = coda::effectiveSize(x)

	# autocorrelation puts the ESS below n
	expect_lt(ess, n)
	# and well above n / 5 at rho = 0.5
	expect_gt(ess, n / 5)
})

test_that("array_to_list conversion works", {
	n = 10
	Tn = 4

	# an array with actor names
	Y_array = array(rnorm(n * n * Tn), dim = c(n, n, Tn))
	rownames(Y_array) = colnames(Y_array) = paste0("Actor", 1:n)

	# every actor present at every time
	actorList = lapply(1:Tn, function(t) paste0("Actor", 1:n))

	# array_to_list() takes all three arguments
	Y_list = array_to_list(Y_array, actorList, paste0("T", 1:Tn))

	expect_true(is.list(Y_list))
	expect_equal(length(Y_list), Tn)

	# each element is a matrix
	for (t in 1:Tn) {
		expect_true(is.matrix(Y_list[[t]]))
		expect_equal(dim(Y_list[[t]]), c(n, n))
		# array_to_list() sets the diagonal to NA
		Y_expected = Y_array[, , t]
		diag(Y_expected) = NA
		expect_equal(Y_list[[t]], Y_expected)
	}
})

test_that("sampsonmonks liking layers share the roster actor order (regression)", {
	# the like layers follow the roster order; these nominee sets match ergm::samplk1/2/3 cell for cell
	data("sampsonmonks", package = "lame")
	S = sampsonmonks
	expect_equal(dimnames(S)[[1]], dimnames(S)[[2]])
	nominees = function(layer, monk) sort(names(which(S[monk, , layer] > 0)))
	expect_equal(nominees("like_m2", "ROMUL"), sort(c("PETER", "VICTOR", "HUGH")))
	expect_equal(nominees("like_m2", "JOHN"), sort(c("BONAVEN", "HUGH", "BASIL")))
	expect_equal(nominees("like_m1", "ROMUL"), sort(c("BONAVEN", "PETER", "AMAND")))
	expect_equal(nominees("like_m1", "BASIL"), sort(c("JOHN", "ELIAS", "SIMP")))
	expect_equal(nominees("like", "ROMUL"), sort(c("BONAVEN", "AMBROSE", "PETER", "AMAND")))
	expect_equal(nominees("like", "JOHN"), sort(c("WINF", "HUGH", "BASIL")))
	# liking agrees with esteem well above chance
	off = row(S[, , 1]) != col(S[, , 1])
	agree = mean(((S[, , "like"] > 0) == (S[, , "esteem"] > 0))[off])
	expect_gt(agree, 0.8)
})

test_that("get_start_vals passes dynamic hyperparameter scalars through (regression)", {
	# get_start_vals() passes the dynamic scalars of als_start_vals() through, so a dynamic MCMC run can warm-start from an ALS fit
	set.seed(1)
	Y = array(rnorm(6 * 6 * 2), c(6, 6, 2))
	for (t in 1:2) diag(Y[, , t]) = NA
	sv = list(rho_uv = 0.55, sigma_uv = 0.4, rho_ab = 0.65, sigma_ab = 0.3)
	out = get_start_vals(
		sv,
		Y,
		family = "normal",
		xP = 0,
		rvar = TRUE,
		cvar = TRUE,
		R = 2
	)
	expect_equal(out$rho_uv, 0.55)
	expect_equal(out$sigma_uv, 0.4)
	expect_equal(out$rho_ab, 0.65)
	expect_equal(out$sigma_ab, 0.3)
	# absent scalars stay absent rather than appearing as NULL-filled slots
	out2 = get_start_vals(
		list(rho_uv = 0.5),
		Y,
		family = "normal",
		xP = 0,
		rvar = TRUE,
		cvar = TRUE,
		R = 2
	)
	expect_null(out2$sigma_uv)
	# and junk is rejected loudly
	expect_error(
		get_start_vals(
			list(rho_uv = "high"),
			Y,
			family = "normal",
			xP = 0,
			rvar = TRUE,
			cvar = TRUE,
			R = 2
		),
		"finite"
	)
})

test_that("get_start_vals takes a partial warm start without cross-slot matching", {
	set.seed(6886)
	n = 8
	Y = array(rnorm(n * n * 2), c(n, n, 2))
	for (t in 1:2) diag(Y[, , t]) = NA
	out = get_start_vals(
		list(beta = c(0.3, 0.1)),
		Y,
		family = "normal",
		xP = 2,
		rvar = TRUE,
		cvar = TRUE,
		R = 0
	)
	expect_equal(out$beta, c(0.3, 0.1))
	# b keeps its own default (one value per actor), not the beta vector
	expect_length(out$b, n)
	expect_length(out$a, n)
	# rho stays scalar when only the dynamic hyperparameters are supplied
	out2 = get_start_vals(
		list(rho_uv = 0.5, rho_ab = 0.7),
		Y,
		family = "normal",
		xP = 0,
		rvar = TRUE,
		cvar = TRUE,
		R = 0
	)
	expect_length(out2$rho, 1L)
	expect_equal(out2$rho_uv, 0.5)
})

test_that("the directed effects kernel reproduces the exact symmetric posterior with a tied prior", {
	set.seed(3)
	n = 20
	v = 0.3
	a = rnorm(n, 0, sqrt(v))
	E = matrix(0, n, n)
	E[upper.tri(E)] = rnorm(n * (n - 1) / 2)
	E = E + t(E)
	mu = -0.2
	Z = mu + outer(a, a, "+") + E
	diag(Z) = NA
	# exact posterior of (mu, a) under z_ij = mu + a_i + a_j + e, one observation per unordered pair, a ~ N(0, v)
	ut = which(upper.tri(Z), arr.ind = TRUE)
	D = matrix(0, nrow(ut), n + 1)
	D[, 1] = 1
	for (r in seq_len(nrow(ut))) {
		D[r, 1 + ut[r, 1]] = 1
		D[r, 1 + ut[r, 2]] = 1
	}
	P = crossprod(D) + diag(c(1e-8, rep(1 / v, n)))
	Vp = solve(P)
	mp = drop(Vp %*% crossprod(D, Z[upper.tri(Z)]))
	exact_a2 = mean(mp[-1]^2 + diag(Vp)[-1])
	X = precomputeX(design_array(intercept = TRUE, n = n))
	draw_a2 = function(Sab) {
		d = replicate(300, {
			Zd = Z
			diag(Zd) = mu + 2 * a + rnorm(n, 0, sqrt(1.9999))
			tmp = rbeta_ab_fc(Zd, Sab, 0.9999, X, 1)
			(tmp$a + tmp$b) / 2
		})
		mean(d^2)
	}
	# the tied (rank-one) prior the symmetric fits use matches the exact
	# posterior; the diagonal prior under-disperses the averaged effects
	expect_lt(abs(draw_a2(matrix(v, 2, 2)) - exact_a2), 0.03)
	expect_lt(draw_a2(diag(2) * v), 0.85 * exact_a2)
})

# the binary start reproduces the observed tie rate given the starting effects
test_that("a binary fit starts at the level of its tie rate", {
	set.seed(6886)
	n = 30; TT = 3
	Y = array(rbinom(n * n * TT, 1, 0.02), c(n, n, TT))
	for (t in seq_len(TT)) diag(Y[, , t]) = NA
	sv = lame:::get_start_vals(NULL, Y, "binary", xP = 2, rvar = TRUE, cvar = TRUE, R = 0, odmax = NULL, bip = FALSE)
	p = mean(Y > 0, na.rm = TRUE)
	implied = mean(pnorm(sv$beta[1] + outer(sv$a, sv$b, "+")), na.rm = TRUE)
	expect_lt(abs(implied - p) / p, 0.05)
})

# above the threshold, and on a family whose cells carry counts, nothing changes
test_that("a denser binary fit and a sparse poisson fit keep their effect variances", {
	set.seed(6886)
	n = 25; TT = 3
	nm = paste0("a", seq_len(n))
	Yd = lapply(seq_len(TT), function(t) {
		m = matrix(rbinom(n * n, 1, 0.2), n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	fd = suppressMessages(suppressWarnings(lame(
		Yd,
		family = "binary",
		R = 0,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_gt(fd$events_per_actor_period, 1)
	expect_false(lame:::.lame_ab_var_unidentified(fd))
	sd_ = summary(fd)
	expect_false(any(is.na(sd_$variance[c("va", "cab", "vb"), ])))
	expect_false(isTRUE(sd_$ab_var_unidentified))
	expect_false(any(is.na(confint(fd, parm = "va"))))
	# a poisson cell carries a count, so the same rate is not the same
	# information and the variances are still reported
	Yp = lapply(seq_len(TT), function(t) {
		m = matrix(rpois(n * n, 0.5 / (n - 1)) * 4L, n, n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	fp = suppressMessages(suppressWarnings(lame(
		Yp,
		family = "poisson",
		R = 0,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_lt(fp$events_per_actor_period, 1)
	expect_false(lame:::.lame_ab_var_unidentified(fp))
	expect_false(any(is.na(summary(fp)$variance[c("va", "vb"), ])))
	# a sparse bipartite panel keeps them as well: its row and column actors see
	# different numbers of cells, so one rate does not describe both sides
	nA = 40; nB = 10
	Yb = lapply(seq_len(TT), function(t) {
		m = matrix(
			rbinom(nA * nB, 1, 0.9 / nB),
			nA,
			nB,
			dimnames = list(paste0("r", seq_len(nA)), paste0("c", seq_len(nB)))
		)
		m
	})
	fb = suppressMessages(suppressWarnings(lame(
		Yb,
		family = "binary",
		mode = "bipartite",
		R_row = 0,
		R_col = 0,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_lt(fb$events_per_actor_period, 1)
	expect_false(lame:::.lame_ab_var_unidentified(fb))
	expect_false(any(is.na(summary(fb)$variance[c("va", "vb"), ])))
	# a fit without the recorded rate reports its variances
	old = fd; old$events_per_actor_period = NULL
	expect_false(lame:::.lame_ab_var_unidentified(old))
})

# a cbin or frn nomination cap sitting at (almost) every row's outdegree erases
# the cross-row spread the additive effects rest on
.sat_cbin_panel = function(n = 20, TT = 3, k = 3, seed = 6886) {
	set.seed(seed)
	nm = paste0("a", seq_len(n))
	lapply(seq_len(TT), function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) {
			m[i, sample(setdiff(seq_len(n), i), k)] = 1
		}
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
}

test_that("a saturating odmax in a multi-period fit reports every parameter", {
	Y = .sat_cbin_panel()
	msgs = character()
	fit = withCallingHandlers(
		suppressWarnings(lame(
			Y,
			family = "cbin",
			odmax = 3,
			R = 0,
			nscan = 200,
			burn = 100,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	# every row nominates exactly the cap
	expect_equal(fit$odmax_saturated_share, 1)
	expect_gt(fit$events_per_actor_period, 1)
	expect_false(lame:::.lame_ab_var_unidentified(fit))
	expect_true(is.na(lame:::.lame_ab_coef_reason(fit)))
	expect_false(grepl("not estimated", gsub("\\s+", " ", paste(msgs, collapse = " "))))
	s = summary(fit)
	expect_false(any(is.na(s$variance[c("va", "vb", "cab", "ve"), ])))
	expect_false(any(is.na(s$coefficients["intercept", ])))
	expect_equal(s$odmax_saturated_share, 1)
	expect_false(any(is.na(confint(fit, parm = c("intercept", "va", "vb", "ve")))))
	td = tidy(fit)
	expect_false(is.na(td$estimate[td$term == "intercept"]))
	expect_false(grepl("not estimated", gsub("\\s+", " ", paste(capture.output(print(s)), collapse = " "))))
})

test_that("hidden cells do not hide a saturated cap", {
	# every row nominates the full cap; a fifth of the cells is then hidden,
	# nominations included, so the observed outdegrees fall below the cap
	Y = .sat_cbin_panel()
	set.seed(9)
	Ym = lapply(Y, function(m) {
		off = which(!is.na(m))
		m[sample(off, round(0.2 * length(off)))] = NA
		m
	})
	full = lame:::.lame_odmax_saturated_share(Y, 3, "cbin")
	hid = lame:::.lame_odmax_saturated_share(Ym, 3, "cbin")
	expect_equal(full, 1)
	# the share is estimated from the thinned outdegrees: it reaches the
	# level at which a one-period design withholds its intercept
	expect_gte(as.numeric(hid), lame:::.lame_cap_level_share_hidden)
	one = list(family = "cbin", mode = "unipartite", odmax_saturated_share = hid, cells_hidden = TRUE, n_time = 1)
	expect_identical(lame:::.lame_ab_coef_reason(one), "saturated")
	expect_true(is.na(lame:::.lame_ab_coef_reason(modifyList(one, list(cells_hidden = FALSE)))))
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(one)))
	expect_match(lame:::.lame_sat_phrase(hid), "nominate the full odmax")
	# a cap far above the outdegrees stays unsaturated with the same holes
	far = lame:::.lame_odmax_saturated_share(Ym, 15, "cbin")
	expect_lt(far, 0.2)
	expect_false(isTRUE(attr(far, "not_ruled_out")))
	# outdegrees spread below a small cap are not read as saturated
	set.seed(12)
	Ys = lapply(1:4, function(t) {
		m = matrix(0, 30, 30)
		for (i in 1:30) m[i, sample(setdiff(1:30, i), sample(1:5, 1))] = 1
		diag(m) = NA
		off = which(!is.na(m)); m[sample(off, round(0.2 * length(off)))] = NA
		m
	})
	spread = lame:::.lame_odmax_saturated_share(Ys, 5, "cbin")
	expect_lt(spread, 0.5)
	expect_false(isTRUE(attr(spread, "not_ruled_out")))
	# an absent actor (a whole row and column missing) is not a hidden cell
	Ya = lapply(Y, function(m) { m[1, ] = NA; m[, 1] = NA; m })
	expect_equal(
		lame:::.lame_odmax_saturated_share(Ya, 3, "cbin"),
		lame:::.lame_odmax_saturated_share(lapply(Y, function(m) m[-1, -1]), 3, "cbin")
	)
})

test_that("a cap above the typical outdegree reports everything", {
	set.seed(4)
	n = 20; TT = 3
	nm = paste0("a", seq_len(n))
	Y = lapply(seq_len(TT), function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) {
			k = sample(0:5, 1)
			if (k > 0) m[i, sample(setdiff(seq_len(n), i), k)] = 1
		}
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	fit = suppressMessages(suppressWarnings(lame(
		Y,
		family = "cbin",
		odmax = 9,
		R = 0,
		nscan = 200,
		burn = 100,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_lt(fit$odmax_saturated_share, lame:::.lame_odmax_sat_threshold)
	expect_false(lame:::.lame_ab_var_unidentified(fit))
	expect_false(any(is.na(summary(fit)$variance[c("va", "vb"), ])))
	expect_false(any(is.na(confint(fit, parm = "va"))))
	# a one-period fit at a saturated share withholds its intercept
	one = fit; one$n_time = 1; one$n.periods = 1; one$odmax_saturated_share = 1
	expect_identical(lame:::.lame_ab_coef_reason(one), "saturated")
	expect_false(lame:::.lame_ab_var_unidentified(one))
	# a fit without the recorded share reports its intercept
	one$odmax_saturated_share = NULL
	expect_true(is.na(lame:::.lame_ab_coef_reason(one)))
})

test_that("ame() withholds the intercept of a saturated cap and reports the variances", {
	Y = .sat_cbin_panel(n = 20, TT = 1, k = 3, seed = 7)[[1]]
	msgs = character()
	fit = withCallingHandlers(
		suppressWarnings(ame(
			Y,
			family = "cbin",
			odmax = 3,
			R = 0,
			nscan = 200,
			burn = 100,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)),
		message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
	)
	expect_equal(fit$odmax_saturated_share, 1)
	expect_false(lame:::.lame_ab_var_unidentified(fit))
	expect_identical(lame:::.lame_ab_coef_reason(fit), "saturated")
	expect_match(
		gsub("\\s+", " ", paste(msgs, collapse = " ")),
		"The intercept and any row-covariate coefficient are not estimated"
	)
	s = summary(fit)
	expect_false(any(is.na(s$variance[c("va", "vb", "cab", "ve"), c("Estimate", "StdError")])))
	expect_true(all(is.na(s$coefficients["intercept", ])))
	expect_match(gsub("\\s+", " ", paste(capture.output(print(s)), collapse = " ")), "not estimated")
	out = gsub("\\s+", " ", paste(capture.output(print(fit)), collapse = " "))
	expect_match(out, "not estimated")
	expect_true(all(is.na(confint(fit, parm = "intercept"))))
	expect_false(any(is.na(confint(fit, parm = c("va", "vb")))))
	expect_true(all(is.finite(fit$BETA[, "intercept"])))
})

test_that("the saturation share is computed only where a cap applies", {
	Y = .sat_cbin_panel(n = 12, TT = 2, k = 2, seed = 3)
	f = lame:::.lame_odmax_saturated_share
	expect_equal(f(Y, 2, "cbin"), 1)
	expect_equal(f(Y, 2, "frn"), 1)
	# a cap above every row's outdegree saturates nothing
	expect_equal(f(Y, 5, "cbin"), 0)
	# families without a cap carry no share at all
	expect_true(is.na(f(Y, 2, "binary")))
	expect_true(is.na(f(Y, 2, "normal")))
	# a missing cap resolves to the largest observed outdegree, as lame() does
	expect_equal(f(Y, NULL, "cbin"), 1)
	# a per-row cap is read row by row
	expect_equal(f(Y, c(rep(2, 6), rep(5, 6)), "cbin"), 0.5)
	# a 3-d array is the same panel as a list of slices
	arr = array(unlist(Y), c(12, 12, 2))
	expect_equal(f(arr, 2, "cbin"), 1)
})

test_that("the named constants set the saturated-cap setting and the reporting levels", {
	thr = lame:::.lame_odmax_sat_threshold
	expect_true(is.numeric(thr) && length(thr) == 1L && thr > 0.5 && thr < 1)
	expect_equal(lame:::.lame_cap_level_periods, 1L)
	expect_true(lame:::.lame_cap_level_share_hidden < lame:::.lame_cap_level_share)
	reason = function(s, ...) lame:::.lame_ab_coef_reason(list(
		family = "cbin",
		mode = "unipartite",
		odmax_saturated_share = s,
		n_time = 1,
		...
	))
	expect_identical(reason(1), "saturated")
	expect_true(is.na(reason(lame:::.lame_cap_level_share - 0.01)))
	expect_identical(reason(lame:::.lame_cap_level_share_hidden, cells_hidden = TRUE), "saturated")
	expect_true(is.na(reason(lame:::.lame_cap_level_share_hidden - 0.01, cells_hidden = TRUE)))
	expect_match(lame:::.lame_sat_phrase(0.97), "97 percent of the observed rows nominate the full odmax")
	# the share is computed on the observed outdegrees
	n = 20
	k_at = ceiling(thr * n - 1e-9)
	Y = lapply(1:2, function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, setdiff(seq_len(n), i)[seq_len(if (i <= k_at) 3 else 2)]] = 1
		diag(m) = NA; m
	})
	sh = lame:::.lame_odmax_saturated_share(Y, 3, "cbin")
	expect_equal(as.numeric(sh), k_at / n)
	# the pre-fit message of a one-period fit fires at a full cap and not under it
	msg_of = function(Yp) {
		msgs = character()
		withCallingHandlers(
			suppressWarnings(ame(
				Yp,
				family = "cbin",
				odmax = 3,
				R = 0,
				nscan = 20,
				burn = 10,
				odens = 1,
				verbose = FALSE,
				gof = FALSE
			)),
			message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
		)
		gsub("\\s+", " ", paste(msgs, collapse = " "))
	}
	nm = paste0("a", seq_len(n))
	full = .sat_cbin_panel(n = n, TT = 1, k = 3, seed = 5)[[1]]
	under = full; under[1:5, ] = 0; under[cbind(1:5, 6:10)] = 1; under[cbind(1:5, 11:15)] = 1; diag(under) = NA
	expect_match(msg_of(full), "nominate the full odmax")
	expect_false(grepl("nominate the full odmax", msg_of(under)))
	rd = tools::Rd_db("lame")
	for (f in c("lame.Rd", "ame.Rd")) {
		txt = gsub(
			"\\s+",
			" ",
			paste(
				capture.output(tools::Rd2txt(
					rd[[f]],
					options = list(underline_titles = FALSE)
				)),
				collapse = " "
			)
		)
		expect_match(txt, paste0(round(100 * thr), " percent of the observed rows"))
	}
})

# a panel of cbin nominations with a fifth of the cells hidden
.hidden_cbin_panel = function(n = 12, cap = 3, mu = -0.2, seed = 41, bip = FALSE) {
	set.seed(seed)
	nb = if (bip) 9 else n
	z = mu + outer(rnorm(n, 0, 0.6), rnorm(nb, 0, 0.4), "+") + matrix(rnorm(n * nb), n, nb)
	if (!bip) diag(z) = -Inf
	y = matrix(0, n, nb)
	for (i in seq_len(n)) {
		top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
		y[i, top[z[i, top] > 0]] = 1
	}
	if (!bip) diag(y) = NA
	off = which(!is.na(y))
	y[sample(off, round(0.2 * length(off)))] = NA
	y
}
# the nomination rule applied to a latent row set: the top min(cap, #positive)
# positive cells of each row
.nominate = function(Z, cap, sq) {
	if (sq) diag(Z) = -Inf
	t(apply(Z, 1, function(z) {
		k = min(cap, sum(z > 0)); y = numeric(length(z))
		if (k > 0) y[order(z, decreasing = TRUE)[seq_len(k)]] = 1
		y
	}))
}

test_that("cbin latent draws stay consistent with the observed cells when cells are hidden", {
	for (bip in c(FALSE, TRUE)) {
		Y = .hidden_cbin_panel(bip = bip)
		n = nrow(Y); nb = ncol(Y); cap = 3
		EZ = matrix(-0.2, n, nb)
		Z = matrix(rnorm(n * nb), n, nb)
		na_off = is.na(Y); if (!bip) diag(na_off) = FALSE
		odobs = rowSums(Y > 0, na.rm = TRUE)
		obs = !is.na(Y)
		ok = logical(200)
		for (s in 1:200) {
			room = odobs + rowSums(na_off & Z > 0)
			Z = if (bip) lame:::rZ_cbin_bip_fc(Z, EZ, Y, rep(cap, n), room)
			else lame:::rZ_cbin_fc(Z, EZ, 0.2, Y, rep(cap, n), room)
			ok[s] = all(.nominate(Z, cap, !bip)[obs] == Y[obs])
		}
		# after every sweep the nominations the latent values imply match
		# every observed cell
		expect_true(all(ok))
	}
})

test_that("a hidden cell of a row whose observed nominations fill the cap stays below them", {
	n = 10; cap = 3
	Y = matrix(0, n, n); diag(Y) = NA
	Y[1, 2:4] = 1; Y[1, 8:10] = NA
	EZ = matrix(0.5, n, n)
	Z = matrix(rnorm(n * n), n, n)
	set.seed(3)
	hid_max = nom_min = numeric(300)
	for (s in 1:300) {
		Z = lame:::rZ_cbin_fc(
			Z,
			EZ,
			0,
			Y,
			rep(cap, n),
			rowSums(Y > 0, na.rm = TRUE) + rowSums(is.na(Y) & row(Y) != col(Y) & Z > 0)
		)
		hid_max[s] = max(Z[1, 8:10]); nom_min[s] = min(Z[1, 2:4])
	}
	# the cap is spent on the observed nominations, so no hidden cell of the
	# row can outrank them
	expect_true(all(hid_max < nom_min))
})

test_that("complete cbin panels take the unchanged draw", {
	set.seed(5); n = 12
	Y = matrix(0, n, n); for (i in 1:n) Y[i, sample(setdiff(1:n, i), sample(1:4, 1))] = 1; diag(Y) = NA
	Z = matrix(rnorm(n * n), n, n); EZ = matrix(rnorm(n * n, -0.5), n, n)
	# the sweep first repairs a start outside the nomination constraints, so the reference starts from the repaired state
	Z = lame:::.lame_nom_feasible(Z, Y, NULL, rep(3, n), rowSums(Y > 0, na.rm = TRUE), FALSE)
	set.seed(9); z_new = lame:::rZ_cbin_fc(Z, EZ, 0.3, Y, rep(3, n), rowSums(Y > 0, na.rm = TRUE))
	set.seed(9); z_old = lame:::.rZ_cbin_sweep(Z, EZ, 0.3, Y, rep(3, n), rowSums(Y > 0, na.rm = TRUE))
	expect_identical(z_new, z_old)
})

test_that("the likelihood-ratio saturation flag needs an estimate near the threshold", {
	thr = lame:::.lame_odmax_sat_threshold
	# half the rows at a cap of 5, a fifth of the cells hidden: the estimate
	# sits well below the threshold, and the flag must not be set whatever
	# the likelihood-ratio test says
	n_checked = 0
	for (seed in 1:6) {
		set.seed(seed)
		Ys = lapply(1:4, function(t) {
			m = matrix(0, 30, 30)
			for (i in 1:30) m[i, sample(setdiff(1:30, i), sample(c(2, 3, 4, 5, 5, 5), 1))] = 1
			diag(m) = NA
			off = which(!is.na(m)); m[sample(off, round(0.2 * length(off)))] = NA
			m
		})
		sh = lame:::.lame_odmax_saturated_share(Ys, 5, "cbin")
		if (sh < thr - 0.08) {
			n_checked = n_checked + 1
			expect_false(isTRUE(attr(sh, "not_ruled_out")))
			expect_true(is.na(lame:::.lame_ab_var_unident_reason(list(
				family = "cbin",
				mode = "unipartite",
				odmax_saturated_share = sh
			))))
		}
	}
	expect_gt(n_checked, 3)
})

test_that("the compiled hidden-cell draw gives the R reference's draws", {
	for (bip in c(FALSE, TRUE)) {
		Y = .hidden_cbin_panel(bip = bip, seed = 43)
		n = nrow(Y); nb = ncol(Y)
		hid = lame:::.lame_nom_hidden_cells(Y, bip)
		set.seed(8); Z = matrix(rnorm(n * nb), n, nb); EZ = matrix(rnorm(n * nb, -0.3), n, nb)
		expect_identical(
			lame:::.lame_cbin_implied_labels(Z, Y, 3, hid),
			lame:::.lame_cbin_implied_labels_r(Z, Y, 3, hid)
		)
		rho = if (bip) 0 else 0.3
		set.seed(9); z_c = lame:::.lame_cbin_hidden_moves(Z, EZ, rho, Y, 3, hid, bip)
		set.seed(9); z_r = lame:::.lame_cbin_hidden_moves_r(Z, EZ, rho, Y, 3, hid, bip)
		expect_equal(z_c, z_r, tolerance = 1e-10)
	}
})

test_that("frn latent draws keep the observed order and membership when cells are hidden", {
	set.seed(44)
	n = 12; cap = 4
	z = -0.2 + outer(rnorm(n, 0, 0.6), rnorm(n, 0, 0.4), "+") + matrix(rnorm(n * n), n, n)
	diag(z) = -Inf
	Y = matrix(0, n, n)
	for (i in seq_len(n)) {
		top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
		top = top[z[i, top] > 0]
		if (length(top)) Y[i, top] = rank(z[i, top])
	}
	diag(Y) = NA
	off = which(!is.na(Y)); Y[sample(off, round(0.2 * length(off)))] = NA
	YL = t(apply(Y, 1, function(y) match(seq_len(cap), y)))
	EZ = matrix(-0.2, n, n); Z = matrix(rnorm(n * n), n, n)
	na_off = is.na(Y); diag(na_off) = FALSE
	obs = !is.na(Y)
	# from a random start the order constraints take a few sweeps to be met
	# (as with complete panels); from then on every sweep must meet them
	ok = logical(150)
	for (s in 1:160) {
		Z = lame:::rZ_frn_fc(Z, EZ, 0.2, Y, YL, rep(cap, n), rowSums(Y > 0, na.rm = TRUE) + rowSums(na_off & Z > 0))
		if (s <= 10) next
		nom = .nominate(Z, cap, TRUE)
		ord_ok = all(vapply(seq_len(n), function(i) {
			j = which(obs[i, ] & Y[i, ] > 0)
			length(j) < 2 || all(diff(Z[i, j[order(Y[i, j])]]) > 0)
		}, logical(1)))
		ok[s - 10] = all((nom[obs] > 0) == (Y[obs] > 0)) && ord_ok
	}
	# after every sweep the nominations the latent values imply match every
	# observed cell and the observed nominations keep their observed order
	expect_true(all(ok))
})

test_that("the posterior-predictive nomination draw lets hidden cells hold nominations", {
	Y = .hidden_cbin_panel()
	n = nrow(Y)
	Y[5, ] = NA; Y[, 5] = NA
	hid = lame:::.lame_nom_hidden_cells(Y, FALSE)
	YO = lame:::.lame_nom_predictive_yo(Y, FALSE)
	# hidden cells stay in contention; the diagonal and the absent actor do not
	expect_true(all(YO[hid] == 0))
	expect_true(all(is.na(diag(YO))) && all(is.na(YO[5, ])) && all(is.na(YO[, 5])))
	expect_identical(is.na(YO), !(!is.na(Y) | hid))
	EZ = matrix(0.3, n, n)
	set.seed(2); d1 = simY_frn(EZ, 0, 3, YO = YO)
	YO2 = matrix(0, n, n); YO2[5, ] = NA; YO2[, 5] = NA
	set.seed(2); d2 = simY_frn(EZ, 0, 3, YO = YO2)
	# the same draw as a panel with no hidden cell: every row fills its cap
	# over the cells a nomination can sit in
	expect_identical(d1, d2)
	expect_true(all(rowSums(d1 > 0, na.rm = TRUE)[-5] <= 3))
	expect_true(sum(d1[hid] > 0) > 0)
})

test_that("YPM of a cbin fit with hidden cells matches the observed nomination rate", {
	set.seed(61)
	n = 20; TT = 2; cap = 4
	a = rnorm(n, 0, 0.5); b = rnorm(n, 0, 0.5)
	Y = lapply(seq_len(TT), function(t) {
		z = 0.3 + outer(a, b, "+") + matrix(rnorm(n * n), n, n)
		diag(z) = -Inf
		y = matrix(0, n, n)
		for (i in seq_len(n)) {
			top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
			y[i, top[z[i, top] > 0]] = 1
		}
		diag(y) = NA
		off = which(!is.na(y)); y[sample(off, round(0.2 * length(off)))] = NA
		dimnames(y) = list(sprintf("a%02d", 1:n), sprintf("a%02d", 1:n))
		y
	})
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		family = "cbin",
		odmax = cap,
		R = 0,
		nscan = 400,
		burn = 100,
		odens = 4,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 3
	)))
	obs = lapply(Y, function(y) !is.na(y))
	yp = unlist(lapply(seq_len(TT), function(t) fit$YPM[[t]][obs[[t]]]))
	yo = unlist(lapply(seq_len(TT), function(t) Y[[t]][obs[[t]]]))
	# on the observed cells the predictive mean matches the observed mean
	expect_lt(abs(mean(yp) / mean(yo) - 1), 0.08)
	hid = lame:::.lame_nom_hidden_cells(Y[[1]], FALSE)
	expect_true(all(is.finite(fit$YPM[[1]][hid])))
	expect_true(all(is.na(diag(fit$YPM[[1]]))))
	# the settling check of hidden-cell fits is stored
	expect_named(fit$hidden_cell_ess, c("intercept", "va", "draws"))
})

# brute-force check of the shifts that keep a latent row consistent with its
# observed cells, over a grid of shifts
.shift_feasible = function(Z, Y, cap, i, deltas) {
	TT = dim(Z)[3]
	rows = lapply(seq_len(TT), function(t) {
		y = Y[, , t]
		hid = lame:::.lame_nom_hidden_cells(y, nrow(y) != ncol(y))
		elig = !is.na(y) | hid
		if (nrow(y) == ncol(y)) diag(elig) = FALSE
		o = !is.na(y[i, ]); if (nrow(y) == ncol(y)) o[i] = FALSE
		zi = Z[i, , t]; zi[!elig[i, ]] = -Inf
		list(z = zi, o = o, y = y[i, ])
	})
	vapply(deltas, function(d) {
		all(vapply(rows, function(r) {
			zi = r$z + d
			k = min(cap, sum(zi > 0)); nom = numeric(length(zi))
			if (k > 0) nom[order(zi, decreasing = TRUE)[seq_len(k)]] = 1
			all((nom[r$o] > 0) == (r$y[r$o] > 0))
		}, logical(1)))
	}, logical(1))
}

test_that("the row-shift bounds are exactly the shifts that keep a row consistent", {
	for (bip in c(FALSE, TRUE)) {
		Y1 = .hidden_cbin_panel(bip = bip, seed = 41); Y2 = .hidden_cbin_panel(bip = bip, seed = 42)
		if (!bip) { Y2[7, ] = NA; Y2[, 7] = NA }
		Y = array(c(Y1, Y2), c(dim(Y1), 2))
		n = nrow(Y1); nb = ncol(Y1); cap = 3
		Z = array(rnorm(n * nb * 2), c(n, nb, 2))
		# a consistent state from the latent draw
		for (s in 1:40) for (t in 1:2) {
			na_off = is.na(Y[, , t]); if (!bip) diag(na_off) = FALSE
			room = rowSums(Y[, , t] > 0, na.rm = TRUE) + rowSums(na_off & Z[, , t] > 0)
			Z[, , t] = if (bip) lame:::rZ_cbin_bip_fc(Z[, , t], matrix(-0.2, n, nb), Y[, , t], rep(cap, n), room)
			else lame:::rZ_cbin_fc(Z[, , t], matrix(-0.2, n, nb), 0, Y[, , t], rep(cap, n), room)
		}
		bd = lame:::.lame_nom_shift_bounds(Z, cap, lame:::.lame_nom_shift_masks(Y, bip))
		deltas = seq(-3, 3, by = 0.005)
		for (i in seq_len(n)) {
			f = .shift_feasible(Z, Y, cap, i, deltas)
			expect_true(f[deltas == 0])
			# feasible exactly on (lo, hi), up to the grid
			inside = deltas > bd$lo[i] + 0.005 & deltas < bd$hi[i] - 0.005
			outside = deltas < bd$lo[i] - 0.005 | deltas > bd$hi[i] + 0.005
			expect_true(all(f[inside]) && !any(f[outside]))
		}
	}
})

test_that("the row-shift move keeps the panel consistent and draws the row effect from its truncated prior", {
	Y1 = .hidden_cbin_panel(seed = 41); Y2 = .hidden_cbin_panel(seed = 43)
	Y = array(c(Y1, Y2), c(dim(Y1), 2))
	n = nrow(Y1); cap = 3
	masks = lame:::.lame_nom_shift_masks(Y, FALSE)
	Z = array(rnorm(n * n * 2), c(n, n, 2))
	a = rnorm(n, 0, 0.3); b = rnorm(n, 0, 0.3)
	Sab = matrix(c(0.3, 0.05, 0.05, 0.2), 2, 2)
	obs = !is.na(Y)
	for (t in 1:2) { o = obs[, , t]; diag(o) = FALSE; obs[, , t] = o }
	ok = logical(100)
	for (s in 1:140) {
		for (t in 1:2) {
			na_off = is.na(Y[, , t]); diag(na_off) = FALSE
			room = rowSums(Y[, , t] > 0, na.rm = TRUE) + rowSums(na_off & Z[, , t] > 0)
			Z[, , t] = lame:::rZ_cbin_fc(Z[, , t], matrix(-0.2, n, n) + a, 0, Y[, , t], rep(cap, n), room)
		}
		sh = lame:::.lame_nom_row_shift(Z, a, b, Sab, cap, masks, FALSE, rep(TRUE, n))
		Z = sh$Z; a = sh$a
		if (s > 40) ok[s - 40] = all(vapply(1:2, function(t) {
			all((.nominate(Z[, , t], cap, TRUE)[obs[, , t]] > 0) == (Y[, , t][obs[, , t]] > 0))
		}, logical(1)))
	}
	expect_true(all(ok))
	# from one state, many draws of the move: the row effect follows its
	# conditional prior truncated to (a + lo, a + hi)
	bd = lame:::.lame_nom_shift_bounds(Z, cap, masks)
	r = Sab[1, 2] / Sab[2, 2]; s = sqrt(Sab[1, 1] - r * Sab[1, 2])
	m = b * r
	draws = replicate(4000, lame:::.lame_nom_row_shift(Z, a, b, Sab, cap, masks, FALSE, rep(TRUE, n))$a)
	for (i in seq_len(n)) {
		lo = (a[i] + bd$lo[i] - m[i]) / s; hi = (a[i] + bd$hi[i] - m[i]) / s
		mu = m[i] + s * (dnorm(lo) - dnorm(hi)) / (pnorm(hi) - pnorm(lo))
		expect_lt(abs(mean(draws[i, ]) - mu), 4 * sd(draws[i, ]) / sqrt(4000) + 1e-8)
	}
})

test_that("simulate() leaves the actors absent from a period out of the nomination draw", {
	set.seed(71)
	n = 14; cap = 3
	nm = sprintf("a%02d", 1:n)
	Y = lapply(1:2, function(t) {
		keep = if (t == 1) nm else nm[-(1:3)]
		k = length(keep)
		z = 1.5 + matrix(rnorm(k * k), k, k); diag(z) = -Inf
		y = matrix(0, k, k)
		for (i in seq_len(k)) {
			top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
			y[i, top[z[i, top] > 0]] = 1
		}
		diag(y) = NA; dimnames(y) = list(keep, keep)
		y
	})
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		family = "cbin",
		odmax = cap,
		R = 0,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		seed = 5
	)))
	s = simulate(fit, nsim = 20, seed = 3)
	y2 = s$Y[[1]][[2]]
	expect_identical(dim(y2), c(11L, 11L))
	rs = unlist(lapply(s$Y, function(traj) rowSums(traj[[2]] > 0, na.rm = TRUE)))
	# the rows of the second period fill their cap over the actors present
	expect_gt(mean(rs == cap), 0.85)
})

test_that("a bipartite ame() cbin fit with hidden cells predicts the observed nomination rate", {
	set.seed(63)
	nA = 24; nB = 18; cap = 6
	a = rnorm(nA, 0, 0.5); b = rnorm(nB, 0, 0.5)
	z = 0.2 + outer(a, b, "+") + matrix(rnorm(nA * nB), nA, nB)
	Y = matrix(0, nA, nB)
	for (i in seq_len(nA)) {
		top = order(z[i, ], decreasing = TRUE)[seq_len(cap)]
		Y[i, top[z[i, top] > 0]] = 1
	}
	dimnames(Y) = list(sprintf("r%02d", 1:nA), sprintf("c%02d", 1:nB))
	Y[sample(length(Y), round(0.2 * length(Y)))] = NA
	fit = suppressWarnings(suppressMessages(ame(
		Y,
		family = "cbin",
		mode = "bipartite",
		odmax = cap,
		nscan = 400,
		burn = 100,
		odens = 4,
		verbose = FALSE,
		gof = TRUE,
		seed = 4
	)))
	obs = !is.na(Y)
	expect_lt(abs(mean(fit$YPM[obs]) / mean(Y[obs]) - 1), 0.08)
	expect_true(all(is.finite(fit$YPM[!obs])))
	expect_true(all(is.finite(fit$VC[, "va"])))
	expect_true(all(is.finite(fit$GOF)))
})

test_that("a hidden-cell cbin / frn fit records the effective sample size of its intercept and va", {
	set.seed(9)
	slow = cumsum(rnorm(200, 0, 0.05)); fast = rnorm(200)
	fake = list(
		family = "cbin",
		BETA = cbind(intercept = slow, x1_dyad = fast),
		VC = cbind(va = fast + 3, vb = fast + 3)
	)
	he = lame:::.lame_nom_hidden_ess(fake, TRUE)
	expect_lt(he[["intercept"]], 20)
	expect_gt(he[["va"]], 100)
	expect_equal(he[["draws"]], 200)
	expect_null(lame:::.lame_nom_hidden_ess(fake, FALSE))
	expect_null(lame:::.lame_nom_hidden_ess(modifyList(fake, list(family = "binary")), TRUE))
	# a fit with hidden cells carries the record and reports the intercept and va whatever the chain's effective size
	Y = .hidden_cbin_panel(n = 14, seed = 47)
	fit = suppressWarnings(suppressMessages(ame(
		Y,
		family = "cbin",
		odmax = 3,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	)))
	expect_named(fit$hidden_cell_ess, c("intercept", "va", "draws"))
	fit$hidden_cell_ess[["intercept"]] = 5
	expect_false(grepl("not usable", paste(capture.output(print(summary(fit))), collapse = " "), ignore.case = TRUE))
	expect_false(any(is.na(summary(fit)$variance["va", c("Estimate", "StdError")])))
	# a complete panel carries none
	Yc = Y; Yc[is.na(Yc)] = 0; diag(Yc) = NA
	fit_c = suppressWarnings(suppressMessages(ame(
		Yc,
		family = "cbin",
		odmax = 3,
		nscan = 200,
		burn = 50,
		odens = 5,
		verbose = FALSE,
		gof = FALSE,
		seed = 2
	)))
	expect_null(fit_c$hidden_cell_ess)
})

test_that("a single period is judged against its own saturated-cap setting", {
	thr1 = lame:::.lame_odmax_sat_threshold_single
	thr = lame:::.lame_odmax_sat_threshold
	expect_true(is.numeric(thr1) && length(thr1) == 1L && thr1 > 0.3 && thr1 < thr)
	expect_equal(lame:::.lame_odmax_threshold_for(1L), thr1)
	expect_equal(lame:::.lame_odmax_threshold_for(4L), thr)
	reason = function(s, cls = NULL, ...) lame:::.lame_ab_coef_reason(structure(
		list(
			family = "cbin",
			mode = "unipartite",
			odmax_saturated_share = s,
			...
		),
		class = cls
	))
	# an ame() fit has one period, a lame() fit the periods it carries, and an unclassed list the panel rule
	expect_identical(reason(1, "ame"), "saturated")
	expect_true(is.na(reason(1, c("lame", "ame"), n_time = 4L)))
	expect_identical(reason(1, c("lame", "ame"), n_time = 1L), "saturated")
	expect_true(is.na(reason(1)))
	# a slice with about seven in ten rows at the cap reports everything
	n = 20; k_at = 14
	nm = paste0("a", seq_len(n))
	mk = function(seed) {
		set.seed(seed)
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, sample(setdiff(seq_len(n), i), if (i <= k_at) 3 else 1)] = 1
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	}
	Y1 = mk(1)
	expect_equal(as.numeric(lame:::.lame_odmax_saturated_share(Y1, 3, "cbin")), k_at / n)
	fa = suppressMessages(suppressWarnings(ame(
		Y1,
		family = "cbin",
		odmax = 3,
		R = 0,
		nscan = 100,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		gof = FALSE
	)))
	expect_equal(fa$odmax_sat_threshold, thr1)
	expect_true(is.na(lame:::.lame_ab_coef_reason(fa)))
	expect_false(any(is.na(summary(fa)$variance["va", c("Estimate", "StdError")])))
	fl = suppressMessages(suppressWarnings(lame(
		list(Y1, mk(2)),
		family = "cbin",
		odmax = 3,
		R = 0,
		nscan = 100,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(fl$odmax_sat_threshold, thr)
	# the one-period rule does not apply to two periods
	expect_false(identical(lame:::.lame_ab_coef_reason(fl), "saturated"))
	f1 = suppressMessages(suppressWarnings(lame(
		list(Y1),
		family = "cbin",
		odmax = 3,
		R = 0,
		nscan = 100,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(f1$odmax_sat_threshold, thr1)
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(f1)))
	rd = tools::Rd_db("lame")
	for (f in c("lame.Rd", "ame.Rd")) {
		txt = gsub(
			"\\s+",
			" ",
			paste(
				capture.output(tools::Rd2txt(
					rd[[f]],
					options = list(underline_titles = FALSE)
				)),
				collapse = " "
			)
		)
		expect_match(txt, paste0(round(100 * thr1), " percent in a single period"))
	}
})

test_that("check_format refuses a NULL Xdyad period", {
	Y_list = lapply(1:2, function(t) matrix(rnorm(25), 5, 5))
	X_list = list(array(rnorm(25), c(5, 5, 1)), NULL)
	expect_error(check_format(Y_list, Xdyad = X_list), "Xdyad\\[\\[2\\]\\] is NULL")
	expect_no_error(check_format(Y_list, Xdyad = list(X_list[[1]], X_list[[1]])))
})

test_that("prediction_draws_long handles a fit with no coefficients", {
	set.seed(43)
	n = 8
	nm = paste0("a", seq_len(n))
	Yl = lapply(1:2, function(t) {
		y = matrix(rnorm(n * n), n, n, dimnames = list(nm, nm))
		diag(y) = NA
		y
	})
	fit = suppressWarnings(lame(
		Yl,
		family = "normal",
		R = 0,
		intercept = FALSE,
		burn = 5,
		nscan = 20,
		odens = 1,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		chains = 1
	))
	expect_equal(ncol(fit$BETA), 0L)
	pd = prediction_draws_long(fit, n_draws = 3, seed = 1)
	expect_equal(nrow(pd), 3L * 2L * n * n)
	# with no coefficients the linear predictor is the additive effects alone
	row1 = pd[pd$.draw == 1L & pd$period == 1L, ]
	expect_equal(row1$.value, fit$APM[row1$i] + fit$BPM[row1$j], ignore_attr = TRUE)
})
