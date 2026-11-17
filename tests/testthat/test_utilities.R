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
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(fd)))
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
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(fp)))
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
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(fb)))
	expect_false(any(is.na(summary(fb)$variance[c("va", "vb"), ])))
	# a fit without the recorded rate reports its variances
	old = fd; old$events_per_actor_period = NULL
	expect_true(is.na(lame:::.lame_ab_var_unident_reason(old)))
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
