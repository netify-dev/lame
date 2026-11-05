# dynamic_G parses, FALSE runs the default sampler, and TRUE stops on every route

test_that("dynamic_G = FALSE is byte-identical to default", {
	skip_on_cran()
	set.seed(2026)
	nA = 8; nB = 10; Tn = 3
	Y_list = lapply(seq_len(Tn), function(t) {
		mat = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
		rownames(mat) = paste0("a", seq_len(nA))
		colnames(mat) = paste0("b", seq_len(nB))
		mat
	})
	f1 = lame(
		Y_list,
		family = "binary",
		mode = "bipartite",
		R_row = 1,
		R_col = 1,
		nscan = 20,
		burn = 10,
		odens = 5,
		rvar = FALSE,
		cvar = FALSE,
		verbose = FALSE,
		seed = 7
	)
	f2 = lame(
		Y_list,
		family = "binary",
		mode = "bipartite",
		R_row = 1,
		R_col = 1,
		dynamic_G = FALSE,
		nscan = 20,
		burn = 10,
		odens = 5,
		rvar = FALSE,
		cvar = FALSE,
		verbose = FALSE,
		seed = 7
	)
	expect_identical(f1$BETA, f2$BETA)
	expect_identical(f1$YPM, f2$YPM)
	expect_null(f1$dynamic_G)
	expect_null(f2$dynamic_G)
})

test_that("dynamic_G = TRUE is refused on every route", {
	skip_on_cran()
	set.seed(2027)
	nA = 8; nB = 10; Tn = 4
	Y_list = lapply(seq_len(Tn), function(t) {
		mat = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
		rownames(mat) = paste0("a", seq_len(nA))
		colnames(mat) = paste0("b", seq_len(nB))
		mat
	})
	expect_error(
		lame(
			Y_list,
			family = "binary",
			mode = "bipartite",
			R_row = 2,
			R_col = 2,
			dynamic_G = TRUE,
			nscan = 30,
			burn = 15,
			odens = 5,
			rvar = FALSE,
			cvar = FALSE,
			verbose = FALSE,
			seed = 13
		),
		"not supported"
	)
	expect_error(
		lame(
			Y_list,
			family = "binary",
			mode = "bipartite",
			R_row = 2,
			R_col = 2,
			dynamic_G = TRUE,
			method = "als",
			verbose = FALSE,
			seed = 13
		),
		"not supported"
	)
	# the message tells moving positions (not identified) apart from static ones (not implemented)
	expect_error(
		lame(
			Y_list,
			family = "binary",
			mode = "bipartite",
			R_row = 2,
			R_col = 2,
			dynamic_G = TRUE,
			verbose = FALSE,
			seed = 13
		),
		"with static positions it is identified but not implemented"
	)
	# the argument stays in the signature so existing calls parse
	expect_true("dynamic_G" %in% names(formals(lame)))
	expect_true("dynamic_G" %in% names(formals(lame:::lame_dynamic_als)))
})

test_that("a rank-one dynamic_G fit is refused", {
	set.seed(1)
	nA = 14
	nB = 10
	TT = 3
	rn = sprintf("r%02d", seq_len(nA))
	cn = sprintf("c%02d", seq_len(nB))
	Y = lapply(seq_len(TT), function(t) {
		m = 1 * (matrix(rnorm(nA * nB), nA, nB) > 0)
		dimnames(m) = list(rn, cn)
		m
	})
	expect_error(
		lame(
			Y,
			family = "binary",
			mode = "bipartite",
			R_row = 1,
			R_col = 1,
			dynamic_uv = TRUE,
			dynamic_G = TRUE,
			nscan = 60,
			burn = 20,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		),
		"not supported"
	)
})
