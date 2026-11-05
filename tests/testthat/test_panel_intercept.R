# panel designs retain an intercept unless their observed covariates span it
skip_on_cran()

.panel_constant_design = function(C, n = 5L, m = n) {
	if (is.null(dim(C))) C = matrix(C, ncol = 1L)
	lapply(seq_len(nrow(C)), function(t) {
		a = array(rep(C[t, ], each = n * m), c(n, m, ncol(C)))
		dimnames(a) = list(NULL, NULL, paste0("x", seq_len(ncol(C))))
		a
	})
}

test_that("panel rank distinguishes reference indicators from a constant span", {
	Y = array(0, c(5, 5, 4))
	for (t in 1:4) diag(Y[, , t]) = NA
	cases = list(
		reference = list(C = rbind(0, diag(3)), redundant = FALSE),
		all_periods = list(C = diag(4), redundant = TRUE),
		continuous = list(C = c(-2, 0, 1, 4), redundant = FALSE),
		constant = list(C = rep(2, 4), redundant = TRUE),
		zero = list(C = rep(0, 4), redundant = FALSE),
		combination = list(C = cbind(1:4, 1 - (1:4)), redundant = TRUE),
		zero_plus_trend = list(C = cbind(0, 1:4), redundant = FALSE)
	)
	for (case in cases) {
		expect_identical(
			lame:::.lame_panel_intercept_redundant(
				.panel_constant_design(case$C),
				Y
			),
			case$redundant
		)
	}
	# an unobserved reference period supplies no identifying information
	Y[, , 1] = NA
	expect_true(lame:::.lame_panel_intercept_redundant(
		.panel_constant_design(rbind(0, diag(3))),
		Y
	))
})

test_that("panel rank matches an independent stacked design with missing cells", {
	set.seed(20261001)
	Y = array(0, c(5, 3, 4))
	Y[sample(length(Y), 10)] = NA
	for (redundant in c(FALSE, TRUE)) {
		X = replicate(4, array(rnorm(5 * 3 * 3), c(5, 3, 3)), simplify = FALSE)
		for (t in 1:4) {
			if (redundant) X[[t]][, , 3] = 1 - X[[t]][, , 1]
			X[[t]][, , 1] = X[[t]][, , 1] * 1e-8
			X[[t]][, , 2] = X[[t]][, , 2] * 1e8
		}
		X[[2]][1, 1, 2] = NA
		stacked = do.call(rbind, lapply(1:4, function(t) {
			M = matrix(X[[t]], ncol = 3)
			M[!is.na(c(Y[, , t])) & complete.cases(M), , drop = FALSE]
		}))
		stacked = sweep(stacked, 2, sqrt(colSums(stacked^2)), "/")
		expected = qr(cbind(stacked, 1))$rank == qr(stacked)$rank
		expect_identical(expected, redundant)
		expect_identical(lame:::.lame_panel_intercept_redundant(X, Y), expected)
	}
})

test_that("row, column and dyadic period covariates use the same panel decision", {
	n = 5L; TT = 4L
	Y = array(0, c(n, n, TT), dimnames = list(letters[1:n], letters[1:n], 1:TT))
	for (t in 1:TT) diag(Y[, , t]) = NA
	for (type in c("row", "col", "dyad")) {
		for (C in list(rbind(0, diag(3)), matrix(1:4, ncol = 1), matrix(0, 4, 1))) {
			p = ncol(C)
			node = array(
				rep(t(C), each = n),
				c(n, p, TT),
				dimnames = list(NULL, paste0("x", 1:p), NULL)
			)
			dyad = array(
				rep(t(C), each = n * n),
				c(n, n, p, TT),
				dimnames = list(NULL, NULL, paste0("x", 1:p), NULL)
			)
			args = list(
				Y = Y,
				Xrow = NULL,
				Xcol = NULL,
				Xdyad = NULL,
				actorSet = letters[1:n],
				intercept = TRUE,
				n = n,
				N = TT,
				pr = 0,
				pc = 0,
				pd = 0
			)
			args[[paste0("X", type)]] = if (type == "dyad") dyad else node
			args[[switch(type, row = "pr", col = "pc", dyad = "pd")]] = p
			design = do.call(get_design_rep, args)
			expect_true(design$intercept)
			expect_identical(dimnames(design$X)[[3]][1], "intercept")
			expect_equal(unname(design$X[1, 2, -1, ]), if (p == 1) c(C) else t(C))
			args$intercept = FALSE
			expect_false(do.call(get_design_rep, args)$intercept)
		}
	}
})

test_that("public panel fits keep reference-period intercepts in both modes", {
	for (bip in c(FALSE, TRUE)) {
		n = 8L; m = if (bip) 6L else n; TT = 4L
		set.seed(59)
		Y = lapply(1:TT, function(t) {
			M = matrix(
				rbinom(n * m, 1, .4),
				n,
				m,
				dimnames = list(paste0("r", 1:n), paste0(if (bip) "c" else "r", 1:m))
			)
			if (!bip) diag(M) = NA
			M
		})
		X = .panel_constant_design(rbind(0, diag(TT - 1L)), n, m)
		for (dynamic in c(FALSE, TRUE)) {
			fit = suppressMessages(lame(
				Y,
				Xdyad = X,
				family = "binary",
				R = 0,
				mode = if (bip) "bipartite" else "unipartite",
				dynamic_ab = dynamic,
				burn = 20,
				nscan = 40,
				odens = 5,
				seed = 59,
				gof = FALSE,
				plot = FALSE,
				verbose = FALSE
			))
			expect_identical(colnames(fit$BETA)[1], "intercept")
			expect_equal(ncol(fit$BETA), TT)
			expect_true(all(is.finite(fit$BETA)))
		}
		# all period indicators, in contrast, already span the constant
		X = .panel_constant_design(diag(TT), n, m)
		expect_warning({fit = suppressMessages(lame(
			Y,
			Xdyad = X,
			family = "binary",
			R = 0,
			mode = if (bip) "bipartite" else "unipartite",
			burn = 20,
			nscan = 40,
			odens = 5,
			seed = 59,
			gof = FALSE,
			plot = FALSE,
			verbose = FALSE
		))}, "collinear with the intercept")
		expect_false("intercept" %in% colnames(fit$BETA))
		expect_equal(ncol(fit$BETA), TT)
	}
})
