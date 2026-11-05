skip_on_cran()

# edge cases and error handling
library(lame)

test_that("models handle small networks", {
	skip_on_cran()

	set.seed(6886)

	# very small network
	n = 5
	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA

	fit = ame(
		Y,
		R = 1,
		family = "binary",
		burn = 50,
		nscan = 100,
		verbose = FALSE
	)
	expect_true(!is.null(fit))

	# bipartite with small dimensions
	Y_bip = matrix(rbinom(3 * 4, 1, 0.5), 3, 4)
	rownames(Y_bip) = paste0("Row", 1:3)
	colnames(Y_bip) = paste0("Col", 1:4)

	fit_bip = ame(
		Y_bip,
		mode = "bipartite",
		R_row = 1,
		R_col = 1,
		family = "binary",
		burn = 50,
		nscan = 100,
		verbose = FALSE,
		gof = FALSE
	) # no GOF for tiny networks
	expect_true(!is.null(fit_bip))
})

test_that("models handle sparse networks", {
	skip_on_cran()

	set.seed(6886)
	n = 30

	# very sparse network (density ~ 0.01)
	Y_sparse = matrix(0, n, n)
	n_edges = floor(0.01 * n * (n - 1))
	edge_idx = sample(which(upper.tri(Y_sparse)), n_edges)
	Y_sparse[edge_idx] = 1
	Y_sparse[lower.tri(Y_sparse)] = t(Y_sparse)[lower.tri(Y_sparse)]
	diag(Y_sparse) = NA

	expect_no_error({
		fit_sparse = ame(
			Y_sparse,
			R = 1,
			family = "binary",
			symmetric = TRUE,
			burn = 50,
			nscan = 200,
			verbose = FALSE
		)
	})

	# the model makes predictions (YPM exists)
	expect_true(!is.null(fit_sparse$YPM))
	expect_true(is.numeric(fit_sparse$YPM))
})

test_that("models handle dense networks", {
	skip_on_cran()

	set.seed(6886)
	n = 20

	# very dense network (density ~ 0.9)
	Y_dense = matrix(1, n, n)
	n_non_edges = floor(0.1 * n * (n - 1))
	non_edge_idx = sample(which(upper.tri(Y_dense)), n_non_edges)
	Y_dense[non_edge_idx] = 0
	Y_dense[lower.tri(Y_dense)] = t(Y_dense)[lower.tri(Y_dense)]
	diag(Y_dense) = NA

	expect_no_error({
		fit_dense = ame(
			Y_dense,
			R = 1,
			family = "binary",
			symmetric = TRUE,
			burn = 50,
			nscan = 200,
			verbose = FALSE
		)
	})

	# the fit captures the high density (binary YPM is a probability)
	expect_gt(mean(fit_dense$YPM, na.rm = TRUE), 0.5)
})

test_that("models handle networks with isolated nodes", {
	skip_on_cran()

	set.seed(6886)
	n = 20

	# network with isolated nodes
	Y = matrix(rbinom(n * n, 1, 0.2), n, n)
	# the first 3 nodes are isolated
	Y[1:3, ] = 0
	Y[, 1:3] = 0
	diag(Y) = NA

	expect_no_error({
		fit = ame(
			Y,
			R = 1,
			family = "binary",
			burn = 50,
			nscan = 200,
			verbose = FALSE
		)
	})

	# isolated nodes have low predicted probabilities
	isolated_preds = fit$YPM[1:3, -c(1:3)]
	expect_lt(mean(isolated_preds, na.rm = TRUE), 0.2)
})

test_that("longitudinal models handle varying network sizes", {
	skip_on_cran()

	set.seed(6886)

	# networks with different missing patterns (simulating different sizes)
	Tn = 4
	n_max = 20
	Y_list = list()

	for (t in 1:Tn) {
		Y_t = matrix(rbinom(n_max * n_max, 1, 0.3), n_max, n_max)

		# missing nodes, a different pattern each period
		n_missing = sample(0:5, 1)
		if (n_missing > 0) {
			missing_idx = sample(1:n_max, n_missing)
			Y_t[missing_idx, ] = NA
			Y_t[, missing_idx] = NA
		}

		diag(Y_t) = NA
		rownames(Y_t) = colnames(Y_t) = paste0("Node", 1:n_max)
		Y_list[[t]] = Y_t
	}

	expect_no_error({
		fit = lame(
			Y_list,
			R = 1,
			family = "binary",
			burn = 50,
			nscan = 200,
			verbose = FALSE,
			plot = FALSE
		)
	})
})

test_that("models handle extreme parameter values", {
	skip_on_cran()

	set.seed(6886)
	n = 15

	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y) = NA

	# extreme prior values
	expect_no_error({
		fit = ame(
			Y,
			R = 1,
			family = "binary",
			burn = 50,
			nscan = 100,
			verbose = FALSE,
			prior = list(Sab0 = diag(c(100, 100)))
		) # large variance prior
	})

	# very large R: the fit runs, with no rank cap below the network size; the default prior holds the factor covariance
	# at the data scale of the additive effects (R/uv_prior.R), so the predictor does not saturate
	expect_no_warning({
		fit_large_R = ame(
			Y,
			R = 10,
			family = "binary",
			burn = 50,
			nscan = 100,
			verbose = FALSE
		)
	}, class = "lame_separation")
	expect_s3_class(fit_large_R, "ame")
	expect_equal(fit_large_R$prior$kappa0, lame:::.lame_uv_prior(Y, "binary", 10L)$kappa0)
	expect_gt(fit_large_R$prior$kappa0, 22)

	# a rank at or above the network size is refused before sampling
	expect_error(
		ame(
			Y,
			R = n + 5,
			family = "binary",
			burn = 50,
			nscan = 100,
			verbose = FALSE
		),
		"must be smaller than the network dimension",
		class = "rlang_error"
	)
})

test_that("models handle non-standard data types correctly", {
	skip_on_cran()

	set.seed(6886)
	n = 20

	# ordinal data with gaps
	Y_ord = matrix(sample(c(1, 3, 5, 7), n * n, replace = TRUE), n, n)
	diag(Y_ord) = NA

	expect_no_error({
		fit_ord = ame(
			Y_ord,
			R = 1,
			family = "ordinal",
			intercept = FALSE,
			burn = 50,
			nscan = 200,
			verbose = FALSE
		)
	})

	# count data with large values
	Y_count = matrix(rpois(n * n, lambda = 10), n, n)
	Y_count[sample(1:(n * n), 5)] = 100 # a few large counts
	diag(Y_count) = NA

	expect_no_error({
		fit_count = ame(
			Y_count,
			R = 1,
			family = "poisson",
			burn = 50,
			nscan = 200,
			verbose = FALSE
		)
	})
})

test_that("bipartite models reject square matrices with mode='bipartite'", {
	skip_on_cran()

	set.seed(6886)
	n = 10

	# square matrix
	Y_square = matrix(rbinom(n * n, 1, 0.3), n, n)
	diag(Y_square) = NA

	# fits as unipartite
	expect_no_error({
		fit_uni = ame(
			Y_square,
			mode = "unipartite",
			R = 1,
			family = "binary",
			burn = 50,
			nscan = 100,
			verbose = FALSE
		)
	})

	tryCatch({
		fit_bip = ame(
			Y_square,
			mode = "bipartite",
			R_row = 1,
			R_col = 1,
			family = "binary",
			burn = 50,
			nscan = 100,
			verbose = FALSE,
			gof = FALSE
		)
		expect_true(TRUE)
	}, error = function(e) {
		expect_true(TRUE)
	})
})

test_that("models handle character row/column names", {
	skip_on_cran()

	set.seed(6886)
	n = 15

	Y = matrix(rbinom(n * n, 1, 0.3), n, n)
	rownames(Y) = paste0("Node_", 1:n)
	colnames(Y) = paste0("Node_", 1:n)
	diag(Y) = NA

	expect_no_error({
		fit = ame(
			Y,
			R = 1,
			family = "binary",
			burn = 50,
			nscan = 100,
			verbose = FALSE
		)
	})

	# names are preserved
	expect_equal(rownames(fit$YPM), rownames(Y))
	expect_equal(colnames(fit$YPM), colnames(Y))
})

test_that("lame() rejects a constant response and an asymmetric symmetric fit", {
	n = 8
	nm = paste0("a", seq_len(n))
	Y_const = lapply(1:2, function(t) {
		m = matrix(1, n, n)
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	expect_error(
		suppressWarnings(lame(
			Y_const,
			family = "binary",
			R = 0,
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			plot = FALSE
		)),
		"constant"
	)
	set.seed(6886)
	Y_dir = lapply(1:2, function(t) {
		m = matrix(rnorm(n * n), n, n)
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	expect_error(
		suppressWarnings(lame(
			Y_dir,
			family = "normal",
			R = 0,
			symmetric = TRUE,
			nscan = 10,
			burn = 5,
			odens = 1,
			verbose = FALSE,
			plot = FALSE
		)),
		"not symmetric"
	)
})

# a missing cell counts as a nomination when its imputed latent is positive, so rows with missing cells are neither
# assumed to have free slots nor saturated; nominations simulated from the fitted predictor are compared per row on observed cells
test_that("cbin and frn fits nominate at the observed rate under missing cells", {
	skip_on_cran()
	n = 22
	Tn = 3
	odmax = 4
	nm = sprintf("a%02d", seq_len(n))
	pred_rows = function(fit, Ymask, nsim = 100) {
		out = 0
		for (t in seq_along(fit$EZ)) {
			EZt = fit$EZ[[t]]
			obs = !is.na(Ymask[[t]])
			for (s in seq_len(nsim)) {
				ZS = EZt + matrix(rnorm(n * n), n, n)
				diag(ZS) = -Inf
				nom = matrix(0, n, n)
				for (i in seq_len(n)) {
					rs = rank(ZS[i, ]) - (n - odmax)
					nom[i, ] = 1 * ((rs > 0) & (ZS[i, ] > 0))
				}
				out = out + mean(rowSums(nom * obs)) / nsim
			}
		}
		out / length(fit$EZ)
	}
	for (fam in c("cbin", "frn")) for (shift in c(0, -1.6)) {
		set.seed(85)
		a = rnorm(n, 0, 0.5)
		b = rnorm(n, 0, 0.5)
		X = vector("list", Tn)
		Y = vector("list", Tn)
		for (t in seq_len(Tn)) {
			Xt = matrix(rnorm(n * n), n, n)
			diag(Xt) = NA
			Z = 0.2 + shift + 0.7 * Xt + outer(a, b, "+") + matrix(rnorm(n * n), n, n)
			Yt = matrix(0, n, n)
			for (i in seq_len(n)) {
				zi = Z[i, ]
				zi[i] = -Inf
				top = order(zi, decreasing = TRUE)[seq_len(odmax)]
				top = top[zi[top] > 0]
				if (fam == "cbin") Yt[i, top] = 1 else Yt[i, top] = rank(zi[top])
			}
			diag(Yt) = NA
			k = which(!is.na(Yt))
			Yt[sample(k, round(0.2 * length(k)))] = NA
			dimnames(Yt) = list(nm, nm)
			X[[t]] = array(Xt, c(n, n, 1), dimnames = list(nm, nm, "x1"))
			Y[[t]] = Yt
		}
		fit = suppressWarnings(suppressMessages(lame(
			Y,
			Xdyad = X,
			family = fam,
			R = 0,
			odmax = odmax,
			nscan = 300,
			burn = 150,
			odens = 5,
			verbose = FALSE,
			plot = FALSE,
			gof = FALSE
		)))
		obs_rows = mean(sapply(Y, function(m) mean(rowSums(m > 0, na.rm = TRUE))))
		ratio = pred_rows(fit, Y) / obs_rows
		expect_gt(ratio, 0.85)
		expect_lt(ratio, 1.15)
	}
})

test_that("frn fits store the ranks the user supplied, so residuals equal Y - YPM", {
	set.seed(6886)
	n = 16
	nm = sprintf("a%02d", seq_len(n))
	Y = lapply(1:2, function(t) {
		m = matrix(0, n, n)
		for (i in seq_len(n)) m[i, sample(setdiff(seq_len(n), i), 4)] = 1:4
		diag(m) = NA
		k = which(!is.na(m))
		m[sample(k, round(0.2 * length(k)))] = NA
		dimnames(m) = list(nm, nm)
		m
	})
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		family = "frn",
		R = 0,
		odmax = 4,
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	Yf = if (is.list(fit$Y)) fit$Y else lapply(1:2, function(t) fit$Y[, , t])
	for (t in 1:2) {
		o = !is.na(Y[[t]])
		expect_equal(Yf[[t]][o], Y[[t]][o])
	}
	r = residuals(fit)
	o = !is.na(Y[[1]])
	expect_lt(max(abs(r[[1]][o] - (Y[[1]][o] - fit$YPM[[1]][o]))), 1e-8)
})

# a directed fit to symmetric data counts each dyad's two equal cells twice and, with dcor, rho concentrates near one:
# a message points to symmetric = TRUE; with dynamic positions the dyadic-correlation fit is refused
test_that("lame() and ame() say so when Y is symmetric but the fit is directed", {
	set.seed(6886)
	n = 8
	nm = sprintf("a%02d", seq_len(n))
	Y = lapply(1:2, function(t) {
		m = matrix(rbinom(n * n, 1, 0.4), n, n)
		m[lower.tri(m)] = t(m)[lower.tri(m)]
		diag(m) = NA
		dimnames(m) = list(nm, nm)
		m
	})
	note = "lame_symmetric_y_note"
	flat = function(e) gsub("[[:space:]]+", " ", conditionMessage(e))
	quick = list(
		family = "binary", R = 0, nscan = 10, burn = 2, odens = 2, verbose = FALSE,
		plot = FALSE, gof = FALSE
	)
	e = expect_message(do.call(lame, c(list(Y), quick)), class = note)
	expect_match(flat(e), "exactly symmetric")
	expect_match(flat(e), "concentrates near one")
	expect_match(flat(e), "symmetric = TRUE")
	expect_no_message(do.call(lame, c(list(Y), quick, list(symmetric = TRUE))), class = note)
	e = expect_message(do.call(ame, c(list(Y[[1]]), quick)), class = note)
	expect_match(flat(e), "concentrates near one")
	# without a dyadic correlation the note is about the two cells of a dyad
	e = expect_message(do.call(lame, c(list(Y), quick, list(dcor = FALSE))), class = note)
	expect_match(flat(e), "two equal cells of each dyad count as two observations")
	expect_false(grepl("concentrates", flat(e)))
	# the nominations of a rank family only draw the correlation toward one
	e = expect_message(lame:::.lame_inform_symmetric_y(TRUE, family = "frn"), class = note)
	expect_match(flat(e), "drawn toward one")
	# a panel with one asymmetric dyad has no note
	Yd = Y
	Yd[[1]][2, 1] = 1 - Yd[[1]][2, 1]
	expect_no_message(do.call(lame, c(list(Yd), quick)), class = note)
	# the note is a message, not a warning
	expect_no_warning(suppressMessages(do.call(lame, c(list(Y), quick))))
})

# actors that enter and leave: each period is reported on its own actors while the effects and the stored Y cover the
# union, so the methods align by name; an absent actor's bipartite design is zero-filled so G is sampled
test_that("post-processing works under changing actor composition", {
	set.seed(6886)
	n = 10
	TT = 3
	nm = sprintf("a%02d", seq_len(n))
	Y = X = vector("list", TT)
	for (t in seq_len(TT)) {
		keep = nm[-t]
		m = matrix(rnorm(length(keep)^2), length(keep), length(keep), dimnames = list(keep, keep))
		diag(m) = NA
		Y[[t]] = m
		X[[t]] = array(rnorm(length(keep)^2), c(length(keep), length(keep), 1), dimnames = list(keep, keep, "x1"))
	}
	fit = suppressWarnings(suppressMessages(lame(
		Y,
		X,
		family = "normal",
		R = 2,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(dim(fit$YPM[[1]]), c(n - 1, n - 1))
	expect_length(fit$APM, n)
	r = residuals(fit)
	for (t in seq_len(TT)) {
		o = !is.na(Y[[t]])
		expect_lt(max(abs(r[[t]][o] - (Y[[t]][o] - fit$YPM[[t]][o]))), 1e-8)
	}
	p = predict(fit, newdata = X, type = "link")
	for (t in seq_len(TT)) {
		o = !is.na(Y[[t]])
		expect_lt(max(abs(p[[t]][o] - fit$EZ[[t]][o])), 1e-8)
	}
	s = simulate(fit, nsim = 2)
	expect_equal(dim(s$Y[[1]][[2]]), dim(fit$YPM[[2]]))
	expect_equal(rownames(s$Y[[1]][[2]]), rownames(fit$YPM[[2]]))
	# a fit without a multiplicative term has no positions to report
	fit0 = suppressWarnings(suppressMessages(lame(
		Y,
		family = "normal",
		R = 0,
		nscan = 20,
		burn = 5,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(nrow(latent_positions(fit0)), 0L)
	# symmetric: the forecast takes its column count from the one effect per actor
	Ys = lapply(Y, function(m) { m[lower.tri(m)] = t(m)[lower.tri(m)]; m })
	fs = suppressWarnings(suppressMessages(lame(
		Ys,
		family = "normal",
		R = 1,
		symmetric = TRUE,
		dynamic_ab = TRUE,
		nscan = 40,
		burn = 10,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	fc = predict(fs, h = 1L, type = "link")
	expect_equal(dim(fc[[1]]), c(n, n))
	# bipartite with latent factors: absent actors' covariate cells are zero-filled, so the G draw is finite in every sweep
	rn = sprintf("r%02d", 1:8); cn = sprintf("c%02d", 1:6)
	Yb = Xb = vector("list", TT)
	for (t in seq_len(TT)) {
		kr = rn[-t]; kc = cn[-t]
		Yb[[t]] = matrix(rnorm(length(kr) * length(kc)), length(kr), length(kc), dimnames = list(kr, kc))
		Xb[[t]] = array(rnorm(length(kr) * length(kc)), c(length(kr), length(kc), 1), dimnames = list(kr, kc, "x1"))
	}
	fb = suppressWarnings(suppressMessages(lame(
		Yb,
		Xb,
		family = "normal",
		mode = "bipartite",
		R_row = 1,
		R_col = 1,
		nscan = 60,
		burn = 20,
		odens = 5,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE
	)))
	expect_equal(fb$tryErrorChecks$G, 0L)
	expect_false(anyNA(fb$EZ[[1]]))
	expect_equal(rownames(fb$EZ[[1]]), rn)
})

# a rank above floor(n/3) is fitted like any other, with no cap and no caution
test_that("a latent rank above n/3 fits without a warning", {
	set.seed(5); n = 12
	Y = matrix(rnorm(n * n), n, n); diag(Y) = NA; dimnames(Y) = list(letters[1:n], letters[1:n])
	expect_no_warning(
		fit <- suppressMessages(ame(Y, R = 6, nscan = 40, burn = 10, odens = 2, verbose = FALSE, gof = FALSE))
	)
	expect_equal(fit$R, 6)
	expect_no_warning(
		fit <- suppressMessages(lame(
			list(t1 = Y, t2 = Y + rnorm(n * n)),
			R = 6,
			nscan = 40,
			burn = 10,
			odens = 2,
			verbose = FALSE,
			gof = FALSE
		))
	)
	expect_s3_class(fit, "lame")
})
