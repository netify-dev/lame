.reporting_private288 = function(name) get(name, environment(lame::lame), inherits = TRUE)

# capture the sampler coordinates just before reporting
.reporting_capture288 = function(original) {
	records = new.env(parent = emptyenv()); records$draws = list()
	hook = function(frame) {
		before = .Random.seed
		field = function(name, default = NULL) if (exists(name, frame, inherits = FALSE)) get(name, frame, inherits = FALSE) else default
		X = field("Xlist"); if (is.null(X)) X = list(field("X"))
		a = field("a"); b = field("b"); beta = field("beta"); U = field("U"); V = field("V")
		bip = field("bip", identical(field("mode"), "bipartite")); G = if (bip) field("G") else NULL
		eta = array(0, c(length(a), length(b), length(X)))
		for (t in seq_along(X)) {
			z = outer(a, b, "+");p = dim(X[[t]])[3L]
			for (k in seq_len(p)) z = z + X[[t]][, , k] * beta[k]
			u = if (length(dim(U)) == 3L) matrix(U[, , 1L], dim(U)[1L], dim(U)[2L]) else U
			v = if (length(dim(V)) == 3L) matrix(V[, , 1L], dim(V)[1L], dim(V)[2L]) else V
			if (!is.null(u) && NCOL(u) > 0L) z = z + if (bip) u %*% G %*% t(v) else tcrossprod(u, v)
			if (isTRUE(field("symmetric"))) z = (z + t(z)) / 2
			eta[, , t] = z
		}
		records$draws[[length(records$draws) + 1L]] = list(
			eta = eta,
			beta = beta,
			a = a,
			b = b,
			U = U,
			V = V,
			G = G,
			rng = before
		)
		stopifnot(identical(.Random.seed, before))
		invisible(NULL)
	}
	insertions = 0L
	transform = function(x) {
		if (!is.call(x)) return(x)
		if (identical(x, quote(beta_rep <- beta))) {
			insertions <<- insertions + 1L
			return(quote({ .capture_reporting288(environment()); beta_rep = beta }))
		}
		for (j in seq_along(x)[-1L]) {
			if (identical(x[[j]], quote(expr = ))) next
			x[j] = list(transform(x[[j]]))
		}
		x
	}
	f = original;body(f) = transform(body(original));stopifnot(insertions == 1L)
	e = new.env(parent = environment(original));e$.capture_reporting288 = hook;environment(f) = e
	list(fitter = f, records = records)
}

.reporting_cases288 = function() {
	grid = expand.grid(
		fitter = c("ame", "lame"),
		mode = c("directed", "symmetric", "bipartite"),
		intercept = c(FALSE, TRUE),
		margin = c("both", "row", "column", "none"),
		stringsAsFactors = FALSE
	)
	grid = grid[grid$mode != "symmetric" | grid$margin %in% c("both", "none"), ]
	grid$id = paste(grid$fitter, grid$mode, grid$margin, grid$intercept, sep = "_")
	grid$default_margin_flags = FALSE
	defaults = grid[grid$mode == "symmetric" & grid$margin == "none" & !grid$intercept, ]
	defaults$default_margin_flags = TRUE
	defaults$id = paste0(defaults$id, "_default_margin_flags")
	rbind(grid, defaults)
}

.reporting_run288 = function(case, save_log_lik = TRUE) {
	set.seed(601288);n = 6L;m = if (case$mode == "bipartite") 5L else n
	X = array(rnorm(n * m), c(n, m, 1L));Y = matrix(rnorm(n * m), n, m)
	if (case$mode == "symmetric") { X[, , 1] = (X[, , 1] + t(X[, , 1])) / 2;Y = (Y + t(Y)) / 2 }
	if (case$mode != "bipartite") diag(Y) = NA
	args = list(
		Y = Y,
		Xdyad = X,
		family = "normal",
		R = 1L,
		intercept = case$intercept,
		rvar = case$margin %in% c("both", "row"),
		cvar = case$margin %in% c("both", "column"),
		nvar = case$margin == "both",
		symmetric = case$mode == "symmetric",
		dcor = FALSE,
		mode = if (case$mode == "bipartite") "bipartite" else "unipartite",
		nscan = 4L,
		burn = 2L,
		odens = 1L,
		seed = 601288L,
		verbose = FALSE,
		gof = FALSE,
		posterior_opts = list(save_ab = TRUE, save_UV = TRUE),
		save_log_lik = save_log_lik
	)
	if (isTRUE(case$default_margin_flags)) { args$rvar = NULL;args$cvar = NULL }
	if (case$fitter == "lame") { args$Y = list(Y, Y);args$Xdyad = list(X, X);args$plot = FALSE }
	fit = suppressWarnings(do.call(getExportedValue("lame", case$fitter), args))
	list(fit = fit, Y = if (case$fitter == "lame") array(rep(Y, 2), c(n, m, 2L)) else array(Y, c(n, m, 1L)))
}

test_that("stored joint predictors retain the actual restricted-model sampler state", {
	skip_on_cran()
	ame_hook = .reporting_capture288(.reporting_private288("ame_unipartite"))
	lame_hook = .reporting_capture288(lame::lame)
	local_mocked_bindings(ame_unipartite = ame_hook$fitter, lame = lame_hook$fitter, .package = "lame")
	cases = subset(.reporting_cases288(), !(fitter == "ame" & mode == "bipartite"));results = list()
	for (i in seq_len(nrow(cases))) {
		case = cases[i, ];collector = if (case$fitter == "ame") ame_hook$records else lame_hook$records
		collector$draws = list();out = .reporting_run288(case);fit = out$fit
		# bipartite ame has no centering block; its u/g/v draw is checked separately below as the reference route
		if (case$fitter == "ame" && case$mode == "bipartite") next
		captured = collector$draws;expect_length(captured, nrow(fit$BETA))
		state = .reporting_private288(".joint_predictive_state")(fit)
		maximum = 0
		for (s in seq_along(captured)) for (t in seq_len(dim(captured[[s]]$eta)[3L])) {
			eta = .reporting_private288(".joint_predictor")(fit, state, s, t)
			maximum = max(maximum, max(abs(unname(eta) - captured[[s]]$eta[, , t])))
			expect_lt(max(abs(unname(eta) - captured[[s]]$eta[, , t])), 1e-10)
			index = fit$log_lik_index;cols = which(index[, 3L] == t);pairs = index[cols, 1:2, drop = FALSE]
			ref = dnorm(out$Y[, , t][pairs], captured[[s]]$eta[, , t][pairs], sqrt(fit$VC[s, "ve"]), log = TRUE)
			expect_equal(fit$log_lik[s, cols], ref, tolerance = 1e-11)
		}
		if (case$mode == "symmetric" && case$margin == "none") {
			expect_identical(fit$rvar, FALSE);expect_identical(fit$cvar, FALSE)
		}
		if (case$mode == "symmetric" && case$margin == "both" && !case$intercept)
			expect_equal(fit$a_samples, if (is.null(fit$b_samples)) fit$a_samples else fit$b_samples, tolerance = 1e-12)
		results[[case$id]] = list(maximum_eta_error = maximum, saved_draws = length(captured))
	}
	if (exists(
		".reporting_results288",
		.GlobalEnv,
		inherits = FALSE
	)) assign(".reporting_results288", results, .GlobalEnv)
})

test_that("bipartite AME preserves raw per-draw interaction coordinates", {
	cases = subset(.reporting_cases288(), fitter == "ame" & mode == "bipartite")
	for (i in seq_len(nrow(cases))) {
		out = .reporting_run288(cases[i, ]);fit = out$fit;state = .reporting_private288(".joint_predictive_state")(fit)
		for (s in seq_len(nrow(fit$BETA))) {
			eta = matrix(0, nrow(out$Y), dim(out$Y)[2L]);beta = fit$BETA[s, ]
			for (k in seq_along(beta)) eta = eta + beta[k] * fit$X[, , k]
			if (isTRUE(fit$rvar)) eta = eta + fit$a_samples[, s]
			if (isTRUE(fit$cvar)) eta = sweep(eta, 2, fit$b_samples[, s], "+")
			u = matrix(fit$U_samples[, , s], nrow(eta), 1);v = matrix(fit$V_samples[, , s], ncol(eta), 1)
			eta = eta + u %*% matrix(fit$G_samples[, , s], 1, 1) %*% t(v)
			expect_equal(
				unname(.reporting_private288(".joint_predictor")(fit, state, s)),
				unname(eta),
				tolerance = 1e-12
			)
		}
	}
})
