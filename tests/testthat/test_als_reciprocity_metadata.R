# small saved-shape fixtures exercise reporting and public dispatch without fitting
als_rho_fixture = function(legacy = FALSE, family = "normal") {
	n = 4L; TT = 3L; nm = paste0("actor", seq_len(n)); times = c("1999", "2000", "2001")
	a = matrix(c(-.3, .1, .2, 0), n, TT, dimnames = list(nm, times))
	b = matrix(c(.2, -.2, .1, -.1), n, TT, dimnames = list(nm, times))
	eta = lapply(seq_len(TT), function(t) {
		x = .2 + outer(a[, t], b[, t], "+");diag(x) = NA_real_;dimnames(x) = list(nm, nm);x
	});names(eta) = times
	y = array(NA_real_, c(n, n, TT), dimnames = list(nm, nm, times))
	for (t in seq_len(TT)) {
		x = if (family == "poisson") matrix(seq_len(n * n) %% 5L, n) else eta[[t]] + matrix(seq_len(n * n) / 40, n)
		diag(x) = NA_real_;y[, , t] = x
	}
	path = setNames(c(.7, -.2, .6), times)
	fit = list(
		family = family,
		mode = "unipartite",
		symmetric = FALSE,
		R = 0L,
		n_time = TT,
		dims = list(n_row = n, n_col = n, p = 1L),
		row_names = nm,
		col_names = nm,
		time_names = times,
		a = rowMeans(a),
		b = rowMeans(b),
		APM = rowMeans(a),
		BPM = rowMeans(b),
		a_dynamic = a,
		b_dynamic = b,
		coefficients = c(intercept = .2),
		coef_path = matrix(.2, 1L, TT, dimnames = list("intercept", times)),
		BETA = matrix(.2, 1L, 1L, dimnames = list("point", "intercept")),
		U = NULL,
		V = NULL,
		L = NULL,
		Omat = matrix(0, n, n),
		Oarr = array(0, c(n, n, TT)),
		EZ = eta,
		fitted = eta,
		Y = y,
		VC = c(va = .5, cab = 0, vb = .6, rho = .3, ve = .4),
		s2 = .4,
		dynamic_ab = TRUE,
		dynamic_uv = FALSE,
		dynamic_beta = FALSE,
		dynamic_G = FALSE,
		dynamic_rho = legacy,
		rho_path = path,
		rho_ab = .8,
		sigma_ab = .2,
		converged = TRUE,
		iterations = 3L,
		deviance = 1,
		objective_label = "penalized least squares",
		diagnostics = list(stationarity = list(pass = TRUE)),
		call = quote(lame(method = "als"))
	)
	if (!legacy) {fit$residual_rho_path = path;fit$rho_path_role = "residual_diagnostic"}
	class(fit) = c("lame_dynamic_als", "lame_als", "ame_als")
	fit
}

with_absent_rng_als = function(code) {
	had = exists(".Random.seed", .GlobalEnv, inherits = FALSE)
	old = if (had) get(".Random.seed", .GlobalEnv) else NULL
	on.exit(
		{if (had) assign(".Random.seed", old, .GlobalEnv) else if (exists(
			".Random.seed",
			.GlobalEnv,
			inherits = FALSE
		)) rm(".Random.seed", envir = .GlobalEnv)})
	if (had) rm(".Random.seed", envir = .GlobalEnv)
	force(code)
}

test_that("ALS warm starts use modeled pooled rho for current and serialized legacy paths", {
	for (legacy in c(FALSE, TRUE)) for (family in c("normal", "poisson")) {
		f = als_rho_fixture(legacy, family);p = tempfile(fileext = ".rds");saveRDS(f, p);restored = readRDS(p);unlink(p)
		before = restored
		with_absent_rng_als({
			s = als_start_vals(restored, jitter = 0)
			expect_identical(s$rho, unname(restored$VC[["rho"]]))
			expect_length(s$rho, 1L)
			expect_identical(s$a, restored$a_dynamic)
			expect_identical(s$b, restored$b_dynamic)
			expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
		})
		set.seed(811903);rng = .Random.seed
		s = als_start_vals(restored, jitter = .01, seed = 812903)
		expect_identical(.Random.seed, rng)
		expect_identical(s$rho, .3)
		expect_identical(restored, before)
	}
	# a scalar-only ALS object has the same covariance contract
	f = als_rho_fixture(TRUE);f$rho_path = NULL
	expect_identical(als_start_vals(f)$rho, .3)
	# this helper does not convert posterior MCMC paths
	m = f;class(m) = "lame"
	expect_error(als_start_vals(m), "must be an ALS fit")
})

test_that("legacy diagnostic flags do not change public ALS simulation or fitted-period prediction", {
	f = als_rho_fixture(FALSE);old = als_rho_fixture(TRUE);f_before = f;old_before = old
	set.seed(811904);rng = .Random.seed
	method = getS3method("simulate", "ame_als")
	expect_identical(method, get("simulate.ame_als", asNamespace("lame")))
	public = stats::simulate(f, nsim = 2L, seed = 812904)
	expect_identical(.Random.seed, rng)
	expect_identical(public, stats::simulate(old, nsim = 2L, seed = 812904))
	expect_identical(public, method(f, nsim = 2L, seed = 812904))
	expect_identical(.Random.seed, rng)
	expect_identical(stats::predict(f, type = "link"), f$EZ)
	expect_identical(stats::predict(old, type = "link"), f$EZ)
	expect_identical(summary(f)$VC, summary(old)$VC)
	expect_error(stats::predict(f, h = 1L), "does not support h-step forecasting")
	expect_error(stats::predict(old, h = 1L), "does not support h-step forecasting")
	expect_identical(f, f_before);expect_identical(old, old_before);expect_identical(.Random.seed, rng)
})

test_that("unsupported public ALS reciprocity paths fail before fitting and RNG", {
	Y = lapply(
		seq_len(3L),
		function(t) {x = matrix(seq_len(16L) / 10, 4L);diag(x) = NA_real_;dimnames(x) = list(letters[1:4], letters[1:4]);x}
	)
	names(Y) = as.character(1999:2001)
	for (ab in c(FALSE, TRUE)) {
		with_absent_rng_als({
			expect_error(
				lame(
					Y,
					family = "normal",
					R = 0L,
					method = "als",
					dynamic_rho = TRUE,
					dynamic_ab = ab,
					verbose = FALSE
				),
				class = "lame_als_dynamic_rho"
			)
			expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
		})
		set.seed(811905);rng = .Random.seed
		expect_error(
			lame(Y, family = "normal", R = 0L, method = "als", dynamic_rho = TRUE, dynamic_ab = ab, verbose = FALSE),
			"Use.*mcmc",
			class = "lame_als_dynamic_rho"
		)
		expect_identical(.Random.seed, rng)
	}
})

test_that("Poisson latent moment diagnostics and count correction preserve the path alias", {
	f = als_rho_fixture(TRUE, "poisson");f$link = "log"
	set.seed(811906);rng = .Random.seed
	updated = lame:::.ae_update_response_rho(f, "log")
	expect_false(updated$dynamic_rho)
	expect_identical(updated$rho_path_role, "residual_diagnostic")
	expect_identical(updated$residual_rho_path, updated$rho_path)
	expect_identical(updated$rho_path_working, f$rho_path)
	expect_identical(updated$EZ, f$EZ)
	expect_identical(.Random.seed, rng)
	moments = list(
		variance = .2,
		rho = .4,
		rho_path = setNames(c(.1, .2, .3), f$time_names),
		boundary = FALSE,
		path_boundary = rep(FALSE, 3L)
	)
	updated$mu = unname(updated$coefficients["intercept"])
	corrected = lame:::.lda_apply_count_moments(updated, moments)
	expect_false(corrected$dynamic_rho)
	expect_identical(corrected$residual_rho_path, moments$rho_path)
	expect_identical(corrected$rho_path, moments$rho_path)
	expect_identical(corrected$rho_path_role, "residual_diagnostic")
	expect_identical(corrected$raw$rho_path, updated$rho_path)
	expect_identical(corrected$raw$residual_rho_path, updated$residual_rho_path)
	expect_identical(unname(corrected$VC["rho"]), .4)
	expect_identical(unname(corrected$VC["ve"]), .2)
	expect_identical(als_start_vals(corrected)$rho, .4)
	expect_identical(.Random.seed, rng)
})
