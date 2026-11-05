checkpoint_panel = function(bip = FALSE, family = "normal", symmetric = FALSE) {
	set.seed(406)
	nr = 8L; nc = if (bip) 7L else nr; nt = 4L
	Y = array(
		rnorm(nr * nc * nt),
		c(nr, nc, nt),
		dimnames = list(paste0("a", seq_len(nr)), paste0(if (bip) "b" else "a", seq_len(nc)), paste0("t", seq_len(nt)))
	)
	X = array(
		rnorm(nr * nc * nt),
		c(nr, nc, 1L, nt),
		dimnames = list(dimnames(Y)[[1]], dimnames(Y)[[2]], "x", dimnames(Y)[[3]])
	)
	if (family == "ordinal") Y = array(as.numeric(cut(Y, c(-Inf, -.4, .5, Inf))), dim(Y), dimnames(Y))
	for (tt in seq_len(nt)) {
		if (symmetric) { Y[, , tt] = (Y[, , tt] + t(Y[, , tt])) / 2; X[, , 1, tt] = (X[, , 1, tt] + t(X[, , 1, tt])) / 2 }
		if (!bip) diag(Y[, , tt]) = NA_real_
	}
	list(
		Y = Y,
		Xdyad = X,
		family = family,
		mode = if (bip) "bipartite" else "unipartite",
		symmetric = symmetric,
		R = 1L,
		burn = 20L,
		nscan = 30L,
		odens = 5L,
		seed = 406L,
		gof = FALSE,
		plot = FALSE,
		verbose = FALSE
	)
}

checkpoint_without_sweeps = function(args, state) {
	args$start_vals = state; args$burn = 0L; args$nscan = 10L
	args$checkpoint_path = tempfile(fileext = ".rds"); args$max_seconds = 1e-12
	suppressWarnings(suppressMessages(do.call(lame, args)))
	out = readRDS(args$checkpoint_path)
	unlink(args$checkpoint_path)
	out
}

test_that("resuming preserves dynamic coefficient paths and actor-level decomposition", {
	skip_on_cran()
	for (bip in c(FALSE, TRUE)) {
		a = checkpoint_panel(bip)
		a$dynamic_beta = c("intercept", "x_dyad")
		a$dynamic_beta_kind = "rw1"; a$dynamic_beta_pool = "both"
		a$dynamic_ab = TRUE
		f = suppressWarnings(suppressMessages(do.call(lame, a)))
		state = f$start_vals$checkpoint_state$values
		expect_true(any(abs(diff(state$beta_dyn_path[, 1])) > 1e-8))
		expect_true(any(abs(state$a_lev - rowMeans(state$a_mat)) > 1e-8))
		ck = checkpoint_without_sweeps(a, f$start_vals)
		got = ck$start_vals$checkpoint_state$values
		expect_equal(got$beta_dyn_path, state$beta_dyn_path, tolerance = 0)
		expect_equal(got$beta_static, state$beta_static, tolerance = 0)
		expect_equal(got$sigma_beta_by_group, state$sigma_beta_by_group, tolerance = 0)
		expect_equal(got$rho_beta_by_group, state$rho_beta_by_group, tolerance = 0)
		expect_equal(got$a_lev, state$a_lev, tolerance = 0)
		expect_equal(got$b_lev, state$b_lev, tolerance = 0)
		expect_equal(got$a_mat, state$a_mat, tolerance = 0)
		expect_equal(got$b_mat, state$b_mat, tolerance = 0)
		expect_equal(got$Z, state$Z, tolerance = 0)
	}
})

test_that("resuming preserves transition auxiliaries and the symmetric gauge", {
	skip_on_cran()
	for (kind in c("ar1", "snap", "t", "symmetric")) {
		a = checkpoint_panel(symmetric = kind == "symmetric")
		a$dynamic_ab = TRUE; a$dynamic_uv = TRUE
		a$dynamic_uv_kind = if (kind == "symmetric") "ar1" else kind
		f = suppressWarnings(suppressMessages(do.call(lame, a)))
		ck = checkpoint_without_sweeps(a, f$start_vals)
		fields = c(
			"U_cube", "V_cube", "rho_uv", "sigma_uv", "L_dyn", "lambda_u", "lambda_v", "pi_snap", "delta_uv", "delta_uv_v"
		)
		for (nm in fields) expect_equal(
			ck$start_vals$checkpoint_state$values[[nm]],
			f$start_vals$checkpoint_state$values[[nm]],
			tolerance = 0,
			info = paste(kind, nm)
		)
	}
})

test_that("resuming preserves explicit ordinal cutpoints", {
	skip_on_cran()
	a = checkpoint_panel(family = "ordinal"); a$ordinal_cutpoints = "explicit"
	f = suppressWarnings(suppressMessages(do.call(lame, a)))
	ck = checkpoint_without_sweeps(a, f$start_vals)
	expect_equal(ck$start_vals$checkpoint_state$values$alpha, f$start_vals$checkpoint_state$values$alpha, tolerance = 0)
	expect_equal(
		ck$start_vals$checkpoint_state$values$tau_prop,
		f$start_vals$checkpoint_state$values$tau_prop,
		tolerance = 0
	)
})

test_that("resume completes unfinished burn-in and rejects incompatible layouts", {
	skip_on_cran()
	a = checkpoint_panel(); a$R = 0L; a$burn = 40L; a$max_seconds = 1e-12
	ck = tempfile(fileext = ".rds"); a$checkpoint_path = ck
	suppressWarnings(suppressMessages(do.call(lame, a)))
	saved = readRDS(ck)
	expect_identical(saved$iter, 0L)
	expect_equal(saved$burn_remaining, 40L)
	f = suppressWarnings(suppressMessages(lame_resume(ck, nscan_more = 50L)))
	expect_identical(f$call$burn, 40L)
	expect_equal(nrow(f$BETA), 10L)
	expect_error(
		suppressWarnings(lame_resume(ck, nscan_more = 50L, dynamic_ab = TRUE)),
		class = "lame_checkpoint_layout"
	)
	saved$start_vals$checkpoint_state = NULL; saved$burn_remaining = NULL
	saveRDS(saved, ck)
	expect_warning(f2 <- suppressMessages(lame_resume(ck, nscan_more = 50L)), class = "lame_checkpoint_legacy")
	expect_equal(f2$call$burn, 40L)
	unlink(ck)
})


test_that("pooled persistence retains its hierarchical state across resume", {
	skip_on_cran()
	a = checkpoint_panel()
	a$dynamic_beta = c("intercept", "x_dyad")
	a$dynamic_beta_pool = "rho"
	a$R = 0L
	f = suppressMessages(do.call(lame, a))
	state = f$start_vals$checkpoint_state$values
	expect_true(is.finite(state$pool_rho_mu))
	expect_true(is.finite(state$pool_rho_tau) && state$pool_rho_tau > 0)
	expect_false(identical(state$pool_rho_tau, 1))
	ck = checkpoint_without_sweeps(a, f$start_vals)
	expect_identical(ck$start_vals$checkpoint_state$values$pool_rho_mu, state$pool_rho_mu)
	expect_identical(ck$start_vals$checkpoint_state$values$pool_rho_tau, state$pool_rho_tau)
})
