# dynamic_ab makes the intercept a dynamic_beta coefficient, a random walk next to ar1 blocks
skip_on_cran()

.abi_panel = function(n = 10L, nt = 4L, seed = 5, family = "normal", mu0 = 0) {
	set.seed(seed)
	nm = sprintf("a%02d", seq_len(n))
	X = lapply(seq_len(nt), function(t) array(rnorm(n * n), c(n, n, 1L), dimnames = list(nm, nm, "x")))
	shift = seq(-0.6, 0.6, length.out = nt)
	Y = lapply(seq_len(nt), function(t) {
		eta = mu0 + shift[t] + 0.5 * X[[t]][, , 1] + outer(rnorm(n, 0, 0.3), rnorm(n, 0, 0.3), "+")
		m = if (family == "binary") 1 * (eta + matrix(rnorm(n * n), n) > 0) else eta + matrix(rnorm(n * n), n)
		diag(m) = NA; dimnames(m) = list(nm, nm); m
	})
	list(Y = Y, X = X)
}
.abi_fit = function(d, ...) suppressWarnings(suppressMessages(lame(
	d$Y,
	Xdyad = d$X,
	R = 0,
	dynamic_ab = TRUE,
	burn = 30,
	nscan = 200,
	odens = 5,
	chains = 1,
	verbose = FALSE,
	gof = FALSE,
	...
)))

test_that("the rule adds the intercept to every form of dynamic_beta", {
	f = function(db, ...) .lame_ab_dynamic_intercept(db, TRUE, TRUE, "normal", ...)
	expect_identical(f(FALSE), "intercept")
	expect_identical(f(NULL), "intercept")
	expect_identical(f(TRUE), TRUE)
	expect_identical(f("dyad"), c("intercept", "dyad"))
	expect_identical(f(c("intercept", "dyad")), c("intercept", "dyad"))
	expect_equal(f(2), c(1, 2))
	expect_identical(f(c(FALSE, TRUE)), c(TRUE, TRUE))
	# untouched: no dynamic_ab, no intercept, ordinal
	expect_identical(.lame_ab_dynamic_intercept(FALSE, FALSE, TRUE, "normal"), FALSE)
	expect_identical(.lame_ab_dynamic_intercept(FALSE, TRUE, FALSE, "normal"), FALSE)
	expect_identical(.lame_ab_dynamic_intercept("dyad", TRUE, TRUE, "ordinal"), "dyad")
})

test_that("the added intercept is rw1 and the blocks the user names keep ar1", {
	d = .abi_panel()
	fa = .abi_fit(d)
	expect_true(isTRUE(fa$dynamic_intercept_from_ab))
	expect_identical(fa$dynamic_beta_kind, "rw1")
	expect_identical(fa$dynamic_beta_kind_requested, "ar1")
	expect_identical(fa$dynamic_beta_kind_by_block, c(intercept = "rw1"))
	expect_true(all(fa$RHO_BETA[, "intercept"] == 1))
	fm = .abi_fit(d, dynamic_beta = "dyad")
	expect_identical(fm$dynamic_beta_kind, "ar1")
	expect_identical(
		fm$dynamic_beta_kind_by_block[c("dyad", "intercept")],
		c(dyad = "ar1", intercept = "rw1")
	)
	expect_true(all(fm$RHO_BETA[, "intercept"] == 1))
	expect_true(all(fm$RHO_BETA[, "dyad"] < 1) && stats::sd(fm$RHO_BETA[, "dyad"]) > 0)
	# an intercept the user names keeps the user's kind
	fu = .abi_fit(d, dynamic_beta = "intercept")
	expect_false(isTRUE(fu$dynamic_intercept_from_ab))
	expect_identical(fu$dynamic_beta_kind_by_block, c(intercept = "ar1"))
	# the rw1 block raises neither the weak-persistence note nor the near-unit-root note
	out = paste(capture.output(print(summary(fa))), collapse = "\n")
	expect_false(grepl("weakly informed by the data for block", out))
	expect_false(grepl("concentrated near 1", out))
	expect_match(out, "kind = rw1")
	# an intercept = FALSE fit keeps static coefficients
	expect_length(dim(.abi_fit(d, intercept = FALSE)$BETA), 2L)
})

test_that("the forecast holds the last intercept level of a binary panel", {
	d = .abi_panel(n = 12L, family = "binary", mu0 = -1.5, seed = 9)
	fit = .abi_fit(d, family = "binary")
	# freeze every transition but the intercept's, so an rw1 intercept carries the last period forward
	fz = fit
	fz$rho_ab[] = 1; fz$sigma_ab[] = 0; fz$SIGMA_BETA[] = 0
	fc = suppressWarnings(predict(fz, h = 2L, type = "link"))
	off = function(m) m[row(m) != col(m)]
	expect_equal(off(fc[[1L]]), off(fit$EZ[[4L]]), tolerance = 1e-10)
	expect_equal(off(fc[[2L]]), off(fit$EZ[[4L]]), tolerance = 1e-10)
	# an ar1 intercept the user names reverts toward its long-run mean (fit$MU_BETA), not toward zero
	fu = .abi_fit(d, family = "binary", dynamic_beta = "intercept")
	fu$rho_ab[] = 1; fu$sigma_ab[] = 0; fu$SIGMA_BETA[] = 0
	fcu = suppressWarnings(predict(fu, h = 2L, type = "link"))
	step = mean(off(fcu[[2L]])) - mean(off(fu$EZ[[4L]]))
	gap = mean(fu$MU_BETA[, "intercept"]) - mean(fu$BETA[, "intercept", 4L])
	expect_gt(abs(step), 0.01)
	expect_identical(sign(step), sign(gap))
	# no unit-root warning comes from the rw1 block
	expect_no_warning(predict(fit, h = 1L, type = "link"))
})

test_that("a checkpoint of a fit with mixed kinds resumes", {
	d = .abi_panel()
	ck = tempfile(fileext = ".rds")
	on.exit(unlink(ck), add = TRUE)
	f1 = .abi_fit(d, dynamic_beta = "dyad", checkpoint_path = ck, checkpoint_every = 50L)
	expect_true(file.exists(ck))
	f2 = suppressWarnings(suppressMessages(lame_resume(ck, nscan = 100L)))
	expect_identical(f2$dynamic_beta_kind_by_block, f1$dynamic_beta_kind_by_block)
	expect_true(all(f2$RHO_BETA[, "intercept"] == 1))
	expect_equal(dim(f2$BETA)[3], 4L)
})

test_that("lame(method = \"als\") resolves the same intercept and kinds", {
	d = .abi_panel(n = 12L)
	a0 = suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		R = 0,
		dynamic_ab = TRUE,
		method = "als",
		verbose = FALSE
	)))
	expect_identical(a0$dynamic$beta_kind, "rw1")
	a1 = suppressWarnings(suppressMessages(lame(
		d$Y,
		Xdyad = d$X,
		R = 0,
		dynamic_ab = TRUE,
		dynamic_beta = "dyad",
		method = "als",
		verbose = FALSE
	)))
	m1 = .abi_fit(d, dynamic_beta = "dyad")
	expect_setequal(names(a1$dynamic_beta_kind_by_block), names(m1$dynamic_beta_kind_by_block))
	expect_identical(
		a1$dynamic_beta_kind_by_block[c("dyad", "intercept")],
		m1$dynamic_beta_kind_by_block[c("dyad", "intercept")]
	)
	expect_equal(unname(a1$rho_beta[["intercept"]]), 1)
	expect_lt(a1$rho_beta[["dyad"]], 1)
	# penalised fit with the rw1 intercept path, called as lame() calls it
	a2 = suppressWarnings(suppressMessages(lame:::lame_dynamic_als(
		d$Y,
		Xdyad = d$X,
		R = 0,
		dynamic_ab = TRUE,
		dynamic_beta = c("intercept", "dyad"),
		beta_kind_intercept = "rw1",
		bias_correct = "none",
		verbose = FALSE
	)))
	expect_s3_class(a2, "lame_dynamic_als")
	expect_true(all(is.finite(a2$coef_path)))
})

test_that("the gate monitors a static coefficient of a dynamic fit once", {
	d = .abi_panel()
	fit = .abi_fit(d)
	nm = colnames(.lame_gate_series(fit))
	expect_true("x_dyad" %in% nm)
	expect_false(any(grepl("^x_dyad\\[", nm)))
	expect_true(all(paste0("intercept[t", 1:4, "]") %in% nm))
	# the persistence of an rw1 block is not monitored, that of an ar1 block is
	expect_false(any(grepl("^rho_beta", nm)))
	expect_true("sigma_beta[intercept]" %in% nm)
	nm2 = colnames(.lame_gate_series(.abi_fit(d, dynamic_beta = "dyad")))
	expect_true("rho_beta[dyad]" %in% nm2)
	expect_false("rho_beta[intercept]" %in% nm2)
})
