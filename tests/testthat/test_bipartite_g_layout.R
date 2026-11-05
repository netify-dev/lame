bipartite_g_layout_fixture = function() {
	periods = c("P3", "P1", "P2")
	Y = X = vector("list", length(periods))
	for (t in seq_along(periods)) {
		y = outer(seq_len(5L), seq_len(6L), function(i, j) ((i + 2L * j + t) %% 4L < 2L) * 1)
		dimnames(y) = list(paste0("r", seq_len(5L)), paste0("c", seq_len(6L)))
		x = array(0, c(5L, 6L, 3L), dimnames = list(rownames(y), colnames(y), c("x", "z", "w")))
		for (k in seq_len(3L)) x[, , k] = outer(
			seq_len(5L),
			seq_len(6L),
			function(i, j) ((i + 3L * j + 2L * t + k) %% 7L - 3L) / 4
		)
		Y[[t]] = y; X[[t]] = x
	}
	names(Y) = names(X) = periods
	list(
		Y = Y,
		Xdyad = X,
		family = "binary",
		mode = "bipartite",
		R = 0L,
		dynamic_ab = TRUE,
		burn = 2L,
		nscan = 8L,
		odens = 2L,
		seed = 53081L,
		gof = FALSE,
		verbose = FALSE,
		posterior_opts = list(save_ab = TRUE)
	)
}

test_that("bipartite scalar g follows the fitted coefficient layout", {
	withr::local_preserve_seed()
	args = bipartite_g_layout_fixture()
	fit = do.call(lame, args)
	expect_identical(fit$g, rep(25, 4L))
	expect_identical(length(fit$g), dim(fit$Xlist[[1L]])[3L])
	expect_identical(colnames(fit$BETA), c("intercept", "x_dyad", "z_dyad", "w_dyad"))
	args$g = 7
	fit = do.call(lame, args)
	expect_identical(fit$g, rep(7, 4L))
	args$dynamic_ab = FALSE
	fit = do.call(lame, args)
	expect_identical(fit$g, rep(7, 4L))
	args$intercept = FALSE
	fit = do.call(lame, args)
	expect_identical(fit$g, rep(7, 3L))
	expect_identical(colnames(fit$BETA), c("x_dyad", "z_dyad", "w_dyad"))
})

test_that("mixed bipartite coefficients expand scalar g for the final static column", {
	withr::local_preserve_seed()
	args = bipartite_g_layout_fixture()
	args$dynamic_beta = "intercept"
	args$g = 7
	fit = do.call(lame, args)
	expect_identical(fit$g, rep(7, 4L))
	expect_equal(fit$tryErrorChecks$beta, 0L)
	args$g = rep(7, 4L)
	explicit = do.call(lame, args)
	expect_identical(fit$refit_options$g, 7)
	expect_identical(explicit$refit_options$g, rep(7, 4L))
	# the call is kept and differs only in the requested g
	expect_identical(fit$call[[1L]], explicit$call[[1L]])
	expect_identical(eval(fit$call[["g"]], envir = baseenv()), 7)
	expect_identical(eval(explicit$call[["g"]], envir = baseenv()), rep(7, 4L))
	call_scalar = as.list(fit$call); call_vector = as.list(explicit$call)
	expect_identical(call_scalar[names(call_scalar) != "g"], call_vector[names(call_vector) != "g"])
	fit$call[["g"]] = explicit$call[["g"]]
	fit$refit_options$g = explicit$refit_options$g
	expect_identical(fit, explicit)
})

test_that("bipartite supplied g vectors retain their exact values", {
	withr::local_preserve_seed()
	args = bipartite_g_layout_fixture()
	args$g = c(2, 3, 5, 7)
	fit = do.call(lame, args)
	expect_identical(fit$g, args$g)
	expect_identical(length(fit$g), ncol(fit$BETA))
})
