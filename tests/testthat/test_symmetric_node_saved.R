skip_on_cran()

# symmetric fits saved by lame 1.4.2 (fixtures/make_symmetric_142.R): the sampler spelled the node
# covariate `growth` as `gnodeth_node` and the fast route kept Xcol apart as W_col
old = readRDS(test_path("fixtures", "symmetric_142.rds"))
off = function(p) { diag(p) = NA; p }

test_that("a saved sampler fit with the 1.4.2 node name predicts and simulates", {
	m = old$mcmc
	expect_identical(colnames(m$BETA), c("intercept", "gnodeth_node", "d_dyad"))
	expect_equal(predict(m), old$mcmc_pred)
	expect_equal(off(predict(m, newdata = old$Xd)), predict(m))
	expect_false(isTRUE(all.equal(off(predict(m, newdata = old$Xd * 0)), predict(m))))
	expect_s3_class(simulate(m, nsim = 1), "ame.sim")
})

test_that("a saved dynamic_beta sampler fit with the 1.4.2 node name finds change points and forecasts", {
	f = old$lame_db
	expect_identical(dimnames(f$BETA)[[2L]], c("intercept", "gnodeth_node"))
	expect_no_error(detect_change_point(f))
	# the stored link took the pooled coefficient positionally; the rebuilt one matches newdata
	nd = lapply(f$Xlist, function(X) X[, , integer(0), drop = FALSE])
	gap = Map(function(a, b) max(abs(a - b), na.rm = TRUE), predict(f, type = "link"),
		predict(f, newdata = nd, type = "link"))
	expect_lt(max(unlist(gap)), 1e-10)
	p = predict(f, h = 1)
	expect_true(all(is.finite(unlist(p)[!is.na(unlist(p))])))
})

test_that("a saved symmetric fast fit with W_col refits to its coefficients and predicts", {
	f = old$als
	expect_false(is.null(f$W_col))
	expect_equal(coef(suppressWarnings(ame_als_refit(f))), old$als_refit, tolerance = 1e-6)
	expect_equal(predict(f), old$als_pred)
})
