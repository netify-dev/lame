test_that("binary bipartite AR1 fits expose the exact scale move and opt-out", {
	rn = paste0("row", 1:5); cn = paste0("col", 1:6)
	Y = lapply(1:3, function(t) matrix(as.numeric(sin(seq_len(30) / 4 + t) > .1), 5, 6, dimnames = list(rn, cn)))
	names(Y) = paste0("year", 1:3)
	Y[[1]][2, 4] = NA; Y[[3]][4, 1] = NA
	args = list(
		Y = Y,
		mode = "bipartite",
		family = "binary",
		R = 1L,
		dynamic_uv = TRUE,
		dynamic_uv_kind = "ar1",
		rvar = FALSE,
		cvar = FALSE,
		burn = 20L,
		nscan = 40L,
		odens = 5L,
		seed = 47335623L,
		verbose = FALSE,
		plot = FALSE,
		gof = FALSE,
		posterior_opts = list(save_UV = TRUE, thin_UV = 1L)
	)
	set.seed(47335601L); caller = .Random.seed
	expect_no_warning({ fit = do.call(lame, args) })
	expect_identical(caller, .Random.seed)
	expect_equal(unname(fit$bipartite_scale_move["attempts"]), 50)
	expect_gt(unname(fit$bipartite_scale_move["accepts"]), 0)
	expect_equal(unname(fit$bipartite_scale_move["numerical_rejections"]), 0)
	expect_true(all(fit$sampler_health$numerical_failures == 0))
	expect_equal(nrow(fit$BETA), 8L)
	args$prior = list(bipartite_scale_move = FALSE)
	expect_no_warning({ off = do.call(lame, args) })
	expect_identical(caller, .Random.seed)
	expect_null(off$bipartite_scale_move)
	off_prior = off$prior; off_prior[["bipartite_scale_move"]] = NULL
	expect_identical(off_prior, fit$prior)
	expect_identical(off$g, fit$g)
})
