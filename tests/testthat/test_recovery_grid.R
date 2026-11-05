skip_on_cran()

# fast slice of the recovery grid (helper_recovery_grid.R): one cell per estimator and family group at a short chain
# (the full grid is described in inst/recovery/README.md)

quick = function(cell) {
	cell$nscan = 150L
	cell$burn = 50L
	cell$odens = 5L
	cell$n = 16L
	cell$nA = 14L
	cell$nB = 10L
	cell
}

test_that("the recovery grid's fast slice runs clean across estimators", {
	cells = list(
		quick(grid_cell("q_ame_normal_R2_dcor", "normal", "ame", R = 2, dcor = TRUE, group = "quick")),
		quick(grid_cell("q_ame_binary_sym_R1", "binary", "ame", symmetric = TRUE, R = 1, group = "quick")),
		quick(
			grid_cell(
				"q_ame_bip_poisson_R11", "poisson", "ame", mode = "bipartite", R_row = 1, R_col = 1, group = "quick"
			)
		),
		quick(grid_cell("q_ame_bip_ordinal_R00", "ordinal", "ame", mode = "bipartite", group = "quick")),
		quick(grid_cell("q_lame_binary_dynab", "binary", "lame", T = 3, dynamic_ab = TRUE, group = "quick")),
		quick(
			grid_cell(
				"q_lame_poisson_sym_dynuv", "poisson", "lame", T = 3, symmetric = TRUE, R = 1, dynamic_uv = TRUE, group = "quick"
			)
		),
		quick(grid_cell("q_lame_normal_dynbeta", "normal", "lame", T = 3, dynamic_beta = TRUE, group = "quick")),
		quick(grid_cell("q_lame_frn_miss", "frn", "lame", T = 3, missing = 0.2, group = "quick")),
		quick(
			grid_cell(
				"q_lame_bip_binary_dynuv", "binary", "lame", T = 3, mode = "bipartite", R_row = 1, R_col = 1, dynamic_uv = TRUE, group = "quick"
			)
		),
		grid_cell("q_ame_als_normal_R2", "normal", "ame_als", R = 2, group = "quick"),
		grid_cell("q_lame_als_binary", "binary", "lame_als", T = 3, group = "quick"),
		grid_cell(
			"q_lame_dynamic_als_poisson", "poisson", "lame_dynamic_als", T = 3, dynamic_ab = TRUE, group = "quick"
		),
		grid_cell(
			"q_ame_als_ordinal_expected", "ordinal", "ame_als", expect_error = "ordinal|not supported|ALS", group = "quick"
		)
	)
	csv = tempfile(fileext = ".csv")
	df = grid_run(cells, seeds = 1L, out_csv = csv, verbose = FALSE)
	expect_equal(nrow(df), length(cells))
	expect_true(all(df$status %in% c("ok", "expected_error")))
	ok = df[df$status == "ok", ]
	# exact post-processing identities on every MCMC fit
	mc = ok[!grepl("als", ok$estimator), ]
	expect_true(all(mc$inv_pred < 1e-6, na.rm = TRUE))
	expect_true(all(mc$inv_fit == 0, na.rm = TRUE))
	expect_true(all(mc$inv_res == 0, na.rm = TRUE))
	expect_true(all(mc$inv_coef == 0, na.rm = TRUE))
	expect_true(all(mc$inv_ci == 1, na.rm = TRUE))
	expect_true(all(mc$inv_summary == 1, na.rm = TRUE))
	# sampler fallbacks stay rare, no errors recorded by the post-processing calls
	expect_true(all(mc$tryerr_frac <= 0.05, na.rm = TRUE))
	expect_false(any(grepl("predict_newdata:|fitted:|residuals:|simulate:|coef:|summary:|forecast", ok$notes)))
	# the fitted predictor tracks the truth even at this chain length
	expect_true(all(ok$ez_cor > 0.5, na.rm = TRUE))
})
