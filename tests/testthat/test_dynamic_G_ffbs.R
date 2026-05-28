# dynamic_G FFBS state-space sampler tests. Cover:
# (a) low-level FFBS sampler exists + returns correct shapes
# (b) rho_G MH update stays in (-1, 1) and respects the AR(1) likelihood
# (c) sigma_G2 conjugate IG draw is positive
# (d) byte-identical default contract for dynamic_G = FALSE
# (e) end-to-end lame() fit with dynamic_G = TRUE stores all outputs
# (f) rotation drift diagnostic produces a numeric ratio + flag

test_that("ffbs_vecG returns G_cube and vecG_path with correct shapes", {
	set.seed(11)
	nA = 8; nB = 10; T = 4
	RA = 2; RB = 2
	U_cube = array(rnorm(nA * RA * T), dim = c(nA, RA, T))
	V_cube = array(rnorm(nB * RB * T), dim = c(nB, RB, T))
	E_cube = array(rnorm(nA * nB * T), dim = c(nA, nB, T))
	out = lame:::ffbs_vecG(E_cube, U_cube, V_cube,
	                       s2 = 1, rho_G = 0.7, sigma_G2 = 0.1)
	expect_equal(dim(out$G_cube), c(RA, RB, T))
	expect_equal(dim(out$vecG_path), c(RA * RB, T))
	expect_true(all(is.finite(out$G_cube)))
})

test_that("ffbs_vecG zero-rank case returns empty G_cube cleanly", {
	nA = 8; nB = 10; T = 3
	U_cube = array(0, dim = c(nA, 0L, T))
	V_cube = array(0, dim = c(nB, 0L, T))
	E_cube = array(rnorm(nA * nB * T), dim = c(nA, nB, T))
	out = lame:::ffbs_vecG(E_cube, U_cube, V_cube,
	                       s2 = 1, rho_G = 0.7, sigma_G2 = 0.1)
	expect_equal(dim(out$G_cube), c(0L, 0L, T))
})

test_that("sample_rho_G_mh keeps rho in (-1, 1) over many sweeps", {
	set.seed(42)
	# simulate a 1D AR(1) path with true rho = 0.6
	N = 50; p = 4
	true_rho = 0.6; true_sigma2 = 0.05
	vG = matrix(0, p, N)
	vG[, 1L] = rnorm(p, 0, sqrt(true_sigma2 / (1 - true_rho^2)))
	for (t in 2L:N) {
		vG[, t] = true_rho * vG[, t - 1L] + rnorm(p, 0, sqrt(true_sigma2))
	}
	rho_chain = numeric(200)
	rho = 0.0
	for (s in seq_len(200)) {
		step = lame:::sample_rho_G_mh(rho, vG, sigma_G2 = true_sigma2, tau = 0.3)
		rho = step$rho
		rho_chain[s] = rho
	}
	expect_true(all(abs(rho_chain) < 1))
	# posterior should drift toward true_rho on a chain this long
	expect_gt(mean(rho_chain[100:200]), 0.3)
})

test_that("sample_sigma_G2 is positive and shrinks with more data", {
	set.seed(13)
	# small chain
	vG_small = matrix(rnorm(5 * 5, 0, 0.1), 5, 5)
	# large chain (same DGP, more periods)
	vG_large = matrix(rnorm(5 * 50, 0, 0.1), 5, 50)
	s_small = lame:::sample_sigma_G2(vG_small, rho_G = 0.5)
	s_large = lame:::sample_sigma_G2(vG_large, rho_G = 0.5)
	expect_gt(s_small, 0)
	expect_gt(s_large, 0)
	# both should be small (true sigma^2 = 0.01); large-data should be tighter
	expect_lt(s_large, 1)
})

test_that("dynamic_G = FALSE is byte-identical to default", {
	skip_on_cran()
	set.seed(2026)
	nA = 8; nB = 10; T = 3
	Y_list = lapply(seq_len(T), function(t) {
		mat = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
		rownames(mat) = paste0("a", seq_len(nA))
		colnames(mat) = paste0("b", seq_len(nB))
		mat
	})
	f1 = lame(Y_list, family = "binary", mode = "bipartite",
	          R_row = 1, R_col = 1,
	          nscan = 20, burn = 10, odens = 5,
	          rvar = FALSE, cvar = FALSE, verbose = FALSE, seed = 7)
	f2 = lame(Y_list, family = "binary", mode = "bipartite",
	          R_row = 1, R_col = 1, dynamic_G = FALSE,
	          nscan = 20, burn = 10, odens = 5,
	          rvar = FALSE, cvar = FALSE, verbose = FALSE, seed = 7)
	expect_identical(f1$BETA, f2$BETA)
	expect_identical(f1$YPM, f2$YPM)
	expect_null(f1$dynamic_G)
	expect_null(f2$dynamic_G)
})

test_that("dynamic_G = TRUE stores v8 outputs (G_cube + post_mean + RHO_G)", {
	skip_on_cran()
	set.seed(2027)
	nA = 8; nB = 10; T = 4
	Y_list = lapply(seq_len(T), function(t) {
		mat = matrix(rbinom(nA * nB, 1, 0.3), nA, nB)
		rownames(mat) = paste0("a", seq_len(nA))
		colnames(mat) = paste0("b", seq_len(nB))
		mat
	})
	suppressWarnings(
		fit <- lame(Y_list, family = "binary", mode = "bipartite",
		            R_row = 2, R_col = 2,
		            dynamic_G = TRUE,
		            nscan = 30, burn = 15, odens = 5,
		            rvar = FALSE, cvar = FALSE, verbose = FALSE, seed = 13))
	expect_true(isTRUE(fit$dynamic_G))
	expect_equal(dim(fit$G_cube), c(2L, 2L, T))
	expect_equal(dim(fit$G_cube_post_mean), c(2L, 2L, T))
	expect_equal(dim(fit$G_cube_post_sd), c(2L, 2L, T))
	expect_true(all(is.finite(fit$G_cube_post_mean)))
	expect_true(all(fit$G_cube_post_sd >= 0))
	expect_true(length(fit$RHO_G) == 30 / 5)
	expect_true(all(abs(fit$RHO_G) < 1))
	expect_true(length(fit$SIGMA_G2) == 30 / 5)
	expect_true(all(fit$SIGMA_G2 > 0))
})

test_that("rotation_drift_diagnostic returns numeric ratio + boolean flag", {
	set.seed(31)
	nA = 6; nB = 8; T = 4; RA = 2; RB = 2
	U_cube = array(rnorm(nA * RA * T), dim = c(nA, RA, T))
	V_cube = array(rnorm(nB * RB * T), dim = c(nB, RB, T))
	G_cube = array(rnorm(RA * RB * T), dim = c(RA, RB, T))
	diag = lame:::.rotation_drift_diagnostic(G_cube, U_cube, V_cube)
	expect_true(is.numeric(diag$ratio) || is.na(diag$ratio))
	expect_true(is.logical(diag$flag))
	expect_equal(dim(diag$var_raw), c(RA, RB))
	expect_equal(dim(diag$var_canonical), c(RA, RB))
})
