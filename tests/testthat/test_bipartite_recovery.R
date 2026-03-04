# bipartite parameter recovery

# align true params with the alphabetically sorted actor order
align_true_params <- function(true_vals, orig_names, sorted_names) {
	names(true_vals) <- orig_names
	true_vals[sorted_names]
}

test_that("ame() bipartite normal recovers beta and additive effects", {
	skip_on_cran()

	nA <- 15
	nB <- 10
	R <- 2
	true_intercept <- 1.0
	true_beta <- 0.8

	dat <- simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 1,
		family = "normal", R = R,
		beta_intercept = true_intercept,
		beta_dyad = true_beta,
		sd_a = 0.5, sd_b = 0.5, s2 = 1,
		mode = "bipartite", seed = 123
	)

	fit <- ame(
		dat$Y, Xdyad = dat$Xdyad,
		R = R, family = "normal", mode = "bipartite",
		burn = 2000, nscan = 8000, odens = 2,
		print = FALSE, gof = FALSE, seed = 42
	)

	# beta recovery
	beta_hat <- colMeans(fit$BETA)
	expect_true(abs(beta_hat[1] - true_intercept) < 0.5,
		info = paste("Intercept:", beta_hat[1], "vs true:", true_intercept))
	expect_true(abs(beta_hat[2] - true_beta) < 0.5,
		info = paste("Beta dyad:", beta_hat[2], "vs true:", true_beta))

	# additive effects recovery (align by name)
	true_a <- align_true_params(dat$true_params$a[[1]],
		rownames(dat$Y), names(fit$APM))
	true_b <- align_true_params(dat$true_params$b[[1]],
		colnames(dat$Y), names(fit$BPM))
	cor_a <- cor(fit$APM, true_a)
	cor_b <- cor(fit$BPM, true_b)
	expect_true(cor_a > 0.4,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.4,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))

	# latent position recovery (UV product correlation)
	# true UV is in original order, est UV is in sorted order
	# correlation of vectorized products is order-independent if we
	# align both to the same actor ordering
	true_UV <- tcrossprod(dat$true_params$U[[1]], dat$true_params$V[[1]])
	est_UV <- fit$U %*% fit$G %*% t(fit$V)
	# reorder true_UV to match sorted actor order
	row_order <- match(names(fit$APM), rownames(dat$Y))
	col_order <- match(names(fit$BPM), colnames(dat$Y))
	true_UV_aligned <- true_UV[row_order, col_order]
	cor_UV <- cor(c(true_UV_aligned), c(est_UV))
	expect_true(cor_UV > 0.2,
		info = paste("UV product correlation:", round(cor_UV, 3)))

	# variance s2 should be near 1
	s2_hat <- mean(fit$VC[, "ve"])
	expect_true(abs(s2_hat - 1) < 1,
		info = paste("s2:", round(s2_hat, 3), "vs true: 1"))
})


test_that("ame() bipartite binary recovers beta and additive effects", {
	skip_on_cran()

	nA <- 20
	nB <- 12
	R <- 2

	dat <- simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 1,
		family = "binary", R = R,
		beta_intercept = 0, beta_dyad = 0.8,
		sd_a = 0.4, sd_b = 0.4, s2 = 1,
		mode = "bipartite", seed = 456
	)

	fit <- ame(
		dat$Y, Xdyad = dat$Xdyad,
		R = R, family = "binary", mode = "bipartite",
		burn = 2000, nscan = 8000, odens = 2,
		print = FALSE, gof = FALSE, seed = 42
	)

	# beta_dyad should be positive
	beta_hat <- colMeans(fit$BETA)
	expect_true(beta_hat[2] > 0,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "should be positive"))

	# additive effects correlation (align by name)
	true_a <- align_true_params(dat$true_params$a[[1]],
		rownames(dat$Y), names(fit$APM))
	true_b <- align_true_params(dat$true_params$b[[1]],
		colnames(dat$Y), names(fit$BPM))
	cor_a <- cor(fit$APM, true_a)
	cor_b <- cor(fit$BPM, true_b)
	expect_true(cor_a > 0.3,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.3,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))

	# UV product correlation
	true_UV <- tcrossprod(dat$true_params$U[[1]], dat$true_params$V[[1]])
	est_UV <- fit$U %*% fit$G %*% t(fit$V)
	row_order <- match(names(fit$APM), rownames(dat$Y))
	col_order <- match(names(fit$BPM), colnames(dat$Y))
	true_UV_aligned <- true_UV[row_order, col_order]
	cor_UV <- cor(c(true_UV_aligned), c(est_UV))
	expect_true(cor_UV > 0.15,
		info = paste("UV product correlation:", round(cor_UV, 3)))
})


test_that("lame() bipartite normal recovers beta and additive effects", {
	skip_on_cran()

	nA <- 12
	nB <- 8
	n_time <- 3
	true_intercept <- 0.5
	true_beta <- 0.6

	dat <- simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = n_time,
		family = "normal", R = 0,
		beta_intercept = true_intercept,
		beta_dyad = true_beta,
		sd_a = 0.5, sd_b = 0.5, s2 = 1,
		mode = "bipartite", seed = 789
	)

	# lame() expects Xdyad as list of matrices
	Xdyad_list <- lapply(dat$Xdyad, function(x) x[,,1])

	fit <- lame(
		dat$Y, Xdyad = Xdyad_list,
		R = 0, family = "normal", mode = "bipartite",
		burn = 1000, nscan = 5000, odens = 2,
		print = FALSE, gof = FALSE, seed = 42
	)

	# beta recovery
	beta_hat <- colMeans(fit$BETA)
	expect_true(abs(beta_hat[1] - true_intercept) < 0.6,
		info = paste("Intercept:", round(beta_hat[1], 3), "vs true:", true_intercept))
	expect_true(abs(beta_hat[2] - true_beta) < 0.6,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "vs true:", true_beta))

	# additive effects (align by name: lame sorts actors alphabetically)
	orig_row_names <- rownames(dat$Y[[1]])
	orig_col_names <- colnames(dat$Y[[1]])
	sorted_row_names <- names(fit$APM)
	sorted_col_names <- names(fit$BPM)

	true_a_avg <- rowMeans(do.call(cbind, dat$true_params$a))
	true_b_avg <- rowMeans(do.call(cbind, dat$true_params$b))
	true_a_aligned <- align_true_params(true_a_avg, orig_row_names, sorted_row_names)
	true_b_aligned <- align_true_params(true_b_avg, orig_col_names, sorted_col_names)

	cor_a <- cor(fit$APM, true_a_aligned)
	cor_b <- cor(fit$BPM, true_b_aligned)
	expect_true(cor_a > 0.3,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.3,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))

	# s2 should be near 1
	s2_hat <- mean(fit$VC[, "ve"])
	expect_true(abs(s2_hat - 1) < 1.5,
		info = paste("s2:", round(s2_hat, 3), "vs true: 1"))
})


test_that("lame() bipartite binary recovers beta and effects", {
	skip_on_cran()

	nA <- 15
	nB <- 10
	n_time <- 3

	dat <- simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = n_time,
		family = "binary", R = 0,
		beta_intercept = 0, beta_dyad = 0.8,
		sd_a = 0.4, sd_b = 0.4, s2 = 1,
		mode = "bipartite", seed = 101
	)

	Xdyad_list <- lapply(dat$Xdyad, function(x) x[,,1])

	fit <- lame(
		dat$Y, Xdyad = Xdyad_list,
		R = 0, family = "binary", mode = "bipartite",
		burn = 1000, nscan = 5000, odens = 2,
		print = FALSE, gof = FALSE, seed = 42
	)

	# beta_dyad should be positive
	beta_hat <- colMeans(fit$BETA)
	expect_true(beta_hat[2] > 0,
		info = paste("Beta dyad:", round(beta_hat[2], 3), "should be positive"))

	# additive effects correlation (align by name)
	orig_row_names <- rownames(dat$Y[[1]])
	orig_col_names <- colnames(dat$Y[[1]])

	true_a_avg <- rowMeans(do.call(cbind, dat$true_params$a))
	true_b_avg <- rowMeans(do.call(cbind, dat$true_params$b))
	true_a_aligned <- align_true_params(true_a_avg, orig_row_names, names(fit$APM))
	true_b_aligned <- align_true_params(true_b_avg, orig_col_names, names(fit$BPM))

	cor_a <- cor(fit$APM, true_a_aligned)
	cor_b <- cor(fit$BPM, true_b_aligned)
	expect_true(cor_a > 0.2,
		info = paste("Sender effects correlation:", round(cor_a, 3)))
	expect_true(cor_b > 0.2,
		info = paste("Receiver effects correlation:", round(cor_b, 3)))
})


test_that("lame() rejects unsupported bipartite families", {
	nA <- 8
	nB <- 6
	dat <- simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 2,
		family = "normal", R = 0,
		mode = "bipartite", seed = 999
	)
	Xdyad_list <- lapply(dat$Xdyad, function(x) x[,,1])

	for (fam in c("ordinal", "cbin", "frn", "poisson")) {
		expect_error(
			lame(dat$Y, Xdyad = Xdyad_list, R = 0, family = fam,
				mode = "bipartite", burn = 10, nscan = 10, odens = 1,
				print = FALSE, gof = FALSE, seed = 42),
			"not supported for bipartite",
			info = paste("Family:", fam, "should be rejected for bipartite")
		)
	}
})


test_that("ame() bipartite dimensions are correct", {
	nA <- 10
	nB <- 7
	R <- 2

	dat <- simulate_test_network(
		n = nA, nA = nA, nB = nB, n_time = 1,
		family = "normal", R = R,
		mode = "bipartite", seed = 222
	)

	fit <- ame(
		dat$Y, Xdyad = dat$Xdyad,
		R = R, family = "normal", mode = "bipartite",
		burn = 100, nscan = 200, odens = 2,
		print = FALSE, gof = FALSE, seed = 42
	)

	# Check output dimensions
	expect_equal(length(fit$APM), nA)
	expect_equal(length(fit$BPM), nB)
	expect_equal(nrow(fit$U), nA)
	expect_equal(ncol(fit$U), R)
	expect_equal(nrow(fit$V), nB)
	expect_equal(ncol(fit$V), R)
	expect_equal(nrow(fit$YPM), nA)
	expect_equal(ncol(fit$YPM), nB)
	expect_equal(nrow(fit$G), R)
	expect_equal(ncol(fit$G), R)
})
