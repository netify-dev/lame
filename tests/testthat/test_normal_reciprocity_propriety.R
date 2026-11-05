skip_on_cran()

test_that("normal reciprocity refusal uses an exact conservative observed-data certificate", {
	check = getFromNamespace(".ame_check_normal_reciprocity", "lame")
	n = 8L
	Y = outer(seq_len(n), seq_len(n), "+") / 7
	diag(Y) = NA_real_
	X = array(1, c(n, n, 1L))
	args = list(
		Y = Y,
		X = X,
		rvar = TRUE,
		cvar = TRUE,
		dcor = TRUE,
		symmetric = FALSE,
		R = 0,
		family = "normal"
	)
	expect_error(do.call(check, args), class = "lame_improper_normal_reciprocity")
	# one unequal response invalidates the all-pairs boundary proof
	unequal = args
	unequal$Y[1, 2] = unequal$Y[1, 2] + 1e-14
	expect_null(do.call(check, unequal))
	# an exact finite inequality counts even when a symmetry warning uses all.equal
	expect_true(isTRUE(all.equal(unequal$Y, t(unequal$Y))))
	for (field in c("dcor", "symmetric", "family")) {
		outside = args
		outside[[field]] = switch(field, dcor = FALSE, symmetric = TRUE, family = "bin")
		expect_null(do.call(check, outside))
	}
	# one-sided observations are not reciprocal constraints
	one_side = args
	one_side$Y[lower.tri(Y)] = NA_real_
	expect_null(do.call(check, one_side))
	# a forest has full-row-rank skew incidence and is not certified improper
	forest = args
	forest$Y[, ] = NA_real_
	for (i in seq_len(n - 1L)) forest$Y[i, i + 1L] = forest$Y[i + 1L, i] = i
	certificate = do.call(check, forest)
	expect_equal(certificate$reciprocal_pairs, n - 1L)
	expect_equal(certificate$antisymmetric_rank_bound, n - 1L)
	no_actor = forest
	no_actor$rvar = no_actor$cvar = FALSE
	expect_error(do.call(check, no_actor), class = "lame_improper_normal_reciprocity")
	# enough asymmetric columns make the conservative bound inconclusive
	full_columns = args
	full_columns$X = array(0, c(n, n, 28L))
	pairs = which(upper.tri(Y), arr.ind = TRUE)
	for (k in seq_len(nrow(pairs))) full_columns$X[pairs[k, 1L], pairs[k, 2L], k] = 1
	expect_equal(do.call(check, full_columns)$antisymmetric_rank_bound, 35L)
	# disconnected components are counted exactly, isolated actors included; the posterior is improper once the pairs
	# exceed the rank bound by four, the excess the 2 Beta(2, 2) - 1 prior absorbs
	expect_identical(lame:::.lame_rho_improper_dyads(), 4L)
	graph = function(edges) {
		g = args
		g$Y[, ] = NA_real_
		for (ij in edges) g$Y[ij[1L], ij[2L]] = g$Y[ij[2L], ij[1L]] = 1
		g
	}
	clique = function(nodes) combn(nodes, 2L, simplify = FALSE)
	# a triangle: 3 pairs, rank 2, one pair too many
	certificate = do.call(check, graph(clique(1:3)))
	expect_equal(certificate$reciprocal_pairs, 3L)
	expect_equal(certificate$actor_rank_bound, 2L)
	# four actors all paired: 6 pairs, rank 3, three too many
	expect_equal(do.call(check, graph(clique(1:4)))$reciprocal_pairs, 6L)
	# four too many: the same four actors plus a fifth tied to two of them (8 pairs, rank 4)
	expect_error(
		do.call(check, graph(c(clique(1:4), list(c(1L, 5L), c(2L, 5L))))),
		class = "lame_improper_normal_reciprocity"
	)
	# separate components add their ranks: a four-clique and a triangle (9 pairs, rank 5), the eighth actor isolated
	expect_error(do.call(check, graph(c(clique(1:4), clique(5:7)))), class = "lame_improper_normal_reciprocity")
	# a four-clique and a path of two pairs (8 pairs, rank 6) is three too many
	expect_equal(do.call(check, graph(c(clique(1:4), list(c(5L, 6L), c(6L, 7L)))))$reciprocal_pairs, 8L)
})

test_that("proper forest initialization uses one unchanged Gaussian conditional draw", {
	check = getFromNamespace(".ame_check_normal_reciprocity", "lame")
	refresh = getFromNamespace(".ame_normal_reciprocity_initial_draw", "lame")
	Y = matrix(1.3, 2L, 2L)
	diag(Y) = NA_real_
	X = precomputeX(design_array(n = 2L, intercept = TRUE))
	certificate = check(Y, X, TRUE, TRUE, TRUE, FALSE, 0, "normal")
	Z = matrix(c(.1, 1.3, 1.3, -.2), 2L, 2L)
	EZ = matrix(0, 2L, 2L)
	Sab = matrix(c(.8, .23, .23, 1.1), 2L, 2L)
	set.seed(47333511)
	actual = refresh(Z, EZ, Sab, .2, X, 1.4, EZ, 4, certificate, TRUE, TRUE)
	actual_seed = .Random.seed
	set.seed(47333511)
	expected = rbeta_ab_fc(Z, Sab, .2, X, 1.4, offset = EZ, g = 4)
	expect_identical(actual, expected)
	expect_identical(.Random.seed, actual_seed)
	before = .Random.seed
	non_degenerate = Z
	non_degenerate[1L, 2L] = non_degenerate[1L, 2L] + .1
	expect_null(refresh(non_degenerate, EZ, Sab, .2, X, 1.4, EZ, 4, certificate, TRUE, TRUE))
	expect_identical(.Random.seed, before)
	expect_error(
		refresh(Z, EZ, Sab, .2, X, 1.4, EZ, 0, certificate, TRUE, TRUE),
		class = "lame_normal_reciprocity_initialization"
	)
	expect_identical(.Random.seed, before)
	zero_contrast = matrix(1, 2L, 2L)
	expect_error(
		refresh(Z, EZ, zero_contrast, .2, X, 1.4, EZ, 4, certificate, TRUE, TRUE),
		class = "lame_normal_reciprocity_initialization"
	)
	expect_identical(.Random.seed, before)
	# the helper tests every complete latent reciprocal pair, not only observed graph edges
	Y3 = matrix(NA_real_, 3L, 3L)
	Y3[1L, 2L] = Y3[2L, 1L] = .4
	X3 = precomputeX(design_array(n = 3L, intercept = TRUE))
	certificate3 = check(Y3, X3, TRUE, TRUE, TRUE, FALSE, 0, "normal")
	Z3 = matrix(.4, 3L, 3L)
	Z3[1L, 3L] = .7
	expect_null(refresh(
		Z3,
		matrix(0, 3L, 3L),
		Sab,
		.2,
		X3,
		1.4,
		matrix(0, 3L, 3L),
		4,
		certificate3,
		TRUE,
		TRUE
	))
	expect_identical(.Random.seed, before)
})

test_that("a design with cycles the prior lets stand gets the same exact refresh when the next conditional needs it", {
	check = getFromNamespace(".ame_check_normal_reciprocity", "lame")
	refresh = getFromNamespace(".ame_normal_reciprocity_initial_draw", "lame")
	Sab = matrix(c(.8, .23, .23, 1.1), 2L, 2L)
	for (n in 3:4) {
		Y = matrix(.4, n, n)
		diag(Y) = NA_real_
		X = precomputeX(design_array(n = n, intercept = TRUE))
		certificate = check(Y, X, TRUE, TRUE, TRUE, FALSE, 0, "normal")
		expect_equal(certificate$reciprocal_pairs, n * (n - 1L) / 2L)
		# not a forest: the actor rank is n - 1
		expect_lt(certificate$actor_rank_bound, certificate$reciprocal_pairs)
		Z = matrix(.4, n, n)
		EZ = matrix(0, n, n)
		set.seed(47334001)
		actual = refresh(Z, EZ, Sab, .2, X, 1.4, EZ, 4, certificate, TRUE, TRUE)
		if (n == 3L) {
			# three dyads: the conditional of rho is proper without it
			expect_null(actual)
		} else {
			set.seed(47334001)
			expect_identical(actual, rbeta_ab_fc(Z, Sab, .2, X, 1.4, offset = EZ, g = 4))
			# without actor effects there is nothing to draw and the rho update returns the state itself
			expect_null(refresh(Z, EZ, Sab, .2, X, 1.4, EZ, 4, certificate, FALSE, FALSE))
			expect_null(refresh(Z, EZ, Sab, .2, X, 1.4, EZ, 0, certificate, TRUE, TRUE))
		}
	}
})

test_that("ame() fits a complete symmetric four-actor panel that the arc-sine prior had to refuse", {
	Y = matrix(c(NA, .4, .8, 1.1, .4, NA, 1.2, .3, .8, 1.2, NA, .6, 1.1, .3, .6, NA), 4L, 4L)
	dimnames(Y) = list(letters[1:4], letters[1:4])
	set.seed(47334002)
	before = .Random.seed
	expect_message({fit = ame(
		Y,
		R = 0,
		g = 4,
		nscan = 60,
		burn = 20,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)}, "exactly symmetric")
	expect_s3_class(fit, "ame")
	expect_true(all(is.finite(fit$VC)))
	expect_true(all(abs(fit$VC[, "rho"]) < 1))
	# the same panel with five actors is improper under the prior and refused
	Y5 = matrix(0, 5L, 5L)
	Y5[upper.tri(Y5)] = seq(.2, 1.1, length.out = 10L)
	Y5 = Y5 + t(Y5)
	diag(Y5) = NA_real_
	dimnames(Y5) = list(letters[1:5], letters[1:5])
	expect_message(
		expect_error(
			ame(
				Y5,
				R = 0,
				g = 4,
				nscan = 6,
				burn = 2,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE
			),
			class = "lame_improper_normal_reciprocity"
		),
		"exactly symmetric"
	)
})

test_that("proper observed forest fits with an explicit proper regression prior", {
	Y = matrix(NA_real_, 3L, 3L)
	Y[1L, 2L] = Y[2L, 1L] = .4
	Y[2L, 3L] = Y[3L, 2L] = .8
	dimnames(Y) = list(c("a", "b", "c"), c("a", "b", "c"))
	set.seed(47333512)
	before = .Random.seed
	expect_message({fit = ame(
		Y,
		R = 0,
		g = 4,
		nscan = 6,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)}, "exactly symmetric")
	expect_s3_class(fit, "ame")
	expect_true(all(is.finite(fit$VC)))
	expect_identical(.Random.seed, before)
})

test_that("normal coefficient priors refuse invalid g before RNG and exempt unused g", {
	check = getFromNamespace(".ame_check_normal_g", "lame")
	set.seed(47333515)
	before = .Random.seed
	for (g in list(0, -1, Inf, NaN, NA_real_, numeric(0), c(1, 2), "4", TRUE)) {
		expect_error(check(g, "normal", 1L), class = "lame_invalid_normal_g")
		expect_identical(.Random.seed, before)
		expect_null(check(g, "normal", 0L))
	}
	expect_null(check(4, "normal", 1L))
	expect_identical(.Random.seed, before)
	Y = matrix(c(NA, .4, .8, .4, NA, 1.2, .8, 1.2, NA), 3L, 3L)
	dimnames(Y) = list(c("a", "b", "c"), c("a", "b", "c"))
	expect_message(
		expect_error(
			ame(
				Y,
				R = 0,
				g = 0,
				dcor = FALSE,
				nscan = 6,
				burn = 2,
				odens = 1,
				verbose = FALSE,
				gof = FALSE,
				plot = FALSE
			),
			class = "lame_invalid_normal_g"
		),
		"exactly symmetric"
	)
	expect_identical(.Random.seed, before)
})

test_that("higher-rank normal refusal uses only its sufficient complete-pair bound", {
	check = getFromNamespace(".ame_check_normal_reciprocity", "lame")
	n = 8L
	Y = outer(seq_len(n), seq_len(n), "+") / 7
	diag(Y) = NA_real_
	X = array(1, c(n, n, 1L))
	set.seed(47334711)
	before = .Random.seed
	for (rank in c(1L, 2L)) {
		expect_error(
			check(Y, X, TRUE, TRUE, TRUE, FALSE, rank, "normal"),
			class = "lame_improper_normal_reciprocity"
		)
		expect_identical(.Random.seed, before)
	}
	# failing this sufficient inequality is unclassified, never certified proper
	expect_null(check(Y, X, TRUE, TRUE, TRUE, FALSE, 3L, "normal"))
	unequal = Y
	unequal[1L, 2L] = unequal[1L, 2L] + 1e-14
	expect_null(check(unequal, X, TRUE, TRUE, TRUE, FALSE, 1L, "normal"))
	expect_true(isTRUE(all.equal(unequal, t(unequal))))
	expect_null(check(Y, X, TRUE, TRUE, FALSE, FALSE, 1L, "normal"))
	expect_null(check(Y, X, TRUE, TRUE, TRUE, TRUE, 1L, "normal"))
	expect_null(check(Y, X, TRUE, TRUE, TRUE, FALSE, 1L, "binary"))
	# one-way observations add no reciprocal determinant factors
	one_side = Y
	one_side[lower.tri(Y)] = NA_real_
	expect_null(check(one_side, X, TRUE, TRUE, TRUE, FALSE, 1L, "normal"))
	# the conservative asymmetric-column count can make the result inconclusive
	full_columns = array(0, c(n, n, 28L))
	pairs = which(upper.tri(Y), arr.ind = TRUE)
	for (k in seq_len(nrow(pairs))) full_columns[pairs[k, 1L], pairs[k, 2L], k] = 1
	expect_null(check(Y, full_columns, TRUE, TRUE, TRUE, FALSE, 1L, "normal"))
	# refusal needs the pairs to exceed the actor-plus-factor bound (7 + 8) by the four dyads the prior absorbs
	threshold = matrix(NA_real_, n, n)
	tree = cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)
	remaining = pairs[!(pairs[, 2L] == pairs[, 1L] + 1L), , drop = FALSE]
	connected = rbind(tree, remaining)
	for (k in seq_len(15L)) threshold[connected[k, 1L], connected[k, 2L]] = threshold[connected[k, 2L], connected[k, 1L]] = 1
	expect_null(check(threshold, X, TRUE, TRUE, TRUE, FALSE, 1L, "normal"))
	for (k in 16:18) threshold[connected[k, 1L], connected[k, 2L]] = threshold[connected[k, 2L], connected[k, 1L]] = 1
	expect_null(check(threshold, X, TRUE, TRUE, TRUE, FALSE, 1L, "normal"))
	threshold[connected[19L, 1L], connected[19L, 2L]] = threshold[connected[19L, 2L], connected[19L, 1L]] = 1
	expect_error(
		check(threshold, X, TRUE, TRUE, TRUE, FALSE, 1L, "normal"),
		class = "lame_improper_normal_reciprocity"
	)
	# higher-rank checks never pass an R0 forest certificate to initialization
	forest = matrix(NA_real_, n, n)
	for (i in seq_len(n - 1L)) forest[i, i + 1L] = forest[i + 1L, i] = i
	expect_null(check(forest, X, TRUE, TRUE, TRUE, FALSE, 1L, "normal"))
	expect_identical(.Random.seed, before)
})

test_that("public higher-rank normal refusal preserves RNG and fixed-rho scope", {
	Y = outer(seq_len(8L), seq_len(8L), "+") / 7
	diag(Y) = NA_real_
	set.seed(47334712)
	before = .Random.seed
	for (rank in c(1L, 2L)) {
		expect_message(
			expect_error(
				ame(
					Y,
					R = rank,
					family = "normal",
					g = 4,
					nscan = 6,
					burn = 2,
					odens = 1,
					verbose = FALSE,
					gof = FALSE,
					plot = FALSE
				),
				class = "lame_improper_normal_reciprocity"
			),
			"exactly symmetric"
		)
		expect_identical(.Random.seed, before)
	}
	expect_message({fit = ame(
		Y,
		R = 1L,
		family = "normal",
		g = 4,
		dcor = FALSE,
		nscan = 6,
		burn = 2,
		odens = 1,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE
	)}, "exactly symmetric")
	expect_s3_class(fit, "ame")
	expect_true(all(is.finite(fit$VC)))
	expect_identical(.Random.seed, before)
})
