# uv_plot dims, nodes, large sides, heatmap and latent_similarity on small
# hand-built fits (no sampling)

uvd_fit = function(n_r = 12, n_c = 10, d = c(5, 3, 1.5), bipartite = TRUE, seed = 1) {
	set.seed(seed)
	k = length(d)
	A = qr.Q(qr(matrix(rnorm(n_r * k), n_r, k)))
	B = qr.Q(qr(matrix(rnorm(n_c * k), n_c, k)))
	P = A %*% diag(d, k) %*% t(B)
	rn = sprintf("r%02d", seq_len(n_r))
	cn = sprintf("c%02d", seq_len(n_c))
	if (!bipartite) cn = rn
	dimnames(P) = list(rn, cn)
	fit = list(
		U = matrix(0, n_r, k, dimnames = list(rn, NULL)),
		V = matrix(0, n_c, k, dimnames = list(cn, NULL)),
		UVPM = P,
		mode = if (bipartite) "bipartite" else "unipartite",
		symmetric = FALSE
	)
	if (bipartite) fit$G = diag(k)
	class(fit) = "ame"
	fit
}

uvd_dyn = function(n = 8, k = 3, n_t = 3) {
	set.seed(2)
	nm = sprintf("a%02d", 1:n)
	fit = list(
		U = array(rnorm(n * k * n_t), c(n, k, n_t), dimnames = list(nm, NULL, paste0("t", 1:n_t))),
		V = array(rnorm(n * k * n_t), c(n, k, n_t), dimnames = list(nm, NULL, paste0("t", 1:n_t))),
		dynamic_uv = TRUE,
		uv_positions = "posterior_mean_product",
		mode = "unipartite",
		symmetric = FALSE
	)
	class(fit) = "lame"
	fit
}

test_that("dims picks the dimensions and labels them with their shares", {
	fit = uvd_fit()
	st = lame:::.lame_dim_strength(fit)
	p = uv_plot(fit, layout = "biplot")
	expect_s3_class(p, "ggplot")
	expect_equal(p$labels$x, sprintf("Dimension 1 (%s%% of latent product)", formatC(100 * st$share[1], format = "f", digits = 0)))
	expect_match(p$labels$caption, sprintf("Dimensions 1 and 2 of 3 shown: %s%% of the latent product", formatC(100 * sum(st$share[1:2]), format = "f", digits = 0)), fixed = TRUE)
	p13 = uv_plot(fit, layout = "biplot", dims = c(1, 3))
	expect_match(p13$labels$y, "^Dimension 3 \\(")
	# three dimensions: one ggplot with a panel per pair
	p3 = uv_plot(fit, dims = 1:3) + ggplot2::ggtitle("x")
	expect_s3_class(p3, "ggplot")
	expect_equal(nlevels(p3$data$panel), 3L)
	expect_equal(p3$labels$title, "x")
	expect_length(ggplot2::ggplot_build(p3)$layout$layout$PANEL, 3L)
	expect_equal(nlevels(uv_plot(fit, dims = 1:3, layout = "biplot")$data$panel), 3L)
})

test_that("dims refuses values outside the fitted rank, repeats and more than four", {
	fit = uvd_fit()
	expect_error(uv_plot(fit, dims = c(1, 4)), "rank 3")
	expect_error(uv_plot(fit, dims = c(1, 1)), "distinct")
	expect_error(uv_plot(fit, dims = 1.5), "whole numbers")
	U6 = matrix(rnorm(60), 10, 6)
	expect_error(uv_plot(U = U6, V = U6, dims = 1:5), "at most four")
	expect_s3_class(uv_plot(U = U6, V = U6, dims = c(2, 5)), "ggplot")
})

test_that("a fit with one or two dimensions keeps its labels and has no caption", {
	fit = uvd_fit(n_c = 12, d = c(4, 2), bipartite = FALSE)
	p = uv_plot(fit, layout = "biplot")
	expect_equal(p$labels$x, "Latent dimension 1")
	expect_null(p$labels$caption)
	expect_null(uv_plot(fit)$labels$caption)
	expect_s3_class(uv_plot(uvd_fit(d = 3)), "ggplot")
})

test_that("matrix input names the dimensions without shares", {
	set.seed(3)
	U = matrix(rnorm(30), 10, 3)
	p = uv_plot(U = U, V = U, layout = "biplot", dims = c(1, 3))
	expect_equal(p$labels$x, "Dimension 1")
	expect_equal(p$labels$caption, "Dimensions 1 and 3 of 3 shown.")
})

test_that("the caption carries a tie with a dimension that is not shown", {
	st = data.frame(dimension = 1:3, strength = c(5, 3, 2.9), share = c(0.6, 0.2, 0.2), tied_with_next = c(FALSE, TRUE, FALSE))
	cap = lame:::.uv_dims_caption(1:2, 3L, st)
	expect_length(cap, 2L)
	expect_equal(cap[2], "Dimensions 2 and 3 have nearly equal strength, so their two axes are not well separated; read them together as a plane.")
	expect_length(lame:::.uv_dims_caption(2:3, 3L, st), 1L)
	expect_length(lame:::.uv_dims_caption(1:3, 3L, st), 1L)
})

test_that("dynamic fits take exactly two dims and no longer say which they use", {
	fit = uvd_dyn()
	expect_error(uv_plot(fit, dims = 1:3), "exactly two")
	expect_error(uv_plot(fit, dims = 2, plot_type = "trajectory"), "exactly two")
	expect_no_message(p <- uv_plot(fit, dims = c(1, 3), plot_type = "trajectory"))
	expect_equal(p$labels$y, "Dimension 3")
	expect_equal(p$labels$caption, "Dimensions 1 and 3 of 3 shown.")
	expect_s3_class(uv_plot(fit, dims = c(2, 3), plot_type = "faceted", nodes = "row"), "ggplot")
	expect_equal(unique(ggplot2::ggplot_build(uv_plot(fit, plot_type = "faceted", nodes = "col"))$plot$data$type), "Receiver")
})

test_that("nodes draws one side, and a symmetric fit refuses it", {
	fit = uvd_fit()
	p = uv_plot(fit, nodes = "row")
	expect_equal(unique(as.character(p$data$type)), "Row actor")
	expect_equal(nrow(p$data), 12L)
	p = uv_plot(fit, layout = "biplot", nodes = "col")
	expect_equal(unique(as.character(p$data$type)), "Column actor")
	sym = uvd_fit(n_c = 12, bipartite = FALSE)
	sym$symmetric = TRUE
	sym$UVPM = (sym$UVPM + t(sym$UVPM)) / 2
	sym$V = NULL
	expect_error(uv_plot(sym, nodes = "col"), "one set of positions")
	expect_error(uv_plot(fit, nodes = "row", show.edges = TRUE, Y = fit$UVPM > 0), "both sides")
})

test_that("a symmetric fit draws the similarity heatmap with nodes row or col", {
	sym = uvd_fit(n_c = 12, bipartite = FALSE)
	sym$symmetric = TRUE
	sym$UVPM = (sym$UVPM + t(sym$UVPM)) / 2
	sym$V = NULL
	p_r = uv_plot(sym, layout = "heatmap", nodes = "row")
	p_c = uv_plot(sym, layout = "heatmap", nodes = "col")
	expect_s3_class(p_r, "ggplot")
	expect_equal(p_r$data, p_c$data)
	expect_equal(nrow(p_r$data), 144L)
	expect_s3_class(uv_plot(sym, layout = "heatmap"), "ggplot")
	expect_error(uv_plot(sym, layout = "biplot", nodes = "row"), "layout = \"heatmap\"", fixed = TRUE)
})

test_that("a side over 50 nodes labels its 20 strongest; show.usernames wins", {
	fit = uvd_fit(n_r = 80, n_c = 10)
	p = uv_plot(fit, layout = "biplot")
	lab = p$data[p$data$label != "", ]
	expect_equal(sum(lab$type == "Row actor"), 20L)
	expect_equal(sum(lab$type == "Column actor"), 10L)
	b = ggplot2::ggplot_build(p)$data
	expect_true(any(vapply(b, nrow, 1L) == 60L))
	expect_match(p$labels$caption, "Row positions rescaled for display", fixed = TRUE)
	# the rescaled row cloud has the column cloud's rms radius
	rms = function(z) sqrt(mean(z$x^2 + z$y^2))
	ps = uv_plot(fit, layout = "biplot", dims = 1:2, vscale = 1)$data
	expect_equal(rms(ps[ps$type == "Row actor", ]), rms(ps[ps$type == "Column actor", ]), tolerance = 0.05)
	p_all = uv_plot(fit, show.usernames = TRUE)
	expect_equal(sum(p_all$data$label != ""), 90L)
	expect_false(any(grepl("rescaled", uv_plot(fit, layout = "biplot", nodes = "row")$labels$caption)))
})

test_that("the heatmap uses every dimension and refuses dims", {
	fit = uvd_fit()
	p = uv_plot(fit, layout = "heatmap", nodes = "col")
	expect_equal(nrow(p$data), 100L)
	b = ggplot2::ggplot_build(p)
	expect_equal(b$plot$scales$get_scales("fill")$name, "Similarity (all 3 dimensions)")
	pb = uv_plot(fit, layout = "heatmap")
	expect_equal(nrow(pb$data), 120L)
	expect_equal(range(pb$data$fill_value), range(fit$UVPM), tolerance = 1e-10)
	expect_error(uv_plot(fit, layout = "heatmap", dims = 1:2), "takes no")
	U = matrix(rnorm(2001 * 2), 2001, 2)
	V = matrix(rnorm(10 * 2), 10, 2)
	expect_error(uv_plot(U = U, V = V, layout = "heatmap", nodes = "row"), "at most 2,000")
	expect_equal(nrow(uv_plot(U = U, V = V, layout = "heatmap", nodes = "col")$data), 100L)
	dyn = uvd_dyn()
	expect_s3_class(uv_plot(dyn, layout = "heatmap", nodes = "row", time_point = 2), "ggplot")
	expect_error(uv_plot(dyn, layout = "heatmap", time_point = "average"), "one period")
})

test_that("latent_similarity is the cosine between the product's rows or columns", {
	fit = uvd_fit()
	P = fit$UVPM
	s = latent_similarity(fit, "col")
	Pn = sweep(P, 2, sqrt(colSums(P^2)), "/")
	expect_equal(unclass(s)[, ], crossprod(Pn), tolerance = 1e-8, ignore_attr = TRUE)
	expect_equal(dimnames(s), list(colnames(P), colnames(P)))
	expect_setequal(attr(s, "order"), seq_len(10))
	expect_identical(attr(s, "weak"), character(0))
	r = latent_similarity(fit, "row")
	Pr = P / sqrt(rowSums(P^2))
	expect_equal(unclass(r)[, ], tcrossprod(Pr), tolerance = 1e-8, ignore_attr = TRUE)
	expect_error(latent_similarity(fit), "required")
	expect_error(latent_similarity(fit, "col", time_point = 1), "dynamic positions")
	expect_equal(dim(latent_similarity(uvd_dyn(), "row", time_point = "t2")), c(8L, 8L))
})

test_that("nodes not distinguishable from zero are weak and placed last", {
	fit = lame:::.lame_positions(uvd_fit())
	fit$V[2, ] = 0.01
	fit$V_sd = matrix(0.1, nrow(fit$V), ncol(fit$V))
	s = latent_similarity(fit, "col")
	expect_identical(attr(s, "weak"), "c02")
	expect_equal(utils::tail(attr(s, "order"), 1L), 2L)
	p = uv_plot(fit, layout = "heatmap", nodes = "col")
	expect_match(p$labels$caption, "c02")
	expect_true(all(is.na(p$data$fill_value[p$data$x == 10])))
})
