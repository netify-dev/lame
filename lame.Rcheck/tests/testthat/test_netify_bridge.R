test_that("netify cross-sectional inputs work with ame and ame_als", {
	skip_if_not_installed("netify", minimum_version = "1.5.3")
	set.seed(11)
	df = expand.grid(
		i = paste0("a", seq_len(5)),
		j = paste0("a", seq_len(5)),
		stringsAsFactors = FALSE
	)
	df = df[df$i != df$j, , drop = FALSE]
	df$y <- rbinom(nrow(df), 1, 0.35)

	net = netify::netify(
		df,
		actor1 = "i", actor2 = "j",
		weight = "y",
		symmetric = FALSE,
		mode = "unipartite",
		missing_to_zero = FALSE
	)

	fit_als = suppressWarnings(ame_als(net, R = 0, max_iter = 20, verbose = FALSE))
	expect_s3_class(fit_als, "ame_als")
	expect_equal(fit_als$family, "binary")
	expect_equal(fit_als$mode, "unipartite")
	expect_equal(dim(fit_als$Y), c(5L, 5L, 1L))

	fit_mcmc = ame(
		net, R = 0,
		nscan = 12, burn = 6, odens = 3,
		gof = FALSE, verbose = FALSE
	)
	expect_s3_class(fit_mcmc, "ame")
	expect_equal(fit_mcmc$family, "binary")
	expect_equal(fit_mcmc$mode, "unipartite")
})

test_that("netify bridge rejects mismatched entry points and duplicate covariates", {
	skip_if_not_installed("netify", minimum_version = "1.5.3")
	df = data.frame(
		i = c("a", "b", "a", "b"),
		j = c("b", "a", "b", "a"),
		year = c(2001, 2001, 2002, 2002),
		y = c(1, 0, 0, 1)
	)
	net_long = netify::netify(
		df,
		actor1 = "i", actor2 = "j", time = "year",
		weight = "y", symmetric = FALSE,
		output_format = "longit_list",
		missing_to_zero = FALSE
	)
	expect_error(ame(net_long), "Use .*lame")

	net_xs = netify::netify(
		df[df$year == 2001, ],
		actor1 = "i", actor2 = "j",
		weight = "y", symmetric = FALSE,
		missing_to_zero = FALSE
	)
	expect_error(
		ame_als(net_xs, Xdyad = array(0, c(2, 2, 1))),
		"pass covariates through"
	)
})

.netify_bipartite_panel = function(add_nodes = TRUE) {
	df = data.frame(
		row = c("r1", "r2", "r1", "r3", "r2", "r3"),
		col = c("c1", "c1", "c2", "c2", "c3", "c3"),
		year = c(2001, 2001, 2002, 2002, 2003, 2003),
		y = c(1, 0, 2, 1, 0, 3),
		z = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
	)
	roster = data.frame(
		actor = c("r1", "r2", "r3", "c1", "c2", "c3"),
		min_time = c(2001, 2001, 2002, 2001, 2002, 2003),
		max_time = c(2002, 2003, 2003, 2001, 2002, 2003)
	)
	net = netify::netify(
		df,
		actor1 = "row", actor2 = "col", time = "year",
		weight = "y", dyad_vars = "z",
		mode = "bipartite", symmetric = FALSE,
		output_format = "longit_list",
		actor_time_uniform = FALSE,
		actor_pds = roster,
		missing_to_zero = FALSE
	)
	if (isTRUE(add_nodes)) {
		node_df = data.frame(
			actor = rep(c("r1", "r2", "r3", "c1", "c2", "c3"), each = 3),
			time = rep(2001:2003, 6),
			score = seq_len(18)
		)
		net = netify::add_node_vars(net, node_df, actor = "actor", time = "time")
	}
	net
}

test_that("netify to_lame ragged bipartite output matches lame padding", {
	skip_if_not_installed("netify", minimum_version = "1.5.3")
	net = .netify_bipartite_panel(add_nodes = TRUE)
	nl = netify::to_lame(net, lame = TRUE, family = "normal", pad = FALSE)

	row_actors = sort(unique(unlist(lapply(nl$Y, rownames), use.names = FALSE)))
	col_actors = sort(unique(unlist(lapply(nl$Y, colnames), use.names = FALSE)))
	padded = list_to_array_bipartite(
		row_actors, col_actors,
		nl$Y, nl$Xdyad, nl$Xrow, nl$Xcol
	)

	expect_equal(dim(padded$Y), c(3L, 3L, 3L))
	expect_equal(dimnames(padded$Y), list(row_actors, col_actors, names(nl$Y)))
	expect_equal(padded$Y["r1", "c1", "2001"], 1)
	expect_true(is.na(padded$Y["r3", "c1", "2001"]))
	expect_true(is.na(padded$Y["r1", "c3", "2001"]))

	expect_equal(dim(padded$Xrow), c(3L, 1L, 3L))
	expect_equal(dim(padded$Xcol), c(3L, 1L, 3L))
	expect_equal(dimnames(padded$Xrow)[[1]], row_actors)
	expect_equal(dimnames(padded$Xcol)[[1]], col_actors)
	expect_equal(padded$Xrow["r1", "score_row", "2001"], 1)
	expect_equal(padded$Xcol["c1", "score_col", "2001"], 10)
})

test_that("ragged bipartite netify panels fit through lame entry points", {
	skip_if_not_installed("netify", minimum_version = "1.5.3")
	net = .netify_bipartite_panel(add_nodes = FALSE)

	fit_als = suppressWarnings(lame_als(
		net, R = 0, family = "normal",
		max_iter = 20, verbose = FALSE
	))
	expect_s3_class(fit_als, "ame_als")
	expect_equal(fit_als$mode, "bipartite")
	expect_equal(dim(fit_als$Y), c(3L, 3L, 3L))

	fit_dispatch = suppressWarnings(lame(
		net, R = 0, family = "normal", method = "als",
		als_max_iter = 20, verbose = FALSE
	))
	expect_s3_class(fit_dispatch, "ame_als")
	expect_equal(fit_dispatch$mode, "bipartite")
	expect_equal(dim(fit_dispatch$Y), c(3L, 3L, 3L))

	fit_mcmc = suppressWarnings(lame(
		net, R = 0, family = "normal",
		nscan = 12, burn = 6, odens = 3,
		gof = FALSE, verbose = FALSE
	))
	expect_s3_class(fit_mcmc, "lame")
	expect_equal(fit_mcmc$mode, "bipartite")
	expect_equal(dim(fit_mcmc$Y), c(3L, 3L, 3L))
})

test_that("snap ALS accepts normal ragged bipartite netify panels", {
	skip_if_not_installed("netify", minimum_version = "1.5.3")
	set.seed(12)
	df = expand.grid(
		row = paste0("r", seq_len(4)),
		col = paste0("c", seq_len(4)),
		year = 2001:2003,
		stringsAsFactors = FALSE
	)
	df = df[!(df$row == "r4" & df$year == 2001) &
		!(df$col == "c4" & df$year < 2003), , drop = FALSE]
	df$y <- rnorm(nrow(df))
	roster = data.frame(
		actor = c(paste0("r", seq_len(4)), paste0("c", seq_len(4))),
		min_time = c(2001, 2001, 2001, 2002, 2001, 2001, 2001, 2003),
		max_time = 2003
	)
	net = netify::netify(
		df,
		actor1 = "row", actor2 = "col", time = "year",
		weight = "y", mode = "bipartite", symmetric = FALSE,
		output_format = "longit_list",
		actor_time_uniform = FALSE,
		actor_pds = roster,
		missing_to_zero = FALSE
	)

	fit = suppressWarnings(lame_snap_als(
		net, R = 1, family = "normal",
		max_iter = 3, stability = "none",
		verbose = FALSE
	))
	expect_s3_class(fit, "lame_snap_als")
	expect_equal(fit$mode, "bipartite")
	expect_equal(dim(fit$Y), c(4L, 4L, 3L))
})
