# ab_plot(effect = "both") for ame / lame fits: sender and receiver effects
# as two panels, each sorted on its own. fake fit objects keep this fast; one
# real fit at the end checks the actual object layout.

fake_static = function(a, b, mode = "unipartite", cls = "ame") {
	fit = list(APM = a, BPM = b, mode = mode)
	class(fit) = cls
	fit
}

# the panel key is "<panel>\r<actor>"; strip the panel prefix
key_actor = function(keys) sub(".*\r", "", keys)

a_eff = c(x = 0.5, y = -1, z = 0.2)
b_eff = c(x = -0.3, y = 0.9, z = 0.1)

test_that("ab_plot both facets sender and receiver for a static ame fit", {
	fit = fake_static(a_eff, b_eff)
	p = ab_plot(fit, effect = "both")
	expect_s3_class(p, "ggplot")
	expect_equal(levels(p$data$margin),
	             c("Sender Effects (a)", "Receiver Effects (b)"))
	expect_equal(nrow(p$data), 6L)
	# each panel is ordered by its own values
	lv = levels(p$data$key)
	expect_equal(key_actor(lv[startsWith(lv, "Sender")]), names(sort(a_eff)))
	expect_equal(key_actor(lv[startsWith(lv, "Receiver")]), names(sort(b_eff)))
	expect_setequal(p$data$id, c("x", "y", "z"))
	expect_match(p$labels$title, "Sender and Receiver")
})

test_that("ab_plot both keeps input order when sorted = FALSE", {
	fit = fake_static(a_eff, b_eff)
	p = ab_plot(fit, effect = "both", sorted = FALSE)
	lv = levels(p$data$key)
	expect_equal(key_actor(lv[startsWith(lv, "Sender")]), names(a_eff))
	expect_equal(key_actor(lv[startsWith(lv, "Receiver")]), names(b_eff))
	expect_false(grepl("sorted", p$labels$title))
})

test_that("ab_plot both handles bipartite fits with different actor sets", {
	a = c(r1 = 0.1, r2 = -0.4, r3 = 0.8)
	b = c(c1 = 0.3, c2 = -0.2, c3 = 0.5, c4 = -0.9)
	fit = fake_static(a, b, mode = "bipartite")
	p = ab_plot(fit, effect = "both")
	expect_equal(levels(p$data$margin),
	             c("Row Actor Effects (a)", "Column Actor Effects (b)"))
	expect_equal(nrow(p$data), 7L)
	expect_match(p$labels$title, "Row and Column Actor")
})

test_that("ab_plot both aborts cleanly for symmetric fits", {
	fit = fake_static(c(x = 0.5, y = -1), NULL)
	expect_error(ab_plot(fit, effect = "both"), "no receiver")
	expect_error(ab_plot(fit, effect = "receiver"), "no receiver")
	expect_s3_class(ab_plot(fit, effect = "sender"), "ggplot")
})

test_that("ab_plot single-margin static plots are unchanged", {
	fit = fake_static(a_eff, b_eff)
	p = ab_plot(fit, effect = "sender")
	expect_false("margin" %in% names(p$data))
	expect_equal(levels(p$data$key), names(sort(a_eff)))
	expect_equal(p$labels$y, "Sender Effects (a)")
	expect_equal(p$labels$title, "Additive Sender Effects (sorted)")
	p2 = ab_plot(fit, effect = "receiver", sorted = FALSE, title = "custom")
	expect_equal(levels(p2$data$key), names(b_eff))
	expect_equal(p2$labels$title, "custom")
	expect_equal(p2$labels$y, "Receiver Effects (b)")
})

test_that("ab_plot both works for dynamic snapshots and refuses the other views", {
	am = matrix(c(0, 1, 2, 0, -1, -2), nrow = 2, byrow = TRUE,
	            dimnames = list(c("A", "B"), c("2001", "2002", "2003")))
	fit = list(a_dynamic = am, b_dynamic = -am, mode = "unipartite")
	class(fit) = c("lame", "ame")
	p = ab_plot(fit, effect = "both")
	expect_s3_class(p, "ggplot")
	expect_equal(nrow(p$data), 4L)
	expect_match(p$labels$title, "at 2003")
	expect_equal(p$data$mu[p$data$margin == "Sender Effects (a)"], c(2, -2))
	expect_equal(p$data$mu[p$data$margin == "Receiver Effects (b)"], c(-2, 2))
	p_avg = ab_plot(fit, effect = "both", time_point = "average")
	expect_match(p_avg$labels$title, "Average")
	expect_equal(p_avg$data$mu[p_avg$data$margin == "Sender Effects (a)"], c(1, -1))
	p1 = ab_plot(fit, effect = "both", time_point = 1)
	expect_equal(p1$data$mu, c(0, 0, 0, 0))
	expect_error(ab_plot(fit, effect = "both", plot_type = "trajectory"), "separately")
	expect_error(ab_plot(fit, effect = "both", plot_type = "ribbon"), "separately")
	expect_error(ab_plot(fit, effect = "both", time_point = "all"), "separately")
	expect_error(ab_plot(fit, effect = "both", time_point = 7), "between 1 and 3")
	# single-margin dynamic plots still go through the dynamic path
	p_traj = ab_plot(fit, effect = "sender", plot_type = "trajectory")
	expect_equal(p_traj$data$effect[p_traj$data$actor == "A"], c(0, 1, 2))
})

test_that("ab_plot both still delegates static ame_als fits to their own chart", {
	fit = list(a = a_eff, b = b_eff)
	class(fit) = "ame_als"
	p = ab_plot(fit, effect = "both")
	expect_s3_class(p, "ggplot")
	expect_setequal(unique(p$data$effect), c("sender", "receiver"))
})

test_that("ab_plot both runs on a real ame fit", {
	skip_on_cran()
	set.seed(3); n = 12
	Y = matrix(rnorm(n * n), n, n); diag(Y) = NA
	dimnames(Y) = list(paste0("a", 1:n), paste0("a", 1:n))
	fit = ame(Y, R = 0, burn = 5, nscan = 10, odens = 1, verbose = FALSE)
	p = ab_plot(fit, effect = "both")
	expect_s3_class(p, "ggplot")
	expect_equal(nrow(p$data), 2L * n)
	expect_setequal(p$data$id, paste0("a", 1:n))
})
