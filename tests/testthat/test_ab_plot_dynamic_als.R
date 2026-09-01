test_that("ab_plot honors dynamic ALS trajectories", {
	fit = list(
		a = c(a1 = 0, a2 = 0),
		b = c(a1 = 0, a2 = 0),
		a_dynamic = matrix(
			c(0, 1, 2, 0, -1, -2),
			nrow = 2,
			byrow = TRUE,
			dimnames = list(c("State A", "State B"), c("2001", "2002", "2003"))
		),
		b_dynamic = matrix(
			c(0, -1, -2, 0, 1, 2),
			nrow = 2,
			byrow = TRUE,
			dimnames = list(c("State A", "State B"), c("2001", "2002", "2003"))
		)
	)
	class(fit) = c("lame_dynamic_als", "lame_als", "ame_als")

	p = ab_plot(
		fit,
		effect = "sender",
		plot_type = "trajectory",
		show_actors = "State A"
	)

	expect_s3_class(p, "ggplot")
	expect_equal(unique(p$data$actor), "State A")
	expect_equal(p$data$time_label, c("2001", "2002", "2003"))
	expect_equal(p$data$effect, c(0, 1, 2))
})
