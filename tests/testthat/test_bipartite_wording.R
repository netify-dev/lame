# bipartite output speaks of rows and columns, and prints only the variance components
# the model has and estimates

wording_fit = function(Y, ...) {
	suppressWarnings(ame(Y, R = 0, burn = 20, nscan = 100, odens = 5, verbose = FALSE, chains = 1, seed = 1, ...))
}

printed = function(x) paste(utils::capture.output(print(x)), collapse = "\n")

# printed text on one line, for legends that wrap
flat = function(txt) gsub("\\s+", " ", txt)

test_that("gof labels name rows and columns for a bipartite network", {
	expect_equal(
		unname(gof_label(c("sd.row", "sd.colmean", "four.cycles"), bipartite = TRUE)),
		c("Row degree heterogeneity", "Column degree heterogeneity", "Four cycles")
	)
	expect_equal(
		unname(gof_label(c("sd.rowmean", "sd.col"))),
		c("Sender degree heterogeneity", "Receiver degree heterogeneity")
	)
})

test_that("a bipartite summary prints the row and column variances only", {
	skip_on_cran()
	set.seed(2)
	Yb = matrix(rbinom(40 * 8, 1, 0.4), 40, 8, dimnames = list(paste0("p", 1:40), paste0("s", 1:8)))
	fit = wording_fit(Yb, mode = "bipartite", family = "binary")
	txt = printed(summary(fit))
	expect_match(flat(txt), "va = row variance, vb = column variance", fixed = TRUE)
	expect_no_match(txt, "sender|receiver|fixed to 0|\\ncab |\\nrho |\\nve ")
	# the summary object keeps every component
	expect_equal(rownames(summary(fit)$variance), colnames(fit$VC))
	expect_match(printed(fit), "Row effects: enabled", fixed = TRUE)
	Yn = Yb + matrix(rnorm(length(Yb)), nrow(Yb))
	txt_n = printed(summary(wording_fit(Yn, mode = "bipartite", family = "normal")))
	expect_match(flat(txt_n), "va = row variance, vb = column variance, ve = residual variance", fixed = TRUE)
	expect_match(txt_n, "\nve ")
})

test_that("a one-mode normal summary keeps every component and its legend", {
	skip_on_cran()
	set.seed(3)
	Y = matrix(rnorm(15 * 15), 15, 15)
	diag(Y) = NA
	txt = printed(summary(wording_fit(Y)))
	expect_match(
		flat(txt),
		"va = sender variance, cab = sender-receiver covariance, vb = receiver variance, rho = dyadic correlation, ve = residual variance",
		fixed = TRUE
	)
	expect_match(txt, "\ncab ")
	expect_match(txt, "\nrho ")
	expect_match(txt, "\nve ")
})
