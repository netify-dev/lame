# one table has no slices to align, so names on one margin only raise no warning about slices; a panel whose
# slices differ in having names keeps the warning

test_that("a single two-mode table with column names only is fitted without a warning", {
	skip_on_cran()
	set.seed(4101)
	Y = matrix(rbinom(60 * 8, 1, 0.4), 60, dimnames = list(NULL, sprintf("s%d", 1:8)))
	expect_no_warning(suppressMessages(ame_als(Y, mode = "bipartite", family = "binary", R = 0, verbose = FALSE)))
})

test_that("a panel whose slices differ in having actor names warns that matching is by position", {
	skip_on_cran()
	set.seed(4102)
	nm = sprintf("a%02d", 1:15)
	Y1 = matrix(rnorm(225), 15, dimnames = list(nm, nm)); Y2 = matrix(rnorm(225), 15)
	diag(Y1) = NA; diag(Y2) = NA
	expect_warning(suppressMessages(lame_als(list(Y1, Y2), family = "normal", R = 0, verbose = FALSE)), "actor names")
})
