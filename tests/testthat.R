library(testthat)
library(lame)

# cran runs the essential file only; the full suite runs with NOT_CRAN set to "true"
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
	test_check("lame")
} else {
	test_check("lame", filter = "cran_essential")
}
