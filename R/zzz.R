#' @importFrom grDevices devAskNewPage
#' @importFrom stats aggregate dnorm rbinom rchisq rWishart setNames
#' @importFrom utils head object.size
NULL

# register S3 methods (posterior, loo, ggplot2::autoplot, generics::tidy) at load so
# package::generic(fit) dispatches to them
.onLoad <- function(libname, pkgname) {
	.register_posterior_methods()
	.register_loo_methods()
	.register_autoplot_method()
	.register_tidy_methods()
}

# warn at attach when amen is also loaded: both export ame()
.onAttach <- function(libname, pkgname) {
	if ("amen" %in% loadedNamespaces()) {
		packageStartupMessage(
			"`lame` and `amen` both export `ame()`. ",
			"Unqualified calls dispatch to whichever package was attached last. ",
			"Use `lame::ame(...)` or `amen::ame(...)` to be explicit. ",
			"See `?lame::ame` for migration notes from amen."
		)
	}
}
