#' @importFrom grDevices devAskNewPage
#' @importFrom stats aggregate dnorm rbinom rWishart setNames
#' @importFrom utils head object.size
NULL

# register s3 methods (posterior / loo / ggplot2::autoplot / generics::tidy)
# at package load so qualified package::generic(fit) dispatches to the
# methods defined here. fires once per session.
.onLoad <- function(libname, pkgname) {
	.register_posterior_methods()
	.register_loo_methods()
	.register_autoplot_method()
	.register_tidy_methods()
}

# warn at attach time when amen is also loaded. both packages export
# ame() and unqualified calls dispatch to whichever was attached last.
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
