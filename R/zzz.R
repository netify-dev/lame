#' @importFrom grDevices devAskNewPage
#' @importFrom stats aggregate dnorm rbinom rchisq rWishart setNames
#' @importFrom utils head object.size
NULL

# methods for posterior, loo, generics and ggplot2 generics are registered by
# NAMESPACE when those packages load, so loading lame loads none of them
.onLoad <- function(libname, pkgname) {
	# exported autoplot is ggplot2's generic, bound lazily so loading lame
	# does not load ggplot2
	delayedAssign(
		"autoplot",
		ggplot2::autoplot,
		eval.env = baseenv(),
		assign.env = asNamespace(pkgname)
	)
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
