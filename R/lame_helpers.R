####
#' Convert a graph object to a lame-ready adjacency matrix
#'
#' Convenience wrapper that takes an \pkg{igraph} or \pkg{network} object and
#' returns the kind of matrix that \code{\link{ame}} / \code{\link{lame}}
#' expect: a numeric adjacency matrix with \code{NA} on the diagonal (self-ties
#' not modelled) and (if available) actor names preserved on the row/column
#' dimnames.
#'
#' Plain matrices and data.frames are accepted too; data.frames are coerced
#' to a numeric matrix and only the diagonal is rewritten to \code{NA}.
#'
#' @param x an \code{igraph}, \code{network}, \code{matrix}, or
#'   \code{data.frame} representing an adjacency / sociomatrix.
#' @param na_diag logical: replace the diagonal with \code{NA}? Defaults to
#'   \code{TRUE} for square (unipartite) matrices and is ignored for
#'   rectangular (bipartite) matrices.
#' @return A numeric matrix suitable to pass as \code{Y} to
#'   \code{\link{ame}} / \code{\link{lame}}.
#' @examples
#' \donttest{
#' # mat <- as_lame_y(igraph::sample_smallworld(1, 12, 2, 0.1))
#' # ame(mat, family = "binary", R = 1, burn = 50, nscan = 200,
#' #     odens = 5, verbose = FALSE, symmetric = TRUE)
#' }
#' @export
as_lame_y <- function(x, na_diag = TRUE) {
	if (inherits(x, "igraph")) {
		if (!requireNamespace("igraph", quietly = TRUE)) {
			cli::cli_abort("Package {.pkg igraph} is required to convert an {.cls igraph} object.")
		}
		m <- as.matrix(igraph::as_adjacency_matrix(x, sparse = FALSE))
	} else if (inherits(x, "network")) {
		if (!requireNamespace("network", quietly = TRUE)) {
			cli::cli_abort("Package {.pkg network} is required to convert a {.cls network} object.")
		}
		m <- as.matrix(x)
	} else if (is.data.frame(x)) {
		m <- as.matrix(x)
	} else if (is.matrix(x)) {
		m <- x
	} else {
		cli::cli_abort(c(
			"{.arg x} must be an {.cls igraph}, {.cls network}, {.cls matrix}, or {.cls data.frame}.",
			"i" = "Got {.cls {class(x)[1]}}."))
	}
	storage.mode(m) <- "double"
	if (na_diag && nrow(m) == ncol(m)) diag(m) <- NA_real_
	m
}

####
#' ERGM-style covariate helpers for ame() / lame()
#'
#' \code{nodematch(x)} returns an \code{n x n} matrix with \code{1} where
#' \code{x[i] == x[j]} and \code{0} otherwise -- the AME analogue of ERGM's
#' \code{nodematch} term, suitable for \code{Xdyad} as a single
#' homophily covariate. \code{absdiff(x)} returns the absolute difference
#' \code{|x[i] - x[j]|}, the AME analogue of ERGM's \code{absdiff}.
#' \code{nodefactor(x)} returns the dyadic matrix of \code{x[i] + x[j]} for
#' a numeric \code{x} (or a list of per-level binary dyadic indicators for
#' a factor / character \code{x}).
#'
#' All helpers return a matrix \emph{or} list of matrices that is ready to
#' wrap into an \code{n x n x p} \code{Xdyad} array via \code{simplify2array}
#' or \code{array()}.
#'
#' @param x a numeric, factor, or character vector of length \code{n}.
#' @param na_diag logical: set the diagonal to \code{NA}? Defaults to
#'   \code{TRUE}.
#' @return An \code{n x n} numeric matrix (or, for a factor \code{x} passed
#'   to \code{nodefactor}, a named list of such matrices).
#' @examples
#' \donttest{
#' n <- 12
#' grp <- sample(letters[1:3], n, replace = TRUE)
#' age <- rnorm(n, 40, 10)
#' # build a single same-group homophily covariate
#' Xdyad <- array(nodematch(grp), dim = c(n, n, 1),
#'                dimnames = list(NULL, NULL, "same_group"))
#' # combine homophily + age difference
#' Xdyad <- array(c(nodematch(grp), absdiff(age)),
#'                dim = c(n, n, 2),
#'                dimnames = list(NULL, NULL, c("same_group", "age_diff")))
#' }
#' @export
nodematch <- function(x, na_diag = TRUE) {
	if (length(x) < 2L) cli::cli_abort("{.arg x} must have length >= 2.")
	m <- outer(x, x, FUN = function(a, b) as.integer(a == b))
	storage.mode(m) <- "double"
	if (na_diag) diag(m) <- NA_real_
	m
}

#' @rdname nodematch
#' @export
absdiff <- function(x, na_diag = TRUE) {
	if (!is.numeric(x)) cli::cli_abort("{.arg x} must be numeric for {.fn absdiff}.")
	m <- abs(outer(x, x, "-"))
	if (na_diag) diag(m) <- NA_real_
	m
}

#' @rdname nodematch
#' @export
nodefactor <- function(x, na_diag = TRUE) {
	if (is.numeric(x)) {
		m <- outer(x, x, "+")
		if (na_diag) diag(m) <- NA_real_
		return(m)
	}
	if (is.factor(x) || is.character(x)) {
		levels_x <- if (is.factor(x)) levels(x) else sort(unique(x))
		out <- lapply(levels_x, function(lv) {
			ind <- as.integer(x == lv)
			m <- outer(ind, ind, "+")
			storage.mode(m) <- "double"
			if (na_diag) diag(m) <- NA_real_
			m
		})
		names(out) <- levels_x
		return(out)
	}
	cli::cli_abort("{.arg x} must be numeric, factor, or character for {.fn nodefactor}.")
}
