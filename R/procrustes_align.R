#' Procrustes alignment of latent positions across time
#'
#' @description
#' Aligns dynamic latent positions across time periods to remove arbitrary
#' rotational indeterminacy. This is essential for interpreting temporal
#' trajectories of latent positions, since the latent space is only identified
#' up to rotation at each time point.
#'
#' Uses Procrustes rotation: at each time step t, finds the orthogonal rotation
#' matrix that best aligns U_t to U_{t-1} (minimizing Frobenius norm), then
#' applies it sequentially from t=2 to T.
#'
#' @param object A fitted \code{ame} or \code{lame} object, or NULL if raw
#'   arrays are provided via \code{U}/\code{V}/\code{G}.
#' @param U Optional 3D array \code{[n, R, T]} of sender latent positions.
#'   If \code{object} is provided and \code{U} is NULL, extracted from
#'   \code{object$U}.
#' @param V Optional 3D array \code{[n, R, T]} of receiver latent positions.
#'   For unipartite asymmetric models, aligned independently. For bipartite
#'   models, aligned jointly with U via the G interaction matrix.
#' @param G Optional \code{[R_row, R_col]} interaction matrix for bipartite
#'   models. Updated to maintain the invariant U G V'.
#' @param return_fit Logical. If TRUE and \code{object} is provided, returns
#'   a modified copy of the fit object with aligned latent positions.
#'   Default FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return If \code{return_fit = FALSE} (default): a list with components
#'   \code{U} (aligned sender positions), \code{V} (aligned receiver positions,
#'   if applicable), and \code{G} (updated interaction matrix, for bipartite).
#'
#'   If \code{return_fit = TRUE}: a copy of \code{object} with aligned latent
#'   positions replacing the originals.
#'
#' @details
#' For unipartite networks, U and V are aligned independently using separate
#' Procrustes rotations. For symmetric networks, only U is present and aligned.
#'
#' For bipartite networks, U and V are aligned jointly: separate rotation
#' matrices are computed for U and V, and the G interaction matrix is updated
#' as \code{G_aligned = t(R_U) \%*\% G \%*\% R_V} to preserve the product
#' \code{U \%*\% G \%*\% t(V)}.
#'
#' If U is a 2D matrix (static model with a single time point), it is returned
#' unchanged with an informational message.
#'
#' @seealso \code{\link{latent_positions}} for extracting aligned positions as
#'   a tidy data frame, \code{\link{uv_plot}} for visualizing latent positions
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit <- lame(YX_bin_list$Y, Xdyad = YX_bin_list$X, R = 2,
#'             family = "binary", dynamic_uv = TRUE,
#'             burn = 5, nscan = 5, odens = 1,
#'             verbose = FALSE)
#' aligned <- procrustes_align(fit)
#' str(aligned$U)  # aligned 3D array [n, R, T]
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
procrustes_align <- function(object = NULL, U = NULL, V = NULL, G = NULL,
                              return_fit = FALSE, ...) {

	# extract from fit object if not provided directly
	if (!is.null(object)) {
		if (is.null(U)) U <- object$U
		if (is.null(V)) V <- object$V
		if (is.null(G)) G <- object$G
	}

	if (is.null(U)) {
		stop("No latent positions found. Either provide a fitted object with R > 0 or pass U directly.")
	}

	# static model (2D matrix) - nothing to align
	if (is.matrix(U) || length(dim(U)) < 3L) {
		cli::cli_alert_info("Latent positions are static (single time point). No alignment needed.")
		result <- list(U = U, V = V, G = G)
		if (return_fit && !is.null(object)) return(object)
		return(result)
	}

	# detect bipartite
	is_bip <- FALSE
	if (!is.null(object) && !is.null(object$mode)) {
		is_bip <- identical(object$mode, "bipartite")
	} else if (!is.null(V) && !is.null(U) && length(dim(V)) == 3L) {
		# infer bipartite from differing row counts
		is_bip <- nrow(U) != nrow(V)
	}

	if (is_bip && !is.null(V) && !is.null(G)) {
		# bipartite: joint alignment of U, V, G
		aligned <- align_over_time_bip(U, V, G)
		U_al <- aligned$U
		V_al <- aligned$V
		G_al <- aligned$G
	} else {
		# unipartite: align U (and V independently if present)
		U_al <- align_over_time_unip(U)
		V_al <- if (!is.null(V) && length(dim(V)) == 3L) {
			align_over_time_unip(V)
		} else {
			V
		}
		G_al <- G
	}

	if (return_fit && !is.null(object)) {
		fit_copy <- object
		fit_copy$U <- U_al
		fit_copy$V <- V_al
		if (!is.null(G_al)) fit_copy$G <- G_al
		return(fit_copy)
	}

	list(U = U_al, V = V_al, G = G_al)
}
