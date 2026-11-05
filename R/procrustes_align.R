#' Procrustes alignment of latent positions across time
#'
#' @description
#' Aligns dynamic latent positions across time periods to remove arbitrary
#' rotational indeterminacy. This is essential for interpreting temporal
#' trajectories of latent positions, since the latent space is only identified
#' up to rotation at each time point.
#'
#' Uses Procrustes rotation: at each time step, finds the orthogonal rotation
#' matrix that best aligns the current sender coordinates to the preceding
#' period, then applies it sequentially through the series.
#'
#' @param object A fitted \code{ame} or \code{lame} object, or NULL if raw
#'   arrays are provided via \code{U}/\code{V}/\code{G}.
#' @param U Optional 3D array \code{[n, R, T]} of sender latent positions.
#'   If \code{object} is provided and \code{U} is NULL, extracted from
#'   \code{object$U}.
#' @param V Optional 3D array \code{[n, R, T]} of receiver latent positions.
#'   For unipartite asymmetric models, U and V share one rotation per period so U V' is unchanged. For bipartite
#'   models, aligned jointly with U through the G interaction matrix.
#' @param G Optional \code{[R_row, R_col]} interaction matrix for bipartite
#'   models; it is returned unchanged, and each period's U G V' is preserved.
#' @param return_fit Logical. If TRUE and \code{object} is provided, returns
#'   a modified copy of the fit object with aligned latent positions.
#'   Default FALSE.
#' @param per_draw Logical. When TRUE, run Procrustes alignment per
#'   posterior draw rather than on the posterior-mean trajectory. This
#'   uses \code{object$U_draws} / \code{object$V_draws} on a dynamic fit,
#'   or \code{object$U_samples} on a static one (populated by
#'   \code{posterior_opts = list(save_UV_draws = TRUE)}) when present;
#'   otherwise falls back to mean-trajectory
#'   alignment and emits an informational note. On a bipartite fit each draw
#'   is aligned with its own interaction-matrix draw from
#'   \code{object$G_samples}, so every draw's U G V' is unchanged; without
#'   those draws per-draw alignment is an error. Per-draw alignment is the
#'   methodologically correct treatment -- rotation indeterminacy is a
#'   per-draw property, not a per-mean one -- but is more memory-hungry.
#'   Default FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return If \code{return_fit = FALSE} (default): a list with components
#'   \code{U} (aligned sender positions), \code{V} (aligned receiver positions,
#'   if applicable), and \code{G} (the interaction matrix, for bipartite).
#'
#'   If \code{return_fit = TRUE}: a copy of \code{object} with aligned latent
#'   positions replacing the originals.
#'
#' @details
#' For directed unipartite networks, U and V share one Procrustes rotation per
#' period, found from the stacked positions, so the product U V' is unchanged.
#' For symmetric networks the product is \eqn{U L U'} with a diagonal L, which
#' only a sign change of a column of U leaves unchanged; each period takes the
#' column signs that best match the period before it.
#'
#' For bipartite networks, one G serves every period, so a period's product
#' \code{U \%*\% G \%*\% t(V)} is unchanged only by transformations
#' \code{R_U}, \code{R_V} with \code{t(R_U) \%*\% G \%*\% R_V = G}. In the
#' singular-vector bases of G these are sign changes of a singular direction,
#' applied to U and V together; each period takes the sign changes that best
#' match the period before it, and G is returned unchanged.
#'
#' If U is a 2D matrix (static model with a single time point), it is returned
#' unchanged with an informational message.
#'
#' Alignment is within one fit: each period of a single dynamic trajectory is
#' rotated onto the period before it. Two different fits are never aligned to
#' each other automatically. Because raw arrays are accepted, you can do that
#' explicitly: stack the two fits' \code{U} matrices into an \code{[n, R, 2]}
#' array (same actors in the same order, same \code{R}) and call
#' \code{procrustes_align(U = that_array)}; the second slice comes back
#' rotated onto the first.
#'
#' @seealso \code{\link{latent_positions}} for extracting aligned positions as
#'   a tidy data frame, \code{\link{uv_plot}} for visualizing latent positions
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' # YX_bin_list$Y stores latent-scale values; threshold to 0/1 first
#' Y_bin <- lapply(YX_bin_list$Y, function(y) 1 * (y > 0))
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(Y_bin, Xdyad = YX_bin_list$X, R = 2,
#'             family = "binary", dynamic_uv = TRUE,
#'             burn = 5, nscan = 5, odens = 1,
#'             verbose = FALSE, chains = 1)
#' aligned <- procrustes_align(fit)
#' str(aligned$U)  # aligned 3D array [n, R, T]
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
procrustes_align <- function(
	object = NULL,
	U = NULL,
	V = NULL,
	G = NULL,
	return_fit = FALSE,
	per_draw = FALSE,
	...
) {

	# per-draw alignment: align each posterior draw's U/V trajectory separately
	if (isTRUE(per_draw)) {
		U_full <- object$U_full %||% object$U_draws %||% object$U_samples
		V_full <- object$V_full %||% object$V_draws
		if (is.null(U_full) || length(dim(U_full)) != 4L) {
			cli::cli_inform(c(
				"i" = "Per-draw alignment needs the 4-D per-draw trajectory cube {.code object$U_draws} ([n, R, T, S]); refit with {.code posterior_opts = list(save_UV_draws = TRUE)} to store one.",
				"i" = "Falling back to posterior-mean (mean-trajectory) alignment."
			))
			per_draw <- FALSE
		} else {
			S <- dim(U_full)[4L]
			U_aligned <- U_full
			V_aligned <- V_full
			# a bipartite draw keeps its own U G V', so it is aligned with its own G draw
			G_full <- object$G_samples
			bip_draws <- identical(object$mode, "bipartite") && length(dim(V_full)) == 4L
			if (bip_draws && (length(dim(G_full)) != 3L || dim(G_full)[3L] != S)) {
				cli::cli_abort(c(
					"Per-draw alignment of a bipartite fit needs one interaction-matrix draw per latent draw in {.code object$G_samples}.",
					"i" = "Refit with {.code posterior_opts = list(save_UV_draws = TRUE)}, or use {.code per_draw = FALSE}."
				))
			}
			for (s in seq_len(S)) {
				if (bip_draws) {
					a <- align_over_time_bip(
						array(U_full[, , , s], dim(U_full)[1:3]),
						array(V_full[, , , s], dim(V_full)[1:3]),
						matrix(G_full[, , s], dim(G_full)[1L])
					)
					U_aligned[, , , s] <- a$U
					V_aligned[, , , s] <- a[["V"]]
				} else if (length(dim(V_full)) == 4L) {
					a <- procrustes_align(
						NULL,
						U = U_full[, , , s],
						V = V_full[, , , s],
						G = G
					)
					U_aligned[, , , s] <- a$U
					V_aligned[, , , s] <- a[["V"]]
				} else if (isTRUE(object$symmetric)) {
					U_aligned[, , , s] <- align_over_time_sym(array(U_full[, , , s], dim(U_full)[1:3]))
				} else {
					a <- procrustes_align(NULL, U = U_full[, , , s])
					U_aligned[, , , s] <- a$U
				}
			}
			result <- list(
				U = U_aligned,
				V = V_aligned,
				G = G,
				per_draw = TRUE,
				n_draws = S
			)
			if (return_fit && !is.null(object)) {
				fit_copy <- object
				fit_copy$U_full_aligned <- U_aligned
				fit_copy$V_full_aligned <- V_aligned
				# fitted U and V become the posterior means of the aligned cube
				fit_copy$U <- apply(U_aligned, c(1, 2, 3), mean)
				if (length(dim(V_full)) == 4L) {
					fit_copy[["V"]] <- apply(V_aligned, c(1, 2, 3), mean)
				}
				return(fit_copy)
			}
			return(result)
		}
	}

	# take positions from the fit; exact [[ ]] access because object$G would
	# partial-match GOF
	if (!is.null(object)) {
		if (is.null(U)) U <- object[["U"]]
		if (is.null(V)) V <- object[["V"]]
		if (is.null(G)) G <- object[["G"]]
	}

	if (is.null(U)) {
		cli::cli_abort("No latent positions found. Either provide a fitted object with R > 0 or pass U directly.")
	}

	# static model has nothing to align
	if (is.matrix(U) || length(dim(U)) < 3L) {
		cli::cli_alert_info("Latent positions are static (single time point). No alignment needed.")
		result <- list(U = U, V = V, G = G)
		if (return_fit && !is.null(object)) return(object)
		return(result)
	}

	is_bip <- FALSE
	if (!is.null(object) && !is.null(object$mode)) {
		is_bip <- identical(object$mode, "bipartite")
	} else if (!is.null(V) && !is.null(U) && length(dim(V)) == 3L) {
		# infer bipartite from differing row counts
		is_bip <- nrow(U) != nrow(V)
	}

	if (is_bip && !is.null(V) && !is.null(G)) {
		# align row positions, column positions and G jointly
		aligned <- align_over_time_bip(U, V, G)
		U_al <- aligned$U
		V_al <- aligned[["V"]]
		G_al <- aligned$G
	} else {
		# directed fits: one rotation per period for U and V together, so U V' and
		# every fitted value are unchanged; symmetric fits have U only
		if (!is.null(V) && length(dim(V)) == 3L) {
			pair <- align_over_time_unip_pair(U, V)
			U_al <- pair$U
			V_al <- pair[["V"]]
		} else if (!is.null(object) && isTRUE(object$symmetric)) {
			U_al <- align_over_time_sym(U)
			V_al <- V
		} else {
			U_al <- align_over_time_unip(U)
			V_al <- V
		}
		G_al <- G
	}

	if (return_fit && !is.null(object)) {
		fit_copy <- object
		fit_copy$U <- U_al
		fit_copy[["V"]] <- V_al
		if (!is.null(G_al)) fit_copy$G <- G_al
		return(fit_copy)
	}

	list(U = U_al, V = V_al, G = G_al)
}
