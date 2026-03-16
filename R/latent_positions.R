#' Extract latent positions as a tidy data frame
#'
#' @description
#' Extracts multiplicative latent factor positions (U and V) from a fitted
#' \code{ame} or \code{lame} model and returns them as a tidy data frame
#' suitable for plotting and analysis. Optionally applies Procrustes alignment
#' for dynamic models and includes posterior standard deviations when
#' posterior samples are available.
#'
#' @param object A fitted \code{ame} or \code{lame} model object with R > 0.
#' @param align Logical. For dynamic models (\code{dynamic_uv = TRUE}), apply
#'   Procrustes alignment across time to remove rotational indeterminacy.
#'   Default is \code{FALSE} for \code{ame} objects and \code{TRUE} for
#'   \code{lame} objects.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{actor}{Character. Actor name (from rownames of U or V).}
#'   \item{dimension}{Integer. Latent dimension index (1 to R).}
#'   \item{time}{Character or NA. Time period label for dynamic models;
#'     \code{NA} for static (cross-sectional) models.}
#'   \item{value}{Numeric. The posterior mean latent position.}
#'   \item{posterior_sd}{Numeric. Posterior standard deviation of the latent
#'     position, or \code{NA} if posterior samples are not available.
#'     To enable, fit the model with
#'     \code{posterior_opts = posterior_options(save_UV = TRUE)}.}
#'   \item{type}{Character. \code{"U"} for sender/row positions,
#'     \code{"V"} for receiver/column positions. Symmetric models have
#'     only \code{"U"}.}
#' }
#'
#' Returns a zero-row data frame with correct column names if R = 0.
#'
#' @seealso \code{\link{procrustes_align}} for standalone Procrustes alignment,
#'   \code{\link{uv_plot}} for visualizing latent positions,
#'   \code{\link{posterior_options}} for enabling posterior sampling of U/V
#'
#' @examples
#' \donttest{
#' data(YX_nrm)
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            burn = 5, nscan = 5, odens = 1, verbose = FALSE)
#' lp <- latent_positions(fit)
#' head(lp)
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
latent_positions <- function(object, ...) {
	UseMethod("latent_positions")
}

#' @rdname latent_positions
#' @method latent_positions ame
#' @export
latent_positions.ame <- function(object, align = FALSE, ...) {
	.extract_latent_positions(object, align = align)
}

#' @rdname latent_positions
#' @method latent_positions lame
#' @export
latent_positions.lame <- function(object, align = TRUE, ...) {
	.extract_latent_positions(object, align = align)
}

# internal workhorse
.extract_latent_positions <- function(object, align = FALSE) {

	empty_df <- data.frame(
		actor = character(0),
		dimension = integer(0),
		time = character(0),
		value = numeric(0),
		posterior_sd = numeric(0),
		type = character(0),
		stringsAsFactors = FALSE
	)

	U <- object$U
	V <- object$V
	G <- object$G

	if (is.null(U)) return(empty_df)

	is_symmetric <- isTRUE(object$symmetric)
	is_bip <- identical(object$mode, "bipartite")

	# apply Procrustes alignment if requested
	if (align && length(dim(U)) == 3L) {
		aligned <- procrustes_align(object = object)
		U <- aligned$U
		V <- aligned$V
		G <- aligned$G
	}

	# compute posterior SDs from samples if available
	U_sd <- .compute_posterior_sd(object$U_samples)
	V_sd <- .compute_posterior_sd(object$V_samples)

	# build data frames from arrays
	u_df <- .array_to_df(U, type = "U", sd_array = U_sd)

	if (is_symmetric || is.null(V)) {
		return(u_df)
	}

	v_df <- .array_to_df(V, type = "V", sd_array = V_sd)
	rbind(u_df, v_df)
}

# compute element-wise posterior SD from a 3D samples array [n, R, n_samples]
.compute_posterior_sd <- function(samples) {
	if (is.null(samples)) return(NULL)
	if (length(dim(samples)) != 3L) return(NULL)
	# apply SD across the samples dimension (3rd dim)
	apply(samples, c(1, 2), sd, na.rm = TRUE)
}

# convert a 2D matrix or 3D array to tidy data frame
# sd_array: optional matrix [n, R] of posterior SDs (same shape as 2D x)
.array_to_df <- function(x, type, sd_array = NULL) {
	if (is.matrix(x) || length(dim(x)) == 2L) {
		# static: [n, R]
		actors <- rownames(x)
		if (is.null(actors)) actors <- paste0("node", seq_len(nrow(x)))
		R <- ncol(x)
		dfs <- lapply(seq_len(R), function(r) {
			sd_vals <- if (!is.null(sd_array) && ncol(sd_array) >= r) {
				sd_array[, r]
			} else {
				rep(NA_real_, nrow(x))
			}
			data.frame(
				actor = actors,
				dimension = r,
				time = NA_character_,
				value = x[, r],
				posterior_sd = sd_vals,
				type = type,
				stringsAsFactors = FALSE,
				row.names = NULL
			)
		})
		do.call(rbind, dfs)
	} else if (length(dim(x)) == 3L) {
		# dynamic: [n, R, T]
		actors <- dimnames(x)[[1]]
		if (is.null(actors)) actors <- paste0("node", seq_len(dim(x)[1]))
		time_labels <- dimnames(x)[[3]]
		if (is.null(time_labels)) time_labels <- as.character(seq_len(dim(x)[3]))
		R <- dim(x)[2]
		TT <- dim(x)[3]
		dfs <- vector("list", R * TT)
		idx <- 1L
		for (r in seq_len(R)) {
			for (t in seq_len(TT)) {
				# For dynamic models, posterior_sd from samples is per-actor per-dim
				# (not per-time), so use the same SD for all time points
				sd_vals <- if (!is.null(sd_array) && ncol(sd_array) >= r) {
					sd_array[, r]
				} else {
					rep(NA_real_, dim(x)[1])
				}
				dfs[[idx]] <- data.frame(
					actor = actors,
					dimension = r,
					time = time_labels[t],
					value = x[, r, t],
					posterior_sd = sd_vals,
					type = type,
					stringsAsFactors = FALSE,
					row.names = NULL
				)
				idx <- idx + 1L
			}
		}
		do.call(rbind, dfs)
	} else {
		data.frame(
			actor = character(0), dimension = integer(0),
			time = character(0), value = numeric(0),
			posterior_sd = numeric(0),
			type = character(0), stringsAsFactors = FALSE
		)
	}
}
