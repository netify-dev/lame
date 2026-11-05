#' Extract latent positions as a tidy data frame
#'
#' @description
#' Extracts multiplicative latent factor positions (U and V) from a fitted
#' \code{ame}, \code{lame} or \code{ame_als} model and returns them as a
#' tidy data frame suitable for plotting and analysis, with posterior standard
#' deviations when posterior draws of the positions are available.
#'
#' @param object A fitted \code{ame}, \code{lame} or \code{ame_als} model
#'   object with R > 0.
#' @param align Logical. For dynamic \code{ame_als} fits, apply Procrustes
#'   alignment across the time periods of the fit. It has no effect on
#'   \code{ame} fits (a single period) or on \code{lame} fits, whose dynamic
#'   positions are aligned when the model is fitted (see Details).
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The latent product of each period is identified, the positions only up to a
#' rotation or sign change that leaves the product unchanged. The positions of
#' a \code{lame} fit with \code{dynamic_uv = TRUE} therefore factor the
#' posterior mean product of every period at the fitted rank: \eqn{U_t V_t'}
#' with the scale split evenly between \eqn{U_t} and \eqn{V_t};
#' \eqn{U_t L U_t'} for symmetric fits, with \code{fit$L} fixed across periods;
#' and \eqn{U_t G V_t'} for bipartite fits, with a diagonal \code{fit$G} fixed
#' across periods. Each period is matched to the one before by moves that keep
#' its product. The posterior standard deviations come from the saved draws
#' (\code{posterior_opts = list(save_UV_draws = TRUE)}), each draw factored the
#' same way and matched onto the reported positions.
#'
#' Alignment only ever happens across the periods of one fit. Nothing here aligns two
#' different fits (two chains, or the same model with and without a
#' covariate) to a common orientation, so their latent positions are not
#' directly comparable. To compare fits that way, stack their \code{U}
#' matrices into an \code{[n, R, 2]} array (same actors in the same order,
#' same \code{R}) and call \code{procrustes_align(U = that_array)} on it.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{actor}{Character. Actor name (from rownames of U or V).}
#'   \item{dimension}{Integer. Latent dimension index (1 to R).}
#'   \item{time}{Character. Time period label. Dynamic fits use the time
#'     labels from the input; static (cross-sectional) fits return
#'     \code{"1"} for every row so downstream filtering by \code{time}
#'     behaves the same in both cases.}
#'   \item{value}{Numeric. The posterior mean latent position.}
#'   \item{posterior_sd}{Numeric. Posterior standard deviation of the latent
#'     position, or \code{NA} if posterior draws are not available. To enable,
#'     fit the model with \code{posterior_opts = posterior_options(save_UV = TRUE)}
#'     (static positions) or \code{posterior_opts = list(save_UV_draws = TRUE)}
#'     (dynamic positions).}
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
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            burn = 5, nscan = 5, odens = 1, verbose = FALSE, chains = 1)
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
	.lame_require_draws(object, "latent_positions")
	.extract_latent_positions(object, align = align)
}

#' @rdname latent_positions
#' @method latent_positions lame
#' @export
latent_positions.lame <- function(object, align = TRUE, ...) {
	.lame_require_draws(object, "latent_positions")
	if (isTRUE(object$dynamic_uv)) {
		# dynamic positions are factored and aligned when the model is fitted
		if (!identical(object$uv_positions, "posterior_mean_product")) {
			object <- .lame_dynamic_positions(object)
		}
		align <- FALSE
	}
	.extract_latent_positions(object, align = align)
}

#' @rdname latent_positions
#' @method latent_positions ame_als
#' @export
latent_positions.ame_als <- function(object, align = FALSE, ...) {
	# ALS positions are static with no draws, so posterior_sd is NA; a withheld
	# product (.ae_withhold_product) gives NA positions
	if (.ae_product_is_withheld(object)) {
		pw <- object$eb$product_withheld
		cli::cli_inform(
			c(
				"i" = "The latent product of this fit is withheld: it is not identified by the data at rank {pw$rank} ({signif(pw$per_loading, 2)} noise variances of signal per loading, below {pw$threshold}); the positions are NA."
			)
		)
	}
	.extract_latent_positions(object, align = align)
}

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
	V <- object[["V"]]
	G <- object[["G"]]

	# a fit without a multiplicative term carries an empty (n x 0) U
	if (is.null(U) || length(U) == 0L || NCOL(U) == 0L) {
		cli::cli_inform("This fit has no latent positions ({.code R = 0}); returning an empty table.")
		return(empty_df)
	}

	is_symmetric <- isTRUE(object$symmetric)

	# Procrustes alignment if requested
	if (align && length(dim(U)) == 3L) {
		aligned <- procrustes_align(object = object)
		U <- aligned$U
		V <- aligned[["V"]]
		G <- aligned$G
	}

	# posterior sds: static fits keep [n, R, S] draws in U_samples / V_samples;
	# dynamic fits carry the sds of the reported coordinates in U_sd / V_sd
	U_sd <- .compute_posterior_sd(object$U_samples)
	V_sd <- .compute_posterior_sd(object$V_samples)
	if (is.null(U_sd) && !is.null(object$U_sd) &&
			identical(dim(object$U_sd), dim(U))) {
		U_sd <- object$U_sd
	}
	if (is.null(V_sd) && !is.null(object$V_sd) && !is.null(V) &&
			identical(dim(object$V_sd), dim(V))) {
		V_sd <- object$V_sd
	}

	u_df <- .array_to_df(U, type = "U", sd_array = U_sd)

	if (is_symmetric || is.null(V)) {
		out <- u_df
	} else {
		v_df <- .array_to_df(V, type = "V", sd_array = V_sd)
		out <- rbind(u_df, v_df)
	}

	# explain why posterior_sd is NA
	if (nrow(out) > 0L && all(is.na(out$posterior_sd)) && inherits(object, "ame_als")) {
		cli::cli_inform(
			c(
				"i" = "{.code posterior_sd} is {.val NA} because an ALS fit holds point estimates of the positions and no posterior draws.",
				"i" = "For posterior SDs, fit by MCMC with {.fn ame} or {.fn lame} and {.code posterior_opts = posterior_options(save_UV = TRUE)}."
			),
			.frequency = "once",
			.frequency_id = "latent_positions_no_sd_als"
		)
	} else if (nrow(out) > 0L && all(is.na(out$posterior_sd))) {
		cli::cli_inform(
			c(
				"i" = "{.code posterior_sd} is {.val NA} because U/V draws were not saved.",
				"i" = if (length(dim(U)) == 3L) {
					"To get posterior SDs, refit with {.code posterior_opts = list(save_UV_draws = TRUE)}."
				} else {
					"To get posterior SDs, refit with {.code posterior_opts = posterior_options(save_UV = TRUE)}."
				}
			),
			.frequency = "once",
			.frequency_id = "latent_positions_no_sd"
		)
	}

	out
}

# element-wise posterior sd from a 3-d draws array [n, R, n_samples]
.compute_posterior_sd <- function(samples) {
	if (is.null(samples)) return(NULL)
	if (length(dim(samples)) != 3L) return(NULL)
	apply(samples, c(1, 2), sd, na.rm = TRUE)
}

# 2-d matrix or 3-d array to a tidy data frame; sd_array an optional array of
# posterior sds matching x
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
				# static fits get time = "1" so they filter like dynamic fits
				time = "1",
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
				# a 3-d sd_array holds per-period sds; a 2-d one is constant over time
				sd_vals <- if (!is.null(sd_array) && length(dim(sd_array)) == 3L &&
						dim(sd_array)[2] >= r && dim(sd_array)[3] >= t) {
					sd_array[, r, t]
				} else if (!is.null(sd_array) && length(dim(sd_array)) == 2L &&
						ncol(sd_array) >= r) {
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
			actor = character(0),
			dimension = integer(0),
			time = character(0),
			value = numeric(0),
			posterior_sd = numeric(0),
			type = character(0),
			stringsAsFactors = FALSE
		)
	}
}
