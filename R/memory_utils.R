# memory utilities for ame fits; each n x n matrix takes 8 n^2 bytes

#' Optimize AME model output for memory efficiency
#'
#' @description
#' Optimizes AME model storage by converting matrices to sparse format.
#'
#' @param fit Fitted AME model
#' @param use_sparse_matrices Logical; whether to use sparse matrix storage (default: FALSE)
#'
#' @return Memory-optimized AME model object
#'
#' @details
#' When \code{use_sparse_matrices = TRUE}, the posterior-mean matrices
#' \code{YPM} and \code{EZ} are converted to sparse storage via the
#' \pkg{Matrix} package, but only when a matrix's nonzero density is below
#' 0.5 -- converting a dense posterior-mean matrix would \emph{grow} memory,
#' so denser matrices are left as ordinary dense matrices. For \code{lame}
#' fits, where \code{YPM}/\code{EZ} are per-period lists, the check and
#' conversion are applied per element. The additive-effect summaries
#' \code{APM}/\code{BPM} are named numeric vectors and are never converted.
#' Names and dimnames are preserved.
#'
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' fit_small <- compact_ame(fit)
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
compact_ame <- function(fit, use_sparse_matrices = FALSE) {

	# drop NULL entries; sparse posterior-mean matrices optionally become Matrix objects
	compact_fit <- fit
	compact_fit <- compact_fit[!sapply(compact_fit, is.null)]

	if (!use_sparse_matrices) {
		# without sparse conversion, return the NULL-stripped fit
		class(compact_fit) <- class(fit)
		return(compact_fit)
	}

	if (requireNamespace("Matrix", quietly = TRUE)) {
		# YPM and EZ (per-period lists for lame) become sparse only when their nonzero
		# density, counting NA cells, is below 0.5; APM and BPM stay numeric vectors
		sparsify_matrix <- function(m) {
			if (!is.matrix(m) || !is.numeric(m) || length(m) == 0L) return(m)
			density <- sum(m != 0 | is.na(m)) / length(m)
			if (isTRUE(density < 0.5)) Matrix::Matrix(m, sparse = TRUE) else m
		}
		sparsify <- function(x) {
			if (is.list(x) && !is.matrix(x)) lapply(x, sparsify_matrix) else sparsify_matrix(x)
		}
		for (comp in c("YPM", "EZ")) {
			if (!is.null(compact_fit[[comp]])) {
				compact_fit[[comp]] <- sparsify(compact_fit[[comp]])
			}
		}
		cli::cli_inform("Sparse storage applied to posterior-mean matrices with nonzero density below {.val {0.5}}; denser matrices kept dense.")
	} else {
		cli::cli_warn("Matrix package not available; using dense storage")
	}

	class(compact_fit) <- class(fit)
	return(compact_fit)
}

#' Calculate memory usage of AME model components
#'
#' @param fit Fitted AME model
#' @param detailed Logical; show detailed breakdown (default TRUE)
#'
#' @return Invisibly returns data.frame with memory usage
#'
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' usage <- ame_memory_usage(fit, detailed = FALSE)
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
ame_memory_usage <- function(fit = NULL, detailed = TRUE) {

	if (is.null(fit)) {
		cli::cli_abort(c(
			"{.arg fit} is required.",
			"i" = "Pass a fitted {.cls ame} / {.cls lame} object: {.code ame_memory_usage(fit)}."
		))
	}

	# size of each component
	components <- names(fit)
	sizes <- sapply(components, function(x) {
		as.numeric(object.size(fit[[x]]))
	})

	memory_df <- data.frame(
		Component = components,
		Bytes = sizes,
		MB = round(sizes / 1024^2, 2),
		stringsAsFactors = FALSE
	)

	memory_df <- memory_df[order(memory_df$Bytes, decreasing = TRUE), ]

	if (detailed) {
		cli::cli_h3("AME Model Memory Usage")

		# component groups
		param_components <- c("BETA", "VC")
		pred_components <- c("YPM", "EZ", "UVPM")
		factor_components <- c("U", "V", "G", "L")
		sample_components <- c("U_samples", "V_samples", "a_samples", "b_samples")

		total_mb <- sum(memory_df$MB)

		cli::cli_text("Total: {.val {round(total_mb, 2)}} MB")
		cli::cli_text("")

		for (category in list(
			list(name = "Predictions", components = pred_components),
			list(name = "Parameters", components = param_components),
			list(name = "Factors", components = factor_components),
			list(name = "Samples", components = sample_components)
		)) {
			cat_rows <- memory_df[memory_df$Component %in% category$components, ]
			if (nrow(cat_rows) > 0) {
				cat_total <- sum(cat_rows$MB)
				cli::cli_text("{.strong {category$name}}: {.val {round(cat_total, 2)}} MB")
				for (i in 1:nrow(cat_rows)) {
					cli::cli_text("  {cat_rows$Component[i]}: {.val {cat_rows$MB[i]}} MB")
				}
			}
		}
	}

	invisible(memory_df)
}

#' Display memory usage information for AME models
#'
#' @description
#' Shows estimated memory usage for networks of given size.
#' Memory optimization is automatic, so this is informational only.
#'
#' @param n_nodes Number of nodes in network
#' @param R Rank of multiplicative effects (default: 2)
#'
#' @return Invisibly returns memory estimates
#'
#' @details
#' Memory levers available in the package:
#' - Run \code{compact_ame()} on a fitted model to drop empty slots and,
#'   for genuinely sparse posterior means, use sparse storage via
#'   \code{use_sparse_matrices = TRUE}
#' - Increase \code{odens} in \code{ame()}/\code{lame()} to store fewer
#'   posterior draws
#' - Pass \code{posterior_opts = list(thin_UV = ..., thin_ab = ...)} to thin
#'   the stored latent-factor and additive-effect draws
#'
#' @examples
#' ame_memory_settings(n_nodes = 12, R = 1)
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
ame_memory_settings <- function(n_nodes, R = 2) {

	# expected memory usage
	matrix_size_mb <- (8 * n_nodes^2) / 1024^2
	factor_size_mb <- (8 * n_nodes * R * 2) / 1024^2

	cli::cli_h3("Memory Usage Estimates for n={.val {n_nodes}} Network")

	cli::cli_text("Estimated memory usage:")
	cli::cli_text("  - Each n x n matrix: {.val {round(matrix_size_mb, 1)}} MB")
	cli::cli_text("  - Factor matrices (U,V): {.val {round(factor_size_mb, 1)}} MB")

	# totals for each storage configuration; networks of more than 100 nodes keep YPM only
	full_size <- matrix_size_mb * 3 + factor_size_mb # YPM, EZ, UVPM and factors
	ypm_only <- n_nodes > 100
	optimized_size <- if (ypm_only) matrix_size_mb + factor_size_mb else full_size

	cli::cli_text("")
	cli::cli_text("Memory usage by configuration:")
	cli::cli_text("  - All posterior-mean matrices (YPM, EZ, UVPM): {.val {round(full_size, 1)}} MB")
	if (ypm_only) {
		cli::cli_text("  - Optimized (YPM only, for n > 100): {.val {round(optimized_size, 1)}} MB")
	} else {
		cli::cli_text("  - Optimized (all matrices kept for n <= 100): {.val {round(optimized_size, 1)}} MB")
	}

	cli::cli_text("")
	cli::cli_text("Recommendations:")
	cli::cli_text("  - Run {.fn compact_ame} on the fitted model to drop empty slots.")
	cli::cli_text("  - Increase {.arg odens} to store fewer posterior draws.")
	cli::cli_text("  - Pass {.code posterior_opts = list(thin_UV = ..., thin_ab = ...)} to thin stored latent-factor and additive-effect draws.")
	if (n_nodes > 1000) {
		cli::cli_alert_info("For n > 1000, {.code compact_ame(fit, use_sparse_matrices = TRUE)} helps if the posterior means are genuinely sparse.")
	}

	invisible(list(
		matrix_mb = matrix_size_mb,
		factor_mb = factor_size_mb,
		total_mb = full_size,
		optimized_mb = optimized_size
	))
}
