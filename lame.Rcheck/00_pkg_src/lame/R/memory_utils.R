# memory optimization utilities for ame models
#
# manages memory usage for large networks
# each n x n matrix uses 8n^2 bytes (e.g., ~8mb for n=1000)

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
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
compact_ame <- function(fit, use_sparse_matrices = FALSE) {

	# drop null entries unconditionally. when use_sparse_matrices = true we
	# additionally cast dense posterior-mean matrices to matrix::matrix.
	compact_fit <- fit
	compact_fit <- compact_fit[!sapply(compact_fit, is.null)]

	if(!use_sparse_matrices) {
		# preserve class and return the null-stripped fit. without sparse
		# conversion the size reduction is modest and depends on how many
		# slots were left empty by the fit assembly.
		class(compact_fit) <- class(fit)
		return(compact_fit)
	}

	# sparse matrices
	if(use_sparse_matrices && requireNamespace("Matrix", quietly = TRUE)) {
		# only the posterior-mean matrices ypm/ez are candidates; apm/bpm are
		# named numeric vectors and must stay that way. lame fits store ypm/ez
		# as per-period lists, so convert element-wise. a matrix is cast to
		# sparse only when its nonzero density (counting na cells, which sparse
		# storage must also hold) is below 0.5 -- converting a dense
		# posterior-mean matrix grows memory instead of shrinking it.
		sparsify_matrix <- function(m) {
			if(!is.matrix(m) || !is.numeric(m) || length(m) == 0L) return(m)
			density <- sum(m != 0 | is.na(m)) / length(m)
			if(isTRUE(density < 0.5)) Matrix::Matrix(m, sparse = TRUE) else m
		}
		sparsify <- function(x) {
			if(is.list(x) && !is.matrix(x)) lapply(x, sparsify_matrix) else sparsify_matrix(x)
		}
		for(comp in c("YPM", "EZ")) {
			if(!is.null(compact_fit[[comp]])) {
				compact_fit[[comp]] <- sparsify(compact_fit[[comp]])
			}
		}
		cli::cli_inform("Sparse storage applied to posterior-mean matrices with nonzero density below {.val {0.5}}; denser matrices kept dense.")
	} else if(use_sparse_matrices) {
		cli::cli_warn("Matrix package not available; using dense storage")
	}
	
	# posterior samples
	if(!is.null(fit$U_samples)) compact_fit$U_samples <- fit$U_samples
	if(!is.null(fit$V_samples)) compact_fit$V_samples <- fit$V_samples
	if(!is.null(fit$a_samples)) compact_fit$a_samples <- fit$a_samples
	if(!is.null(fit$b_samples)) compact_fit$b_samples <- fit$b_samples
	
	# preserve class
	class(compact_fit) <- class(fit)
	return(compact_fit)
}

####
# ame_memory_usage
####

#' Calculate memory usage of AME model components
#' 
#' @param fit Fitted AME model
#' @param detailed Logical; show detailed breakdown (default TRUE)
#' 
#' @return Invisibly returns data.frame with memory usage
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
ame_memory_usage <- function(fit = NULL, detailed = TRUE) {

	if (is.null(fit)) {
		cli::cli_abort(c(
			"{.arg fit} is required.",
			"i" = "Pass a fitted {.cls ame} / {.cls lame} object: {.code ame_memory_usage(fit)}."))
	}

	# object names and sizes
	components <- names(fit)
	sizes <- sapply(components, function(x) {
		as.numeric(object.size(fit[[x]]))
	})
	
	# summary table
	memory_df <- data.frame(
		Component = components,
		Bytes = sizes,
		MB = round(sizes / 1024^2, 2),
		stringsAsFactors = FALSE
	)
	
	# sort by size
	memory_df <- memory_df[order(memory_df$Bytes, decreasing = TRUE), ]
	
	if(detailed) {
		cli::cli_h3("AME Model Memory Usage")
		
		# group by type
		param_components <- c("BETA", "VC")
		pred_components <- c("YPM", "EZ", "UVPM")
		factor_components <- c("U", "V", "G", "L")
		sample_components <- c("U_samples", "V_samples", "a_samples", "b_samples")
		
		total_mb <- sum(memory_df$MB)
		
		cli::cli_text("Total: {.val {round(total_mb, 2)}} MB")
		cli::cli_text("")
		
		# by category
		for(category in list(
			list(name = "Predictions", components = pred_components),
			list(name = "Parameters", components = param_components),
			list(name = "Factors", components = factor_components),
			list(name = "Samples", components = sample_components)
		)) {
			cat_rows <- memory_df[memory_df$Component %in% category$components, ]
			if(nrow(cat_rows) > 0) {
				cat_total <- sum(cat_rows$MB)
				cli::cli_text("{.strong {category$name}}: {.val {round(cat_total, 2)}} MB")
				for(i in 1:nrow(cat_rows)) {
					cli::cli_text("  {cat_rows$Component[i]}: {.val {cat_rows$MB[i]}} MB")
				}
			}
		}
	}
	
	invisible(memory_df)
}

####
# ame_memory_settings
####

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
	
	# memory for different options
	full_size <- matrix_size_mb * 3 + factor_size_mb  # ypm, ez, uvpm + factors
	optimized_size <- matrix_size_mb + factor_size_mb  # only ypm + factors
	
	cli::cli_text("")
	cli::cli_text("Memory usage by configuration:")
	cli::cli_text("  - All posterior-mean matrices (YPM, EZ, UVPM): {.val {round(full_size, 1)}} MB")
	cli::cli_text("  - YPM only: {.val {round(optimized_size, 1)}} MB")

	cli::cli_text("")
	cli::cli_text("Recommendations:")
	cli::cli_text("  - Run {.fn compact_ame} on the fitted model to drop empty slots.")
	cli::cli_text("  - Increase {.arg odens} to store fewer posterior draws.")
	cli::cli_text("  - Pass {.code posterior_opts = list(thin_UV = ..., thin_ab = ...)} to thin stored latent-factor and additive-effect draws.")
	if(n_nodes > 1000) {
		cli::cli_alert_info("For n > 1000, {.code compact_ame(fit, use_sparse_matrices = TRUE)} helps if the posterior means are genuinely sparse.")
	}
	
	invisible(list(
		matrix_mb = matrix_size_mb,
		factor_mb = factor_size_mb,
		total_mb = matrix_size_mb * 3 + factor_size_mb,
		optimized_mb = if(n_nodes > 100) matrix_size_mb + factor_size_mb else matrix_size_mb * 3 + factor_size_mb
	))
}
