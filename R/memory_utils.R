#' Memory optimization utilities for AME models
#' 
#' @description
#' Functions to manage memory usage in AME models, especially for large networks.
#' For a network of size n x n with R latent dimensions:
#' - YPM: n x n matrix (8n^2 bytes for doubles)
#' - EZ: n x n matrix (8n^2 bytes)
#' - UVPM: n x n matrix (8n^2 bytes)
#' - U: n x R matrix (8nR bytes)
#' - V: n x R matrix (8nR bytes)
#' 
#' For n=1000, each n x n matrix uses ~8MB. For n=10000, each uses ~800MB.
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau

#' Optimize AME model output for memory efficiency
#' 
#' @description
#' Optimizes AME model storage by optionally removing redundant components
#' and/or using sparse matrices based on user preferences.
#' 
#' @param fit Fitted AME model
#' @param use_sparse_matrices Logical; whether to use sparse matrix storage (default: FALSE)
#' 
#' @return Memory-optimized AME model object
#' 
#' @details
#' This function provides two optimization strategies:
#' 
#' 1. remove_redundant = TRUE: Removes redundant matrices (EZ, UVPM) that can be 
#'    reconstructed if needed, keeping only essential components (BETA, VC, YPM).
#'    Recommended for networks with > 100 nodes.
#' 
#' 2. use_sparse_matrices = TRUE: Converts matrices to sparse format. 
#'    Only recommended when your network is actually sparse (< 10\% non-zero entries).
#'    For dense networks, sparse storage will be slower and may use more memory.
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
compact_ame <- function(fit, use_sparse_matrices = FALSE) {
  
  # Only apply sparse matrices if requested
  if(!use_sparse_matrices) {
    return(fit)  # Return unchanged if no sparse matrices requested
  }
  
  # Start with the current fit (already has redundant matrices removed)
  compact_fit <- fit
  
  # Remove NULL entries
  compact_fit <- compact_fit[!sapply(compact_fit, is.null)]
  
  # Use sparse matrices if requested by user
  if(use_sparse_matrices && requireNamespace("Matrix", quietly = TRUE)) {
    if(!is.null(compact_fit$YPM)) {
      compact_fit$YPM <- Matrix::Matrix(fit$YPM, sparse = TRUE)
    }
    if(!is.null(compact_fit$APM)) {
      compact_fit$APM <- Matrix::Matrix(fit$APM, sparse = TRUE)
    }
    if(!is.null(compact_fit$BPM)) {
      compact_fit$BPM <- Matrix::Matrix(fit$BPM, sparse = TRUE)
    }
    cli::cli_inform("Using sparse matrix storage as requested")
  } else if(use_sparse_matrices) {
    cli::cli_warn("Matrix package not available; using dense storage")
  }
  
  # Add posterior samples if they exist (already thinned)
  if(!is.null(fit$U_samples)) compact_fit$U_samples <- fit$U_samples
  if(!is.null(fit$V_samples)) compact_fit$V_samples <- fit$V_samples
  if(!is.null(fit$a_samples)) compact_fit$a_samples <- fit$a_samples
  if(!is.null(fit$b_samples)) compact_fit$b_samples <- fit$b_samples
  
  # Preserve class
  class(compact_fit) <- class(fit)
  
  # Silent optimization - no need to report every time
  return(compact_fit)
}

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
ame_memory_usage <- function(fit, detailed = TRUE) {
  
  # Get all object names and sizes
  components <- names(fit)
  sizes <- sapply(components, function(x) {
    as.numeric(object.size(fit[[x]]))
  })
  
  # Create summary table
  memory_df <- data.frame(
    Component = components,
    Bytes = sizes,
    MB = round(sizes / 1024^2, 2),
    stringsAsFactors = FALSE
  )
  
  # Sort by size
  memory_df <- memory_df[order(memory_df$Bytes, decreasing = TRUE), ]
  
  if(detailed) {
    cli::cli_h3("AME Model Memory Usage")
    
    # Group by type
    param_components <- c("BETA", "VC")
    pred_components <- c("YPM", "EZ", "UVPM")
    factor_components <- c("U", "V", "G", "L")
    sample_components <- c("U_samples", "V_samples", "a_samples", "b_samples")
    
    total_mb <- sum(memory_df$MB)
    
    cli::cli_text("Total: {.val {round(total_mb, 2)}} MB")
    cli::cli_text("")
    
    # Show by category
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

#' Display memory usage information for AME models
#' 
#' @description
#' Shows estimated memory usage for networks of given size.
#' Memory optimization is now automatic, so this is informational only.
#' 
#' @param n_nodes Number of nodes in network
#' @param R Rank of multiplicative effects (default: 2)
#' 
#' @return Invisibly returns memory estimates
#' 
#' @details
#' Memory optimization is now user-controlled:
#' - Use remove_redundant=TRUE to remove EZ and UVPM matrices
#' - Use use_sparse_matrices=TRUE for sparse networks
#' - Both options can be combined for maximum memory savings
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
ame_memory_settings <- function(n_nodes, R = 2) {
  
  # Calculate expected memory usage
  matrix_size_mb <- (8 * n_nodes^2) / 1024^2
  factor_size_mb <- (8 * n_nodes * R * 2) / 1024^2
  
  cli::cli_h3("Memory Usage Estimates for n={.val {n_nodes}} Network")
  
  cli::cli_text("Estimated memory usage:")
  cli::cli_text("  • Each n×n matrix: {.val {round(matrix_size_mb, 1)}} MB")
  cli::cli_text("  • Factor matrices (U,V): {.val {round(factor_size_mb, 1)}} MB")
  
  # Show memory for different options
  full_size <- matrix_size_mb * 3 + factor_size_mb  # YPM, EZ, UVPM + factors
  optimized_size <- matrix_size_mb + factor_size_mb  # Only YPM + factors
  
  cli::cli_text("")
  cli::cli_text("Memory usage by configuration:")
  cli::cli_text("  • Default (all matrices): {.val {round(full_size, 1)}} MB")
  cli::cli_text("  • With remove_redundant=TRUE: {.val {round(optimized_size, 1)}} MB")
  cli::cli_text("    (Saves {.val {round((1 - optimized_size/full_size) * 100, 0)}}% memory)")
  
  cli::cli_text("")
  cli::cli_text("Recommendations:")
  if(n_nodes > 100) {
    cli::cli_alert_info("For n > 100, consider using remove_redundant=TRUE")
  }
  if(n_nodes > 1000) {
    cli::cli_alert_info("For n > 1000, consider use_sparse_matrices=TRUE if network is sparse")
  }
  
  invisible(list(
    matrix_mb = matrix_size_mb,
    factor_mb = factor_size_mb,
    total_mb = matrix_size_mb * 3 + factor_size_mb,
    optimized_mb = if(n_nodes > 100) matrix_size_mb + factor_size_mb else matrix_size_mb * 3 + factor_size_mb
  ))
}