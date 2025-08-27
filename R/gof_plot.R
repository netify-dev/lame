#' Visualize goodness-of-fit statistics for AME and LAME models
#' 
#' @description
#' Creates comprehensive diagnostic plots comparing observed network statistics to their 
#' posterior predictive distributions. This function is essential for assessing whether 
#' the fitted AME/LAME model adequately captures the structural properties and dependencies
#' present in the observed network data.
#' 
#' @details
#' \strong{Overview:}
#' 
#' Goodness-of-fit (GOF) assessment is crucial for network models because standard
#' residual diagnostics are often inadequate for capturing network dependencies.
#' This function implements posterior predictive checking by:
#' 
#' 1. Computing key network statistics from the observed data
#' 2. Generating multiple networks from the model's posterior predictive distribution
#' 3. Computing the same statistics on simulated networks
#' 4. Visualizing the comparison to identify model inadequacies
#' 
#' \strong{Network Statistics Evaluated:}
#' 
#' \emph{For unipartite (square) networks:}
#' \describe{
#'   \item{\code{sd.row} (Out-degree heterogeneity)}{
#'     Standard deviation of row means. High values indicate substantial variation
#'     in how active nodes are as senders/initiators. If the model underestimates
#'     this, it may be missing important sender effects or covariates.
#'   }
#'   \item{\code{sd.col} (In-degree heterogeneity)}{
#'     Standard deviation of column means. High values indicate substantial variation
#'     in node popularity as receivers. Underestimation suggests missing receiver
#'     effects or popularity-related covariates.
#'   }
#'   \item{\code{dyad.dep} (Reciprocity/Mutuality)}{
#'     Correlation between Y[i,j] and Y[j,i]. Positive values indicate reciprocity
#'     (mutual ties are more likely). The AME model captures this through the
#'     dyadic correlation parameter rho. Poor fit here suggests the need to adjust
#'     the dcor parameter.
#'   }
#'   \item{\code{triad.dep} (Transitivity/Clustering)}{
#'     Measures tendency for triadic closure (friend of a friend is a friend).
#'     Calculated as correlation between Y[i,j] and sum(Y[i,k]*Y[k,j])/sqrt(n-2).
#'     AME captures this through multiplicative effects (U,V). Poor fit suggests
#'     increasing the latent dimension R.
#'   }
#' }
#' 
#' \emph{For bipartite (rectangular) networks:}
#' \describe{
#'   \item{\code{sd.row} (Type A activity variation)}{
#'     Standard deviation of row means for Type A nodes. Indicates heterogeneity
#'     in how actively Type A nodes connect to Type B nodes.
#'   }
#'   \item{\code{sd.col} (Type B popularity variation)}{
#'     Standard deviation of column means for Type B nodes. Indicates heterogeneity
#'     in how popular Type B nodes are with Type A nodes.
#'   }
#'   \item{\code{four.cycles} (Bipartite clustering)}{
#'     Count of 4-cycles (rectangular paths A1→B1→A2→B2→A1). High values indicate
#'     that pairs of Type A nodes tend to connect to the same Type B nodes.
#'     Captured through bipartite multiplicative effects with appropriate R_row, R_col.
#'   }
#' }
#' 
#' \strong{Interpretation Guide:}
#' 
#' \emph{Histogram plots (static models):}
#' - Red vertical line: observed statistic value
#' - Blue histogram: distribution from posterior predictive simulations
#' - Good fit: red line falls within the bulk of the histogram
#' - Poor fit: red line in the tail or outside the distribution
#' 
#' \emph{Common model inadequacies and solutions:}
#' \describe{
#'   \item{Observed sd.row/sd.col too high}{
#'     Model underestimates degree heterogeneity.
#'     Solutions: Add row/column covariates (Xrow, Xcol), enable random effects
#'     (rvar=TRUE, cvar=TRUE), or increase their variance.
#'   }
#'   \item{Observed dyad.dep outside distribution}{
#'     Reciprocity not captured well.
#'     Solutions: For positive reciprocity, ensure dcor=TRUE. For negative,
#'     consider transformation or different family.
#'   }
#'   \item{Observed triad.dep too high}{
#'     Clustering/transitivity underestimated.
#'     Solutions: Increase latent dimension R, add network covariates that
#'     capture homophily, or consider including community structure covariates.
#'   }
#'   \item{Multiple statistics showing poor fit}{
#'     Fundamental model misspecification.
#'     Solutions: Change family, add missing covariates, or consider
#'     different model class.
#'   }
#' }
#' 
#' \strong{Mathematical Details:}
#' 
#' The posterior predictive p-value for statistic s is:
#' \deqn{p = P(s(Y_{rep}) \geq s(Y_{obs}) | Y_{obs})}
#' 
#' where Y_rep is drawn from the posterior predictive distribution.
#' Values near 0 or 1 indicate poor fit. The function visualizes the full
#' distribution rather than just p-values for richer diagnostics.
#' 
#' \strong{Customization:}
#' 
#' The function accepts custom GOF statistics through the model's custom_gof
#' argument. These are automatically included in the plot. Custom statistics
#' should capture network features important for your specific application.
#' 
#' \strong{Computational Notes:}
#' 
#' GOF computation involves simulating multiple networks, which can be
#' computationally intensive. The model uses the networks simulated during
#' MCMC if gof=TRUE was specified. For post-hoc GOF, use the gof() function
#' which generates new simulations.
#' 
#' Good model fit is indicated when:
#' - Observed statistics fall within the central 95% of posterior predictive distributions
#' - No systematic patterns of over/underestimation across statistics
#' - Custom statistics (if provided) also show adequate fit
#' - For longitudinal models, observed trajectories track within credible bands
#' 
#' @param fit An object of class "ame" or "lame" containing GOF statistics
#' @param type Character string: "auto" (default), "static", or "longitudinal".
#'        If "auto", determined by model class.
#' @param statistics Character vector specifying which statistics to plot.
#'        Default is all four: c("sd.row", "sd.col", "dyad.dep", "triad.dep")
#' @param credible.level Numeric between 0 and 1; credible interval level for
#'        longitudinal plots (default 0.95)
#' @param ncol Number of columns for faceted plot layout (default 2)
#' @param point.size Size of points in longitudinal plots (default 2)
#' @param line.size Width of lines in plots (default 1)
#' @param title Optional title for the plot
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \dontrun{
#' # Fit an AME model
#' fit_ame <- ame(Y, X, gof = TRUE)
#' 
#' # Basic GOF plot
#' gof_plot(fit_ame)
#' 
#' # Plot only degree-related statistics
#' gof_plot(fit_ame, statistics = c("sd.row", "sd.col"))
#' 
#' # Fit a LAME model
#' fit_lame <- lame(Y_list, X_list)
#' 
#' # Longitudinal GOF plot with 90% credible intervals
#' gof_plot(fit_lame, credible.level = 0.90)
#' }
#' @export
#' @import ggplot2
#' @import reshape2
gof_plot <- function(
  fit, 
  type = c("auto", "static", "longitudinal"),
  statistics = c("sd.row", "sd.col", "dyad.dep", "triad.dep"),
  credible.level = 0.95,
  ncol = 2,
  point.size = 2,
  line.size = 1,
  title = NULL) {
  
  # Check input
  if (!inherits(fit, c("ame", "lame"))) {
    stop("fit must be an object of class 'ame' or 'lame'")
  }
  
  # Check for GOF statistics
  if (is.null(fit$GOF)) {
    stop("No GOF statistics found. Re-run model with gof = TRUE")
  }
  
  # Determine plot type
  type <- match.arg(type)
  if (type == "auto") {
    type <- ifelse(inherits(fit, "lame"), "longitudinal", "static")
  }
  
  # Detect if network is bipartite based on GOF columns
  is_bipartite <- FALSE
  if(!is.null(colnames(fit$GOF))) {
    is_bipartite <- "four.cycles" %in% colnames(fit$GOF)
  } else if(!is.null(fit$mode)) {
    is_bipartite <- fit$mode == "bipartite"
  }
  
  # Validate statistics argument based on network type
  if(is_bipartite) {
    # For bipartite networks
    stat.names <- c("sd.row" = "sd.rowmean", 
                    "sd.col" = "sd.colmean",
                    "four.cycles" = "four.cycles")
    
    # Filter requested statistics to only those applicable to bipartite
    valid_stats <- intersect(statistics, c("sd.row", "sd.col", "four.cycles"))
    if(length(valid_stats) == 0) {
      # If no valid stats specified, use all bipartite stats
      statistics <- c("sd.row", "sd.col", "four.cycles")
    } else {
      statistics <- valid_stats
    }
  } else {
    # For unipartite networks
    stat.names <- c("sd.row" = "sd.rowmean", 
                    "sd.col" = "sd.colmean",
                    "dyad.dep" = "dyad.dep", 
                    "triad.dep" = "cycle.dep",
                    "trans.dep" = "trans.dep")
    
    statistics <- match.arg(statistics, several.ok = TRUE,
                           choices = names(stat.names))
  }
  
  # Add any custom GOF columns to available statistics
  if (!is.null(colnames(fit$GOF))) {
    custom_cols <- setdiff(colnames(fit$GOF), unname(stat.names))
    if (length(custom_cols) > 0) {
      # Add custom columns to stat.names mapping (identity mapping)
      for (col in custom_cols) {
        stat.names[col] <- col
      }
      # If no specific statistics requested, include custom ones
      if (missing(statistics)) {
        statistics <- c(statistics, custom_cols)
      }
    }
  }
  
  # Create appropriate plot
  if (type == "static") {
    p <- gof_plot_static(fit, statistics, stat.names, ncol, line.size, title)
  } else {
    p <- gof_plot_longitudinal(fit, statistics, stat.names, credible.level,
                              ncol, point.size, line.size, title)
  }
  
  return(p)
}

# Helper function for static GOF plots
gof_plot_static <- function(fit, statistics, stat.names, ncol, line.size, title) {
  
  # Extract GOF data
  gof_data <- fit$GOF
  
  # Handle different GOF formats
  if (is.list(gof_data) && !is.data.frame(gof_data)) {
    # GOF is a named list (converted format)
    # Convert to matrix format [iterations, statistics]
    stat_names_full <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    n_iter <- nrow(gof_data[[1]])
    gof_matrix <- matrix(NA, n_iter, length(stat_names_full))
    colnames(gof_matrix) <- stat_names_full
    
    for (i in seq_along(stat_names_full)) {
      if (stat_names_full[i] %in% names(gof_data)) {
        # Average across time if there are multiple columns
        if (ncol(gof_data[[stat_names_full[i]]]) > 1) {
          gof_matrix[, i] <- rowMeans(gof_data[[stat_names_full[i]]], na.rm = TRUE)
        } else {
          gof_matrix[, i] <- gof_data[[stat_names_full[i]]][, 1]
        }
      }
    }
    gof_data <- gof_matrix
  } else if (length(dim(gof_data)) == 3) {
    # GOF is [statistics, time, iterations]
    # Average across time to get [statistics, iterations]
    gof_data_avg <- apply(gof_data, c(1, 3), mean, na.rm = TRUE)
    # Transpose to get [iterations, statistics]
    gof_data <- t(gof_data_avg)
  }
  
  obs_vals <- gof_data[1, ]
  pred_vals <- gof_data[-1, , drop = FALSE]
  
  # Prepare data for plotting
  plot_list <- list()
  
  for (stat in statistics) {
    # Check if this is a custom statistic or predefined
    if (stat %in% names(stat.names)) {
      stat_col <- stat.names[stat]
    } else {
      # Custom statistic - use as-is if it exists in columns
      stat_col <- stat
    }
    
    # Skip if column doesn't exist or no data
    if (!stat_col %in% colnames(gof_data)) {
      next
    }
    
    # Check if we have predictive samples
    if (nrow(pred_vals) == 0) {
      warning("No posterior predictive samples available for GOF plots")
      next
    }
    
    # Create data frame for this statistic
    df <- data.frame(
      value = as.vector(pred_vals[, stat_col]),
      statistic = stat,
      stringsAsFactors = FALSE
    )
    
    # Create individual plot
    p_stat <- ggplot(df, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), 
                    bins = 30, fill = "lightblue", 
                    color = "white", alpha = 0.7) +
      geom_vline(xintercept = obs_vals[stat_col], 
                color = "red", linewidth = line.size) +
      labs(
        x = "",
        y = "Density",
        title = stat
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5)
      )
    
    plot_list[[stat]] <- p_stat
  }
  
  # Combine plots
  if (length(plot_list) == 1) {
    p <- plot_list[[1]]
  } else {
    # Use patchwork or gridExtra to combine plots
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- patchwork::wrap_plots(plot_list, ncol = ncol)
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      p <- gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)
    } else {
      warning("Install 'patchwork' or 'gridExtra' for better plot layout")
      p <- plot_list[[1]]  # Return first plot only
    }
  }
  
  # Add overall title
  if (!is.null(title)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- p + patchwork::plot_annotation(title = title)
    }
  } else {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- p + patchwork::plot_annotation(
        title = "Goodness-of-Fit: Observed (red) vs Posterior Predictive"
      )
    }
  }
  
  return(p)
}

# Helper function for longitudinal GOF plots
gof_plot_longitudinal <- function(
  fit, statistics, stat.names, credible.level,
  ncol, point.size, line.size, title) {
  
  # Extract GOF data (should be time-indexed for LAME)
  if (!is.null(fit$GOF_T)) {
    gof_data <- fit$GOF_T
  } else if (!is.null(fit$GOF) && is.list(fit$GOF) && !is.data.frame(fit$GOF)) {
    # Handle GOF as named list format (new format from get_fit_object)
    stat_names_full <- c("sd.rowmean", "sd.colmean", "dyad.dep", "cycle.dep", "trans.dep")
    
    # Check if we have time-indexed data (multiple columns)
    if (ncol(fit$GOF[[1]]) > 1) {
      # Convert to list of matrices for each time point
      gof_data <- list()
      n_time <- ncol(fit$GOF[[1]])
      
      for (t in 1:n_time) {
        # Create matrix for this time point [iterations, statistics]
        gof_t <- matrix(NA, nrow(fit$GOF[[1]]), length(stat_names_full))
        colnames(gof_t) <- stat_names_full
        
        for (i in seq_along(stat_names_full)) {
          if (stat_names_full[i] %in% names(fit$GOF)) {
            gof_t[, i] <- fit$GOF[[stat_names_full[i]]][, t]
          }
        }
        gof_data[[t]] <- gof_t
      }
    } else {
      # Single time point, fall back to static
      warning("No time-indexed GOF found, using static GOF")
      return(gof_plot_static(fit, statistics, stat.names, ncol, line.size, title))
    }
  } else if (!is.null(fit$GOF) && length(dim(fit$GOF)) == 3) {
    # Handle LAME 3D GOF array [statistics, time, iterations]
    # Convert to list of matrices for each time point
    gof_data <- list()
    n_time <- dim(fit$GOF)[2]
    stat_names <- dimnames(fit$GOF)[[1]]
    for (t in 1:n_time) {
      # Extract and transpose to get [iterations, statistics]
      gof_t <- t(fit$GOF[, t, ])
      colnames(gof_t) <- stat_names
      gof_data[[t]] <- gof_t
    }
  } else {
    # Fall back to regular GOF if no time-indexed version
    warning("No time-indexed GOF found, using static GOF")
    return(gof_plot_static(fit, statistics, stat.names, ncol, line.size, title))
  }
  
  # Calculate quantiles
  alpha <- (1 - credible.level) / 2
  
  # Prepare data for plotting
  n_time <- length(gof_data)
  plot_data <- data.frame()
  
  for (t in 1:n_time) {
    gof_t <- gof_data[[t]]
    obs_t <- gof_t[1, ]
    pred_t <- gof_t[-1, , drop = FALSE]
    
    for (stat in statistics) {
      stat_col <- stat.names[stat]
      
      # Calculate statistics
      plot_data <- rbind(plot_data, data.frame(
        time = t,
        statistic = stat,
        observed = obs_t[stat_col],
        median = median(pred_t[, stat_col]),
        lower = quantile(pred_t[, stat_col], alpha),
        upper = quantile(pred_t[, stat_col], 1 - alpha),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = time)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
               alpha = 0.3, fill = "blue") +
    geom_line(aes(y = median), color = "blue", 
             linewidth = line.size, linetype = "dashed") +
    geom_line(aes(y = observed), color = "red", 
             linewidth = line.size) +
    geom_point(aes(y = observed), color = "red", 
              size = point.size) +
    facet_wrap(~ statistic, scales = "free_y", ncol = ncol) +
    labs(
      x = "Time Period",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      panel.spacing = unit(1, "lines")
    )
  
  # Add title
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  } else {
    p <- p + ggtitle(sprintf(
      "Longitudinal GOF: Observed (red) vs %d%% Credible Intervals",
      round(credible.level * 100)
    ))
  }
  
  return(p)
}