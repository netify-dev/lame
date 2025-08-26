#' MCMC trace plots and density plots for AME/LAME model parameters
#' 
#' Creates diagnostic plots for Markov Chain Monte Carlo (MCMC) samples from 
#' AME or LAME models. Displays trace plots to assess convergence and mixing,
#' alongside density plots to visualize posterior distributions.
#' 
#' @details
#' This function produces two types of diagnostic plots:
#' \describe{
#'   \item{Trace plots}{Show the evolution of parameter values across MCMC 
#'         iterations. Good mixing is indicated by rapid exploration of the 
#'         parameter space with no trends or stuck periods.}
#'   \item{Density plots}{Show the posterior distribution of parameters.
#'         Multiple modes may indicate identification issues or convergence problems.}
#' }
#' 
#' The plots help diagnose:
#' - Convergence: Has the chain reached the stationary distribution?
#' - Mixing: Is the chain exploring the parameter space efficiently?
#' - Autocorrelation: Are successive samples highly correlated?
#' 
#' Parameters displayed include:
#' - Regression coefficients (beta)
#' - Variance components (va, vb, cab, rho, ve)
#' 
#' @param fit An object of class "ame" or "lame" containing MCMC samples
#' @param params Character vector specifying which parameters to plot:
#'        "beta" for regression coefficients, "variance" for variance components,
#'        or "all" (default) for both
#' @param include Character vector of specific parameter names to include
#' @param exclude Character vector of specific parameter names to exclude
#' @param ncol Number of columns for plot layout (default 3)
#' @param nrow Number of rows for plot layout (default NULL, determined automatically)
#' @param burn.in Number of initial iterations to exclude as burn-in when 
#'        calculating statistics (default 0, assumes burn-in already removed)
#' @param thin Thinning interval for display (default 1, no thinning)
#' @param title Optional title for the plot
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \dontrun{
#' # Fit an AME model
#' fit <- ame(Y, X, nscan = 10000, burn = 1000)
#' 
#' # Basic trace plots for all parameters
#' trace_plot(fit)
#' 
#' # Only regression coefficients
#' trace_plot(fit, params = "beta")
#' 
#' # Only variance components
#' trace_plot(fit, params = "variance")
#' 
#' # Exclude intercept from plot
#' trace_plot(fit, exclude = "intercept")
#' 
#' # Thin the display for clearer visualization
#' trace_plot(fit, thin = 10)
#' }
#' @export
#' @import ggplot2
#' @import reshape2
trace_plot <- function(
  fit,
  params = c("all", "beta", "variance"),
  include = NULL,
  exclude = NULL,
  ncol = 3,
  nrow = NULL,
  burn.in = 0,
  thin = 1,
  title = NULL
  ){
  
  # Check input
  if (!inherits(fit, c("ame", "lame"))) {
    stop("fit must be an object of class 'ame' or 'lame'")
  }
  
  # Match params argument
  params <- match.arg(params)
  
  # Prepare data based on parameter selection
  plot_data <- data.frame()
  
  # Add regression coefficients if requested
  if (params %in% c("all", "beta") && !is.null(fit$BETA) && ncol(fit$BETA) > 0) {
    beta_data <- fit$BETA
    if (burn.in > 0 && burn.in < nrow(beta_data)) {
      beta_data <- beta_data[-(1:burn.in), , drop = FALSE]
    }
    
    # Apply thinning
    if (thin > 1) {
      keep_idx <- seq(1, nrow(beta_data), by = thin)
      beta_data <- beta_data[keep_idx, , drop = FALSE]
    }
    
    # Reshape for plotting
    for (j in 1:ncol(beta_data)) {
      param_name <- colnames(beta_data)[j]
      if (is.null(param_name)) param_name <- paste0("beta[", j, "]")
      
      plot_data <- rbind(plot_data, data.frame(
        iteration = 1:nrow(beta_data),
        value = beta_data[, j],
        parameter = param_name,
        type = "Regression",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Add variance components if requested
  if (params %in% c("all", "variance") && !is.null(fit$VC)) {
    vc_data <- fit$VC
    if (burn.in > 0 && burn.in < nrow(vc_data)) {
      vc_data <- vc_data[-(1:burn.in), , drop = FALSE]
    }
    
    # Apply thinning
    if (thin > 1) {
      keep_idx <- seq(1, nrow(vc_data), by = thin)
      vc_data <- vc_data[keep_idx, , drop = FALSE]
    }
    
    # Reshape for plotting
    for (j in 1:ncol(vc_data)) {
      param_name <- colnames(vc_data)[j]
      if (is.null(param_name)) param_name <- paste0("vc[", j, "]")
      
      plot_data <- rbind(plot_data, data.frame(
        iteration = 1:nrow(vc_data),
        value = vc_data[, j],
        parameter = param_name,
        type = "Variance",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Apply include/exclude filters
  if (!is.null(include)) {
    plot_data <- plot_data[plot_data$parameter %in% include, ]
  }
  if (!is.null(exclude)) {
    plot_data <- plot_data[!(plot_data$parameter %in% exclude), ]
  }
  
  # Check if any data remains
  if (nrow(plot_data) == 0) {
    stop("No parameters to plot after filtering")
  }
  
  # Create trace plots
  trace_plots <- ggplot(plot_data, aes(x = iteration, y = value)) +
    geom_line(alpha = 0.7, color = "steelblue") +
    facet_wrap(~ parameter, scales = "free_y", 
              ncol = ncol, nrow = nrow) +
    labs(
      x = "Iteration",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      panel.spacing = unit(0.5, "lines")
    )
  
  # Create density plots
  density_plots <- ggplot(plot_data, aes(x = value)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    facet_wrap(~ parameter, scales = "free", 
              ncol = ncol, nrow = nrow) +
    labs(
      x = "Value",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 9, face = "bold"),
      panel.spacing = unit(0.5, "lines")
    )
  
  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p <- trace_plots / density_plots
    
    # Add title
    if (!is.null(title)) {
      p <- p + patchwork::plot_annotation(title = title)
    } else {
      p <- p + patchwork::plot_annotation(
        title = "MCMC Diagnostics: Trace Plots (top) and Posterior Densities (bottom)"
      )
    }
  } else {
    warning("Install 'patchwork' for combined trace and density plots. Returning trace plots only.")
    p <- trace_plots
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle("MCMC Trace Plots")
    }
  }
  
  return(p)
}