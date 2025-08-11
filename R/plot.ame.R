#' Plot comprehensive diagnostics for an AME model fit
#' 
#' Creates a comprehensive set of diagnostic plots for an AME (Additive and 
#' Multiplicative Effects) model, including MCMC diagnostics, parameter estimates,
#' and goodness-of-fit checks. This is the default plot method for AME objects.
#' 
#' @details
#' The function produces a multi-panel plot containing:
#' \describe{
#'   \item{MCMC trace plots}{Shows mixing and convergence of key parameters}
#'   \item{Posterior distributions}{Density plots of regression coefficients 
#'         and variance components}
#'   \item{Goodness-of-fit}{Comparison of observed network statistics to 
#'         posterior predictive distributions (if gof=TRUE in model fit)}
#'   \item{Effect visualizations}{Plots of additive sender/receiver effects
#'         and multiplicative effects if present}
#' }
#' 
#' The plot adapts to the model specification:
#' - Shows only relevant variance components (e.g., omits cab for symmetric networks)
#' - Includes multiplicative effects plots only if R > 0
#' - Includes GOF plots only if computed during model fitting
#' 
#' @param x an object of class "ame" from fitting an AME model
#' @param which numeric or character vector specifying which plots to produce:
#'        1 or "trace" = MCMC trace plots,
#'        2 or "density" = posterior density plots,
#'        3 or "gof" = goodness-of-fit plots,
#'        4 or "effects" = additive and multiplicative effects.
#'        Default is c(1,2,3,4) to show all plots.
#' @param ask logical; if TRUE, user is prompted before each plot page
#' @param pages character string specifying how to arrange plots:
#'        "single" = one comprehensive page (default),
#'        "multiple" = separate pages for each plot type
#' @param ... additional arguments (currently not used)
#' @return NULL (invisibly). Plots are displayed as side effects.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @seealso \code{\link{ame}}, \code{\link{trace_plot}}, \code{\link{gof_plot}}, 
#'          \code{\link{ab_plot}}, \code{\link{uv_plot}}
#' @examples
#' \dontrun{
#' # Fit an AME model
#' fit <- ame(Y, X, R = 2, gof = TRUE)
#' 
#' # Default comprehensive plot
#' plot(fit)
#' 
#' # Only MCMC diagnostics
#' plot(fit, which = c("trace", "density"))
#' 
#' # Only effects plots
#' plot(fit, which = "effects")
#' 
#' # Separate pages for each plot type
#' plot(fit, pages = "multiple")
#' }
#' @method plot ame
#' @export
#' @import ggplot2
#' @import patchwork
plot.ame <- function(x, 
                     which = c(1, 2, 3, 4),
                     ask = FALSE,
                     pages = c("single", "multiple"),
                     ...) {
  
  fit <- x
  pages <- match.arg(pages)
  
  # Parse which argument
  plot_types <- c()
  if (is.numeric(which)) {
    if (1 %in% which) plot_types <- c(plot_types, "trace")
    if (2 %in% which) plot_types <- c(plot_types, "density")
    if (3 %in% which) plot_types <- c(plot_types, "gof")
    if (4 %in% which) plot_types <- c(plot_types, "effects")
  } else {
    plot_types <- which
  }
  
  # Check for patchwork
  has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
  if (!has_patchwork && pages == "single") {
    warning("Package 'patchwork' not installed. Switching to multiple pages.")
    pages <- "multiple"
  }
  
  # Store plots
  plot_list <- list()
  
  # 1. MCMC Trace plots
  if ("trace" %in% plot_types) {
    # Prepare data for trace plots
    trace_data <- data.frame()
    
    # Add variance components
    if (!is.null(fit$VC)) {
      vc_names <- colnames(fit$VC)
      for (i in 1:ncol(fit$VC)) {
        trace_data <- rbind(trace_data, data.frame(
          iteration = 1:nrow(fit$VC),
          value = fit$VC[, i],
          parameter = vc_names[i],
          type = "Variance",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add key regression coefficients (first 6 or all if fewer)
    if (!is.null(fit$BETA) && ncol(fit$BETA) > 0) {
      n_beta <- min(6, ncol(fit$BETA))
      beta_names <- colnames(fit$BETA)[1:n_beta]
      for (i in 1:n_beta) {
        trace_data <- rbind(trace_data, data.frame(
          iteration = 1:nrow(fit$BETA),
          value = fit$BETA[, i],
          parameter = beta_names[i],
          type = "Regression",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if (nrow(trace_data) > 0) {
      p_trace <- ggplot(trace_data, aes(x = iteration, y = value)) +
        geom_line(alpha = 0.7, color = "steelblue") +
        facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
        labs(title = "MCMC Trace Plots", x = "Iteration", y = "Value") +
        theme_minimal() +
        theme(strip.text = element_text(size = 8))
      
      plot_list[["trace"]] <- p_trace
    }
  }
  
  # 2. Posterior Density plots
  if ("density" %in% plot_types) {
    density_data <- data.frame()
    
    # Add variance components
    if (!is.null(fit$VC)) {
      vc_names <- colnames(fit$VC)
      for (i in 1:ncol(fit$VC)) {
        density_data <- rbind(density_data, data.frame(
          value = fit$VC[, i],
          parameter = vc_names[i],
          type = "Variance",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Add regression coefficients
    if (!is.null(fit$BETA) && ncol(fit$BETA) > 0) {
      n_beta <- min(6, ncol(fit$BETA))
      beta_names <- colnames(fit$BETA)[1:n_beta]
      for (i in 1:n_beta) {
        density_data <- rbind(density_data, data.frame(
          value = fit$BETA[, i],
          parameter = beta_names[i],
          type = "Regression",
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if (nrow(density_data) > 0) {
      p_density <- ggplot(density_data, aes(x = value)) +
        geom_density(fill = "steelblue", alpha = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
        facet_wrap(~ parameter, scales = "free", ncol = 3) +
        labs(title = "Posterior Distributions", x = "Value", y = "Density") +
        theme_minimal() +
        theme(strip.text = element_text(size = 8))
      
      plot_list[["density"]] <- p_density
    }
  }
  
  # 3. Goodness-of-fit plots
  if ("gof" %in% plot_types && !is.null(fit$GOF) && nrow(fit$GOF) > 1) {
    gof_data <- fit$GOF
    obs_vals <- gof_data[1, ]
    pred_vals <- gof_data[-1, , drop = FALSE]
    
    # Create GOF plot data
    gof_plot_data <- data.frame()
    for (j in 1:ncol(gof_data)) {
      stat_name <- colnames(gof_data)[j]
      gof_plot_data <- rbind(gof_plot_data, data.frame(
        value = pred_vals[, j],
        statistic = stat_name,
        observed = obs_vals[j],
        stringsAsFactors = FALSE
      ))
    }
    
    p_gof <- ggplot(gof_plot_data, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), 
                    bins = 30, fill = "lightblue", 
                    color = "white", alpha = 0.7) +
      geom_vline(aes(xintercept = observed), 
                color = "red", linewidth = 1) +
      facet_wrap(~ statistic, scales = "free", ncol = 2) +
      labs(title = "Goodness-of-Fit: Observed (red) vs Posterior Predictive",
           x = "Value", y = "Density") +
      theme_minimal() +
      theme(strip.text = element_text(size = 9))
    
    plot_list[["gof"]] <- p_gof
  }
  
  # 4. Effects plots
  if ("effects" %in% plot_types) {
    effects_plots <- list()
    
    # Additive effects
    if (!is.null(fit$APM) && !is.null(fit$BPM)) {
      # Sender effects
      sender_data <- data.frame(
        effect = fit$APM,
        name = names(fit$APM),
        type = "Sender",
        stringsAsFactors = FALSE
      )
      
      # Receiver effects
      receiver_data <- data.frame(
        effect = fit$BPM,
        name = names(fit$BPM),
        type = "Receiver",
        stringsAsFactors = FALSE
      )
      
      # Combine
      effects_data <- rbind(sender_data, receiver_data)
      
      # Sort by effect size within each type
      effects_data <- effects_data[order(effects_data$type, effects_data$effect), ]
      effects_data$index <- 1:nrow(effects_data)
      
      p_additive <- ggplot(effects_data, aes(x = index, y = effect, color = type)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_point(size = 1.5, alpha = 0.7) +
        geom_segment(aes(xend = index, yend = 0), alpha = 0.5) +
        scale_color_manual(values = c("Sender" = "steelblue", "Receiver" = "darkred")) +
        labs(title = "Additive Effects", 
             x = "Actor (sorted)", 
             y = "Effect",
             color = "Type") +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom")
      
      effects_plots[["additive"]] <- p_additive
    }
    
    # Multiplicative effects
    if (!is.null(fit$U) && ncol(fit$U) >= 2) {
      mult_data <- data.frame(
        U1 = fit$U[, 1],
        U2 = fit$U[, 2],
        V1 = fit$V[, 1],
        V2 = fit$V[, 2],
        name = rownames(fit$U)
      )
      
      p_mult <- ggplot(mult_data) +
        geom_point(aes(x = U1, y = U2), color = "steelblue", size = 2, alpha = 0.7) +
        geom_point(aes(x = V1, y = V2), color = "darkred", size = 2, alpha = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        labs(title = "Multiplicative Effects (Blue=U, Red=V)",
             x = "Dimension 1", y = "Dimension 2") +
        theme_minimal() +
        coord_fixed()
      
      effects_plots[["multiplicative"]] <- p_mult
    }
    
    # Combine effects plots
    if (length(effects_plots) > 0) {
      if (has_patchwork) {
        plot_list[["effects"]] <- patchwork::wrap_plots(effects_plots, ncol = 1)
      } else {
        plot_list[["effects"]] <- effects_plots[[1]]
      }
    }
  }
  
  # Display plots
  if (pages == "single" && has_patchwork) {
    # Combine all plots into one page
    if (length(plot_list) > 0) {
      combined_plot <- patchwork::wrap_plots(plot_list, ncol = 1) +
        patchwork::plot_annotation(
          title = "AME Model Diagnostics",
          theme = theme(plot.title = element_text(size = 14, face = "bold"))
        )
      print(combined_plot)
    }
  } else {
    # Display plots on separate pages
    if (ask && length(plot_list) > 1) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    
    for (p in plot_list) {
      print(p)
      if (ask && length(plot_list) > 1) {
        readline(prompt = "Press [enter] to continue")
      }
    }
  }
  
  invisible(NULL)
}