#' Plot comprehensive diagnostics for a LAME model fit
#' 
#' Creates a comprehensive set of diagnostic plots for a LAME (Longitudinal 
#' Additive and Multiplicative Effects) model, including MCMC diagnostics, 
#' parameter evolution over time, and longitudinal goodness-of-fit checks. 
#' This is the default plot method for LAME objects.
#' 
#' @details
#' The function produces a multi-panel plot containing:
#' \describe{
#'   \item{MCMC trace plots}{Shows mixing and convergence of key parameters}
#'   \item{Posterior distributions}{Density plots of regression coefficients 
#'         and variance components}
#'   \item{Longitudinal GOF}{Time series of observed network statistics with
#'         posterior predictive intervals}
#'   \item{Effects over time}{Evolution of additive effects across time periods
#'         (if applicable)}
#'   \item{Network snapshots}{Visualization of network at selected time points}
#' }
#' 
#' The plot adapts to the longitudinal structure:
#' - Shows temporal trends in network statistics
#' - Highlights composition changes if actors enter/exit
#' - Displays credible intervals for time-varying statistics
#' 
#' @param x an object of class "lame" from fitting a LAME model
#' @param which numeric or character vector specifying which plots to produce:
#'        1 or "trace" = MCMC trace plots,
#'        2 or "density" = posterior density plots,
#'        3 or "gof" = longitudinal goodness-of-fit plots,
#'        4 or "effects" = additive and multiplicative effects,
#'        5 or "network" = network snapshots at selected times.
#'        Default is c(1,2,3,4) to show main diagnostic plots.
#' @param time.points numeric vector of time points for network snapshots
#'        (only used if "network" in which). Default is c(1, middle, last).
#' @param ask logical; if TRUE, user is prompted before each plot page
#' @param pages character string specifying how to arrange plots:
#'        "single" = one comprehensive page (default),
#'        "multiple" = separate pages for each plot type
#' @param ... additional arguments (currently not used)
#' @return NULL (invisibly). Plots are displayed as side effects.
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @seealso \code{\link{lame}}, \code{\link{trace_plot}}, \code{\link{gof_plot}}, 
#'          \code{\link{ab_plot}}, \code{\link{uv_plot}}
#' @examples
#' \dontrun{
#' # Fit a LAME model
#' fit <- lame(Y_list, X_list, R = 2)
#' 
#' # Default comprehensive plot
#' plot(fit)
#' 
#' # Only MCMC diagnostics
#' plot(fit, which = c("trace", "density"))
#' 
#' # Include network snapshots at specific times
#' plot(fit, which = c(3, 4, 5), time.points = c(1, 5, 10))
#' 
#' # Separate pages for each plot type
#' plot(fit, pages = "multiple")
#' }
#' @method plot lame
#' @export
#' @import ggplot2
#' @import patchwork
plot.lame <- function(x, 
                      which = c(1, 2, 3, 4),
                      time.points = NULL,
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
    if (5 %in% which) plot_types <- c(plot_types, "network")
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
    
    # Add key regression coefficients
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
  
  # 3. Longitudinal Goodness-of-fit plots
  if ("gof" %in% plot_types) {
    if (!is.null(fit$GOF_T)) {
      # Time-indexed GOF data
      gof_data_list <- fit$GOF_T
      n_time <- length(gof_data_list)
      
      # Prepare longitudinal GOF data
      gof_long_data <- data.frame()
      
      for (t in 1:n_time) {
        gof_t <- gof_data_list[[t]]
        if (!is.null(gof_t) && nrow(gof_t) > 1) {
          obs_t <- gof_t[1, ]
          pred_t <- gof_t[-1, , drop = FALSE]
          
          for (j in 1:ncol(gof_t)) {
            stat_name <- colnames(gof_t)[j]
            gof_long_data <- rbind(gof_long_data, data.frame(
              time = t,
              statistic = stat_name,
              observed = obs_t[j],
              median = median(pred_t[, j]),
              lower = quantile(pred_t[, j], 0.025),
              upper = quantile(pred_t[, j], 0.975),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      if (nrow(gof_long_data) > 0) {
        p_gof <- ggplot(gof_long_data, aes(x = time)) +
          geom_ribbon(aes(ymin = lower, ymax = upper), 
                     alpha = 0.3, fill = "blue") +
          geom_line(aes(y = median), color = "blue", 
                   linewidth = 1, linetype = "dashed") +
          geom_line(aes(y = observed), color = "red", linewidth = 1) +
          geom_point(aes(y = observed), color = "red", size = 2) +
          facet_wrap(~ statistic, scales = "free_y", ncol = 2) +
          labs(title = "Longitudinal GOF: Observed (red) vs 95% Credible Intervals",
               x = "Time Period", y = "Value") +
          theme_minimal() +
          theme(strip.text = element_text(size = 9))
        
        plot_list[["gof"]] <- p_gof
      }
    } else if (!is.null(fit$GOF) && is.array(fit$GOF) && length(dim(fit$GOF)) == 3) {
      gof_long_data <- data.frame()
      n_stats <- dim(fit$GOF)[1]
      n_time <- dim(fit$GOF)[2]
      n_iter <- dim(fit$GOF)[3]
      
      stat_names <- dimnames(fit$GOF)[[1]]
      time_names <- dimnames(fit$GOF)[[2]]
      if (is.null(time_names)) time_names <- 1:n_time
      
      for (t in 1:n_time) {
        obs_t <- fit$GOF[, t, 1]
        
        if (n_iter > 1) {
          for (s in 1:n_stats) {
            stat_name <- stat_names[s]
            pred_vals <- fit$GOF[s, t, 2:n_iter]
            
            gof_long_data <- rbind(gof_long_data, data.frame(
              time = time_names[t],
              statistic = stat_name,
              observed = obs_t[s],
              median = median(pred_vals),
              lower = quantile(pred_vals, 0.025),
              upper = quantile(pred_vals, 0.975),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      if (nrow(gof_long_data) > 0) {
        p_gof <- ggplot(gof_long_data, aes(x = as.numeric(factor(time)))) +
          geom_ribbon(aes(ymin = lower, ymax = upper), 
                     alpha = 0.3, fill = "blue") +
          geom_line(aes(y = median), color = "blue", 
                   linewidth = 1, linetype = "dashed") +
          geom_line(aes(y = observed), color = "red", linewidth = 1) +
          geom_point(aes(y = observed), color = "red", size = 2) +
          facet_wrap(~ statistic, scales = "free_y", ncol = 2) +
          labs(title = "Longitudinal GOF: Observed (red) vs 95% Credible Intervals",
               x = "Time Period", y = "Value") +
          scale_x_continuous(breaks = 1:n_time, labels = time_names) +
          theme_minimal() +
          theme(strip.text = element_text(size = 9))
        
        plot_list[["gof"]] <- p_gof
      }
    } else if (!is.null(fit$GOF) && is.list(fit$GOF) && !is.data.frame(fit$GOF)) {
      # Handle GOF as named list (new format)
      gof_long_data <- data.frame()
      stat_names <- names(fit$GOF)
      
      # Check if we have time-indexed data (multiple columns)
      if (ncol(fit$GOF[[1]]) > 1) {
        n_time <- ncol(fit$GOF[[1]])
        
        for (t in 1:n_time) {
          for (stat_name in stat_names) {
            vals <- fit$GOF[[stat_name]][, t]
            if (length(vals) > 1) {
              obs_val <- vals[1]
              pred_vals <- vals[-1]
              
              gof_long_data <- rbind(gof_long_data, data.frame(
                time = t,
                statistic = stat_name,
                observed = obs_val,
                median = median(pred_vals),
                lower = quantile(pred_vals, 0.025),
                upper = quantile(pred_vals, 0.975),
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        if (nrow(gof_long_data) > 0) {
          p_gof <- ggplot(gof_long_data, aes(x = time)) +
            geom_ribbon(aes(ymin = lower, ymax = upper), 
                       alpha = 0.3, fill = "blue") +
            geom_line(aes(y = median), color = "blue", 
                     linewidth = 1, linetype = "dashed") +
            geom_line(aes(y = observed), color = "red", linewidth = 1) +
            geom_point(aes(y = observed), color = "red", size = 2) +
            facet_wrap(~ statistic, scales = "free_y", ncol = 2) +
            labs(title = "Longitudinal GOF: Observed (red) vs 95% Credible Intervals",
                 x = "Time Period", y = "Value") +
            theme_minimal() +
            theme(strip.text = element_text(size = 9))
          
          plot_list[["gof"]] <- p_gof
        }
      } else {
        # Single time point - create static GOF plot
        gof_plot_data <- data.frame()
        for (stat_name in stat_names) {
          vals <- fit$GOF[[stat_name]][, 1]
          if (length(vals) > 1) {
            obs_val <- vals[1]
            pred_vals <- vals[-1]
            
            gof_plot_data <- rbind(gof_plot_data, data.frame(
              value = pred_vals,
              statistic = stat_name,
              observed = obs_val,
              stringsAsFactors = FALSE
            ))
          }
        }
        
        if (nrow(gof_plot_data) > 0) {
          p_gof <- ggplot(gof_plot_data, aes(x = value)) +
            geom_histogram(aes(y = after_stat(density)), 
                          bins = 30, fill = "lightblue", 
                          color = "white", alpha = 0.7) +
            geom_vline(aes(xintercept = observed), 
                      color = "red", linewidth = 1) +
            facet_wrap(~ statistic, scales = "free", ncol = 2) +
            labs(title = "GOF: Observed (red) vs Posterior Predictive",
                 x = "Value", y = "Density") +
            theme_minimal() +
            theme(strip.text = element_text(size = 9))
          
          plot_list[["gof"]] <- p_gof
        }
      }
    } else if (!is.null(fit$GOF) && is.matrix(fit$GOF) && nrow(fit$GOF) > 1) {
      # Fall back to static GOF if no time-indexed version (2D array)
      gof_data <- fit$GOF
      obs_vals <- gof_data[1, ]
      pred_vals <- gof_data[-1, , drop = FALSE]
      
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
  }
  
  # 4. Effects plots
  if ("effects" %in% plot_types) {
    effects_plots <- list()
    
    # Additive effects
    if (!is.null(fit$APM) && !is.null(fit$BPM)) {
      # Check if we have time-varying effects
      if (!is.null(fit$a_dynamic) && !is.null(fit$b_dynamic)) {
        # Time-varying effects using a_dynamic/b_dynamic
        A <- fit$a_dynamic
        B <- fit$b_dynamic
        
        # Convert to long format using as.table for efficiency
        df_a <- as.data.frame(as.table(A))
        names(df_a) <- c("actor", "time", "effect")
        df_a$type <- "Sender"
        df_a$time <- as.integer(df_a$time)
        
        df_b <- as.data.frame(as.table(B))
        names(df_b) <- c("actor", "time", "effect")
        df_b$type <- "Receiver"
        df_b$time <- as.integer(df_b$time)
        
        effects_time_data <- rbind(df_a, df_b)
        names(effects_time_data)[1] <- "name"  # Rename actor to name for consistency
      } else if (!is.null(fit$APM_T)) {
        # Legacy support for APM_T/BPM_T format
        n_time <- length(fit$APM_T)
        effects_time_data <- data.frame()
        
        for (t in 1:n_time) {
          if (!is.null(fit$APM_T[[t]])) {
            for (i in 1:length(fit$APM_T[[t]])) {
              effects_time_data <- rbind(effects_time_data, data.frame(
                time = t,
                effect = fit$APM_T[[t]][i],
                name = names(fit$APM_T[[t]])[i],
                type = "Sender",
                stringsAsFactors = FALSE
              ))
            }
          }
          if (!is.null(fit$BPM_T[[t]])) {
            for (i in 1:length(fit$BPM_T[[t]])) {
              effects_time_data <- rbind(effects_time_data, data.frame(
                time = t,
                effect = fit$BPM_T[[t]][i],
                name = names(fit$BPM_T[[t]])[i],
                type = "Receiver",
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        if (nrow(effects_time_data) > 0) {
          # Select top actors by variance across time
          actor_var <- aggregate(effect ~ name + type, effects_time_data, var)
          top_actors <- head(actor_var[order(-actor_var$effect), ], 20)
          effects_time_data <- effects_time_data[
            effects_time_data$name %in% top_actors$name, 
          ]
          
          p_effects_time <- ggplot(effects_time_data, 
                                   aes(x = time, y = effect, 
                                       color = name, linetype = type)) +
            geom_line(alpha = 0.7) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
            labs(title = "Top 20 Actors: Effects Over Time",
                 x = "Time Period", y = "Effect",
                 color = "Actor", linetype = "Type") +
            theme_minimal() +
            theme(legend.position = "right")
          
          effects_plots[["time_varying"]] <- p_effects_time
        }
      } else {
        # Static effects (averaged over time)
        sender_data <- data.frame(
          effect = fit$APM,
          name = names(fit$APM),
          type = "Sender",
          stringsAsFactors = FALSE
        )
        
        receiver_data <- data.frame(
          effect = fit$BPM,
          name = names(fit$BPM),
          type = "Receiver",
          stringsAsFactors = FALSE
        )
        
        effects_data <- rbind(sender_data, receiver_data)
        effects_data <- effects_data[order(effects_data$type, effects_data$effect), ]
        effects_data$index <- 1:nrow(effects_data)
        
        p_additive <- ggplot(effects_data, aes(x = index, y = effect, color = type)) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
          geom_point(size = 1.5, alpha = 0.7) +
          geom_segment(aes(xend = index, yend = 0), alpha = 0.5) +
          scale_color_manual(values = c("Sender" = "steelblue", "Receiver" = "darkred")) +
          labs(title = "Additive Effects (Time-Averaged)", 
               x = "Actor (sorted)", y = "Effect", color = "Type") +
          theme_minimal() +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.position = "bottom")
        
        effects_plots[["additive"]] <- p_additive
      }
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
  
  # 5. Network snapshots (optional)
  if ("network" %in% plot_types && !is.null(fit$Y_T)) {
    n_time <- length(fit$Y_T)
    
    # Determine time points to plot
    if (is.null(time.points)) {
      time.points <- c(1, ceiling(n_time/2), n_time)
      time.points <- unique(time.points[time.points <= n_time])
    }
    time.points <- time.points[time.points >= 1 & time.points <= n_time]
    
    network_plots <- list()
    
    for (t in time.points) {
      Y_t <- fit$Y_T[[t]]
      if (!is.null(Y_t)) {
        # Create network data
        net_data <- data.frame()
        n <- nrow(Y_t)
        for (i in 1:n) {
          for (j in 1:n) {
            if (i != j && !is.na(Y_t[i,j]) && Y_t[i,j] > 0) {
              net_data <- rbind(net_data, data.frame(
                from = i,
                to = j,
                weight = Y_t[i,j],
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        if (nrow(net_data) > 0) {
          # Simple circular layout
          angles <- seq(0, 2*pi, length.out = n+1)[1:n]
          node_positions <- data.frame(
            node = 1:n,
            x = cos(angles),
            y = sin(angles)
          )
          
          # Merge positions
          net_data <- merge(net_data, node_positions, by.x = "from", by.y = "node")
          names(net_data)[names(net_data) %in% c("x", "y")] <- c("x_from", "y_from")
          net_data <- merge(net_data, node_positions, by.x = "to", by.y = "node")
          names(net_data)[names(net_data) %in% c("x", "y")] <- c("x_to", "y_to")
          
          p_net <- ggplot() +
            geom_segment(data = net_data, 
                        aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
                        alpha = 0.3, color = "gray50") +
            geom_point(data = node_positions, 
                      aes(x = x, y = y), 
                      size = 3, color = "steelblue") +
            labs(title = paste("Network at Time", t)) +
            theme_void() +
            coord_fixed()
          
          network_plots[[paste0("t", t)]] <- p_net
        }
      }
    }
    
    if (length(network_plots) > 0) {
      if (has_patchwork) {
        plot_list[["network"]] <- patchwork::wrap_plots(network_plots, ncol = length(network_plots))
      } else {
        plot_list[["network"]] <- network_plots[[1]]
      }
    }
  }
  
  # Display plots
  if (pages == "single" && has_patchwork) {
    # Combine all plots into one page
    if (length(plot_list) > 0) {
      combined_plot <- patchwork::wrap_plots(plot_list, ncol = 1) +
        patchwork::plot_annotation(
          title = "LAME Model Diagnostics",
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