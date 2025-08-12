#' Print methods for AME and LAME simulation objects
#' 
#' @param x simulation object of class "ame.sim" or "lame.sim"
#' @param ... additional arguments (not used)
#' @return the simulation object invisibly
#' @author Shahryar Minhas
#' @importFrom cli cli_h1 cli_h2 cli_text cli_ul cli_li cli_end cli_alert_info
#' @importFrom cli cli_alert_success col_blue col_green col_red style_bold
#' @method print ame.sim
#' @export
print.ame.sim <- function(x, ...) {
  cli::cli_h1("AME Network Simulations")
  
  # Basic info box
  cli::cli_div(theme = list(span.field = list(color = "blue", "font-weight" = "bold")))
  
  cli::cli_ul()
  cli::cli_li("{.field Networks simulated}: {.val {length(x$Y)}}")
  cli::cli_li("{.field Network mode}: {.emph {x$mode}}")
  cli::cli_li("{.field Model family}: {.emph {x$family}}")
  
  # Get dimensions
  if (length(x$Y) > 0) {
    dims <- dim(x$Y[[1]])
    if (x$mode == "bipartite") {
      cli::cli_li("{.field Dimensions}: {.val {dims[1]}} {cli::symbol$times} {.val {dims[2]}} (bipartite)")
    } else {
      cli::cli_li("{.field Dimensions}: {.val {dims[1]}} {cli::symbol$times} {.val {dims[2]}}")
    }
  }
  
  if (!is.null(x$Z)) {
    cli::cli_li("{.field Latent Z matrices}: {col_green(cli::symbol$tick)} included")
  }
  cli::cli_end()
  
  # Usage hint
  cli::cli_alert_info("Access networks: {.code x$Y[[i]]}")
  
  # Quick summary stats if binary
  if (x$family == "binary" && length(x$Y) > 0) {
    densities <- sapply(x$Y[1:min(5, length(x$Y))], 
                       function(y) mean(y == 1, na.rm = TRUE))
    cli::cli_text("")
    cli::cli_text("{.field Density} (first {min(5, length(x$Y))} networks): {.val {round(densities, 3)}}")
  }
  
  invisible(x)
}

#' @method print lame.sim
#' @export
print.lame.sim <- function(x, ...) {
  cli::cli_h1("LAME Longitudinal Network Simulations")
  
  # Basic info
  cli::cli_div(theme = list(span.field = list(color = "blue", "font-weight" = "bold")))
  
  cli::cli_ul()
  cli::cli_li("{.field Trajectories simulated}: {.val {length(x$Y)}}")
  cli::cli_li("{.field Time periods}: {.val {x$n_time}} per trajectory")
  cli::cli_li("{.field Network mode}: {.emph {x$mode}}")
  cli::cli_li("{.field Model family}: {.emph {x$family}}")
  
  # Get dimensions
  if (length(x$Y) > 0 && length(x$Y[[1]]) > 0) {
    dims <- dim(x$Y[[1]][[1]])
    if (x$mode == "bipartite") {
      cli::cli_li("{.field Dimensions}: {.val {dims[1]}} {cli::symbol$times} {.val {dims[2]}} (bipartite)")
    } else {
      cli::cli_li("{.field Dimensions}: {.val {dims[1]}} {cli::symbol$times} {.val {dims[2]}}")
    }
  }
  
  if (!is.null(x$Z)) {
    cli::cli_li("{.field Latent Z matrices}: {col_green(cli::symbol$tick)} included")
  }
  cli::cli_end()
  
  # Total networks
  total_nets <- length(x$Y) * x$n_time
  cli::cli_text("")
  cli::cli_alert_success("Total networks generated: {.strong {total_nets}}")
  
  # Usage hints
  cli::cli_text("")
  cli::cli_alert_info("Access trajectories: {.code x$Y[[i]][[t]]}")
  cli::cli_text("  {cli::symbol$arrow_right} {.emph i} = trajectory index (1-{length(x$Y)})")
  cli::cli_text("  {cli::symbol$arrow_right} {.emph t} = time index (1-{x$n_time})")
  
  invisible(x)
}

#' Summary method for AME simulations
#' 
#' @param object simulation object of class "ame.sim"
#' @param ... additional arguments (not used)
#' @return summary statistics invisibly
#' @importFrom cli cli_h2 cli_alert_info format_inline
#' @importFrom stats quantile
#' @method summary ame.sim
#' @export
summary.ame.sim <- function(object, ...) {
  x <- object
  
  cli::cli_h2("Summary of AME Network Simulations")
  
  # Basic info
  cli::cli_alert_info("Simulations: {.val {length(x$Y)}}, Family: {.emph {x$family}}, Mode: {.emph {x$mode}}")
  
  # Calculate summary statistics
  if (x$family == "binary") {
    densities <- sapply(x$Y, function(y) {
      sum(y == 1, na.rm = TRUE) / sum(!is.na(y))
    })
    
    cli::cli_text("")
    cli::cli_text("{style_bold('Network Density Statistics:')}")
    
    # Create a nice summary table
    q <- quantile(densities, c(0.025, 0.25, 0.5, 0.75, 0.975))
    
    cli::cli_ul()
    cli::cli_li("Mean: {.val {round(mean(densities), 4)}} (SD: {.val {round(sd(densities), 4)}})")
    cli::cli_li("Range: [{.val {round(min(densities), 4)}}, {.val {round(max(densities), 4)}}]")
    cli::cli_li("95% CI: [{.val {round(q[1], 4)}}, {.val {round(q[5], 4)}}]")
    cli::cli_li("IQR: [{.val {round(q[2], 4)}}, {.val {round(q[4], 4)}}]")
    cli::cli_end()
    
    # Visual density bar
    cli::cli_text("")
    cli::cli_text("{style_bold('Density Distribution:')}")
    
    # Create a simple ASCII histogram
    hist_data <- cut(densities, breaks = 10)
    hist_table <- table(hist_data)
    max_count <- max(hist_table)
    
    for(i in 1:length(hist_table)) {
      bar_length <- round(30 * hist_table[i] / max_count)
      bar <- paste(rep("\u2588", bar_length), collapse = "")
      range_label <- names(hist_table)[i]
      cli::cli_text("  {.field {sprintf('%-15s', range_label)}} {bar} {.val {hist_table[i]}}")
    }
  }
  
  # For continuous networks
  if (x$family %in% c("normal", "tobit")) {
    means <- sapply(x$Y, mean, na.rm = TRUE)
    vars <- sapply(x$Y, var, na.rm = TRUE)
    
    cli::cli_text("")
    cli::cli_text("{style_bold('Network Value Statistics:')}")
    
    cli::cli_ul()
    cli::cli_li("Mean across networks: {.val {round(mean(means), 4)}} (SD: {.val {round(sd(means), 4)}})")
    cli::cli_li("Variance across networks: {.val {round(mean(vars), 4)}} (SD: {.val {round(sd(vars), 4)}})")
    cli::cli_end()
  }
  
  # Count statistics
  if (x$family == "poisson") {
    edge_counts <- sapply(x$Y, sum, na.rm = TRUE)
    cli::cli_text("")
    cli::cli_text("{style_bold('Edge Count Statistics:')}")
    cli::cli_ul()
    cli::cli_li("Mean edges: {.val {round(mean(edge_counts), 1)}}")
    cli::cli_li("SD: {.val {round(sd(edge_counts), 1)}}")
    cli::cli_li("Range: [{.val {min(edge_counts)}}, {.val {max(edge_counts)}}]")
    cli::cli_end()
  }
  
  invisible(list(
    n_sims = length(x$Y),
    mode = x$mode,
    family = x$family,
    densities = if(x$family == "binary") densities else NULL
  ))
}

#' Summary method for LAME simulations
#' 
#' @param object simulation object of class "lame.sim"
#' @param ... additional arguments (not used)
#' @return summary statistics
#' @method summary lame.sim
#' @export
summary.lame.sim <- function(object, ...) {
  x <- object
  
  cli::cli_h2("Summary of LAME Longitudinal Network Simulations")
  
  # Calculate summary statistics across simulations and time
  n_sims <- length(x$Y)
  n_time <- x$n_time
  
  # Basic info
  cli::cli_alert_info("Trajectories: {.val {n_sims}}, Time periods: {.val {n_time}}, Family: {.emph {x$family}}, Mode: {.emph {x$mode}}")
  
  # For binary networks, calculate density over time
  if (x$family == "binary") {
    # Matrix: rows are simulations, columns are time points
    density_matrix <- matrix(NA, n_sims, n_time)
    
    for (i in 1:n_sims) {
      for (t in 1:n_time) {
        y <- x$Y[[i]][[t]]
        density_matrix[i, t] <- sum(y == 1, na.rm = TRUE) / sum(!is.na(y))
      }
    }
    
    # Average density at each time point
    mean_density <- colMeans(density_matrix)
    sd_density <- apply(density_matrix, 2, sd)
    
    cli::cli_text("")
    cli::cli_text("{style_bold('Network Density Over Time:')}")
    
    cli::cli_ul()
    cli::cli_li("Time 1: {.val {round(mean_density[1], 4)}} (SD: {.val {round(sd_density[1], 4)}})")
    if (n_time > 2) {
      mid_t <- floor(n_time/2)
      cli::cli_li("Time {mid_t}: {.val {round(mean_density[mid_t], 4)}} (SD: {.val {round(sd_density[mid_t], 4)}})")
    }
    cli::cli_li("Time {n_time}: {.val {round(mean_density[n_time], 4)}} (SD: {.val {round(sd_density[n_time], 4)}})")
    cli::cli_end()
    
    # Test for temporal trend
    time_vec <- 1:n_time
    trends <- apply(density_matrix, 1, function(d) {
      if (any(is.na(d))) return(NA)
      cor(d, time_vec)
    })
    
    mean_trend <- mean(trends, na.rm = TRUE)
    prop_increasing <- mean(trends > 0.1, na.rm = TRUE)
    prop_decreasing <- mean(trends < -0.1, na.rm = TRUE)
    
    cli::cli_text("")
    cli::cli_text("{style_bold('Temporal Trends:')}")
    
    cli::cli_ul()
    cli::cli_li("Mean correlation with time: {.val {round(mean_trend, 4)}}")
    cli::cli_li("Proportion increasing: {.val {round(prop_increasing, 3)}} {col_green(paste0('(', round(prop_increasing*100, 1), '%)'))}")
    cli::cli_li("Proportion decreasing: {.val {round(prop_decreasing, 3)}} {col_red(paste0('(', round(prop_decreasing*100, 1), '%)'))}")
    cli::cli_li("Proportion stable: {.val {round(1 - prop_increasing - prop_decreasing, 3)}}")
    cli::cli_end()
    
    # Visualize trend with sparkline-like display
    if (n_time <= 20) {
      cli::cli_text("")
      cli::cli_text("{style_bold('Density Trajectory (mean \u00b1 SD):')}")
      
      # Create a simple text-based plot
      for (t in 1:n_time) {
        bar_val <- round(mean_density[t] * 20)  # Scale to ~20 chars
        bar <- paste(rep("\u2593", bar_val), collapse = "")
        time_label <- sprintf("T%02d", t)
        density_label <- sprintf("%.3f", mean_density[t])
        cli::cli_text("  {.field {time_label}}: {bar} {.val {density_label}}")
      }
    }
  }
  
  # For continuous networks
  if (x$family %in% c("normal", "tobit")) {
    cli::cli_text("")
    cli::cli_text("{style_bold('Network Value Statistics Over Time:')}")
    
    # Calculate mean values over time
    mean_matrix <- matrix(NA, n_sims, n_time)
    for (i in 1:n_sims) {
      for (t in 1:n_time) {
        mean_matrix[i, t] <- mean(x$Y[[i]][[t]], na.rm = TRUE)
      }
    }
    
    mean_by_time <- colMeans(mean_matrix, na.rm = TRUE)
    
    cli::cli_ul()
    cli::cli_li("Mean at T1: {.val {round(mean_by_time[1], 4)}}")
    cli::cli_li("Mean at T{n_time}: {.val {round(mean_by_time[n_time], 4)}}")
    cli::cli_li("Overall mean: {.val {round(mean(mean_matrix, na.rm = TRUE), 4)}}")
    cli::cli_end()
  }
  
  # Summary of variation
  cli::cli_text("")
  cli::cli_alert_success("Total networks generated: {.strong {n_sims * n_time}}")
  
  invisible(list(
    n_sims = n_sims,
    n_time = n_time,
    mode = x$mode,
    family = x$family,
    density_matrix = if(x$family == "binary") density_matrix else NULL
  ))
}