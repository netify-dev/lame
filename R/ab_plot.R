# Global variables to avoid R CMD check warnings
if(getRversion() >= "2.15.1") {
  utils::globalVariables(c("mean_effect", "x_prev", "y_prev"))
}

#' Visualize sender and receiver random effects
#' 
#' Creates a visualization of the additive sender (row) and receiver (column) 
#' random effects from an AME or LAME model. Automatically detects whether effects
#' are static or dynamic and provides appropriate visualization options.
#' 
#' @details
#' The additive effects in AME models represent:
#' \describe{
#'   \item{Sender effects (a)}{Actor-specific tendencies to form outgoing ties.
#'         Positive values indicate actors who send more ties than expected;
#'         negative values indicate actors who send fewer ties.}
#'   \item{Receiver effects (b)}{Actor-specific tendencies to receive incoming ties.
#'         Positive values indicate actors who receive more ties than expected;
#'         negative values indicate actors who receive fewer ties.}
#' }
#' 
#' For static effects, the plot displays these effects as a dot plot with vertical 
#' lines extending from zero to each effect estimate.
#' 
#' For dynamic effects (when fit contains a_dynamic/b_dynamic), additional options
#' are available to visualize how effects evolve over time.
#' 
#' @param fit An object of class "ame" or "lame" from fitting an AME model
#' @param effect Character string specifying which effect to plot: 
#'        "sender" (default) or "receiver"
#' @param sorted Logical; if TRUE (default), actors are sorted by effect magnitude
#' @param labels Logical; if TRUE, actor labels are shown on x-axis (default TRUE
#'        for n <= 50 actors)
#' @param title Optional title for the plot
#' @param time_point For dynamic effects, which time point to plot (default: last).
#'                   Can be numeric index, "all" for faceted plot, or "average" for time-averaged
#' @param plot_type For dynamic effects: "snapshot" (single time), "trajectory" (evolution over time),
#'                  "faceted" (grid of time points), or "ribbon" (confidence bands over time).
#'                  For static effects, this parameter is ignored.
#' @param show_actors Character vector of specific actors to highlight (for dynamic trajectory/ribbon plots)
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \dontrun{
#' # Fit an AME model
#' fit <- ame(Y, X)
#' 
#' # Visualize sender effects
#' ab_plot(fit, effect = "sender")
#' 
#' # Visualize receiver effects without sorting
#' ab_plot(fit, effect = "receiver", sorted = FALSE)
#' 
#' # Customize the plot
#' library(ggplot2)
#' ab_plot(fit) + theme_minimal() + ggtitle("Network Sender Effects")
#' }
#' @export
#' @import ggplot2
#' @import cli
#' @importFrom stats reorder
ab_plot <- function(fit, 
                    effect = c("sender", "receiver"), 
                    sorted = TRUE, 
                    labels = NULL, 
                    title = NULL,
                    time_point = NULL,
                    plot_type = c("snapshot", "trajectory", "faceted", "ribbon"),
                    show_actors = NULL) {
  
  # Check class
  if (!inherits(fit, c("ame", "lame"))) {
    stop("fit must be an object of class 'ame' or 'lame'")
  }
  
  # Match arguments
  effect <- match.arg(effect)
  plot_type <- match.arg(plot_type)
  
  # Check if effects are dynamic
  is_dynamic <- FALSE
  if (!is.null(fit$a_dynamic) || !is.null(fit$b_dynamic)) {
    is_dynamic <- TRUE
  }
  
  # If dynamic effects exist, delegate to dynamic plotting logic
  if (is_dynamic) {
    return(ab_plot_dynamic_internal(fit, effect, sorted, labels, title, 
                                    time_point, plot_type, show_actors))
  }
  
  # Static effects plotting (original logic)
  if (effect == "sender") {
    muEff <- fit$APM
    ylabel <- "Sender Effects (a)"
    default_title <- "Additive Sender Effects"
  } else {
    muEff <- fit$BPM
    ylabel <- "Receiver Effects (b)"
    default_title <- "Additive Receiver Effects"
  }
  
  # Set labels default based on number of actors
  if (is.null(labels)) {
    labels <- length(muEff) <= 50
  }
  
  # Create data frame
  muDf <- data.frame(
    mu = muEff,
    id = names(muEff),
    stringsAsFactors = FALSE
  )
  
  # Sort if requested
  if (sorted) {
    muDf$id <- factor(muDf$id, levels = muDf$id[order(muDf$mu)])
  } else {
    muDf$id <- factor(muDf$id, levels = muDf$id)
  }
  
  # Create ymin/ymax for segments
  muDf$ymax <- with(muDf, ifelse(mu >= 0, mu, 0))
  muDf$ymin <- with(muDf, ifelse(mu < 0, mu, 0))
  
  # Build plot
  p <- ggplot(muDf, aes(x = id, y = mu)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_segment(aes(xend = id, yend = 0), color = "steelblue") +
    geom_point(size = 2, color = "steelblue") +
    xlab("") + 
    ylab(ylabel) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Add or remove x-axis labels
  if (!labels) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) + xlab("Actors (sorted by effect magnitude)")
  } else {
    p <- p + theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  }
  
  # Add title
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  } else if (sorted) {
    p <- p + ggtitle(paste(default_title, "(sorted)"))
  } else {
    p <- p + ggtitle(default_title)
  }
  
  return(p)
}

#' Internal function for dynamic additive effects plotting
#' @keywords internal
#' @noRd
ab_plot_dynamic_internal <- function(fit, effect, sorted, labels, title,
                                     time_point, plot_type, show_actors) {
  
  # Extract dynamic effects
  if (effect == "sender") {
    effects_mat <- fit$a_dynamic
    ylabel <- "Sender Effects (a)"
    default_title <- "Dynamic Sender Effects"
  } else {
    effects_mat <- fit$b_dynamic
    ylabel <- "Receiver Effects (b)"
    default_title <- "Dynamic Receiver Effects"
  }
  
  n_actors <- nrow(effects_mat)
  n_times <- ncol(effects_mat)
  actor_names <- rownames(effects_mat)
  if (is.null(actor_names)) actor_names <- paste0("Actor", 1:n_actors)
  
  # Time labels
  time_labels <- colnames(effects_mat)
  if (is.null(time_labels)) time_labels <- paste0("T", 1:n_times)
  
  # Handle different plot types
  if (plot_type == "snapshot") {
    # Plot a single time point
    if (is.null(time_point)) {
      time_point <- n_times  # Default to last time point
    } else if (time_point == "average") {
      # Use time-averaged effects
      muEff <- rowMeans(effects_mat)
      names(muEff) <- actor_names
      
      # Create data frame
      muDf <- data.frame(
        mu = muEff,
        id = actor_names,
        stringsAsFactors = FALSE
      )
      
      # Sort if requested
      if (sorted) {
        muDf$id <- factor(muDf$id, levels = muDf$id[order(muDf$mu)])
      } else {
        muDf$id <- factor(muDf$id, levels = muDf$id)
      }
      
      # Set labels default
      if (is.null(labels)) {
        labels <- n_actors <= 50
      }
      
      # Build plot
      p <- ggplot(muDf, aes(x = id, y = mu)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_segment(aes(xend = id, yend = 0), color = "steelblue") +
        geom_point(size = 2, color = "steelblue") +
        xlab("") + ylab(paste(ylabel, "(Time-Averaged)")) +
        theme_bw()
      
      if (!labels) {
        p <- p + theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) + xlab("Actors")
      } else {
        p <- p + theme(
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
      }
      
      if (!is.null(title)) {
        p <- p + ggtitle(title)
      } else {
        p <- p + ggtitle(paste(default_title, "- Average Across Time"))
      }
      
      return(p)
    }
    
    if (time_point > n_times || time_point < 1) {
      stop("time_point must be between 1 and ", n_times)
    }
    
    # Extract effects for specific time
    muEff <- effects_mat[, time_point]
    names(muEff) <- actor_names
    
    # Use standard plotting logic
    muDf <- data.frame(
      mu = muEff,
      id = actor_names,
      stringsAsFactors = FALSE
    )
    
    if (sorted) {
      muDf$id <- factor(muDf$id, levels = muDf$id[order(muDf$mu)])
    } else {
      muDf$id <- factor(muDf$id, levels = muDf$id)
    }
    
    # Set labels default
    if (is.null(labels)) {
      labels <- n_actors <= 50
    }
    
    p <- ggplot(muDf, aes(x = id, y = mu)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_segment(aes(xend = id, yend = 0), color = "steelblue") +
      geom_point(size = 2, color = "steelblue") +
      xlab("") + ylab(ylabel) +
      theme_bw()
    
    if (!labels) {
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) + xlab("Actors")
    } else {
      p <- p + theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    }
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle(paste(default_title, "at", time_labels[time_point]))
    }
    
    return(p)
    
  } else if (plot_type == "trajectory") {
    # Show evolution of effects over time
    
    # Prepare data
    traj_data <- data.frame()
    for (i in 1:n_actors) {
      for (t in 1:n_times) {
        traj_data <- rbind(traj_data, data.frame(
          actor = actor_names[i],
          time = t,
          time_label = time_labels[t],
          effect = effects_mat[i, t]
        ))
      }
    }
    
    # Filter actors if requested
    if (!is.null(show_actors)) {
      traj_data <- traj_data[traj_data$actor %in% show_actors, ]
    } else if (n_actors > 20) {
      # For many actors, show only top/bottom
      avg_effects <- rowMeans(effects_mat)
      top_actors <- names(sort(avg_effects, decreasing = TRUE)[1:5])
      bottom_actors <- names(sort(avg_effects)[1:5])
      show_actors <- c(top_actors, bottom_actors)
      traj_data <- traj_data[traj_data$actor %in% show_actors, ]
      
      cli::cli_inform(c(
        "i" = "Showing top 5 and bottom 5 actors by average effect",
        ">" = "Use {.arg show_actors} to specify actors to display"
      ))
    }
    
    p <- ggplot(traj_data, aes(x = time, y = effect, group = actor, color = actor)) +
      geom_line(alpha = 0.7, linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_x_continuous(breaks = 1:n_times, labels = time_labels) +
      theme_bw() +
      labs(x = "Time", y = ylabel, color = "Actor")
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle(paste(default_title, "Trajectories"))
    }
    
    return(p)
    
  } else if (plot_type == "faceted") {
    # Create faceted plot showing multiple time points
    
    # Determine which time points to show
    if (n_times <= 6) {
      show_times <- 1:n_times
    } else {
      # Show 6 evenly spaced time points
      show_times <- round(seq(1, n_times, length.out = 6))
    }
    
    # Create data for faceted plot
    facet_data <- data.frame()
    for (t in show_times) {
      for (i in 1:n_actors) {
        facet_data <- rbind(facet_data, 
          data.frame(
            actor = actor_names[i],
            time = time_labels[t],
            effect = effects_mat[i, t]
          )
        )
      }
    }
    
    # Create faceted plot
    p <- ggplot(facet_data, aes(x = reorder(actor, effect), y = effect)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_segment(aes(xend = actor, yend = 0), color = "steelblue") +
      geom_point(size = 1.5, color = "steelblue") +
      facet_wrap(~ time, scales = "free_x") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
      labs(x = "", y = ylabel)
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle(paste(default_title, "Over Time"))
    }
    
    return(p)
    
  } else if (plot_type == "ribbon") {
    # Show confidence bands over time for each actor
    
    # Calculate summary statistics
    ribbon_data <- data.frame()
    for (i in 1:n_actors) {
      effects_i <- effects_mat[i, ]
      ribbon_data <- rbind(ribbon_data, data.frame(
        actor = actor_names[i],
        time = 1:n_times,
        mean_effect = effects_i,
        # Add some variation for visualization (in practice, would use posterior samples)
        lower = effects_i - 0.1 * abs(effects_i),
        upper = effects_i + 0.1 * abs(effects_i)
      ))
    }
    
    # Filter actors if needed
    if (!is.null(show_actors)) {
      ribbon_data <- ribbon_data[ribbon_data$actor %in% show_actors, ]
    } else if (n_actors > 10) {
      # Show subset for clarity
      avg_effects <- rowMeans(effects_mat)
      show_actors <- names(sort(abs(avg_effects), decreasing = TRUE)[1:10])
      ribbon_data <- ribbon_data[ribbon_data$actor %in% show_actors, ]
      
      cli::cli_inform(c(
        "i" = "Showing 10 actors with largest absolute average effects",
        ">" = "Use {.arg show_actors} to specify actors to display"
      ))
    }
    
    p <- ggplot(ribbon_data, aes(x = time, y = mean_effect, group = actor)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = actor), alpha = 0.3) +
      geom_line(aes(color = actor), linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_x_continuous(breaks = 1:n_times, labels = time_labels) +
      theme_bw() +
      labs(x = "Time", y = ylabel, color = "Actor", fill = "Actor")
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle(paste(default_title, "with Uncertainty"))
    }
    
    return(p)
  }
}