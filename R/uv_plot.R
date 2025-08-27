#' Visualize multiplicative effects (latent factors) from AME models
#' 
#' Creates a two-dimensional visualization of the multiplicative effects (U and V)
#' from an AME or LAME model. These latent factors capture network structure 
#' beyond what is explained by covariates and additive effects, including 
#' clustering, transitivity, and other higher-order dependencies.
#' 
#' @details
#' The multiplicative effects in AME models provide a low-rank representation
#' of network structure through latent factors:
#' \describe{
#'   \item{U matrix}{Sender-specific latent positions (row factors)}
#'   \item{V matrix}{Receiver-specific latent positions (column factors)}
#'   \item{UV' product}{Captures dyad-specific effects beyond additive terms}
#' }
#' 
#' The visualization can show:
#' \describe{
#'   \item{Circular layout}{Default layout placing nodes on a circle with
#'         latent positions shown as deviations}
#'   \item{Biplot layout}{Shows U and V positions directly in latent space}
#'   \item{Network overlay}{Optional display of actual network ties}
#' }
#' 
#' Interpretation:
#' - Nodes close together in latent space tend to have similar connection patterns
#' - The distance between sender position (U) and receiver position (V) relates
#'   to the likelihood of a tie
#' - Clustering in the latent space indicates community structure
#' 
#' @param fit An object of class "ame" or "lame" containing multiplicative effects,
#'        or a network matrix Y if U and V are provided separately
#' @param Y Network matrix (only needed if fit is not provided)
#' @param U Matrix of sender latent factors (extracted from fit if not provided)
#' @param V Matrix of receiver latent factors (extracted from fit if not provided)
#' @param row.names Names for row nodes (defaults to rownames of Y or U)
#' @param col.names Names for column nodes (defaults to colnames of Y or V)
#' @param layout Character string specifying layout: "circle" (default) or "biplot"
#' @param vscale Scaling factor for V positions relative to U (default 0.8)
#' @param show.edges Logical; if TRUE, show network edges (default FALSE)
#' @param edge.alpha Transparency for edges (default 0.3)
#' @param node.size Size of nodes, or "degree" to scale by degree (default 3)
#' @param label.nodes Logical; if TRUE, show node labels (default TRUE)
#' @param label.size Size of node labels (default 3)
#' @param colors Optional vector of colors for nodes (e.g., for communities)
#' @param title Optional title for the plot
#' @param time_point For dynamic UV, which time point to plot (default: last).
#'                   Can be numeric index or "average" for time-averaged positions
#' @param plot_type For dynamic UV: "snapshot" (single time), "trajectory" (evolution),
#'                  "faceted" (grid of time points). For static UV, this is ignored.
#' @param show_arrows For trajectory plots, whether to show directional arrows
#' @return A ggplot2 object that can be further customized
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \dontrun{
#' # Fit an AME model with multiplicative effects
#' fit <- ame(Y, X, R = 2)
#' 
#' # Basic visualization
#' uv_plot(fit)
#' 
#' # Show network edges
#' uv_plot(fit, show.edges = TRUE)
#' 
#' # Use biplot layout
#' uv_plot(fit, layout = "biplot")
#' 
#' # Color by known communities
#' uv_plot(fit, colors = community_membership)
#' }
#' @export
#' @import ggplot2
#' @import ggrepel
#' @import cli
uv_plot <- function(
  fit = NULL, Y = NULL, U = NULL, V = NULL,
  row.names = NULL, col.names = NULL,
  layout = c("circle", "biplot"),
  vscale = 0.8,
  show.edges = FALSE,
  edge.alpha = 0.3,
  node.size = 3,
  label.nodes = TRUE,
  label.size = 3,
  colors = NULL,
  title = NULL,
  time_point = NULL,
  plot_type = c("snapshot", "trajectory", "faceted"),
  show_arrows = TRUE
  ){
  
  # Match arguments
  layout <- match.arg(layout)
  plot_type <- match.arg(plot_type)
  
  # Extract components from fit if provided
  if (!is.null(fit)) {
    if (!inherits(fit, c("ame", "lame"))) {
      stop("fit must be an object of class 'ame' or 'lame'")
    }
    if (is.null(U)) U <- fit$U
    if (is.null(V)) V <- fit$V
    
    # Check for dynamic UV (3D arrays)
    is_dynamic <- FALSE
    if (!is.null(U) && length(dim(U)) == 3) {
      is_dynamic <- TRUE
      # Delegate to dynamic plotting if appropriate plot_type is requested
      if (plot_type != "snapshot" || (!is.null(time_point) && time_point != "average")) {
        return(uv_plot_dynamic_internal(fit, U, V, time_point, plot_type, 
                                       show_arrows, title, label.nodes, 
                                       label.size, colors))
      }
      # For snapshot at specific time or average
      n_times <- dim(U)[3]
      if (is.null(time_point)) {
        time_point <- n_times  # Default to last time point
        cli::cli_inform(c(
          "i" = "Detected dynamic UV effects",
          "*" = "Using last time point for visualization"
        ))
      }
      if (time_point == "average") {
        U <- apply(U, c(1, 2), mean)
        if (!is.null(V) && length(dim(V)) == 3) V <- apply(V, c(1, 2), mean)
        cli::cli_inform("Using time-averaged positions")
      } else {
        U <- U[, , time_point]
        if (!is.null(V) && length(dim(V)) == 3) V <- V[, , time_point]
      }
    }
    if (is.null(Y)) {
      if (!is.null(fit$Y)) {
        Y <- fit$Y
      } else if (!is.null(fit$YPM)) {
        Y <- fit$YPM
      }
      # For LAME, Y might be a list - use the mean across time
      if (is.list(Y)) {
        # Convert list to array and average
        if (length(Y) > 0) {
          Y_avg <- Y[[1]]
          if (length(Y) > 1) {
            for (t in 2:length(Y)) {
              Y_avg <- Y_avg + Y[[t]]
            }
            Y_avg <- Y_avg / length(Y)
          }
          Y <- Y_avg
        }
      }
    }
  }
  
  # Check that we have U and V
  if (is.null(U) || is.null(V)) {
    if (is.null(Y)) {
      stop("Must provide either fit object with U and V, or Y matrix to compute factors")
    }
    # Compute factors from Y using SVD
    Y_centered <- Y - mean(Y, na.rm = TRUE)
    Y_centered[is.na(Y_centered)] <- 0
    svd_Y <- svd(Y_centered)
    U <- svd_Y$u[, 1:2]
    V <- svd_Y$v[, 1:2] * vscale
  }
  
  # Ensure we have 2D factors
  if (ncol(U) > 2) U <- U[, 1:2]
  if (ncol(V) > 2) V <- V[, 1:2]
  if (ncol(U) == 1) U <- cbind(U, 0)
  if (ncol(V) == 1) V <- cbind(V, 0)
  
  # Get node names
  n_row <- nrow(U)
  n_col <- nrow(V)
  if (is.null(row.names)) {
    row.names <- rownames(U)
    if (is.null(row.names)) row.names <- paste0("S", 1:n_row)
  }
  if (is.null(col.names)) {
    col.names <- rownames(V)
    if (is.null(col.names)) col.names <- paste0("R", 1:n_col)
  }
  
  # Standardize U and V to have comparable scales
  # This prevents one dimension from dominating
  U_scaled <- U
  V_scaled <- V
  
  # Scale each dimension to unit variance if variance > 0
  for(d in 1:ncol(U)) {
    u_sd <- sd(U[,d])
    v_sd <- sd(V[,d])
    combined_sd <- sqrt((u_sd^2 * n_row + v_sd^2 * n_col) / (n_row + n_col))
    
    if(combined_sd > 1e-10) {
      U_scaled[,d] <- U[,d] / combined_sd
      V_scaled[,d] <- V[,d] / combined_sd
    }
  }
  
  # Create node data frame
  if (layout == "circle") {
    # Circular layout - order by first principal component
    n_total <- n_row + n_col
    
    # Use first dimension to determine angular position
    combined_scores <- c(U_scaled[,1], V_scaled[,1])
    # Map scores to angles
    score_range <- range(combined_scores)
    if(diff(score_range) > 0) {
      angles <- 2 * pi * (combined_scores - score_range[1]) / diff(score_range)
    } else {
      angles <- seq(0, 2*pi, length.out = n_total + 1)[1:n_total]
    }
    
    # Use second dimension for radial displacement
    radii <- 1 + 0.3 * c(U_scaled[,2], V_scaled[,2])
    
    # Compute positions
    node_data <- data.frame(
      x = radii * cos(angles),
      y = radii * sin(angles),
      name = c(row.names, col.names),
      type = factor(c(rep("Sender", n_row), rep("Receiver", n_col))),
      label = "",  # Initialize label column
      stringsAsFactors = FALSE
    )
  } else {
    # Biplot layout - show scaled positions directly
    node_data <- data.frame(
      x = c(U_scaled[,1], V_scaled[,1] * vscale),
      y = c(U_scaled[,2], V_scaled[,2] * vscale),
      name = c(row.names, col.names),
      type = factor(c(rep("Sender", n_row), rep("Receiver", n_col))),
      label = "",  # Initialize label column
      stringsAsFactors = FALSE
    )
  }
  
  # Add colors if provided
  if (!is.null(colors)) {
    if (length(colors) != nrow(node_data)) {
      warning("Length of colors doesn't match number of nodes, ignoring colors")
      colors <- NULL
    } else {
      node_data$color <- colors
    }
  }
  
  # Calculate node sizes if based on degree
  if (is.character(node.size) && node.size == "degree" && !is.null(Y)) {
    out_degree <- rowSums(Y > 0, na.rm = TRUE)
    in_degree <- colSums(Y > 0, na.rm = TRUE)
    node_data$size <- c(out_degree, in_degree)
    node_data$size <- 1 + 4 * (node_data$size - min(node_data$size)) / 
                     (max(node_data$size) - min(node_data$size) + 1)
  } else {
    node_data$size <- node.size
  }
  
  # Start building plot
  p <- ggplot(node_data, aes(x = x, y = y))
  
  # Add edges if requested
  if (show.edges && !is.null(Y)) {
    edge_data <- data.frame()
    for (i in 1:n_row) {
      for (j in 1:n_col) {
        if (!is.na(Y[i,j]) && Y[i,j] > 0) {
          # Use the actual node positions from node_data
          edge_data <- rbind(edge_data, data.frame(
            x = node_data$x[i],  # Sender position from node_data
            y = node_data$y[i],
            xend = node_data$x[n_row + j],  # Receiver position from node_data
            yend = node_data$y[n_row + j],
            weight = Y[i,j]
          ))
        }
      }
    }
    if (nrow(edge_data) > 0) {
      p <- p + geom_segment(data = edge_data,
                           aes(x = x, y = y, xend = xend, yend = yend),
                           alpha = edge.alpha, color = "gray50")
    }
  }
  
  # Add nodes
  if (!is.null(colors)) {
    p <- p + geom_point(aes(color = color, shape = type, size = size))
  } else {
    # Use black points with different shapes for types
    p <- p + geom_point(aes(shape = type, size = size), color = "black")
  }
  
  # Add labels
  if (label.nodes) {
    if (nrow(node_data) <= 50) {
      # Use ggrepel for small networks
      node_data$label <- node_data$name
    } else {
      # Just show labels for high-degree nodes in large networks
      if (!is.null(Y)) {
        degree <- c(rowSums(Y > 0, na.rm = TRUE), colSums(Y > 0, na.rm = TRUE))
        top_nodes <- degree >= quantile(degree, 0.9)
        node_data$label <- ifelse(top_nodes, node_data$name, "")
      } else {
        # If no Y matrix, show no labels for large networks
        node_data$label <- ""
      }
    }
    
    # Add the labels using the label column
    if (any(node_data$label != "")) {
      p <- p + geom_text_repel(
        data = node_data,
        aes(x = x, y = y, label = label, color = type),
        size = label.size,
        max.overlaps = 20
      )
    }
  }
  
  # Customize appearance
  p <- p + 
    scale_size_continuous(range = c(2, 6), guide = "none") +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom"
    ) +
    coord_fixed()
  
  # Add title
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  } else if (layout == "circle") {
    p <- p + ggtitle("Multiplicative Effects (Circular Layout)")
  } else {
    p <- p + ggtitle("Multiplicative Effects (Biplot)")
  }
  
  return(p)
}

#' Internal function for dynamic UV plotting
#' @keywords internal
#' @noRd
uv_plot_dynamic_internal <- function(fit, U, V, time_point, plot_type, 
                                     show_arrows, title, label.nodes,
                                     label.size, colors) {
  
  # Get dimensions
  n <- dim(U)[1]
  R <- dim(U)[2]
  n_times <- dim(U)[3]
  
  # Get node names
  node_names <- rownames(U)
  if (is.null(node_names)) node_names <- paste0("Node", 1:n)
  
  # Time labels
  time_labels <- dimnames(U)[[3]]
  if (is.null(time_labels)) time_labels <- paste0("T", 1:n_times)
  
  # Ensure we're using first 2 dimensions for visualization
  if (R > 2) {
    cli::cli_inform("Using first 2 dimensions of {R}-dimensional latent space")
    U <- U[, 1:2, , drop = FALSE]
    if (!is.null(V)) V <- V[, 1:2, , drop = FALSE]
  }
  
  if (plot_type == "snapshot") {
    # Single time point - already handled in main function
    # This shouldn't be reached but included for completeness
    if (is.null(time_point)) time_point <- n_times
    
    U_t <- U[, , time_point]
    V_t <- if (!is.null(V)) V[, , time_point] else U_t
    
    # Create data frame
    positions <- data.frame(
      x = U_t[, 1],
      y = U_t[, 2],
      node = node_names,
      type = "U"
    )
    
    if (!is.null(V)) {
      positions <- rbind(positions,
        data.frame(
          x = V_t[, 1],
          y = V_t[, 2],
          node = node_names,
          type = "V"
        )
      )
    }
    
    p <- ggplot(positions, aes(x = x, y = y)) +
      geom_point(aes(color = if (!is.null(colors)) colors else type), size = 3) +
      theme_minimal() +
      coord_fixed() +
      labs(x = "Dimension 1", y = "Dimension 2")
    
    if (label.nodes) {
      p <- p + geom_text_repel(aes(label = node), size = label.size)
    }
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle(paste("UV Positions at", time_labels[time_point]))
    }
    
    return(p)
    
  } else if (plot_type == "trajectory") {
    # Show trajectories over time
    
    # Prepare trajectory data
    traj_data <- data.frame()
    for (i in 1:n) {
      for (t in 1:n_times) {
        traj_data <- rbind(traj_data, data.frame(
          x = U[i, 1, t],
          y = U[i, 2, t],
          time = t,
          time_label = time_labels[t],
          node = node_names[i],
          type = "U"
        ))
        
        if (!is.null(V)) {
          traj_data <- rbind(traj_data, data.frame(
            x = V[i, 1, t],
            y = V[i, 2, t],
            time = t,
            time_label = time_labels[t],
            node = node_names[i],
            type = "V"
          ))
        }
      }
    }
    
    # Create plot
    p <- ggplot(traj_data, aes(x = x, y = y, group = interaction(node, type))) +
      geom_path(aes(color = node, linetype = type), alpha = 0.5) +
      theme_minimal() +
      coord_fixed() +
      labs(x = "Dimension 1", y = "Dimension 2")
    
    if (show_arrows) {
      # Add arrows to show direction
      arrow_data <- traj_data[traj_data$time == n_times, ]
      prev_data <- traj_data[traj_data$time == n_times - 1, ]
      arrow_data$x_prev <- prev_data$x
      arrow_data$y_prev <- prev_data$y
      
      p <- p + geom_segment(data = arrow_data,
                           aes(x = x_prev, y = y_prev, xend = x, yend = y,
                               color = node),
                           arrow = arrow(length = unit(0.2, "cm")),
                           size = 1)
    }
    
    # Add start and end points
    start_data <- traj_data[traj_data$time == 1, ]
    end_data <- traj_data[traj_data$time == n_times, ]
    
    p <- p + 
      geom_point(data = start_data, aes(color = node, shape = type), 
                size = 3, alpha = 0.5) +
      geom_point(data = end_data, aes(color = node, shape = type), 
                size = 4)
    
    if (label.nodes && n <= 20) {
      p <- p + geom_text_repel(data = end_data,
                               aes(label = node, color = node),
                               size = label.size)
    }
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle("UV Trajectories Over Time")
    }
    
    return(p)
    
  } else if (plot_type == "faceted") {
    # Faceted plot showing multiple time points
    
    # Select time points to show
    if (n_times <= 6) {
      show_times <- 1:n_times
    } else {
      show_times <- round(seq(1, n_times, length.out = 6))
    }
    
    # Prepare data
    facet_data <- data.frame()
    for (t in show_times) {
      for (i in 1:n) {
        facet_data <- rbind(facet_data, data.frame(
          x = U[i, 1, t],
          y = U[i, 2, t],
          node = node_names[i],
          time = time_labels[t],
          type = "U"
        ))
        
        if (!is.null(V)) {
          facet_data <- rbind(facet_data, data.frame(
            x = V[i, 1, t],
            y = V[i, 2, t],
            node = node_names[i],
            time = time_labels[t],
            type = "V"
          ))
        }
      }
    }
    
    # Create plot
    p <- ggplot(facet_data, aes(x = x, y = y)) +
      geom_point(aes(color = if (!is.null(colors)) colors else type, 
                    shape = type), size = 2) +
      facet_wrap(~ time, scales = "free") +
      theme_minimal() +
      coord_fixed() +
      labs(x = "Dimension 1", y = "Dimension 2")
    
    if (label.nodes && n <= 10) {
      p <- p + geom_text_repel(aes(label = node), size = label.size * 0.8)
    }
    
    if (!is.null(title)) {
      p <- p + ggtitle(title)
    } else {
      p <- p + ggtitle("UV Positions Over Time")
    }
    
    return(p)
  }
}