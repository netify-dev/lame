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
uv_plot <- function(fit = NULL, Y = NULL, U = NULL, V = NULL,
                    row.names = NULL, col.names = NULL,
                    layout = c("circle", "biplot"),
                    vscale = 0.8,
                    show.edges = FALSE,
                    edge.alpha = 0.3,
                    node.size = 3,
                    label.nodes = TRUE,
                    label.size = 3,
                    colors = NULL,
                    title = NULL) {
  
  # Match layout argument
  layout <- match.arg(layout)
  
  # Extract components from fit if provided
  if (!is.null(fit)) {
    if (!inherits(fit, c("ame", "lame"))) {
      stop("fit must be an object of class 'ame' or 'lame'")
    }
    if (is.null(U)) U <- fit$U
    if (is.null(V)) V <- fit$V
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
  
  # Create node data frame
  if (layout == "circle") {
    # Circular layout
    n_total <- n_row + n_col
    angles <- seq(0, 2*pi, length.out = n_total + 1)[1:n_total]
    
    # Base positions on circle
    x_base <- cos(angles)
    y_base <- sin(angles)
    
    # Combine U and V with offsets from circle
    node_data <- data.frame(
      x = c(x_base[1:n_row] + U[,1]*0.2, 
            x_base[(n_row+1):n_total] + V[,1]*0.2),
      y = c(y_base[1:n_row] + U[,2]*0.2, 
            y_base[(n_row+1):n_total] + V[,2]*0.2),
      name = c(row.names, col.names),
      type = factor(c(rep("Sender", n_row), rep("Receiver", n_col))),
      label = "",  # Initialize label column
      stringsAsFactors = FALSE
    )
  } else {
    # Biplot layout
    node_data <- data.frame(
      x = c(U[,1], V[,1]),
      y = c(U[,2], V[,2]),
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
          edge_data <- rbind(edge_data, data.frame(
            x = U[i,1],
            y = U[i,2],
            xend = V[j,1],
            yend = V[j,2],
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
    p <- p + geom_point(aes(color = type, shape = type, size = size))
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