#' Visualize sender and receiver random effects
#' 
#' Creates a visualization of the additive sender (row) and receiver (column) 
#' random effects from an AME or LAME model. These effects capture heterogeneity 
#' in actors' propensities to send and receive ties that is not explained by 
#' covariates.
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
#' The plot displays these effects as a dot plot with vertical lines extending 
#' from zero to each effect estimate, making it easy to identify:
#' - Which actors are most active/inactive as senders
#' - Which actors are most popular/unpopular as receivers
#' - The overall distribution and range of effects
#' 
#' @param fit An object of class "ame" or "lame" from fitting an AME model
#' @param effect Character string specifying which effect to plot: 
#'        "sender" (default) or "receiver"
#' @param sorted Logical; if TRUE (default), actors are sorted by effect magnitude
#' @param labels Logical; if TRUE, actor labels are shown on x-axis (default TRUE
#'        for n <= 50 actors)
#' @param title Optional title for the plot
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
ab_plot <- function(fit, effect = c("sender", "receiver"), 
                    sorted = TRUE, labels = NULL, title = NULL) {
  
  # Check class
  if (!inherits(fit, c("ame", "lame"))) {
    stop("fit must be an object of class 'ame' or 'lame'")
  }
  
  # Match effect argument
  effect <- match.arg(effect)
  
  # Extract appropriate effects
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