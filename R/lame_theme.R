# internal ggplot theme used across the package. follows the project
# style rules: no panel border, no axis ticks, gridlines kept, legend
# on top, facet strips in black with left-aligned white text.

#' @noRd
.lame_theme <- function() {
	ggplot2::theme_bw() +
		ggplot2::theme(
			panel.border     = ggplot2::element_blank(),
			axis.ticks       = ggplot2::element_blank(),
			legend.position  = "top",
			strip.background = ggplot2::element_rect(fill = "black", color = "black"),
			strip.text.x     = ggplot2::element_text(color = "white", hjust = 0),
			strip.text.y     = ggplot2::element_text(color = "white", hjust = 0)
		)
}
