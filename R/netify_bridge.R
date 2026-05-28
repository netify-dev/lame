#' Using `lame` with `netify` data
#'
#' For data prepared with the \code{netify} package, load \code{netify}
#' and convert to \code{lame}-ready input with
#' \code{netify::to_amen(netlet, lame = TRUE)}. Pass the result's
#' \code{Y}, \code{Xdyad}, \code{Xrow}, \code{Xcol} components straight
#' to \code{lame()} or \code{ame()}.
#'
#' \preformatted{
#' library(netify); library(lame)
#' amen_in <- netify::to_amen(netlet, lame = TRUE)
#' fit <- lame(
#'     Y      = amen_in$Y,
#'     Xdyad  = amen_in$Xdyad,
#'     family = "binary", R = 2,
#'     nscan  = 1000, burn = 250, odens = 5
#' )
#' }
#'
#' For longitudinal models, set \code{lame = TRUE} on the converter so
#' the dyadic covariates carry the period axis.
#'
#' @name netify_to_lame
#' @keywords internal
NULL
