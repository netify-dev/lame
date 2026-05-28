#' TIES sanctions data for vignettes
#'
#' Longitudinal directed binary network of international economic sanctions
#' among 35 countries across four years (1993, 1994, 1995, 2000), derived from
#' the Threat and Imposition of Sanctions (TIES) dataset (Morgan et al. 2014).
#' Entry \eqn{y_{ij,t} = 1} means country \eqn{i} imposed sanctions on country
#' \eqn{j} in year \eqn{t}. Network density is approximately 2\%.
#'
#' When loaded via \code{data(vignette_data)}, the following objects are placed
#' in the calling environment: \code{Y}, \code{Xdyad}, \code{Xrow}, \code{Xcol}.
#'
#' @name vignette_data
#' @docType data
#' @usage data(vignette_data)
#' @format Four objects:
#' \describe{
#'   \item{Y}{List of 4 binary adjacency matrices (35 x 35), one per year.}
#'   \item{Xdyad}{List of 4 arrays (35 x 35 x 2) of dyadic covariates:
#'     \code{distance} (geographic) and \code{shared_igos} (shared IGO memberships).}
#'   \item{Xrow}{List of 4 matrices (35 x 2) of sender covariates:
#'     \code{log_gdp} and \code{log_pop}.}
#'   \item{Xcol}{List of 4 matrices (35 x 2) of receiver covariates:
#'     \code{log_gdp} and \code{log_pop}.}
#' }
#' @references Morgan, T. Clifton, Bapat, N., & Kobayashi, Y. (2014). Threat
#'   and Imposition of Economic Sanctions 1945--2005. \emph{Conflict Management
#'   and Peace Science}, 31(5), 541--558.
#' @keywords datasets
#' @examples
#' data(vignette_data)
#' cat("Countries:", nrow(Y[[1]]), "\n")
#' cat("Time periods:", length(Y), "\n")
NULL

#' Synthetic longitudinal binary relational data, list-form (latent-scale)
#'
#' The list-form sibling of \code{\link{YX_bin_long}}: same synthetic
#' 4-period 50-actor panel reshaped from arrays to lists. The shipped
#' \code{Y} entries hold the \emph{latent probit-scale predictor}
#' (range roughly -20 to 17), not 0/1 ties. To recover the binary tie
#' indicator the dataset name implies, threshold at zero, e.g.
#' \code{Y_bin <- lapply(YX_bin_list$Y, function(Yt) 1 * (Yt > 0))}.
#' If passed unthresholded to \code{lame(..., family = "binary")} the
#' fit will warn and silently apply the same \code{Y > 0} threshold.
#'
#' @name YX_bin_list
#' @docType data
#' @format A list with two elements:
#'   \describe{
#'     \item{Y}{List of 4 numeric matrices, each \code{[50, 50]}, with
#'       \code{NA} on the diagonal. Latent probit-scale predictor;
#'       threshold at 0 for the binary indicator.}
#'     \item{X}{List of 4 numeric arrays, each \code{[50, 50, 3]} of
#'       dyadic covariates.}
#'   }
#' @keywords datasets
NULL

#' Column covariates
#'
#' Column covariates loaded by \code{data(vignette_data)}.
#' See \code{\link{vignette_data}} for full details.
#'
#' @name Xcol
#' @docType data
#' @format A list of matrices of column covariates
#' @seealso \code{\link{vignette_data}}
#' @keywords datasets
NULL

#' Dyadic covariates
#'
#' Dyadic covariates loaded by \code{data(vignette_data)}.
#' See \code{\link{vignette_data}} for full details.
#'
#' @name Xdyad
#' @docType data
#' @format A list of arrays of dyadic covariates
#' @seealso \code{\link{vignette_data}}
#' @keywords datasets
NULL

#' Row covariates
#'
#' Row covariates loaded by \code{data(vignette_data)}.
#' See \code{\link{vignette_data}} for full details.
#'
#' @name Xrow
#' @docType data
#' @format A list of matrices of row covariates
#' @seealso \code{\link{vignette_data}}
#' @keywords datasets
NULL

#' Relational matrix
#'
#' Relational matrix loaded by \code{data(vignette_data)}.
#' See \code{\link{vignette_data}} for full details.
#'
#' @name Y
#' @docType data
#' @format A list of square binary adjacency matrices
#' @seealso \code{\link{vignette_data}}
#' @keywords datasets
NULL