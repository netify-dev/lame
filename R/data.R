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

#' Binary relational data list
#' 
#' A list containing binary relational data in list format with associated covariates
#' 
#' @name YX_bin_list
#' @docType data
#' @format A list containing Y and X matrices for binary relational data
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