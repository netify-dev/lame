#' Validate input data format for lame function
#'
#' Internal validation function that checks the format and consistency of input
#' data for longitudinal AME models. Ensures that network data and covariates
#' are properly formatted as lists with consistent dimensions across time periods.
#'
#' @usage check_format(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL)
#' @param Y a list of T network matrices, where T is the number of time periods.
#' Each element should be an n x m matrix representing the network at time t.
#' @param Xdyad an optional list of T dyadic covariates. Each element can be
#' either an n x m matrix (single covariate) or an n x m x p array (p covariates).
#' Must have the same length as Y if provided.
#' @param Xrow an optional list of T matrices of row/sender covariates.
#' Each element should be an n x pr matrix where pr is the number of row covariates.
#' Must have the same length as Y if provided.
#' @param Xcol an optional list of T matrices of column/receiver covariates.
#' Each element should be an m x pc matrix where pc is the number of column covariates.
#' Must have the same length as Y if provided.
#'
#' @details
#' This function performs comprehensive validation of input data for longitudinal
#' network analysis, including:
#' \itemize{
#'   \item Verifying Y is a non-empty list of matrices
#'   \item Checking that all covariates (if provided) are lists of appropriate length
#'   \item Validating data types for all elements
#'   \item Warning about dimension inconsistencies across time periods
#' }
#'
#' The function uses informative error messages via the cli package to help
#' users identify and correct data formatting issues.
#'
#' @return Invisible TRUE if all checks pass. Throws an error with an informative
#' message if any validation fails.
#'
#' @note This is an internal function primarily used by \code{lame()} but exported
#' for advanced users who want to validate their data before model fitting.
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#'
#' @export check_format

check_format <- function(Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL) {
  # Check Y is a list
  if (!is.list(Y)) {
    cli::cli_abort(c(
      "Y must be a list of matrices.",
      "i" = "Received object of class {class(Y)[1]}."
    ))
  }

  # Check Y is not empty
  if (length(Y) == 0) {
    cli::cli_abort("Y must contain at least one matrix.")
  }

  # Check each element of Y is a matrix
  for (t in seq_along(Y)) {
    if (!is.matrix(Y[[t]]) && !is.array(Y[[t]])) {
      cli::cli_abort(c(
        "All elements of Y must be matrices.",
        "i" = "Y[[{t}]] is of class {class(Y[[t]])[1]}."
      ))
    }
  }

  # Get the number of time periods
  N <- length(Y)

  # Check Xdyad if provided
  if (!is.null(Xdyad)) {
    if (!is.list(Xdyad)) {
      cli::cli_abort(c(
        "Xdyad must be NULL or a list of matrices/arrays.",
        "i" = "Received object of class {class(Xdyad)[1]}."
      ))
    }
    if (length(Xdyad) != N) {
      cli::cli_abort(c(
        "Xdyad must have the same length as Y.",
        "i" = "Y has {N} time periods, but Xdyad has {length(Xdyad)}."
      ))
    }
    # Check each element is a matrix or array
    for (t in seq_along(Xdyad)) {
      if (!is.null(Xdyad[[t]]) && !is.matrix(Xdyad[[t]]) && !is.array(Xdyad[[t]])) {
        cli::cli_abort(c(
          "All non-NULL elements of Xdyad must be matrices or arrays.",
          "i" = "Xdyad[[{t}]] is of class {class(Xdyad[[t]])[1]}."
        ))
      }
    }
  }

  # Check Xrow if provided
  if (!is.null(Xrow)) {
    if (!is.list(Xrow)) {
      cli::cli_abort(c(
        "Xrow must be NULL or a list of matrices.",
        "i" = "Received object of class {class(Xrow)[1]}."
      ))
    }
    if (length(Xrow) != N) {
      cli::cli_abort(c(
        "Xrow must have the same length as Y.",
        "i" = "Y has {N} time periods, but Xrow has {length(Xrow)}."
      ))
    }
    # Check each element is a matrix
    for (t in seq_along(Xrow)) {
      if (!is.null(Xrow[[t]]) && !is.matrix(Xrow[[t]]) && !is.data.frame(Xrow[[t]])) {
        cli::cli_abort(c(
          "All non-NULL elements of Xrow must be matrices or data frames.",
          "i" = "Xrow[[{t}]] is of class {class(Xrow[[t]])[1]}."
        ))
      }
    }
  }

  # Check Xcol if provided
  if (!is.null(Xcol)) {
    if (!is.list(Xcol)) {
      cli::cli_abort(c(
        "Xcol must be NULL or a list of matrices.",
        "i" = "Received object of class {class(Xcol)[1]}."
      ))
    }
    if (length(Xcol) != N) {
      cli::cli_abort(c(
        "Xcol must have the same length as Y.",
        "i" = "Y has {N} time periods, but Xcol has {length(Xcol)}."
      ))
    }
    # Check each element is a matrix
    for (t in seq_along(Xcol)) {
      if (!is.null(Xcol[[t]]) && !is.matrix(Xcol[[t]]) && !is.data.frame(Xcol[[t]])) {
        cli::cli_abort(c(
          "All non-NULL elements of Xcol must be matrices or data frames.",
          "i" = "Xcol[[{t}]] is of class {class(Xcol[[t]])[1]}."
        ))
      }
    }
  }

  # Check dimension consistency within Y
  if (N > 1) {
    first_nrow <- nrow(Y[[1]])
    first_ncol <- ncol(Y[[1]])
    for (t in 2:N) {
      if (nrow(Y[[t]]) != first_nrow || ncol(Y[[t]]) != first_ncol) {
        cli::cli_warn(c(
          "Network dimensions vary across time periods.",
          "i" = "Y[[1]] has dimensions {first_nrow} x {first_ncol}.",
          "i" = "Y[[{t}]] has dimensions {nrow(Y[[t]])} x {ncol(Y[[t]])}.",
          "i" = "This is allowed but actors will be aligned by names."
        ))
        break  # Only warn once
      }
    }
  }

  invisible(TRUE)
}