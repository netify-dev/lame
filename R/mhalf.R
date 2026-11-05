#' Symmetric square root of a matrix
#'
#' Computes the symmetric square root of a positive semidefinite matrix
#'
#'
#' @usage mhalf(M)
#' @param M a positive semidefinite matrix
#' @return a matrix \code{H} such that \code{H^2} equals \code{M}
#' @examples
#' M <- matrix(c(4, 1, 1, 4), 2, 2)
#' H <- mhalf(M)
#' H %*% H  # recovers M
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export mhalf
mhalf <-
	function(M) {
		if (!is.matrix(M) || !is.numeric(M) || nrow(M) == 0L ||
				nrow(M) != ncol(M) || any(!is.finite(M))) {
			cli::cli_abort("M must be a finite nonempty square matrix.")
		}
		M <- (M + t(M)) / 2
		tmp <- eigen(M, symmetric = TRUE)
		# tolerate eigensolver roundoff around zero
		tolerance <- 8 * nrow(M) * .Machine$double.eps * max(abs(tmp$values))
		if (any(tmp$values < -tolerance)) {
			cli::cli_abort("M must be positive semidefinite.")
		}
		D <- pmax(tmp$values, 0)
		tmp$vectors %*% diag(sqrt(D), nrow(M)) %*% t(tmp$vectors)
	}
