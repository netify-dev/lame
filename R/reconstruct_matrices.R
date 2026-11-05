#' Reconstruct EZ and UVPM matrices from AME model output
#'
#' @description
#' Helper functions to recover EZ (linear predictor) and UVPM (posterior
#' mean of the multiplicative product) matrices. When the fit stores the
#' quantity (\code{UVPM} on asymmetric fits, \code{ULUPM} on symmetric
#' fits), the stored posterior mean is returned directly; otherwise it is
#' reconstructed from the saved factors.
#'
#' Note: EZ returns the linear predictor (\eqn{\eta}), not the response:
#' - For Gaussian: \eqn{EZ = \eta = \mu} (identity link)
#' - For Poisson: \eqn{EZ = \eta = \log(\lambda)} (can be negative)
#' - For Binary: \eqn{EZ = \eta} = probit inverse of p (can be any real value)
#' Use YPM for predictions on the response scale.
#'
#' @param fit Fitted AME model object
#' @param X Covariate array (optional, will use \code{fit$X} if available)
#'
#' @return Reconstructed matrix
#'
#' @details
#' Fits that carry the stored posterior mean (\code{UVPM} or, for
#' symmetric fits, \code{ULUPM}) have it returned as-is; older fit objects
#' without one are rebuilt from the posterior factor means.
#'
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' EZ <- reconstruct_EZ(fit)
#' UVPM <- reconstruct_UVPM(fit)
#' dim(EZ)
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
reconstruct_EZ <- function(fit, X = NULL) {
	# lame fits: YPM is a list, so return the stored EZ
	if (!is.null(fit$YPM) && is.list(fit$YPM)) {
		if (!is.null(fit$EZ)) return(fit$EZ)
		cli::cli_abort("reconstruct_EZ for lame objects requires fit$EZ to be stored")
	}
	# amen fits store EZ directly and carry no mode or symmetric metadata
	if (is.null(fit[["mode", exact = TRUE]]) && !is.null(fit[["EZ", exact = TRUE]]) &&
			is.matrix(fit[["EZ", exact = TRUE]])) {
		return(fit[["EZ", exact = TRUE]])
	}
	mode <- fit[["mode", exact = TRUE]]
	symmetric <- isTRUE(fit[["symmetric", exact = TRUE]])
	if (!is.null(fit$YPM)) {
		n <- nrow(fit$YPM)
		m <- ncol(fit$YPM)
	} else {
		cli::cli_abort("Cannot determine dimensions from fit object")
	}

	if (!is.null(fit$BETA)) {
		beta_mean <- colMeans(fit$BETA)

		if (!is.null(fit[["X"]])) {
			X <- fit[["X"]]
		}

		if (!is.null(X)) {
			# symmetric fits keep one merged nodal coefficient per covariate, so map
			# each design slice to its coefficient by name
			beta_slice <- .beta_on_design(beta_mean, dimnames(X)[[3]])
			EZ <- matrix(0, n, m)
			for (k in seq_len(dim(X)[3])) {
				if (k <= length(beta_slice)) {
					EZ <- EZ + X[, , k] * beta_slice[k]
				}
			}
		} else {
			EZ <- matrix(beta_mean[1], n, m)
		}
	} else {
		EZ <- matrix(0, n, m)
	}

	# additive effects
	if (!is.null(fit$APM)) {
		if (identical(mode, "bipartite")) {
			for (i in 1:n) EZ[i, ] <- EZ[i, ] + fit$APM[i]
			if (!is.null(fit$BPM)) {
				for (j in 1:m) EZ[, j] <- EZ[, j] + fit$BPM[j]
			}
		} else if (symmetric) {
			EZ <- EZ + outer(fit$APM, fit$APM, "+")
		} else {
			if (!is.null(fit$BPM)) {
				EZ <- EZ + outer(fit$APM, fit$BPM, "+")
			}
		}
	}

	# multiplicative effects: prefer the stored posterior-mean product (ULUPM on
	# symmetric fits, UVPM otherwise), as predict() does with newdata
	if (!is.null(fit[["ULUPM"]])) {
		EZ <- EZ + fit[["ULUPM"]]
	} else if (!is.null(fit[["UVPM"]])) {
		EZ <- EZ + fit[["UVPM"]]
	} else if (!is.null(fit$U)) {
		if (identical(mode, "bipartite") && !is.null(fit[["G"]]) && !is.null(fit[["V"]])) {
			EZ <- EZ + fit$U %*% fit[["G"]] %*% t(fit[["V"]])
		} else if (symmetric && !is.null(fit$L)) {
			if (is.matrix(fit$L)) {
				EZ <- EZ + fit$U %*% fit$L %*% t(fit$U)
			} else if (is.numeric(fit$L)) {
				EZ <- EZ + fit$U %*% diag(fit$L, nrow = length(fit$L)) %*% t(fit$U)
			}
		} else if (!is.null(fit[["V"]])) {
			EZ <- EZ + fit$U %*% t(fit[["V"]])
		}
	}

	return(EZ)
}

#' @rdname reconstruct_EZ
#' @export
reconstruct_UVPM <- function(fit) {
	# prefer the stored posterior mean (ULUPM on symmetric fits)
	if (!is.null(fit$UVPM)) {
		return(fit$UVPM)
	}
	if (!is.null(fit$ULUPM)) {
		return(fit$ULUPM)
	}
	if (is.null(fit$U)) {
		return(NULL)
	}

	# bipartite lame stores U/V as n x R x T cubes even for static fits; take the
	# last slice and keep the matrix shape for rank-1 sides
	.uv_last_slice <- function(A) {
		if (is.null(A) || length(dim(A)) != 3L) return(A)
		matrix(
			A[, , dim(A)[3L]],
			nrow = dim(A)[1L],
			ncol = dim(A)[2L],
			dimnames = list(dimnames(A)[[1]], NULL)
		)
	}
	fit$U <- .uv_last_slice(fit$U)
	fit[["V"]] <- .uv_last_slice(fit[["V"]])

	if (isTRUE(fit[["symmetric", exact = TRUE]])) {
		if (!is.null(fit$L) && !is.null(fit$U)) {
			if (is.matrix(fit$L)) {
				return(fit$U %*% fit$L %*% t(fit$U))
			} else if (is.numeric(fit$L)) {
				return(fit$U %*% diag(fit$L, nrow = length(fit$L)) %*% t(fit$U))
			}
		}
	} else {
		if (!is.null(fit[["V"]])) {
			if (!is.null(fit[["G", exact = TRUE]]) &&
					identical(fit[["mode", exact = TRUE]], "bipartite")) {
				if (ncol(fit$U) == nrow(fit[["G"]]) && ncol(fit[["G"]]) == ncol(fit[["V"]])) {
					return(fit$U %*% fit[["G"]] %*% t(fit[["V"]]))
				} else {
					return(fit$U %*% t(fit[["V"]]))
				}
			} else {
				return(fit$U %*% t(fit[["V"]]))
			}
		}
	}

	return(NULL)
}
