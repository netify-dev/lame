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
		if (!is.null(fit$EZ)) return(.lame_stored_link(fit))
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

# TRUE for a lame fit whose symmetric design keeps <v>_row and <v>_col slices for one stored
# <v>_node coefficient
.lame_pooled_node_fit <- function(fit) {
	fd <- fit[["Xlist", exact = TRUE]]
	B <- fit[["BETA", exact = TRUE]]
	bn <- if (length(dim(B)) == 3L) dimnames(B)[[2L]] else colnames(B)
	isTRUE(fit[["symmetric", exact = TRUE]]) && is.list(fit[["EZ", exact = TRUE]]) && is.list(fd) &&
		length(fd) == length(fit[["EZ", exact = TRUE]]) && !is.null(bn) &&
		length(dim(fd[[1L]])) == 3L && length(bn) < dim(fd[[1L]])[3L]
}

# the stored link fit$EZ of a lame fit. fits saved by lame 1.4.2 and earlier contracted the pooled
# <v>_node coefficients of a symmetric fit against its design positionally, so the link of such a
# fit is rebuilt from its own dyadic slices, as predict(newdata) builds it, keeping the cells the
# stored EZ leaves empty
.lame_stored_link <- function(fit) {
	ez <- fit[["EZ", exact = TRUE]]
	if (!.lame_pooled_node_fit(fit)) return(ez)
	own <- lapply(fit$Xlist, function(X) X[, , grepl("_dyad$|\\.dyad$", dimnames(X)[[3L]]), drop = FALSE])
	link <- predict.lame(fit, newdata = own, type = "link")
	for (t in seq_along(ez)) {
		e <- ez[[t]]
		l <- link[[t]]
		if (!is.null(rownames(e)) && !is.null(rownames(l))) l <- l[rownames(e), colnames(e), drop = FALSE]
		l[is.na(e)] <- NA_real_
		dimnames(l) <- dimnames(e)
		ez[[t]] <- l
	}
	ez
}

# the stored response mean fit$YPM of a lame fit; a Poisson fit forms it from the link, so a
# symmetric Poisson fit with pooled node coefficients takes it from the rebuilt link
.lame_stored_fitted <- function(fit) {
	ypm <- fit[["YPM", exact = TRUE]]
	if (!identical(fit[["family", exact = TRUE]], "poisson") || !.lame_pooled_node_fit(fit)) return(ypm)
	.poisson_fitted_from_fit(fit, EZ = .lame_stored_link(fit))
}
