#' Format AME call for display
#'
#' Creates a readable formula-like representation of the AME model
#'
#' @param fit An AME fit object
#' @return A character string representing the model formula
#' @keywords internal
#' @noRd
format_ame_call <- function(fit) {
	call_str <- "Y ~ "

	####
	# add covariates if present. bipartite fits don't set fit$x_names but do set
	# colnames(fit$beta); read that as a fallback so the call line shows the
	# real covariates rather than only "a[i] + b[j]".
	x_names <- if (!is.null(fit$X_names) && length(fit$X_names) > 0) {
		fit$X_names
	} else if (!is.null(colnames(fit$BETA))) {
		setdiff(colnames(fit$BETA), "intercept")
	} else {
		character(0)
	}
	covars <- c()
	if(length(x_names) > 0) {
		dyad_vars <- grep("_dyad$", x_names, value = TRUE)
		row_vars <- grep("_row$", x_names, value = TRUE)
		col_vars <- grep("_col$", x_names, value = TRUE)
		other_vars <- setdiff(x_names, c(dyad_vars, row_vars, col_vars))

		if(length(other_vars) > 0) {
			covars <- c(covars, other_vars)
		}
		if(length(dyad_vars) > 0) {
			covars <- c(covars, paste0("dyad(", paste(gsub("_dyad$", "", dyad_vars), collapse=", "), ")"))
		}
		if(length(row_vars) > 0) {
			covars <- c(covars, paste0("row(", paste(gsub("_row$", "", row_vars), collapse=", "), ")"))
		}
		if(length(col_vars) > 0) {
			covars <- c(covars, paste0("col(", paste(gsub("_col$", "", col_vars), collapse=", "), ")"))
		}
	}
	####

	####
	# add random effects
	effects <- c()
	if(!is.null(fit$rvar) && fit$rvar) effects <- c(effects, "a[i]")
	if(!is.null(fit$cvar) && fit$cvar) effects <- c(effects, "b[j]")
	if(!is.null(fit$dcor) && fit$dcor) effects <- c(effects, "rho*e[ji]")
	# bipartite fits store r_row / r_col separately and may leave r unset; surface
	# the multiplicative term in either case so the printed call matches the fit
	is_bipartite <- identical(fit$mode, "bipartite")
	if (is_bipartite) {
		rr <- if (!is.null(fit$R_row)) fit$R_row else if (!is.null(fit$R)) fit$R else 0
		rc <- if (!is.null(fit$R_col)) fit$R_col else if (!is.null(fit$R)) fit$R else 0
		if (rr > 0 && rc > 0) {
			effects <- c(effects, paste0("U[i,1:", rr, "] %*% G %*% V[j,1:", rc, "]'"))
		}
	} else if (!is.null(fit$R) && fit$R > 0) {
		if (!is.null(fit$symmetric) && fit$symmetric) {
			effects <- c(effects, paste0("U[i,1:", fit$R, "] %*% L %*% U[j,1:", fit$R, "]"))
		} else {
			effects <- c(effects, paste0("U[i,1:", fit$R, "] %*% V[j,1:", fit$R, "]"))
		}
	}
	####

	####
	# combine into formula string
	if(length(covars) > 0) {
		call_str <- paste0(call_str, paste(covars, collapse = " + "))
	}

	if(length(effects) > 0) {
		if(length(covars) > 0) call_str <- paste0(call_str, " + ")
		call_str <- paste0(call_str, paste(effects, collapse = " + "))
	}

	if(!is.null(fit$family)) {
		call_str <- paste0(call_str, ", family = '", fit$family, "'")
	}
	####

	return(call_str)
}
