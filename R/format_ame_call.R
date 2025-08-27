#' Format AME call for display
#' 
#' Creates a readable formula-like representation of the AME model
#' 
#' @param fit An AME fit object
#' @return A character string representing the model formula
#' @keywords internal
#' @noRd
format_ame_call <- function(fit) {
  # Start with network name
  call_str <- "Y ~ "
  
  # Add covariates if present
  covars <- c()
  if(!is.null(fit$X_names) && length(fit$X_names) > 0) {
    # Group by type
    dyad_vars <- grep("_dyad$", fit$X_names, value = TRUE)
    row_vars <- grep("_row$", fit$X_names, value = TRUE)
    col_vars <- grep("_col$", fit$X_names, value = TRUE)
    other_vars <- setdiff(fit$X_names, c(dyad_vars, row_vars, col_vars))
    
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
  
  # Add random effects
  effects <- c()
  if(!is.null(fit$rvar) && fit$rvar) effects <- c(effects, "a[i]")
  if(!is.null(fit$cvar) && fit$cvar) effects <- c(effects, "b[j]")
  if(!is.null(fit$dcor) && fit$dcor) effects <- c(effects, "rho*e[ji]")
  if(!is.null(fit$R) && fit$R > 0) {
    if(!is.null(fit$symmetric) && fit$symmetric) {
      effects <- c(effects, paste0("U[i,1:", fit$R, "] %*% L %*% U[j,1:", fit$R, "]"))
    } else {
      effects <- c(effects, paste0("U[i,1:", fit$R, "] %*% V[j,1:", fit$R, "]"))
    }
  }
  
  # Combine
  if(length(covars) > 0) {
    call_str <- paste0(call_str, paste(covars, collapse = " + "))
  }
  
  if(length(effects) > 0) {
    if(length(covars) > 0) call_str <- paste0(call_str, " + ")
    call_str <- paste0(call_str, paste(effects, collapse = " + "))
  }
  
  # Add family info
  if(!is.null(fit$family)) {
    call_str <- paste0(call_str, ", family = '", fit$family, "'")
  }
  
  return(call_str)
}