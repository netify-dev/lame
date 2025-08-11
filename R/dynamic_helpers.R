#' Helper functions for dynamic effects in LAME
#' 
#' Internal functions to support dynamic additive and multiplicative effects

#' Compute expected values with dynamic additive effects
#' @keywords internal
#' @noRd
get_EZ_dynamic_ab <- function(Xlist, beta, a_mat, b_mat, U, V, N) {
  EZ <- array(0, dim=c(nrow(a_mat), nrow(a_mat), N))
  
  for(t in 1:N) {
    # Use time-specific additive effects
    ab_t <- outer(a_mat[,t], b_mat[,t], "+")
    EZ[,,t] <- get_EZ_cpp(list(Xlist[[t]]), beta, ab_t, U, V)[,,1]
  }
  
  return(EZ)
}

#' Extract time-specific or averaged additive effects
#' @keywords internal
#' @noRd
get_ab_for_time <- function(a_mat, b_mat, t=NULL, average=FALSE) {
  if(is.null(a_mat) || is.null(b_mat)) {
    return(list(a=0, b=0))
  }
  
  if(average || is.null(t)) {
    # Return time-averaged effects
    return(list(a=rowMeans(a_mat), b=rowMeans(b_mat)))
  } else {
    # Return time-specific effects
    return(list(a=a_mat[,t], b=b_mat[,t]))
  }
}