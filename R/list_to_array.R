# internal: coerce a 2D dyadic covariate into a 3D array with 1 slice
.as_3d_array <- function(M, name = "X1") {
  if (is.null(M)) return(NULL)
  if (length(dim(M)) == 3L) return(M)
  if (is.matrix(M)) {
    A <- array(M, dim = c(nrow(M), ncol(M), 1L))
    dn <- dimnames(M)
    dimnames(A) <- list(
      if (!is.null(dn)) dn[[1]] else NULL,
      if (!is.null(dn)) dn[[2]] else NULL,
      name
    )
    return(A)
  }
  stop("Xdyad element must be matrix or 3D array")
}

# internal: coerce every time-slice to a 3D array, naming 3rd dim if missing
.coerce_Xdyad_3d <- function(Xdyad) {
  if (is.null(Xdyad)) return(NULL)
  lapply(seq_along(Xdyad), function(t) {
    x <- Xdyad[[t]]
    if (is.null(x)) return(NULL)
    if (length(dim(x)) < 3L) {
      .as_3d_array(x, name = "X1")
    } else {
      # ensure third dimnames exists
      dn <- dimnames(x)
      if (length(dn) < 3L || is.null(dn[[3]])) {
        dimnames(x) <- list(dn[[1]], dn[[2]], paste0("X", seq_len(dim(x)[3])))
      }
      x
    }
  })
}

#' Convert list to array
#' 
#' @usage list_to_array(actors, Y, Xdyad, Xrow, Xcol)
#' @param actors vector of actors
#' @param Y dv in list format
#' @param Xdyad dyadic covariates in list format
#' @param Xrow sender covariates in list format
#' @param Xcol receiver covariates in list format
#' @return transforms Y, Xdyad, Xrow, and Xcol to arrays
#' @author Shahryar Minhas
#' 
#' @export list_to_array

list_to_array <- function(actors, Y, Xdyad, Xrow, Xcol){
  
  # dims
  N <- length(Y)
  n <- length(actors)
  pdLabs <- names(Y)
  
  # convert into large array format
  tmp <- array(NA, dim=c(n,n,N),
               dimnames=list( actors, actors, names(Y) ) )
  for(t in 1:N){ tmp[rownames(Y[[t]]),rownames(Y[[t]]),t] <- Y[[t]] }
  Y <- tmp
  
  # Coerce Xdyad to 3D arrays if needed
  Xdyad <- .coerce_Xdyad_3d(Xdyad)
  
  if(!is.null(Xdyad)){
    # Guard dimnames access
    get_third_names <- function(x) {
      dn <- dimnames(x)
      if (length(dn) >= 3L && !is.null(dn[[3]])) dn[[3]] else paste0("X", seq_len(dim(x)[3]))
    }
    
    var_names <- get_third_names(Xdyad[[1]])
    tmp <- array(NA, dim=c(n,n,dim(Xdyad[[1]])[3],N),
                 dimnames=list( actors, actors, var_names, pdLabs ) )
    for(t in 1:N){
      tmp[rownames(Xdyad[[t]]),rownames(Xdyad[[t]]),,t] <- Xdyad[[t]] }
    Xdyad <- tmp
  }

  if(!is.null(Xrow)){
    tmp <- array(NA, dim=c(n, dim(Xrow[[1]])[2], N),
      dimnames=list( actors, colnames(Xrow[[1]]), pdLabs) )
    for(t in 1:N){
      tmp[rownames(Xrow[[t]]),,t] <- Xrow[[t]] }
    Xrow <- tmp
  }

  if(!is.null(Xcol)){
    tmp <- array(NA, dim=c(n, dim(Xcol[[1]])[2], N),
      dimnames=list( actors, colnames(Xcol[[1]]), pdLabs) )
    for(t in 1:N){
      tmp[rownames(Xcol[[t]]),,t] <- Xcol[[t]] }
    Xcol <- tmp
  }

  return( list(Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol) )
}