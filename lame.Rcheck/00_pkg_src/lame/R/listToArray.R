#' Convert list to array
#' 
#' @usage listToArray(actors, Y, Xdyad, Xrow, Xcol)
#' @param actors vector of actors
#' @param Y dv in list format
#' @param Xdyad dyadic covariates in list format
#' @param Xrow sender covariates in list format
#' @param Xcol receiver covariates in list format
#' @return transforms Y, Xdyad, Xrow, and Xcol to arrays
#' @author Shahryar Minhas
#' 
#' @export listToArray

listToArray <- function(actors, Y, Xdyad, Xrow, Xcol){
  
  # dims
  N <- length(Y)
  n <- length(actors)
  pdLabs <- names(Y)
  
  # convert into large array format
  tmp <- array(NA, dim=c(n,n,N),
               dimnames=list( actors, actors, names(Y) ) )
  for(t in 1:N){ tmp[rownames(Y[[t]]),rownames(Y[[t]]),t] <- Y[[t]] }
  Y <- tmp
  
  if(!is.null(Xdyad)){
    tmp <- array(NA, dim=c(n,n,dim(Xdyad[[1]])[3],N),
                 dimnames=list( actors, actors, dimnames(Xdyad[[1]])[[3]], pdLabs ) )
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