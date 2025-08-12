#' Check formatting of input objects into ame_repL function
#' 
#' Function used within ame_repL to ensure that input objects are formatted 
#' correctly. Will stop function run if discrepancies found.
#' 
#' @usage check_format(Y, Xdyad, Xrow, Xcol)
#' @param Y a T length list of n x n relational matrices, where T corresponds to the number of replicates (over time, for example). See
#' family below for various data types.
#' @param Xdyad a T length list of n x n x pd arrays of covariates
#' @param Xrow a T length list of n x pr matrices of nodal row covariates
#' @param Xcol a T length list of n x pc matrices of nodal column covariates
#' @author Shahryar Minhas
#' @export
#' 

check_format <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL){
  # make sure inputted data is formatted correctly
  if( !is.list(Y)
  ){ stop('Y needs to be inputted as a list.') }
  
  if(!is.null(Xdyad)){
    if( !is.list(Xdyad)
    ){
      stop('Xdyad needs to be inputted as a list.') } }
  
  if(!is.null(Xrow)){
    if( !is.list(Xrow)
    ){
      stop('Xrow needs to be inputted as a list.') } }
  
  if(!is.null(Xcol)){
    if( !is.list(Xcol)
    ){ stop('Xcol needs to be inputted as a list.') } }
  
  # make sure actor labels are provided
  if(
    do.call('sum', lapply(Y, function(y) !is.null( rownames(y) ) ) )!=length(Y) | 
    do.call('sum', lapply(Y, function(y) !is.null( colnames(y) ) ) )!=length(Y) 
  ){
    stop('Actor labels need to be provided for all Y list objects.') }
  
  if(!is.null(Xdyad)){ 
    if(
      do.call('sum', lapply(Xdyad, function(x) !is.null( dimnames(x)[[1]] ) ) )!=length(Y) | 
      do.call('sum', lapply(Xdyad, function(x) !is.null( dimnames(x)[[2]] ) ) )!=length(Y)
    ){
      stop('Actor labels need to be provided for all Xdyad list objects.') } }
  
  if(!is.null(Xrow)){
    if(
      do.call('sum', lapply(Xrow, function(x) !is.null(rownames(x)) )) != length(Y)
    ){
      stop('Actor labels need to be provided for all Xrow list objects.') } }
  
  if(!is.null(Xcol)){
    if(
      do.call('sum', lapply(Xcol, function(x) !is.null(rownames(x)) )) != length(Y)
    ){
      stop('Actor labels need to be provided for all Xcol list objects.') } }
  
  # make sure actor labels appear in same order across objects
  N <- length(Y)
  
  # Check if this is bipartite (rectangular matrices)
  is_bipartite <- nrow(Y[[1]]) != ncol(Y[[1]])
  
  for(t in 1:N){
    tRowNames <- rownames(Y[[t]])
    tColNames <- colnames(Y[[t]])
    
    if(is_bipartite) {
      # For bipartite, row and column names can be different
      check <- TRUE  # Start with TRUE
      if(!is.null(Xdyad)) {
        check <- c(check, identical(tRowNames, dimnames(Xdyad[[t]])[[1]]))
        check <- c(check, identical(tColNames, dimnames(Xdyad[[t]])[[2]]))
      }
      if(!is.null(Xrow)) {
        check <- c(check, identical(tRowNames, rownames(Xrow[[t]])))
      }
      if(!is.null(Xcol)) {
        check <- c(check, identical(tColNames, rownames(Xcol[[t]])))
      }
    } else {
      # For unipartite, row and column names must be identical
      check <- identical(tRowNames, tColNames)
      if(!is.null(Xdyad)) {
        check <- c(check, identical(tRowNames, dimnames(Xdyad[[t]])[[1]]))
        check <- c(check, identical(tRowNames, dimnames(Xdyad[[t]])[[2]]))
      }
      if(!is.null(Xrow)) {
        check <- c(check, identical(tRowNames, rownames(Xrow[[t]])))
      }
      if(!is.null(Xcol)) {
        check <- c(check, identical(tRowNames, rownames(Xcol[[t]])))
      }
    }
    
    if(sum(check)/length(check) != 1) {
      stop('Actor labels are not identical across inputted data within time periods.')
    }
  }
}