#' Convert array to list.
#' 
#' @usage arrayToList(arrayObj, actorList, sliceLabel)
#' @param arrayObj 3d array object
#' @param actorList list of actor names
#' @param sliceLabel labels for array slices
#' @return array in list format
#' @author Shahryar Minhas
#' 
#' @export arrayToList
arrayToList <- function(arrayObj, actorList, sliceLabel){
  listObj <- lapply(1:dim(arrayObj)[3], function(t){
    actorT <- actorList[[t]]
    mat <- arrayObj[actorT,actorT,t]
    diag(mat) <- NA
    return(mat) })
  names(listObj) <- sliceLabel
  return(listObj)
}