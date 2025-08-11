#' Convert array to list.
#' 
#' @usage array_to_list(arrayObj, actorList, sliceLabel)
#' @param arrayObj 3d array object
#' @param actorList list of actor names
#' @param sliceLabel labels for array slices
#' @return array in list format
#' @author Shahryar Minhas
#' 
#' @export array_to_list
array_to_list <- function(arrayObj, actorList, sliceLabel){
  listObj <- lapply(1:dim(arrayObj)[3], function(t){
    actorT <- actorList[[t]]
    mat <- arrayObj[actorT,actorT,t]
    diag(mat) <- NA
    return(mat) })
  names(listObj) <- sliceLabel
  return(listObj)
}