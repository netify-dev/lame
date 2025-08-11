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
    # Get indices for actors present at time t
    all_actors <- rownames(arrayObj)
    if(is.null(all_actors)) {
      all_actors <- 1:nrow(arrayObj)
    }
    actor_indices <- match(actorT, all_actors)
    # Remove NAs from indices (actors not in full set)
    actor_indices <- actor_indices[!is.na(actor_indices)]
    
    # Extract submatrix for actors present at time t
    mat <- arrayObj[actor_indices, actor_indices, t, drop=FALSE]
    if(length(dim(mat)) == 3) {
      mat <- mat[,,1]
    }
    rownames(mat) <- colnames(mat) <- actorT[actorT %in% all_actors]
    diag(mat) <- NA
    return(mat) })
  names(listObj) <- sliceLabel
  return(listObj)
}