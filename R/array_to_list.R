#' Convert array to list.
#'
#' @usage array_to_list(arrayObj, actorList, sliceLabel)
#' @param arrayObj 3d array object
#' @param actorList list of actor names
#' @param sliceLabel labels for array slices
#' @return array in list format
#' @examples
#' n <- 12
#' arr <- array(rnorm(n * n * 3), dim = c(n, n, 3))
#' rownames(arr) <- colnames(arr) <- paste0("a", 1:n)
#' actors <- rep(list(paste0("a", 1:n)), 3)
#' Y_list <- array_to_list(arr, actors, sliceLabel = paste0("t", 1:3))
#' length(Y_list)
#' @author Shahryar Minhas
#'
#' @export array_to_list
array_to_list <- function(arrayObj, actorList, sliceLabel) {
	list_obj <- lapply(1:dim(arrayObj)[3], function(t) {
		actor_t <- actorList[[t]]
		all_actors <- rownames(arrayObj)
		if (is.null(all_actors)) {
			all_actors <- 1:nrow(arrayObj)
		}
		actor_indices <- match(actor_t, all_actors)
		actor_indices <- actor_indices[!is.na(actor_indices)]

		mat <- arrayObj[actor_indices, actor_indices, t, drop = FALSE]
		if (length(dim(mat)) == 3) {
			mat <- mat[, , 1]
		}
		rownames(mat) <- colnames(mat) <- actor_t[actor_t %in% all_actors]
		diag(mat) <- NA
		return(mat) })
	names(list_obj) <- sliceLabel
	return(list_obj)
}
