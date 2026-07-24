#' rank-based z-scores
#'
#' Maps a numeric vector to normal quantiles determined by the ranks of its
#' observed entries, leaving missing values as \code{NA}.
#'
#' @usage zscores(y)
#' @param y a numeric vector
#' @return a numeric vector of the same length as \code{y}
#' @keywords internal
#' @author lame authors
#' @export zscores
zscores <- function(y) {
	scores <- rep(NA_real_, length(y))
	observed <- which(!is.na(y))
	n_obs <- length(observed)
	if (n_obs > 0L) {
		avg_ranks <- rank(y[observed], ties.method = "average")
		scores[observed] <- qnorm(avg_ranks / (n_obs + 1))
	}
	scores
}
