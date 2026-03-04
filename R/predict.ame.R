####
#' Predict method for AME models
#' 
#' @description
#' Generate predictions from fitted AME models, including point estimates
#' and predictive distributions.
#' 
#' @param object Fitted AME model object
#' @param newdata Optional new data (covariates) for prediction
#' @param type Character; type of prediction:
#'   - "response": predicted values on response scale (default)
#'   - "link": predicted values on link scale
#'   - "distribution": full posterior predictive distribution
#' @param n_samples For type="distribution", number of posterior samples
#' @param include_uncertainty Logical; include parameter uncertainty (default TRUE)
#' @param ... Additional arguments (not used)
#' 
#' @return Depending on type:
#'   - "response"/"link": Matrix of predictions
#'   - "distribution": Array of posterior predictive samples
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @examples
#' \donttest{
#' # Fit model
#' data(YX_nrm)
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            nscan = 100, burn = 10, odens = 1, print = FALSE)
#'
#' # Point predictions
#' Y_pred <- predict(fit)
#'
#' # Predictions on link scale
#' Y_link <- predict(fit, type = "link")
#' }
#' 
#' @export
predict.ame <- function(
	object, 
	newdata = NULL,
	type = c("response", "link", "distribution"),
	n_samples = 100,
	include_uncertainty = TRUE, ...
	){
	
	type <- match.arg(type)
	
	# network dimensions
	n <- ifelse(object$mode == "bipartite", object$nA, nrow(object$YPM))
	m <- ifelse(object$mode == "bipartite", object$nB, ncol(object$YPM))
	
	# validate new covariates
	if(!is.null(newdata)) {
		if(!is.array(newdata) || length(dim(newdata)) != 3) {
			cli::cli_abort("newdata must be a 3-dimensional array")
		}
		if(dim(newdata)[1] != n || dim(newdata)[2] != m) {
			cli::cli_abort("newdata dimensions don't match model ({.val {n}} x {.val {m}})")
		}
		X <- newdata
	} else {
		X <- NULL
	}
	
	if(type == "distribution") {
		return(predict_distribution(object, X, n_samples, include_uncertainty))
	} else {
		if(!is.null(X)) {
			# linear predictor from new covariates
			beta_mean <- colMeans(object$BETA)
			beta_names <- colnames(object$BETA)
			has_intercept <- !is.null(beta_names) && beta_names[1] == "intercept"
			XB <- array(0, dim = c(n, m))

			if(has_intercept) {
				# intercept at index 1, covariates start at index 2
				for(k in 1:dim(X)[3]) {
					XB <- XB + X[,,k] * beta_mean[k+1]
				}
				XB <- XB + beta_mean[1]
			} else {
				# no intercept: covariates start at index 1
				for(k in 1:dim(X)[3]) {
					XB <- XB + X[,,k] * beta_mean[k]
				}
			}
			
			# add sender/receiver random effects
			if(!is.null(object$APM)) {
				for(i in 1:n) XB[i,] <- XB[i,] + object$APM[i]
			}
			if(!is.null(object$BPM)) {
				for(j in 1:m) XB[,j] <- XB[,j] + object$BPM[j]
			}
			
			# multiplicative effects
			UVPM <- reconstruct_UVPM(object)
			if(!is.null(UVPM)) {
				XB <- XB + UVPM
			}
			
			EZ <- XB
		} else {
			if(type == "link") {
				EZ <- reconstruct_EZ(object)
				return(EZ)
			} else {
				# apply inverse link for response scale
				EZ <- reconstruct_EZ(object)
				if(object$family == "binary") {
					return(pnorm(EZ))  # probit link
				} else if(object$family == "normal") {
					return(EZ)  # identity link
				} else if(object$family == "poisson") {
					return(exp(EZ))  # log link
				} else if(object$family == "tobit") {
					return(pmax(0, EZ))  # censored at 0
				} else {
					return(object$YPM)
				}
			}
		}
	}
}

####
#' Generate posterior predictive distribution
#' @noRd
predict_distribution <- function(object, X, n_samples, include_uncertainty) {

	n <- ifelse(object$mode == "bipartite", object$nA, nrow(object$YPM))
	m <- ifelse(object$mode == "bipartite", object$nB, ncol(object$YPM))

	beta_names <- colnames(object$BETA)
	has_intercept <- !is.null(beta_names) && beta_names[1] == "intercept"

	Y_pred <- array(NA, dim = c(n, m, n_samples))

	if(include_uncertainty && !is.null(object$BETA)) {
		n_mcmc <- nrow(object$BETA)
		if(n_mcmc < n_samples) {
			cli::cli_warn("Only {.val {n_mcmc}} MCMC samples available, using all")
			n_samples <- n_mcmc
		}

		sample_idx <- sample(n_mcmc, n_samples, replace = (n_samples > n_mcmc))

		for(i in 1:n_samples) {
			beta <- object$BETA[sample_idx[i], ]
			if(!is.null(X)) {
				XB <- array(0, dim = c(n, m))
				if(has_intercept) {
					for(k in 1:dim(X)[3]) {
						XB <- XB + X[,,k] * beta[k+1]
					}
					XB <- XB + beta[1]
				} else {
					for(k in 1:dim(X)[3]) {
						XB <- XB + X[,,k] * beta[k]
					}
				}
			} else {
				if(has_intercept) {
					XB <- matrix(beta[1], n, m)
				} else {
					XB <- matrix(0, n, m)
				}
			}
			
			if(!is.null(object$APM)) {
				XB <- XB + object$APM
			}
			if(!is.null(object$BPM)) {
				XB <- sweep(XB, 2, object$BPM, "+")
			}
			
			if(!is.null(object$U) && !is.null(object$V)) {
				if(!is.null(object$U_samples) && !is.null(object$V_samples)) {
					idx <- min(i, dim(object$U_samples)[3])
					UV <- object$U_samples[,,idx] %*% t(object$V_samples[,,idx])
				} else {
					# jitter posterior means
					s2 <- object$VC[sample_idx[i], ncol(object$VC)]
					R <- ncol(object$U)
					U_sim <- object$U + matrix(rnorm(n * R, 0, sqrt(s2/R/10)), n, R)
					V_sim <- object$V + matrix(rnorm(m * R, 0, sqrt(s2/R/10)), m, R)
					UV <- U_sim %*% t(V_sim)
				}
				XB <- XB + UV
			}
			
			# observation noise
			s2 <- object$VC[sample_idx[i], ncol(object$VC)]
			Z <- XB + matrix(rnorm(n * m, 0, sqrt(s2)), n, m)
			
			# apply inverse link
			Y_pred[,,i] <- transform_to_response(Z, object$family)
		}
	} else {
		EZ <- if(!is.null(X)) {
			beta_mean <- colMeans(object$BETA)
			XB <- array(0, dim = c(n, m))
			if(has_intercept) {
				for(k in 1:dim(X)[3]) {
					XB <- XB + X[,,k] * beta_mean[k+1]
				}
				XB + beta_mean[1]
			} else {
				for(k in 1:dim(X)[3]) {
					XB <- XB + X[,,k] * beta_mean[k]
				}
				XB
			}
		} else {
			reconstruct_EZ(object)
		}
		
		s2 <- mean(object$VC[, ncol(object$VC)])
		
		for(i in 1:n_samples) {
			Z <- EZ + matrix(rnorm(n * m, 0, sqrt(s2)), n, m)
			Y_pred[,,i] <- transform_to_response(Z, object$family)
		}
	}
	
	return(Y_pred)
}

####
#' Transform linear predictor to response scale
#' @noRd
transform_to_response <- function(Z, family) {
	switch(family,
		normal = Z,
		binary = pnorm(Z),
		poisson = exp(Z),
		tobit = pmax(0, Z),
		ordinal = Z,  # proper ordinal transform needs cutpoints we don't have here
		cbin = pnorm(Z),
		frn = Z,
		rrl = Z,
		Z  # default passthrough
	)
}

####
#' Extract fitted values from AME model
#' 
#' @param object Fitted AME model
#' @param ... Additional arguments
#' 
#' @return Matrix of fitted values
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
fitted.ame <- function(object, ...) {
	return(object$YPM)
}

####
#' Extract residuals from AME model
#' 
#' @param object Fitted AME model
#' @param type Type of residuals ("response" or "pearson")
#' @param ... Additional arguments
#' 
#' @return Matrix of residuals
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
residuals.ame <- function(object, type = c("response", "pearson"), ...) {
	type <- match.arg(type)

	if(!is.null(object$Y)) {
		Y <- object$Y
	} else if(!is.null(object$Y_obs)) {
		Y <- object$Y_obs
	} else {
		cli::cli_warn("original Y not stored in model object, residuals may be incorrect")
		Y <- object$YPM
	}

	if(type == "response") {
		return(Y - object$YPM)
	} else {
		# pearson residuals scaled by variance
		s2 <- mean(object$VC[, ncol(object$VC)])
		return((Y - object$YPM) / sqrt(s2))
	}
}

####
#' Predict method for LAME models
#'
#' @description
#' Generate predictions from fitted longitudinal AME models.
#' Returns a list of matrices (one per time point) on the requested scale.
#'
#' @param object Fitted LAME model object
#' @param type Character; "response" or "link"
#' @param ... Additional arguments (not used)
#'
#' @return List of prediction matrices (one per time point)
#'
#' @method predict lame
#' @export
predict.lame <- function(object, type = c("response", "link"), ...) {
	type <- match.arg(type)

	EZ <- object$EZ
	if(is.null(EZ)) {
		cli::cli_warn("EZ not available in model object")
		return(object$YPM)
	}

	if(type == "link") {
		return(EZ)
	}

	# apply inverse link per time point
	family <- object$family %||% "normal"
	lapply(EZ, function(ez) {
		switch(family,
			binary = pnorm(ez),
			poisson = exp(ez),
			tobit = pmax(0, ez),
			ez
		)
	})
}

####
#' Extract fitted values from LAME model
#'
#' @param object Fitted LAME model
#' @param ... Additional arguments
#'
#' @return List of fitted value matrices (one per time point)
#'
#' @method fitted lame
#' @export
fitted.lame <- function(object, ...) {
	return(object$YPM)
}

####
#' Extract residuals from LAME model
#'
#' @param object Fitted LAME model
#' @param type Type of residuals ("response" or "pearson")
#' @param ... Additional arguments
#'
#' @return List of residual matrices (one per time point)
#'
#' @method residuals lame
#' @export
residuals.lame <- function(object, type = c("response", "pearson"), ...) {
	type <- match.arg(type)

	if(!is.null(object$Y)) {
		Y <- object$Y
	} else {
		cli::cli_warn("original Y not stored in model object, residuals may be incorrect")
		Y <- object$YPM
	}

	# Y may be 3D array; convert to list of matrices
	if(is.array(Y) && length(dim(Y)) == 3) {
		n_time <- dim(Y)[3]
		Y <- lapply(seq_len(n_time), function(t) Y[,,t])
	}

	# Y and YPM are lists for lame
	s2 <- mean(object$VC[, ncol(object$VC)])

	mapply(function(y, ypm) {
		if(type == "response") {
			y - ypm
		} else {
			(y - ypm) / sqrt(s2)
		}
	}, Y, object$YPM, SIMPLIFY = FALSE)
}
####