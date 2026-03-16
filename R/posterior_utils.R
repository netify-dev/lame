#' Options for saving posterior samples during MCMC
#'
#' @param save_UV Logical; whether to save samples of U and V matrices (default FALSE)
#' @param save_ab Logical; whether to save samples of additive effects a and b (default FALSE)
#' @param thin_UV Integer; thinning interval for U/V samples (default 10)
#' @param thin_ab Integer; thinning interval for a/b samples (default 10)
#' 
#' @return List of posterior saving options
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
posterior_options <- function(
	save_UV = FALSE, save_ab = FALSE,
	thin_UV = 10, thin_ab = 10
	){
	list(
		save_UV = save_UV,
		save_ab = save_ab,
		thin_UV = thin_UV,
		thin_ab = thin_ab
	)
}

####
# simulate_posterior
####

#' Simulate posterior distributions from fitted AME model
#'
#' @param fit Fitted ame model object
#' @param component Character; which component to simulate: "UV", "ab", "beta", "Y"
#' @param n_samples Number of posterior samples to generate
#' @param use_full_cov For UV component, whether to use empirical variance scaling (default FALSE)
#' @param seed Random seed for reproducibility
#' 
#' @return Array or matrix of posterior samples
#' 
#' @details
#' This function can simulate posterior distributions even when they weren't
#' saved during MCMC, by using the posterior means and variance components.
#' 
#' For more accurate posteriors, use posterior_options() during model fitting
#' to save the actual MCMC samples.
#' 
#' @examples
#' \donttest{
#' # Fit a model with multiplicative effects
#' data(YX_nrm)
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE)
#'
#' # Get posterior samples of regression coefficients
#' beta_post <- simulate_posterior(fit, "beta", n_samples = 50)
#' }
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
simulate_posterior <- function(
	fit, 
	component = c("UV", "ab", "beta", "Y"),
	n_samples = 100,
	use_full_cov = FALSE,
	seed = NULL
	){
	
	component <- match.arg(component)
	
	if(!is.null(seed)) set.seed(seed)
	
	has_samples <- switch(component,
		UV = !is.null(fit$U_samples) && !is.null(fit$V_samples),
		ab = !is.null(fit$a_samples) && !is.null(fit$b_samples),
		beta = !is.null(fit$BETA),
		Y = FALSE  # Always simulate Y
	)
	
	if(has_samples && component != "Y") {
		cli::cli_text("Using saved MCMC samples for {.field {component}}")
		
		if(component == "UV") {
			return(sample_from_saved_UV(fit, n_samples))
		} else if(component == "ab") {
			return(sample_from_saved_ab(fit, n_samples))
		} else if(component == "beta") {
			return(sample_from_saved_beta(fit, n_samples))
		}
		
	} else {
		cli::cli_text("Simulating from approximate posterior for {.field {component}}")
		
		if(component == "UV") {
			return(simulate_UV_posterior(fit, n_samples, use_full_cov))
		} else if(component == "ab") {
			return(simulate_ab_posterior(fit, n_samples))
		} else if(component == "beta") {
			return(simulate_beta_posterior(fit, n_samples))
		} else if(component == "Y") {
			return(simulate_Y_posterior(fit, n_samples))
		}
	}
}

####
# sample_from_saved_UV
####

#' Sample from saved U/V posteriors
#' @noRd
sample_from_saved_UV <- function(fit, n_samples) {
	n_saved <- dim(fit$U_samples)[3]
	
	if(n_samples > n_saved) {
		cli::cli_warn("Requested {.val {n_samples}} samples but only {.val {n_saved}} available")
		n_samples <- n_saved
	}
	
	idx <- sample(n_saved, n_samples, replace = (n_samples > n_saved))

	n_row <- nrow(fit$U_samples)
	n_col <- nrow(fit$V_samples)
	
	UV_samples <- array(NA, dim = c(n_row, n_col, n_samples))
	
	for(i in 1:n_samples) {
		UV_samples[,,i] <- fit$U_samples[,,idx[i]] %*% t(fit$V_samples[,,idx[i]])
	}
	
	return(UV_samples)
}

####
# simulate_UV_posterior
####

#' Simulate U/V from approximate posterior
#' @noRd
simulate_UV_posterior <- function(fit, n_samples, use_full_cov = FALSE) {

	if(is.null(fit$U) || is.null(fit$V)) {
		cli::cli_abort("No multiplicative effects in model")
	}

	n_row <- nrow(fit$U)
	n_col <- nrow(fit$V)
	R_u <- ncol(fit$U)
	R_v <- ncol(fit$V)

	if(!is.null(fit$VC)) {
		s2 <- mean(fit$VC[, ncol(fit$VC)])
	} else {
		s2 <- 1
	}

	UV_samples <- array(NA, dim = c(n_row, n_col, n_samples))

	if(use_full_cov && !is.null(fit$BETA) && nrow(fit$BETA) >= 100) {
		n_mcmc <- nrow(fit$BETA)
		sd_scale <- sqrt(s2 / R_u) * sqrt(1 + 50/n_mcmc)
		cli::cli_inform("Using empirical variance scaling based on {.val {n_mcmc}} MCMC samples")
	} else {
		sd_scale <- sqrt(s2 / R_u)
	}

	bip <- !is.null(fit$mode) && fit$mode == "bipartite"

	for(s in 1:n_samples) {
		U_sim <- fit$U + matrix(rnorm(n_row * R_u, 0, sd_scale), n_row, R_u)
		V_sim <- fit$V + matrix(rnorm(n_col * R_v, 0, sd_scale), n_col, R_v)

		if(bip && !is.null(fit$G)) {
			UV_samples[,,s] <- U_sim %*% fit$G %*% t(V_sim)
		} else {
			UV_samples[,,s] <- U_sim %*% t(V_sim)
		}
	}

	return(UV_samples)
}

####
# sample_from_saved_ab
####

#' Sample from saved a/b posteriors
#' @noRd
sample_from_saved_ab <- function(fit, n_samples) {
	n_saved <- ncol(fit$a_samples)
	
	if(n_samples > n_saved) {
		cli::cli_warn("Requested {.val {n_samples}} samples but only {.val {n_saved}} available")
		n_samples <- n_saved
	}
	
	idx <- sample(n_saved, n_samples, replace = (n_samples > n_saved))

	ab_samples <- list(
		a = fit$a_samples[, idx],
		b = fit$b_samples[, idx]
	)
	
	return(ab_samples)
}

####
# simulate_ab_posterior
####

#' Simulate a/b from approximate posterior
#' @noRd
simulate_ab_posterior <- function(fit, n_samples) {

	if(is.null(fit$APM) || is.null(fit$BPM)) {
		cli::cli_abort("No additive effects in model")
	}

	nA <- length(fit$APM)
	nB <- length(fit$BPM)

	if(!is.null(fit$VC) && ncol(fit$VC) >= 2) {
		va <- mean(fit$VC[, 1])
		# separate variance for b if available (col 3 for asymmetric)
		if(ncol(fit$VC) >= 3) {
			vb <- mean(fit$VC[, 3])
		} else {
			vb <- va
		}
	} else {
		va <- 1
		vb <- 1
	}

	ab_samples <- list(
		a = matrix(NA, nA, n_samples),
		b = matrix(NA, nB, n_samples)
	)

	for(s in 1:n_samples) {
		ab_samples$a[, s] <- fit$APM + rnorm(nA, 0, sqrt(va))
		ab_samples$b[, s] <- fit$BPM + rnorm(nB, 0, sqrt(vb))
	}

	return(ab_samples)
}

####
# sample_from_saved_beta
####

#' Sample from saved beta posteriors
#' @noRd
sample_from_saved_beta <- function(fit, n_samples) {
	n_saved <- nrow(fit$BETA)
	
	if(n_samples > n_saved) {
		cli::cli_warn("Requested {.val {n_samples}} samples but only {.val {n_saved}} available")
		n_samples <- n_saved
	}
	
	idx <- sample(n_saved, n_samples, replace = (n_samples > n_saved))

	return(fit$BETA[idx, , drop = FALSE])
}

####
# simulate_beta_posterior
####

#' Simulate beta from approximate posterior
#' @noRd
simulate_beta_posterior <- function(fit, n_samples) {
	
	if(is.null(fit$BETA)) {
		return(matrix(0, n_samples, 0))
	}
	
	beta_mean <- colMeans(fit$BETA)
	beta_cov <- cov(fit$BETA)

	if(requireNamespace("MASS", quietly = TRUE)) {
		beta_samples <- MASS::mvrnorm(n_samples, beta_mean, beta_cov)
	} else {
		p <- length(beta_mean)
		beta_samples <- matrix(NA, n_samples, p)
		for(i in 1:p) {
			beta_samples[, i] <- rnorm(n_samples, beta_mean[i], sqrt(beta_cov[i, i]))
		}
	}
	
	return(beta_samples)
}

####
# simulate_Y_posterior
####

#' Simulate Y from posterior predictive distribution
#' @noRd
simulate_Y_posterior <- function(fit, n_samples) {
	
	n <- nrow(fit$Y)
	family <- fit$family
	
	beta_samples <- simulate_posterior(fit, "beta", n_samples)

	if(fit$symmetric) {
		Y_samples <- array(NA, dim = c(n, n, n_samples))
	} else {
		Y_samples <- array(NA, dim = c(nrow(fit$Y), ncol(fit$Y), n_samples))
	}
	
	for(s in 1:n_samples) {
		EZ <- matrix(0, nrow(fit$Y), ncol(fit$Y))

		if(!is.null(fit$X) && ncol(beta_samples) > 0) {
			for(k in 1:dim(fit$X)[3]) {
				EZ <- EZ + beta_samples[s, k] * fit$X[,,k]
			}
		}
		
		if(!is.null(fit$APM) && !is.null(fit$BPM)) {
			ab_sim <- simulate_posterior(fit, "ab", 1)
			EZ <- EZ + outer(ab_sim$a[,1], ab_sim$b[,1], "+")
		}
		
		if(!is.null(fit$U_postmean) && !is.null(fit$V_postmean)) {
			UV_sim <- simulate_posterior(fit, "UV", 1)
			EZ <- EZ + UV_sim[,,1]
		}
		
		Y_samples[,,s] <- generate_Y_from_EZ(EZ, family, fit)
	}
	
	return(Y_samples)
}

####
# generate_Y_from_EZ
####

#' Generate Y from linear predictor EZ
#' @noRd
generate_Y_from_EZ <- function(EZ, family, fit) {
	
	n <- nrow(EZ)
	m <- ncol(EZ)
	
	Y <- switch(family,
		normal = EZ + matrix(rnorm(n * m), n, m),
		binary = matrix(rbinom(n * m, 1, pnorm(EZ)), n, m),
		poisson = matrix(rpois(n * m, exp(pmin(EZ, 10))), n, m),
		EZ
	)
	
	if(!is.null(fit$Y)) {
		Y[is.na(fit$Y)] <- NA
	}
	
	return(Y)
}

####
# posterior_quantiles
####

#' Extract posterior quantiles for model components
#'
#' @param fit Fitted ame model object
#' @param component Character; which component: "beta", "UV", "ab"
#' @param probs Numeric vector of probabilities for quantiles
#' 
#' @return Matrix or array of posterior quantiles
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
posterior_quantiles <- function(
	fit, 
	component = c("beta", "UV", "ab"),
	probs = c(0.025, 0.5, 0.975)
	){
	
	component <- match.arg(component)
	
	if(component == "beta") {
		if(is.null(fit$BETA)) {
			cli::cli_abort("No regression coefficients in model")
		}
		
		quantiles <- apply(fit$BETA, 2, quantile, probs = probs, na.rm = TRUE)
		rownames(quantiles) <- paste0("q", probs * 100)
		
		return(quantiles)
		
	} else if(component == "UV") {
		if(is.null(fit$U_postmean) || is.null(fit$V_postmean)) {
			cli::cli_abort("No multiplicative effects in model")
		}
		
		UV_samples <- simulate_posterior(fit, "UV", n_samples = 1000)
		n <- dim(UV_samples)[1]
		m <- dim(UV_samples)[2]
		
		UV_quantiles <- array(NA, dim = c(n, m, length(probs)))
		
		for(i in 1:n) {
			for(j in 1:m) {
				UV_quantiles[i, j, ] <- quantile(UV_samples[i, j, ], 
																				probs = probs, na.rm = TRUE)
			}
		}
		
		dimnames(UV_quantiles)[[3]] <- paste0("q", probs * 100)
		
		return(UV_quantiles)
		
	} else if(component == "ab") {
		if(is.null(fit$APM) || is.null(fit$BPM)) {
			cli::cli_abort("No additive effects in model")
		}
		
		ab_samples <- simulate_posterior(fit, "ab", n_samples = 1000)
		
		a_quantiles <- apply(ab_samples$a, 1, quantile, probs = probs, na.rm = TRUE)
		b_quantiles <- apply(ab_samples$b, 1, quantile, probs = probs, na.rm = TRUE)
		
		rownames(a_quantiles) <- paste0("q", probs * 100)
		rownames(b_quantiles) <- paste0("q", probs * 100)
		
		return(list(a = a_quantiles, b = b_quantiles))
	}
}