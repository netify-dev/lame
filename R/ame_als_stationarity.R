# certify the likelihood equations at the fitted factor variances and the empirical variance fixed
# point separately (a flat working-response SSE establishes neither for a non-normal likelihood).
# `eta`, `xmat` (dyadic design as a cell-by-covariate matrix in array order), `working` and, for
# probit, `terms` may be passed by a caller that already has them; `ab_prior` (the EB actor-effect
# prior of ame_als_eb.R) adds that estimator's equations, including the EM fixed point of Sab
.ae_irls_stationarity <- function(
	prep,
	fit,
	family,
	link,
	R,
	symmetric,
	tol = 1e-6,
	eta = NULL,
	xmat = NULL,
	working = NULL,
	terms = NULL,
	ab_prior = NULL,
	ve_fixed = FALSE
) {
	nr <- prep$nr
	nc <- prep$nc
	nt <- prep$Tt
	p <- prep$p_dyad
	y <- prep$Yarr
	if (is.null(xmat) && p > 0L) {
		xmat <- matrix(aperm(prep$X, c(1L, 2L, 4L, 3L)), length(y), p)
	}
	if (p > 0L && !all(is.finite(xmat))) xmat[!is.finite(xmat)] <- 0
	if (is.null(eta)) {
		eta <- array(
			rep(fit$mu + outer(fit$a, fit$b, "+") + fit$Omat, nt),
			dim(y)
		)
		if (p > 0L) eta <- eta + as.vector(xmat %*% fit$beta)
	}
	# a count network under the EB prior: score and information of the Poisson-lognormal cell at the
	# iterate (ame_als_pln.R)
	pln <- !is.null(ab_prior) && identical(family, "poisson")
	# a normal network under the EB prior: score and information of the Gaussian cell (ame_als_gauss.R)
	gauss <- !is.null(ab_prior) && identical(family, "normal")
	if (pln || gauss) {
		if (is.null(terms)) terms <- if (pln) .ae_pln_terms(y, eta, fit$ve) else .ae_gauss_terms(y, eta, fit$ve)
		score <- terms$s
		information <- terms$w
	} else if (family == "poisson") {
		mu <- exp(eta)
		score <- y - mu
		information <- mu
	} else if (link == "probit") {
		sign_y <- 2 * y - 1
		if (is.null(terms)) terms <- .ae_probit_terms(eta, which(y == 1))
		score <- sign_y * exp(terms$ld - terms$outcome)
		information <- exp(2 * terms$ld - terms$lp1 - terms$lp0)
	} else {
		mu <- stats::plogis(eta)
		score <- y * stats::plogis(-eta) - (1 - y) * mu
		information <- mu * stats::plogis(-eta)
	}
	observed <- is.finite(y)
	score[!observed] <- 0
	information[!observed] <- 0
	threshold <- sqrt(tol)
	if (any(!is.finite(score)) || any(!is.finite(information))) {
		return(list(
			pass = FALSE,
			max_score = Inf,
			scale_change = Inf,
			unused_rank_score = Inf,
			threshold = threshold,
			reason = "nonfinite likelihood score"
		))
	}
	score_sum <- rowSums(score, dims = 2L)
	info_sum <- rowSums(information, dims = 2L)
	# working weights enter only through the factor-variance update
	weight_sum <- NULL
	if (R > 0L && symmetric && !is.null(ab_prior)) {
		# the symmetric EB factor scale update uses the factor rows' posterior second moment, which needs
		# the dyad weights
		wt <- if (is.null(working)) {
			if (pln) .ae_irls_pln(eta, terms)$w
			else if (gauss) terms$w
			else if (family == "poisson") .ae_irls_poisson(y, eta)$w
			else if (link == "probit") .ae_irls_binary_probit(y, eta, terms = terms)$w
			else .ae_irls_binary_logit(y, eta)$w
		} else working$w
		wt[!observed] <- 0
		weight_sum <- rowSums(wt, dims = 2L)
	}
	if (R > 0L && !symmetric) {
		if (is.null(working)) {
			working <- if (pln) {
				.ae_irls_pln(eta, terms)
			} else if (gauss) {
				list(z = y, w = terms$w)
			} else if (family == "poisson") {
				.ae_irls_poisson(y, eta)
			} else if (link == "probit") {
				.ae_irls_binary_probit(y, eta, terms = terms)
			} else {
				.ae_irls_binary_logit(y, eta)
			}
		}
		wt <- working$w
		wt[!observed] <- 0
		weight_sum <- rowSums(wt, dims = 2L)
	}
	beta_gradient <- beta_information <- numeric(p)
	if (p > 0L) {
		beta_gradient <- -2 * as.vector(crossprod(xmat, as.vector(score)))
		beta_information <- 2 * as.vector(crossprod(xmat^2, as.vector(information)))
	}
	gradient <- c(-2 * sum(score_sum), beta_gradient)
	curvature <- c(2 * sum(info_sum), beta_information)
	if (symmetric) {
		gradient <- c(gradient, -2 * (rowSums(score_sum) + colSums(score_sum)))
		curvature <- c(curvature, 2 * (rowSums(info_sum) + colSums(info_sum)))
	} else {
		gradient <- c(gradient, -2 * rowSums(score_sum), -2 * colSums(score_sum))
		curvature <- c(curvature, 2 * rowSums(info_sum), 2 * colSums(info_sum))
	}
	eb_extra <- NULL
	eb_img <- NULL
	if (!is.null(ab_prior)) {
		# Laplace correction and EM image of the covariance at the iterate's own weights and predictor
		w_eb <- if (!is.null(working)) working$w else if (pln) .ae_irls_pln(eta, terms)$w else if (gauss) terms$w else if (identical(link, "probit")) {
			.ae_irls_binary_probit(y, eta, terms = terms)$w
		} else {
			.ae_irls_binary_logit(y, eta)$w
		}
		eb_img <- tryCatch(
			.ae_eb_at_iterate(prep, fit, w_eb, eta, terms, link, ab_prior),
			error = function(e) NULL
		)
		if (is.null(eb_img)) {
			return(list(
				pass = FALSE,
				max_score = Inf,
				scale_change = Inf,
				unused_rank_score = Inf,
				threshold = threshold,
				reason = "empirical-bayes correction not computable"
			))
		}
		# a symmetric fit's penalty and certificate are on the ordered-cell scale, counting each dyad twice
		dup <- 1 / ab_prior$nu
		Lam <- .ae_eb_precision(fit$sab, ab_prior)
		M <- .ae_eb_mean_design(ab_prior)
		al <- if (symmetric) fit$a else c(fit$a, fit$b)
		if (!is.null(M) && length(fit$gamma) == ncol(M)) {
			al <- al - as.vector(M %*% fit$gamma)
		}
		pos <- 1L + p + seq_len(ab_prior$q)
		gradient[pos] <- gradient[pos] + 2 * dup * as.vector(Lam %*% al)
		curvature[pos] <- curvature[pos] + 2 * dup * diag(Lam)
		# Laplace-marginal score of (mu, beta)
		corr <- eb_img$corr_th
		if (length(corr) == p + 1L) {
			gradient[seq_len(p + 1L)] <- gradient[seq_len(p + 1L)] + 2 * corr
		}
		if (!is.null(M) && length(fit$gamma) == ncol(M)) {
			corr_g <- if (length(eb_img$corr_g) == ncol(M)) eb_img$corr_g else numeric(ncol(M))
			eb_extra <- list(
				gradient = dup * (-2 * as.vector(crossprod(M, Lam %*% al)) + 2 * corr_g),
				curvature = dup * 2 * pmax(diag(crossprod(M, Lam %*% M)), 0)
			)
		}
	}
	# the factors take the same Laplace cell correction as the coefficients
	score_sum_f <- score_sum
	if (!is.null(eb_img) && identical(dim(eb_img$kd), dim(score_sum))) {
		score_sum_f <- score_sum - eb_img$kd
	}
	scale_change <- unused_rank_score <- 0
	suv0 <- .ae_uv_prior_scale(y, family)
	kappa0 <- 2 * R + 2
	# a directed unipartite EB fit uses the factor prior of ame() and lame() (R/uv_prior.R), as in
	# .ae_irls_run
	if (R > 0L && !symmetric && !isTRUE(prep$bip) && !is.null(ab_prior)) {
		uvp <- .lame_uv_prior(y, family, R)
		kappa0 <- uvp$df
		suv0 <- uvp$psi
	}
	if (R > 0L) {
		if (symmetric) {
			factor <- sweep(fit$U, 2, sqrt(abs(fit$L)), "*")
			signed_factor <- sweep(factor, 2, sign(fit$L), "*")
			sig2 <- fit$sig2_u[1]
			gradient <- c(
				gradient,
				c(
					-2 * (score_sum_f + t(score_sum_f)) %*%
						signed_factor + 4 * factor / sig2
				)
			)
			curvature <- c(
				curvature,
				c(
					2 * (info_sum + t(info_sum)) %*%
						(signed_factor^2) + 4 / sig2
				)
			)
			next_scale <- if (is.null(ab_prior)) {
				(kappa0 * suv0 + 2 * sum(abs(fit$L))) / (kappa0 + 2 * nr * R)
			} else {
				(suv0 + sum(abs(fit$L)) + .ae_sym_trace(factor, weight_sum, sig2)) /
					(nr * R + kappa0)
			}
			scale_change <- abs(log(next_scale / sig2))

			# a zero eigenvalue also needs the unused-rank optimality check (its square-root-coordinate
			# gradient is zero)
			active <- abs(fit$L) > sqrt(.Machine$double.eps) * max(1, abs(fit$L))
			if (any(!active)) {
				projection <- diag(nr)
				if (any(active)) {
					projection <- projection - tcrossprod(fit$U[, active, drop = FALSE])
				}
				residual_score <- projection %*% ((score_sum_f + t(score_sum_f)) / 2) %*%
					projection
				eigen_score <- max(abs(
					eigen(
						residual_score,
						symmetric = TRUE,
						only.values = TRUE
					)$values
				))
				unused_rank_score <- max(0, eigen_score - 1 / sig2) / max(1, 1 / sig2)
			}
		} else {
			sig2_u <- fit$sig2_u
			sig2_v <- fit$sig2_v
			gradient <- c(
				gradient,
				c(-2 * score_sum_f %*% fit[["V"]] + 2 * sweep(fit$U, 2, sig2_u, "/")),
				c(-2 * t(score_sum_f) %*% fit$U + 2 * sweep(fit[["V"]], 2, sig2_v, "/"))
			)
			curvature <- c(
				curvature,
				c(2 * info_sum %*% (fit[["V"]]^2) + matrix(rep(2 / sig2_u, each = nr), nr, R)),
				c(2 * t(info_sum) %*% (fit$U^2) + matrix(rep(2 / sig2_v, each = nc), nc, R))
			)
			# the variance update's posterior second moments, leaving the supplied factors unchanged
			next_variance <- als_factor_update_cpp(
				matrix(0, nr, nc),
				weight_sum,
				fit$U,
				fit[["V"]],
				sig2_u,
				sig2_v,
				suv0,
				kappa0,
				1L,
				1e-6
			)
			scale_change <- max(abs(log(c(
				next_variance$sig2_u / sig2_u,
				next_variance$sig2_v / sig2_v
			))))
		}
	}
	# EM fixed point of the actor-effect covariance, at the iterate itself
	if (!is.null(ab_prior)) {
		old_s <- fit$sab; new_s <- eb_img$sab_new
		sab_gap <- max(
			abs(log(diag(new_s) / diag(old_s))),
			if (identical(ab_prior$kind, "directed"))
				abs(new_s[1, 2] - old_s[1, 2]) / sqrt(prod(diag(old_s))) else 0
		)
		scale_change <- if (is.finite(sab_gap)) max(scale_change, sab_gap) else Inf
		# EM fixed point of the noise variance of a count or normal network (a variance held fixed has
		# none to certify)
		if ((pln || gauss) && !isTRUE(ve_fixed)) {
			ve_gap <- abs(
				log(
					.ae_pln_ve_image(
						sum(terms$ee),
						sum(observed),
						ab_prior$nu,
						ab_prior$ve_prior,
						.ae_ve_floor(ab_prior)
					) / fit$ve
				)
			)
			scale_change <- if (is.finite(ve_gap)) max(scale_change, ve_gap) else Inf
		}
		if (!is.null(eb_extra)) {
			gradient <- c(gradient, eb_extra$gradient)
			curvature <- c(curvature, eb_extra$curvature)
		}
	}
	scaled_score <- abs(gradient) / sqrt(pmax(curvature, .Machine$double.xmin))
	max_score <- max(scaled_score)
	list(
		pass = all(is.finite(scaled_score)) && max_score <= threshold &&
			is.finite(scale_change) && scale_change <= threshold &&
			unused_rank_score <= threshold,
		max_score = max_score,
		scale_change = scale_change,
		unused_rank_score = unused_rank_score,
		threshold = threshold,
		score = scaled_score,
		gradient = gradient,
		information = curvature
	)
}
