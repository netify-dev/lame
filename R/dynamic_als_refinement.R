# refine a stalled Gaussian bipartite fit with the objective and analytic
# gradient of its block solver; additive paths stay on their centring gauge
.lda_refine_bip_normal <- function(prep, state, settings, tol, max_iter) {
	s <- settings
	v <- state
	nr <- prep$nr
	nc <- prep$nc
	nt <- prep$Tt
	p <- nrow(v$coef_path)
	initial <- c(
		v$coef_path[, 1L],
		if (s$dynamic_ab) c(v$a_mat) else v$a_mat[, 1L],
		if (s$dynamic_ab) c(v$b_mat) else v$b_mat[, 1L],
		c(v$U),
		c(v$V),
		c(v$G_static)
	)
	unpack <- function(theta) {
		out <- v
		index <- 0L
		take <- function(n) {
			old <- index
			index <<- index + n
			theta[old + seq_len(n)]
		}
		out$coef_path[] <- take(p)
		a <- if (s$dynamic_ab) matrix(take(nr * nt), nr, nt) else
			matrix(rep(take(nr), nt), nr, nt)
		b <- if (s$dynamic_ab) matrix(take(nc * nt), nc, nt) else
			matrix(rep(take(nc), nt), nc, nt)
		out$a_mat <- sweep(a, 2L, colMeans(a), "-")
		out$b_mat <- sweep(b, 2L, colMeans(b), "-")
		out$U[] <- take(length(v$U))
		out[["V"]][] <- take(length(v$V))
		out$G_static[] <- take(length(v$G_static))
		out$Oarr <- .lda_uvg_to_o_array(
			out$U,
			out[["V"]],
			out$G_static,
			Tt = nt,
			time_names = prep$time_names
		)
		out$Omat <- apply(out$Oarr, c(1L, 2L), mean)
		out
	}
	weights <- array(as.numeric(is.finite(prep$Yarr)), dim(prep$Yarr))
	transition_u <- .lda_transition_from_presence(prep$row_presence)
	transition_v <- .lda_transition_from_presence(prep$col_presence)
	cache_theta <- cache_value <- cache_gradient <- cache_state <- cache_cert <- NULL
	evaluate <- function(theta) {
		if (identical(theta, cache_theta)) return(invisible(NULL))
		z <- unpack(theta)
		value <- .lda_objective(
			prep$Yarr,
			weights,
			prep,
			z$coef_path,
			z$a_mat,
			z$b_mat,
			z$Oarr,
			s$beta_dyn,
			s$lambda_beta,
			s$dynamic_beta_kind,
			s$rho_beta,
			s$dynamic_ab,
			s$lambda_ab,
			s$rho_ab,
			transition_a = transition_u,
			transition_b = transition_v,
			symmetric = FALSE,
			dynamic_uv = s$dynamic_uv,
			U = z$U,
			V = z$V,
			lambda_uv = s$lambda_uv,
			rho_uv = s$rho_uv,
			transition_u = transition_u,
			transition_v = transition_v,
			dynamic_uv_kind = s$dynamic_uv_kind,
			G_static = z$G_static,
			level_prec_ab = s$level_prec_ab,
			level_prec_uv = s$level_prec_uv,
			level_prec_g = s$level_prec_g
		)
		certificate <- .lda_stationarity(prep, z, s, tol)
		cache_theta <<- theta
		cache_value <<- value
		cache_gradient <<- unlist(certificate$gradient, use.names = FALSE)
		cache_state <<- z
		cache_cert <<- certificate
		invisible(NULL)
	}
	objective <- function(theta) {
		evaluate(theta)
		cache_value
	}
	gradient <- function(theta) {
		evaluate(theta)
		cache_gradient
	}
	evaluate(initial)
	initial_value <- cache_value
	information <- unlist(cache_cert$information, use.names = FALSE)
	optimum <- stats::optim(
		initial,
		objective,
		gradient,
		method = "L-BFGS-B",
		control = list(
			maxit = max_iter,
			factr = 10,
			pgtol = min(1e-5, sqrt(tol) / 10),
			parscale = 1 / sqrt(pmax(information, 1e-8)),
			lmm = 20L
		)
	)
	evaluate(optimum$par)
	accepted <- isTRUE(cache_cert$pass) && is.finite(cache_value) &&
		cache_value <= initial_value + 64 * .Machine$double.eps * max(1, abs(initial_value))
	list(
		state = cache_state,
		accepted = accepted,
		diagnostics = list(
			method = "L-BFGS-B",
			accepted = accepted,
			objective_before = initial_value,
			objective_after = cache_value,
			counts = optimum$counts,
			max_iter = max_iter,
			optimizer_convergence = optimum$convergence,
			max_score = cache_cert$max_score,
			threshold = cache_cert$threshold
		)
	)
}
