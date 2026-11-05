# refine a stalled directed count fit on its observed-data objective,
# rescaling curvature within the joint iteration budget
.lda_refine_directed_count <- function(prep, state, settings, tol, max_iter) {
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
		c(v$V)
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
		out$Oarr <- .lda_uv_to_o_array(out$U, out[["V"]])
		out$Omat <- apply(out$Oarr, c(1L, 2L), mean)
		out
	}
	unit_weights <- array(1, dim(prep$Yarr))
	transition_u <- .lda_transition_from_presence(prep$row_presence)
	transition_v <- .lda_transition_from_presence(prep$col_presence)
	cache_theta <- cache_value <- cache_gradient <- cache_state <- cache_cert <- NULL
	evaluate <- function(theta) {
		if (identical(theta, cache_theta)) return(invisible(NULL))
		z <- unpack(theta)
		eta <- .lda_eta_array(.lda_eta(
			prep,
			z$coef_path,
			z$a_mat,
			z$b_mat,
			z$Oarr,
			symmetric = FALSE,
			beta_dyn = s$beta_dyn
		))
		penalty <- .lda_objective(
			eta,
			unit_weights,
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
			level_prec_ab = s$level_prec_ab,
			level_prec_uv = s$level_prec_uv
		)
		value <- .ae_irls_deviance(prep$Yarr, eta, "poisson", "log") + penalty
		certificate <- .lda_stationarity(prep, z, s, tol)
		cache_theta <<- theta
		cache_value <<- value
		cache_gradient <<- unlist(certificate$gradient, use.names = FALSE)
		cache_state <<- z
		cache_cert <<- certificate
		invisible(NULL)
	}
	evaluate(initial)
	initial_value <- cache_value
	function_calls <- gradient_calls <- 0L
	phase_value <- initial_value
	no_worse <- function(value, reference = initial_value) is.finite(value) &&
		value <= reference + 64 * .Machine$double.eps * max(1, abs(reference))
	objective <- function(theta) {
		function_calls <<- function_calls + 1L
		evaluate(theta)
		if (isTRUE(cache_cert$pass) && no_worse(cache_value, phase_value)) {
			stop(structure(
				list(message = "model score reached", call = NULL, par = theta),
				class = c("lame_als_score_reached", "error", "condition")
			))
		}
		cache_value
	}
	gradient <- function(theta) {
		gradient_calls <<- gradient_calls + 1L
		evaluate(theta)
		cache_gradient
	}
	theta <- initial
	remaining <- as.integer(max_iter)
	cycles <- list()
	while (remaining > 0L) {
		evaluate(theta)
		information <- unlist(cache_cert$information, use.names = FALSE)
		budget <- min(50L, remaining)
		phase_value <- cache_value
		before_calls <- c("function" = function_calls, "gradient" = gradient_calls)
		optimum <- tryCatch(
			stats::optim(
				theta,
				objective,
				gradient,
				method = "L-BFGS-B",
				control = list(
					maxit = budget,
					factr = 10,
					pgtol = min(1e-5, sqrt(tol) / 10),
					parscale = 1 / sqrt(pmax(information, 1e-8)),
					lmm = 20L
				)
			),
			lame_als_score_reached = function(e) list(
				par = e$par,
				convergence = 0L,
				message = "model score reached"
			)
		)
		evaluate(optimum$par)
		step_ok <- no_worse(cache_value, phase_value)
		cycles[[length(cycles) + 1L]] <- list(
			budget = budget,
			objective_before = phase_value,
			objective_after = cache_value,
			max_score = cache_cert$max_score,
			accepted_step = step_ok,
			counts = c("function" = function_calls, "gradient" = gradient_calls) - before_calls,
			optimizer_convergence = optimum$convergence
		)
		remaining <- remaining - budget
		if (!step_ok) break
		theta <- optimum$par
		if (isTRUE(cache_cert$pass)) break
	}
	optimum$par <- theta
	optimum$counts <- c("function" = function_calls, "gradient" = gradient_calls)
	evaluate(optimum$par)
	accepted <- isTRUE(cache_cert$pass) && no_worse(cache_value)
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
			cycles = cycles,
			optimizer_convergence = optimum$convergence,
			termination = optimum$message,
			max_score = cache_cert$max_score,
			threshold = cache_cert$threshold
		)
	)
}
