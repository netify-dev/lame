# observed-probit replica exchange for static symmetric factors and AR actors
.lame_tempered_integer <- function(x, name, lower = 1L) {
	if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
			x < lower || x != floor(x) || x > .Machine$integer.max)
		cli::cli_abort("{.arg {name}} must be an integer at least {lower}.")
	as.integer(x)
}

.lame_tempered_input <- function(Y, Xdyad, R) {
	if (is.array(Y) && length(dim(Y)) == 3L) {
		y <- Y
		Y <- lapply(seq_len(dim(y)[3L]), function(t) y[, , t])
		names(Y) <- dimnames(y)[[3L]]
	}
	if (is.array(Xdyad) && length(dim(Xdyad)) == 4L) {
		x <- Xdyad
		Xdyad <- lapply(seq_len(dim(x)[4L]), function(t) {
			y <- x[, , , t, drop = FALSE]
			dim(y) <- dim(x)[1:3]
			dimnames(y) <- dimnames(x)[1:3]
			y
		})
		names(Xdyad) <- dimnames(x)[[4L]]
	}
	input <- .lame_robust_input(Y, Xdyad, R)
	if (nrow(input$Y[[1L]]) < 4L || input$R < 1L)
		cli::cli_abort("The tempered sampler requires at least four actors and positive rank.")
	for (y in input$Y) {
		if (!isTRUE(all.equal(y, t(y), check.attributes = FALSE)) ||
				any(!y[is.finite(y)] %in% c(0, 1)))
			cli::cli_abort("Tempered binary outcomes must be symmetric, including missingness, and contain only zero or one.")
	}
	for (x in input$Xdyad) for (j in seq_len(dim(x)[3L])) {
		if (!isTRUE(all.equal(x[, , j], t(x[, , j]), check.attributes = FALSE)))
			cli::cli_abort("The tempered sampler requires symmetric dyadic covariates.")
	}
	input
}

.lame_tempered_prediction_design <- function(object, X) {
	if (!identical(object[["sampler"]], "tempered") || is.null(X)) return(X)
	ok <- is.numeric(X) && !is.complex(X) && length(dim(X)) == 3L &&
		dim(X)[1L] == dim(X)[2L]
	if (ok) {
		rn <- dimnames(X)[[1L]]; cn <- dimnames(X)[[2L]]
		ok <- is.null(rn) || is.null(cn) || identical(rn, cn)
		for (j in seq_len(dim(X)[3L])) {
			x <- X[, , j]
			ok <- ok && all(is.finite(x[row(x) != col(x)])) &&
				isTRUE(all.equal(x, t(x), check.attributes = FALSE))
		}
	}
	if (!ok) cli::cli_abort("Tempered predictions require symmetric finite dyadic covariates on the same row and column actors.")
	X
}

.lame_tempered_prior <- function(prior, Y) {
	allowed <- c(
		"uv_prior", "Sab0", "eta0", "etaab", "rho_ab_mean", "rho_ab_sd",
		"sigma_ab_shape", "sigma_ab_scale", "ab_min_observed", "ab_var_cap"
	)
	if (!is.list(prior) || (length(prior) && (is.null(names(prior)) ||
			anyNA(names(prior)) || anyDuplicated(names(prior)) ||
			any(!names(prior) %in% allowed) || any(vapply(prior, is.null, logical(1))))))
		cli::cli_abort("For the tempered sampler, {.arg prior} must be a named list of documented tempered prior settings.")
	prior$uv_prior <- prior$uv_prior %||% "hierarchical"
	if (!is.character(prior$uv_prior) || length(prior$uv_prior) != 1L ||
			is.na(prior$uv_prior) || !prior$uv_prior %in% c("unit", "hierarchical"))
		cli::cli_abort("{.code prior$uv_prior} must be {.val unit} or {.val hierarchical}.")
	# same first-period empirical prior as the other binary path
	y <- Y[[1L]]
	ydist <- table(y)
	ymode <- if (length(ydist)) as.numeric(names(ydist)[ydist == max(ydist)])[1L] else 0
	yb <- 1 * (y != ymode)
	ybar <- mean(yb, na.rm = TRUE)
	vscale <- 1
	eta <- round(4 + 3 * nrow(y) / 100)
	if (is.finite(ybar) && ybar > 0 && ybar < 1) {
		e <- (yb - ybar) / dnorm(qnorm(ybar))
		diag(e) <- 0
		a <- rowMeans(e, na.rm = TRUE); a[is.na(a)] <- 0
		b <- colMeans(e, na.rm = TRUE); b[is.na(b)] <- 0
		vscale <- mean(c(stats::var(a), stats::var(b)))
		if (!is.finite(vscale) || vscale <= 0) vscale <- 1
		phat <- pnorm(qnorm(ybar) + outer(a, b, "+"))
		eta <- round(4 * 0.25 / mean(phat * (1 - phat)))
	}
	prior$Sab0 <- prior$Sab0 %||% (diag(2L) * vscale)
	sab <- prior$Sab0
	if (!is.matrix(sab) || !is.numeric(sab) || is.complex(sab) ||
			!identical(dim(sab), c(2L, 2L)) || any(!is.finite(sab)) ||
			any(diag(sab) <= 0) || any(sab[row(sab) != col(sab)] != 0) ||
			sab[1L, 1L] != sab[2L, 2L])
		cli::cli_abort("For shared symmetric actors, {.code prior$Sab0} must be a positive multiple of the two by two identity matrix.")
	prior$eta0 <- prior$eta0 %||% eta
	prior$etaab <- prior$etaab %||% prior$eta0
	prior$rho_ab_mean <- prior$rho_ab_mean %||% 0.8
	prior$rho_ab_sd <- prior$rho_ab_sd %||% 0.15
	prior$sigma_ab_shape <- prior$sigma_ab_shape %||% 2
	for (nm in c("eta0", "etaab", "rho_ab_sd", "sigma_ab_shape")) {
		x <- prior[[nm]]
		if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) || x <= 0)
			cli::cli_abort("{.code prior${nm}} must be a finite positive number.")
	}
	rho <- prior$rho_ab_mean
	if (!is.numeric(rho) || is.complex(rho) || length(rho) != 1L ||
			!is.finite(rho) || abs(rho) >= 0.99)
		cli::cli_abort("{.code prior$rho_ab_mean} must lie strictly between -0.99 and 0.99.")
	prior$sigma_ab_scale <- prior$sigma_ab_scale %||%
		.lame_sigma_ab_prior_scale(prior$Sab0, prior$etaab, rho, TRUE, TRUE)
	x <- prior$sigma_ab_scale
	if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) || x <= 0)
		cli::cli_abort("{.code prior$sigma_ab_scale} must be a finite positive number.")
	prior$ab_min_observed <- .lame_tempered_integer(
		prior$ab_min_observed %||% max(10L, length(Y)),
		"prior$ab_min_observed",
		0L
	)
	prior$ab_var_cap <- prior$ab_var_cap %||% Inf
	if (!identical(as.double(prior$ab_var_cap), Inf))
		cli::cli_abort("The tempered sampler currently requires {.code prior$ab_var_cap = Inf}.")
	prior
}

.lame_tempered_design <- function(input, g) {
	Y <- input$Y; Xdyad <- input$Xdyad
	n <- nrow(Y[[1L]]); nt <- length(Y)
	y <- array(unlist(Y), c(n, n, nt), dimnames = c(dimnames(Y[[1L]]), list(names(Y))))
	pd <- if (is.null(Xdyad)) 0L else dim(Xdyad[[1L]])[3L]
	x <- if (is.null(Xdyad)) NULL else array(
		unlist(Xdyad),
		c(n, n, pd, nt),
		dimnames = c(dimnames(Xdyad[[1L]]), list(names(Y)))
	)
	design <- get_design_rep(y, x, NULL, NULL, rownames(Y[[1L]]), TRUE, n, nt, pd, 0L, 0L)
	if (identical(g, NA)) g <- NA_real_
	if (!is.numeric(g) || is.complex(g) || length(g) != 1L ||
			(!is.na(g) && (!is.finite(g) || g <= 0)))
		cli::cli_abort("For the tempered sampler, {.arg g} must be NA or a finite positive scalar.")
	if (is.na(g)) g <- n^2
	Q <- apply(design$xxLong, c(1L, 2L), sum) / (nt * g)
	if (any(!is.finite(Q)) || inherits(try(chol(Q), silent = TRUE), "try-error"))
		cli::cli_abort("The tempered sampler requires a full-rank coefficient design for its proper Gaussian prior.")
	list(Y = y, Xlist = design$Xlist, Q = Q, g = g)
}

.lame_tempered_context <- function(Y, X, rank, Q, prior) {
	n <- dim(Y)[1L]; nt <- dim(Y)[3L]; p <- dim(X[[1L]])[3L]
	ij <- do.call(rbind, lapply(seq_len(nt), function(t) {
		a <- which(upper.tri(Y[, , t]) & is.finite(Y[, , t]), arr.ind = TRUE)
		cbind(a, rep.int(t, nrow(a)))
	}))
	x <- do.call(rbind, lapply(seq_len(nt), function(t) {
		matrix(X[[t]], n * n, p)[which(upper.tri(Y[, , t]) & is.finite(Y[, , t])), , drop = FALSE]
	}))
	list(
		n = n,
		tt = nt,
		r = rank,
		p = p,
		H = qr.Q(qr(stats::contr.helmert(n))),
		B = solve(chol(Q)),
		X = x,
		ii = as.integer(ij[, 1L] - 1L),
		jj = as.integer(ij[, 2L] - 1L),
		time = as.integer(ij[, 3L] - 1L),
		sy = 2 * Y[ij] - 1,
		av = prior$etaab / 2,
		bv = prior$etaab * prior$Sab0[1L, 1L] / 2,
		as = prior$sigma_ab_shape,
		bs = prior$sigma_ab_scale,
		rm = prior$rho_ab_mean,
		rs = prior$rho_ab_sd,
		hierarchical = FALSE,
		collapsed = identical(prior$uv_prior, "hierarchical")
	)
}

.lame_tempered_decode <- function(q, cx, Q) {
	d <- cx$n * cx$r + cx$r + cx$p + (cx$n - 1L) * (cx$tt + 1L)
	va <- exp(q[d + 1L]); sigma <- exp(q[d + 2L]); rho <- 0.99 * tanh(q[d + 3L])
	k <- 0L
	take <- function(n) { i <- k + seq_len(n); k <<- k + n; q[i] }
	U <- matrix(take(cx$n * cx$r), cx$n, cx$r)
	lambda <- take(cx$r)
	beta <- drop(backsolve(chol(Q), take(cx$p)))
	level <- drop(cx$H %*% take(cx$n - 1L)) * sqrt(va)
	innov <- matrix(take((cx$n - 1L) * cx$tt), cx$n - 1L, cx$tt)
	dev <- innov
	dev[, 1L] <- sigma / sqrt(1 - rho^2) * innov[, 1L]
	for (t in 2:cx$tt) dev[, t] <- rho * dev[, t - 1L] + sigma * innov[, t]
	list(
		U = U,
		V = sweep(U, 2L, lambda, "*"),
		lambda = lambda,
		beta = beta,
		level = level,
		path = level + cx$H %*% dev,
		va = va,
		sigma = sigma,
		rho = rho
	)
}

.lame_tempered_conditional_scale <- function(theta) {
	b <- 0.5
	repeat {
		if (abs(theta) <= b) {
			x <- rgamma(1L, shape = 2, rate = b)
			lr <- -0.5 * (theta * x)^2
		} else {
			x <- sqrt(-2 * log(runif(1L))) / abs(theta)
			lr <- -b * x
		}
		if (log(runif(1L)) < lr) return(1 / x)
	}
}

.lame_tempered_storage <- function(options) {
	if (is.null(options)) options <- list()
	allowed <- c("save_UV", "save_UV_draws", "save_ab", "thin_UV", "thin_ab")
	if (!is.list(options) || (length(options) && (is.null(names(options)) ||
			anyNA(names(options)) || anyDuplicated(names(options)) || any(!names(options) %in% allowed))))
		cli::cli_abort("{.arg posterior_opts} must be a named list of documented storage options.")
	for (nm in intersect(names(options), c("save_UV", "save_UV_draws", "save_ab"))) {
		if (!is.logical(options[[nm]]) || length(options[[nm]]) != 1L || is.na(options[[nm]]))
			cli::cli_abort("{.code posterior_opts${nm}} must be TRUE or FALSE.")
	}
	for (nm in c("thin_UV", "thin_ab"))
		options[[nm]] <- .lame_tempered_integer(options[[nm]] %||% 1L, paste0("posterior_opts$", nm))
	options$save_UV_draws <- isTRUE(options$save_UV_draws) || isTRUE(options$save_UV)
	options$save_ab <- isTRUE(options$save_ab)
	options
}

.lame_tempered_output <- function(a, input, design, options, save_log_lik = FALSE, gof = FALSE) {
	cx <- a$context
	n <- cx$n; nt <- cx$tt; S <- ncol(a$q)
	Y <- input$Y; X <- design$Xlist
	actors <- rownames(Y[[1L]]); periods <- names(Y)
	f <- list(
		family = "binary",
		symmetric = TRUE,
		mode = "unipartite",
		n_time = nt,
		dynamic_ab = TRUE,
		dynamic_uv = FALSE,
		dynamic_beta = FALSE,
		dynamic_rho = FALSE,
		intercept = TRUE,
		rvar = TRUE,
		cvar = TRUE,
		dcor = FALSE,
		R = cx$r,
		g = design$g,
		prior = a$prior,
		Y = Y,
		Xlist = X,
		method = "mcmc",
		sampler = "tempered",
		likelihood_unit = "unordered_dyad_period",
		observed_dyads = length(cx$sy)
	)
	f$BETA <- matrix(0, S, cx$p, dimnames = list(NULL, dimnames(X[[1L]])[[3L]]))
	f$VC <- matrix(0, S, 2L, dimnames = list(NULL, c("va", "ve")))
	f$tempering <- list(
		temperatures = a$temperatures,
		stepsizes = a$stepsizes,
		counts = a$counts,
		warmcounts = a$warmcounts,
		modes = a$modes,
		kernel = a$kernel,
		start_iteration = a$start_iteration,
		end_iteration = a$end_iteration
	)
	bad <- sum(a$counts$nonfinite) + sum(a$warmcounts$nonfinite)
	f$sampler_health <- data.frame(
		step = "tempered",
		numerical_failures = bad,
		prior_bound_rejections = 0L
	)
	f$sampler_health_iterations <- a$end_iteration - a$start_iteration
	f$tryErrorChecks <- f$mh_counters <- list(tempered = bad)
	f$refusal_share <- c(UV = 0, G = 0)
	if (S == 0L) {
		f$no_draws_stored <- TRUE
		class(f) <- c("lame", "ame")
		return(f)
	}
	f$VC[, "ve"] <- 1
	f$a_draws <- array(0, c(n, nt, S), dimnames = list(actors, periods, NULL))
	f$a_level_samples <- matrix(0, n, S, dimnames = list(actors, NULL))
	f$U_samples <- f$V_samples <- array(
		0,
		c(n, cx$r, S),
		dimnames = list(actors, paste0("factor", seq_len(cx$r)), NULL)
	)
	f$rho_ab <- f$sigma_ab <- numeric(S)
	product_mean <- matrix(0, n, n, dimnames = list(actors, actors))
	ez <- lapply(Y, function(y) matrix(0, n, n, dimnames = dimnames(y)))
	ypm <- ez
	obs <- cbind(row = cx$ii + 1L, col = cx$jj + 1L, time = cx$time + 1L)
	if (isTRUE(save_log_lik)) {
		f$log_lik <- matrix(0, S, nrow(obs))
		f$log_lik_index <- obs
	}
	gof_values <- NULL
	if (isTRUE(gof)) {
		observed <- vapply(Y, gof_stats, numeric(5L), mode = "unipartite")
		gof_values <- array(
			NA_real_,
			c(5L, nt, S + 1L),
			dimnames = list(rownames(observed), periods, c("obs", as.character(seq_len(S))))
		)
		gof_values[, , 1L] <- observed
	}
	for (s in seq_len(S)) {
		z <- .lame_tempered_decode(a$q[, s], cx, a$Q)
		if (isTRUE(cx$collapsed)) {
			vscale <- if (is.null(a$factor_variances)) {
				vapply(z$lambda, .lame_tempered_conditional_scale, numeric(1))
			} else a$factor_variances[, s]
			z$U <- sweep(z$U, 2L, sqrt(vscale), "*")
			z$lambda <- z$lambda / vscale
		}
		mu <- colMeans(z$U)
		u <- sweep(z$U, 2L, mu)
		shift <- drop(u %*% (z$lambda * mu))
		beta <- z$beta
		beta[1L] <- beta[1L] + sum(z$lambda * mu^2)
		path <- z$path + shift
		v <- sweep(u, 2L, z$lambda, "*")
		f$BETA[s, ] <- beta
		f$VC[s, "va"] <- z$va
		f$a_draws[, , s] <- path
		f$a_level_samples[, s] <- z$level + shift
		f$U_samples[, , s] <- u
		f$V_samples[, , s] <- v
		f$rho_ab[s] <- z$rho
		f$sigma_ab[s] <- z$sigma
		product <- tcrossprod(u, v)
		product_mean <- product_mean + product / S
		for (t in seq_len(nt)) {
			eta <- matrix(matrix(X[[t]], n * n, cx$p) %*% beta, n, n) +
				outer(path[, t], path[, t], "+") + product
			ez[[t]] <- ez[[t]] + eta / S
			ypm[[t]] <- ypm[[t]] + pnorm(eta) / S
			if (isTRUE(save_log_lik)) {
				rows <- which(obs[, "time"] == t)
				cells <- obs[rows, c("row", "col"), drop = FALSE]
				f$log_lik[s, rows] <- pnorm(cx$sy[rows] * eta[cells], log.p = TRUE)
			}
			if (isTRUE(gof)) {
				y <- matrix(0, n, n, dimnames = dimnames(Y[[t]]))
				y[upper.tri(y)] <- stats::rbinom(n * (n - 1L) / 2L, 1L, pnorm(eta[upper.tri(eta)]))
				y <- y + t(y)
				y[is.na(Y[[t]])] <- NA_real_
				gof_values[, t, s + 1L] <- gof_stats(y, mode = "unipartite")
			}
		}
	}
	f$b_draws <- f$a_draws
	f$b_level_samples <- f$a_level_samples
	f$a_dynamic <- f$b_dynamic <- apply(f$a_draws, c(1L, 2L), mean)
	f$a_dynamic_sd <- f$b_dynamic_sd <- apply(f$a_draws, c(1L, 2L), stats::sd)
	f$APM <- f$BPM <- rowMeans(f$a_level_samples)
	f$a_level <- f$b_level <- f$APM
	f$ab_draw_index <- f$uv_draw_index <- seq_len(S)
	f$EZ <- ez; f$YPM <- ypm; f$ULUPM <- product_mean
	ev <- eigen(product_mean, symmetric = TRUE)
	take <- order(abs(ev$values), decreasing = TRUE)[seq_len(cx$r)]
	f$U <- ev$vectors[, take, drop = FALSE]
	rownames(f$U) <- actors
	f$L <- ev$values[take]
	f[["V"]] <- sweep(f$U, 2L, f$L, "*")
	f$forecast_state <- list(
		a = matrix(f$a_draws[, nt, ], n, S, dimnames = list(actors, NULL)),
		b = matrix(f$b_draws[, nt, ], n, S, dimnames = list(actors, NULL)),
		a_level = f$a_level_samples,
		b_level = f$b_level_samples,
		U = f$U_samples,
		V = f$V_samples,
		draw_index = seq_len(S)
	)
	if (!isTRUE(options$save_ab)) {
		f$a_draws <- f$b_draws <- f$a_level_samples <- f$b_level_samples <- NULL
		f$ab_draw_index <- NULL
	} else {
		idx <- seq.int(1L, S, by = options$thin_ab)
		f$a_draws <- f$a_draws[, , idx, drop = FALSE]
		f$b_draws <- f$b_draws[, , idx, drop = FALSE]
		f$a_level_samples <- f$a_level_samples[, idx, drop = FALSE]
		f$b_level_samples <- f$b_level_samples[, idx, drop = FALSE]
		f$ab_draw_index <- idx
	}
	# as in lame()'s other MCMC path, latent draws align with every coefficient draw
	if (!isTRUE(options$save_UV_draws)) {
		f$U_samples <- f$V_samples <- NULL
		f$uv_draw_index <- NULL
	}
	if (!is.null(gof_values)) {
		f$GOF <- lapply(
			seq_len(5L),
			function(k) matrix(
				gof_values[k, , ],
				nt,
				S + 1L,
				dimnames = dimnames(gof_values)[2:3]
			)
		)
		names(f$GOF) <- dimnames(gof_values)[[1L]]
	}
	f$likelihood_unit <- "unordered_dyad_period"
	f$observed_dyads <- nrow(obs)
	f$log_lik_method <- "observed_exact"
	f$coefficient_prior <- list(
		kind = "panel_design",
		description = "Regression coefficients have zero-mean Gaussian prior precision sum_t(X_t'X_t) / (T * g), independent of residual variance."
	)
	f$static_uv_prior <- list(
		distribution = a$prior$uv_prior,
		description = if (isTRUE(cx$collapsed)) {
			"Signed weights are independent N(0,1); each coordinate column has variance v ~ inverse-gamma(1,0.5)."
		} else "Signed weights and raw latent coordinates are independent N(0,1); coordinate variance is fixed at one."
	)
	f$sigma_ab_prior_scale <- a$prior$sigma_ab_scale
	f$dynamic <- list(uv = FALSE, ab = TRUE, beta = FALSE, G = FALSE)
	class(f) <- c("lame", "ame")
	f
}

.lame_tempered_temperatures <- function(kind) {
	# preset temperature ladders; the two proper priors need different ones
	if (identical(kind, "unit")) {
		c(
			1, 0.83293052762751074, 0.69376546452775989, 0.57727773755775169,
			0.48020608894183497, 0.39601412666417241, 0.3255210900759884,
			0.26257388413370109, 0.21303530620497704, 0.174050371129475,
			0.14488550518048893, 0.12179415390241212, 0.1005882595034364,
			0.081939200184798736, 0.059297344944402877, 0.038204922023243111,
			0.021535445555524498, 0.01011524252308655, 0.0037440366708393043,
			0.0010000000000000002, 0
		)
	} else {
		c(
			1, 0.81660208791927114, 0.67577955428889891, 0.58878682375107649,
			0.51299261959422671, 0.4359242437035557, 0.36343784534420809,
			0.29795827489591942, 0.24117643806532418, 0.19966454733678568,
			0.16596438284466739, 0.14069545360295471, 0.11958901945740609,
			0.095852344808021156, 0.07194199271292033, 0.045800843055654836,
			0.024973254045624158, 0.011817315149644198, 0.0040257258030339579,
			0.0010000000000000002, 0
		)
	}
}

.lame_tempered_rng <- function(seed, fun) {
	if (!is.null(seed)) seed <- .lame_tempered_integer(seed, "seed", 0L)
	had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	if (had) old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	kind <- RNGkind()
	if (kind[2L] == "Box-Muller" || any(kind == "user-supplied"))
		cli::cli_abort("Exact tempered checkpoints require a fully serializable RNG; use a built-in generator with the Inversion normal method.")
	on.exit({
		do.call(RNGkind, as.list(kind))
		if (had) assign(".Random.seed", old, envir = .GlobalEnv)
		else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
			rm(".Random.seed", envir = .GlobalEnv)
	}, add = TRUE)
	set.seed(seed)
	fun()
}

.lame_tempered_checkpoint_write <- function(checkpoint, path) {
	if (is.null(path)) return(invisible(NULL))
	tmp <- tempfile(pattern = ".lame-checkpoint-", tmpdir = dirname(path))
	on.exit(unlink(tmp), add = TRUE)
	saveRDS(checkpoint, tmp)
	# renaming over a file is atomic on unix; on windows keep the last valid
	# checkpoint in a sibling backup while replacing it
	if (.Platform$OS.type != "windows" || !file.exists(path)) {
		if (!file.rename(tmp, path)) cli::cli_abort("Could not save checkpoint {.file {path}}.")
	} else {
		backup <- paste0(path, ".previous")
		if (file.exists(backup)) unlink(backup)
		if (!file.rename(path, backup)) cli::cli_abort("Could not preserve checkpoint {.file {path}}.")
		if (!file.rename(tmp, path)) {
			file.rename(backup, path)
			cli::cli_abort("Could not replace checkpoint {.file {path}}; its previous version was retained.")
		}
		unlink(backup)
	}
	invisible(path)
}

.lame_tempered_validate_rng <- function(seed, kind) {
	generators <- c(
		"Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
		"Mersenne-Twister", "Knuth-TAOCP", "Knuth-TAOCP-2002", "L'Ecuyer-CMRG"
	)
	codes <- c(0L, 1L, 2L, 3L, 4L, 6L, 7L)
	lengths <- c(4L, 3L, 3L, 626L, 102L, 102L, 7L)
	normals <- c("Buggy Kinderman-Ramage", "Ahrens-Dieter", "Inversion", "Kinderman-Ramage")
	normal_codes <- c(0L, 1L, 4L, 5L)
	# an RNGkind() with binom.kind serialises a fourth generator; a three-element
	# kind encodes the legacy binomial generator as zero
	kind_lengths <- if ("binom.kind" %in% names(formals(RNGkind))) c(3L, 4L) else 3L
	ok <- is.character(kind) && length(kind) %in% kind_lengths && !anyNA(kind)
	if (ok) {
		g <- match(kind[1L], generators); n <- match(kind[2L], normals)
		s <- match(kind[3L], c("Rounding", "Rejection"))
		b <- if (length(kind) == 4L) match(kind[4L], c("Buggy BTPE", "BTPE")) - 1L else 0L
		ok <- !anyNA(c(g, n, s, b))
		# only the kind code is checked for NA: the remaining signed integers hold
		# raw generator bits, and INT_MIN among them reads as NA in R
		if (ok) ok <- is.integer(seed) && length(seed) == lengths[g] &&
			!is.na(seed[1L]) && seed[1L] == codes[g] + 100L * normal_codes[n] +
			10000L * (s - 1L) + 100000L * b
	}
	if (!ok) cli::cli_abort("The tempered checkpoint must contain a complete compatible RNG state and generator description.")
	invisible(TRUE)
}

.lame_tempered_run <- function(args, input, design, cx, call, requested_prior, resume = NULL) {
	kernel <- "tempered_binary_137_v1"
	temperatures <- .lame_tempered_temperatures(args$prior$uv_prior)
	K <- length(temperatures)
	base <- cx$n * cx$r + cx$r + cx$p + (cx$n - 1L) * (cx$tt + 1L)
	counts <- list(
		attempts = numeric(K),
		accepts = numeric(K),
		nonfinite = numeric(K),
		swap_attempts = numeric(K - 1L),
		swap_accepts = numeric(K - 1L),
		reflection_attempts = numeric(K),
		reflection_accepts = numeric(K),
		prior_refreshes = 0L
	)
	warmcounts <- counts
	if (is.null(resume)) {
		q <- matrix(0, base + 3L, K)
		beta <- numeric(cx$p)
		beta[1L] <- qnorm(min(.99, max(.01, mean(design$Y, na.rm = TRUE))))
		for (j in seq_len(K)) {
			q[, j] <- rnorm(base + 3L, sd = .1)
			q[cx$n * cx$r + seq_len(cx$r), j] <- rnorm(cx$r)
			q[cx$n * cx$r + cx$r + seq_len(cx$p), j] <- drop(chol(design$Q) %*% beta)
			q[base + 1L, j] <- log(args$prior$etaab * args$prior$Sab0[1L, 1L] / (args$prior$etaab + 2))
			q[base + 2L, j] <- .5 * log(args$prior$sigma_ab_scale / (args$prior$sigma_ab_shape + 1))
			q[base + 3L, j] <- atanh(args$prior$rho_ab_mean / .99)
		}
		if (!is.null(args$start_vals)) q <- args$start_vals$tempered
		state <- list(
			q = q,
			stepsizes = rep(.15, K),
			iteration = 0L,
			burn = args$burn,
			rng = get(".Random.seed", envir = .GlobalEnv),
			rng_kind = RNGkind()
		)
		if (isTRUE(cx$collapsed)) {
			auxiliary_seed <- sample.int(.Machine$integer.max, 1L)
			state$scale_rng <- .lame_tempered_scale_draws(
				q[, FALSE, drop = FALSE],
				cx,
				seed = auxiliary_seed
			)$rng
			state$rng <- get(".Random.seed", envir = .GlobalEnv)
		}
	} else {
		if (!identical(resume$kernel, kernel) || !identical(resume$context, cx) ||
				!identical(resume$Q, design$Q) || !identical(resume$temperatures, temperatures) ||
				!identical(resume$resolved_prior, args$prior))
			cli::cli_abort("The tempered checkpoint does not match this model, prior, or kernel.")
		state <- resume$state
		.lame_tempered_validate_rng(state$rng, state$rng_kind)
		do.call(RNGkind, as.list(state$rng_kind))
		assign(".Random.seed", state$rng, envir = .GlobalEnv)
	}
	if (!is.matrix(state$q) || !is.numeric(state$q) || is.complex(state$q) ||
			!identical(dim(state$q), c(base + 3L, K)) || any(!is.finite(state$q)) ||
			!is.numeric(state$stepsizes) || length(state$stepsizes) != K ||
			any(!is.finite(state$stepsizes) | state$stepsizes <= 0))
		cli::cli_abort("Invalid tempered replica state or step sizes.")
	if (isTRUE(cx$collapsed)) .lame_tempered_validate_rng(
		state$scale_rng,
		c("L'Ecuyer-CMRG", "Inversion", state$rng_kind[-c(1L, 2L)])
	)
	state$iteration <- .lame_tempered_integer(state$iteration, "checkpoint iteration", 0L)
	state$burn <- .lame_tempered_integer(state$burn, "checkpoint burn", 0L)
	start_iteration <- state$iteration
	end <- max(state$iteration, state$burn) + as.double(args$nscan)
	if (end > .Machine$integer.max) cli::cli_abort("The continued sampling schedule exceeds the integer limit.")
	end <- as.integer(end)
	checkpoint <- function() {
		parameters <- args
		parameters$start_vals <- NULL
		list(
			type = "lame_tempered",
			version = 1L,
			kernel = kernel,
			call_args = call,
			parameters = parameters,
			state = state,
			context = cx,
			Q = design$Q,
			temperatures = temperatures,
			resolved_prior = args$prior,
			requested_prior = requested_prior,
			counts = counts,
			warmcounts = warmcounts
		)
	}
	# state$rng belongs to the last complete block, so an interrupt inside C++
	# rolls back that block
	on.exit(.lame_tempered_checkpoint_write(checkpoint(), args$checkpoint_path), add = TRUE)
	draws <- modes <- scales <- list()
	index <- 0L
	started <- proc.time()[["elapsed"]]
	while (state$iteration < end) {
		if (proc.time()[["elapsed"]] - started >= args$max_seconds) break
		start <- state$iteration
		count <- min(100L, end - start, args$checkpoint_every - start %% args$checkpoint_every)
		if (start < state$burn) count <- min(count, state$burn - start)
		z <- tempered_binary_block(
			state$q,
			cx,
			temperatures,
			state$stepsizes,
			start,
			count,
			state$burn,
			args$odens
		)
		next_state <- state
		next_state$q <- z$q
		next_state$stepsizes <- z$stepsizes
		next_state$iteration <- start + count
		if (isTRUE(cx$collapsed) && ncol(z$draws)) {
			aux <- .lame_tempered_scale_draws(z$draws, cx, rng = state$scale_rng)
			next_state$scale_rng <- aux$rng
			scales[[index + 1L]] <- aux$values
		}
		next_state$rng <- get(".Random.seed", envir = .GlobalEnv)
		next_state$rng_kind <- RNGkind()
		state <- next_state
		if (start < state$burn) {
			for (nm in names(counts)) warmcounts[[nm]] <- warmcounts[[nm]] + z[[nm]]
		} else for (nm in names(counts)) counts[[nm]] <- counts[[nm]] + z[[nm]]
		if (ncol(z$draws)) {
			index <- index + 1L
			draws[[index]] <- z$draws
			modes[[index]] <- z$mode_trace
		}
		if (state$iteration %% args$checkpoint_every == 0L)
			.lame_tempered_checkpoint_write(checkpoint(), args$checkpoint_path)
		if (isTRUE(args$verbose) && (state$iteration %% 5000L == 0L || state$iteration == end))
			cli::cli_inform("Tempered sampling: {state$iteration}/{end} sweeps.")
	}
	list(
		q = if (length(draws)) do.call(cbind, draws) else matrix(numeric(), base + 3L, 0L),
		modes = if (length(modes)) do.call(cbind, modes) else matrix(numeric(), K, 0L),
		factor_variances = if (length(scales)) do.call(cbind, scales) else NULL,
		context = cx,
		Q = design$Q,
		prior = args$prior,
		temperatures = temperatures,
		stepsizes = state$stepsizes,
		counts = counts,
		warmcounts = warmcounts,
		kernel = kernel,
		start_iteration = start_iteration,
		end_iteration = state$iteration,
		terminated_early = state$iteration < end,
		terminal_state = state$q,
		checkpoint = checkpoint()
	)
}

.lame_tempered_scope <- function(e) {
	if (!identical(e$family, "binary") || !identical(e$mode, "unipartite") ||
			!isTRUE(e$symmetric) || !isTRUE(e$dynamic_ab) || !isFALSE(e$dynamic_uv) ||
			!isFALSE(e$dynamic_beta) || !isFALSE(e$dynamic_rho) || !isFALSE(e$dcor) ||
			!isTRUE(e$intercept) || !isTRUE(e$rvar) || !isTRUE(e$cvar) ||
			!identical(e$method, "mcmc"))
		cli::cli_abort(c(
			"The tempered sampler requires binary symmetric unipartite networks, static factors and coefficients, shared dynamic additive effects, and an intercept.",
			"i" = "Use {.code family = \"binary\", symmetric = TRUE, dynamic_ab = TRUE, dynamic_uv = FALSE, R > 0}."
		))
	for (nm in c(
		"Xrow", "Xcol", "R_row", "R_col", "odmax", "dynamic_beta_per_actor",
		"time_index", "period_exposure", "log_lik_path", "out_file"
	)) {
		if (!is.null(get(nm, envir = e, inherits = FALSE)))
			cli::cli_abort("{.arg {nm}} is not supported by the tempered sampler.")
	}
	if (isTRUE(e$periodic_save))
		cli::cli_abort("Use {.arg checkpoint_path} to save tempered sampler checkpoints.")
	if (!is.null(e$bias_correct) || !isTRUE(e$bootstrap == 0L))
		cli::cli_abort("Bias correction and bootstrap settings do not apply to the tempered MCMC sampler.")
	if (!identical(e$log_lik_method[1L], "observed_exact"))
		cli::cli_abort("The tempered sampler uses {.code log_lik_method = \"observed_exact\"}.")
	if (!identical(e$keep_snap_draws[1L], "none") || !identical(e$dynamic_beta_pool[1L], "none"))
		cli::cli_abort("Dynamic factor and coefficient storage/pooling settings do not apply to the tempered sampler.")
}

.lame_tempered_fit <- function(e) {
	.lame_tempered_scope(e)
	input <- .lame_tempered_input(e$Y, e$Xdyad, e$R)
	requested_prior <- e$.lame_input_prior
	prior <- .lame_tempered_prior(e$prior, input$Y)
	design <- .lame_tempered_design(input, e$g)
	cx <- .lame_tempered_context(design$Y, design$Xlist, input$R, design$Q, prior)
	options <- .lame_tempered_storage(e$posterior_opts)
	args <- mget(
		c(
			"seed", "nscan", "burn", "odens", "verbose", "plot", "gof",
			"save_log_lik", "max_seconds", "checkpoint_path", "checkpoint_every", "start_vals",
			"posterior_opts", "freeze_call", "model.name"
		),
		envir = e,
		inherits = FALSE
	)
	for (nm in c("nscan", "burn", "odens", "checkpoint_every"))
		args[[nm]] <- .lame_tempered_integer(args[[nm]], nm, if (nm == "burn") 0L else 1L)
	if (args$nscan %% args$odens != 0L)
		cli::cli_abort("For the tempered sampler, {.arg nscan} must be a multiple of {.arg odens}.")
	if (as.double(args$nscan) + args$burn > .Machine$integer.max)
		cli::cli_abort("{.arg burn} plus {.arg nscan} exceeds the integer sampling limit.")
	for (nm in c("verbose", "plot", "gof", "save_log_lik", "freeze_call")) {
		if (!is.logical(args[[nm]]) || length(args[[nm]]) != 1L || is.na(args[[nm]]))
			cli::cli_abort("For the tempered sampler, {.arg {nm}} must be TRUE or FALSE.")
	}
	if (!is.numeric(args$max_seconds) || is.complex(args$max_seconds) ||
			length(args$max_seconds) != 1L || is.na(args$max_seconds) || args$max_seconds <= 0)
		cli::cli_abort("{.arg max_seconds} must be a positive number or Inf.")
	if (!is.null(args$checkpoint_path)) {
		p <- args$checkpoint_path
		if (!is.character(p) || length(p) != 1L || is.na(p) || !nzchar(p) || dir.exists(p))
			cli::cli_abort("{.arg checkpoint_path} must name one file.")
		if (!dir.exists(dirname(p)) && !dir.create(dirname(p), recursive = TRUE))
			cli::cli_abort("Could not create checkpoint directory {.file {dirname(p)}}.")
	}
	resume <- NULL
	if (!is.null(args$start_vals)) {
		st <- args$start_vals
		if (!is.list(st) || length(st) != 1L || is.null(names(st)) ||
				!names(st) %in% c("tempered", "tempered_checkpoint"))
			cli::cli_abort("Tempered starting values must contain only {.code tempered}, a full replica-state matrix; use {.fn lame_resume} for exact continuation.")
		if (identical(names(st), "tempered_checkpoint")) {
			resume <- st$tempered_checkpoint
			if (!is.list(resume) || !identical(resume$type, "lame_tempered") ||
					!identical(resume$version, 1L))
				cli::cli_abort("Invalid tempered checkpoint.")
		}
	}
	args$prior <- prior
	args$Y <- input$Y; args$Xdyad <- input$Xdyad; args$R <- input$R; args$g <- design$g
	fit <- .lame_tempered_rng(args$seed, function() {
		a <- .lame_tempered_run(args, input, design, cx, e$mc, requested_prior, resume)
		f <- .lame_tempered_output(a, input, design, options, args$save_log_lik, args$gof)
		f$terminated_early <- a$terminated_early
		f$start_vals <- list(tempered = a$terminal_state)
		f$tempered_checkpoint <- a$checkpoint
		f
	})
	fit$call <- e$mc
	fit$model.name <- e$model.name
	fit$prior_supplied <- .lame_ab_prior_given(requested_prior)
	fit$events_per_actor_period <- mean(vapply(input$Y, function(y) sum(y > 0, na.rm = TRUE) / nrow(y), numeric(1)))
	fit$events_per_actor_total <- sum(
		vapply(input$Y, function(y) sum(y > 0, na.rm = TRUE), numeric(1))
	) / nrow(input$Y[[1L]])
	fit$refit_options <- e$.lame_refit_options
	fit$refit_options$sampler <- "tempered"
	fit$refit_options$prior <- requested_prior
	fit$refit_data <- e$.lame_refit_data
	fit$seed <- args$seed
	fit$checkpoint_path <- args$checkpoint_path
	fit$burn <- args$burn; fit$nscan <- args$nscan; fit$odens <- args$odens
	if (isTRUE(args$freeze_call)) {
		fit$freeze_call <- TRUE
		fit$data_snapshot <- e$.lame_freeze_snapshot
	}
	.lame_report_sampler_health(fit$sampler_health, fit$sampler_health_iterations)
	if (isTRUE(fit$no_draws_stored)) {
		cli::cli_warn("No posterior draw was retained before the time limit. Resume the checkpoint or increase {.arg max_seconds}.")
	} else {
		fit <- .lame_chain_policy_tempered(fit, e, args)
		.lame_separation_check(fit)
		if (isTRUE(args$plot)) plot(fit)
	}
	fit
}

.lame_tempered_resume <- function(checkpoint, nscan_more, overrides, block = NULL) {
	if (identical(overrides[["sampler"]], "tempered")) overrides$sampler <- NULL
	allowed <- c(
		"nscan", "checkpoint_path", "checkpoint_every", "max_seconds", "verbose",
		"plot", "gof", "posterior_opts", "save_log_lik", "freeze_call", "model.name"
	)
	if (length(overrides) && (is.null(names(overrides)) || anyNA(names(overrides)) || anyDuplicated(names(overrides)) ||
			any(!names(overrides) %in% allowed)))
		cli::cli_abort("Exact tempered continuation permits only sampling length, output, checkpoint, and time-limit overrides. Start a new fit to change data, priors, seed, or model.")
	args <- checkpoint$parameters
	if (!is.list(args) || is.null(checkpoint$state) || is.null(checkpoint$resolved_prior))
		cli::cli_abort("The tempered checkpoint is incomplete.")
	args$checkpoint_path <- NULL
	args$max_seconds <- Inf
	args$start_vals <- list(tempered_checkpoint = checkpoint)
	args$prior <- checkpoint$resolved_prior
	attr(args$prior, "lame_checkpoint_input_prior") <- checkpoint$requested_prior %||% list()
	if (!is.null(block)) attr(args$prior, "lame_chain_block") <- block
	args$family <- "binary"; args$symmetric <- TRUE; args$dynamic_ab <- TRUE
	args$sampler <- "tempered"
	args[names(overrides)] <- overrides
	if (!is.null(nscan_more)) args$nscan <- nscan_more
	# remaining warmup; the kernel keeps the absolute adaptation boundary and
	# thinning phase in checkpoint$state
	args$burn <- max(0L, checkpoint$state$burn - checkpoint$state$iteration)
	do.call(lame, args)
}

# separate saved stream for the hierarchical coordinate scales, so they and the
# primary chain resume exactly without altering the marginalised chain
.lame_tempered_scale_draws <- function(q, cx, rng = NULL, seed = NULL) {
	main_rng <- get(".Random.seed", envir = .GlobalEnv)
	main_kind <- RNGkind()
	on.exit({
		do.call(RNGkind, as.list(main_kind))
		assign(".Random.seed", main_rng, envir = .GlobalEnv)
	}, add = TRUE)
	RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion")
	if (is.null(rng)) set.seed(seed) else assign(".Random.seed", rng, envir = .GlobalEnv)
	values <- matrix(0, cx$r, ncol(q))
	for (s in seq_len(ncol(q))) values[, s] <- vapply(
		q[cx$n * cx$r + seq_len(cx$r), s],
		.lame_tempered_conditional_scale,
		numeric(1)
	)
	if (any(!is.finite(values) | values <= 0))
		cli::cli_abort("Nonfinite conditional factor scale in the tempered sampler.")
	list(values = values, rng = get(".Random.seed", envir = .GlobalEnv))
}
