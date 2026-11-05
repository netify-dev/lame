# chain policy of ame() and lame(): continue a chain that is still moving at the end of the
# post-burn run, and refuse with a classed error if it still moves after the extensions

# thresholds of the check (see .lame_chain_moving)
.lame_chain_z_max <- 3
.lame_chain_sd_min <- 0.5
.lame_chain_n_batch <- 10L
.lame_chain_default_extensions <- 3L

# scalar series the check monitors: intercept (or first coefficient), va, vb, ve, rho and the dynamic
# hyperparameters, leaving out quantities the fit reports as not estimated
.lame_chain_series <- function(fit) {
	reason <- .lame_ab_var_unident_reason(fit)
	hidden <- character(0)
	if (!is.na(reason)) hidden <- .lame_ab_var_hidden_rows(reason)
	coef_reason <- .lame_ab_coef_reason(fit)
	if (.lame_uv_hyper_unidentified(fit)) hidden <- c(hidden, "rho_uv", "sigma_uv")
	if (.lame_snap_rate_unidentified(fit)) hidden <- c(hidden, "pi_snap")
	series <- list()
	B <- fit$BETA
	if (!is.null(B) && length(B) > 0L) {
		if (length(dim(B)) == 3L) B <- apply(B, c(1L, 2L), mean)
		if (NCOL(B) > 0L) {
			nm <- colnames(B) %||% paste0("beta", seq_len(NCOL(B)))
			keep <- setdiff(nm, .lame_ab_coef_hidden_rows(coef_reason, nm))
			if (length(keep) > 0L) {
				k <- which(nm == (if ("intercept" %in% keep) "intercept" else keep[1L]))[1L]
				series[[nm[k]]] <- B[, k]
			}
		}
	}
	if (!is.null(fit$VC) && !is.null(colnames(fit$VC))) {
		for (v in setdiff(intersect(colnames(fit$VC), c("va", "vb", "ve", "rho")), hidden)) series[[v]] <- fit$VC[, v]
	}
	for (v in setdiff(c("rho_ab", "sigma_ab", "rho_uv", "sigma_uv", "pi_snap"), hidden)) {
		if (is.numeric(fit[[v]]) && length(fit[[v]]) > 1L) series[[v]] <- as.numeric(fit[[v]])
	}
	series
}

# first third against last third of each series with batch-means standard errors (n_batch batches
# per third): both means, the standardised difference z and the shift in posterior sds
.lame_chain_stats <- function(series, n_batch = .lame_chain_n_batch) {
	empty <- data.frame(
		parameter = character(0),
		mean_first = numeric(0),
		mean_last = numeric(0),
		z = numeric(0),
		shift = numeric(0),
		stringsAsFactors = FALSE
	)
	if (length(series) == 0L) return(empty)
	n <- min(lengths(series))
	if (n < 6L * n_batch) return(empty)
	batch_se <- function(s) {
		id <- rep(seq_len(n_batch), each = ceiling(length(s) / n_batch))[seq_along(s)]
		bm <- tapply(s, id, mean)
		stats::sd(bm) / sqrt(length(bm))
	}
	third <- n %/% 3L
	rows <- lapply(names(series), function(v) {
		s <- series[[v]][seq_len(n)]
		# skip a constant series; the floor is relative so a normal fit's check does not depend on the
		# response unit
		if (!all(is.finite(s)) || stats::sd(s) <= 1e-10 * max(abs(s))) return(NULL)
		a <- s[seq_len(third)]
		b <- s[(n - third + 1L):n]
		z <- (mean(b) - mean(a)) / sqrt(batch_se(a)^2 + batch_se(b)^2)
		data.frame(
			parameter = v,
			mean_first = mean(a),
			mean_last = mean(b),
			z = z,
			shift = (mean(b) - mean(a)) / stats::sd(s),
			stringsAsFactors = FALSE
		)
	})
	rows <- rows[!vapply(rows, is.null, logical(1))]
	if (length(rows) == 0L) return(empty)
	do.call(rbind, rows)
}

# series whose first and last thirds differ by more than z_max standard errors and sd_min posterior
# sds; zero rows when the chain has settled, with a `text` column "name (a to b)"
.lame_chain_moving <- function(
	fit,
	z_max = .lame_chain_z_max,
	sd_min = .lame_chain_sd_min,
	n_batch = .lame_chain_n_batch
) {
	st <- .lame_chain_stats(.lame_chain_series(fit), n_batch)
	flag <- is.finite(st$z) & abs(st$z) > z_max & abs(st$shift) > sd_min
	st <- st[flag, , drop = FALSE]
	st$text <- sprintf("%s (%.3g to %.3g)", st$parameter, st$mean_first, st$mean_last)
	st
}

# prior$chain_extensions: continuations of a flagged chain before refusal (0 refuses at the first flag)
.lame_chain_limit <- function(prior) {
	k <- if (is.list(prior)) prior[["chain_extensions"]] else NULL
	if (is.null(k)) return(.lame_chain_default_extensions)
	k
}

.lame_chain_check_prior <- function(prior) {
	k <- if (is.list(prior)) prior[["chain_extensions"]] else NULL
	if (!is.null(k) && (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 0 || k != floor(k))) {
		cli::cli_abort("prior$chain_extensions must be a non-negative whole number.")
	}
	invisible(TRUE)
}

# record carried as fit$auto_extend by every checked fit
.lame_chain_record <- function(
	K,
	extensions = 0L,
	nscan = 0L,
	triggers = list(),
	next_seed = NA_real_,
	skipped = NULL
) {
	rec <- list(
		max_extensions = as.integer(K),
		extensions = as.integer(extensions),
		discarded_sweeps = as.integer(extensions) * as.integer(nscan),
		triggers = triggers,
		z_max = .lame_chain_z_max,
		sd_min = .lame_chain_sd_min,
		next_seed = next_seed
	)
	if (!is.null(skipped)) rec$skipped <- skipped
	rec
}

# chain policy for lame()

# captured at the top of lame(): either a continuation block of an extension, or the user's call
# with the evaluated arguments a continuation reruns
.lame_chain_enter <- function(env, mc, caller) {
	prior <- get("prior", envir = env, inherits = FALSE)
	blk <- attr(prior, "lame_chain_block", exact = TRUE)
	if (!is.null(blk)) return(
		list(block = TRUE, mc = blk$mc, record = isTRUE(blk$record), nopolicy = isTRUE(blk$nopolicy))
	)
	supplied <- setdiff(names(mc)[-1L], c("print", "resume_from"))
	nin <- get0("netify_input", envir = env, inherits = FALSE)
	if (is.list(nin) && isTRUE(nin$converted)) {
		supplied <- union(supplied, c("Y", "Xdyad", "Xrow", "Xcol", "mode", "family", "symmetric"))
	}
	supplied <- union(supplied, "verbose")
	ok <- TRUE
	args <- list()
	for (nm in supplied) {
		v <- tryCatch(list(get(nm, envir = env, inherits = FALSE)), error = function(e) NULL)
		if (is.null(v)) { ok <- FALSE; next }
		args[nm] <- v
	}
	K <- tryCatch(.lame_chain_limit(prior), error = function(e) .lame_chain_default_extensions)
	list(block = FALSE, mc = mc, args = args, ok = ok, K = K, caller = caller, record = FALSE)
}

# policy at the end of lame(): done = TRUE means `fit` is the value of the extended call; otherwise
# the caller attaches `record` to its own fit
.lame_chain_policy <- function(fit, ctx, env) {
	if (isTRUE(ctx$block)) return(list(done = FALSE, record = NULL))
	K <- as.integer(ctx$K)
	nscan <- get("nscan", envir = env, inherits = FALSE)
	burn <- get("burn", envir = env, inherits = FALSE)
	seed <- get("seed", envir = env, inherits = FALSE)
	next_seed <- as.numeric(seed) + burn + nscan
	if (isTRUE(fit$terminated_early)) {
		return(list(done = FALSE, record = .lame_chain_record(K, next_seed = next_seed, skipped = "max_seconds")))
	}
	moving <- .lame_chain_moving(fit)
	if (nrow(moving) == 0L) return(list(done = FALSE, record = .lame_chain_record(K, next_seed = next_seed)))
	input_prior <- get(".lame_input_prior", envir = env, inherits = FALSE)
	verbose <- isTRUE(get("verbose", envir = env, inherits = FALSE))
	max_seconds <- get("max_seconds", envir = env, inherits = FALSE)
	elapsed <- as.numeric(difftime(Sys.time(), get(".lame_start_time", envir = env, inherits = FALSE), units = "secs"))
	call_args <- ctx$args
	call_args$verbose <- verbose
	meta <- fit[intersect(c("call", "refit_options", "refit_data", "als_start"), names(fit))]
	resolved_prior <- fit$prior
	patch <- function(f) {
		for (nm in names(meta)) f[[nm]] <- meta[[nm]]
		f
	}
	spec <- list(
		K = K,
		nscan = nscan,
		seed = as.numeric(seed),
		burn = burn,
		verbose = verbose,
		call_env = env,
		resume_args = call_args,
		sampler = "lame",
		max_seconds = max_seconds,
		elapsed = elapsed,
		t0 = Sys.time()
	)
	if (K == 0L || !isTRUE(ctx$ok)) {
		.lame_chain_refuse(patch(fit), spec, list(moving$text), moving$text, 0L)
	}
	prev <- fit["start_vals"]
	first_moving <- moving$text
	# free the sampler's frame; the flagged block's draws are discarded and only its final state continues
	rm(list = setdiff(ls(env, all.names = TRUE), c(".had_seed", ".old_seed")), envir = env)
	rm(fit, moving)
	run_block <- function(prev, k, seed_k, max_seconds_k) {
		a <- call_args
		pr <- resolved_prior
		attr(pr, "lame_checkpoint_input_prior") <- input_prior
		attr(pr, "lame_chain_block") <- list(mc = ctx$mc, k = k)
		a$prior <- pr
		a$start_vals <- prev$start_vals
		a$burn <- 0L
		a$nscan <- nscan
		a$seed <- seed_k
		if (is.finite(max_seconds_k)) a$max_seconds <- max_seconds_k
		do.call(lame, a, envir = ctx$caller)
	}
	list(done = TRUE, fit = .lame_chain_extend(prev, list(first_moving), spec, run_block, patch))
}

# continue the chain from `prev` until a block settles or K extensions are made; `run_block` runs
# nscan sweeps from the previous state and `patch` restores the user's call metadata on its fit
.lame_chain_extend <- function(prev, triggers, spec, run_block, patch) {
	K <- spec$K
	for (k in seq_len(K)) {
		if (isTRUE(spec$verbose)) {
			now <- paste(utils::tail(triggers, 1L)[[1L]], collapse = ", ")
			cli::cli_inform(c("i" = "The chain is still moving ({now}): extending ({k} of {K})."))
		}
		remaining <- Inf
		if (is.finite(spec$max_seconds)) {
			remaining <- max(
				1e-9,
				spec$max_seconds - (spec$elapsed + as.numeric(difftime(Sys.time(), spec$t0, units = "secs")))
			)
		}
		seed_k <- spec$seed + spec$burn + k * spec$nscan
		fit_k <- patch(run_block(prev, k, seed_k, remaining))
		if (isTRUE(fit_k$terminated_early)) {
			return(.lame_chain_refuse(fit_k, spec, triggers, character(0), k, reason = "max_seconds"))
		}
		mv <- .lame_chain_moving(fit_k)
		if (nrow(mv) == 0L) {
			fit_k$auto_extend <- .lame_chain_record(
				spec$K,
				k,
				spec$nscan,
				triggers,
				spec$seed + spec$burn + (k + 1L) * spec$nscan
			)
			return(fit_k)
		}
		if (k == K) return(.lame_chain_refuse(fit_k, spec, triggers, mv$text, k))
		triggers[[length(triggers) + 1L]] <- mv$text
		prev <- fit_k[intersect(c("start_vals", "tempered_checkpoint"), names(fit_k))]
		rm(fit_k, mv)
	}
}

# classed error carrying the fit of the last block
.lame_chain_refuse <- function(fit, spec, triggers, moving, extensions, reason = "moving") {
	nscan <- spec$nscan
	rec <- .lame_chain_record(spec$K, extensions, nscan, triggers, spec$seed + spec$burn + (extensions + 1L) * nscan)
	rec$refused <- TRUE
	rec$moving <- moving
	rec$call_args <- if (is.function(spec$resume_args)) spec$resume_args() else spec$resume_args
	rec$sampler <- spec$sampler
	fit$auto_extend <- rec
	sweeps <- format(extensions * nscan, big.mark = ",")
	nsc <- format(nscan, big.mark = ",")
	made <- if (identical(reason, "max_seconds")) {
		"The chain was still moving and the time under {.arg max_seconds} ran out during extension {extensions}"
	} else if (extensions == 0L && spec$K == 0L) {
		"The chain was still moving when it stopped and {.code prior$chain_extensions} is 0, so it was not extended"
	} else if (extensions == 0L) {
		"The chain was still moving when it stopped"
	} else {
		"The chain was still moving after {extensions} extension{?s} ({sweeps} sweeps discarded)"
	}
	moving_txt <- if (length(moving)) ": {moving}." else "."
	fix <- paste0(
		"Continue the same chain with {.code lame_resume(e$fit, nscan = ",
		nscan,
		")}, where {.code e} is the caught condition, or refit with a larger {.arg nscan} and {.arg burn}."
	)
	cli::cli_abort(
		c(
			paste0(made, moving_txt),
			"i" = paste0(
				"The fit with the draws of its last block of {nsc} sweeps is in the {.field fit} field of the condition. ",
				fix
			)
		),
		class = "lame_chain_unsettled",
		fit = fit,
		moving = moving,
		extensions = as.integer(extensions),
		call = spec$call_env
	)
}

# chain policy for ame()

# policy at the end of ame(): `fun` is the sampler that produced `fit`, `env` the frame of ame()
.ame_chain_policy <- function(fit, env, fun) {
	prior <- get("prior", envir = env, inherits = FALSE)
	K <- as.integer(.lame_chain_limit(prior))
	nscan <- get("nscan", envir = env, inherits = FALSE)
	burn <- get("burn", envir = env, inherits = FALSE)
	seed <- get("seed", envir = env, inherits = FALSE)
	next_seed <- as.numeric(seed) + burn + nscan
	moving <- .lame_chain_moving(fit)
	if (nrow(moving) == 0L) {
		fit$auto_extend <- .lame_chain_record(K, next_seed = next_seed)
		return(fit)
	}
	nm <- intersect(names(formals(fun)), ls(env, all.names = TRUE))
	call_args <- mget(nm, envir = env)
	mc <- get("mc", envir = env, inherits = FALSE)
	resume_args <- function() {
		nm <- setdiff(intersect(names(formals(ame)), ls(env, all.names = TRUE)), c("print", "..."))
		mget(nm, envir = env)
	}
	spec <- list(
		K = K,
		nscan = nscan,
		seed = as.numeric(seed),
		burn = burn,
		verbose = isTRUE(call_args$verbose),
		call_env = env,
		max_seconds = Inf,
		resume_args = resume_args,
		sampler = "ame"
	)
	patch <- function(f) { f$call <- mc; f }
	if (K == 0L) .lame_chain_refuse(patch(fit), spec, list(moving$text), moving$text, 0L)
	prev <- fit["start_vals"]
	first_moving <- moving$text
	rm(list = "fit", envir = env)
	rm(fit, moving)
	run_block <- function(prev, k, seed_k, max_seconds_k) {
		a <- call_args
		a$start_vals <- prev$start_vals
		a$burn <- 0L
		a$seed <- seed_k
		do.call(fun, a)
	}
	.lame_chain_extend(prev, list(first_moving), spec, run_block, patch)
}

# print() line for a fit whose chain was extended or pooled
.lame_chain_print_line <- function(rec, fit = NULL) {
	out <- character(0)
	if (!is.null(fit) && isTRUE(fit$is_combined) && !is.null(fit$chains)) {
		rh <- fit$rhat
		top <- if (is.data.frame(rh) && any(is.finite(rh$rhat))) max(rh$rhat, na.rm = TRUE) else NA_real_
		out <- paste0(
			fit$chains,
			" chains (seeds ",
			paste(fit$chain_seeds, collapse = ", "),
			")",
			if (is.finite(top)) paste0(
				"; largest R-hat of the monitored set ",
				formatC(top, digits = 3, format = "f")
			) else "",
			"."
		)
	}
	if (is.list(rec) && isTRUE(rec$extensions > 0L)) {
		out <- c(
			out,
			paste0(
				"The chains were extended ",
				rec$extensions,
				if (rec$extensions == 1L) " time" else " times",
				" (",
				format(rec$discarded_sweeps, big.mark = ","),
				" sweeps of draws discarded); not settled: ",
				paste(unique(sub(" .*", "", sub("^chain [0-9]+: ", "", unlist(rec$triggers)))), collapse = ", "),
				"."
			)
		)
	}
	if (length(out) == 0L) NULL else out
}

# chain policy for sampler = "tempered"

# the tempered sampler continues exactly from its checkpoint, so a flagged chain resumes it for
# another nscan sweeps
.lame_chain_policy_tempered <- function(fit, e, args) {
	ctx <- e$.chain_ctx
	if (isTRUE(ctx$block)) return(fit)
	K <- as.integer(ctx$K)
	nscan <- args$nscan
	moving <- .lame_chain_moving(fit)
	if (nrow(moving) == 0L) {
		fit$auto_extend <- .lame_chain_record(K)
		return(fit)
	}
	meta <- fit[intersect(
		c("call", "refit_options", "refit_data", "seed", "burn", "nscan", "odens", "checkpoint_path"),
		names(fit)
	)]
	patch <- function(f) {
		for (nm in names(meta)) f[[nm]] <- meta[[nm]]
		f
	}
	spec <- list(
		K = K,
		nscan = nscan,
		seed = 0,
		burn = 0L,
		verbose = isTRUE(args$verbose),
		call_env = e,
		max_seconds = Inf,
		resume_args = ctx$args,
		sampler = "tempered"
	)
	if (K == 0L || !isTRUE(ctx$ok) || is.null(fit$tempered_checkpoint)) {
		.lame_chain_refuse(patch(fit), spec, list(moving$text), moving$text, 0L)
	}
	prev <- fit["tempered_checkpoint"]
	first_moving <- moving$text
	rm(fit, moving)
	run_block <- function(prev, k, seed_k, max_seconds_k) {
		.lame_tempered_resume(prev$tempered_checkpoint, nscan, list(), block = list(mc = ctx$mc, k = k))
	}
	.lame_chain_extend(prev, list(first_moving), spec, run_block, patch)
}

# global summaries recorded at every stored draw for the between-chain gate: rms of a and b (per
# period when dynamic), rms of U G V', norm of G (bipartite), log U/V scale split (directed) and the
# largest absolute actor effect; fit$chain_summary holds the draws-by-summaries matrix

.lame_chain_draw <- function(a, b, U = NULL, V = NULL, G = NULL, symmetric = FALSE) {
	rms <- function(x) sqrt(mean(x^2))
	per <- function(x) if (is.matrix(x)) x else matrix(x, ncol = 1L)
	A <- if (is.null(a)) NULL else per(a)
	B <- if (is.null(b)) NULL else per(b)
	nm <- function(base, TT) if (TT > 1L) paste0(base, "[p", seq_len(TT), "]") else base
	out <- numeric(0)
	if (!is.null(A)) {
		out <- c(
			out,
			stats::setNames(apply(A, 2L, rms), nm("a_rms", ncol(A))),
			stats::setNames(apply(abs(A), 2L, max), nm("a_max", ncol(A)))
		)
	}
	if (!is.null(B)) {
		out <- c(
			out,
			stats::setNames(apply(B, 2L, rms), nm("b_rms", ncol(B))),
			stats::setNames(apply(abs(B), 2L, max), nm("b_max", ncol(B)))
		)
	}
	if (!is.null(U) && !is.null(V) && length(dim(U)) >= 2L && ncol(U) > 0L) {
		dyn <- length(dim(U)) == 3L
		Tu <- if (dyn) dim(U)[3L] else 1L
		uv <- split <- numeric(Tu)
		for (t in seq_len(Tu)) {
			Ut <- if (dyn) matrix(U[, , t], nrow = dim(U)[1L]) else as.matrix(U)
			Vt <- if (dyn) matrix(V[, , t], nrow = dim(V)[1L]) else as.matrix(V)
			M <- if (is.null(G)) Ut else Ut %*% G
			uv[t] <- sqrt(sum(crossprod(M) * crossprod(Vt)) / (nrow(Ut) * nrow(Vt)))
			split[t] <- log(sqrt(sum(Ut^2)) / sqrt(sum(Vt^2)))
		}
		out <- c(out, stats::setNames(uv, nm("uv_rms", Tu)))
		if (!isTRUE(symmetric)) out <- c(out, stats::setNames(split, nm("uv_split", Tu)))
	}
	if (!is.null(G)) out <- c(out, g_norm = sqrt(sum(G^2)))
	out
}

# sampler hook at every stored draw; creates the record on the first call (`acc` NULL)
.lame_chain_acc_update <- function(
	acc,
	iter,
	n_draws,
	a,
	b,
	U = NULL,
	V = NULL,
	G = NULL,
	symmetric = FALSE,
	actors = NULL
) {
	if (n_draws < 1L) return(acc)
	v <- tryCatch(.lame_chain_draw(a, b, U, V, G, symmetric), error = function(e) NULL)
	if (is.null(v)) return(acc)
	if (is.null(acc)) {
		acc <- new.env(parent = emptyenv())
		acc$draws <- matrix(NA_real_, n_draws, length(v), dimnames = list(NULL, names(v)))
		acc$n <- 0L
	}
	if (iter <= nrow(acc$draws) && length(v) == ncol(acc$draws)) {
		acc$draws[iter, ] <- v
		acc$n <- max(acc$n, iter)
	}
	acc
}

# matrix of draws by summaries, kept only for complete runs
.lame_chain_acc_finish <- function(acc) {
	if (is.null(acc) || acc$n < nrow(acc$draws)) return(NULL)
	acc$draws
}
