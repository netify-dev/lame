# independently seeded chains per fit (seed + 1000 * (k - 1)) and the between-chain gate:
# every chain is extended from its final state while one drifts or a monitored split
# R-hat exceeds .lame_chain_rhat_max, up to prior$chain_extensions extensions

# threshold of the between-chain gate
.lame_chain_rhat_max <- 1.05
# stored draws per chain the gate needs (as for the drift check)
.lame_chain_gate_min_draws <- 60L

# number of chains when the user does not say: options(lame.chains), else 2
.lame_chains_default <- function() getOption("lame.chains", 2L)

#' @noRd
.lame_chains_check <- function(chains) {
	if (!is.numeric(chains) || length(chains) != 1L || !is.finite(chains) || chains < 1 || chains != floor(chains)) {
		cli::cli_abort("{.arg chains} must be a whole number of at least 1.", class = "lame_chains_arg")
	}
	as.integer(chains)
}

# the path of chain k when a single path was given: file_chain<k>.ext
.lame_chain_path <- function(path, k) {
	if (is.null(path)) return(NULL)
	if (grepl("\\.[^./\\\\]+$", path)) sub(
		"(\\.[^./\\\\]+)$",
		paste0("_chain", k, "\\1"),
		path
	) else paste0(path, "_chain", k)
}

# the operating system, as a function so tests can stand in for Windows
.lame_chains_os <- function() .Platform$OS.type

# the number of cores for the chains: the argument, else options(lame.cores), else all but one
.lame_chains_cores <- function(cores, chains) {
	if (is.null(cores)) cores <- getOption("lame.cores", NULL)
	if (is.null(cores)) cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
	if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores < 1) {
		cli::cli_abort("{.arg cores} must be a whole number of at least 1.", class = "lame_chains_arg")
	}
	if (identical(.lame_chains_os(), "windows")) return(1L)
	as.integer(min(chains, cores))
}

# run `fn(j)` for j in seq_len(n) on `cores` forked workers (sequentially for one)
.lame_chains_map <- function(n, fn, cores) {
	if (cores > 1L && n > 1L) {
		out <- parallel::mclapply(seq_len(n), fn, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
		bad <- vapply(out, inherits, logical(1), "try-error")
		if (any(bad)) cli::cli_abort(
			c("A chain worker stopped.", "i" = "{as.character(out[[which(bad)[1L]]])}"),
			class = "lame_chains_worker"
		)
		out
	} else lapply(seq_len(n), fn)
}

# one chain call; its warnings and messages are kept for the parent to re-signal,
# and a chain that is still moving comes back as its condition
.lame_chain_worker <- function(fun, args) {
	conds <- list()
	res <- withCallingHandlers(
		tryCatch(do.call(fun, args), lame_chain_unsettled = function(e) e, error = function(e) e),
		warning = function(w) { conds[[length(conds) + 1L]] <<- w; invokeRestart("muffleWarning") },
		message = function(m) { conds[[length(conds) + 1L]] <<- m; invokeRestart("muffleMessage") }
	)
	list(res = res, conds = conds)
}

# re-signal the warnings and messages of the chains, each distinct one once
.lame_chains_resignal <- function(rs, seen) {
	for (r in rs) for (cnd in r$conds) {
		key <- paste(class(cnd)[1L], conditionMessage(cnd))
		if (key %in% seen$keys) next
		seen$keys <- c(seen$keys, key)
		if (inherits(cnd, "warning")) warning(cnd) else message(cnd)
	}
	invisible(seen)
}

# the monitored series of one chain, a matrix of stored draws by series
.lame_gate_series <- function(fit) {
	reason <- .lame_ab_var_unident_reason(fit)
	hidden <- character(0)
	if (!is.na(reason)) hidden <- .lame_ab_var_hidden_rows(reason)
	coef_reason <- .lame_ab_coef_reason(fit)
	if (.lame_uv_hyper_unidentified(fit)) hidden <- c(hidden, "rho_uv", "sigma_uv")
	if (.lame_snap_rate_unidentified(fit)) hidden <- c(hidden, "pi_snap")
	cols <- list()
	B <- fit$BETA
	if (!is.null(B) && length(B) > 0L) {
		nm <- dimnames(B)[[2L]] %||% paste0("beta", seq_len(dim(B)[2L]))
		keep <- which(!nm %in% .lame_ab_coef_hidden_rows(coef_reason, nm))
		if (length(dim(B)) == 3L) {
			per <- dimnames(B)[[3L]] %||% paste0("p", seq_len(dim(B)[3L]))
			for (j in keep) for (t in seq_len(dim(B)[3L])) cols[[paste0(nm[j], "[", per[t], "]")]] <- B[, j, t]
		} else for (j in keep) cols[[nm[j]]] <- B[, j]
	}
	if (!is.null(fit$VC) && !is.null(colnames(fit$VC))) {
		for (v in setdiff(
			intersect(colnames(fit$VC), c("va", "vb", "cab", "ve", "rho")),
			hidden
		)) cols[[v]] <- fit$VC[, v]
	}
	for (v in setdiff(c("rho_ab", "sigma_ab", "rho_uv", "sigma_uv", "pi_snap"), hidden)) {
		if (is.numeric(fit[[v]]) && length(fit[[v]]) > 1L) cols[[v]] <- as.numeric(fit[[v]])
	}
	for (v in c("RHO_BETA", "SIGMA_BETA", "RHO")) {
		M <- fit[[v]]
		if (is.numeric(M) && length(M) > 1L) {
			M <- as.matrix(M)
			nmv <- colnames(M) %||% if (ncol(M) > 1L) paste0("g", seq_len(ncol(M))) else ""
			for (j in seq_len(ncol(M))) cols[[if (nzchar(nmv[j])) paste0(tolower(v), "[", nmv[j], "]") else tolower(v)]] <- M[, j]
		}
	}
	S <- fit$chain_summary
	if (is.matrix(S)) for (j in seq_len(ncol(S))) cols[[colnames(S)[j]]] <- S[, j]
	if (length(cols) == 0L) return(NULL)
	n <- min(lengths(cols))
	out <- do.call(cbind, lapply(cols, function(x) as.numeric(x)[seq_len(n)]))
	out
}

# split rank-normalised R-hat of every monitored series across chains (empty when the runs are too short)
.lame_gate_table <- function(fits) {
	empty <- data.frame(parameter = character(0), rhat = numeric(0), stringsAsFactors = FALSE)
	if (length(fits) < 2L) return(empty)
	ser <- lapply(fits, .lame_gate_series)
	if (any(vapply(ser, is.null, logical(1)))) return(empty)
	n <- min(vapply(ser, nrow, integer(1)))
	if (n < .lame_chain_gate_min_draws) return(empty)
	common <- Reduce(intersect, lapply(ser, colnames))
	rh <- vapply(
		common,
		function(v) .lame_rhat_chains(vapply(ser, function(m) m[seq_len(n), v], numeric(n))),
		numeric(1)
	)
	data.frame(parameter = common, rhat = unname(rh), stringsAsFactors = FALSE)
}

# the R-hat every reported surface uses, so summary(), compute_mcmc_diagnostics() and
# rhat_dynamic_beta() match fit$rhat; NA under 8 draws per chain or for non-finite or constant series
.lame_rhat_chains <- function(x) {
	x <- as.matrix(x)
	if (nrow(x) < 8L || ncol(x) < 2L) return(NA_real_)
	if (!all(is.finite(x)) || all(apply(x, 2L, stats::sd) < 1e-10 * max(1, max(abs(x))))) return(NA_real_)
	.lame_rhat(x)
}

# split rank-normalised R-hat (Vehtari et al. 2021) of a draws-by-chains matrix: the larger
# of the bulk and folded (tail) values, as posterior::rhat() reports it
.lame_rhat <- function(x) {
	x <- as.matrix(x)
	half <- nrow(x) %/% 2L
	if (half < 4L) return(NA_real_)
	split <- cbind(x[seq_len(half), , drop = FALSE], x[nrow(x) - half + seq_len(half), , drop = FALSE])
	basic <- function(m) {
		n <- nrow(m)
		w <- mean(apply(m, 2L, stats::var))
		b <- n * stats::var(colMeans(m))
		if (!is.finite(w) || w <= 0) return(NA_real_)
		sqrt(((n - 1) / n * w + b / n) / w)
	}
	z <- function(m) {
		r <- rank(m, ties.method = "average")
		matrix(stats::qnorm((r - 3 / 8) / (length(m) + 1 / 4)), nrow(m), ncol(m))
	}
	bulk <- basic(z(split))
	tail <- basic(z(abs(split - stats::median(split))))
	max(bulk, tail, na.rm = TRUE)
}

# the series above the threshold
.lame_gate_flagged <- function(tab, rhat_max = .lame_chain_rhat_max) {
	tab[is.finite(tab$rhat) & tab$rhat > rhat_max, , drop = FALSE]
}

# run `chains` chains of `fun` ("ame" or "lame") on the evaluated `args`; returns the pooled
# fit, or the chain list when `as_list`; `start_list` holds per-chain start values
.lame_chains_run <- function(
	fun,
	args,
	chains,
	seeds,
	call,
	cores = NULL,
	caller = parent.frame(),
	as_list = FALSE,
	start_list = NULL,
	verbose = FALSE,
	plot = FALSE,
	resolved = NULL,
	input_prior = NULL
) {
	fn <- get(fun, envir = asNamespace("lame"), mode = "function")
	prior <- args$prior %||% list()
	K <- as.integer(.lame_chain_limit(prior))
	nscan <- args$nscan %||% eval(formals(fn)$nscan)
	burn <- args$burn %||% eval(formals(fn)$burn)
	cores <- .lame_chains_cores(cores, chains)
	seen <- new.env(); seen$keys <- character(0)
	if (is.null(start_list)) start_list <- vector("list", chains)
	if (!is.null(start_list[[1L]])) burn <- 0L
	if (is.null(input_prior)) input_prior <- prior
	# the arguments of chain j in round r (0: the first run)
	chain_args <- function(j, r, prev) {
		a <- args
		a$chains <- 1L
		a$cores <- NULL
		a$seed <- seeds[j] + (if (r > 0L) burn + r * nscan else 0)
		a$verbose <- FALSE
		a$plot <- FALSE
		for (nm in c("out_file", "checkpoint_path")) if (!is.null(a[[nm]])) a[[nm]] <- .lame_chain_path(a[[nm]], j)
		cont <- r > 0L || !is.null(start_list[[j]])
		use_resolved <- cont && identical(fun, "lame") && !is.null(resolved)
		pr <- if (use_resolved) resolved else prior
		attr(pr, "lame_chain_block") <- list(mc = call, k = r, record = TRUE, nopolicy = TRUE)
		if (use_resolved) attr(pr, "lame_checkpoint_input_prior") <- input_prior
		a$prior <- pr
		if (r > 0L) { a$start_vals <- prev[[j]]; a$burn <- 0L }
		else if (!is.null(start_list[[j]])) { a$start_vals <- start_list[[j]]; a$burn <- 0L }
		a
	}
	run_round <- function(r, prev) {
		rs <- .lame_chains_map(chains, function(j) .lame_chain_worker(fn, chain_args(j, r, prev)), cores)
		.lame_chains_resignal(rs, seen)
		res <- lapply(rs, `[[`, "res")
		errs <- Filter(function(x) inherits(x, "error") && !inherits(x, "lame_chain_unsettled"), res)
		if (length(errs) > 0L) stop(errs[[1L]])
		res
	}
	say <- function(r, res) {
		if (!isTRUE(verbose)) return(invisible())
		for (j in seq_len(chains)) cli::cli_inform(
			c(
				"i" = "Chain {j} of {chains} (seed {seeds[j] + (if (r > 0L) burn + r * nscan else 0)}): {if (r == 0L) 'run' else paste0('extension ', r)} complete."
			)
		)
	}
	triggers <- list()
	r <- 0L
	prev <- NULL
	res <- run_round(0L, NULL)
	say(0L, res)
	repeat {
		own <- Filter(function(x) inherits(x, "lame_chain_unsettled"), res)
		fits <- lapply(res, function(x) if (inherits(x, "lame_chain_unsettled")) x$fit else x)
		if (identical(fun, "lame") && is.null(resolved)) {
			resolved <- fits[[1L]]$prior
			if (!is.null(fits[[1L]]$refit_options$prior)) input_prior <- fits[[1L]]$refit_options$prior
		}
		if (identical(fun, "lame") && is.null(resolved)) stop("internal: no resolved prior")
		early <- any(vapply(fits, function(f) isTRUE(f$terminated_early), logical(1)))
		moving <- lapply(fits, function(f) if (early) NULL else .lame_chain_moving(f))
		gate <- if (early) .lame_gate_table(list()) else .lame_gate_table(fits)
		bad <- .lame_gate_flagged(gate)
		drift_txt <- unlist(
			lapply(
				seq_along(moving),
				function(j) if (!is.null(moving[[j]]) && nrow(moving[[j]]) > 0L) paste0(
					"chain ",
					j,
					": ",
					moving[[j]]$text
				)
			)
		)
		rhat_txt <- if (nrow(bad) > 0L) sprintf("%s (R-hat %.3f)", bad$parameter, bad$rhat) else character(0)
		if (length(drift_txt) == 0L && length(rhat_txt) == 0L) break
		if (isTRUE(verbose)) cli::cli_inform(
			c("i" = "The chains are not settled ({paste(c(drift_txt, rhat_txt), collapse = ', ')}); round {r} of {K}.")
		)
		if (r >= K) {
			pooled <- .lame_chains_pool(
				fits,
				seeds,
				burn,
				nscan,
				r,
				triggers,
				call,
				args,
				gate,
				input_prior,
				as_list = FALSE,
				refused = TRUE
			)
			cli::cli_abort(
				c(
					if (r == 0L && K == 0L) "The chains were not settled when they stopped and {.code prior$chain_extensions} is 0, so they were not extended: {c(drift_txt, rhat_txt)}."
					else if (r == 0L) "The chains were not settled when they stopped: {c(drift_txt, rhat_txt)}."
					else "The chains were still not settled after {r} extension{?s} ({format(r * nscan, big.mark = ',')} sweeps discarded): {c(drift_txt, rhat_txt)}.",
					"i" = "The pooled fit of the last block is in the {.field fit} field of the condition and the per-chain fits in {.field fits}; the offending series are in {.field rhat}. Continue every chain with {.code lame_resume(e$fit, nscan = {nscan})}, or one chain with {.code lame_resume(e$fits[[k]], nscan = {nscan})}, where {.code e} is the caught condition."
				),
				class = "lame_chain_unsettled",
				fit = pooled,
				fits = fits,
				rhat = gate,
				moving = c(drift_txt, rhat_txt),
				extensions = as.integer(r),
				call = caller
			)
		}
		triggers[[length(triggers) + 1L]] <- c(drift_txt, rhat_txt)
		prev <- lapply(fits, `[[`, "start_vals")
		rm(fits, res, own)
		r <- r + 1L
		res <- run_round(r, prev)
		say(r, res)
	}
	pooled <- .lame_chains_pool(
		fits,
		seeds,
		burn,
		nscan,
		r,
		triggers,
		call,
		args,
		gate,
		input_prior,
		as_list = as_list,
		refused = FALSE
	)
	pooled
}

# the pooled fit (or the list of chain fits) with the chain records
.lame_chains_pool <- function(
	fits,
	seeds,
	burn,
	nscan,
	extensions,
	triggers,
	call,
	args,
	gate,
	input_prior,
	as_list = FALSE,
	refused = FALSE
) {
	chains <- length(fits)
	K <- as.integer(.lame_chain_limit(args$prior %||% list()))
	next_seeds <- seeds + burn + (extensions + 1L) * nscan
	rec <- .lame_chain_record(K, extensions, nscan, triggers, NA_real_)
	rec$next_seed <- NULL
	rec$rhat_max <- .lame_chain_rhat_max
	rec$rhat <- gate
	rec$chains <- chains
	rec$next_seeds <- next_seeds
	if (refused) { rec$refused <- TRUE; rec$call_args <- args }
	for (j in seq_len(chains)) fits[[j]]$chain_id <- j
	if (as_list) {
		for (j in seq_len(chains)) fits[[j]]$auto_extend <- rec
		return(fits)
	}
	start_list <- lapply(fits, `[[`, "start_vals")
	diag <- tryCatch(suppressMessages(suppressWarnings(compute_mcmc_diagnostics(fits))), error = function(e) NULL)
	fits_slim <- fits
	for (j in seq_len(chains)) fits_slim[[j]]$chain_summary <- NULL
	empty <- all(vapply(fits, .lame_has_no_draws, logical(1)))
	pooled <- if (empty) {
		# no chain stored a draw: the fit of the first chain, with the chain records
		f <- fits_slim[[1L]]
		f$n_chains <- chains
		f$is_combined <- TRUE
		f$chain_indicator <- integer(0)
		f
	} else suppressWarnings(suppressMessages(combine_ame_chains(fits_slim, diagnostics = FALSE)))
	pooled$chain_summary <- NULL
	pooled$diagnostics <- diag
	pooled$call <- call
	pooled$chains <- chains
	pooled$chain_seeds <- seeds
	pooled$rhat <- gate
	pooled$chain_start_vals <- start_list
	pooled$auto_extend <- rec
	pooled$terminated_early <- any(vapply(fits, function(f) isTRUE(f$terminated_early), logical(1)))
	if (empty) pooled <- .lame_empty_draws(pooled)
	pooled
}

# the chains of a user's ame() call; `env` is the frame of ame() after its early
# argument handling
.lame_chains_from_ame <- function(env, mc, chains) {
	fm <- setdiff(names(formals(ame)), c("...", "print", "chains", "n_chains", "cores"))
	args <- mget(intersect(fm, ls(env, all.names = TRUE)), envir = env)
	args$n_chains <- 1
	cores <- get("cores", envir = env, inherits = FALSE)
	.lame_chains_entry("ame", args, chains, mc, cores, parent.frame(2L))
}

# the chains of a user's lame() call: `ctx` holds the evaluated arguments the user
# supplied (after any conversion of the input)
.lame_chains_from_lame <- function(env, ctx, chains) {
	args <- ctx$args
	for (nm in c("seed", "burn", "nscan", "odens", "prior", "verbose", "plot")) {
		args[nm] <- list(get(nm, envir = env, inherits = FALSE))
	}
	args$print <- NULL
	args$chains <- NULL
	cores <- get("cores", envir = env, inherits = FALSE)
	args$cores <- NULL
	.lame_chains_entry("lame", args, chains, ctx$mc, cores, ctx$caller)
}

.lame_chains_entry <- function(fun, args, chains, call, cores, caller) {
	prior <- args$prior %||% list()
	as_list <- isTRUE(attr(prior, "lame_chain_list", exact = TRUE))
	if (as_list) { attr(prior, "lame_chain_list") <- NULL; args$prior <- prior }
	if (!is.numeric(args$seed) || length(args$seed) != 1L || !is.finite(args$seed)) args$seed <- 6886
	seeds <- as.numeric(args$seed) + 1000 * (seq_len(chains) - 1L)
	verbose <- isTRUE(args$verbose)
	plot <- isTRUE(args$plot)
	out <- .lame_chains_run(
		fun,
		args,
		chains,
		seeds,
		call,
		cores = cores,
		caller = caller,
		as_list = as_list,
		verbose = verbose,
		plot = plot
	)
	if (plot && !as_list) try(plot(out), silent = TRUE)
	out
}
