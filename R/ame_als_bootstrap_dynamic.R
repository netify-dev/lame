# replicate loop of the bootstrap of a penalised dynamic ALS fit, under the rules of
# ame_als_bootstrap(): a replicate that does not converge is retried at four times the
# budget and otherwise dropped, and dropping more than .ae_boot_max_fail of them aborts

# replicate panels of a dynamic fit, drawn as ame_als_bootstrap() draws them: new
# additive effects (whole paths under dynamic_ab) around the fitted rest of the
# predictor, then the residual
.lda_boot_sims <- function(fit, R_boot, seed) {
	if (!is.null(seed)) {
		if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			old_seed <- get(".Random.seed", envir = globalenv())
			on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
		} else {
			on.exit(
				if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
					rm(".Random.seed", envir = globalenv()),
				add = TRUE
			)
		}
		set.seed(seed)
	}
	Tt <- fit$n_time
	nr <- fit$dims$n_row
	nc <- fit$dims$n_col
	rho_hat <- as.numeric(fit$VC["rho"])
	if (!is.finite(rho_hat)) rho_hat <- 0
	dyn_ab <- isTRUE(fit$dynamic_ab) && is.matrix(fit$a_dynamic) && ncol(fit$a_dynamic) == Tt
	paths <- if (dyn_ab) .ae_boot_ab_paths(fit)
	ab_orig <- if (!dyn_ab) outer(as.numeric(fit$a), as.numeric(fit$b), "+")
	Y <- lapply(seq_len(R_boot), function(b) {
		ab_new <- if (dyn_ab) paths$draw() else .ae_boot_draw_ab(fit)
		out <- array(NA_real_, dim = c(nr, nc, Tt))
		for (t in seq_len(Tt)) {
			ez <- if (dyn_ab) fit$EZ[[t]] - paths$orig[[t]] + ab_new[[t]] else fit$EZ[[t]] - ab_orig + ab_new
			ez[!is.finite(ez)] <- 0
			yt <- as.matrix(.ae_boot_draw_response(fit, ez, rho_hat, fit$s2))
			storage.mode(yt) <- "double"
			if (isTRUE(fit$symmetric)) {
				yt[lower.tri(yt)] <- 0
				yt <- yt + t(yt)
			}
			yt[is.na(fit$Y[, , t])] <- NA_real_
			out[, , t] <- yt
		}
		out
	})
	list(Y = Y)
}

# one replicate's refit at the fit's budget, then four times it; a convergence refusal
# counts as not converged, any other error as a refit error. `fit` is NULL when dropped
.lda_boot_refit <- function(args, budget) {
	last <- NULL
	for (it in c(budget, 4 * budget)) {
		args$max_iter <- it
		r <- tryCatch(
			suppressWarnings(do.call(lame_dynamic_als, args)),
			error = function(e) e
		)
		last <- r
		if (inherits(r, "error")) {
			if (!inherits(r, "lame_als_convergence")) {
				return(list(fit = NULL, reason = "error", last = r))
			}
			next
		}
		if (isTRUE(r$converged)) return(list(fit = r, reason = "", last = r))
	}
	list(fit = NULL, reason = "nonconverged", last = last)
}

# refits of the replicates `sims` of a dynamic fit; `mask_sim` puts a simulated panel
# in the form the fitter takes (observed cells, changing composition)
.lda_boot_run <- function(
	fit,
	sims,
	mask_sim,
	call_args,
	seed,
	R_boot,
	verbose = FALSE
) {
	budget <- call_args$max_iter
	max_fail <- floor(.ae_boot_max_fail * R_boot + 1e-9)
	refits <- vector("list", R_boot)
	fail_reason <- character(R_boot)
	failures <- vector("list", R_boot)
	# refit every replicate from the fitted trajectories so replicates stay at the same
	# optimum and rotation
	boot_start <- .ae_static_seed(fit)
	first_failed <- NULL
	for (b in seq_len(R_boot)) {
		args <- call_args
		args$Y <- mask_sim(sims$Y[[b]])
		args$seed <- (if (is.null(seed)) 6886 else seed) + b
		args$verbose <- FALSE
		args$bootstrap <- 0L
		args$stability <- "none"
		args$start <- boot_start
		args$start_jitter <- 0
		args$call <- NULL
		res <- .lda_boot_refit(args, budget)
		if (is.null(res$fit)) {
			fail_reason[b] <- res$reason
			failures[[b]] <- res$last
			if (is.null(first_failed)) {
				first_failed <- list(replicate = b, fit = res$last)
			}
			# stop as soon as the dropped share exceeds the cap
			if (sum(nzchar(fail_reason)) > max_fail) {
				.ae_boot_abort_failed(
					fail_reason,
					R_boot,
					max_fail,
					first_failed,
					budget,
					2L,
					b
				)
			}
		} else {
			refits[[b]] <- res$fit
		}
		if (isTRUE(verbose) && b %% 10L == 0L) {
			cli::cli_text("Dynamic bootstrap refit {b}/{R_boot}")
		}
	}
	kept <- !nzchar(fail_reason)
	valid <- refits[kept]
	coef_paths <- if (length(valid) > 0L) {
		array(
			unlist(lapply(valid, function(x) x$coef_path), use.names = FALSE),
			dim = c(nrow(fit$coef_path), ncol(fit$coef_path), length(valid)),
			dimnames = list(rownames(fit$coef_path), colnames(fit$coef_path), NULL)
		)
	} else NULL
	list(
		type = "parametric",
		R = R_boot,
		n_success = sum(kept),
		n_failed = sum(!kept),
		fail_reason = fail_reason,
		max_fail_share = .ae_boot_max_fail,
		coef_paths = coef_paths,
		converged = rep(TRUE, length(valid)),
		refits = valid,
		failures = failures[!kept]
	)
}
