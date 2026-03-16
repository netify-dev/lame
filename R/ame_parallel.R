#' Run AME model with multiple parallel chains
#'
#' @param Y Network data matrix
#' @param n_chains Number of parallel chains to run (default = 4)
#' @param cores Number of CPU cores to use for parallel processing.
#'   Default is n_chains. Use 1 for sequential processing.
#' @param combine_method Method for combining chains: "pool" (default) or "list"
#' @param ... Additional arguments passed to ame()
#'
#' @return If combine_method = "pool": A single ame object with pooled chains
#'         If combine_method = "list": A list of ame objects, one per chain
#'
#' @examples
#' \donttest{
#' # Run 2 chains sequentially
#' data(YX_nrm)
#' fit_parallel <- ame_parallel(YX_nrm$Y, Xdyad = YX_nrm$X,
#'                              n_chains = 2, cores = 1,
#'                              nscan = 100, burn = 10, odens = 1,
#'                              verbose = FALSE)
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
ame_parallel <- function(
	Y, n_chains = 4, cores = n_chains,
	combine_method = c("pool", "list"), ...
	){

	combine_method <- match.arg(combine_method)

	if(cores > 1 && !requireNamespace("parallel", quietly = TRUE)) {
		cli::cli_warn("Package 'parallel' not available. Running chains sequentially.")
		cores <- 1
	}

	ame_args <- list(...)

	####
	# generate seeds and define chain runner
	base_seed <- ame_args$seed %||% 6886
	chain_seeds <- base_seed + seq_len(n_chains) * 1000

	run_chain <- function(chain_id) {
		cli::cli_text("Starting chain {.val {chain_id}}")

		ame_args$seed <- chain_seeds[chain_id]

		if(!is.null(ame_args$out_file)) {
			ame_args$out_file <- gsub("\\.RData$",
																paste0("_chain", chain_id, ".RData"),
																ame_args$out_file)
		}

		fit <- do.call(ame, c(list(Y = Y), ame_args))
		fit$chain_id <- chain_id

		cli::cli_text("Completed chain {.val {chain_id}}")

		return(fit)
	}
	####

	####
	# run chains
	if(cores > 1) {
		cli::cli_h2("Running {.val {n_chains}} chains in parallel on {.val {cores}} cores")

		cl <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cl))

		parallel::clusterEvalQ(cl, {
			library(lame)
		})

		parallel::clusterExport(cl, c("Y", "ame_args", "run_chain"), envir = environment())

		chain_results <- parallel::parLapply(cl, seq_len(n_chains), run_chain)

	} else {
		cli::cli_h2("Running {.val {n_chains}} chains sequentially")

		chain_results <- lapply(seq_len(n_chains), run_chain)
	}
	####

	####
	# combine results
	if(combine_method == "pool") {
		cli::cli_text("Combining chains...")
		combined_fit <- combine_ame_chains(chain_results)
		return(combined_fit)
	} else {
		return(chain_results)
	}
	####
}

#' Combine multiple AME chains
#'
#' @param chain_list List of ame fit objects from multiple chains
#' @param diagnostics Logical; whether to compute convergence diagnostics (default TRUE)
#'
#' @return A single ame object with combined chains
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
combine_ame_chains <- function(chain_list, diagnostics = TRUE) {

	n_chains <- length(chain_list)

	if(n_chains == 0) {
		cli::cli_abort("No chains provided")
	}

	if(n_chains == 1) {
		cli::cli_warn("Only one chain provided, returning as-is")
		return(chain_list[[1]])
	}

	combined <- chain_list[[1]]

	####
	# combine BETA and VC matrices
	if(!is.null(combined$BETA)) {
		BETA_list <- lapply(chain_list, function(x) x$BETA)
		combined$BETA <- do.call(rbind, BETA_list)

		n_samples <- sapply(BETA_list, nrow)
		combined$chain_indicator <- rep(seq_len(n_chains), n_samples)
	}

	if(!is.null(combined$VC)) {
		VC_list <- lapply(chain_list, function(x) x$VC)
		combined$VC <- do.call(rbind, VC_list)
	}

	if(!is.null(combined$GOF)) {
		GOF_list <- lapply(chain_list, function(x) x$GOF)
		combined$GOF <- do.call(rbind, GOF_list)
	}
	####

	####
	# combine posterior samples if saved
	if(!is.null(combined$U_samples)) {
		U_samples_list <- lapply(chain_list, function(x) x$U_samples)
		combined$U_samples <- abind::abind(U_samples_list, along = 3)
	}

	if(!is.null(combined$V_samples)) {
		V_samples_list <- lapply(chain_list, function(x) x$V_samples)
		combined$V_samples <- abind::abind(V_samples_list, along = 3)
	}

	if(!is.null(combined$a_samples)) {
		a_samples_list <- lapply(chain_list, function(x) x$a_samples)
		combined$a_samples <- do.call(cbind, a_samples_list)
	}

	if(!is.null(combined$b_samples)) {
		b_samples_list <- lapply(chain_list, function(x) x$b_samples)
		combined$b_samples <- do.call(cbind, b_samples_list)
	}
	####

	####
	# average posterior means across chains
	if(!is.null(combined$U)) {
		U_list <- lapply(chain_list, function(x) x$U)
		combined$U <- Reduce("+", U_list) / n_chains
	}

	if(!is.null(combined$V)) {
		V_list <- lapply(chain_list, function(x) x$V)
		combined$V <- Reduce("+", V_list) / n_chains
	}

	if(!is.null(combined$UVPM)) {
		UVPM_list <- lapply(chain_list, function(x) x$UVPM)
		combined$UVPM <- Reduce("+", UVPM_list) / n_chains
	}

	if(!is.null(combined$APM)) {
		APM_list <- lapply(chain_list, function(x) x$APM)
		combined$APM <- Reduce("+", APM_list) / n_chains
	}

	if(!is.null(combined$BPM)) {
		BPM_list <- lapply(chain_list, function(x) x$BPM)
		combined$BPM <- Reduce("+", BPM_list) / n_chains
	}

	if(!is.null(combined$YPM)) {
		YPM_list <- lapply(chain_list, function(x) x$YPM)
		combined$YPM <- Reduce("+", YPM_list) / n_chains
	}

	if(!is.null(combined$EZ)) {
		EZ_list <- lapply(chain_list, function(x) x$EZ)
		combined$EZ <- Reduce("+", EZ_list) / n_chains
	}
	####

	####
	# convergence diagnostics
	if(diagnostics && !is.null(combined$BETA)) {
		combined$diagnostics <- compute_mcmc_diagnostics(chain_list)
	}

	combined$n_chains <- n_chains
	combined$is_combined <- TRUE
	####

	return(combined)
}

#' Compute MCMC convergence diagnostics for multiple chains
#'
#' @param chain_list List of ame fit objects from multiple chains
#'
#' @return List containing convergence diagnostics
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
compute_mcmc_diagnostics <- function(chain_list) {

	n_chains <- length(chain_list)

	if(n_chains < 2) {
		cli::cli_warn("Need at least 2 chains for convergence diagnostics")
		return(NULL)
	}

	BETA_list <- lapply(chain_list, function(x) x$BETA)

	n_params <- ncol(BETA_list[[1]])
	n_samples <- sapply(BETA_list, nrow)

	####
	# compute R-hat and effective sample size
	rhat <- numeric(n_params)
	ess <- numeric(n_params)

	for(p in 1:n_params) {
		chains_p <- lapply(BETA_list, function(x) x[, p])

		# within-chain variance
		W <- mean(sapply(chains_p, var))

		# between-chain variance
		chain_means <- sapply(chains_p, mean)
		B <- var(chain_means) * min(n_samples)

		# r-hat
		var_plus <- ((min(n_samples) - 1) * W + B) / min(n_samples)
		rhat[p] <- sqrt(var_plus / W)

		# effective sample size
		if(requireNamespace("coda", quietly = TRUE)) {
			mcmc_list <- coda::mcmc.list(lapply(chains_p, coda::mcmc))
			ess[p] <- coda::effectiveSize(mcmc_list)
		}
	}
	####

	####
	# assemble and report diagnostics
	param_names <- colnames(BETA_list[[1]])
	if(is.null(param_names)) {
		param_names <- paste0("beta[", 1:n_params, "]")
	}

	diagnostics <- list(
		rhat = rhat,
		ess = ess,
		param_names = param_names,
		n_chains = n_chains,
		n_samples = n_samples
	)

	cli::cli_h3("MCMC Convergence Diagnostics")
	cli::cli_text("Number of chains: {.val {n_chains}}")
	cli::cli_text("Samples per chain: {.val {paste(n_samples, collapse = ', ')}}")

	converged <- all(rhat < 1.1, na.rm = TRUE)
	if(converged) {
		cli::cli_alert_success("All parameters converged (R-hat < 1.1)")
	} else {
		n_not_converged <- sum(rhat >= 1.1, na.rm = TRUE)
		cli::cli_alert_warning("{.val {n_not_converged}} parameters have R-hat >= 1.1")

		problem_params <- which(rhat >= 1.1)
		for(i in problem_params[1:min(5, length(problem_params))]) {
			cli::cli_text("  {param_names[i]}: R-hat = {.val {round(rhat[i], 3)}}")
		}
		if(length(problem_params) > 5) {
			cli::cli_text("  ... and {.val {length(problem_params) - 5}} more")
		}
	}
	####

	return(diagnostics)
}

#' Helper function to handle NULL with default
#' @noRd
`%||%` <- function(x, y) {
	if(is.null(x)) y else x
}
