# transition state that survives a checkpoint, excluding posterior accumulators
.lame_checkpoint_fields <- function() {
	c(
		"Z", "beta", "a", "b", "a_mat", "b_mat", "a_lev", "b_lev",
		"U", "V", "U_cube", "V_cube", "L_dyn", "G", "rho", "s2", "Sab",
		"rho_ab", "sigma_ab", "rho_uv", "sigma_uv", "sigma2_a", "sigma2_b",
		"beta_dyn_path", "beta_static", "rho_beta_by_group", "sigma_beta_by_group",
		"rho_by_coef", "sigma_by_coef", "pool_rho_mu", "pool_rho_tau", "pi_snap", "delta_uv", "delta_uv_v",
		"lambda_u", "lambda_v", "alpha", "tau_prop",
		"theta_actor", "theta_actor_prev", "rho_actor", "sigma_actor_sq",
		"lev_sd", "dev_sd", "marginal_sd", "dhmc_sd", "mh_sd_U", "mh_sd_V",
		"uv_actor_path_version", "ab_hyper_par"
	)
}

.lame_checkpoint_capture <- function(envir, layout) {
	fields <- .lame_checkpoint_fields()
	fields <- fields[vapply(
		fields,
		exists,
		logical(1),
		envir = envir,
		inherits = FALSE
	)]
	list(
		version = 1L,
		layout = layout,
		values = mget(fields, envir = envir, inherits = FALSE)
	)
}

.lame_checkpoint_restore <- function(state, envir, layout) {
	if (!is.list(state) || !identical(state$version, 1L) ||
			!is.list(state$values) || is.null(names(state$values)) ||
			anyDuplicated(names(state$values)) ||
			!all(names(state$values) %in% .lame_checkpoint_fields())) {
		cli::cli_abort(
			"The checkpoint carries an invalid sampler state.",
			class = "lame_checkpoint_state"
		)
	}
	if (!identical(state$layout, layout)) {
		cli::cli_abort(
			c(
				"The checkpoint's model layout differs from this call.",
				"i" = "Resume with the original actors, periods and model options. For a different model, remove {.code start_vals$checkpoint_state} and supply a new burn-in."
			),
			class = "lame_checkpoint_layout"
		)
	}
	# tuning fields a checkpoint may lack keep their initialized values
	required <- setdiff(
		names(.lame_checkpoint_capture(envir, layout)$values),
		c("marginal_sd", "dhmc_sd", "uv_actor_path_version", "ab_hyper_par")
	)
	if (!all(required %in% names(state$values)) ||
			!all(vapply(state$values, function(x) is.null(x) || is.numeric(x), logical(1)))) {
		cli::cli_abort(
			"The checkpoint carries incomplete or invalid transition variables.",
			class = "lame_checkpoint_state"
		)
	}
	for (nm in names(state$values)) assign(nm, state$values[[nm]], envir = envir)
	invisible(NULL)
}
