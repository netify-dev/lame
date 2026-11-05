# sum per-period likelihood statistics for the actor levels
.lame_ab_level_stats_period <- function(resid, dev, other, s2, margin, rho, skip_diag) {
	nt <- dim(resid)[3L]
	if (!is.numeric(rho) || length(rho) != nt || any(!is.finite(rho)) || any(abs(rho) >= 1))
		cli::cli_abort("Residual correlations must have one finite value in (-1, 1) per period.")
	if (all(rho == rho[1L])) return(.lame_ab_level_stats(resid, dev, other, s2, margin, rho[1L], skip_diag))
	if (any(!is.finite(dev)) || any(!is.finite(other))) cli::cli_abort("Dynamic actor paths must be finite.")
	states <- lapply(
		seq_len(nt),
		function(t)
			.lame_ab_level_stats(
				resid[, , t, drop = FALSE],
				dev[, t, drop = FALSE],
				other[, t, drop = FALSE],
				s2,
				margin,
				rho[t],
				skip_diag
			)
	)
	num <- Reduce(`+`, lapply(states, `[[`, "num"))
	prec <- Reduce(`+`, lapply(states, `[[`, "prec_d"))
	off <- Filter(Negate(is.null), lapply(states, `[[`, "prec_off"))
	off <- if (length(off)) Reduce(`+`, off) else NULL
	if (any(!is.finite(c(num, prec, off)))) cli::cli_abort("Nonfinite dynamic actor level precision.")
	list(num = num, prec_d = prec, prec_off = off)
}

# constant correlations fall back to the scalar transition
.lame_sample_dynamic_ab_period <- function(
	a_current,
	b_current,
	Z_array,
	EZ_array,
	rho_ab,
	sigma_ab,
	s2,
	symmetric,
	rho_dyad = 0,
	center = FALSE,
	rvar = TRUE,
	cvar = TRUE
) {
	nt <- ncol(a_current)
	if (!is.numeric(rho_dyad) || !length(rho_dyad) ||
			!(length(rho_dyad) %in% c(1L, nt)) || any(!is.finite(rho_dyad)) || any(abs(rho_dyad) >= 1))
		cli::cli_abort("Residual correlations must be scalar or one finite value in (-1, 1) per period.")
	if (length(rho_dyad) == 1L || all(rho_dyad == rho_dyad[1L]))
		return(sample_dynamic_ab_cpp(
			a_current,
			b_current,
			Z_array,
			EZ_array,
			rho_ab,
			sigma_ab,
			s2,
			symmetric,
			rho_dyad[1L],
			center,
			rvar,
			cvar
		))
	sample_dynamic_ab_rho_cpp(
		a_current,
		b_current,
		Z_array,
		EZ_array,
		rho_ab,
		sigma_ab,
		s2,
		symmetric,
		rho_dyad,
		center,
		rvar,
		cvar
	)
}
