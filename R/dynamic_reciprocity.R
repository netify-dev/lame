# period-specific reciprocity correlations within the shared coefficient and actor kernel
.lame_dynamic_rho_beta_ab_draw <- function(
	ZT,
	mX,
	mXt,
	XX,
	Sab,
	rho,
	s2,
	g,
	center = FALSE
) {
	n <- dim(ZT)[1L]; nt <- dim(ZT)[3L]; p <- dim(mX)[2L]
	if (is.null(nt) || dim(ZT)[2L] != n ||
			!identical(dim(mX), c(n * n, p, nt)) ||
			!identical(dim(XX), c(p, p, nt)) ||
			any(!is.finite(ZT)) || any(!is.finite(XX)) ||
			length(g) != 1L || !is.finite(g) || g <= 0)
		cli::cli_abort("Invalid complete-data replicated Gaussian design or prior scale.")
	.lame_replicated_gaussian_draw(ZT, mX, Sab, rho, s2, g, center, XX)
}
