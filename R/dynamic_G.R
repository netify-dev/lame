# refuse dynamic_G = TRUE on every route (MCMC and ALS): with moving positions G_t is not
# identified apart from them, since U_t G_t is unconstrained; with static positions it is
# identified but has no sampler or estimator
.lame_refuse_dynamic_G <- function(dynamic_G) {
	if (!isTRUE(dynamic_G)) return(invisible(FALSE))
	cli::cli_abort(c(
		"{.arg dynamic_G = TRUE} is not supported.",
		"x" = "lame does not fit a per-period interaction matrix. With moving positions ({.code dynamic_uv = TRUE}) it is not identified apart from them, since {.code U_t G_t V_t'} and {.code (U_t G_t) V_t'} describe the same products; with static positions it is identified but not implemented.",
		"i" = "Use {.code dynamic_uv = TRUE} for a time-varying multiplicative term and read the per-period product {.code fit$UVPM}, which is the identified quantity.",
		"i" = "A single interaction matrix ({.code dynamic_G = FALSE}, the default) is still estimated for every bipartite fit and returned on it as {.code G}."
	))
}
