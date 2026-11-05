# units of a response whose scale is free (the normal family) and the default priors written in them.
# the other families' links fix the latent scale; for the normal family every default prior and starting value
# is scaled by the response, so a prior means the same for y and 1000 y. ame(), lame(), the ALS fits and the
# dynamic EB engines read the same numbers from these helpers and R/uv_prior.R.
# unit_var: variance of y, the unit of s2, coefficients, interaction weights and dynamic innovation scales;
# ab_scale: variance scale of the additive effects, the unit of Sab0 and (as its square root) of the factor rows.
# a new prior is (default for a unit response) x (unit helper of the same degree); floors are relative too

# variance of the observed response, the unit of s2 and of everything in units of y^2; 1 where the link
# fixes the scale. falls back to the mean square, then to 1 for an identically zero response
.lame_unit_var <- function(Y, family) {
	if (!identical(family, "normal")) return(1)
	y <- as.numeric(unlist(Y, use.names = FALSE))
	y <- y[is.finite(y)]
	v <- if (length(y) > 1L) stats::var(y) else NA_real_
	if (!is.finite(v) || v <= 0) v <- mean(y^2)
	if (!is.finite(v) || v <= 0) v <- 1
	v
}

# sd of the observed response, the unit of y itself
.lame_unit_sd <- function(Y, family) sqrt(.lame_unit_var(Y, family))

# variance scale of the normal-family additive effects: mean variance of the row and column means of the
# centred response (3-d slices pooled), or the response variance when those means do not vary
.lame_ab_scale <- function(Y) {
	Ec <- Y - mean(Y, na.rm = TRUE)
	if (length(dim(Ec)) == 2L) {
		a_tmp <- rowMeans(Ec, na.rm = TRUE); b_tmp <- colMeans(Ec, na.rm = TRUE)
	} else {
		a_tmp <- apply(Ec, 1, mean, na.rm = TRUE); b_tmp <- apply(Ec, 2, mean, na.rm = TRUE)
	}
	a_tmp[is.na(a_tmp)] <- 0; b_tmp[is.na(b_tmp)] <- 0
	s <- mean(c(stats::var(a_tmp), stats::var(b_tmp)), na.rm = TRUE)
	if (!is.finite(s) || s <= 0) s <- .lame_unit_var(Y, "normal")
	s
}

# prior variance of one entry of the bipartite G or one symmetric eigen-weight (the positions are unit free):
# the response variance for the normal family, 1 otherwise
.lame_weight_prior_var <- function(Y, family) .lame_unit_var(Y, family)

# the default prior variance of the initial dynamic coefficients, 10 in the units of the response
.lame_beta0_var_default <- function(Y, family) 10 * .lame_unit_var(Y, family)

# default inverse-gamma scale of the dynamic position innovation variance: the factor-row scale (R/uv_prior.R)
# for a unipartite normal network, whose u_i'v_j carries the units of y; 1 otherwise
.lame_sigma_uv_scale_default <- function(Y, family, bip = FALSE, symmetric = FALSE) {
	if (!identical(family, "normal") || isTRUE(bip) || isTRUE(symmetric)) return(1)
	.lame_uv_prior_scale(Y, family)
}

# multipliers of the unit-response safety rails for positions and G: unipartite coordinates scale with
# sqrt(sd(y)), bipartite and symmetric coordinates are unit free and G carries the unit
.lame_rail_units <- function(Y, family, bip = FALSE, symmetric = FALSE) {
	if (!identical(family, "normal")) return(list(uv = 1, G = 1))
	u <- .lame_unit_sd(Y, family)
	list(uv = if (isTRUE(bip) || isTRUE(symmetric)) 1 else sqrt(u), G = u)
}
