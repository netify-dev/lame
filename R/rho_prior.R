# prior of the dyadic correlation rho, density proportional to (1 - rho^2)^power
# on (-1, 1); every update of rho in ame() and lame() reads it from here

# the exponent of (1 - rho^2) in the density of the prior
.lame_rho_prior_power <- function() 1

# note for panels with few mutual dyads: the absence of mutual ties is information too
.lame_rho_few_mutual_text <- paste(
	"Its posterior can stay near the prior, proportional to {.code 1 - rho^2} on (-1, 1), or move well away from it,",
	"toward -1 when the fitted means predict more mutual ties than the panel has; compare {.code fit$VC[, \"rho\"]} with that prior.",
	"{.code dcor = FALSE} fits a model without reciprocal residual dependence."
)

# one draw from the prior
.lame_rho_prior_draw <- function(power = .lame_rho_prior_power()) {
	2 * stats::rbeta(1L, power + 1, power + 1) - 1
}

# the log density of the prior, up to a constant
.lame_rho_prior_log <- function(rho, power = .lame_rho_prior_power()) {
	power * (log1p(rho) + log1p(-rho))
}

# smallest number of zero-contrast dyads from which the conditional of rho is improper
# (they contribute (1 - rho^2)^(-m / 2), integrable only when m / 2 - power < 1)
.lame_rho_improper_dyads <- function(power = .lame_rho_prior_power()) {
	max(1L, as.integer(ceiling(2 + 2 * power - 1e-9)))
}
