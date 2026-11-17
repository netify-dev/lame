# prior of the innovation variance sigma_uv^2 of the dynamic latent positions (?lame, prior).
# the default is an inverse gamma on the stationary variance w = sigma_uv^2 / (1 - rho_uv^2) of a
# coordinate path, independent of rho_uv, so that the persistence keeps its own prior at every w.
# given rho_uv it is the inverse gamma IG(a, b (1 - rho_uv^2)) on sigma_uv^2, which the updates at a
# fixed persistence read through .lame_uv_ig_scale(); the updates that move rho_uv take `stationary`.
# a prior that sets sigma_uv_shape or sigma_uv_scale is the inverse gamma on sigma_uv^2 itself

# resolve the prior entries; `scale_default` is the unit of the positions (.lame_sigma_uv_scale_default).
# every kind takes the stationary prior unless sigma_uv_shape or sigma_uv_scale is set
.lame_uv_var_prior_resolve <- function(prior, scale_default, kind = "ar1") {
	innovation <- !is.null(prior$sigma_uv_shape) || !is.null(prior$sigma_uv_scale)
	stationary <- !is.null(prior$uv_var_shape) || !is.null(prior$uv_var_scale)
	if (innovation && stationary) {
		cli::cli_abort(
			c(
				"{.code prior} sets both the innovation prior and the stationary-variance prior of the positions.",
				"i" = "Set {.code sigma_uv_shape} and {.code sigma_uv_scale} or {.code uv_var_shape} and {.code uv_var_scale}, not both."
			),
			class = "lame_uv_var_prior"
		)
	}
	if (innovation) {
		if (is.null(prior$sigma_uv_shape)) prior$sigma_uv_shape <- 2
		if (is.null(prior$sigma_uv_scale)) prior$sigma_uv_scale <- scale_default
		return(prior)
	}
	if (is.null(prior$uv_var_shape)) prior$uv_var_shape <- 2
	if (is.null(prior$uv_var_scale)) prior$uv_var_scale <- 0.5 * scale_default
	for (nm in c("uv_var_shape", "uv_var_scale")) {
		x <- prior[[nm]]
		if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
			cli::cli_abort(
				"{.code prior${nm}} must be a single positive number.",
				class = "lame_uv_var_prior"
			)
		}
	}
	prior
}

# TRUE when the prior is on the stationary variance
.lame_uv_var_stationary <- function(prior) !is.null(prior$uv_var_shape)

# shape of the inverse gamma on sigma_uv^2 given rho_uv
.lame_uv_ig_shape <- function(prior) {
	if (.lame_uv_var_stationary(prior)) prior$uv_var_shape else prior$sigma_uv_shape
}

# scale of the inverse gamma on sigma_uv^2 given rho_uv
.lame_uv_ig_scale <- function(prior, rho) {
	if (.lame_uv_var_stationary(prior)) prior$uv_var_scale * (1 - rho * rho) else prior$sigma_uv_scale
}

# scale of the prior itself: on w when stationary, on sigma_uv^2 otherwise
.lame_uv_var_scale <- function(prior) {
	if (.lame_uv_var_stationary(prior)) prior$uv_var_scale else prior$sigma_uv_scale
}
