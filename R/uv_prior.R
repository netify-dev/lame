# the inverse-Wishart prior of the multiplicative effects of a unipartite network, read
# by every Suv update, every move that integrates Suv out and the empirical-Bayes ALS
# estimators. directed rows (u_i, v_i) ~ N(0, Suv) with Suv ~ IW(df, psi I), written in
# ame() as kappa0 = df and Suv0 = (psi / df) I (Hoff's convention in amen). Suv0 = s I
# with s the additive-effect scale of the data (probit moments for binary, ordinal,
# cbin and frn, log rates for poisson, a square root for normal so that u_i'v_j has the
# units of y), and df = round((2R + 2) m) with m the inverse Bernoulli information per
# cell (m = 1 for ordinal, poisson and normal). symmetric and bipartite factors have
# their own priors and are not covered here

# the smallest degrees of freedom of the inverse-Wishart prior of Suv with a finite mean, the df of the normal family
.lame_uv_prior_df <- function(R) 2L * as.integer(R) + 2L

# the unit of the factor rows: 1 when the link fixes the scale, the data scale for normal
# (in the units of y, as the product of the rows is)
.lame_uv_prior_scale <- function(Y, family) {
	if (!identical(family, "normal")) return(1)
	sqrt(.lame_ab_scale(Y))
}

# probit moments of a binary response (a matrix, or an array of periods pooled): vscale
# the variance of the actor means of the standardised residuals of the minority class,
# and vdfmlt the inverse Bernoulli information 0.25 / mean(p (1 - p)) of the probit fit
.lame_binary_probit_scale <- function(Y) {
	fallback <- list(vscale = 1, vdfmlt = 1)
	d <- dim(Y)
	if (length(d) == 2L) d <- c(d, 1L)
	if (length(d) != 3L) return(fallback)
	ydist <- table(Y)
	if (length(ydist) < 2L) return(fallback)
	ymode <- as.numeric(names(ydist)[ydist == max(ydist)])[1L]
	YB <- array(1 * (Y != ymode), d)
	ybar <- mean(YB, na.rm = TRUE)
	if (!is.finite(ybar) || ybar <= 0 || ybar >= 1) return(fallback)
	mu_bin <- stats::qnorm(ybar)
	E <- (YB - ybar) / stats::dnorm(stats::qnorm(ybar))
	if (d[1L] == d[2L]) for (t in seq_len(d[3L])) diag(E[, , t]) <- 0
	if (d[3L] == 1L) {
		a_tmp <- rowMeans(E[, , 1L], na.rm = TRUE); b_tmp <- colMeans(E[, , 1L], na.rm = TRUE)
	} else {
		a_tmp <- apply(E, 1L, mean, na.rm = TRUE); b_tmp <- apply(E, 2L, mean, na.rm = TRUE)
	}
	a_tmp[is.na(a_tmp)] <- 0; b_tmp[is.na(b_tmp)] <- 0
	vscale <- mean(diag(stats::cov(cbind(a_tmp, b_tmp))))
	if (!is.finite(vscale) || vscale <= 0) vscale <- 1
	PHAT <- stats::pnorm(mu_bin + outer(a_tmp, b_tmp, "+"))
	vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
	if (!is.finite(vdfmlt) || vdfmlt < 1) vdfmlt <- 1
	list(vscale = vscale, vdfmlt = vdfmlt)
}

# the data scale s and the multiplier m of the degrees of freedom of the factor prior of a family whose link fixes the
# scale (see the header); odmax the nomination cap of cbin and frn
.lame_uv_prior_moments <- function(Y, family, odmax = NULL) {
	switch(
		family,
		binary = .lame_binary_probit_scale(Y),
		ordinal = ,
		cbin = ,
		frn = {
			sc <- .lame_probit_ab_scale(Y, family, odmax)
			list(vscale = sc$vscale, vdfmlt = sc$vdfmlt)
		},
		poisson = list(vscale = .ame_pois_effect_scale(Y), vdfmlt = 1),
		list(vscale = 1, vdfmlt = 1)
	)
}

# the default prior of a directed network of rank R: df, psi = df * s, and ame()'s
# kappa0 = df and Suv0 = s I (NULL without a latent rank); `directed = FALSE` gives the
# unit prior the moves of a symmetric fit read
.lame_uv_prior <- function(Y, family, R, odmax = NULL, directed = TRUE) {
	if (is.null(R) || R < 1L) return(NULL)
	base <- .lame_uv_prior_df(R)
	if (identical(family, "normal") || !isTRUE(directed)) {
		df <- base
		s <- .lame_uv_prior_scale(Y, family)
		m <- 1
	} else {
		mo <- .lame_uv_prior_moments(Y, family, odmax)
		s <- mo$vscale
		m <- mo$vdfmlt
		df <- as.integer(round(base * m))
	}
	list(df = df, scale = s, psi = df * s, Suv0 = diag(s, 2L * R), kappa0 = df, df_multiplier = m)
}

# the prior in the form the moves take, list(df, scale) with scale the inverse-Wishart scale matrix
.lame_uv_prior_move <- function(prior, R) list(df = prior$df, scale = diag(prior$psi, 2L * R))

# the unit prior IW(2R + 2, (2R + 2) I) the moves use when the host passes none
.lame_uv_prior_unit <- function(R) {
	df <- .lame_uv_prior_df(R)
	list(df = df, scale = diag(df, 2L * R))
}

# a bipartite network has no row covariance to put an inverse-Wishart prior on (G carries
# the scale), so Suv0 and kappa0 are refused there, not ignored
.lame_refuse_bipartite_suv0 <- function(prior) {
	given <- intersect(c("Suv0", "kappa0"), names(prior))
	if (length(given)) {
		cli::cli_abort(
			c(
				"{.field {given}} of {.arg prior} {?does/do} not apply to a bipartite network.",
				"i" = "The rows of the factors of a bipartite network have a unit-free N(0, I) prior and the interaction matrix {.code G} carries the scale.",
				"i" = "Remove {.field {given}} from {.arg prior}."
			),
			class = "lame_prior_unused",
			given = given
		)
	}
	invisible(NULL)
}
