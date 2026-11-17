# long-run means of the ar1 coefficient paths of lame(). an ar1 path reverts to its own mean,
#     beta_kt - mu_k = rho (beta_k,t-1 - mu_k) + e_kt,  e_kt ~ N(0, sigma^2 lam_k),
#     beta_k1 - mu_k ~ N(0, sigma^2 lam_k / (1 - rho^2)),  mu_k ~ N(0, v_k),
# with v_k the prior variance the coefficient has when static: in a unipartite fit the marginal of the
# panel prior on its column, free of s2 (.lame_static_panel_var()), in a bipartite one g_k s2 with g_k
# the scale of its g-prior (.lame_beta_mean_scale()). the path is drawn given mu_k as its
# deviation from mu_k, and mu_k given the path from its gaussian conditional. rw1, rw2 and matern32
# paths carry no mean, so fits without an ar1 block never call these helpers

# which dynamic coefficients follow a mean-reverting ar1 path
.lame_beta_ar1_mask <- function(group_id, kind_by_group) {
	unname(kind_by_group[group_id] == "ar1")
}

# prior scale g_k of each long-run mean of a bipartite fit, the scale of the static g-prior of that coefficient: g over
# the mean square of the design column for a default g (one for the intercept), else g as given
.lame_beta_mean_scale <- function(g, dyn_idx, XtX_per_t, n_dyads_per_t, g_default, intercept) {
	out <- .lame_static_g_scale(g, dyn_idx)
	if (!isTRUE(g_default)) return(out)
	n_cells <- max(sum(unlist(n_dyads_per_t)), 1)
	for (j in seq_along(dyn_idx)) {
		k <- dyn_idx[j]
		if (isTRUE(intercept) && k == 1L) next
		ms_k <- sum(vapply(XtX_per_t, function(m) m[k, k], numeric(1))) / n_cells
		if (is.finite(ms_k) && ms_k > 0) out[j] <- out[j] / ms_k
	}
	out
}

# covariance of the first state: stationary for an ar1 coefficient, the fixed prior otherwise
.lame_beta_init_cov <- function(beta0_cov, ar1, rho, sigma, lam) {
	if (!any(ar1)) return(beta0_cov)
	if (any(abs(rho[ar1]) >= 1)) stop("Mean-reverting coefficient paths need a persistence below one.")
	v <- diag(beta0_cov)
	v[ar1] <- sigma[ar1]^2 * lam[ar1] / (1 - rho[ar1]^2)
	diag(v, length(v))
}

# mean of the first state of the deviation path: zero for an ar1 coefficient
.lame_beta_init_mean <- function(beta0_mean, ar1) {
	beta0_mean[ar1] <- 0
	beta0_mean
}

# paths as deviations from their long-run means (columns without a mean have mu = 0), and back
.lame_beta_deviation <- function(path, mu) {
	path - rep(mu, each = nrow(path))
}
.lame_beta_add_mean <- function(dev, mu) {
	dev + rep(mu, each = nrow(dev))
}

# offsets that carry the long-run means, so the path sampler draws the deviations
.lame_beta_mean_offsets <- function(offset_list, Xdyn_list, mu) {
	lapply(seq_along(offset_list), function(t) {
		off <- offset_list[[t]]
		off + matrix(drop(Xdyn_list[[t]] %*% mu), nrow(off), ncol(off))
	})
}

# per-period statistics of the deviation from those of the path: h_t - Omega_t mu
.lame_beta_mean_stats <- function(st, mu) {
	p <- length(mu)
	h <- st$h
	for (t in seq_len(nrow(h))) h[t, ] <- h[t, ] - drop(matrix(st$Omega[, , t], p, p) %*% mu)
	list(Omega = st$Omega, h = h)
}

# gaussian conditional draw of each long-run mean given its path, (rho, sigma) and its prior
# variance; `gaps` the time gaps between periods
.lame_beta_mean_draw <- function(path, mu, ar1, rho, sigma, lam, prior_var, gaps) {
	nt <- nrow(path)
	for (k in which(ar1)) {
		z <- .lame_beta_gap_transition(rho[k], gaps)
		q <- sigma[k]^2 * lam[k]
		x <- path[, k]
		e <- x[-1L] - z$f * x[-nt]
		prec <- ((1 - rho[k]^2) + sum((1 - z$f)^2 / z$q)) / q + 1 / prior_var[k]
		lin <- ((1 - rho[k]^2) * x[1L] + sum((1 - z$f) * e / z$q)) / q
		if (!is.finite(prec) || prec <= 0 || !is.finite(lin)) {
			stop("Invalid conditional of a coefficient long-run mean.")
		}
		mu[k] <- stats::rnorm(1L, lin / prec, 1 / sqrt(prec))
	}
	mu
}

# persistence of each ar1 block given its deviation paths, on equally spaced periods. the transitions
# and the gaussian part of the stationary first state give a truncated normal, which is the proposal of an
# independence step that corrects for the remaining (1 - rho^2)^(n / 2) of the first states
.lame_sample_rho_beta_stationary <- function(
	dev,
	group_id,
	n_groups,
	lam,
	sigma_by_coef,
	rho_current,
	rho_prior_mean,
	rho_prior_sd,
	rho_lower,
	rho_upper,
	stationary
) {
	nt <- nrow(dev)
	for (g in which(stationary)) {
		idx <- which(group_id == g)
		if (!length(idx) || nt < 2L) next
		x <- dev[, idx, drop = FALSE]
		w <- 1 / (sigma_by_coef[idx]^2 * lam[idx])
		inner <- x[-c(1L, nt), , drop = FALSE]
		A <- sum(colSums(inner^2) * w)
		B <- sum(colSums(x[-nt, , drop = FALSE] * x[-1L, , drop = FALSE]) * w)
		prec <- A + 1 / rho_prior_sd[g]^2
		m <- (B + rho_prior_mean[g] / rho_prior_sd[g]^2) / prec
		prop <- .lame_ar_truncnorm_draw(m, 1 / sqrt(prec), rho_lower, rho_upper)
		lr <- length(idx) / 2 * (log1p(-prop^2) - log1p(-rho_current[g]^2))
		if (is.finite(lr) && log(stats::runif(1L)) < lr) rho_current[g] <- prop
	}
	rho_current
}

# log prior terms of the long-run means for the scale moves, which scale them with the paths:
# quadratic form and dimension
.lame_beta_mean_terms <- function(mu, ar1, prior_var) {
	c(q = sum(mu[ar1]^2 / prior_var[ar1]), d = sum(ar1))
}

# long-run means in the stored coefficient layout: `to_store` maps each design column to its stored
# column (a symmetric fit stores one coefficient for the row and column slices of a nodal covariate)
.lame_beta_mean_store <- function(mu, ar1, dyn_idx, to_store, store_cols) {
	full <- rep(NA_real_, length(to_store))
	full[dyn_idx[ar1]] <- mu[ar1]
	vapply(store_cols, function(j) mean(full[to_store == j], na.rm = TRUE), numeric(1))
}

# long-run mean of every coefficient of a fit for draw i: the stored draw for the ar1 coefficients,
# zero for the others (which follow their path without a mean)
.lame_beta_mean_row <- function(object, i, coef_names) {
	out <- stats::setNames(numeric(length(coef_names)), coef_names)
	M <- object$MU_BETA
	if (is.null(M) || !ncol(M)) return(out)
	out[colnames(M)] <- M[i, ]
	out
}
