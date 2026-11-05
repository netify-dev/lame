# dyadic correlation of the dynamic binary EB estimator (directed, probit, dcor = TRUE).
#
# the two latent errors of a dyad are bivariate normal with correlation rho, so
#   P(y_ij, y_ji | eta_ij, eta_ji, rho) = Phi2(s_ij eta_ij, s_ji eta_ji; s_ij s_ji rho) ,  s = 2 y - 1 .
# a dyad with one cell missing keeps the probit probability of the observed cell.
#
# this file is the cell problem of such a pair under the noise of the product of
# the positions (the pair analogue of .lde_work_gh()): the expected score, the
# expected Fisher information (I_11, I_22, cross weight I_12) with its
# derivatives, and the expected log probability. pairs are integrated on the
# shared grid (lame_dynamic_als_eb_pairgrid.R) where the noise is wide and by the
# tensor Gauss-Hermite rule of R/ame_als_pair.R otherwise. at rho = 0 the pair
# reduces to two independent cells.
#
# rho is a hyperparameter of the Laplace-marginal likelihood, stepped by
# .ae_pair_rho_next() with the Laplace-marginal score from .lde_laplace().

# the pairs of a panel

# pairs of a directed panel (na x na x Tt arrays of the reduced actor set): dyads
# i < j of every period with both cells observed, as linear indices of their two
# cells, with actors, period and outcome signs
.lde_pair_context <- function(y) {
	n <- dim(y)[1L]
	Tt <- dim(y)[3L]
	ut <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
	i <- ut[, 1L]
	j <- ut[, 2L]
	off <- (seq_len(Tt) - 1L) * n * n
	one <- as.vector(outer(i + (j - 1L) * n, off, "+"))
	two <- as.vector(outer(j + (i - 1L) * n, off, "+"))
	keep <- is.finite(y[one]) & is.finite(y[two])
	list(
		one = one[keep],
		two = two[keep],
		i = rep(i, Tt)[keep],
		j = rep(j, Tt)[keep],
		t = rep(seq_len(Tt), each = length(i))[keep],
		s1 = 2 * y[one[keep]] - 1,
		s2 = 2 * y[two[keep]] - 1,
		n = n,
		Tt = Tt,
		n_pairs = sum(keep)
	)
}

# expectation under the noise of the product

# nodes a side for a pair with product variances v1, v2: one node without noise,
# else the smallest rung with at least `per_sd` sqrt(max v) nodes
.lde_pair_rule <- function(v1, v2, per_sd = 9) {
	vm <- pmax(v1, v2)
	rungs <- c(5L, 9L, 15L, 21L, 31L, 41L, 61L)
	k <- rungs[pmin(findInterval(sqrt(vm) * per_sd, rungs, left.open = TRUE) + 1L, length(rungs))]
	k[!(vm > 1e-12)] <- 1L
	k
}

# expected pair problem under q for independent predictors of the two cells: per
# pair the expected log probability ell, score (psi1, psi2), Fisher information
# (I11, I22, I12), with `deriv` its derivatives in the two means and with
# `rho_deriv` the derivatives in rho, and the nodes a side `K`. pairs the shared
# grid resolves are integrated on it, the rest by the direct kernel, and with the
# grid off (options(lame.pair_grid = FALSE)) or all sds below 0.08 by the tensor
# Gauss-Hermite rule of .ae_pair_quad_fisher(). `tier` 1 is the coarse rule of
# the inner fits of the position-scale search, 2 the rule of the fixed point
.lde_pair_quad <- function(m1, m2, v1, v2, s1, s2, rho, deriv = TRUE, rho_deriv = FALSE, tier = 2L) {
	n <- length(m1)
	rho_deriv <- deriv && rho_deriv
	opt <- getOption("lame.pair_grid", TRUE)
	plan <- NULL
	if (!identical(opt, FALSE) && n > 0L) {
		sd1 <- sqrt(pmax(v1, 0))
		sd2 <- sqrt(pmax(v2, 0))
		variant <- 2L * (s1 < 0) + (s2 < 0)
		ncol <- 6L + (if (deriv) 6L else 0L) + (if (rho_deriv) 7L else 0L)
		plan <- .lde_pair_plan(m1, m2, sd1, sd2, variant, rho, ncol, tier, force = identical(opt, "force"))
	}
	if (is.null(plan)) {
		kk <- .lde_pair_rule(v1, v2)
		out <- .ae_pair_quad_fisher(m1, m2, v1, v2, 0, s1, s2, rho, kk, deriv, rho_deriv)
		out$K <- kk
		return(out)
	}
	nm <- c(
		"ell",
		"psi1",
		"psi2",
		"I11",
		"I22",
		"I12",
		if (deriv) c("I11a", "I11b", "I22a", "I22b", "I12a", "I12b"),
		if (rho_deriv) c("s0", "rc0", "psi1r", "psi2r", "I11r", "I22r", "I12r")
	)
	out <- stats::setNames(lapply(nm, function(f) numeric(n)), nm)
	K <- integer(n)
	r_engine <- identical(.ae_engine(), "r")
	threads <- .ae_threads()
	todo <- plan$direct
	ix <- which(plan$tab)
	if (length(ix)) {
		tab <- .lde_pair_table(plan, rho, plan$vmask, threads, tier)
		ap <- if (r_engine) {
			.lde_pair_apply_r(tab, m1[ix], m2[ix], sd1[ix], sd2[ix], variant[ix], plan$cutoff, deriv, rho_deriv)
		} else {
			lde_pair_apply_cpp(
				tab,
				m1[ix],
				m2[ix],
				sd1[ix],
				sd2[ix],
				variant[ix],
				plan$cutoff,
				deriv,
				rho_deriv,
				threads
			)
		}
		good <- ap$ok == 1L
		for (f in nm) out[[f]][ix[good]] <- ap[[f]][good]
		K[ix[good]] <- ap$K[good]
		# a pair whose window the table did not cover is integrated directly
		if (!all(good)) {
			bad <- ix[!good]
			todo <- c(todo, bad)
			plan$type1 <- c(plan$type1, rep(0L, length(bad))); plan$type2 <- c(plan$type2, rep(0L, length(bad)))
			plan$K1 <- c(plan$K1, .lde_gh_rung(sd1[bad])); plan$K2 <- c(plan$K2, .lde_gh_rung(sd2[bad]))
		}
	}
	if (length(todo)) {
		dr <- if (r_engine) {
			.lde_pair_direct_r(
				m1[todo],
				m2[todo],
				sd1[todo],
				sd2[todo],
				variant[todo],
				plan$type1,
				plan$type2,
				plan$K1,
				plan$K2,
				rho,
				plan$cc,
				plan$du,
				plan$cutoff,
				deriv,
				rho_deriv
			)
		} else {
			lde_pair_direct_cpp(
				m1[todo],
				m2[todo],
				sd1[todo],
				sd2[todo],
				variant[todo],
				plan$type1,
				plan$type2,
				plan$K1,
				plan$K2,
				rho,
				plan$cc,
				plan$du,
				plan$cutoff,
				deriv,
				rho_deriv,
				threads
			)
		}
		okd <- dr$ok == 1L
		for (f in nm) out[[f]][todo[okd]] <- dr[[f]][okd]
		K[todo[okd]] <- dr$K[okd]
		if (!all(okd)) {
			# unreachable in practice: the shared kernel is the last resort
			bad <- todo[!okd]
			kk <- .lde_pair_rule(v1[bad], v2[bad])
			gh <- .ae_pair_quad_fisher(
				m1[bad],
				m2[bad],
				v1[bad],
				v2[bad],
				0,
				s1[bad],
				s2[bad],
				rho,
				kk,
				deriv,
				rho_deriv
			)
			for (f in nm) out[[f]][bad] <- gh[[f]]
			K[bad] <- kk
		}
	}
	out$K <- K
	out
}

# twice the negative log probability of the observed pairs with the noise
# integrated out of the likelihood: Phi2 of the means over sqrt(1 + v) with
# correlation r / sqrt((1 + v1)(1 + v2)) (used by the position-scale profile)
.lde_pair_dev_marg <- function(m1, m2, v1, v2, s1, s2, rho) {
	sg1 <- sqrt(1 + v1)
	sg2 <- sqrt(1 + v2)
	-2 * .ae_pair_logphi2(s1 * m1 / sg1, s2 * m2 / sg2, s1 * s2 * rho / (sg1 * sg2))
}

# working problem of the cells of a panel with the dyadic correlation

# expected working problem of every cell when the two cells of a dyad are a
# pair: the arrays of .lde_work_gh() / .lde_work(), the cross weight `cw` (equal
# in both cells of a dyad, zero without a partner) and, for the Laplace
# correction, the pair record `pair`. `vnoise` NULL: no noise. a cell without an
# observed partner is an independent cell.
.lde_work_pair <- function(ctx, y, eta, vnoise, rho, what = 1L, link = "probit") {
	pc <- ctx$pc
	one <- pc$one
	two <- pc$two
	yy <- y
	yy[c(one, two)] <- NA_real_
	wk <- if (is.null(vnoise)) .lde_work(yy, eta, link) else .lde_work_gh(yy, eta, vnoise, what = what, link = link)
	if (!length(one)) {
		wk$cw <- array(0, dim(y))
		wk$pair <- NULL
		return(wk)
	}
	deriv <- what >= 1L
	m1 <- eta[one]
	m2 <- eta[two]
	v1 <- if (is.null(vnoise)) numeric(length(one)) else vnoise[one]
	v2 <- if (is.null(vnoise)) numeric(length(one)) else vnoise[two]
	# derivatives in rho are needed only where rho takes its step
	q <- .lde_pair_quad(
		m1,
		m2,
		v1,
		v2,
		pc$s1,
		pc$s2,
		rho,
		deriv,
		rho_deriv = deriv && is.null(ctx$rho_hold),
		tier = if (is.null(ctx$pair_tier)) 2L else ctx$pair_tier
	)
	# the working-response increment is capped at 30 as for a cell
	cap <- function(psi, w) {
		w <- pmax(w, 1e-8)
		big <- abs(psi / w) > 30
		if (any(big)) w[big] <- abs(psi[big]) / 30
		w
	}
	w1 <- cap(q$psi1, q$I11)
	w2 <- cap(q$psi2, q$I22)
	wk$z[one] <- m1 + q$psi1 / w1
	wk$z[two] <- m2 + q$psi2 / w2
	wk$w[one] <- w1
	wk$w[two] <- w2
	wk$score[one] <- q$psi1
	wk$score[two] <- q$psi2
	cw <- array(0, dim(y))
	cw[one] <- q$I12
	cw[two] <- q$I12
	wk$cw <- cw
	if (deriv) {
		wk$hprime[one] <- q$I11a
		wk$hprime[two] <- q$I22b
		# the curvature of the weight only damps the coefficient step
		# (.lde_step_factor()), so the independent cell's curvature stands in
		yp <- y
		yp[-c(one, two)] <- NA_real_
		wi <- if (is.null(vnoise)) .lde_work(yp, eta, link) else .lde_work_gh(yp, eta, vnoise, what = what, link = link)
		wk$hpp[one] <- wi$hpp[one]
		wk$hpp[two] <- wi$hpp[two]
		wk$hprime[!is.finite(wk$hprime)] <- 0
		if (!is.null(wk$dev_marg) && is.finite(wk$dev_marg)) wk$dev_marg <- wk$dev_marg + sum(
			.lde_pair_dev_marg(m1, m2, v1, v2, pc$s1, pc$s2, rho)
		)
	}
	# data deviance: expected log probability (with noise) or log probability at the means
	if (!is.null(wk$dev) && is.finite(wk$dev)) wk$dev <- wk$dev - 2 * sum(q$ell)
	wk$pair <- c(
		pc[c("one", "two", "i", "j", "t")],
		q[c("ell", "K")],
		if (deriv) q[c("I11a", "I11b", "I22a", "I22b", "I12a", "I12b")],
		if (deriv && is.null(ctx$rho_hold)) q[c("s0", "rc0", "psi1r", "psi2r", "I11r", "I22r", "I12r")]
	)
	wk
}

# rho (update and posterior sd are the static estimator's .ae_pair_rho_next()
# and .ae_pair_rho_sd(), with score `rs` and curvature `rc` from the quadrature)

# working problem of the additive predictor alone: the cells, or the pairs of
# cells under a dyadic correlation; used on the R = 0 path and for reported deviances
.lde_work_plain <- function(ctx, eta, hyp) {
	if (isTRUE(ctx$nrm)) return(.lde_work_nrm(ctx, eta, hyp))
	if (isTRUE(ctx$pln)) return(.lde_work_pln(ctx, eta, hyp$ve))
	if (isTRUE(ctx$dcor)) return(.lde_work_pair(ctx, ctx$y, eta, NULL, hyp$rho_d, 1L, ctx$link))
	.lde_work(ctx$y, eta, ctx$link)
}

# deviance of the additive predictor alone: the cells' deviance, and for pairs
# -2 times the log bivariate probit probability of the observed pair
.lde_dev_plain <- function(ctx, eta, hyp) {
	if (isTRUE(ctx$nrm)) return(.lde_dev_nrm(ctx, eta, hyp))
	if (isTRUE(ctx$pln)) return(.lde_dev_pln(ctx, eta, hyp$ve))
	if (!isTRUE(ctx$dcor)) return(.lde_work(ctx$y, eta, ctx$link)$dev)
	pc <- ctx$pc
	yy <- ctx$y
	yy[c(pc$one, pc$two)] <- NA_real_
	dev <- .lde_work(yy, eta, ctx$link)$dev
	if (pc$n_pairs) dev <- dev - 2 * sum(
		.ae_pair_logphi2(pc$s1 * eta[pc$one], pc$s2 * eta[pc$two], pc$s1 * pc$s2 * hyp$rho_d)
	)
	dev
}

# scoring step of rho at a state: the Laplace-marginal score, the curvature at
# the mode (.lde_laplace()) and the next value; NULL without dcor, without a
# pair, or when rho is held
.lde_pair_rho_step <- function(ctx, hyp, lap, eta = NULL, vextra = NULL) {
	if (isTRUE(ctx$nrm)) return(if (is.null(eta)) NULL else .lde_nrm_rho_step(ctx, hyp, eta, lap, vextra))
	if (!isTRUE(ctx$dcor) || !is.null(ctx$rho_hold) || is.null(lap$rho_score)) return(NULL)
	q <- lap$rho_score
	if (!is.finite(q$rs)) return(NULL)
	# the step uses twice the plug-in curvature, close to the curvature of the
	# Laplace-marginal profile, so the iteration does not oscillate; the reported
	# curvature and sd are the plug-in's
	rho <- hyp$rho_d
	if (-q$rc > 0.05) {
		new <- .ae_pair_rho_next(list(rs = q$rs, rc = 2 * q$rc), rho)
	} else {
		# uninformative or non-concave pairs: a damped Newton step on the log
		# posterior with the prior's curvature, which returns rho to the prior mean
		g <- q$rs + 2 * .lame_rho_prior_power() * (-rho) / (1 - rho^2)
		h <- 2 * .lame_rho_prior_power() * (1 + rho^2) / (1 - rho^2)^2
		new <- max(min(rho + max(min(g / h, 0.2), -0.2), 0.98), -0.98)
	}
	list(rho = new, rs = q$rs, rc = q$rc)
}

# cell problem of a state under the noise of the product of the positions:
# independent cells (.lde_work_gh()) or pairs of cells under a dyadic correlation
.lde_work_noise <- function(ctx, eta, V, hyp, what) {
	if (isTRUE(ctx$nrm)) return(.lde_work_nrm(ctx, eta, hyp, V))
	if (isTRUE(ctx$pln)) return(.lde_work_gh_pln(ctx, eta, V, hyp$ve))
	if (isTRUE(ctx$dcor)) return(.lde_work_pair(ctx, ctx$y, eta, V, hyp$rho_d, what, ctx$link))
	.lde_work_gh(ctx$y, eta, V, what = what, link = ctx$link)
}

# cavity weights of the cells (.ae_mf_info_w()) or of the pairs of cells
# (.ae_pair_info_w(), over the total predictor variance and the additive-part
# covariance `lap$cvpair`). `w` diagonal weights, `cw` cross weights (NULL
# without a dyadic correlation)
.lde_cavity_w <- function(ctx, eta, V3, lap, hyp) {
	if (isTRUE(ctx$nrm)) return(.lde_cavity_w_nrm(ctx, eta, V3, lap, hyp))
	if (isTRUE(ctx$pln)) return(.lde_cavity_w_pln(ctx, eta, V3, hyp))
	if (!isTRUE(ctx$dcor) || is.null(lap$cvpair)) return(list(w = .ae_mf_info_w(ctx$y, eta, V3, ctx$link), cw = NULL))
	.ae_pair_info_w(list(pair = ctx$pc, y = ctx$y, link = ctx$link), eta, V3, lap$cvpair, hyp$rho_d)
}
