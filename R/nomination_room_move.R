# collapsed translation of a period's intercept for periods that have rows with
# room (cbin and frn dynamic_beta fits, lame()). when nearly every row is at its
# cap the data barely inform the period's level and the augmentation sweep moves
# it as a slow random walk.
# the move shifts the whole latent slice of the period by s with its intercept:
# residuals are unchanged, rows at their cap translate rigidly, and the
# non-nominations of rows with room are integrated out and redrawn afterwards.
# the shift takes a Metropolis step with a symmetric random-walk proposal on its
# log-concave marginal density, so (shift, cells) is a block update. periods with
# a row whose hidden cells leave its room status open are left to the exact
# translations.

# what the data decide about the rows, once per fit: the complete rows with room
# the move integrates, the "unsure" rows whose period is skipped, the nomination
# cells (which bound the shift below) and each row's non-nomination columns
.lame_nom_room_context <- function(Y, odmax, bip) {
	n_r <- dim(Y)[1L]
	N <- dim(Y)[3L]
	cap <- rep_len(as.numeric(odmax), n_r)
	lapply(seq_len(N), function(t) {
		y <- Y[, , t]
		obs <- !is.na(y)
		if (!bip) diag(obs) <- FALSE
		hid <- .lame_nom_hidden_cells(y, bip)
		o1 <- obs & y > 0
		o0 <- obs & y == 0
		room <- rowSums(o1) < cap
		hidrow <- rowSums(hid) > 0L
		pinned <- room & rowSums(o0) > 0L & !hidrow
		list(
			pinned = pinned,
			unsure = room & hidrow & rowSums(o0) > 0L,
			nom_all = which(o1),
			Wcols = lapply(seq_len(n_r), function(i) if (pinned[i]) which(o0[i, ]) else integer(0))
		)
	})
}

# cells the move integrates out in rows `rows` of one period
.lame_nom_room_cells <- function(cx, Zt, EZt, rows, bip, rho) {
	nc <- lengths(cx$Wcols[rows])
	ii <- rep(rows, nc)
	jj <- unlist(cx$Wcols[rows], use.names = FALSE)
	k <- length(ii)
	out <- list(
		n = k,
		ii = ii,
		jj = jj,
		slack = if (k) -EZt[cbind(ii, jj)] else numeric(0),
		e = numeric(k),
		pw = logical(k),
		p1 = integer(0),
		p2 = integer(0),
		rho = rep(rho, k),
		sr = rep(sqrt(1 - rho^2), k)
	)
	if (k == 0L || bip || rho == 0) return(out)
	pi <- cbind(jj, ii)
	out$e <- Zt[pi] - EZt[pi]
	# a dyad of two integrated cells is treated jointly; the partner of any other is retained
	n_r <- nrow(Zt)
	key <- ii + (jj - 1L) * n_r
	pkey <- jj + (ii - 1L) * n_r
	kp <- match(pkey, key)
	out$pw <- !is.na(kp)
	first <- which(out$pw & seq_len(k) < kp)
	out$p1 <- first
	out$p2 <- kp[first]
	out$e[out$pw] <- 0
	out
}

# log density of the integrated cells at a shift s of their predictors
.lame_nom_room_logf <- function(cl, s) {
	if (cl$n == 0L) return(0)
	keep <- !cl$pw
	v <- 0
	if (any(keep)) v <- v + sum(
		stats::pnorm((cl$slack[keep] - s - cl$rho[keep] * cl$e[keep]) / cl$sr[keep], log.p = TRUE)
	)
	if (length(cl$p1)) {
		p <- pbivnorm::pbivnorm(cl$slack[cl$p1] - s, cl$slack[cl$p2] - s, cl$rho[cl$p1])
		v <- v + sum(log(pmax(p, 1e-300)))
	}
	v
}

# integrated cells redrawn at shift s from their conditional given the retained
# cells (a pair jointly, by rejection); NULL when the rejection fails
.lame_nom_room_redraw <- function(cl, s) {
	e_new <- numeric(cl$n)
	if (cl$n == 0L) return(e_new)
	sl <- cl$slack - s
	keep <- !cl$pw
	if (any(keep)) {
		w <- rtnorm_interval_logp(rep(0, sum(keep)), -Inf, (sl[keep] - cl$rho[keep] * cl$e[keep]) / cl$sr[keep])
		e_new[keep] <- cl$rho[keep] * cl$e[keep] + cl$sr[keep] * w
	}
	if (length(cl$p1)) {
		p1 <- cl$p1
		p2 <- cl$p2
		r <- cl$rho[p1]
		sr <- cl$sr[p1]
		pend <- seq_along(p1)
		e1 <- numeric(length(p1))
		e2 <- numeric(length(p1))
		for (pass in 1:5000) {
			if (!length(pend)) break
			x1 <- rtnorm_interval_logp(rep(0, length(pend)), -Inf, sl[p1[pend]])
			b2 <- (sl[p2[pend]] - r[pend] * x1) / sr[pend]
			ok <- stats::runif(length(pend)) < stats::pnorm(b2)
			if (any(ok)) {
				w <- rtnorm_interval_logp(rep(0, sum(ok)), -Inf, b2[ok])
				e1[pend[ok]] <- x1[ok]
				e2[pend[ok]] <- r[pend[ok]] * x1[ok] + sr[pend[ok]] * w
				pend <- pend[!ok]
			}
		}
		# a pair the rejection did not reach makes the whole move a rejection
		if (length(pend)) return(NULL)
		e_new[p1] <- e1
		e_new[p2] <- e2
	}
	e_new
}

# Metropolis steps on the shift: `target(s)` is the log density (-Inf outside its
# support), `lower` the lowest allowed shift; returns the final shift and the counts
.lame_nom_room_walk <- function(target, lower, tries, ladder) {
	s <- 0
	cur <- target(0)
	acc <- 0L
	if (!is.finite(cur)) return(list(s = 0, att = 0L, acc = 0L))
	for (k in seq_len(tries)) {
		sc <- ladder[sample.int(length(ladder), 1L)]
		prop <- s + sc * stats::rnorm(1L)
		if (!(prop > lower)) next
		new <- target(prop)
		if (is.finite(new) && log(stats::runif(1L)) < new - cur) {
			s <- prop
			cur <- new
			acc <- acc + 1L
		}
	}
	list(s = s, att = tries, acc = acc)
}

# prior terms of a shift of the intercept path by delta u: the log prior changes
# by -(Q delta^2 / 2 + lin delta)
.lame_nom_room_path_terms <- function(ic, u) {
	N <- length(u)
	du <- u[-1L] - ic$rho * u[-N]
	path <- ic$path
	c(
		Q = unname(u[1L]^2 / ic$v0 + sum(du^2) / ic$sd2),
		lin = unname(u[1L] * (path[1L] - ic$m0) / ic$v0 + sum(du * (path[-1L] - ic$rho * path[-N])) / ic$sd2)
	)
}

# one sweep of the move over the periods with rows with room, in random order;
# returns the updated Z, EZ and ic and the counts
.lame_nom_room_period <- function(Z, EZ, cx, bip, rho_t, ic, tries = 8L, ladder = c(0.02, 0.05, 0.1, 0.2, 0.4)) {
	n_r <- dim(Z)[1L]
	n_c <- dim(Z)[2L]
	N <- dim(Z)[3L]
	stats_out <- c(attempts = 0, accepts = 0, stuck = 0)
	res <- list(Z = Z, EZ = EZ, ic = ic, stats = stats_out)
	if (is.null(ic) || !identical(ic$kind, "path") || N < 2L || length(ic$path) != N) return(res)
	rho_t <- rep_len(rho_t, N)
	# slices as matrices, since the move changes one slice at a time
	Zl <- lapply(seq_len(N), function(t) Z[, , t])
	Ezl <- lapply(seq_len(N), function(t) EZ[, , t])
	moved <- FALSE
	for (t in sample.int(N)) {
		ct <- cx[[t]]
		pr <- which(ct$pinned)
		if (!length(pr) || any(ct$unsure)) next
		cl <- .lame_nom_room_cells(ct, Zl[[t]], Ezl[[t]], pr, bip, rho_t[t])
		if (cl$n == 0L) next
		# every nomination of the period translates with its row and must stay positive
		lower <- if (length(ct$nom_all)) -min(Zl[[t]][ct$nom_all]) else -Inf
		ut <- numeric(N)
		ut[t] <- 1
		pt <- .lame_nom_room_path_terms(ic, ut)
		target <- function(s) .lame_nom_room_logf(cl, s) - 0.5 * pt[["Q"]] * s^2 - pt[["lin"]] * s
		wk <- .lame_nom_room_walk(target, lower, tries, ladder)
		res$stats[c("attempts", "accepts")] <- res$stats[c("attempts", "accepts")] + c(wk$att, wk$acc)
		s <- wk$s
		if (s == 0) next
		e_new <- .lame_nom_room_redraw(cl, s)
		if (is.null(e_new)) {
			res$stats[["stuck"]] <- res$stats[["stuck"]] + 1
			next
		}
		ij <- cbind(cl$ii, cl$jj)
		new_w <- Ezl[[t]][ij] + s + e_new
		Zl[[t]] <- Zl[[t]] + s
		Zl[[t]][ij] <- new_w
		Ezl[[t]] <- Ezl[[t]] + s
		ic$path[t] <- ic$path[t] + s
		moved <- TRUE
	}
	if (moved) {
		res$Z <- array(unlist(Zl, use.names = FALSE), dim = c(n_r, n_c, N), dimnames = dimnames(Z))
		res$EZ <- array(unlist(Ezl, use.names = FALSE), dim = c(n_r, n_c, N), dimnames = dimnames(EZ))
	}
	res$ic <- ic
	res
}
