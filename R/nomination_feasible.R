# rebuild the rows of a latent slice that break the cbin/frn nomination constraints, so the
# Gibbs sweep starts inside its support; rows already inside it are returned untouched

# `Z` the latent slice, `Y` the observed nominations (0, positive ranks or counts, NA missing),
# `YL` the column of each rank (frn; NULL for cbin), `odmax` and `odobs` each row's cap and its
# nominations made
.lame_nom_feasible <- function(Z, Y, YL, odmax, odobs, bip) {
	nr <- nrow(Z)
	room <- rep_len(odobs < odmax, nr)
	room[is.na(room)] <- FALSE
	# hidden cells also use up a row's room, so rows with hidden cells skip the room rule
	bad <- nom_bad_rows_cpp(Z, Y, if (is.null(YL)) matrix(0L, 0L, 0L) else matrix(as.integer(YL), nrow(YL)), room, bip)
	if (!any(bad)) return(Z)
	nom <- !is.na(Y) & Y > 0
	non <- !is.na(Y) & Y == 0
	Zn <- Z
	Zn[!nom] <- Inf
	Zo <- Z
	Zo[!non] <- -Inf
	mx <- apply(Zo, 1L, max)
	for (i in which(bad)) {
		cols <- if (!is.null(YL)) {
			cc <- YL[i, ]
			cc[!is.na(cc)]
		} else which(nom[i, ])
		if (!length(cols)) next
		v <- sort(Z[i, cols])
		floor0 <- max(0.05, if (is.finite(mx[i]) && !room[i]) mx[i] + 0.05 else 0.05)
		v <- v - min(v) + floor0
		if (length(v) > 1L) for (k in 2:length(v)) v[k] <- max(v[k], v[k - 1L] + 1e-3)
		Z[i, cols] <- v
		o <- which(non[i, ])
		if (length(o)) Z[i, o] <- pmin(Z[i, o], min(v) - 0.05, if (room[i]) -0.05 else Inf)
	}
	Z
}
