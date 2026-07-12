// fused c++ kernel for the binary-family probit latent update.
//
// drop-in replacement for the r kernel rZ_bin_fc() in R/rZ_bin_fc.R:
// it samples from exactly the same full conditional, with the same
// sweep structure, but in a single cell-wise pass per triangle with
// no full-matrix temporaries (no t(Z)/t(EZ) copies, no full-matrix
// pnorm/qnorm). all randomness comes from r's rng (unif_rand /
// norm_rand / R::rnorm), so this is a valid drop-in gibbs kernel.
//
// the update has three stages, mirroring the r code:
//   1. truncated-normal gibbs by triangle: for y in {-1 (missing),
//      0, 1}, update the upper triangle then the lower triangle;
//      each cell z_ij is drawn from
//        N(ez_ij + rho * (z_ji - ez_ji), 1 - rho^2)
//      truncated to (-Inf, 0] for y = 0, [0, Inf) for y = 1, and
//      unconstrained for missing cells, conditioning on the current
//      value of the reciprocal cell z_ji.
//   2. joint proposal for consistent dyads: zp = ez + c*e + d*t(e)
//      with e iid standard normal, accepted only where both cells of
//      the dyad are sign-consistent with y (missing cells always
//      consistent).
//   3. diagonal refresh: z_ii ~ N(ez_ii, 1 + rho).
//
// numerical strategy: one-sided truncated draws use the log.p-scale
// inverse cdf, z = ez + sz * qnorm(log(u) + log-tail-mass, log.p),
// with pnorm/qnorm evaluated on the log scale in the tail that
// contains the truncation interval. this is exact for arbitrary |ez|
// (log-scale pnorm/qnorm stay accurate out to |ez| of several
// hundred), so no probability clamp is needed and no draw can land
// outside the truncation region: the earlier clamp + keep-current
// fallback froze constraint-contradicting cells (|ez| beyond ~8.2)
// at their initial, possibly sign-violating, values for the entire
// chain. a non-finite backstop keeps the current value in the
// (theoretically unreachable) case unif_rand() returns an endpoint.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

static inline double sgn_dbl(double x) {
	return (double)((x > 0.0) - (x < 0.0));
}

// one draw from the truncated-normal full conditional of a single cell
static inline double draw_tn_cell(double z_old, double ez, double sz,
                                  double lb, double ub) {
	double z_new;
	if (lb == R_NegInf && ub == R_PosInf) {
		// unconstrained (missing cell)
		z_new = ez + sz * norm_rand();
	} else if (ub == R_PosInf) {
		// truncation [lb, inf): upper-tail sampling on the log scale
		double logS = R::pnorm((lb - ez) / sz, 0.0, 1.0, 0, 1);
		z_new = ez + sz * R::qnorm(std::log(unif_rand()) + logS,
		                           0.0, 1.0, 0, 1);
	} else if (lb == R_NegInf) {
		// truncation (-inf, ub]: lower-tail sampling on the log scale
		double logF = R::pnorm((ub - ez) / sz, 0.0, 1.0, 1, 1);
		z_new = ez + sz * R::qnorm(std::log(unif_rand()) + logF,
		                           0.0, 1.0, 1, 1);
	} else {
		// two-sided (lb, ub): not used by the binary kernel, kept exact
		// for safety. work in the tail containing the interval.
		double u = unif_rand();
		if (ub - ez <= 0.0) {
			double la = R::pnorm((lb - ez) / sz, 0.0, 1.0, 1, 1);
			double lc = R::pnorm((ub - ez) / sz, 0.0, 1.0, 1, 1);
			double r = std::exp(la - lc);
			z_new = ez + sz * R::qnorm(lc + std::log(u + (1.0 - u) * r),
			                           0.0, 1.0, 1, 1);
		} else {
			double lqa = R::pnorm((lb - ez) / sz, 0.0, 1.0, 0, 1);
			double lqb = R::pnorm((ub - ez) / sz, 0.0, 1.0, 0, 1);
			double r = std::exp(lqb - lqa);
			z_new = ez + sz * R::qnorm(lqa + std::log((1.0 - u) + u * r),
			                           0.0, 1.0, 0, 1);
		}
	}
	// backstop only: with the log-scale inverse cdf a non-finite draw
	// requires unif_rand() to return an exact endpoint, which r's rng
	// never does
	if (!std::isfinite(z_new)) {
		return z_old;
	}
	return z_new;
}

// [[Rcpp::export]]
arma::mat rZ_bin_fused_cpp(arma::mat Z, const arma::mat& EZ, double rho,
                           const arma::mat& Y) {

	// this kernel assumes square (unipartite) networks
	if (Z.n_rows != Z.n_cols) {
		Rcpp::stop("rZ_bin_fused_cpp requires square matrices (unipartite). "
		           "Got %d x %d.", Z.n_rows, Z.n_cols);
	}
	if (EZ.n_rows != Z.n_rows || EZ.n_cols != Z.n_cols ||
	    Y.n_rows != Z.n_rows || Y.n_cols != Z.n_cols) {
		Rcpp::stop("rZ_bin_fused_cpp: Z, EZ and Y must have identical dimensions.");
	}

	const int n = (int)Z.n_rows;
	const double sz = std::sqrt(1.0 - rho * rho);

	// stage 1: truncated-normal gibbs by triangle, in the same sweep
	// order as the r kernel: y in {-1, 0, 1}, upper triangle then
	// lower triangle. within a block each cell conditions only on its
	// reciprocal cell (the other triangle), so sequential cell-wise
	// updates reproduce the r kernel's batch semantics exactly.
	const double lbs[3] = { R_NegInf, R_NegInf, 0.0 };
	const double ubs[3] = { R_PosInf, 0.0, R_PosInf };

	for (int yi = 0; yi < 3; ++yi) {
		const double yv = (double)(yi - 1); // -1 (missing), 0, 1
		const double lb = lbs[yi];
		const double ub = ubs[yi];
		for (int tri = 0; tri < 2; ++tri) {
			for (int j = 0; j < n; ++j) {
				const int i_lo = (tri == 0) ? 0 : j + 1;
				const int i_hi = (tri == 0) ? j : n;
				for (int i = i_lo; i < i_hi; ++i) {
					const double yc = Y(i, j);
					const double yw = std::isnan(yc) ? -1.0 : yc;
					if (yw != yv) {
						continue;
					}
					const double ez = EZ(i, j) + rho * (Z(j, i) - EZ(j, i));
					Z(i, j) = draw_tn_cell(Z(i, j), ez, sz, lb, ub);
				}
			}
		}
	}

	// stage 2: joint proposal for consistent dyads,
	// zp = ez + c*e + d*t(e), accepted only where both cells of the
	// dyad are sign-consistent with y (missing always consistent)
	const double cc = (std::sqrt(1.0 + rho) + std::sqrt(1.0 - rho)) / 2.0;
	const double dd = (std::sqrt(1.0 + rho) - std::sqrt(1.0 - rho)) / 2.0;

	for (int j = 1; j < n; ++j) {
		for (int i = 0; i < j; ++i) {
			const double e_ij = norm_rand();
			const double e_ji = norm_rand();
			const double zp_ij = EZ(i, j) + cc * e_ij + dd * e_ji;
			const double zp_ji = EZ(j, i) + cc * e_ji + dd * e_ij;
			const double yw_ij = std::isnan(Y(i, j)) ? -1.0 : Y(i, j);
			const double yw_ji = std::isnan(Y(j, i)) ? -1.0 : Y(j, i);
			const bool acc_ij = (yw_ij == -1.0) ||
				(sgn_dbl(zp_ij) == sgn_dbl(yw_ij - 0.5));
			const bool acc_ji = (yw_ji == -1.0) ||
				(sgn_dbl(zp_ji) == sgn_dbl(yw_ji - 0.5));
			if (acc_ij && acc_ji) {
				Z(i, j) = zp_ij;
				Z(j, i) = zp_ji;
			}
		}
	}

	// stage 3: diagonal refresh (the r kernel's diag(A) <- TRUE
	// proposal is immediately overwritten by this same draw, so the
	// diagonal law is exactly N(ez_ii, 1 + rho))
	const double sd_diag = std::sqrt(1.0 + rho);
	for (int i = 0; i < n; ++i) {
		Z(i, i) = R::rnorm(EZ(i, i), sd_diag);
	}

	return Z;
}
