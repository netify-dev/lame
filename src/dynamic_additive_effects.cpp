// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// validate all periods before sampling so numerical overflow cannot masquerade as missingness.
static void ab_validate_data(const arma::mat& a_current, const arma::mat& b_current,
                             const arma::cube& Z_array, const arma::cube& EZ_array) {
	if (!a_current.is_finite() || !b_current.is_finite()) Rcpp::stop("Dynamic actor paths must be finite.");
	if (b_current.n_cols != a_current.n_cols || Z_array.n_slices != a_current.n_cols || EZ_array.n_slices != a_current.n_cols) Rcpp::stop("Dynamic actor period dimensions do not match.");
	if (Z_array.n_rows != a_current.n_rows || Z_array.n_cols != b_current.n_rows || EZ_array.n_rows != Z_array.n_rows || EZ_array.n_cols != Z_array.n_cols) Rcpp::stop("Dynamic actor observation dimensions do not match.");
	for (arma::uword i = 0; i < Z_array.n_elem; ++i) {
		if (std::isinf(Z_array[i]) || std::isinf(EZ_array[i]))
			Rcpp::stop("Dynamic actor observations cannot be infinite.");
		// intentional missing observations remain allowed.
		if (std::isfinite(Z_array[i]) && std::isfinite(EZ_array[i]) &&
		    !std::isfinite(Z_array[i] - EZ_array[i]))
			Rcpp::stop("Finite dynamic actor residual cannot be represented at this scale.");
	}
}

// joint draw of one margin's period-t effects when their conditionals are
// coupled (a symmetric network, or a directed one with a dyadic residual
// correlation): x ~ N(Q^-1 rhs, Q^-1), conditioned on sum(x) = 0 when asked.
// returns false when Q has no cholesky factor so the caller can fall back
static bool ab_joint_draw(const arma::mat& Q, const arma::vec& rhs, bool center,
                          arma::vec& x) {
	arma::mat L;
	if (!arma::chol(L, Q, "lower")) return false;
	arma::mat Lt = L.t();
	arma::vec m = arma::solve(arma::trimatu(Lt), arma::solve(arma::trimatl(L), rhs));
	arma::vec z(rhs.n_elem);
	for (arma::uword k = 0; k < z.n_elem; k++) z(k) = R::norm_rand();
	x = m + arma::solve(arma::trimatu(Lt), z);
	if (center) {
		arma::vec one(rhs.n_elem, arma::fill::ones);
		arma::vec u = arma::solve(arma::trimatu(Lt), arma::solve(arma::trimatl(L), one));
		double su = arma::accu(u);
		if (std::isfinite(su) && su > 0.0) x -= u * (arma::accu(x) / su);
	}
	return x.is_finite();
}

// conditional of independent draws x_i ~ N(m_i, v_i) given sum(x) = 0:
// subtract v_i sum(x) / sum(v). entries with v_i = 0 were not drawn and
// stay out of the constraint
static void ab_center_diag(arma::mat& M, int t, const arma::vec& v) {
	double sv = 0.0, sx = 0.0;
	for (arma::uword i = 0; i < v.n_elem; i++) {
		if (v(i) > 0.0 && std::isfinite(M(i, t))) { sv += v(i); sx += M(i, t); }
	}
	if (!(sv > 0.0) || !std::isfinite(sx)) return;
	for (arma::uword i = 0; i < v.n_elem; i++) {
		if (v(i) > 0.0 && std::isfinite(M(i, t))) M(i, t) -= v(i) * sx / sv;
	}
}

// ar(1) bridge prior of one actor's period-t state given its neighbours in
// time, with the stationary initial condition
static void ab_bridge_prior(const arma::mat& M, int i, int t, int T, double rho_ab,
                            double rho2, double var_innov, double& pp, double& pm) {
	if (T == 1) {
		pp = (1.0 - rho2) / var_innov; pm = 0.0;
	} else if (t == 0) {
		pp = 1.0 / var_innov; pm = rho_ab * M(i, 1);
	} else if (t == T - 1) {
		pp = 1.0 / var_innov; pm = rho_ab * M(i, T - 2);
	} else {
		pp = (1.0 + rho2) / var_innov;
		pm = rho_ab * (M(i, t - 1) + M(i, t + 1)) / (1.0 + rho2);
	}
	if (!std::isfinite(pm)) Rcpp::stop("Nonfinite additive AR conditional mean.");
}

// one margin's period-t step by single-site scans: each effect is drawn from
// its conditional given the freshest values of every other effect. `own` is
// the margin being drawn (rows when a_step, columns otherwise), `other` the
// reciprocal margin. v receives the conditional variances (0 where the draw
// fell back to the prior mean)
static void ab_step_seq(const arma::mat& resid_t, arma::mat& own, const arma::mat& other,
                        int t, int T, bool a_step, bool symmetric, bool use_dyad_rho,
                        double rho_ab, double rho2, double var_innov, double s2_safe,
                        double q_w, double p_w, arma::vec& v) {
	const int n = own.n_rows;
	const int m = a_step ? resid_t.n_cols : resid_t.n_rows;
	for (int i = 0; i < n; i++) {
		double pp, pm;
		ab_bridge_prior(own, i, t, T, rho_ab, rho2, var_innov, pp, pm);
		double dsum = 0.0, dsum_pair = 0.0, csum_pair = 0.0, dsum_iso = 0.0;
		int nobs = 0, npair = 0, niso = 0;
		double rec_i = 0.0;
		if (use_dyad_rho) {
			rec_i = other(i, t);
			if (!std::isfinite(rec_i)) rec_i = 0.0;
		}
		for (int j = 0; j < m; j++) {
			if (j == i && resid_t.n_rows == resid_t.n_cols) continue;
			double r = a_step ? resid_t(i, j) : resid_t(j, i);
			if (!std::isfinite(r)) continue;
			// the reciprocal effect on the same cell: the other margin's
			// effect, or this margin's own other actor on a symmetric network
			double rec = symmetric ? own(j, t) : other(j, t);
			if (!std::isfinite(rec)) rec = 0.0;
			dsum += r - rec;
			nobs++;
			if (use_dyad_rho) {
				double r_rev = a_step ? resid_t(j, i) : resid_t(i, j);
				if (std::isfinite(r_rev)) {
					// the reverse cell carries none of this effect and enters
					// through the cross weight only
					double own_j = own(j, t);
					if (!std::isfinite(own_j)) own_j = 0.0;
					dsum_pair += r - rec;
					csum_pair += r_rev - own_j - rec_i;
					npair++;
				} else {
					dsum_iso += r - rec;
					niso++;
				}
			}
		}
		double post_prec, post_mean;
		if (use_dyad_rho) {
			post_prec = pp + q_w * npair + niso / s2_safe;
			post_mean = (pp * pm + q_w * dsum_pair + p_w * csum_pair + dsum_iso / s2_safe) / post_prec;
		} else {
			post_prec = pp + nobs / s2_safe;
			post_mean = (pp * pm + dsum / s2_safe) / post_prec;
		}
		if (std::isfinite(post_prec) && post_prec > 0.0 && std::isfinite(post_mean)) {
			own(i, t) = R::rnorm(post_mean, 1.0 / sqrt(post_prec));
			v(i) = 1.0 / post_prec;
		} else {
			Rcpp::stop("Unable to resolve the additive Gaussian conditional.");
		}
	}
}

// the coupled case as one joint gaussian: precision Q and right-hand side of
// one margin's period-t conditional, with the reciprocal margin held fixed
static void ab_step_joint_system(const arma::mat& resid_t, const arma::mat& own,
                                 const arma::mat& other, int t, int T, bool a_step,
                                 bool symmetric, bool use_dyad_rho, double rho_ab,
                                 double rho2, double var_innov, double s2_safe,
                                 double q_w, double p_w, arma::mat& Q, arma::vec& rhs) {
	const int n = own.n_rows;
	const int m = a_step ? resid_t.n_cols : resid_t.n_rows;
	Q.zeros(n, n);
	rhs.zeros(n);
	for (int i = 0; i < n; i++) {
		double pp, pm;
		ab_bridge_prior(own, i, t, T, rho_ab, rho2, var_innov, pp, pm);
		Q(i, i) += pp;
		rhs(i) += pp * pm;
		double rec_i = 0.0;
		if (use_dyad_rho) {
			rec_i = other(i, t);
			if (!std::isfinite(rec_i)) rec_i = 0.0;
		}
		for (int j = 0; j < m; j++) {
			if (j == i && resid_t.n_rows == resid_t.n_cols) continue;
			double r = a_step ? resid_t(i, j) : resid_t(j, i);
			if (!std::isfinite(r)) continue;
			if (symmetric) {
				// resid(i, j) = a_i + a_j + e_ij, each dyad once per row
				Q(i, i) += 1.0 / s2_safe;
				Q(i, j) += 1.0 / s2_safe;
				rhs(i) += r / s2_safe;
				continue;
			}
			double rec = other(j, t);
			if (!std::isfinite(rec)) rec = 0.0;
			double r_rev = a_step ? resid_t(j, i) : resid_t(i, j);
			if (use_dyad_rho && std::isfinite(r_rev)) {
				Q(i, i) += q_w;
				Q(i, j) += p_w;
				rhs(i) += q_w * (r - rec) + p_w * (r_rev - rec_i);
			} else {
				Q(i, i) += 1.0 / s2_safe;
				rhs(i) += (r - rec) / s2_safe;
			}
		}
	}
}

// one margin's period-t step: a joint draw when the effects are coupled, the
// independent single-site draws otherwise; centred on request
static void ab_step(const arma::mat& resid_t, arma::mat& own, const arma::mat& other,
                    int t, int T, bool a_step, bool symmetric, bool use_dyad_rho,
                    bool center, double rho_ab, double rho2, double var_innov,
                    double s2_safe, double q_w, double p_w) {
	const int n = own.n_rows;
	const bool coupled = symmetric || use_dyad_rho;
	if (coupled) {
		arma::mat Q; arma::vec rhs, x;
		ab_step_joint_system(resid_t, own, other, t, T, a_step, symmetric, use_dyad_rho,
		                     rho_ab, rho2, var_innov, s2_safe, q_w, p_w, Q, rhs);
		if (ab_joint_draw(Q, rhs, center, x)) {
			own.col(t) = x;
			return;
		}
		Rcpp::stop("Unable to resolve the coupled additive Gaussian conditional.");
	}
	arma::vec v(n, arma::fill::zeros);
	ab_step_seq(resid_t, own, other, t, T, a_step, symmetric, use_dyad_rho,
	            rho_ab, rho2, var_innov, s2_safe, q_w, p_w, v);
	if (center) ab_center_diag(own, t, v);
}

//' Sample dynamic additive effects with AR(1) evolution
//'
//' Gibbs update of the per-period additive effects under the AR(1) state
//' model with stationary initial condition. For each actor and period the
//' full conditional combines the AR(1) bridge prior (stationary init at
//' t = 1) with the dyadic residual likelihood: resid_ij = a_i + b_j + e_ij.
//'
//' With \code{rho_dyad = 0} the dyads are independent, e_ij ~ N(0, s2), and
//' the conditional for a_i uses only row i of the residual. With
//' \code{rho_dyad != 0} the two directed residuals of a dyad are jointly
//' normal, (e_ij, e_ji) ~ N(0, s2 * ((1, rho), (rho, 1))), whose inverse
//' covariance is ((q, p), (p, q)) with q = 1 / (s2 (1 - rho^2)) and
//' p = -rho / (s2 (1 - rho^2)). The effect being drawn enters e_ij directly
//' and enters the pair's cross term through e_ji, so each fully observed
//' dyad contributes q to the precision and q (r_ij - b_j) + p (r_ji - a_j -
//' b_i) to the numerator of the a-step (symmetrically for the b-step). A
//' dyad with only one direction observed falls back to the marginal
//' isotropic weight 1 / s2.
//'
//' The reciprocal effect (b_j for the a-step, the freshly updated a_i for
//' the b-step) is subtracted from each residual, missing residuals are
//' skipped, the diagonal is skipped, and for symmetric networks each dyad
//' contributes exactly once (resid(i,j) = a_i + a_j + e_ij) so the pair term
//' does not apply and \code{rho_dyad} is ignored.
//'
//' @param a_current Current 2D array of row effects (n x T)
//' @param b_current Current 2D array of column effects (n x T)
//' @param Z_array 3D array of latent positions (n x n x T)
//' @param EZ_array 3D array of expected values without additive effects (n x n x T)
//' @param rho_ab AR(1) parameter for additive effects
//' @param sigma_ab Innovation standard deviation
//' @param s2 Dyadic residual variance
//' @param symmetric Whether the network is symmetric
//' @param rho_dyad Within-dyad correlation of the residual pair
//'   (e_ij, e_ji); 0 (the default) treats the dyads as independent. Ignored
//'   for symmetric networks and for non-square slices. Must lie strictly
//'   between -1 and 1 for a directed square panel.
//' @param center Whether each period's effects are drawn conditional on
//'   summing to zero across the actors of their margin, the identification
//'   used when the model has a static intercept. The draw is then the
//'   conditional of the block's joint normal given the constraint, not a
//'   projection of an unconstrained draw.
//' @param rvar Whether the sender margin is active.
//' @param cvar Whether the receiver margin is active.
//' @return List with updated a and b arrays
// [[Rcpp::export]]
List sample_dynamic_ab_cpp(arma::mat a_current, arma::mat b_current,
                           const arma::cube& Z_array, const arma::cube& EZ_array,
                           double rho_ab, double sigma_ab,
                           double s2, bool symmetric, double rho_dyad = 0.0,
                           bool center = false, bool rvar = true, bool cvar = true) {

	const int n = a_current.n_rows;
	const int T = a_current.n_cols;

	arma::mat a_new = a_current;
	arma::mat b_new = b_current;
	if (!rvar) a_new.zeros();
	if (!cvar) b_new.zeros();

	ab_validate_data(a_current, b_current, Z_array, EZ_array);

	const double var_innov = sigma_ab * sigma_ab;
	const double s2_safe = s2;
	const double rho2 = rho_ab * rho_ab;
	if (!std::isfinite(var_innov) || var_innov <= 0.0 ||
	    !std::isfinite(s2) || s2 <= 0.0 || !std::isfinite(rho_ab) || rho2 >= 1.0)
		Rcpp::stop("Invalid additive AR parameters or residual variance.");

	// within-dyad correlation weights: the pair (e_ij, e_ji) has covariance
	// s2 * [[1, rho], [rho, 1]] and inverse covariance [[q, p], [p, q]].
	// a symmetric network holds one observation per dyad, so the pair term
	// does not apply there and the isotropic weights are used instead
	const bool square = (Z_array.n_rows == Z_array.n_cols) &&
		((int)Z_array.n_rows == n) && ((int)b_new.n_rows == n);
	const double rho_d = (!symmetric && square) ? rho_dyad : 0.0;
	if (!std::isfinite(rho_d) || std::abs(rho_d) >= 1.0) {
		Rcpp::stop("Residual correlation must lie strictly between -1 and 1.");
	}
	const bool use_dyad_rho = (!symmetric) && square && (rho_d != 0.0);
	const double q_w = 1.0 / (s2_safe * (1.0 - rho_d * rho_d));
	const double p_w = -rho_d / (s2_safe * (1.0 - rho_d * rho_d));

	// t-outer sweep with in-place updates so each period's conditional sees
	// the freshest neighbouring states and the b-step the fresh a
	for(int t = 0; t < T; t++) {
		arma::mat resid_t = Z_array.slice(t) - EZ_array.slice(t);
		if (rvar) ab_step(resid_t, a_new, b_new, t, T, true, symmetric, use_dyad_rho, center,
		        rho_ab, rho2, var_innov, s2_safe, q_w, p_w);
		if (symmetric) {
			b_new.col(t) = a_new.col(t);
			continue;
		}
		if (cvar) ab_step(resid_t, b_new, a_new, t, T, false, false, use_dyad_rho, center,
		        rho_ab, rho2, var_innov, s2_safe, q_w, p_w);
	}

	return List::create(
		Named("a") = a_new,
		Named("b") = b_new
	);
}

//' Initialize dynamic additive effects with AR(1) structure
//'
//' @param n Number of actors
//' @param Tn Number of time points
//' @param rho_ab AR(1) parameter
//' @param sigma_ab Innovation standard deviation
//' @param mean_a Mean for row effects
//' @param mean_b Mean for column effects
//' @return List with initialized a and b matrices
// [[Rcpp::export]]
List init_dynamic_ab_cpp(int n, int Tn, double rho_ab, double sigma_ab,
                         double mean_a = 0.0, double mean_b = 0.0) {

	arma::mat a_mat(n, Tn);
	arma::mat b_mat(n, Tn);

	double var_stat = sigma_ab * sigma_ab / (1.0 - rho_ab * rho_ab);
	double sd_stat = sqrt(var_stat);

	for(int i = 0; i < n; i++) {
		a_mat(i, 0) = mean_a + R::rnorm(0, sd_stat);
		b_mat(i, 0) = mean_b + R::rnorm(0, sd_stat);

		for(int t = 1; t < Tn; t++) {
			a_mat(i, t) = rho_ab * a_mat(i, t-1) + R::rnorm(0, sigma_ab);
			b_mat(i, t) = rho_ab * b_mat(i, t-1) + R::rnorm(0, sigma_ab);
		}
	}

	return List::create(
		Named("a") = a_mat,
		Named("b") = b_mat
	);
}


// vector-aware internal variant; the scalar boundary remains unchanged.
// [[Rcpp::export]]
List sample_dynamic_ab_rho_cpp(arma::mat a_current, arma::mat b_current,
                           const arma::cube& Z_array, const arma::cube& EZ_array,
                           double rho_ab, double sigma_ab,
                           double s2, bool symmetric, Rcpp::NumericVector rho_dyad,
                           bool center = false, bool rvar = true, bool cvar = true) {

	const int n = a_current.n_rows;
	const int T = a_current.n_cols;

	arma::mat a_new = a_current;
	arma::mat b_new = b_current;
	if (!rvar) a_new.zeros();
	if (!cvar) b_new.zeros();

	ab_validate_data(a_current, b_current, Z_array, EZ_array);
	const double var_innov = sigma_ab * sigma_ab;
	const double s2_safe = s2;
	const double rho2 = rho_ab * rho_ab;
	if (!std::isfinite(var_innov) || var_innov <= 0.0 ||
	    !std::isfinite(s2) || s2 <= 0.0 || !std::isfinite(rho_ab) || rho2 >= 1.0)
		Rcpp::stop("Invalid additive AR parameters or residual variance.");

	// within-dyad correlation weights: the pair (e_ij, e_ji) has covariance
	// s2 * [[1, rho], [rho, 1]] and inverse covariance [[q, p], [p, q]].
	// a symmetric network holds one observation per dyad, so the pair term
	// does not apply there and the isotropic weights are used instead
	const bool square = (Z_array.n_rows == Z_array.n_cols) &&
		((int)Z_array.n_rows == n) && ((int)b_new.n_rows == n);
	if (rho_dyad.size() != T) Rcpp::stop("Residual correlation must have one value per period.");
	arma::vec q_period(T), p_period(T);
	for (int t = 0; t < T; ++t) {
		if (!std::isfinite(rho_dyad[t]) || std::abs(rho_dyad[t]) >= 1.0)
			Rcpp::stop("Residual correlations must be finite and strictly between -1 and 1.");
		const double r = (!symmetric && square) ? rho_dyad[t] : 0.0;
		const double denominator = s2_safe * ((1.0 - r) * (1.0 + r));
		q_period[t] = 1.0 / denominator;
		p_period[t] = -r * q_period[t];
		if (!std::isfinite(denominator) || denominator <= 0.0 || !std::isfinite(q_period[t]) || q_period[t] <= 0.0 || !std::isfinite(p_period[t]))
			Rcpp::stop("Unable to resolve finite dynamic actor residual precision.");
	}

	// t-outer sweep with in-place updates so each period's conditional sees
	// the freshest neighbouring states and the b-step the fresh a
	for(int t = 0; t < T; t++) {
		const double rho_d = (!symmetric && square) ? rho_dyad[t] : 0.0;
		const bool use_dyad_rho = (!symmetric) && square && (rho_d != 0.0);
		const double q_w = q_period[t];
		const double p_w = p_period[t];
		arma::mat resid_t = Z_array.slice(t) - EZ_array.slice(t);
		if (rvar) ab_step(resid_t, a_new, b_new, t, T, true, symmetric, use_dyad_rho, center,
		        rho_ab, rho2, var_innov, s2_safe, q_w, p_w);
		if (symmetric) {
			b_new.col(t) = a_new.col(t);
			continue;
		}
		if (cvar) ab_step(resid_t, b_new, a_new, t, T, false, false, use_dyad_rho, center,
		        rho_ab, rho2, var_innov, s2_safe, q_w, p_w);
	}

	if (!a_new.is_finite() || !b_new.is_finite()) Rcpp::stop("Nonfinite dynamic actor draw.");
	return List::create(
		Named("a") = a_new,
		Named("b") = b_new
	);
}
