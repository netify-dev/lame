#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// orthogonal dyad coordinates avoid subtracting nearly equal quadratics. the
// prior of rho is proportional to (1 - rho^2)^prior_power (the package prior,
// 2 Beta(2, 2) - 1, has power 1: .lame_rho_prior_power() in R/rho_prior.R)
static inline double rrho_log_target_cpp(
	double r, int m, double plus_norm, double minus_norm, double prior_power
) {
	const double log_det = std::log1p(r) + std::log1p(-r);
	const double plus = plus_norm / std::sqrt(1 + r);
	const double minus = minus_norm / std::sqrt(1 - r);
	return -.5 * m * log_det - plus * plus - minus * minus + prior_power * log_det;
}

// metropolis update for the dyadic correlation with independent replicate data
//'
//' @usage rrho_mh_rep_cpp(ET, rho, s2, p_wide = 0.1, prior_power = 1)
//' @param ET Array of square residual relational matrix series. The third
//' dimension of the array is for different replicates. Each slice of the array
//' according to the third dimension is a square residual relational matrix.
//' @param rho current value of rho
//' @param s2 current value of s2
//' @param p_wide probability of proposing from the wide independent
//' component instead of the local random walk. the two components form a
//' mixture of Metropolis kernels for the same target: the random walk moves
//' at the conditional scale, which is a unit-root walk when the likelihood
//' carries little information about rho, while the independent draw from
//' Uniform(-0.99, 0.99) ignores the current value and can cross the whole
//' range in one step. \code{p_wide = 0} gives the random walk alone.
//' @param prior_power exponent of \eqn{1 - \rho^2} in the prior density of
//' rho: 1 is the prior of ame() and lame(), \eqn{2\,\mathrm{Beta}(2, 2) - 1}.
//' The callers read it from \code{.lame_rho_prior_power()}.
//' @return a new value of rho
//' @author Shahryar Minhas
//' @keywords internal
//' @noRd

// [[Rcpp::export]]
double rrho_mh_rep_cpp(
	arma::cube ET, double rho, double s2, double p_wide = 0.1, double prior_power = 1.0
) {
	if (ET.n_rows != ET.n_cols || !std::isfinite(rho) || std::abs(rho) >= 1 ||
		!std::isfinite(s2) || s2 <= 0 || !std::isfinite(p_wide) || p_wide < 0 || p_wide > 1) {
		Rcpp::stop("rrho_mh_rep_cpp: require square residual slices, rho in (-1,1), positive finite s2, and p_wide in [0,1].");
	}
	// the prior (1 - rho^2)^prior_power is integrable only for prior_power > -1
	if (!std::isfinite(prior_power) || prior_power <= -1) {
		Rcpp::stop("rrho_mh_rep_cpp: prior_power must be finite and greater than -1.");
	}
	const arma::uvec index = arma::find(arma::trimatl(arma::ones(ET.n_rows, ET.n_cols)) == 0);
	arma::mat EM(index.n_elem * ET.n_slices, 2);
	bool zero_plus = true, zero_minus = true;
	for (arma::uword t = 0; t < ET.n_slices; ++t) {
		const arma::mat E = ET.slice(t), Et = E.t();
		if (index.n_elem > 0) {
			const arma::vec raw1 = E.elem(index), raw2 = Et.elem(index);
			zero_plus = zero_plus && arma::all(raw1 == -raw2);
			zero_minus = zero_minus && arma::all(raw1 == raw2);
			if (arma::any((raw1 / std::sqrt(s2) == 0) % (raw1 != 0)) ||
				arma::any((raw2 / std::sqrt(s2) == 0) % (raw2 != 0))) {
				Rcpp::stop("rrho_mh_rep_cpp: residuals cannot be represented at this variance.");
			}
			EM.submat(t * index.n_elem, 0, (t + 1) * index.n_elem - 1, 0) = E.elem(index) / std::sqrt(s2);
			EM.submat(t * index.n_elem, 1, (t + 1) * index.n_elem - 1, 1) = Et.elem(index) / std::sqrt(s2);
		}
	}
	if (!EM.is_finite()) {
		Rcpp::stop("rrho_mh_rep_cpp: ET contains non-finite residuals.");
	}
	const arma::vec plus = EM.col(0) / 2 + EM.col(1) / 2;
	const arma::vec minus = EM.col(0) / 2 - EM.col(1) / 2;
	auto stable_norm = [](const arma::vec& x) {
		if (x.n_elem == 0) return 0.;
		const double scale = arma::abs(x).max();
		return scale == 0 ? 0. : scale * std::sqrt(arma::accu(arma::square(x / scale)));
	};
	const double plus_norm = stable_norm(plus), minus_norm = stable_norm(minus);
	if (!std::isfinite(plus_norm) || !std::isfinite(minus_norm)) {
		Rcpp::stop("rrho_mh_rep_cpp: ET contains non-finite sufficient statistics.");
	}
	const int m = EM.n_rows;
	if ((plus_norm == 0 && !zero_plus) || (minus_norm == 0 && !zero_minus)) {
		Rcpp::stop("rrho_mh_rep_cpp: dyadic residual contrasts cannot be represented at this scale.");
	}
	// m dyads with a zero contrast leave (1 - rho^2)^(prior_power - m / 2) at an
	// endpoint, integrable when m / 2 - prior_power < 1 (.lame_rho_improper_dyads)
	const int improper_dyads = std::max(1, static_cast<int>(std::ceil(2 + 2 * prior_power - 1e-9)));
	if (m >= improper_dyads && (zero_plus || zero_minus)) {
		Rcpp::stop("rrho_mh_rep_cpp: correlation conditional is improper for exactly degenerate residual pairs.");
	}
	const double old_target = rrho_log_target_cpp(rho, m, plus_norm, minus_norm, prior_power);
	if (!std::isfinite(old_target)) {
		Rcpp::stop("rrho_mh_rep_cpp: current rho has a non-finite log density.");
	}
	// each mixture component separately preserves the full-support target
	if (p_wide > 0 && R::runif(0, 1) < p_wide) {
		const double rho_wide = 0.99;
		if (std::abs(rho) >= rho_wide) return rho;
		const double rho1 = R::runif(-rho_wide, rho_wide);
		const double ratio = rrho_log_target_cpp(rho1, m, plus_norm, minus_norm, prior_power) - old_target;
		return std::log(R::runif(0, 1)) < ratio ? rho1 : rho;
	}
	// a residual-only scale remains valid for constant or single dyad pairs
	double sr = 1. / std::sqrt(static_cast<double>(std::max(1, m)));
	if (m >= 2 && arma::stddev(EM.col(0)) > 0 && arma::stddev(EM.col(1)) > 0) {
		const arma::mat corr = arma::cor(EM);
		const double scale = 2 * (1 - corr(0, 1) * corr(0, 1)) / std::sqrt(static_cast<double>(m));
		if (std::isfinite(scale) && scale > 0) sr = scale;
	}
	const double lo = R::pnorm((-1 - rho) / sr, 0, 1, 1, 0);
	const double hi = R::pnorm((1 - rho) / sr, 0, 1, 1, 0);
	const double rho1 = rho + sr * R::qnorm(R::runif(lo, hi), 0, 1, 1, 0);
	// inverse-cdf rounding can reach an endpoint outside the open support
	if (!std::isfinite(rho1) || std::abs(rho1) >= 1) return rho;
	const double old_log_z = std::log(hi - lo);
	const double new_log_z = std::log(R::pnorm((1 - rho1) / sr, 0, 1, 1, 0) -
		R::pnorm((-1 - rho1) / sr, 0, 1, 1, 0));
	const double ratio = rrho_log_target_cpp(rho1, m, plus_norm, minus_norm, prior_power) - old_target + old_log_z - new_log_z;
	return std::log(R::runif(0, 1)) < ratio ? rho1 : rho;
}
