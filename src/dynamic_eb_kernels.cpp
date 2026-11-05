// Kernels of the dynamic empirical-Bayes estimator for binary panels
// (R/lame_dynamic_als_eb.R, R/lame_dynamic_als_eb_pos.R).
//
// lde_gh_cells_cpp is the cell problem of a probit (or logit) cell whose predictor is
// N(eta, v): the expected score, Fisher weight and their derivatives by a
// Gauss-Hermite rule, and the expected and marginal log likelihoods. It is the
// arithmetic of .lde_work_gh_ref() (the R version, kept as the reference of the
// tests), cell by cell, and costs a fraction of it: one complementary error
// function per node where the reference takes two log-scale normal tails. The rule
// of a cell is the smallest of a ladder with at least per_sd * sqrt(v) nodes (the
// features of the integrand are of width one and the range is a few sd): the 9-node
// rule of the static estimator where the predictor is nearly known (v up to about
// 1), more where the positions leave it a large variance.
// Column-major storage throughout; cells with a missing outcome are zero.

#include <Rcpp.h>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

// log(sqrt(2 pi)), and 1 / sqrt(2) and 1 / sqrt(2 pi)
const double kLogSqrt2Pi = 0.91893853320467274178;
const double kInvSqrt2 = 0.70710678118654752440;

// beyond this |predictor| the tail probabilities are formed on the log scale by R's own
// routine (as the reference does): the complementary error function would underflow at 37.5,
// and the density divided by a tail probability near the underflow loses digits
const double kTailSwitch = 30.0;

// log Phi(x) for any x
inline double log_phi(double x) {
	if (!(std::fabs(x) <= kTailSwitch)) return R::pnorm(x, 0.0, 1.0, 1, 1);
	const double a = 0.5 * std::erfc(std::fabs(x) * kInvSqrt2);
	return x < 0.0 ? std::log(a) : std::log1p(-a);
}

// the probit terms of one node of predictor e: the inverse Mills ratios of the two outcomes
// (lam1 = phi / Phi(e), lam0 = phi / Phi(-e)) and the log probabilities of the two outcomes
inline void probit_node(double e, double& lam1, double& lam0, double& lp1, double& lp0, bool want_lp1, bool want_lp0) {
	const double ld = -0.5 * e * e - kLogSqrt2Pi;
	if (std::fabs(e) <= kTailSwitch) {
		// one tail by erfc, the other by its complement
		const double a = 0.5 * std::erfc(std::fabs(e) * kInvSqrt2);
		const double b = 1.0 - a;
		const double dens = std::exp(ld);
		if (e < 0.0) {
			lam1 = dens / a;
			lam0 = dens / b;
			if (want_lp1) lp1 = std::log(a);
			if (want_lp0) lp0 = std::log1p(-a);
		} else {
			lam1 = dens / b;
			lam0 = dens / a;
			if (want_lp1) lp1 = std::log1p(-a);
			if (want_lp0) lp0 = std::log(a);
		}
	} else {
		const double l1 = R::pnorm(e, 0.0, 1.0, 1, 1);
		const double l0 = R::pnorm(e, 0.0, 1.0, 0, 1);
		lam1 = std::exp(ld - l1);
		lam0 = std::exp(ld - l0);
		lp1 = l1;
		lp0 = l0;
	}
}

// the same for the logit link: the score of an outcome of 1 is 1 - p1 = p0 (the role of lam1) and minus the
// score of an outcome of 0 is p1 (lam0), so that the Fisher weight is lam1 * lam0 = p0 p1 as for the probit; the
// log probabilities are those of the logistic distribution
inline void logit_node(double e, double& lam1, double& lam0, double& lp1, double& lp0, bool want_lp1, bool want_lp0) {
	// p1 = plogis(e) from the exponential of the smaller side, so that neither tail loses digits
	const double ex = std::exp(-std::fabs(e));
	const double small = ex / (1.0 + ex);
	const double big = 1.0 / (1.0 + ex);
	const double p1 = e < 0.0 ? small : big;
	const double p0 = e < 0.0 ? big : small;
	lam1 = p0;
	lam0 = p1;
	if (want_lp1) lp1 = e < 0.0 ? e - std::log1p(ex) : -std::log1p(ex);
	if (want_lp0) lp0 = e < 0.0 ? -std::log1p(ex) : -e - std::log1p(ex);
}

}  // namespace

// The cell problem at predictors eta with Gaussian noise of variance vnoise (zero for the plain
// problem: one node at zero with weight one), by the rule (znode, wnode) of Gauss-Hermite nodes for
// the standard normal. y, eta, vnoise are equal-length vectors (a cell with a non-finite outcome is
// unobserved: zero everywhere). returns
//   z      the working response eta + E score / E weight, the increment capped at delta_cap
//          (where the cap binds the weight is |E score| / delta_cap: the product of weight and shift
//          stays the cell's expected score)
//   w      the expected Fisher weight (each node floored at wfloor)
//   score  the expected score
//   hprime, hpp  the expected first and second derivatives of the Fisher weight in the predictor
//          (zero where not finite)
//   dev    -2 sum E log p   (the expected log likelihood)
//   dev_marg  -2 sum log Phi(s eta / sqrt(1 + v))   (s = 2 y - 1: the noise integrated out of the
//          likelihood itself)
// `logit` selects the logistic link: the same problem with the scores and weights of the logit (the marginal
// likelihood of a cell, which has no closed form there, is then the quadrature of its probability).
// `what` says how much of it is wanted: 0 the working response and weight (what the position sweeps read),
// 1 also the expected score, the derivatives of the weight and dev_marg, 2 everything (dev, the expected
// log likelihood, needs a logarithm per node more). what is not wanted is not formed: zero, and NA for dev
// and dev_marg
// [[Rcpp::export(rng = false)]]
List lde_gh_cells_cpp(const NumericVector& y, const NumericVector& eta, const NumericVector& vnoise,
                      const List& znodes, const List& wnodes, double wfloor, double delta_cap,
                      int what, double per_sd, bool logit) {
	const R_xlen_t n = y.size();
	if (eta.size() != n || vnoise.size() != n) stop("y, eta and vnoise must have the same length.");
	const int nlev = (int) znodes.size();
	if (nlev < 1 || wnodes.size() != nlev) stop("The quadrature ladder needs matching nodes and weights.");
	std::vector<NumericVector> zl, wl;
	for (int l = 0; l < nlev; ++l) {
		zl.push_back(as<NumericVector>(znodes[l]));
		wl.push_back(as<NumericVector>(wnodes[l]));
		if (zl[l].size() < 1 || wl[l].size() != zl[l].size()) stop("The quadrature rule needs matching nodes and weights.");
	}
	NumericVector z(n), w(n), score(n), hprime(n), hpp(n);
	if (y.hasAttribute("dim")) {
		IntegerVector dm = y.attr("dim");
		z.attr("dim") = dm;
		w.attr("dim") = dm;
		score.attr("dim") = dm;
		hprime.attr("dim") = dm;
		hpp.attr("dim") = dm;
	}
	double sum_el = 0.0, sum_dm = 0.0;
	const bool want_derivs = what >= 1, want_el = what >= 2;
	for (R_xlen_t c = 0; c < n; ++c) {
		const double yc = y[c];
		if (!R_finite(yc)) continue;
		const bool one = (yc == 1.0);
		const double m = eta[c];
		double v = vnoise[c];
		if (v < 0.0) v = 0.0;
		const double sv = std::sqrt(v);
		// the rule: the smallest of the ladder with at least per_sd * sd nodes, the last one beyond it
		int lev = 0;
		const double need = per_sd * sv;
		while (lev + 1 < nlev && (double) zl[lev].size() < need) ++lev;
		const NumericVector& znode = zl[lev];
		const NumericVector& wnode = wl[lev];
		const int K = (int) znode.size();
		double es = 0.0, ew = 0.0, eh = 0.0, ehh = 0.0, el = 0.0, ep = 0.0;
		for (int k = 0; k < K; ++k) {
			const double e = sv * znode[k] + m;
			double lam1, lam0, lp1 = 0.0, lp0 = 0.0;
			if (logit) {
				logit_node(e, lam1, lam0, lp1, lp0, want_el && one, want_el && !one);
			} else {
				probit_node(e, lam1, lam0, lp1, lp0, want_el && one, want_el && !one);
			}
			const double wr = lam1 * lam0;
			const double wk = wnode[k];
			es += wk * (one ? lam1 : lam0);
			ew += wk * ((wr < wfloor) ? wfloor : wr);
			if (want_derivs) {
				if (logit) {
					// w = p0 p1, w' = w (p0 - p1), w'' = w' (p0 - p1) - 2 w^2
					const double mm = lam1 - lam0;
					const double hpr = wr * mm;
					eh += wk * hpr;
					ehh += wk * (hpr * mm - 2.0 * wr * wr);
					// the probability of the outcome, whose expectation is the marginal likelihood
					ep += wk * (one ? lam0 : lam1);
				} else {
					const double mm = lam0 - lam1 - 2.0 * e;
					const double hpr = wr * mm;
					const double hppr = hpr * mm + wr * (lam0 * lam0 + lam1 * lam1 + e * (lam1 - lam0) - 2.0);
					eh += wk * hpr;
					ehh += wk * hppr;
				}
			}
			if (want_el) el += wk * (one ? lp1 : lp0);
		}
		if (!one) es = -es;
		double delta = es / ew;
		if (std::fabs(delta) > delta_cap) {
			delta = delta < 0.0 ? -delta_cap : delta_cap;
			ew = std::fabs(es) / delta_cap;
		}
		z[c] = m + delta;
		w[c] = ew;
		if (want_derivs) {
			score[c] = es;
			hprime[c] = std::isfinite(eh) ? eh : 0.0;
			hpp[c] = std::isfinite(ehh) ? ehh : 0.0;
			sum_dm += logit ? std::log(ep) : log_phi((one ? 1.0 : -1.0) * m / std::sqrt(1.0 + v));
		}
		sum_el += el;
	}
	return List::create(Named("z") = z, Named("w") = w, Named("score") = score,
	                    Named("hprime") = hprime, Named("hpp") = hpp,
	                    Named("dev") = want_el ? -2.0 * sum_el : NA_REAL,
	                    Named("dev_marg") = want_derivs ? -2.0 * sum_dm : NA_REAL);
}
