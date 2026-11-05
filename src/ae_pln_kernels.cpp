// Compiled kernels of the empirical-Bayes Poisson estimator (ame_als(), lame_als() with
// family = "poisson"): the Poisson-lognormal cell of ae_pln.h over vectors of cells, its expectation over
// a Gaussian predictor (the cell quadrature of the variational stage) and the cavity weights. The R
// references are the functions with the suffix _r in R/ame_als_pln.R, which these agree with to rounding
// (the tests compare the two engines). Like the kernels of ae_kernels.cpp, each takes whole vectors,
// allocates its outputs before the loop, calls no R API inside it and writes each output element once, so
// the result does not depend on the number of threads.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "ae_kernels.h"
#include "ae_pln.h"

using namespace Rcpp;
using namespace lame_ae;

namespace {

// a vector that is either of full length or a scalar to recycle
struct PVec {
	const double* p;
	std::ptrdiff_t stride;
	PVec(const NumericVector& x, std::size_t n, const char* name) {
		if (x.size() == 1 && n != 1) {
			stride = 0;
		} else if ((std::size_t) x.size() == n) {
			stride = 1;
		} else {
			stop("Incompatible lengths: '%s' has %d elements, expected %d.", name, (int) x.size(), (int) n);
		}
		p = x.begin();
	}
	double operator()(std::ptrdiff_t i) const { return p[stride * i]; }
};

void check_rule(const NumericVector& zk, const NumericVector& lwk) {
	if (zk.size() < 1 || zk.size() > kPlnMaxNodes || zk.size() != lwk.size()) {
		stop("The inner quadrature rule must have between 1 and %d nodes.", kPlnMaxNodes);
	}
}

}  // namespace

// The cells of the Poisson-lognormal likelihood at the linear predictors `eta` (variance of the noise `ve`,
// recycled): the log likelihood l, the score s, the observed information w, its derivative hp in eta, and the
// second moment ee and the mean me of the posterior of the noise (.ae_pln_node_r in R/ame_als_pln.R). `lgy` is
// lgamma(y + 1); `zk` and `lwk` the nodes and the log weights of the inner Gauss-Hermite rule. Zero off the
// observed cells and where the predictor is not finite.
// [[Rcpp::export(rng = false)]]
List ae_pln_node_cpp(const NumericVector& y, const NumericVector& eta, const NumericVector& ve,
                     const NumericVector& lgy, const NumericVector& zk, const NumericVector& lwk,
                     int threads = 1) {
	const std::size_t n = y.size();
	if ((std::size_t) eta.size() != n || (std::size_t) lgy.size() != n) stop("'y', 'eta' and 'lgy' must have the same length.");
	check_rule(zk, lwk);
	const PVec vv(ve, n, "ve");
	const int K = (int) zk.size();
	NumericVector l(n), s(n), w(n), hp(n), ee(n), me(n);
	const double* yp = y.begin();
	const double* ep = eta.begin();
	const double* gp = lgy.begin();
	const double* zp = zk.begin();
	const double* lp = lwk.begin();
	double *lo = l.begin(), *so = s.begin(), *wo = w.begin(), *ho = hp.begin(), *eo = ee.begin(), *mo = me.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		if (!std::isfinite(yp[i]) || !std::isfinite(ep[i])) continue;
		const PlnNode nd = pln_node(yp[i], ep[i], vv(i), K, zp, lp, gp[i]);
		lo[i] = nd.l;
		so[i] = nd.s;
		wo[i] = nd.w;
		ho[i] = nd.hp;
		eo[i] = nd.ee;
		mo[i] = nd.me;
	}
	return List::create(Named("l") = l, Named("s") = s, Named("w") = w, Named("hp") = hp, Named("ee") = ee,
	                    Named("me") = me);
}

// The expected working problem of the cells whose linear predictor is N(m, v), the outer Gauss-Hermite rules
// (zo, wo: lists of the rules, vcut: the variances above which a cell takes the next rule) of the predictor and
// the inner rule (zk, lwk) of the noise (.ae_pln_cells_r): the working response z = m + E s / E w, the weight
// E w, the expected log likelihood ell and the expected second moment of the noise ee of every observed cell, zero
// off the observed cells. A cell takes the first rule when its predictor variance is at most vcut[0], the second
// when it is at most vcut[1], and so on: few nodes where the predictor is nearly known. The increment of the
// response is capped at 30.
// [[Rcpp::export(rng = false)]]
List ae_pln_cells_cpp(const NumericVector& y, const NumericVector& lgy, const NumericVector& m,
                      const NumericVector& v, const NumericVector& ve, const List& zo,
                      const List& wo, const NumericVector& vcut, const NumericVector& zk, const NumericVector& lwk,
                      int threads = 1) {
	const std::size_t n = y.size();
	if ((std::size_t) m.size() != n || (std::size_t) lgy.size() != n) stop("'y', 'm' and 'lgy' must have the same length.");
	check_rule(zk, lwk);
	const int nrule = zo.size();
	if (nrule < 1 || wo.size() != nrule || vcut.size() != nrule - 1) stop("The outer quadrature rules are not valid.");
	std::vector<const double*> gz(nrule), gw(nrule);
	std::vector<int> gk(nrule);
	std::vector<NumericVector> zs(nrule), ws(nrule);
	for (int r = 0; r < nrule; ++r) {
		zs[r] = as<NumericVector>(zo[r]);
		ws[r] = as<NumericVector>(wo[r]);
		if (zs[r].size() < 1 || zs[r].size() != ws[r].size()) stop("The outer quadrature rules are not valid.");
		gz[r] = zs[r].begin();
		gw[r] = ws[r].begin();
		gk[r] = (int) zs[r].size();
	}
	const double* vc = vcut.begin();
	const PVec vin(v, n, "v"), vv(ve, n, "ve");
	const int K = (int) zk.size();
	NumericVector z(n), w(n), ell(n), ee(n);
	const double* yp = y.begin();
	const double* mp = m.begin();
	const double* gp = lgy.begin();
	const double* zp = zk.begin();
	const double* lp = lwk.begin();
	double *zo_ = z.begin(), *wo_ = w.begin(), *eo_ = ell.begin(), *ee_ = ee.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		if (!std::isfinite(yp[i]) || !std::isfinite(mp[i])) continue;
		const double mm = mp[i];
		const double vi = std::max(vin(i), 0.0);
		const double sdv = std::sqrt(vi);
		const double vei = vv(i);
		int rule = 0;
		while (rule < nrule - 1 && vi > vc[rule]) ++rule;
		const int Ko = gk[rule];
		const double* rz = gz[rule];
		const double* rw = gw[rule];
		double es = 0.0, ew = 0.0, el = 0.0, eev = 0.0;
		for (int j = 0; j < Ko; ++j) {
			const PlnNode nd = pln_node(yp[i], mm + sdv * rz[j], vei, K, zp, lp, gp[i]);
			es += rw[j] * nd.s;
			ew += rw[j] * nd.w;
			el += rw[j] * nd.l;
			eev += rw[j] * nd.ee;
		}
		double delta = es / std::max(ew, 1e-12);
		if (std::fabs(delta) > 30.0) {
			delta = delta > 0.0 ? 30.0 : -30.0;
			ew = std::fabs(es) / 30.0;
		}
		zo_[i] = mm + delta;
		wo_[i] = ew;
		eo_[i] = el;
		ee_[i] = eev;
	}
	return List::create(Named("z") = z, Named("w") = w, Named("ell") = ell, Named("ee") = ee);
}

// The information weight of the cells for the covariance of the factor rows: the observed information of the
// Poisson-lognormal cell with mean m and variance v + ve (.ae_pln_info_w_r); zero off the observed cells.
// [[Rcpp::export(rng = false)]]
NumericVector ae_pln_info_w_cpp(const NumericVector& y, const NumericVector& lgy, const NumericVector& m,
                                const NumericVector& v, const NumericVector& ve, const NumericVector& zk,
                                const NumericVector& lwk, int threads = 1) {
	const std::size_t n = y.size();
	if ((std::size_t) m.size() != n || (std::size_t) lgy.size() != n) stop("'y', 'm' and 'lgy' must have the same length.");
	check_rule(zk, lwk);
	const PVec vin(v, n, "v"), vv(ve, n, "ve");
	const int K = (int) zk.size();
	NumericVector out(n);
	const double* yp = y.begin();
	const double* mp = m.begin();
	const double* gp = lgy.begin();
	const double* zp = zk.begin();
	const double* lp = lwk.begin();
	double* op = out.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		if (!std::isfinite(yp[i]) || !std::isfinite(mp[i])) continue;
		const PlnNode nd = pln_node(yp[i], mp[i], vv(i) + std::max(vin(i), 0.0), K, zp, lp, gp[i]);
		double w = nd.w;
		if (!std::isfinite(w) || w < 0.0) w = 0.0;
		op[i] = w;
	}
	return out;
}

// The pairs of cells of a dyad: for each pair the log probability L of the two counts at the predictor means
// (m1, m2) with the noise covariance (s11, s22, s12), the score (L1, L2) and the information (w11, w22, w12) in the
// means (.ae_pln2_node_r in R/ame_als_pln_pair.R). `zk` and `lwk` are the nodes and log weights of the rule on each
// axis (K nodes, at most kPln2MaxNodes).
// [[Rcpp::export(rng = false)]]
List ae_pln2_node_cpp(const NumericVector& y1, const NumericVector& y2, const NumericVector& lg1, const NumericVector& lg2,
                      const NumericVector& m1, const NumericVector& m2, const NumericVector& s11,
                      const NumericVector& s22, const NumericVector& s12, const NumericVector& zk,
                      const NumericVector& lwk, int threads = 1) {
	const std::size_t n = y1.size();
	if ((std::size_t) y2.size() != n || (std::size_t) m1.size() != n || (std::size_t) m2.size() != n ||
	    (std::size_t) lg1.size() != n || (std::size_t) lg2.size() != n) stop("The vectors of the pairs must have the same length.");
	if (zk.size() < 1 || zk.size() > kPln2MaxNodes || zk.size() != lwk.size()) {
		stop("The rule of the pair must have between 1 and %d nodes.", kPln2MaxNodes);
	}
	const PVec a11(s11, n, "s11"), a22(s22, n, "s22"), a12(s12, n, "s12");
	const int K = (int) zk.size();
	NumericVector L(n), L1(n), L2(n), w11(n), w22(n), w12(n);
	const double *p1 = y1.begin(), *p2 = y2.begin(), *g1 = lg1.begin(), *g2 = lg2.begin(), *q1 = m1.begin(), *q2 = m2.begin();
	const double *zp = zk.begin(), *lp = lwk.begin();
	double *oL = L.begin(), *oL1 = L1.begin(), *oL2 = L2.begin(), *o11 = w11.begin(), *o22 = w22.begin(), *o12 = w12.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		if (!std::isfinite(p1[i]) || !std::isfinite(p2[i]) || !std::isfinite(q1[i]) || !std::isfinite(q2[i])) continue;
		const Pln2Node nd = pln2_node(p1[i], p2[i], q1[i], q2[i], a11(i), a22(i), a12(i), K, zp, lp, g1[i], g2[i]);
		oL[i] = nd.L; oL1[i] = nd.L1; oL2[i] = nd.L2; o11[i] = nd.w11; o22[i] = nd.w22; o12[i] = nd.w12;
	}
	return List::create(Named("L") = L, Named("L1") = L1, Named("L2") = L2, Named("w11") = w11, Named("w22") = w22,
	                    Named("w12") = w12);
}
