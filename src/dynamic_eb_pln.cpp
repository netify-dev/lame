// Compiled kernels of the dynamic empirical-Bayes Poisson estimator (lame_dynamic_als() with family = "poisson"):
// the Poisson-lognormal cell of ae_pln.h with the fourth cumulant of the posterior of the noise, which gives the
// second derivative of the information in the predictor, and the same cell expected over a Gaussian predictor with
// the score, the information and the first two derivatives of the information. The cell and its arithmetic are those
// of pln_node() (ae_pln.h; R/ame_als_pln.R has the derivation): with pi(zeta) ~ exp(y zeta - exp(zeta) - (zeta - eta)^2 /
// (2 ve)) the posterior of zeta = eta + e, an exponential family in eta / ve, the derivatives of the log likelihood are
// its cumulants,
//
//     l' = E[e] / ve ,   w = -l'' = (ve - Var e) / ve^2 ,   h' = dw/deta = -kappa3 / ve^3 ,   h'' = -kappa4 / ve^4 ,
//
// and the pln_node() weights of the nodes carry everything else. Like the kernels of ae_pln_kernels.cpp, each function
// takes whole vectors, allocates its outputs before the loop, calls no R API inside it and writes each output element
// once, so the result does not depend on the number of threads. The R reference is .lde_pln_node_r() in
// R/lame_dynamic_als_eb_pln.R, which these agree with to rounding (the tests compare the two engines).

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

// pln_node() with the fourth cumulant
struct PlnNode4 {
	double l, s, w, hp, hpp, ee, me;
};

inline PlnNode4 pln_node4(double y, double eta, double ve, int K, const double* zk, const double* lwk,
                          double lgamma_y1) {
	eta = std::min(std::max(eta, -700.0), 700.0);
	ve = std::min(std::max(ve, 1e-8), 1e6);
	const double iv = 1.0 / ve;
	const double ly = y > 0.0 ? std::log(y) : -std::numeric_limits<double>::infinity();
	double zeta = std::max(eta, ly);
	const double u = y + eta / ve;
	for (int it = 0; it < 100; ++it) {
		const double ez = std::exp(zeta);
		const double step = (ez + zeta / ve - u) / (ez + iv);
		zeta -= step;
		if (std::fabs(step) <= 1e-13 * std::max(1.0, std::fabs(zeta))) break;
	}
	const double Xs = std::exp(zeta);
	const double s = 1.0 / std::sqrt(Xs + iv);
	const double estar = zeta - eta;
	double lw[kPlnMaxNodes], ex[kPlnMaxNodes];
	double mx = -std::numeric_limits<double>::infinity();
	for (int k = 0; k < K; ++k) {
		const double t = s * zk[k];
		const double e1 = std::expm1(t);
		const double c = Xs > 0.0 ? Xs * pln_remainder(t, e1) : 0.0;
		lw[k] = lwk[k] - c;
		ex[k] = e1;
		if (lw[k] > mx) mx = lw[k];
	}
	double tot = 0.0;
	for (int k = 0; k < K; ++k) {
		lw[k] = std::exp(lw[k] - mx);
		tot += lw[k];
	}
	double ez1 = 0.0;
	for (int k = 0; k < K; ++k) {
		lw[k] /= tot;
		ez1 += lw[k] * zk[k];
	}
	double vz = 0.0, k3 = 0.0, m4 = 0.0, eex = 0.0;
	for (int k = 0; k < K; ++k) {
		if (lw[k] == 0.0) continue;
		const double d = zk[k] - ez1;
		const double d2 = d * d;
		vz += lw[k] * d2;
		k3 += lw[k] * d2 * d;
		m4 += lw[k] * d2 * d2;
		eex += lw[k] * ex[k];
	}
	double cze = 0.0;
	for (int k = 0; k < K; ++k) {
		if (lw[k] == 0.0) continue;
		cze += lw[k] * (zk[k] - ez1) * (ex[k] - eex);
	}
	const double k4 = m4 - 3.0 * vz * vz;
	const double s2 = s * s;
	const double iv2 = iv * iv;
	PlnNode4 o;
	o.me = estar + s * ez1;
	o.s = o.me / ve;
	o.w = Xs * s * cze / ve;
	o.hp = -(s2 * s) * k3 * iv2 * iv;
	o.hpp = -(s2 * s2) * k4 * iv2 * iv2;
	o.ee = s2 * vz + o.me * o.me;
	const double gs = y * estar - Xs - estar * estar / (2.0 * ve);
	o.l = y * eta - lgamma_y1 - 0.5 * std::log(6.283185307179586476925286766559 * ve) + gs + std::log(s) +
	      0.918938533204672741780329736406 + mx + std::log(tot);
	return o;
}

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
// recycled): the log likelihood l, the score s, the observed information w, its first and second derivatives hp and
// hpp in eta, and the second moment ee and the mean me of the posterior of the noise. `lgy` is lgamma(y + 1); `zk`
// and `lwk` the nodes and the log weights of the inner Gauss-Hermite rule. Zero off the observed cells and where the
// predictor is not finite.
// [[Rcpp::export(rng = false)]]
List lde_pln_cells_cpp(const NumericVector& y, const NumericVector& eta, const NumericVector& ve,
                       const NumericVector& lgy, const NumericVector& zk, const NumericVector& lwk,
                       int threads = 1) {
	const std::size_t n = y.size();
	if ((std::size_t) eta.size() != n || (std::size_t) lgy.size() != n) stop("'y', 'eta' and 'lgy' must have the same length.");
	check_rule(zk, lwk);
	const PVec vv(ve, n, "ve");
	const int K = (int) zk.size();
	NumericVector l(n), s(n), w(n), hp(n), hpp(n), ee(n), me(n);
	const double* yp = y.begin();
	const double* ep = eta.begin();
	const double* gp = lgy.begin();
	const double* zp = zk.begin();
	const double* lp = lwk.begin();
	double *lo = l.begin(), *so = s.begin(), *wo = w.begin(), *ho = hp.begin(), *hho = hpp.begin(), *eo = ee.begin(),
	       *mo = me.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		if (!std::isfinite(yp[i]) || !std::isfinite(ep[i])) continue;
		const PlnNode4 nd = pln_node4(yp[i], ep[i], vv(i), K, zp, lp, gp[i]);
		lo[i] = nd.l;
		so[i] = nd.s;
		wo[i] = nd.w;
		ho[i] = nd.hp;
		hho[i] = nd.hpp;
		eo[i] = nd.ee;
		mo[i] = nd.me;
	}
	return List::create(Named("l") = l, Named("s") = s, Named("w") = w, Named("hp") = hp, Named("hpp") = hpp,
	                    Named("ee") = ee, Named("me") = me);
}

// The cells expected over a Gaussian predictor N(m, v) (the outer Gauss-Hermite rules zo, wo of the predictor, the
// variances vcut above which a cell takes the next rule, as in ae_pln_cells_cpp) and the inner rule of the noise: the
// expected score es, information ew, first and second derivatives of the information eh and ehh, the expected log
// likelihood el, the expected second moment of the noise ee, and the log of the expected probability of the count, lm
// = log E_eta p(y | eta), which is the marginal log likelihood of the cell with the predictor integrated out of the
// likelihood (not of its logarithm) and is the Poisson-lognormal log likelihood at the mean m with the variance v + ve,
// formed by its own inner rule (zk2, lwk2) because the variance it integrates over is larger than the noise's.
// Zero off the observed cells.
// [[Rcpp::export(rng = false)]]
List lde_pln_expect_cpp(const NumericVector& y, const NumericVector& lgy, const NumericVector& m,
                        const NumericVector& v, const NumericVector& ve, const List& zo, const List& wo,
                        const NumericVector& vcut, const NumericVector& zk, const NumericVector& lwk,
                        const NumericVector& zk2, const NumericVector& lwk2, int threads = 1) {
	const std::size_t n = y.size();
	if ((std::size_t) m.size() != n || (std::size_t) lgy.size() != n) stop("'y', 'm' and 'lgy' must have the same length.");
	check_rule(zk, lwk);
	check_rule(zk2, lwk2);
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
	const int K2 = (int) zk2.size();
	const double* zp2 = zk2.begin();
	const double* lp2 = lwk2.begin();
	NumericVector es(n), ew(n), eh(n), ehh(n), el(n), ee(n), lm(n);
	const double* yp = y.begin();
	const double* mp = m.begin();
	const double* gp = lgy.begin();
	const double* zp = zk.begin();
	const double* lp = lwk.begin();
	double *so = es.begin(), *wo_ = ew.begin(), *ho = eh.begin(), *hho = ehh.begin(), *lo = el.begin(),
	       *eo = ee.begin(), *mo = lm.begin();
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
		double s1 = 0.0, w1 = 0.0, h1 = 0.0, hh1 = 0.0, l1 = 0.0, e1 = 0.0;
		for (int j = 0; j < Ko; ++j) {
			const PlnNode4 nd = pln_node4(yp[i], mm + sdv * rz[j], vei, K, zp, lp, gp[i]);
			s1 += rw[j] * nd.s;
			w1 += rw[j] * nd.w;
			h1 += rw[j] * nd.hp;
			hh1 += rw[j] * nd.hpp;
			l1 += rw[j] * nd.l;
			e1 += rw[j] * nd.ee;
		}
		so[i] = s1;
		wo_[i] = w1;
		ho[i] = h1;
		hho[i] = hh1;
		lo[i] = l1;
		eo[i] = e1;
		// the marginal log likelihood: the cell at the mean m with the noise variance v + ve
		const PlnNode4 nm = pln_node4(yp[i], mm, vei + vi, K2, zp2, lp2, gp[i]);
		mo[i] = nm.l;
	}
	return List::create(Named("es") = es, Named("ew") = ew, Named("eh") = eh, Named("ehh") = ehh, Named("el") = el,
	                    Named("ee") = ee, Named("lm") = lm);
}
