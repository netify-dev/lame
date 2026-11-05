// Compiled kernels of the binary probit empirical-bayes estimator (ame_als(), lame_als()) and of the
// dynamic estimator that shares them (lame_als()): the cell quadrature, the cavity weights, the
// bivariate normal probability and the pair quadrature. The R references are the functions with the
// suffix _r in R/ame_als_meanfield.R and R/ame_als_pair.R, which these agree with to rounding (typically
// 1e-12 and 1e-9 at worst on realistic predictors; the tests compare the two engines). The scalar
// functions are in ae_kernels.h.
//
// Every kernel takes whole vectors, allocates its outputs before the loop, calls no R API inside
// it and writes each output element once, so the result does not depend on the number of
// threads: `threads` > 1 splits the cells or pairs over OpenMP threads where the compiler has them.

#define USE_FC_LEN_T
#include <Rcpp.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>
#include "ae_kernels.h"

using namespace Rcpp;
using namespace lame_ae;

namespace {

struct GhRule {
	std::vector<double> z, w;
};

// the Gauss Hermite rule of K nodes for the standard normal (nodes z, weights summing to one), from the
// eigenvalues and the first components of the eigenvectors of the Jacobi matrix. It calls the
// LAPACK routine that R's eigen() calls on a symmetric matrix (dsyevr, all eigenvalues, lower
// triangle), so that the rule is the one .ae_mf_gh() gives to the last bit, and the order is the
// decreasing one of eigen()
const GhRule& gh_rule(int K) {
	static std::map<int, GhRule> cache;
	std::map<int, GhRule>::const_iterator hit = cache.find(K);
	if (hit != cache.end()) return hit->second;
	GhRule r;
	r.z.assign(K, 0.0);
	r.w.assign(K, 1.0);
	if (K > 1) {
		std::vector<double> J((std::size_t) K * K, 0.0), val(K), vec((std::size_t) K * K);
		for (int i = 1; i < K; ++i) {
			const double e = std::sqrt(i / 2.0);
			J[(i - 1) + (std::size_t) K * i] = e;
			J[i + (std::size_t) K * (i - 1)] = e;
		}
		char jobv = 'V', range = 'A', uplo = 'L';
		int n = K, lda = K, il = 0, iu = 0, m = 0, ldz = K, lwork = -1, liwork = -1, info = 0;
		double vl = 0.0, vu = 0.0, abstol = 0.0, tmp = 0.0;
		int itmp = 0;
		std::vector<int> isuppz(2 * (std::size_t) K);
		F77_CALL(dsyevr)(&jobv, &range, &uplo, &n, J.data(), &lda, &vl, &vu, &il, &iu, &abstol, &m, val.data(),
                   vec.data(), &ldz, isuppz.data(), &tmp, &lwork, &itmp, &liwork, &info FCONE FCONE FCONE);
		if (info != 0) stop("The Gauss-Hermite rule could not be computed.");
		lwork = (int) tmp;
		liwork = itmp;
		std::vector<double> work((std::size_t) lwork);
		std::vector<int> iwork((std::size_t) liwork);
		F77_CALL(dsyevr)(&jobv, &range, &uplo, &n, J.data(), &lda, &vl, &vu, &il, &iu, &abstol, &m, val.data(),
                   vec.data(), &ldz, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info FCONE FCONE FCONE);
		if (info != 0) stop("The Gauss-Hermite rule could not be computed.");
		// ascending in `val`; eigen() reverses them
		for (int j = 0; j < K; ++j) {
			const int s = K - 1 - j;
			const double first = vec[0 + (std::size_t) K * s];
			r.z[j] = std::sqrt(2.0) * val[s];
			r.w[j] = first * first;
		}
	}
	return cache.insert(std::make_pair(K, r)).first->second;
}

// a vector that is either of full length or a scalar to recycle
struct Vec {
	const double* p;
	std::ptrdiff_t stride;
	Vec(const NumericVector& x, std::size_t n, const char* name) {
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

struct IVec {
	const int* p;
	std::ptrdiff_t stride;
	IVec(const IntegerVector& x, std::size_t n, const char* name) {
		if (x.size() == 1 && n != 1) {
			stride = 0;
		} else if ((std::size_t) x.size() == n) {
			stride = 1;
		} else {
			stop("Incompatible lengths: '%s' has %d elements, expected %d.", name, (int) x.size(), (int) n);
		}
		p = x.begin();
	}
	int operator()(std::ptrdiff_t i) const { return p[stride * i]; }
};

}  // namespace

// --- small exports for the tests ------------------------------------------------------------------

// log Phi(x) (the scalar kernel the others use)
// [[Rcpp::export(rng = false)]]
NumericVector ae_log_pnorm_cpp(const NumericVector& x) {
	const std::size_t n = x.size();
	NumericVector out(n);
	for (std::size_t i = 0; i < n; ++i) out[i] = log_pnorm(x[i]);
	return out;
}

// the Gauss Hermite rule of K nodes for the standard normal: nodes z and weights summing to one
// [[Rcpp::export(rng = false)]]
List ae_gh_rule_cpp(int K) {
	if (K < 1 || K > 200) stop("The number of quadrature nodes must be between 1 and 200.");
	const GhRule& r = gh_rule(K);
	return List::create(Named("z") = NumericVector(r.z.begin(), r.z.end()),
	                    Named("w") = NumericVector(r.w.begin(), r.w.end()));
}

// --- the terms of the probit working response ---------------------------------------------------------

// The standard normal log density and the two log tail probabilities at a probit linear predictor
// (.ae_probit_terms_r in R/ame_als_irls.R): everything the working response, the deviance and the
// likelihood-equation certificate need of it
// [[Rcpp::export(rng = false)]]
List ae_probit_terms_cpp(const NumericVector& eta, int threads = 1) {
	const std::size_t n = eta.size();
	NumericVector ld(n), lp1(n), lp0(n);
	const double* e = eta.begin();
	double *a = ld.begin(), *b = lp1.begin(), *c = lp0.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		a[i] = dnorm_log(e[i]);
		b[i] = log_pnorm(e[i]);
		c[i] = log_pnorm(-e[i]);
	}
	return List::create(Named("ld") = ld, Named("lp1") = lp1, Named("lp0") = lp0);
}

// d h / d eta of the fisher weight h(eta) = lam1 lam0 of a probit cell, h' = h (lam0 - lam1 - 2 eta), lam1 = exp(ld - lp1),
// lam0 = exp(ld - lp0) (.ae_eb_hprime_r)
// [[Rcpp::export(rng = false)]]
NumericVector ae_probit_hprime_cpp(const NumericVector& eta, const NumericVector& ld, const NumericVector& lp1,
                                   const NumericVector& lp0, int threads = 1) {
	const std::size_t n = eta.size();
	if ((std::size_t) ld.size() != n || (std::size_t) lp1.size() != n || (std::size_t) lp0.size() != n) {
		stop("Incompatible lengths for the probit terms.");
	}
	NumericVector out(n);
	const double *e = eta.begin(), *a = ld.begin(), *b = lp1.begin(), *c = lp0.begin();
	double* o = out.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		const double lam1 = std::exp(a[i] - b[i]);
		const double lam0 = std::exp(a[i] - c[i]);
		o[i] = lam1 * lam0 * (lam0 - lam1 - 2.0 * e[i]);
	}
	return out;
}

// --- the cell quadrature ----------------------------------------------------------------------

// The expected working problem of the cells whose linear predictor is N(m, v), K nodes of the Gauss
// Hermite rule (.ae_mf_cells_r in R/ame_als_meanfield.R): the working response z = m + E psi / w,
// the weight w = -E psi' and the expected log likelihood ell of every observed cell (y finite);
// zero off the observed cells. The increment of the response is capped at 30.
// [[Rcpp::export(rng = false)]]
List ae_cells_cpp(const NumericVector& y, const NumericVector& m, const NumericVector& v,
                  bool probit, int K, int threads = 1) {
	const std::size_t n = y.size();
	if (K < 1 || K > 200) stop("The number of quadrature nodes must be between 1 and 200.");
	if ((std::size_t) m.size() != n) stop("'y' and 'm' must have the same length.");
	const Vec vv(v, n, "v");
	const GhRule& rule = gh_rule(K);
	NumericVector z(n), w(n), ell(n);
	const double* yp = y.begin();
	const double* mp = m.begin();
	double* zp = z.begin();
	double* wp = w.begin();
	double* ep = ell.begin();
	const double* gz = rule.z.data();
	const double* gw = rule.w.data();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		const double yi = yp[i];
		if (!std::isfinite(yi)) continue;
		const double sgn = 2.0 * yi - 1.0;
		const double mm = mp[i];
		const double sdv = std::sqrt(std::max(vv(i), 0.0));
		double es = 0.0, ew = 0.0, el = 0.0;
		for (int k = 0; k < K; ++k) {
			const double e = mm + sdv * gz[k];
			const double se = sgn * e;
			double lps, psi, wc;
			if (probit) {
				lps = log_pnorm(se);
				const double lam = std::exp(-0.5 * e * e - kLogSqrt2Pi - lps);
				psi = sgn * lam;
				wc = lam * (lam + se);
			} else {
				lps = log_plogis(se);
				const double p1 = std::exp(lps);
				psi = sgn * (1.0 - p1);
				wc = p1 * (1.0 - p1);
			}
			es += gw[k] * psi;
			ew += gw[k] * wc;
			el += gw[k] * lps;
		}
		double delta = es / std::max(ew, 1e-12);
		if (std::fabs(delta) > 30.0) {
			delta = delta > 0.0 ? 30.0 : -30.0;
			ew = std::fabs(es) / 30.0;
		}
		zp[i] = mm + delta;
		wp[i] = ew;
		ep[i] = el;
	}
	return List::create(Named("z") = z, Named("w") = w, Named("ell") = ell);
}

// The information weight of the cells for the covariance of the factor rows: the curvature
// -d^2 log Z / dm^2 of the log of the likelihood integrated over the Gaussian predictor
// (.ae_mf_info_w_r); zero off the observed cells.
// [[Rcpp::export(rng = false)]]
NumericVector ae_info_w_cpp(const NumericVector& y, const NumericVector& m, const NumericVector& v,
                            bool probit, int K, int threads = 1) {
	const std::size_t n = y.size();
	if (K < 1 || K > 200) stop("The number of quadrature nodes must be between 1 and 200.");
	if ((std::size_t) m.size() != n) stop("'y' and 'm' must have the same length.");
	const Vec vin(v, n, "v");
	const GhRule& rule = gh_rule(K);
	NumericVector out(n);
	const double* yp = y.begin();
	const double* mp = m.begin();
	double* op = out.begin();
	const double* gz = rule.z.data();
	const double* gw = rule.w.data();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (std::ptrdiff_t i = 0; i < nn; ++i) {
		const double yi = yp[i];
		if (!std::isfinite(yi)) continue;
		const double sgn = 2.0 * yi - 1.0;
		const double mm = mp[i];
		const double vv = std::max(vin(i), 0.0);
		double w;
		if (probit) {
			const double sq = std::sqrt(1.0 + vv);
			const double t = sgn * mm / sq;
			const double lam = std::exp(dnorm_log(t) - log_pnorm(t));
			double o = lam * (lam + t);
			if (t < -30.0) {
				const double t2 = t * t;
				o = 1.0 - 1.0 / t2 + 6.0 / (t2 * t2);
			}
			w = o / (1.0 + vv);
		} else {
			const double sdv = std::sqrt(vv);
			double z0 = 0.0, z1 = 0.0, z2 = 0.0;
			for (int k = 0; k < K; ++k) {
				const double x = sgn * (mm + sdv * gz[k]);
				const double p1 = 1.0 / (1.0 + std::exp(-x));
				z0 += gw[k] * p1;
				z1 += gw[k] * sgn * p1 * (1.0 - p1);
				z2 += gw[k] * p1 * (1.0 - p1) * (1.0 - 2.0 * p1);
			}
			w = (z1 / z0) * (z1 / z0) - z2 / z0;
		}
		if (!std::isfinite(w) || w < 0.0) w = 0.0;
		op[i] = w;
	}
	return out;
}

// --- the bivariate normal probability -----------------------------------------------------------

// Phi2(a, b; r) = P(X <= a, Y <= b), the scheme of pbivnorm::pbivnorm() (an infinite limit is
// the exact limit). Not log, and not accurate relative to a probability below 1e-6.
// [[Rcpp::export(rng = false)]]
NumericVector ae_pbivnorm_cpp(const NumericVector& a, const NumericVector& b, const NumericVector& r,
                              int threads = 1) {
	const std::size_t n = a.size();
	const Vec bb(b, n, "b"), rr(r, n, "r");
	NumericVector out(n);
	double* op = out.begin();
	const double* ap = a.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
	const double inf = std::numeric_limits<double>::infinity();
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if (nt > 1)
#endif
	{
		PlanCache cache;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
		for (std::ptrdiff_t i = 0; i < nn; ++i) {
			const double x = ap[i], y = bb(i), rho = rr(i);
			if (x == -inf || y == -inf) op[i] = 0.0;
			else if (x == inf && y == inf) op[i] = 1.0;
			else if (x == inf) op[i] = pnorm(y);
			else if (y == inf) op[i] = pnorm(x);
			else {
				const PairPlan& pp = cache.get(rho);
				const bool k0 = pp.bv.kind == 0;
				op[i] = bvu_plan(-x, -y, pp.bv, k0 ? phid_ref(x) : 0.0, k0 ? phid_ref(y) : 0.0);
			}
		}
	}
	return out;
}

// log Phi2(a, b; r): the log of the probability above 1e-6, the log-domain integral below it
// (.ae_pair_logphi2_r in R/ame_als_pair.R)
// [[Rcpp::export(rng = false)]]
NumericVector ae_pair_logphi2_cpp(const NumericVector& a, const NumericVector& b, const NumericVector& r,
                                  int threads = 1) {
	const std::size_t n = a.size();
	const Vec bb(b, n, "b"), rr(r, n, "r");
	NumericVector out(n);
	double* op = out.begin();
	const double* ap = a.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if (nt > 1)
#endif
	{
		PlanCache cache;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 256)
#endif
		for (std::ptrdiff_t i = 0; i < nn; ++i) op[i] = log_phi2_plan(ap[i], bb(i), cache.get(rr(i)));
	}
	return out;
}

// log P, P = Phi2(a, b; r), with its gradient (La, Lb), hessian (Laa, Lbb, Lab), the score of r
// (Lr = phi2 / P) and its derivative (Lrr) (.ae_pair_logp_r)
// [[Rcpp::export(rng = false)]]
List ae_pair_logp_cpp(const NumericVector& a, const NumericVector& b, const NumericVector& r,
                      int threads = 1) {
	const std::size_t n = a.size();
	const Vec bb(b, n, "b"), rr(r, n, "r");
	NumericVector L(n), La(n), Lb(n), Laa(n), Lbb(n), Lab(n), Lr(n), Lrr(n);
	const double* ap = a.begin();
	double *pL = L.begin(), *pLa = La.begin(), *pLb = Lb.begin(), *pLaa = Laa.begin(), *pLbb = Lbb.begin(),
		*pLab = Lab.begin(), *pLr = Lr.begin(), *pLrr = Lrr.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if (nt > 1)
#endif
	{
		PlanCache cache;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 256)
#endif
		for (std::ptrdiff_t i = 0; i < nn; ++i) {
			const PairPlan& pp = cache.get(rr(i));
			PairTerms t;
			pair_logp_pre(ap[i], bb(i), pp, pp.bv.kind == 0 ? phid_ref(ap[i]) : 0.0, t);
			pL[i] = t.L; pLa[i] = t.La; pLb[i] = t.Lb; pLaa[i] = t.Laa; pLbb[i] = t.Lbb; pLab[i] = t.Lab;
			pLr[i] = t.Lr; pLrr[i] = t.Lrr;
		}
	}
	return List::create(Named("L") = L, Named("La") = La, Named("Lb") = Lb, Named("Laa") = Laa,
	                    Named("Lbb") = Lbb, Named("Lab") = Lab, Named("Lr") = Lr, Named("Lrr") = Lrr);
}

// --- the pair quadrature ------------------------------------------------------------------------------

namespace {

// the four outcomes of a pair in the order (+,+), (+,-), (-,+), (-,-)
const double kSg1[4] = {1.0, 1.0, -1.0, -1.0};
const double kSg2[4] = {1.0, -1.0, 1.0, -1.0};

// the terms of the pair problem at the predictors (a, b) with the correlation r and the observed
// signs (s1, s2) (.ae_pair_node_r): the log probability of the observed pair, its score in the two
// predictors, the Fisher information (I11, I22, I12) and, with `deriv`, its derivatives in a and b
struct FisherNode {
	double lo, psi1, psi2, I11, I22, I12, I11a, I11b, I22a, I22b, I12a, I12b;
	double s0, rc0, psi1r, psi2r, I11r, I22r, I12r;
};

// with `rho_deriv` (which needs `deriv`) also the derivatives in the correlation: the score of the
// log probability in r (s0) and its curvature (rc0), the derivatives of the two scores (psi1r,
// psi2r) and of the Fisher information (I11r, I22r, I12r)
inline void fisher_node(double a, double b, double r, double s1, double s2, bool deriv, bool rho_deriv,
                        PlanCache& cache, FisherNode& o) {
	const double kap = cache.get(r).kap;
	const double ta = (a - r * b) / kap;
	const double tb = (b - r * a) / kap;
	const double fa = dnorm(a);
	const double fb = dnorm(b);
	// the four outcome probabilities: one bivariate normal probability and the margins give the
	// four, and a probability below 1e-5 is recomputed directly, which the subtraction cannot give
	const double pa = pnorm(a);
	const double pb = pnorm(b);
	const double p11 = std::exp(log_phi2_plan(a, b, cache.get(r)));
	double P[4] = {p11, pa - p11, pb - p11, 1.0 - pa - pb + p11};
	for (int j = 1; j < 4; ++j) {
		if (!(P[j] >= 1e-5)) P[j] = std::exp(log_phi2_plan(kSg1[j] * a, kSg2[j] * b, cache.get(kSg1[j] * kSg2[j] * r)));
	}
	double Pc[4];
	for (int j = 0; j < 4; ++j) Pc[j] = std::max(P[j], 1e-150);
	const double ptb = pnorm(tb), ntb = pnorm(-tb), pta = pnorm(ta), nta = pnorm(-ta);
	const double A[4] = {ptb, ntb, ptb, ntb};
	const double B[4] = {pta, pta, nta, nta};
	double F11 = 0.0, F22 = 0.0, F12 = 0.0;
	for (int j = 0; j < 4; ++j) {
		F11 += A[j] * A[j] / Pc[j];
		F22 += B[j] * B[j] / Pc[j];
		F12 += kSg1[j] * kSg2[j] * A[j] * B[j] / Pc[j];
	}
	o.I11 = fa * fa * F11;
	o.I22 = fb * fb * F22;
	o.I12 = fa * fb * F12;
	const int k = 2 * (s1 < 0.0 ? 1 : 0) + (s2 < 0.0 ? 1 : 0);
	const double Po = Pc[k];
	o.psi1 = s1 * fa * A[k] / Po;
	o.psi2 = s2 * fb * B[k] / Po;
	o.lo = std::log(Po);
	if (!deriv) return;
	const double ft = dnorm(tb);
	const double fta = dnorm(ta);
	double Pa[4], Pb[4];
	double F11a = 0.0, F11b = 0.0, F22a = 0.0, F22b = 0.0, F12a = 0.0, F12b = 0.0;
	for (int j = 0; j < 4; ++j) {
		const double Aa = (-r / kap * ft) * kSg2[j];
		const double Ab = (ft / kap) * kSg2[j];
		const double Ba = (fta / kap) * kSg1[j];
		const double Bb = (-r / kap * fta) * kSg1[j];
		Pa[j] = (fa * kSg1[j]) * A[j];
		Pb[j] = (fb * kSg2[j]) * B[j];
		const double P2 = Pc[j] * Pc[j];
		const double ss = kSg1[j] * kSg2[j];
		F11a += 2.0 * A[j] * Aa / Pc[j] - A[j] * A[j] * Pa[j] / P2;
		F11b += 2.0 * A[j] * Ab / Pc[j] - A[j] * A[j] * Pb[j] / P2;
		F22a += 2.0 * B[j] * Ba / Pc[j] - B[j] * B[j] * Pa[j] / P2;
		F22b += 2.0 * B[j] * Bb / Pc[j] - B[j] * B[j] * Pb[j] / P2;
		F12a += ss * ((Aa * B[j] + A[j] * Ba) / Pc[j] - A[j] * B[j] * Pa[j] / P2);
		F12b += ss * ((Ab * B[j] + A[j] * Bb) / Pc[j] - A[j] * B[j] * Pb[j] / P2);
	}
	o.I11a = -2.0 * a * fa * fa * F11 + fa * fa * F11a;
	o.I11b = fa * fa * F11b;
	o.I22a = fb * fb * F22a;
	o.I22b = -2.0 * b * fb * fb * F22 + fb * fb * F22b;
	o.I12a = -a * fa * fb * F12 + fa * fb * F12a;
	o.I12b = -b * fa * fb * F12 + fa * fb * F12b;
	// a pair so far in the tail that a probability is below the floor has information below 1e-40,
	// and the derivatives of the terms of such an outcome are not numbers: they carry no weight
	double* d[6] = {&o.I11a, &o.I11b, &o.I22a, &o.I22b, &o.I12a, &o.I12b};
	for (int j = 0; j < 6; ++j) {
		if (!std::isfinite(*d[j])) *d[j] = 0.0;
	}
	if (!rho_deriv) return;
	// the probability of an outcome has P_r = s phi2(a, b; r) (Plackett, s = s1 s2 of the outcome and
	// phi2 the bivariate normal density phi(a) phi(tb) / kap) and the mixed derivatives
	// P_ar = s phi(a) phi(tb) (r b - a) / kap^3, P_br = s phi(b) phi(ta) (r a - b) / kap^3
	const double kap3 = kap * kap * kap;
	const double phi2 = fa * ft / kap;
	double Pr[4], Par[4], Pbr[4];
	double I11r = 0.0, I22r = 0.0, I12r = 0.0;
	for (int j = 0; j < 4; ++j) {
		const double ss = kSg1[j] * kSg2[j];
		Pr[j] = ss * phi2;
		Par[j] = ss * (fa * ft * (r * b - a) / kap3);
		Pbr[j] = ss * (fb * fta * (r * a - b) / kap3);
		const double P2 = Pc[j] * Pc[j];
		I11r += 2.0 * Pa[j] * Par[j] / Pc[j] - Pa[j] * Pa[j] * Pr[j] / P2;
		I22r += 2.0 * Pb[j] * Pbr[j] / Pc[j] - Pb[j] * Pb[j] * Pr[j] / P2;
		I12r += (Par[j] * Pb[j] + Pa[j] * Pbr[j]) / Pc[j] - Pa[j] * Pb[j] * Pr[j] / P2;
	}
	o.I11r = I11r;
	o.I22r = I22r;
	o.I12r = I12r;
	const double Po_r = Pr[k];
	o.s0 = Po_r / Po;
	o.psi1r = Par[k] / Po - Pa[k] * Po_r / (Po * Po);
	o.psi2r = Pbr[k] / Po - Pb[k] * Po_r / (Po * Po);
	const double Lr = phi2 / Po;
	const double nn = a * a - 2.0 * r * a * b + b * b;
	const double omr = 1.0 - r * r;
	const double dl2 = s1 * s2 * (r * omr + a * b * omr - r * nn) / (omr * omr);
	o.rc0 = Lr * dl2 - Lr * Lr;
	double* e[7] = {&o.I11r, &o.I22r, &o.I12r, &o.s0, &o.psi1r, &o.psi2r, &o.rc0};
	for (int j = 0; j < 7; ++j) {
		if (!std::isfinite(*e[j])) *e[j] = 0.0;
	}
}

}  // namespace

// The expected pair problem of a bivariate probit pair under the Gaussian predictors (eta_1, eta_2)
// with means (m1, m2), variances (v1, v2) and covariance cv, by the tensor Gauss Hermite rule of
// K nodes a side on the Cholesky factor (nodes e1 = m1 + sd1 z_k, e2 = m2 + (cv / sd1) z_k +
// sqrt(v2 - cv^2 / v1) z_l). One row per pair; s1, s2 the signs 2 y - 1 of the cells; rho the
// correlation; cv, rho, s1, s2 and K are recycled when of length one (K may differ by pair: the
// rule of each pair is its own). `mode`:
//   0  the observed information of the static estimator (.ae_pair_quad_r): ell, psi1, psi2,
//      w11, w22, w12, rs, rc  (the expected log probability, the expected score, minus the
//      expected hessian, and the score and curvature of rho)
//   1  the expected Fisher information over the four outcomes, for the dynamic estimator
//      (.ae_pair_quad_fisher_r): ell, psi1, psi2, I11, I22, I12; with `deriv` also its derivatives in
//      the two means, I11a, I11b, I22a, I22b, I12a, I12b; with `deriv` and `rho_deriv` also the
//      expectations of the derivatives in the correlation, s0, rc0, psi1r, psi2r, I11r, I22r, I12r
// [[Rcpp::export(rng = false)]]
List ae_pair_quad_cpp(const NumericVector& m1, const NumericVector& m2, const NumericVector& v1,
                      const NumericVector& v2, const NumericVector& cv, const NumericVector& s1,
                      const NumericVector& s2, const NumericVector& rho, const IntegerVector& K,
                      int mode = 0, bool deriv = false, bool rho_deriv = false, int threads = 1) {
	const std::size_t n = m1.size();
	const Vec m2v(m2, n, "m2"), v1v(v1, n, "v1"), v2v(v2, n, "v2"), cvv(cv, n, "cv"), s1v(s1, n, "s1"),
		s2v(s2, n, "s2"), rhov(rho, n, "rho");
	const IVec Kv(K, n, "K");
	if (mode < 0 || mode > 1) stop("'mode' must be 0 (observed information) or 1 (Fisher information).");
	rho_deriv = rho_deriv && deriv;
	// the rules, formed before the loop
	std::vector<const GhRule*> rules(201, (const GhRule*) 0);
	for (R_xlen_t i = 0; i < K.size(); ++i) {
		const int k = K[i];
		if (k < 1 || k > 200) stop("The number of quadrature nodes must be between 1 and 200.");
		if (!rules[k]) rules[k] = &gh_rule(k);
	}
	const double* m1p = m1.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
	const std::ptrdiff_t nn = (std::ptrdiff_t) n;
	// outputs: mode 0 has 8, mode 1 has 6, 12 or 19
	const int nout = mode == 0 ? 8 : (6 + (deriv ? 6 : 0) + (rho_deriv ? 7 : 0));
	std::vector<NumericVector> out(nout);
	std::vector<double*> op(nout);
	for (int j = 0; j < nout; ++j) {
		out[j] = NumericVector(n);
		op[j] = out[j].begin();
	}
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if (nt > 1)
#endif
	{
		PlanCache cache;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
		for (std::ptrdiff_t i = 0; i < nn; ++i) {
			const GhRule& g = *rules[Kv(i)];
			const int kk = (int) g.z.size();
			const double mm1 = m1p[i], mm2 = m2v(i);
			const double sd1 = std::sqrt(std::max(v1v(i), 0.0));
			const double l21 = sd1 > 0.0 ? cvv(i) / std::max(sd1, 1e-300) : 0.0;
			const double l22 = std::sqrt(std::max(v2v(i) - l21 * l21, 0.0));
			const double sg1 = s1v(i), sg2 = s2v(i);
			const double ss = sg1 * sg2;
			const double rh = rhov(i);
			if (mode == 0) {
				const PairPlan& pp = cache.get(ss * rh);
				const bool k0 = pp.bv.kind == 0;
				double ell = 0.0, psi1 = 0.0, psi2 = 0.0, w11 = 0.0, w22 = 0.0, w12 = 0.0, rs = 0.0, rc = 0.0;
				for (int k = 0; k < kk; ++k) {
					const double e1 = mm1 + sd1 * g.z[k];
					const double a = sg1 * e1;
					const double pa_ref = k0 ? phid_ref(a) : 0.0;
					for (int l = 0; l < kk; ++l) {
						const double wt = g.w[k] * g.w[l];
						const double e2 = mm2 + l21 * g.z[k] + l22 * g.z[l];
						PairTerms d;
						pair_logp_pre(a, sg2 * e2, pp, pa_ref, d);
						ell += wt * d.L;
						psi1 += wt * sg1 * d.La;
						psi2 += wt * sg2 * d.Lb;
						w11 -= wt * d.Laa;
						w22 -= wt * d.Lbb;
						w12 -= wt * ss * d.Lab;
						rs += wt * ss * d.Lr;
						rc += wt * d.Lrr;
					}
				}
				op[0][i] = ell; op[1][i] = psi1; op[2][i] = psi2; op[3][i] = w11; op[4][i] = w22;
				op[5][i] = w12; op[6][i] = rs; op[7][i] = rc;
			} else {
				double acc[19];
				for (int j = 0; j < 19; ++j) acc[j] = 0.0;
				for (int k = 0; k < kk; ++k) {
					const double e1 = mm1 + sd1 * g.z[k];
					for (int l = 0; l < kk; ++l) {
						const double wt = g.w[k] * g.w[l];
						const double e2 = mm2 + l21 * g.z[k] + l22 * g.z[l];
						FisherNode nd;
						fisher_node(e1, e2, rh, sg1, sg2, deriv, rho_deriv, cache, nd);
						acc[0] += wt * nd.lo;
						acc[1] += wt * nd.psi1;
						acc[2] += wt * nd.psi2;
						acc[3] += wt * nd.I11;
						acc[4] += wt * nd.I22;
						acc[5] += wt * nd.I12;
						if (deriv) {
							acc[6] += wt * nd.I11a;
							acc[7] += wt * nd.I11b;
							acc[8] += wt * nd.I22a;
							acc[9] += wt * nd.I22b;
							acc[10] += wt * nd.I12a;
							acc[11] += wt * nd.I12b;
						}
						if (rho_deriv) {
							acc[12] += wt * nd.s0;
							acc[13] += wt * nd.rc0;
							acc[14] += wt * nd.psi1r;
							acc[15] += wt * nd.psi2r;
							acc[16] += wt * nd.I11r;
							acc[17] += wt * nd.I22r;
							acc[18] += wt * nd.I12r;
						}
					}
				}
				for (int j = 0; j < nout; ++j) op[j][i] = acc[j];
			}
		}
	}
	if (mode == 0) {
		return List::create(Named("ell") = out[0], Named("psi1") = out[1], Named("psi2") = out[2],
		                    Named("w11") = out[3], Named("w22") = out[4], Named("w12") = out[5],
		                    Named("rs") = out[6], Named("rc") = out[7]);
	}
	List res = List::create(Named("ell") = out[0], Named("psi1") = out[1], Named("psi2") = out[2],
	                        Named("I11") = out[3], Named("I22") = out[4], Named("I12") = out[5]);
	if (deriv) {
		res["I11a"] = out[6]; res["I11b"] = out[7]; res["I22a"] = out[8];
		res["I22b"] = out[9]; res["I12a"] = out[10]; res["I12b"] = out[11];
	}
	if (rho_deriv) {
		res["s0"] = out[12]; res["rc0"] = out[13]; res["psi1r"] = out[14]; res["psi2r"] = out[15];
		res["I11r"] = out[16]; res["I22r"] = out[17]; res["I12r"] = out[18];
	}
	return res;
}

// --- the cross terms of the additive block ----------------------------------------------------------------

// The sums over time that the cross weights of the pairs add to the normal equations of the additive block
// (.ae_pair_cross_r in R/ame_als_pair.R): with x_1 = 1 and x_{k+1} the dyadic covariates (a non-finite entry counts
// as zero), cw_t the cross weights of period t (symmetric in the two actors of a dyad, zero where there is no pair),
//   Mk[[k]]   = sum_t x_k * cw_t                    (n x n)
//   xx[k, l]  = sum_t sum_ij x_k,ij cw_t,ij x_l,ji   (pa x pa)
// The covariates are an n x n x p x Tt array and the weights an n x n x Tt array.
// [[Rcpp::export(rng = false)]]
List ae_pair_cross_cpp(const NumericVector& cw, const NumericVector& X, int n, int Tt, int p) {
	const int pa = p + 1;
	const std::size_t nn = (std::size_t) n * (std::size_t) n;
	if ((std::size_t) cw.size() != nn * (std::size_t) Tt || (std::size_t) X.size() != nn * (std::size_t) p * (std::size_t) Tt) {
		stop("Incompatible dimensions for the cross sums.");
	}
	NumericVector Mk(nn * (std::size_t) pa);
	NumericMatrix xx(pa, pa);
	double* mk = Mk.begin();
	std::vector<double> xs(nn * (std::size_t) pa), xk(nn);
	for (int t = 0; t < Tt; ++t) {
		const double* c = cw.begin() + nn * (std::size_t) t;
		// the covariates of the period with the non-finite entries zeroed
		for (std::size_t i = 0; i < nn; ++i) xs[i] = 1.0;
		for (int k = 0; k < p; ++k) {
			const double* x = X.begin() + nn * ((std::size_t) k + (std::size_t) p * (std::size_t) t);
			double* d = &xs[nn * (std::size_t) (k + 1)];
			for (std::size_t i = 0; i < nn; ++i) d[i] = std::isfinite(x[i]) ? x[i] : 0.0;
		}
		for (int k = 0; k < pa; ++k) {
			const double* a = &xs[nn * (std::size_t) k];
			double* m = mk + nn * (std::size_t) k;
			for (std::size_t i = 0; i < nn; ++i) {
				const double ci = std::isfinite(c[i]) ? c[i] : 0.0;
				xk[i] = a[i] * ci;
				m[i] += xk[i];
			}
			for (int l = k; l < pa; ++l) {
				const double* b = &xs[nn * (std::size_t) l];
				long double v = 0.0L;
				for (int j = 0; j < n; ++j) {
					for (int i = 0; i < n; ++i) v += (long double) (xk[i + (std::size_t) n * j] * b[j + (std::size_t) n * i]);
				}
				xx(k, l) += (double) v;
				if (l != k) xx(l, k) += (double) v;
			}
		}
	}
	return List::create(Named("Mk") = Mk, Named("xx") = xx);
}

// --- the per-row solves of the factor sweeps ------------------------------------------------------------

// Solve n small systems (G_i + L) x_i = rhs_i, G_i given as the rows of `Gflat` (n x R^2, column-major within a
// row) and L a common R x R prior precision: the solutions x (n x R), the inverses Sig (R x R x n), by the
// cholesky factor of each symmetrised system (.ae_mf_row_solve_r in R/ame_als_meanfield.R). `bad` lists the rows
// whose system is not positive definite, which the caller solves by the eigen route of .ae_safe_inverse()
// [[Rcpp::export(rng = false)]]
List ae_row_solve_cpp(const NumericMatrix& Gflat, const NumericMatrix& L, const NumericMatrix& rhs, int threads = 1) {
	const int n = rhs.nrow(), R = rhs.ncol();
	if (Gflat.nrow() != n || Gflat.ncol() != R * R || L.nrow() != R || L.ncol() != R) {
		stop("Incompatible dimensions for the row solves.");
	}
	NumericMatrix x(n, R);
	NumericVector Sig((std::size_t) R * R * (std::size_t) n);
	std::vector<int> badrow(n, 0);
	const int nt = resolve_threads(threads);
	(void) nt;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if (nt > 1)
#endif
	{
		std::vector<double> g((std::size_t) R * R), U((std::size_t) R * R), Ui((std::size_t) R * R), S((std::size_t) R * R);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
		for (int i = 0; i < n; ++i) {
			for (int k = 0; k < R; ++k) {
				for (int l = 0; l < R; ++l) g[k + (std::size_t) R * l] = Gflat(i, k + R * l) + L(k, l);
			}
			for (int k = 0; k < R; ++k) {
				for (int l = k + 1; l < R; ++l) {
					const double v = (g[k + (std::size_t) R * l] + g[l + (std::size_t) R * k]) / 2.0;
					g[k + (std::size_t) R * l] = v;
					g[l + (std::size_t) R * k] = v;
				}
			}
			// upper cholesky factor U, g = U'U
			bool ok = true;
			for (int j = 0; j < R && ok; ++j) {
				double ajj = g[j + (std::size_t) R * j];
				for (int k = 0; k < j; ++k) ajj -= U[k + (std::size_t) R * j] * U[k + (std::size_t) R * j];
				if (!(ajj > 0.0)) { ok = false; break; }
				ajj = std::sqrt(ajj);
				U[j + (std::size_t) R * j] = ajj;
				for (int l = j + 1; l < R; ++l) {
					double v = g[j + (std::size_t) R * l];
					for (int k = 0; k < j; ++k) v -= U[k + (std::size_t) R * j] * U[k + (std::size_t) R * l];
					U[j + (std::size_t) R * l] = v / ajj;
				}
			}
			if (!ok) {
				badrow[i] = 1;
				continue;
			}
			// the inverse of U (upper triangular), then g^-1 = U^-1 U^-T
			for (int j = 0; j < R; ++j) {
				Ui[j + (std::size_t) R * j] = 1.0 / U[j + (std::size_t) R * j];
				for (int k = 0; k < j; ++k) Ui[j + (std::size_t) R * k] = 0.0;
				for (int l = j + 1; l < R; ++l) {
					double v = 0.0;
					for (int k = j; k < l; ++k) v += Ui[j + (std::size_t) R * k] * U[k + (std::size_t) R * l];
					Ui[j + (std::size_t) R * l] = -v / U[l + (std::size_t) R * l];
				}
			}
			for (int k = 0; k < R; ++k) {
				for (int l = k; l < R; ++l) {
					double v = 0.0;
					for (int m = l; m < R; ++m) v += Ui[k + (std::size_t) R * m] * Ui[l + (std::size_t) R * m];
					S[k + (std::size_t) R * l] = v;
					S[l + (std::size_t) R * k] = v;
				}
			}
			for (int k = 0; k < R; ++k) {
				double v = 0.0;
				for (int l = 0; l < R; ++l) {
					v += S[k + (std::size_t) R * l] * rhs(i, l);
					Sig[(std::size_t) k + (std::size_t) R * l + (std::size_t) R * R * (std::size_t) i] = S[k + (std::size_t) R * l];
				}
				x(i, k) = v;
			}
		}
	}
	IntegerVector bad;
	for (int i = 0; i < n; ++i) {
		if (badrow[i]) bad.push_back(i + 1);
	}
	Sig.attr("dim") = IntegerVector::create(R, R, n);
	return List::create(Named("x") = x, Named("Sig") = Sig, Named("bad") = bad);
}
