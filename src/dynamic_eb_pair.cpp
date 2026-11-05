// The expected pair problem of the dynamic empirical-Bayes estimator under a wide predictor noise
// (R/lame_dynamic_als_eb_pair.R).
//
// A pair of cells (i, j) and (j, i) has the bivariate probit probability Phi2(s1 eta_1, s2 eta_2; s1 s2 rho) and the
// predictors (eta_1, eta_2) are independent normal N(m_k, sd_k^2) under the variational posterior of the positions
// (or of the static factors). The estimator needs the expectations of 19 functions of (eta_1, eta_2): the log
// probability of the observed pair, its score in the two predictors, the Fisher information over the four outcomes
// (I11, I22, I12), the derivatives of the information in the predictors, and the derivatives in the correlation of the
// log probability, of its scores and of the information (R/ame_als_pair.R, .ae_pair_node_fisher_r()).
//
// The tensor Gauss-Hermite rule of the static estimator (ae_pair_quad_cpp, mode 1) needs a number of nodes in
// proportion to the standard deviation of the noise, because the functions have features of width one: 61 nodes a side
// give the information to 3e-3 and its derivatives to 8e-2 at sd 5.5, 3721 evaluations of the bivariate probability a
// pair. This file integrates on a grid that the pairs SHARE. The nodes are x_k = cc sinh(k du), k = klo .. khi (a
// uniform grid in the mapped variable u, so spacing about cc du near the origin and about du |x| far from it), the
// weight of a node for a pair is jac_k N(x_k; m, sd^2) with jac_k = cc cosh(k du) du, and the 19 functions are
// tabulated ONCE on the grid for the correlation of the fit (one correlation for all pairs and periods):
//
//   E[g] = sum_i sum_j w1_i w2_j T[g][i, j] ,   w_k,i = jac_i N(x_i; m_k, sd_k^2) over the nodes within `cutoff` sd of m_k.
//
// The trapezoid rule in u converges exponentially for an integrand that is analytic in a strip, which these are, so the
// error is set by two ratios: the spacing against the standard deviation of the noise (the Gaussian is resolved when the
// spacing is below sd / 1.25) and the spacing against the width of the features of the functions (about 0.3 for 1e-8).
// A pair costs the contraction of its two weight vectors with the table, 19 multiply-adds a node pair, with no
// evaluation of the bivariate probability; the table costs about 1.5 microseconds a node and is shared by all pairs of
// the call (and by the calls of one evaluation).
//
// Four exported kernels: lde_pair_table_cpp() the table for a correlation, lde_pair_table_blend_cpp() the linear combination of
// tables of one grid (a table for a correlation inside a window is the barycentric combination of the tables at its Chebyshev
// nodes, R/lame_dynamic_als_eb_pairgrid.R), lde_pair_apply_cpp() the contraction, and lde_pair_direct_cpp() the expectation of a
// pair the grid does not resolve (the functions at the nodes of one rule a predictor, no table).
//
// The functions that grow in the tails (the log probability, the scores, the derivatives in the correlation) are
// evaluated in the log domain (R/ame_als_pair.R, the terms of the bivariate probit log probability), so that they
// are exact at predictors of any size (a probability of the four outcomes below the smallest double is a logarithm, never a
// subtraction of two underflowed margins); the Fisher information is the expectation over the four outcomes with the
// probabilities floored at 1e-150 as in the shared kernel, which is exact to rounding wherever it is not below 1e-100.
//
// Column-major conventions as elsewhere: the table is a vector of kTabCols blocks of N * N values, index
// [c * N * N + i * N + j] with i the grid index of the predictor of cell 1 and j that of cell 2.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include "ae_kernels.h"

using namespace Rcpp;
using namespace lame_ae;

namespace {

// table columns: the twelve functions of the Fisher information (I11 I22 I12 I11a I11b I22a I22b I12a I12b I11r I22r
// I12r) and, for each of the four outcomes (+,+) (+,-) (-,+) (-,-) of the observed pair, the seven functions of its log
// probability (lo psi1 psi2 s0 rc0 psi1r psi2r)
const int kInfoCols = 12;
const int kObsCols = 7;
const int kTabCols = kInfoCols + 4 * kObsCols;
const double kSg1[4] = {1.0, 1.0, -1.0, -1.0};
const double kSg2[4] = {1.0, -1.0, 1.0, -1.0};

// the output columns of an expectation and the table column each is read from for the outcome v
inline int out_col(int c, int v) {
	static const int info[19] = {-1, -1, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, -1, -1, -1, -1, 9, 10, 11};
	static const int obs[19] = {0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, 4, 5, 6, -1, -1, -1};
	return info[c] >= 0 ? info[c] : kInfoCols + kObsCols * v + obs[c];
}

// log Phi(t) and Phi(t) together
inline void lphi(double t, double& lp, double& p) {
	if (t > -30.0) {
		p = pnorm(t);
		lp = t > 0.0 ? std::log1p(-pnorm(-t)) : std::log(p);
	} else {
		lp = log_pnorm(t);
		p = std::exp(lp);
	}
}

// the table values at the node (a, b) for the correlation r: out[0..11] the information functions, out[12 + 7 v + ..] the
// observed-outcome functions of the outcomes in `vmask` (the others are left as they are)
inline void table_node(double a, double b, double r, PlanCache& cache, int vmask, double* out) {
	const PairPlan& pp = cache.get(r);
	const double kap = pp.kap, omr = pp.omr;
	const double ta = (a - r * b) / kap, tb = (b - r * a) / kap;
	const double la = dnorm_log(a), lb = dnorm_log(b);
	const double fa = std::exp(la), fb = std::exp(lb);
	const double pa = pnorm(a), pb = pnorm(b);
	// the four outcome probabilities and their logs: one bivariate normal probability and the margins, and a
	// probability below 1e-5 recomputed directly (the subtraction cannot give it)
	double LP[4], P[4];
	LP[0] = log_phi2_plan(a, b, pp);
	P[0] = std::exp(LP[0]);
	P[1] = pa - P[0];
	P[2] = pb - P[0];
	P[3] = 1.0 - pa - pb + P[0];
	// the margins give an outcome probability of at least 1e-5 (the series is accurate to 1e-15 absolute); where the
	// probability P[0] is itself small it comes from the log-domain integral, accurate to a relative 1e-13, and the margin of
	// an outcome that is at least a thousandth of the minuend loses three digits at most
	const bool tail0 = P[0] < 1e-6;
	for (int j = 1; j < 4; ++j) {
		const double minuend = j == 1 ? pa : (j == 2 ? pb : 1.0);
		// (a minuend below 1e-250 is no number to subtract from: both terms may have underflowed to zero, which would pass the
		// comparison below with a probability of zero and a log of minus infinity)
		if (P[j] >= 1e-5 || (tail0 && j < 3 && minuend > 1e-250 && P[j] >= 1e-3 * minuend)) {
			LP[j] = std::log(P[j]);
		} else {
			LP[j] = log_phi2_plan(kSg1[j] * a, kSg2[j] * b, cache.get(kSg1[j] * kSg2[j] * r));
			P[j] = std::exp(LP[j]);
		}
	}
	double Pc[4];
	for (int j = 0; j < 4; ++j) Pc[j] = std::max(P[j], 1e-150);
	// A_j = Phi(s2 tb), B_j = Phi(s1 ta)
	double lPta, lMta, lPtb, lMtb, pta, nta, ptb, ntb;
	lphi(ta, lPta, pta);
	lphi(-ta, lMta, nta);
	lphi(tb, lPtb, ptb);
	lphi(-tb, lMtb, ntb);
	const double A[4] = {ptb, ntb, ptb, ntb};
	const double B[4] = {pta, pta, nta, nta};
	const double ft = dnorm(tb), fta = dnorm(ta);
	const double kap3 = kap * kap * kap;
	const double phi2 = fa * ft / kap;
	double F11 = 0.0, F22 = 0.0, F12 = 0.0;
	double F11a = 0.0, F11b = 0.0, F22a = 0.0, F22b = 0.0, F12a = 0.0, F12b = 0.0;
	double I11r = 0.0, I22r = 0.0, I12r = 0.0;
	for (int j = 0; j < 4; ++j) {
		const double s1 = kSg1[j], s2 = kSg2[j], ss = s1 * s2;
		const double Aa = (-r / kap * ft) * s2;
		const double Ab = (ft / kap) * s2;
		const double Ba = (fta / kap) * s1;
		const double Bb = (-r / kap * fta) * s1;
		const double Pa = (fa * s1) * A[j];
		const double Pb = (fb * s2) * B[j];
		const double P2 = Pc[j] * Pc[j];
		F11 += A[j] * A[j] / Pc[j];
		F22 += B[j] * B[j] / Pc[j];
		F12 += ss * A[j] * B[j] / Pc[j];
		F11a += 2.0 * A[j] * Aa / Pc[j] - A[j] * A[j] * Pa / P2;
		F11b += 2.0 * A[j] * Ab / Pc[j] - A[j] * A[j] * Pb / P2;
		F22a += 2.0 * B[j] * Ba / Pc[j] - B[j] * B[j] * Pa / P2;
		F22b += 2.0 * B[j] * Bb / Pc[j] - B[j] * B[j] * Pb / P2;
		F12a += ss * ((Aa * B[j] + A[j] * Ba) / Pc[j] - A[j] * B[j] * Pa / P2);
		F12b += ss * ((Ab * B[j] + A[j] * Bb) / Pc[j] - A[j] * B[j] * Pb / P2);
		const double Pr = ss * phi2;
		const double Par = ss * (fa * ft * (r * b - a) / kap3);
		const double Pbr = ss * (fb * fta * (r * a - b) / kap3);
		I11r += 2.0 * Pa * Par / Pc[j] - Pa * Pa * Pr / P2;
		I22r += 2.0 * Pb * Pbr / Pc[j] - Pb * Pb * Pr / P2;
		I12r += (Par * Pb + Pa * Pbr) / Pc[j] - Pa * Pb * Pr / P2;
	}
	const double I11 = fa * fa * F11, I22 = fb * fb * F22, I12 = fa * fb * F12;
	out[0] = I11;
	out[1] = I22;
	out[2] = I12;
	out[3] = -2.0 * a * fa * fa * F11 + fa * fa * F11a;
	out[4] = fa * fa * F11b;
	out[5] = fb * fb * F22a;
	out[6] = -2.0 * b * fb * fb * F22 + fb * fb * F22b;
	out[7] = -a * fa * fb * F12 + fa * fb * F12a;
	out[8] = -b * fa * fb * F12 + fa * fb * F12b;
	out[9] = I11r;
	out[10] = I22r;
	out[11] = I12r;
	// a pair so far in the tail that a probability is below the floor has information below 1e-40 and derivatives of the
	// terms of such an outcome that are not numbers: they carry no weight
	for (int c = 0; c < kInfoCols; ++c) {
		if (!std::isfinite(out[c])) out[c] = 0.0;
	}
	// the observed-outcome functions from the log probability of the outcome (no floor: they grow in the tails)
	const double lkap = pp.lkap;
	const double ld_tb = dnorm_log(tb);
	const double nn = a * a - 2.0 * r * a * b + b * b;
	for (int j = 0; j < 4; ++j) {
		if (!(vmask & (1 << j))) continue;
		const double s1 = kSg1[j], s2 = kSg2[j], ss = s1 * s2;
		const double lP = LP[j];
		const double lnb = s2 > 0.0 ? lPtb : lMtb;
		const double lna = s1 > 0.0 ? lPta : lMta;
		const double ga = std::exp(la + lnb - lP);
		const double gb = std::exp(lb + lna - lP);
		const double hab = std::exp(la + ld_tb - lkap - lP);
		const double dl2 = ss * (r * omr + a * b * omr - r * nn) / (omr * omr);
		double* o = out + kInfoCols + kObsCols * j;
		o[0] = lP;
		o[1] = s1 * ga;
		o[2] = s2 * gb;
		o[3] = ss * hab;
		o[4] = hab * dl2 - hab * hab;
		o[5] = ss * hab * (r * b - a) / omr - s2 * ga * hab;
		o[6] = ss * hab * (r * a - b) / omr - s1 * gb * hab;
		for (int c = 0; c < kObsCols; ++c) {
			if (!std::isfinite(o[c])) o[c] = 0.0;
		}
	}
}

// the Gauss-Hermite rule of K nodes for the standard normal (nodes z, weights summing to one) by Newton's method on the
// Hermite polynomial (the recurrence of Numerical Recipes' gauher, in the probabilists' scaling); accurate to rounding for
// the K up to 100 that the rules of the engine use
struct GhRule {
	std::vector<double> z, w;
};

GhRule make_gh(int K) {
	GhRule g;
	g.z.assign(K, 0.0);
	g.w.assign(K, 0.0);
	const double pim4 = 0.7511255444649425;
	const int m = (K + 1) / 2;
	double zz = 0.0, pp = 1.0;
	for (int i = 0; i < m; ++i) {
		if (i == 0) zz = std::sqrt((double) (2 * K + 1)) - 1.85575 * std::pow((double) (2 * K + 1), -0.16667);
		else if (i == 1) zz -= 1.14 * std::pow((double) K, 0.426) / zz;
		else if (i == 2) zz = 1.86 * zz - 0.86 * g.z[0];
		else if (i == 3) zz = 1.91 * zz - 0.91 * g.z[1];
		else zz = 2.0 * zz - g.z[i - 2];
		for (int its = 0; its < 100; ++its) {
			double p1 = pim4, p2 = 0.0, p3;
			for (int j = 0; j < K; ++j) {
				p3 = p2;
				p2 = p1;
				p1 = zz * std::sqrt(2.0 / (j + 1)) * p2 - std::sqrt((double) j / (j + 1)) * p3;
			}
			pp = std::sqrt(2.0 * K) * p2;
			const double z1 = zz;
			zz = z1 - p1 / pp;
			if (std::fabs(zz - z1) < 1e-15) break;
		}
		g.z[i] = zz;
		g.z[K - 1 - i] = -zz;
		g.w[i] = 2.0 / (pp * pp);
		g.w[K - 1 - i] = g.w[i];
	}
	// physicists' (weight exp(-t^2), sum sqrt(pi)) to the probabilists' (weight exp(-z^2 / 2), sum 1)
	double sw = 0.0;
	for (int i = 0; i < K; ++i) {
		g.z[i] *= 1.4142135623730951;
		sw += g.w[i];
	}
	for (int i = 0; i < K; ++i) g.w[i] /= sw;
	// ascending order
	std::reverse(g.z.begin(), g.z.end());
	std::reverse(g.w.begin(), g.w.end());
	return g;
}

}  // namespace

// The table of the 19 functions on the grid x_k = cc sinh(k du), k = klo .. khi, for the correlation rho. `vmask` has bit
// j set for the outcomes (+,+) (+,-) (-,+) (-,-) whose observed-outcome functions are wanted. Returns the table, its
// side N and the arguments.
// [[Rcpp::export(rng = false)]]
List lde_pair_table_cpp(double rho, double cc, double du, int klo, int khi, int vmask = 15, int threads = 1) {
	if (!(std::fabs(rho) < 0.999999)) stop("The correlation must be inside (-1, 1).");
	if (!(cc > 0.0 && du > 0.0) || khi <= klo) stop("An empty or invalid grid.");
	const int N = khi - klo + 1;
	const std::size_t NN = (std::size_t) N * (std::size_t) N;
	std::vector<double> x(N);
	for (int i = 0; i < N; ++i) x[i] = cc * std::sinh((klo + i) * du);
	NumericVector tab((R_xlen_t) (NN * kTabCols));
	double* tp = tab.begin();
	const int nt = resolve_threads(threads);
	(void) nt;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if (nt > 1)
#endif
	{
		PlanCache cache;
		double buf[kTabCols];
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 2)
#endif
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < N; ++j) {
				for (int c = 0; c < kTabCols; ++c) buf[c] = 0.0;
				table_node(x[i], x[j], rho, cache, vmask, buf);
				const std::size_t at = (std::size_t) i * (std::size_t) N + (std::size_t) j;
				for (int c = 0; c < kTabCols; ++c) tp[(std::size_t) c * NN + at] = buf[c];
			}
		}
	}
	return List::create(Named("table") = tab, Named("N") = N, Named("cc") = cc, Named("du") = du, Named("klo") = klo,
	                    Named("khi") = khi, Named("rho") = rho, Named("vmask") = vmask);
}

// The linear combination sum_k w_k T_k of tables on one grid (the interpolation of the table in the correlation: the
// tables of a family at the Chebyshev nodes of a window of correlations, the weights those of the barycentric formula at the
// correlation wanted). Each element is the same sum in the same order whatever the number of threads. The tables are the
// `table` vectors that lde_pair_table_cpp() returns.
// [[Rcpp::export(rng = false)]]
NumericVector lde_pair_table_blend_cpp(const List& tabs, const NumericVector& w, int threads = 1) {
	const int nk = (int) tabs.size();
	if (nk < 1 || w.size() != nk) stop("Incompatible lengths.");
	std::vector<NumericVector> keep(nk);
	std::vector<const double*> tp(nk);
	R_xlen_t len = -1;
	for (int k = 0; k < nk; ++k) {
		keep[k] = tabs[k];
		if (len < 0) len = keep[k].size();
		else if (keep[k].size() != len) stop("Tables of different sizes.");
		tp[k] = keep[k].begin();
	}
	NumericVector out(len);
	double* op = out.begin();
	const R_xlen_t chunk = 2048;
	const R_xlen_t nch = (len + chunk - 1) / chunk;
	const int nt = resolve_threads(threads);
	(void) nt;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) if (nt > 1) schedule(static)
#endif
	for (R_xlen_t ch = 0; ch < nch; ++ch) {
		const R_xlen_t a = ch * chunk;
		const R_xlen_t b = std::min(len, a + chunk);
		const double w0 = w[0];
		const double* t0 = tp[0];
		for (R_xlen_t i = a; i < b; ++i) op[i] = w0 * t0[i];
		for (int k = 1; k < nk; ++k) {
			const double wk = w[k];
			const double* tk = tp[k];
			for (R_xlen_t i = a; i < b; ++i) op[i] += wk * tk[i];
		}
	}
	return out;
}

namespace {

// the dot product of the weights w2 (n2 values) with the rows of one column of the table, weighted by w1: sum_i w1[i] sum_j
// w2[j] row_i[j], with four rows at a time so that the four sums run independently (the dot products of a short window are bound
// by the latency of one running sum, not by the arithmetic)
inline double window_sum(const double* T0, std::size_t N, int n1, int n2, const double* w1, const double* w2) {
	double s = 0.0;
	int i = 0;
	for (; i + 3 < n1; i += 4) {
		const double* r0 = T0 + (std::size_t) i * N;
		const double* r1 = r0 + N;
		const double* r2 = r1 + N;
		const double* r3 = r2 + N;
		double d0 = 0.0, d1 = 0.0, d2 = 0.0, d3 = 0.0;
#ifdef _OPENMP
#pragma omp simd reduction(+ : d0, d1, d2, d3)
#endif
		for (int j = 0; j < n2; ++j) {
			const double wj = w2[j];
			d0 += wj * r0[j];
			d1 += wj * r1[j];
			d2 += wj * r2[j];
			d3 += wj * r3[j];
		}
		s += w1[i] * d0 + w1[i + 1] * d1 + w1[i + 2] * d2 + w1[i + 3] * d3;
	}
	for (; i < n1; ++i) {
		const double* r0 = T0 + (std::size_t) i * N;
		double d0 = 0.0;
#ifdef _OPENMP
#pragma omp simd reduction(+ : d0)
#endif
		for (int j = 0; j < n2; ++j) d0 += w2[j] * r0[j];
		s += w1[i] * d0;
	}
	return s;
}

}  // namespace

// The expectations of the functions of a table for pairs with the predictors N(m1, sd1^2) and N(m2, sd2^2) (independent)
// and the outcome v (0..3 for (+,+) (+,-) (-,+) (-,-)): ell, psi1, psi2, I11, I22, I12 and, with `deriv`, I11a, I11b,
// I22a, I22b, I12a, I12b and, with `rho_deriv`, s0, rc0, psi1r, psi2r, I11r, I22r, I12r (as ae_pair_quad_cpp, mode 1).
// A pair whose window (the nodes within `cutoff` sd of the mean) is not inside the grid, or whose sd is not positive, is
// not computed: `ok` is zero and its values are zero. `K` is the larger of the numbers of grid nodes in the two windows.
// [[Rcpp::export(rng = false)]]
List lde_pair_apply_cpp(const List& table, const NumericVector& m1, const NumericVector& m2, const NumericVector& sd1,
                        const NumericVector& sd2, const IntegerVector& variant, double cutoff = 7.5, bool deriv = true,
                        bool rho_deriv = false, int threads = 1) {
	const int N = as<int>(table["N"]);
	const double cc = as<double>(table["cc"]), du = as<double>(table["du"]);
	const int klo = as<int>(table["klo"]), vmask = as<int>(table["vmask"]);
	const NumericVector tab = table["table"];
	const std::size_t NN = (std::size_t) N * (std::size_t) N;
	const double* tp = tab.begin();
	const std::ptrdiff_t n = m1.size();
	if (m2.size() != n || sd1.size() != n || sd2.size() != n || variant.size() != n) stop("Incompatible lengths.");
	if ((std::size_t) tab.size() != NN * (std::size_t) kTabCols) stop("The table does not fit its grid.");
	rho_deriv = rho_deriv && deriv;
	const int nout = 6 + (deriv ? 6 : 0) + (rho_deriv ? 7 : 0);
	std::vector<NumericVector> out(nout);
	std::vector<double*> op(nout);
	for (int j = 0; j < nout; ++j) {
		out[j] = NumericVector((R_xlen_t) n);
		op[j] = out[j].begin();
	}
	IntegerVector ok((R_xlen_t) n), Kn((R_xlen_t) n);
	int* okp = ok.begin();
	int* knp = Kn.begin();
	// the output columns in the order of the result
	static const int order[19] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
	std::vector<double> x(N), jac(N);
	for (int i = 0; i < N; ++i) {
		const double u = (klo + i) * du;
		x[i] = cc * std::sinh(u);
		jac[i] = cc * std::cosh(u) * du;
	}
	const int nt = resolve_threads(threads);
	(void) nt;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if (nt > 1)
#endif
	{
		std::vector<double> w1, w2;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 16)
#endif
		for (std::ptrdiff_t p = 0; p < n; ++p) {
			const double s1d = sd1[p], s2d = sd2[p];
			const int v = variant[p];
			if (!(s1d > 0.0) || !(s2d > 0.0) || v < 0 || v > 3 || !(vmask & (1 << v))) continue;
			// the window of each predictor on the grid
			const double lo1 = m1[p] - cutoff * s1d, hi1 = m1[p] + cutoff * s1d;
			const double lo2 = m2[p] - cutoff * s2d, hi2 = m2[p] + cutoff * s2d;
			if (lo1 < x.front() || hi1 > x.back() || lo2 < x.front() || hi2 > x.back()) continue;
			const int i0 = (int) (std::lower_bound(x.begin(), x.end(), lo1) - x.begin());
			const int i1 = (int) (std::upper_bound(x.begin(), x.end(), hi1) - x.begin()) - 1;
			const int j0 = (int) (std::lower_bound(x.begin(), x.end(), lo2) - x.begin());
			const int j1 = (int) (std::upper_bound(x.begin(), x.end(), hi2) - x.begin()) - 1;
			if (i1 < i0 || j1 < j0) continue;
			const int n1 = i1 - i0 + 1, n2 = j1 - j0 + 1;
			w1.resize(n1);
			w2.resize(n2);
			const double c1 = kInvSqrt2Pi / s1d, c2 = kInvSqrt2Pi / s2d;
			for (int i = 0; i < n1; ++i) {
				const double z = (x[i0 + i] - m1[p]) / s1d;
				w1[i] = jac[i0 + i] * c1 * std::exp(-0.5 * z * z);
			}
			for (int j = 0; j < n2; ++j) {
				const double z = (x[j0 + j] - m2[p]) / s2d;
				w2[j] = jac[j0 + j] * c2 * std::exp(-0.5 * z * z);
			}
			for (int c = 0; c < nout; ++c) {
				const double* T = tp + (std::size_t) out_col(order[c], v) * NN + (std::size_t) i0 * (std::size_t) N + (std::size_t) j0;
				op[c][p] = window_sum(T, (std::size_t) N, n1, n2, w1.data(), w2.data());
			}
			okp[p] = 1;
			knp[p] = std::max(n1, n2);
		}
	}
	List res = List::create(Named("ell") = out[0], Named("psi1") = out[1], Named("psi2") = out[2], Named("I11") = out[3],
	                        Named("I22") = out[4], Named("I12") = out[5]);
	if (deriv) {
		res["I11a"] = out[6]; res["I11b"] = out[7]; res["I22a"] = out[8];
		res["I22b"] = out[9]; res["I12a"] = out[10]; res["I12b"] = out[11];
	}
	if (rho_deriv) {
		res["s0"] = out[12]; res["rc0"] = out[13]; res["psi1r"] = out[14]; res["psi2r"] = out[15];
		res["I11r"] = out[16]; res["I22r"] = out[17]; res["I12r"] = out[18];
	}
	res["ok"] = ok;
	res["K"] = Kn;
	return res;
}

// The expectations of a pair by the product of two one-dimensional rules and the direct evaluation of the functions at the
// product nodes (no table): for the pairs that the shared grid does not resolve, where the standard deviation of one
// predictor is below the grid spacing at its mean while the other is wide, or where both are narrow. Each predictor has its
// own rule: type 0 is the Gauss-Hermite rule of K nodes, type 1 the nodes of the grid x_k = cc sinh(k du) within `cutoff`
// sd of the mean with the trapezoid weights of the grid (the rule of the table, for a wide predictor); a predictor with a
// standard deviation of zero is the single node at its mean. Same output as lde_pair_apply_cpp().
// [[Rcpp::export(rng = false)]]
List lde_pair_direct_cpp(const NumericVector& m1, const NumericVector& m2, const NumericVector& sd1, const NumericVector& sd2,
                         const IntegerVector& variant, const IntegerVector& type1, const IntegerVector& type2,
                         const IntegerVector& K1, const IntegerVector& K2, double rho, double cc, double du,
                         double cutoff = 7.5, bool deriv = true, bool rho_deriv = false, int threads = 1) {
	if (!(std::fabs(rho) < 0.999999)) stop("The correlation must be inside (-1, 1).");
	const std::ptrdiff_t n = m1.size();
	if (m2.size() != n || sd1.size() != n || sd2.size() != n || variant.size() != n || type1.size() != n || type2.size() != n ||
	    K1.size() != n || K2.size() != n) stop("Incompatible lengths.");
	rho_deriv = rho_deriv && deriv;
	const int nout = 6 + (deriv ? 6 : 0) + (rho_deriv ? 7 : 0);
	// the rules, formed before the loop
	std::vector<GhRule> rules(201);
	std::vector<char> have(201, 0);
	for (std::ptrdiff_t p = 0; p < n; ++p) {
		for (int d = 0; d < 2; ++d) {
			const int k = d == 0 ? K1[p] : K2[p];
			const int ty = d == 0 ? type1[p] : type2[p];
			if (ty == 0 && k > 1) {
				if (k > 200) stop("The number of quadrature nodes must be at most 200.");
				if (!have[k]) {
					rules[k] = make_gh(k);
					have[k] = 1;
				}
			}
		}
	}
	std::vector<NumericVector> out(nout);
	std::vector<double*> op(nout);
	for (int j = 0; j < nout; ++j) {
		out[j] = NumericVector((R_xlen_t) n);
		op[j] = out[j].begin();
	}
	IntegerVector ok((R_xlen_t) n), Kn((R_xlen_t) n);
	int* okp = ok.begin();
	int* knp = Kn.begin();
	static const int order[19] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
	const int nt = resolve_threads(threads);
	(void) nt;
#ifdef _OPENMP
#pragma omp parallel num_threads(nt) if (nt > 1)
#endif
	{
		PlanCache cache;
		std::vector<double> x1, x2, w1, w2;
		double buf[kTabCols];
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 4)
#endif
		for (std::ptrdiff_t p = 0; p < n; ++p) {
			const int v = variant[p];
			if (v < 0 || v > 3) continue;
			bool good = true;
			for (int d = 0; d < 2 && good; ++d) {
				const double mm = d == 0 ? m1[p] : m2[p];
				const double sd = d == 0 ? sd1[p] : sd2[p];
				const int ty = d == 0 ? type1[p] : type2[p];
				const int k = d == 0 ? K1[p] : K2[p];
				std::vector<double>& x = d == 0 ? x1 : x2;
				std::vector<double>& w = d == 0 ? w1 : w2;
				x.clear();
				w.clear();
				if (!(sd > 0.0) || (ty == 0 && k <= 1)) {
					x.push_back(mm);
					w.push_back(1.0);
				} else if (ty == 0) {
					const GhRule& g = rules[k];
					for (int i = 0; i < k; ++i) {
						x.push_back(mm + sd * g.z[i]);
						w.push_back(g.w[i]);
					}
				} else {
					const int a0 = (int) std::ceil(std::asinh((mm - cutoff * sd) / cc) / du);
					const int a1 = (int) std::floor(std::asinh((mm + cutoff * sd) / cc) / du);
					if (a1 < a0) { good = false; break; }
					const double c1 = kInvSqrt2Pi / sd;
					for (int i = a0; i <= a1; ++i) {
						const double u = i * du;
						const double xi = cc * std::sinh(u);
						const double z = (xi - mm) / sd;
						x.push_back(xi);
						w.push_back(cc * std::cosh(u) * du * c1 * std::exp(-0.5 * z * z));
					}
				}
			}
			if (!good) continue;
			double acc[19];
			for (int c = 0; c < nout; ++c) acc[c] = 0.0;
			const int n1 = (int) x1.size(), n2 = (int) x2.size();
			for (int i = 0; i < n1; ++i) {
				for (int j = 0; j < n2; ++j) {
					const double wt = w1[i] * w2[j];
					for (int c = 0; c < kTabCols; ++c) buf[c] = 0.0;
					table_node(x1[i], x2[j], rho, cache, 1 << v, buf);
					for (int c = 0; c < nout; ++c) acc[c] += wt * buf[out_col(order[c], v)];
				}
			}
			for (int c = 0; c < nout; ++c) op[c][p] = acc[c];
			okp[p] = 1;
			knp[p] = std::max(n1, n2);
		}
	}
	List res = List::create(Named("ell") = out[0], Named("psi1") = out[1], Named("psi2") = out[2], Named("I11") = out[3],
	                        Named("I22") = out[4], Named("I12") = out[5]);
	if (deriv) {
		res["I11a"] = out[6]; res["I11b"] = out[7]; res["I22a"] = out[8];
		res["I22b"] = out[9]; res["I12a"] = out[10]; res["I12b"] = out[11];
	}
	if (rho_deriv) {
		res["s0"] = out[12]; res["rc0"] = out[13]; res["psi1r"] = out[14]; res["psi2r"] = out[15];
		res["I11r"] = out[16]; res["I22r"] = out[17]; res["I12r"] = out[18];
	}
	res["ok"] = ok;
	res["K"] = Kn;
	return res;
}
