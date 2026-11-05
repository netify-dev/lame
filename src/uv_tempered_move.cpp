// the tempered transition of R/uv_tempered_move.R, compiled whole: the passes up and down the ladder, the elliptical slice steps of the
// coordinate path of each actor and of the weight of one eigen-direction of a symmetric count fit, at each inverse temperature, the latent
// array carried with the predictor. The Poisson factor of the observed cells is kept as the exponentials E = exp(z + offset) of the
// current state, so that a candidate costs one expm1 per cell: sum [ y dz - E expm1(dz) ] is the change of the log factor for the
// change dz of the predictor. The sums are formed in long double in the order of the R reference (.lame_uv_sym_temper_r()), the
// orders of the actors with R's own sampler (R_unif_index), and the random numbers in the same sequence, so that the compiled and
// the R engines give the same draws from the same seed.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace {

struct Chain {
	int n, N, R, r;
	std::vector<double> z, E, U;
	double Lr;
	const double* Y;
	const int* ok;
	const double* off;
	std::vector<double> ui, uj, nu, e, xp, cy, cE, cu;
	std::vector<int> cnt;
	double ll;
	double evA = 0, callsA = 0, evW = 0, callsW = 0;
	inline size_t ix(int i, int j, int t) const { return (size_t)i + (size_t)j * n + (size_t)t * n * n; }
	inline double& Ur(int i, int t) { return U[(size_t)i + (size_t)r * n + (size_t)t * n * R]; }
};

// the log factor of the carried latent array: the upper-triangle observed cells
static double full_ll(const Chain& C) {
	long double s = 0;
	for (int t = 0; t < C.N; ++t) for (int j = 1; j < C.n; ++j) for (int i = 0; i < j; ++i) {
		const size_t k = C.ix(i, j, t);
		if (!C.ok[k]) continue;
		s += C.Y[k] * C.z[k] - std::exp(C.z[k] + C.off[t]);
	}
	return (double)s;
}

static void actor_step(Chain& C, int a, double b, const std::vector<double>& root, double max_abs, bool rail) {
	const int n = C.n, N = C.N;
	for (int t = 0; t < N; ++t) C.ui[t] = C.Ur(a, t);
	for (int t = 0; t < N; ++t) for (int j = 0; j < n; ++j) C.uj[(size_t)j + (size_t)n * t] = C.Ur(j, t);
	const double Lr = C.Lr;
	// the observed partners of the actor in each period (the matrices are symmetric: the column of the actor is the contiguous line)
	C.cnt.assign(N, 0);
	C.cy.clear(); C.cE.clear(); C.cu.clear();
	for (int t = 0; t < N; ++t) for (int j = 0; j < n; ++j) {
		if (j == a) continue;
		const size_t k = C.ix(j, a, t);
		if (!C.ok[k]) continue;
		C.cy.push_back(C.Y[k]); C.cE.push_back(C.E[k]); C.cu.push_back(C.uj[(size_t)j + (size_t)n * t]); C.cnt[t]++;
	}
	// the change of the log factor for the path x of the actor
	auto diff = [&](const std::vector<double>& x) {
		long double s = 0;
		size_t m = 0;
		for (int t = 0; t < N; ++t) {
			const double du = x[t] - C.ui[t];
			for (int q = 0; q < C.cnt[t]; ++q, ++m) {
				const double dz = Lr * C.cu[m] * du;
				s += C.cy[m] * dz - C.cE[m] * std::expm1(dz);
			}
		}
		return (double)s;
	};
	for (int t = 0; t < N; ++t) C.e[t] = norm_rand();
	for (int t = 0; t < N; ++t) { double m = 0; for (int q = 0; q <= t; ++q) m += root[(size_t)t + (size_t)q * N] * C.e[q]; C.nu[t] = m; }
	const double thr = std::log(unif_rand());
	double th = 2 * M_PI * unif_rand(), lo = th - 2 * M_PI, hi = th;
	int ev = 0;
	bool moved = false;
	double val = 0;
	for (;;) {
		const double ct = std::cos(th), st = std::sin(th);
		for (int t = 0; t < N; ++t) C.xp[t] = C.ui[t] * ct + C.nu[t] * st;
		ev++;
		bool feasible = true;
		if (rail) for (int t = 0; t < N; ++t) if (std::fabs(C.xp[t]) > max_abs) { feasible = false; break; }
		if (feasible) {
			val = diff(C.xp);
			if (std::isfinite(val) && b * val > thr) { moved = true; break; }
		}
		if (th < 0) lo = th; else hi = th;
		if (ev >= 60) break;
		th = lo + (hi - lo) * unif_rand();
	}
	C.evA += ev; C.callsA += 1;
	if (!moved) return;
	// the cells of the line move with the predictor, the diagonal too; the exponentials follow
	for (int t = 0; t < N; ++t) {
		const double du = C.xp[t] - C.ui[t];
		for (int j = 0; j < n; ++j) {
			if (j == a) continue;
			const size_t pos = (size_t)j + (size_t)n * t;
			const double dz = Lr * C.uj[pos] * du;
			const size_t k1 = C.ix(a, j, t), k2 = C.ix(j, a, t);
			C.z[k1] += dz;
			C.z[k2] += dz;
			if (C.ok[k1]) { const double f = 1.0 + std::expm1(dz); C.E[k1] *= f; C.E[k2] *= f; }
		}
		C.z[C.ix(a, a, t)] += Lr * (C.xp[t] * C.xp[t] - C.ui[t] * C.ui[t]);
		C.Ur(a, t) = C.xp[t];
	}
	C.ll += val;
}

static void weight_step(Chain& C, double b, double lvar) {
	const int n = C.n, N = C.N;
	const double L0 = C.Lr;
	std::vector<double> base((size_t)n * n * N);
	for (int t = 0; t < N; ++t) for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) base[C.ix(i, j, t)] = C.Ur(i, t) * C.Ur(j, t);
	// the observed upper-triangle cells, in the order of the sums of the R reference
	std::vector<double> wy, wE, wb;
	for (int t = 0; t < N; ++t) for (int j = 1; j < n; ++j) for (int i = 0; i < j; ++i) {
		const size_t k = C.ix(i, j, t);
		if (!C.ok[k]) continue;
		wy.push_back(C.Y[k]); wE.push_back(C.E[k]); wb.push_back(base[k]);
	}
	const size_t nw = wy.size();
	auto diff = [&](double v) {
		long double s = 0;
		const double c = v - L0;
		for (size_t m = 0; m < nw; ++m) {
			const double dz = c * wb[m];
			s += wy[m] * dz - wE[m] * std::expm1(dz);
		}
		return (double)s;
	};
	const double nuL = std::sqrt(lvar) * norm_rand();
	const double thr = std::log(unif_rand());
	double th = 2 * M_PI * unif_rand(), lo = th - 2 * M_PI, hi = th;
	int ev = 0;
	for (;;) {
		const double xp = L0 * std::cos(th) + nuL * std::sin(th);
		ev++;
		const double val = diff(xp);
		if (std::isfinite(val) && b * val > thr) {
			C.evW += ev; C.callsW += 1;
			const double dv = xp - L0;
			for (int t = 0; t < N; ++t) for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i) {
				const size_t k = C.ix(i, j, t);
				const double dz = dv * base[k];
				C.z[k] += dz;
				if (C.ok[k]) C.E[k] *= (1.0 + std::expm1(dz));
			}
			C.Lr = xp;
			C.ll += val;
			return;
		}
		if (th < 0) lo = th; else hi = th;
		if (ev >= 60) return;
		th = lo + (hi - lo) * unif_rand();
	}
}

// the visiting order of one pass: sample.int(n), then the weight (0) when the pass makes the weight step; reversed on the way down
static std::vector<int> pass_order(int n, bool rev, bool with_weight) {
	std::vector<int> x(n), y(n + (with_weight ? 1 : 0));
	for (int i = 0; i < n; ++i) x[i] = i;
	int m = n;
	for (int i = 0; i < n; ++i) {
		const int j = (int)R_unif_index((double)m);
		y[i] = x[j] + 1;
		x[j] = x[--m];
	}
	if (with_weight) y[n] = 0;
	if (rev) std::reverse(y.begin(), y.end());
	return y;
}

}  // namespace

// One transition. z (n x n x N, symmetric), U (n x R x N), L (R), r (0-based direction), Y the counts, ok (observed cells; symmetric),
// off (log exposure of each period), root (N x N lower Cholesky factor of the AR(1) correlation), lvar (prior variance of the weight), bs
// (inverse temperatures, decreasing from 1), cycles (passes per level), max_abs (bound of the coordinates; Inf for none). Returns the
// proposed state and the log acceptance ratio. The weight is stepped in the first pass of each level on the way up and in the last one on the
// way down, so that the kernel of the way down is the reversal of the kernel of the way up; the uniform of the acceptance is drawn by the caller
// [[Rcpp::export]]
List uv_temper_transition_cpp(NumericVector z, NumericVector U, NumericVector L, int r, NumericVector Y, IntegerVector ok, NumericVector off,
		NumericMatrix root, double lvar, NumericVector bs, int cycles, double max_abs) {
	Chain C;
	C.N = off.size();
	C.R = L.size();
	C.n = (int)std::lround(std::sqrt((double)z.size() / C.N));
	C.r = r;
	C.z.assign(z.begin(), z.end());
	C.U.assign(U.begin(), U.end());
	C.Lr = L[r];
	C.Y = Y.begin(); C.ok = ok.begin(); C.off = off.begin();
	const int n = C.n, N = C.N;
	C.E.assign(z.size(), 0.0);
	for (int t = 0; t < N; ++t) for (int k = 0; k < n * n; ++k) { const size_t q = (size_t)k + (size_t)t * n * n; if (ok[q]) C.E[q] = std::exp(C.z[q] + off[t]); }
	C.ui.resize(N); C.uj.resize((size_t)n * N); C.nu.resize(N); C.e.resize(N); C.xp.resize(N);
	std::vector<double> rt((size_t)N * N);
	for (int a = 0; a < N; ++a) for (int c = 0; c < N; ++c) rt[(size_t)a + (size_t)c * N] = root(a, c);
	const bool rail = std::isfinite(max_abs);
	C.ll = full_ll(C);
	const int K = bs.size() - 1;
	double logacc = 0;
	for (int j = 0; j < K; ++j) {
		logacc += (bs[j + 1] - bs[j]) * C.ll;
		for (int cc = 0; cc < cycles; ++cc) {
			std::vector<int> ord = pass_order(n, false, cc == 0);
			for (size_t s = 0; s < ord.size(); ++s) { if (ord[s] == 0) weight_step(C, bs[j + 1], lvar); else actor_step(C, ord[s] - 1, bs[j + 1], rt, max_abs, rail); }
		}
	}
	for (int j = K - 1; j >= 0; --j) {
		logacc += (bs[j] - bs[j + 1]) * C.ll;
		for (int cc = 0; cc < cycles; ++cc) {
			std::vector<int> ord = pass_order(n, true, cc == cycles - 1);
			for (size_t s = 0; s < ord.size(); ++s) { if (ord[s] == 0) weight_step(C, bs[j + 1], lvar); else actor_step(C, ord[s] - 1, bs[j + 1], rt, max_abs, rail); }
		}
	}
	NumericVector zo(C.z.begin(), C.z.end()), Uo(C.U.begin(), C.U.end());
	zo.attr("dim") = z.attr("dim"); Uo.attr("dim") = U.attr("dim");
	NumericVector Lo = clone(L);
	Lo[r] = C.Lr;
	return List::create(Named("z") = zo, Named("U") = Uo, Named("L") = Lo, Named("logacc") = logacc, Named("ll") = C.ll, Named("evals_actor") = C.evA / std::max(1.0, C.callsA), Named("evals_weight") = C.evW / std::max(1.0, C.callsW));
}
