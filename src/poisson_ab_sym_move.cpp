// The sequential actor-pair loop of one period of the symmetric Poisson dynamic_ab move (.lame_pab_sym_period in
// R/poisson_ab_period_move.R). The random numbers are consumed in the order of the R code (rnorm, runif, runif, then one runif
// per rejected angle, for every pair in turn) and the sums are accumulated in long double as sum() does, so the compiled and the
// R versions give the same draws on the same inputs.
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// one elliptical slice step from N(mu, 1 / prec) for exp(delta A - exp(delta) B - exp(-delta) C); evals counted as .lame_pab_slice does
static double pab_slice1(double mu, double prec, double A, double B, double C, int& evals) {
	auto lf = [&](double d) { return d * A - std::exp(d) * B - std::exp(-d) * C; };
	const double twopi = 2.0 * M_PI;
	double nu = mu + (1.0 / std::sqrt(prec)) * norm_rand();
	double thr = lf(0.0) + std::log(unif_rand());
	double th = R::runif(0.0, twopi);
	double lo = th - twopi, hi = th;
	for (int it = 0; it < 200; ++it) {
		double d = mu - mu * std::cos(th) + (nu - mu) * std::sin(th);
		++evals;
		if (lf(d) > thr) return d;
		if (th < 0) lo = th; else hi = th;
		th = R::runif(lo, hi);
	}
	return 0.0;
}

// [[Rcpp::export]]
List pab_sym_period_cpp(NumericMatrix Zt, NumericMatrix Et, NumericMatrix y, double off, NumericVector eff, NumericVector m, NumericVector v,
		bool centred, IntegerVector perm) {
	const int n = Zt.nrow();
	NumericMatrix Z = clone(Zt), E = clone(Et);
	NumericVector ef = clone(eff);
	std::vector<char> obs((size_t)n * n);
	for (size_t q = 0; q < obs.size(); ++q) obs[q] = std::isfinite(y[q]);
	// the counts and the sum of exp(z + offset) over the cells of actor i other than those with i and k
	auto stats_of = [&](int i, int k, double& s1, double& s2) {
		long double a = 0, b = 0;
		for (int j = 0; j < n; ++j) {
			if (j == i || j == k) continue;
			if (obs[i + (size_t)j * n] && std::isfinite(Z(i, j))) { a += y(i, j); b += std::exp(Z(i, j) + off); }
		}
		s1 = (double)a;
		s2 = (double)b;
	};
	auto apply_shift = [&](int i, double d) {
		for (int j = 0; j < n; ++j) Z(i, j) = Z(i, j) + d;
		for (int j = 0; j < n; ++j) Z(j, i) = Z(j, i) + d;
		for (int j = 0; j < n; ++j) E(i, j) = E(i, j) + d;
		for (int j = 0; j < n; ++j) E(j, i) = E(j, i) + d;
	};
	int pairs = 0, evals = 0;
	const int np = centred ? perm.size() / 2 : perm.size();
	for (int q = 0; q < np; ++q) {
		const int i = perm[q] - 1;
		if (centred) {
			const int k = perm[np + q] - 1;
			double si1, si2, sk1, sk2;
			stats_of(i, k, si1, si2);
			stats_of(k, i, sk1, sk2);
			const double prec = 1.0 / v[i] + 1.0 / v[k];
			const double mu = ((m[i] - ef[i]) / v[i] + (ef[k] - m[k]) / v[k]) / prec;
			const double d = pab_slice1(mu, prec, si1 - sk1, si2, sk2, evals);
			if (d != 0) {
				ef[i] = ef[i] + d;
				ef[k] = ef[k] - d;
				apply_shift(i, d);
				apply_shift(k, -d);
			}
		} else {
			double s1, s2;
			stats_of(i, -1, s1, s2);
			const double d = pab_slice1(m[i] - ef[i], 1.0 / v[i], s1, s2, 0.0, evals);
			ef[i] = ef[i] + d;
			apply_shift(i, d);
		}
		++pairs;
	}
	return List::create(_["Z"] = Z, _["E"] = E, _["eff"] = ef, _["pairs"] = pairs, _["evals"] = evals);
}
