// Kernels of the moves of R/poisson_period_move.R: the poisson factor of a translation, and the
// tilted residual carry of the dyadic variance and correlation (.lame_ppm_scale).
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// sum(ez * expm1(d * x)) over the cells of a translation: the part of the log poisson factor of a shift d
// of a coefficient with design x that is not linear in d (ez = exp(z + offset) on the observed cells)
// [[Rcpp::export]]
double ppm_shift_sum_cpp(const NumericVector& ez, const NumericVector& x, double d) {
	const R_xlen_t n = ez.size();
	if (x.size() != n) Rcpp::stop("ppm_shift_sum_cpp: the cells and the design differ in length");
	double s = 0.0;
	for (R_xlen_t i = 0; i < n; ++i) s += ez[i] * std::expm1(d * x[i]);
	return s;
}

// The state of the carried residuals of the dyads at (s2, rho): the poisson term, the gaussian term and the log
// determinant of the conditional covariance of the dyads. eu, el are the residuals of the upper and lower cell of each
// dyad, wu, wl the first-order weights lambda of the observed cells (0 where missing), yu, yl their counts (0 where
// missing). The poisson term drops the part that does not depend on the residuals
// [[Rcpp::export]]
List ppm_dyad_state_cpp(const NumericVector& eu, const NumericVector& el, const NumericVector& wu, const NumericVector& wl,
                        const NumericVector& yu, const NumericVector& yl, double s2, double rho) {
	const R_xlen_t n = eu.size();
	const double q = 1.0 / (s2 * (1.0 - rho * rho));
	const double a12 = -rho * q;
	double ll = 0.0, quad = 0.0, ld = 0.0;
	for (R_xlen_t i = 0; i < n; ++i) {
		ll += yu[i] * eu[i] - wu[i] * std::exp(eu[i]) + yl[i] * el[i] - wl[i] * std::exp(el[i]);
		quad += eu[i] * eu[i] - 2.0 * rho * eu[i] * el[i] + el[i] * el[i];
		const double a11 = q + wu[i], a22 = q + wl[i];
		ld += -0.5 * std::log(a11 * a22 - a12 * a12);
	}
	const double g = -(double) n * (std::log(s2) + 0.5 * std::log1p(-rho * rho)) - quad / (2.0 * s2 * (1.0 - rho * rho));
	return List::create(Named("ll") = ll, Named("g") = g, Named("logdet") = ld);
}

// One proposal of the carried residuals of the dyads: the residuals (eu, el) at (s2c, rhoc) mapped to the residuals at
// (s2n, rhon) by e' = m' + L' L^-1 (e - m), with V = (Sigma^-1 + diag(w))^-1 = L L', m = V b the conditional mean
// (b = y - lambda, w = lambda on the observed cells, 0 where missing). Returns the mapped residuals and the poisson term,
// the gaussian term and the log determinant at the proposal
// [[Rcpp::export]]
List ppm_dyad_carry_cpp(const NumericVector& eu, const NumericVector& el, const NumericVector& wu, const NumericVector& wl,
                        const NumericVector& bu, const NumericVector& bl, const NumericVector& yu, const NumericVector& yl,
                        double s2c, double rhoc, double s2n, double rhon) {
	const R_xlen_t n = eu.size();
	NumericVector eun(n), eln(n);
	const double qc = 1.0 / (s2c * (1.0 - rhoc * rhoc)), a12c = -rhoc * qc;
	const double qn = 1.0 / (s2n * (1.0 - rhon * rhon)), a12n = -rhon * qn;
	double ll = 0.0, quad = 0.0, ld = 0.0;
	for (R_xlen_t i = 0; i < n; ++i) {
		// current conditional law: V = A^-1, A = [[a11, a12], [a12, a22]]; L = [[l11, 0], [l21, l22]] with l22 = 1 / sqrt(a22)
		const double a11c = qc + wu[i], a22c = qc + wl[i];
		const double dtc = a11c * a22c - a12c * a12c;
		const double muc = (a22c * bu[i] - a12c * bl[i]) / dtc;
		const double mlc = (a11c * bl[i] - a12c * bu[i]) / dtc;
		const double l11c = std::sqrt(a22c / dtc), l21c = (-a12c / dtc) / l11c, l22c = 1.0 / std::sqrt(a22c);
		const double z1 = (eu[i] - muc) / l11c;
		const double z2 = (el[i] - mlc - l21c * z1) / l22c;
		// the conditional law at the proposal
		const double a11n = qn + wu[i], a22n = qn + wl[i];
		const double dtn = a11n * a22n - a12n * a12n;
		const double mun = (a22n * bu[i] - a12n * bl[i]) / dtn;
		const double mln = (a11n * bl[i] - a12n * bu[i]) / dtn;
		const double l11n = std::sqrt(a22n / dtn), l21n = (-a12n / dtn) / l11n, l22n = 1.0 / std::sqrt(a22n);
		const double u = mun + l11n * z1;
		const double l = mln + l21n * z1 + l22n * z2;
		eun[i] = u;
		eln[i] = l;
		ll += yu[i] * u - wu[i] * std::exp(u) + yl[i] * l - wl[i] * std::exp(l);
		quad += u * u - 2.0 * rhon * u * l + l * l;
		ld += -0.5 * std::log(dtn);
	}
	const double g = -(double) n * (std::log(s2n) + 0.5 * std::log1p(-rhon * rhon)) - quad / (2.0 * s2n * (1.0 - rhon * rhon));
	return List::create(Named("eu") = eun, Named("el") = eln, Named("ll") = ll, Named("g") = g, Named("logdet") = ld);
}

// The same for one latent value per cell (a bipartite network, a directed network without a dyadic correlation) or per
// unordered pair (a symmetric network): residuals e at variance s2c, conditional law N(m, v) with
// v = 1 / (1 / s2 + w), m = v b
// [[Rcpp::export]]
List ppm_cell_state_cpp(const NumericVector& e, const NumericVector& w, const NumericVector& y, double s2) {
	const R_xlen_t n = e.size();
	const double q = 1.0 / s2;
	double ll = 0.0, ss = 0.0, ld = 0.0;
	for (R_xlen_t i = 0; i < n; ++i) {
		ll += y[i] * e[i] - w[i] * std::exp(e[i]);
		ss += e[i] * e[i];
		ld += -0.5 * std::log(q + w[i]);
	}
	return List::create(Named("ll") = ll, Named("g") = -0.5 * (double) n * std::log(s2) - ss / (2.0 * s2), Named("logdet") = ld);
}

// [[Rcpp::export]]
List ppm_cell_carry_cpp(const NumericVector& e, const NumericVector& w, const NumericVector& b, const NumericVector& y,
                        double s2c, double s2n) {
	const R_xlen_t n = e.size();
	NumericVector en(n);
	const double qc = 1.0 / s2c, qn = 1.0 / s2n;
	double ll = 0.0, ss = 0.0, ld = 0.0;
	for (R_xlen_t i = 0; i < n; ++i) {
		const double ac = qc + w[i], an = qn + w[i];
		const double mc = b[i] / ac, mn = b[i] / an;
		const double v = mn + std::sqrt(ac / an) * (e[i] - mc);
		en[i] = v;
		ll += y[i] * v - w[i] * std::exp(v);
		ss += v * v;
		ld += -0.5 * std::log(an);
	}
	return List::create(Named("e") = en, Named("ll") = ll, Named("g") = -0.5 * (double) n * std::log(s2n) - ss / (2.0 * s2n),
	                    Named("logdet") = ld);
}
