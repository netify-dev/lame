// The cell of the Poisson-lognormal likelihood of the empirical-Bayes Poisson estimator (ame_als(),
// lame_als() with family = "poisson"): a count y with y | z ~ Poisson(exp(z)), z = eta + e,
// e ~ N(0, ve), the model ame() and lame() sample. The marginal log likelihood of the cell at the
// linear predictor eta,
//
//     l(eta) = log int Poisson(y | exp(eta + e)) N(e; 0, ve) de ,
//
// is concave in eta (the integrand is jointly log concave in (eta, e)), and its derivatives are the
// cumulants of the posterior of the noise e given y,
//
//     pi(e) ~ exp(y e - exp(eta + e) - e^2 / (2 ve)) :
//
//     l'   = E[e] / ve                       (the score)
//     -l'' = Cov(e, exp(eta + e)) / ve       (the observed information; it equals (ve - Var e) / ve^2)
//     -l''' = -kappa3(e) / ve^3              (the derivative of the information in eta)
//
// and dl/dve = (E[e^2] - ve) / (2 ve^2), whence the EM image of ve. The posterior is log concave and
// unimodal: its mode is found by a monotone Newton iteration in zeta = eta + e from a point to the right
// of the root, and the integral is taken by an adaptive Gauss-Hermite rule centred at the mode with the
// Laplace scale. Relative to the Gaussian of the Laplace approximation the integrand is
// exp(-Xs (expm1(t) - t - t^2 / 2)), t = s z, Xs = exp(zeta*), so the weights of the nodes carry only the
// part of the integrand that is not Gaussian. Header only, no R API and no state, so every function is
// safe inside an OpenMP loop. The R reference is .ae_pln_node_r() (R/ame_als_pln.R); the two agree to
// rounding.

#ifndef LAME_AE_PLN_H
#define LAME_AE_PLN_H

#include <algorithm>
#include <cmath>
#include <limits>

namespace lame_ae {

// the largest number of nodes of the inner rule (a stack buffer)
static const int kPlnMaxNodes = 64;

// expm1(t) - t - t^2 / 2, the part of the exponential that the Gaussian approximation leaves out: a series for
// a small argument, where the difference of the three terms loses every digit (and a predictor so large that
// exp(zeta) is 1e40 multiplies what is left by 1e40)
inline double pln_remainder(double t, double expm1_t) {
	if (std::fabs(t) < 0.1) {
		return t * t * t / 6.0 * (1.0 + t / 4.0 * (1.0 + t / 5.0 * (1.0 + t / 6.0 * (1.0 + t / 7.0 * (1.0 + t / 8.0 *
			(1.0 + t / 9.0 * (1.0 + t / 10.0 * (1.0 + t / 11.0))))))));
	}
	return expm1_t - t - 0.5 * t * t;
}

struct PlnNode {
	double l;    // the log likelihood of the cell
	double s;    // the score dl/deta
	double w;    // the observed information -d2l/deta2
	double hp;   // d w / d eta = -d3l/deta3
	double ee;   // E[e^2] under the posterior of the noise, the sufficient statistic of the EM step of ve
	double me;   // E[e]
};

// y: the count; eta: the linear predictor; ve: the variance of the noise; zk, lwk: the nodes of the
// standard normal Gauss-Hermite rule and the logs of its weights (K of them, at most kPlnMaxNodes)
inline PlnNode pln_node(double y, double eta, double ve, int K, const double* zk, const double* lwk,
                        double lgamma_y1) {
	// the linear predictor and the variance are held in a range where the exponentials are finite
	eta = std::min(std::max(eta, -700.0), 700.0);
	ve = std::min(std::max(ve, 1e-8), 1e6);
	const double iv = 1.0 / ve;
	const double ly = y > 0.0 ? std::log(y) : -std::numeric_limits<double>::infinity();
	// the mode: exp(zeta) + zeta / ve = y + eta / ve, a monotone Newton iteration from the right of the root
	// (the left side of the equation is increasing and convex in zeta)
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
	// a node whose weight underflowed to zero carries nothing (and its exponential may overflow)
	double vz = 0.0, k3 = 0.0, eex = 0.0;
	for (int k = 0; k < K; ++k) {
		if (lw[k] == 0.0) continue;
		const double d = zk[k] - ez1;
		vz += lw[k] * d * d;
		k3 += lw[k] * d * d * d;
		eex += lw[k] * ex[k];
	}
	double cze = 0.0;
	for (int k = 0; k < K; ++k) {
		if (lw[k] == 0.0) continue;
		cze += lw[k] * (zk[k] - ez1) * (ex[k] - eex);
	}
	PlnNode o;
	o.me = estar + s * ez1;
	o.s = o.me / ve;
	o.w = Xs * s * cze / ve;
	o.hp = -(s * s * s) * k3 / (ve * ve * ve);
	o.ee = s * s * vz + o.me * o.me;
	const double gs = y * estar - Xs - estar * estar / (2.0 * ve);
	o.l = y * eta - lgamma_y1 - 0.5 * std::log(6.283185307179586476925286766559 * ve) + gs + std::log(s) +
	      0.918938533204672741780329736406 + mx + std::log(tot);
	return o;
}

// --- the pair of cells of a dyad ------------------------------------------------------------------------------
//
// The two counts of a dyad with the noises (e1, e2) bivariate normal, covariance Sigma = [[s11, s12], [s12, s22]]
// (R/ame_als_pln_pair.R has the derivation): the log of the pair probability at the predictor means (m1, m2),
// its gradient P E[e] and the information P Cov(e, exp(m + e)) in the means, P = Sigma^-1, the cumulants of the
// posterior of (e1, e2) given the counts. The mode is found by a damped Newton iteration from the modes of the
// two cells alone, and the integral by the product of two adaptive Gauss-Hermite rules on the Cholesky axes of
// the Laplace covariance.

static const int kPln2MaxNodes = 24;

struct Pln2Node {
	double L, L1, L2, w11, w22, w12;
};

// the mode of one cell alone, in the noise: exp(zeta) + zeta / s = y + m / s, e = zeta - m
inline double pln_mode1(double y, double m, double s) {
	const double ly = y > 0.0 ? std::log(y) : -std::numeric_limits<double>::infinity();
	double zeta = std::max(m, ly);
	const double u = y + m / s;
	const double is = 1.0 / s;
	for (int it = 0; it < 100; ++it) {
		const double ez = std::exp(zeta);
		const double step = (ez + zeta / s - u) / (ez + is);
		zeta -= step;
		if (std::fabs(step) <= 1e-13 * std::max(1.0, std::fabs(zeta))) break;
	}
	return zeta - m;
}

inline Pln2Node pln2_node(double y1, double y2, double m1, double m2, double s11, double s22, double s12, int K,
                          const double* zk, const double* lwk, double lg1, double lg2) {
	m1 = std::min(std::max(m1, -700.0), 700.0);
	m2 = std::min(std::max(m2, -700.0), 700.0);
	s11 = std::min(std::max(s11, 1e-8), 1e6);
	s22 = std::min(std::max(s22, 1e-8), 1e6);
	{
		const double lim = 0.999999 * std::sqrt(s11 * s22);
		s12 = std::max(std::min(s12, lim), -lim);
	}
	const double det = s11 * s22 - s12 * s12;
	const double p11 = s22 / det, p22 = s11 / det, p12 = -s12 / det;
	double e1 = pln_mode1(y1, m1, s11), e2 = pln_mode1(y2, m2, s22);
	// the objective of the noise: sum (y e - exp(m + e)) - e' P e / 2
	#define PLN2_PHI(a, b) (y1 * (a) - std::exp(m1 + (a)) + y2 * (b) - std::exp(m2 + (b)) - 0.5 * (p11 * (a) * (a) + 2.0 * p12 * (a) * (b) + p22 * (b) * (b)))
	double f0 = PLN2_PHI(e1, e2);
	for (int it = 0; it < 100; ++it) {
		const double X1 = std::exp(m1 + e1), X2 = std::exp(m2 + e2);
		const double g1 = y1 - X1 - (p11 * e1 + p12 * e2), g2 = y2 - X2 - (p12 * e1 + p22 * e2);
		const double h11 = X1 + p11, h22 = X2 + p22, h12 = p12;
		const double dt = h11 * h22 - h12 * h12;
		const double d1 = (h22 * g1 - h12 * g2) / dt, d2 = (h11 * g2 - h12 * g1) / dt;
		double t = 1.0;
		bool ok = false;
		double n1 = e1, n2 = e2, fn = f0;
		for (int h = 0; h <= 30; ++h) {
			const double c1 = e1 + t * d1, c2 = e2 + t * d2;
			const double fc = PLN2_PHI(c1, c2);
			if (std::isfinite(fc) && fc >= f0 - 1e-12 * (1.0 + std::fabs(f0))) {
				n1 = c1; n2 = c2; fn = fc; ok = true;
				break;
			}
			t /= 2.0;
		}
		e1 = n1; e2 = n2; f0 = fn;
		if (!ok) break;
		if (std::max(std::fabs(t * d1), std::fabs(t * d2)) <= 1e-12 * std::max(1.0, std::max(std::fabs(e1), std::fabs(e2)))) break;
	}
	#undef PLN2_PHI
	const double X1 = std::exp(m1 + e1), X2 = std::exp(m2 + e2);
	const double h11 = X1 + p11, h22 = X2 + p22, h12 = p12;
	const double dh = h11 * h22 - h12 * h12;
	const double a11 = h22 / dh, a22 = h11 / dh, a12 = -h12 / dh;
	const double c11 = std::sqrt(a11);
	const double c21 = a12 / c11;
	const double c22 = std::sqrt(std::max(a22 - c21 * c21, 1e-300));
	double lw[kPln2MaxNodes * kPln2MaxNodes], D1[kPln2MaxNodes * kPln2MaxNodes], D2[kPln2MaxNodes * kPln2MaxNodes];
	double B1[kPln2MaxNodes * kPln2MaxNodes], B2[kPln2MaxNodes * kPln2MaxNodes];
	double mx = -std::numeric_limits<double>::infinity();
	int q = 0;
	for (int i = 0; i < K; ++i) {
		for (int j = 0; j < K; ++j, ++q) {
			const double d1 = c11 * zk[i];
			const double d2 = c21 * zk[i] + c22 * zk[j];
			const double b1 = std::expm1(d1), b2 = std::expm1(d2);
			D1[q] = d1; D2[q] = d2; B1[q] = b1; B2[q] = b2;
			lw[q] = lwk[i] + lwk[j] - X1 * pln_remainder(d1, b1) - X2 * pln_remainder(d2, b2);
			if (lw[q] > mx) mx = lw[q];
		}
	}
	const int nq = K * K;
	double tot = 0.0;
	for (int k = 0; k < nq; ++k) { lw[k] = std::exp(lw[k] - mx); tot += lw[k]; }
	double E1 = 0.0, E2 = 0.0, x1m = 0.0, x2m = 0.0;
	for (int k = 0; k < nq; ++k) {
		lw[k] /= tot;
		E1 += lw[k] * D1[k]; E2 += lw[k] * D2[k];
	}
	for (int k = 0; k < nq; ++k) {
		if (lw[k] == 0.0) continue;
		x1m += lw[k] * B1[k]; x2m += lw[k] * B2[k];
	}
	double s1u1 = 0.0, s1u2 = 0.0, s2u1 = 0.0, s2u2 = 0.0;
	for (int k = 0; k < nq; ++k) {
		if (lw[k] == 0.0) continue;
		const double q1 = D1[k] - E1, q2 = D2[k] - E2;
		const double u1 = B1[k] - x1m, u2 = B2[k] - x2m;
		s1u1 += lw[k] * q1 * u1; s1u2 += lw[k] * q1 * u2;
		s2u1 += lw[k] * q2 * u1; s2u2 += lw[k] * q2 * u2;
	}
	const double M11 = X1 * s1u1, M12 = X2 * s1u2, M21 = X1 * s2u1, M22 = X2 * s2u2;
	const double Ee1 = e1 + E1, Ee2 = e2 + E2;
	const double gs = y1 * e1 - X1 + y2 * e2 - X2 - 0.5 * (p11 * e1 * e1 + 2.0 * p12 * e1 * e2 + p22 * e2 * e2);
	Pln2Node o;
	o.L = y1 * m1 + y2 * m2 - lg1 - lg2 - 0.5 * std::log(det / (a11 * a22 - a12 * a12)) + gs + mx + std::log(tot);
	o.L1 = p11 * Ee1 + p12 * Ee2;
	o.L2 = p12 * Ee1 + p22 * Ee2;
	o.w11 = p11 * M11 + p12 * M21;
	o.w22 = p12 * M12 + p22 * M22;
	o.w12 = 0.5 * ((p11 * M12 + p12 * M22) + (p12 * M11 + p22 * M21));
	return o;
}

}  // namespace lame_ae

#endif
