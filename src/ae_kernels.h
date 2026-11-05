// Scalar kernels of the binary probit empirical-bayes estimator (ame_als(),
// lame_als()): the normal log cdf, the bivariate normal probability, the
// log-domain tail of that probability, and the terms of the bivariate probit
// log probability of a pair of cells. Header only, no R API, no allocation
// after the first call and no state that a call changes, so every function is
// safe inside an OpenMP loop. The R references are .ae_pair_logphi2_r() and
// .ae_pair_logp_r() (R/ame_als_pair.R); the engines agree to rounding.
//
// The bivariate normal probability is the scheme of Drezner and Wesolowsky
// (1989) with the modifications of Genz (2004) for a correlation near one
// (the scheme of Genz's MVBVU, which pbivnorm::pbivnorm() uses): a Gauss
// Legendre rule of 6, 12 or 20 points on the angle for |r| < 0.925, and a
// series about the boundary for the rest. It is accurate to about 1e-15 in
// absolute terms, so a probability below 1e-6 is taken in the log domain.

#ifndef LAME_AE_KERNELS_H
#define LAME_AE_KERNELS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace lame_ae {

// the number of OpenMP threads a kernel uses: the number asked for, at most the processors there are, and
// one where the platform has no OpenMP
inline int resolve_threads(int threads) {
#ifdef _OPENMP
	const int nt = std::max(1, threads);
	return std::max(1, std::min(nt, omp_get_num_procs()));
#else
	(void) threads;
	return 1;
#endif
}

static const double kLogSqrt2Pi = 0.918938533204672741780329736406;
static const double kInvSqrt2 = 0.707106781186547524400844362105;
static const double kTwoPi = 6.283185307179586476925286766559;
static const double kInvSqrt2Pi = 0.398942280401432677939946059934;

// --- the normal distribution -------------------------------------------------------

inline double dnorm_log(double x) { return -0.5 * x * x - kLogSqrt2Pi; }
inline double dnorm(double x) { return kInvSqrt2Pi * std::exp(-0.5 * x * x); }
inline double pnorm(double x) { return 0.5 * std::erfc(-x * kInvSqrt2); }

// log Phi(x) for any x. The upper half is log(1 - Phi(-x)) by log1p, the lower
// half is the log of erfc down to x = -30 and the asymptotic series beyond
inline double log_pnorm(double x) {
	if (x != x) return x;
	if (x > 0.0) return std::log1p(-0.5 * std::erfc(x * kInvSqrt2));
	if (x > -30.0) return std::log(0.5 * std::erfc(-x * kInvSqrt2));
	const double u = 1.0 / (x * x);
	const double s = u * (-1.0 + u * (3.0 + u * (-15.0 + u * (105.0 + u * (-945.0 + u * (10395.0 +
		u * (-135135.0 + u * 2027025.0)))))));
	return -0.5 * x * x - std::log(-x) - kLogSqrt2Pi + std::log1p(s);
}

// log of the logistic cdf, -log(1 + exp(-x)), by the pieces R uses for log1pexp
inline double log_plogis(double x) {
	const double t = -x;
	double l;
	if (t <= 18.0) l = std::log1p(std::exp(t));
	else if (t > 33.3) l = t;
	else l = t + std::exp(-t);
	return -l;
}

// --- Gauss Legendre rules ---------------------------------------------------------------

// the n-point rule on (-1, 1) by Newton's method on the Legendre polynomial; x ascending
inline void gauss_legendre(int n, double* x, double* w) {
	const int m = (n + 1) / 2;
	for (int i = 0; i < m; ++i) {
		double z = std::cos(3.14159265358979323846 * (i + 0.75) / (n + 0.5));
		double pp = 1.0;
		for (int it = 0; it < 100; ++it) {
			double p1 = 1.0, p2 = 0.0;
			for (int j = 1; j <= n; ++j) {
				const double p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
			}
			pp = n * (z * p1 - p2) / (z * z - 1.0);
			const double z1 = z;
			z = z1 - p1 / pp;
			if (std::fabs(z - z1) < 1e-16) break;
		}
		// the derivative at the final abscissa
		double p1 = 1.0, p2 = 0.0;
		for (int j = 1; j <= n; ++j) {
			const double p3 = p2;
			p2 = p1;
			p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
		}
		pp = n * (z * p1 - p2) / (z * z - 1.0);
		x[i] = -z;
		x[n - 1 - i] = z;
		w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * pp * pp);
	}
}

struct GlTables {
	// the half rules (negative abscissae, ascending) of 6, 12 and 20 points, to the 16 digits at which
	// they are tabulated
	double x6[3], w6[3], x12[6], w12[6], x20[10], w20[10];
	// the 64-point rule of the log-domain tail, its weights times 7 (log), and the distances 14 (1 - (x + 1) / 2)
	// of its points below the upper limit
	double x64[64], w64[64], lw64[64], u64[64];
	GlTables() {
		static const double X6[3] = {-0.9324695142031522, -0.6612093864662647, -0.2386191860831970};
		static const double W6[3] = {0.1713244923791705, 0.3607615730481384, 0.4679139345726904};
		static const double X12[6] = {-0.9815606342467191, -0.9041172563704750, -0.7699026741943050,
			-0.5873179542866171, -0.3678314989981802, -0.1252334085114692};
		static const double W12[6] = {0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
			0.2031674267230659, 0.2334925365383547, 0.2491470458134029};
		static const double X20[10] = {-0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
			-0.8391169718222188, -0.7463319064601508, -0.6360536807265150, -0.5108670019508271,
			-0.3737060887154196, -0.2277858511416451, -0.07652652113349733};
		static const double W20[10] = {0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
			0.08327674157670475, 0.1019301198172404, 0.1181945319615184, 0.1316886384491766,
			0.1420961093183821, 0.1491729864726037, 0.1527533871307259};
		for (int i = 0; i < 3; ++i) { x6[i] = X6[i]; w6[i] = W6[i]; }
		for (int i = 0; i < 6; ++i) { x12[i] = X12[i]; w12[i] = W12[i]; }
		for (int i = 0; i < 10; ++i) { x20[i] = X20[i]; w20[i] = W20[i]; }
		gauss_legendre(64, x64, w64);
		for (int i = 0; i < 64; ++i) {
			lw64[i] = std::log(7.0 * w64[i]);
			u64[i] = 14.0 * (1.0 - (x64[i] + 1.0) / 2.0);
		}
	}
};

inline const GlTables& gl_tables() {
	static const GlTables t;
	return t;
}

// --- the bivariate normal probability ---------------------------------------------------------

// the normal cdf of the reference scheme for the bivariate probability: the Chebyshev expansion of
// Schonfelder (Math. Comp. 32, 1978), accurate to 1e-15. Inside bvu() it is used in place of pnorm()
// so that the bivariate probability agrees with pbivnorm::pbivnorm() to the last digits, including
// where the probability is the small difference of terms of order one (a negative correlation)
inline double phid_ref(double z) {
	static const double A[25] = {
		6.10143081923200417926465815756e-1, -4.34841272712577471828182820888e-1,
		1.76351193643605501125840298123e-1, -6.0710795609249414860051215825e-2,
		1.7712068995694114486147141191e-2, -4.321119385567293818599864968e-3,
		8.54216676887098678819832055e-4, -1.27155090609162742628893940e-4,
		1.1248167243671189468847072e-5, 3.13063885421820972630152e-7,
		-2.70988068537762022009086e-7, 3.0737622701407688440959e-8, 2.515620384817622937314e-9,
		-1.028929921320319127590e-9, 2.9944052119949939363e-11, 2.6051789687266936290e-11,
		-2.634839924171969386e-12, -6.43404509890636443e-13, 1.12457401801663447e-13,
		1.7281533389986098e-14, -4.264101694942375e-15, -5.45371977880191e-16,
		1.58697607761671e-16, 2.0899837844334e-17, -5.900526869409e-18};
	const double xa = std::fabs(z) / 1.414213562373095048801688724209;
	double p;
	if (xa > 100.0) {
		p = 0.0;
	} else {
		const double t = (8.0 * xa - 30.0) / (4.0 * xa + 15.0);
		double bm = 0.0, b = 0.0, bp = 0.0;
		for (int i = 24; i >= 0; --i) {
			bp = b;
			b = bm;
			bm = t * b - bp + A[i];
		}
		p = std::exp(-xa * xa) * (bm - bp) / 4.0;
	}
	if (z > 0.0) p = 1.0 - p;
	return p;
}

// What the bivariate probability needs of the correlation alone: the points of the rule on the
// angle (kind 0, |r| below 0.925), or the points of the series about the boundary (kind 1); kind 2
// is |r| >= 1. Formed once, it serves the many probabilities of one correlation: the nodes of a pair
// and the pairs of a fit share |rho|. The cutoffs of the rule are single precision constants in the
// reference code: 0.3 is 0.300000011920929 and 0.925 is 0.925000011920929 there.
struct BvuPlan {
	double r;
	int kind, lg;
	const double *xr, *wr;
	double asr;
	double sn[20], om[20];
	double as, a_full, ah;
	double xs[20], rs[20], twoxs[20], onepr[20], omrs[20], twoopr[20];
};

inline void bvu_make_plan(double r, BvuPlan& p) {
	static const double kCut1 = (double) 0.3f, kCut3 = (double) 0.925f;
	const GlTables& g = gl_tables();
	const double ar = std::fabs(r);
	p.r = r;
	if (ar < kCut1) { p.xr = g.x6; p.wr = g.w6; p.lg = 3; }
	else if (ar < 0.75) { p.xr = g.x12; p.wr = g.w12; p.lg = 6; }
	else { p.xr = g.x20; p.wr = g.w20; p.lg = 10; }
	if (ar < kCut3) {
		p.kind = 0;
		p.asr = std::asin(r);
		for (int i = 0; i < p.lg; ++i) {
			double sn = std::sin(p.asr * (p.xr[i] + 1.0) / 2.0);
			p.sn[2 * i] = sn;
			p.om[2 * i] = 1.0 - sn * sn;
			sn = std::sin(p.asr * (-p.xr[i] + 1.0) / 2.0);
			p.sn[2 * i + 1] = sn;
			p.om[2 * i + 1] = 1.0 - sn * sn;
		}
		return;
	}
	if (ar < 1.0) {
		p.kind = 1;
		p.as = (1.0 - r) * (1.0 + r);
		p.a_full = std::sqrt(p.as);
		p.ah = p.a_full / 2.0;
		for (int i = 0; i < p.lg; ++i) {
			double xs = p.ah * (p.xr[i] + 1.0);
			xs *= xs;
			double rs = std::sqrt(1.0 - xs);
			p.xs[2 * i] = xs;
			p.rs[2 * i] = rs;
			p.twoxs[2 * i] = 2.0 * xs;
			p.onepr[2 * i] = 1.0 + rs;
			const double u = 1.0 - p.xr[i];
			xs = p.as * (u * u) / 4.0;
			rs = std::sqrt(1.0 - xs);
			p.xs[2 * i + 1] = xs;
			p.rs[2 * i + 1] = rs;
			p.omrs[2 * i + 1] = 1.0 - rs;
			p.twoopr[2 * i + 1] = 2.0 * (1.0 + rs);
		}
		return;
	}
	p.kind = 2;
}

// P(X > h, Y > k) for the standard bivariate normal with the correlation of the plan; for kind 0
// the caller passes the marginals pmh = phid_ref(-h) and pmk = phid_ref(-k), which it can share
// among probabilities with a common threshold
inline double bvu_plan(double h, double k, const BvuPlan& p, double pmh, double pmk) {
	const double r = p.r;
	double hk = h * k;
	double bvn = 0.0;
	if (p.kind == 0) {
		const double hs = (h * h + k * k) / 2.0;
		for (int i = 0; i < p.lg; ++i) {
			bvn += p.wr[i] * std::exp((p.sn[2 * i] * hk - hs) / p.om[2 * i]);
			bvn += p.wr[i] * std::exp((p.sn[2 * i + 1] * hk - hs) / p.om[2 * i + 1]);
		}
		return bvn * p.asr / (2.0 * kTwoPi) + pmh * pmk;
	}
	if (r < 0.0) {
		k = -k;
		hk = -hk;
	}
	if (p.kind == 1) {
		const double bs = (h - k) * (h - k);
		const double c = (4.0 - hk) / 8.0;
		const double d = (12.0 - hk) / 16.0;
		bvn = p.a_full * std::exp(-(bs / p.as + hk) / 2.0) *
			(1.0 - c * (bs - p.as) * (1.0 - d * bs / 5.0) / 3.0 + c * d * p.as * p.as / 5.0);
		if (hk > -160.0) {
			const double b = std::sqrt(bs);
			bvn -= std::exp(-hk / 2.0) * std::sqrt(kTwoPi) * phid_ref(-b / p.a_full) * b *
				(1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
		}
		for (int i = 0; i < p.lg; ++i) {
			double xs = p.xs[2 * i];
			bvn += p.ah * p.wr[i] * (std::exp(-bs / p.twoxs[2 * i] - hk / p.onepr[2 * i]) / p.rs[2 * i] -
				std::exp(-(bs / xs + hk) / 2.0) * (1.0 + c * xs * (1.0 + d * xs)));
			xs = p.xs[2 * i + 1];
			bvn += p.ah * p.wr[i] * std::exp(-(bs / xs + hk) / 2.0) *
				(std::exp(-hk * p.omrs[2 * i + 1] / p.twoopr[2 * i + 1]) / p.rs[2 * i + 1] -
				 (1.0 + c * xs * (1.0 + d * xs)));
		}
		bvn = -bvn / kTwoPi;
	}
	if (r > 0.0) bvn += phid_ref(-std::max(h, k));
	else if (r < 0.0) bvn = -bvn + std::max(0.0, phid_ref(-h) - phid_ref(-k));
	return bvn;
}

// P(X > h, Y > k) for the standard bivariate normal with correlation r
inline double bvu(double h, double k, double r) {
	BvuPlan p;
	bvu_make_plan(r, p);
	const double pmh = p.kind == 0 ? phid_ref(-h) : 0.0;
	const double pmk = p.kind == 0 ? phid_ref(-k) : 0.0;
	return bvu_plan(h, k, p, pmh, pmk);
}

// Phi2(a, b; r) = P(X <= a, Y <= b)
inline double phi2_prob(double a, double b, double r) { return bvu(-a, -b, r); }

// log Phi2(a, b; r) of a probability below 1e-6, where the series is not accurate relative to
// the value: the log of int_{a - 14}^{a} phi(x) Phi((b - r x) / sqrt(1 - r^2)) dx by the 64-point rule.
// The terms are positive and the cdf is at most one, so log phi(x_i) + log w_i bounds the log of the
// term of a point from above without the cdf. The sum is formed relative to the largest bound, with
// the cdf itself (no log of it): a point whose bound is below 1e-20 of the term of the point of the
// largest bound is below 1e-20 of the sum and is not evaluated, so that the result is the sum of the
// 64 points to rounding. Where the cdf at the point of the largest bound is itself tiny, the sum
// is formed in the log domain
inline double log_phi2_tail(double a, double b, double r, double kap) {
	const GlTables& g = gl_tables();
	const double inf = std::numeric_limits<double>::infinity();
	double c[64], xs[64];
	double cmax = -inf;
	int imax = 0;
	for (int i = 0; i < 64; ++i) {
		xs[i] = a - g.u64[i];
		c[i] = dnorm_log(xs[i]) + g.lw64[i];
		if (c[i] > cmax || c[i] != c[i]) {
			cmax = c[i];
			imax = i;
		}
	}
	if (cmax != cmax) return cmax;
	const double ik = 1.0 / kap;
	const double zmax = (b - r * xs[imax]) * ik;
	const double pmax = pnorm(zmax);
	if (pmax > 1e-250) {
		const double cut = cmax + std::log(1e-20 * pmax);
		double s = pmax;
		for (int i = 0; i < 64; ++i) {
			if (i == imax || c[i] < cut) continue;
			s += std::exp(c[i] - cmax) * pnorm((b - r * xs[i]) * ik);
		}
		return cmax + std::log(s);
	}
	double lf[64];
	double mx = -inf;
	for (int i = 0; i < 64; ++i) {
		lf[i] = c[i] + log_pnorm((b - r * xs[i]) / kap);
		if (lf[i] > mx || lf[i] != lf[i]) mx = lf[i];
	}
	if (!(mx > -inf) || mx != mx) return mx;
	double s = 0.0;
	for (int i = 0; i < 64; ++i) s += std::exp(lf[i] - mx);
	return mx + std::log(s);
}

// everything that depends on the correlation alone, for the pair functions
struct PairPlan {
	BvuPlan bv;
	double r, omr, kap, lkap, rkap;
};

inline void pair_make_plan(double r, PairPlan& pp) {
	pp.r = r;
	pp.omr = 1.0 - r * r;
	pp.kap = std::sqrt(pp.omr);
	pp.lkap = std::log(pp.kap);
	pp.rkap = r / pp.kap;
	bvu_make_plan(r, pp.bv);
}

// the plans of the two signs of a correlation, kept while a loop runs over pairs whose correlations
// are +r or -r (the sign of the outcomes of the pair flips it)
struct PlanCache {
	PairPlan slot[2];
	bool has[2];
	PlanCache() { has[0] = has[1] = false; }
	const PairPlan& get(double r) {
		const int j = r < 0.0 ? 1 : 0;
		if (!has[j] || slot[j].r != r) {
			pair_make_plan(r, slot[j]);
			has[j] = true;
		}
		return slot[j];
	}
};

// log Phi2(a, b; r): the series where the probability is at least 1e-6, the log-domain integral
// where it is smaller or not finite; `pa_ref` is phid_ref(a) when the plan is of kind 0 (the caller
// shares it among the nodes with a common a)
inline double log_phi2_pre(double a, double b, const PairPlan& pp, double pa_ref) {
	const double pb_ref = pp.bv.kind == 0 ? phid_ref(b) : 0.0;
	const double p = bvu_plan(-a, -b, pp.bv, pa_ref, pb_ref);
	// (the integral is over the smaller of the two limits: Phi2 is symmetric in them, and the mass of the integrand
	// sits below the smaller one, so a window under the larger limit can miss it altogether)
	if (!std::isfinite(p) || p < 1e-6) return a <= b ? log_phi2_tail(a, b, pp.r, pp.kap) : log_phi2_tail(b, a, pp.r, pp.kap);
	return std::log(std::max(p, 1e-300));
}

inline double log_phi2_plan(double a, double b, const PairPlan& pp) {
	return log_phi2_pre(a, b, pp, pp.bv.kind == 0 ? phid_ref(a) : 0.0);
}

inline double log_phi2(double a, double b, double r) {
	PairPlan pp;
	pair_make_plan(r, pp);
	return log_phi2_plan(a, b, pp);
}

// --- the log probability of a pair and its derivatives --------------------------------------------

struct PairTerms {
	// log P, the gradient (La, Lb), the hessian (Laa, Lbb, Lab), the score of r (Lr, which is
	// phi2 / P) and its derivative (Lrr)
	double L, La, Lb, Laa, Lbb, Lab, Lr, Lrr;
};

// log P, P = Phi2(a, b; r), and its derivatives in the thresholds and the correlation, at the
// plan of r (`pa_ref` as in log_phi2_pre)
inline void pair_logp_pre(double a, double b, const PairPlan& pp, double pa_ref, PairTerms& o) {
	const double r = pp.r, omr = pp.omr, kap = pp.kap;
	const double lP = log_phi2_pre(a, b, pp, pa_ref);
	const double ua = (a - r * b) / kap, ub = (b - r * a) / kap;
	const double la = dnorm_log(a), lb = dnorm_log(b);
	const double lna = log_pnorm(ua), lnb = log_pnorm(ub);
	const double lda = dnorm_log(ua), ldb = dnorm_log(ub);
	const double ga = std::exp(la + lnb - lP), gb = std::exp(lb + lna - lP);
	const double hab = std::exp(la + ldb - pp.lkap - lP);
	const double paa = -a * ga - pp.rkap * std::exp(la + ldb - lP);
	const double pbb = -b * gb - pp.rkap * std::exp(lb + lda - lP);
	const double nn = a * a - 2.0 * r * a * b + b * b;
	const double dl2 = (r * omr + a * b * omr - r * nn) / (omr * omr);
	o.L = lP;
	o.La = ga;
	o.Lb = gb;
	o.Laa = paa - ga * ga;
	o.Lbb = pbb - gb * gb;
	o.Lab = hab - ga * gb;
	o.Lr = hab;
	o.Lrr = hab * dl2 - hab * hab;
}

inline void pair_logp(double a, double b, double r, PairTerms& o) {
	PairPlan pp;
	pair_make_plan(r, pp);
	pair_logp_pre(a, b, pp, pp.bv.kind == 0 ? phid_ref(a) : 0.0, o);
}

}  // namespace lame_ae

#endif
