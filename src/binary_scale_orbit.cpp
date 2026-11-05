// Kernel of the scale move of R/binary_scale_orbit.R: the log likelihood of the observed cells of a binary probit panel when the part of the
// linear predictor that the scale carries is multiplied by k. `es` is that part at the observed cells, `mm` the rest of the predictor plus the shared
// dyad factor, `sg` the sign of the cell (+1 event, -1 non-event) and `inv_sd` the reciprocal of the conditional sd of a cell. One pass over the
// cells, no temporaries.
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double bso_loglik_cpp(const NumericVector& es, const NumericVector& mm, const NumericVector& sg, double k, double inv_sd) {
	const R_xlen_t n = es.size();
	if (mm.size() != n || sg.size() != n) Rcpp::stop("bso_loglik_cpp: the predictor, the offset and the signs must have the same length");
	double s = 0.0;
	// a cell classified beyond 9 standard deviations contributes less than 1e-19 to the log likelihood: it is left out of the sum (the cells of a separated
	// panel are mostly of that kind). The log of the normal distribution function is taken from the complementary error function, accurate to a few
	// units of the last place, wherever that does not underflow (z > -30), and from R's pnorm in the far tail
	const double rt2 = 0.70710678118654752440;
	for (R_xlen_t i = 0; i < n; ++i) {
		const double z = sg[i] * (k * es[i] + mm[i]) * inv_sd;
		if (z > 9.0) continue;
		if (z > 0.0) s += std::log1p(-0.5 * std::erfc(z * rt2));
		else if (z > -30.0) s += std::log(0.5 * std::erfc(-z * rt2));
		else s += R::pnorm(z, 0.0, 1.0, 1, 1);
	}
	return s;
}
