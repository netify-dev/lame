#ifndef LAME_AR_TRUNCATED_NORMAL_H
#define LAME_AR_TRUNCATED_NORMAL_H

// Tail-stable finite-interval AR draw. Central probabilities retain the
// historical inverse-CDF arithmetic and consume the same single uniform.
inline double lame_ar_truncnorm_draw(double mean, double sd,
                                    double lower, double upper) {
  if (!std::isfinite(mean) || !std::isfinite(sd) || sd <= 0.0 ||
      !std::isfinite(lower) || !std::isfinite(upper) || lower >= upper)
    Rcpp::stop("Invalid truncated-normal AR proposal parameters.");
  const double flo = R::pnorm(lower, mean, sd, 1, 0);
  const double fhi = R::pnorm(upper, mean, sd, 1, 0);
  for (int attempt = 0; attempt < 1000; ++attempt) {
    const double u = R::runif(0.0, 1.0);
    const double p = flo + u * (fhi - flo);
    double value;
    if (fhi > flo && p >= 1e-12 && p <= 1.0 - 1e-12) {
      value = R::qnorm(p, mean, sd, 1, 0);
    } else if (upper <= mean) {
      const double la = R::pnorm(lower, mean, sd, 1, 1);
      const double lb = R::pnorm(upper, mean, sd, 1, 1);
      const double logp = lb + std::log(u + (1.0 - u) * std::exp(la - lb));
      value = R::qnorm(logp, mean, sd, 1, 1);
    } else {
      const double la = R::pnorm(lower, mean, sd, 0, 1);
      const double lb = R::pnorm(upper, mean, sd, 0, 1);
      const double logp = la + std::log((1.0 - u) + u * std::exp(lb - la));
      value = R::qnorm(logp, mean, sd, 0, 1);
    }
    if (std::isfinite(value) && value > lower && value < upper) return value;
  }
  Rcpp::stop("Truncated-normal AR interval cannot be resolved numerically.");
  return NA_REAL;
}
#endif
