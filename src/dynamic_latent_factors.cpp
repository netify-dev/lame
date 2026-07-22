// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// SPD inverse with an escalating jitter ladder. The AR(1) prior terms supply
// proper prior precision at every t (stationary init at t = 0), but with rho
// near 1 the t = 0 prior precision can be tiny, so a near-singular V'V needs
// a numerical fallback. The jitter is a conditioning safeguard only (1e-8
// scale), not a prior.
static bool safe_inv_sympd_uv(const arma::mat& M, arma::mat& out) {
  arma::mat A = 0.5 * (M + M.t());
  double dscale = std::max(1e-12, A.diag().max());
  double jitters[4] = {0.0, 1e-8, 1e-6, 1e-4};
  for(int k = 0; k < 4; k++) {
    arma::mat A2 = A;
    if(jitters[k] > 0.0) A2.diag() += jitters[k] * dscale;
    bool ok = arma::inv_sympd(out, A2);
    if(ok && out.is_finite()) return true;
  }
  return false;
}

// Gauge-fixing rotation applied identically to every time slice. A single
// COMMON orthogonal rotation of (U_t, V_t) for all t is a symmetry of both
// the likelihood (U_t V_t' invariant) and the isotropic AR(1) prior
// (innovations u_t - rho u_{t-1} are rotated jointly, preserving their law),
// so it leaves every rotation-invariant functional untouched. Per-period
// rotations -- the previous behaviour -- are NOT a symmetry (they change the
// innovations) and were removed.
static void common_rotation_uv(arma::cube& U, arma::cube& V) {
  const int n = U.n_rows;
  const int R = U.n_cols;
  const int T = U.n_slices;
  if(R <= 0 || T <= 0) return;

  arma::mat avg_prod(n, V.n_rows, fill::zeros);
  arma::mat U_avg(n, R, fill::zeros);
  for(int t = 0; t < T; t++) {
    avg_prod += U.slice(t) * V.slice(t).t();
    U_avg += U.slice(t);
  }
  avg_prod /= T;
  U_avg /= T;

  arma::vec s_vals;
  arma::mat U_svd, V_svd;
  if(!svd_econ(U_svd, s_vals, V_svd, avg_prod)) return;
  const int rank = std::min(R, (int)s_vals.n_elem);
  if(rank <= 0) return;

  arma::mat U_ref(n, R, fill::zeros);
  U_ref.head_cols(rank) = U_svd.head_cols(rank) *
    diagmat(sqrt(s_vals.head(rank)));

  // one Procrustes rotation of the time-averaged U onto the reference
  arma::mat M = U_ref.t() * U_avg;
  arma::mat Mu, Mv;
  arma::vec Ms;
  if(!svd_econ(Mu, Ms, Mv, M)) return;
  arma::mat rot = Mv * Mu.t();
  if(!rot.is_finite()) return;

  for(int t = 0; t < T; t++) {
    U.slice(t) = U.slice(t) * rot;
    V.slice(t) = V.slice(t) * rot;
  }
}

//' Update dynamic latent positions using AR(1) process
//'
//' @param U_current Current 3D array of U positions (n x R x T)
//' @param V_current Current 3D array of V positions (n x R x T)
//' @param ET 3D array of residuals (n x n x T)
//' @param rho_uv AR(1) autoregressive parameter for latent positions
//' @param sigma_uv Innovation standard deviation for latent positions
//' @param s2 Dyadic variance
//' @param shrink Whether to apply shrinkage
//' @param symmetric Whether network is symmetric
//' @return List with updated U and V arrays
// [[Rcpp::export]]
List rUV_dynamic_fc_cpp(arma::cube U_current, arma::cube V_current,
                        const arma::cube& ET, double rho_uv,
                        double sigma_uv, double s2, bool shrink,
                        bool symmetric) {

  const int n = U_current.n_rows;
  const int R = U_current.n_cols;
  const int T = U_current.n_slices;

  // ar(1) prior terms need s2 scaling to match the s2-parameterized
  // posterior: VtV_mod = V'V + s2*(ar terms), var = s2 * inv(VtV_mod)
  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2 * sigma2_inv;
  const double s2_rho_sigma2 = s2 * rho_uv * sigma2_inv;
  const double s2_rho2_sigma2 = s2 * rho_uv * rho_uv * sigma2_inv;
  // stationary AR(1) initial condition at t = 0: u_1 ~ N(0, sigma^2/(1-rho^2)),
  // i.e. s2-scaled prior precision s2 * (1 - rho^2) / sigma^2
  const double s2_stat_init = s2 * (1.0 - rho_uv * rho_uv) * sigma2_inv;

  arma::cube& U_new = U_current;
  arma::cube& V_new = V_current;

  arma::mat VtV(R, R);
  arma::mat UtU(R, R);
  arma::mat cholDecomp(R, R);

  for(int t = 0; t < T; t++) {
    const arma::mat& E_t = ET.slice(t);
    arma::mat& U_t = U_new.slice(t);
    arma::mat& V_t = V_new.slice(t);

    VtV = V_t.t() * V_t;

    // batch compute all ei = V_t^T * E_t^T as a single BLAS call.
    // E_all_U = E_t * V_t gives n x R, each row i = E_t.row(i) * V_t
    arma::mat E_all_U = E_t * V_t;  // n x R (one BLAS gemm instead of n gemv calls)

    for(int i = 0; i < n; i++) {
      arma::vec ei = E_all_U.row(i).t();  // Extract row (O(R) not O(nR))
      arma::mat VtV_mod = VtV;

      if(t > 0) {
        VtV_mod.diag() += s2_sigma2_inv;
        ei += s2_rho_sigma2 * U_current.slice(t-1).row(i).t();
      } else {
        // stationary AR(1) prior at t = 0 (the model's marginal init)
        VtV_mod.diag() += s2_stat_init;
      }
      if(t < T-1) {
        VtV_mod.diag() += s2_rho2_sigma2;
        ei += s2_rho_sigma2 * U_current.slice(t+1).row(i).t();
      }

      arma::mat iVtV;
      if(!safe_inv_sympd_uv(VtV_mod, iVtV)) continue;
      arma::vec mu = iVtV * ei;

      cholDecomp = chol(s2 * iVtV, "lower");
      arma::vec randVec(R, fill::randn);
      U_t.row(i) = (mu + cholDecomp * randVec).t();
    }

    if(!symmetric) {
      UtU = U_t.t() * U_t;

      // batch compute all ej = U_t^T * E_t as a single BLAS call
      arma::mat E_all_V = E_t.t() * U_t;  // n x R

      for(int j = 0; j < n; j++) {
        arma::vec ej = E_all_V.row(j).t();  // Extract row
        arma::mat UtU_mod = UtU;

        if(t > 0) {
          UtU_mod.diag() += s2_sigma2_inv;
          ej += s2_rho_sigma2 * V_current.slice(t-1).row(j).t();
        } else {
          UtU_mod.diag() += s2_stat_init;
        }
        if(t < T-1) {
          UtU_mod.diag() += s2_rho2_sigma2;
          ej += s2_rho_sigma2 * V_current.slice(t+1).row(j).t();
        }

        arma::mat iUtU;
        if(!safe_inv_sympd_uv(UtU_mod, iUtU)) continue;
        arma::vec mu = iUtU * ej;

        cholDecomp = chol(s2 * iUtU, "lower");
        arma::vec randVec(R, fill::randn);
        V_t.row(j) = (mu + cholDecomp * randVec).t();
      }
    } else {
      V_new.slice(t) = U_new.slice(t);
    }
  }

  // gauge fixing: one common rotation across all time slices (posterior
  // symmetry); the previous per-period Procrustes and U/V column scale
  // balancing were state-dependent maps that broke the AR(1) prior's law
  // and were removed from the kernel.
  if(shrink) {
    common_rotation_uv(U_new, V_new);
  }

  return List::create(Named("U") = U_new, Named("V") = V_new);
}

//' Initialize dynamic latent positions with AR(1) structure
//'
//' @param n Number of actors
//' @param R Latent dimension
//' @param Tn Number of time points
//' @param rho_uv AR(1) parameter
//' @param sigma_uv Innovation standard deviation
//' @return 3D array of latent positions (n x R x Tn)
// [[Rcpp::export]]
arma::cube init_dynamic_positions(int n, int R, int Tn,
                                  double rho_uv, double sigma_uv) {
  arma::cube positions(n, R, Tn);

  positions.slice(0) = randn(n, R);

  for(int t = 1; t < Tn; t++) {
    for(int i = 0; i < n; i++) {
      for(int r = 0; r < R; r++) {
        double innov = R::rnorm(0.0, sigma_uv);
        positions(i, r, t) = rho_uv * positions(i, r, t-1) + innov;
      }
    }
  }

  return positions;
}

//' Sample AR(1) parameter for dynamic latent factors
//'
//' Independence Metropolis-Hastings update coherent with the stationary
//' AR(1) initial condition used by the dynamic UV kernels. The proposal is
//' the conjugate Normal implied by the transitions and the
//' Normal(prior_mean, prior_sd^2) prior, drawn truncated to (-0.99, 0.99)
//' via the inverse CDF; the acceptance ratio supplies the t = 1 stationary
//' factor proportional to \code{(1 - rho^2)^(nR / 2) * exp(-(1 - rho^2) *
//' S1 / (2 * sigma^2))}
//' that the Gaussian proposal omits. Defaults (prior_mean = 0,
//' prior_sd = 1) preserve the historical behaviour; lame::lame() passes
//' the user-set prior$rho_uv_mean / prior$rho_uv_sd explicitly.
//'
//' @param U_cube 3D array of U positions (n x R x T)
//' @param V_cube 3D array of V positions (n x R x T)
//' @param sigma_uv Innovation standard deviation
//' @param rho_current Current value of rho
//' @param symmetric Whether network is symmetric
//' @param prior_mean Prior mean for rho (default 0)
//' @param prior_sd Prior SD for rho (default 1)
//' @return Updated rho value
// [[Rcpp::export]]
double sample_rho_uv(const arma::cube& U_cube, const arma::cube& V_cube,
                     double sigma_uv, double rho_current, bool symmetric,
                     double prior_mean = 0.0, double prior_sd = 1.0) {

  const int T = U_cube.n_slices;

  double sum_yt_yt1 = 0.0;
  double sum_yt1_yt1 = 0.0;

  for(int t = 1; t < T; t++) {
    const arma::mat& U_t = U_cube.slice(t);
    const arma::mat& U_t1 = U_cube.slice(t-1);

    sum_yt_yt1 += accu(U_t % U_t1);
    sum_yt1_yt1 += accu(U_t1 % U_t1);

    if(!symmetric) {
      const arma::mat& V_t = V_cube.slice(t);
      const arma::mat& V_t1 = V_cube.slice(t-1);

      sum_yt_yt1 += accu(V_t % V_t1);
      sum_yt1_yt1 += accu(V_t1 % V_t1);
    }
  }

  // stationary-init sufficient statistics (t = 1 slice)
  double S1 = accu(U_cube.slice(0) % U_cube.slice(0));
  double n1 = (double)(U_cube.n_rows * U_cube.n_cols);
  if(!symmetric) {
    S1 += accu(V_cube.slice(0) % V_cube.slice(0));
    n1 += (double)(V_cube.n_rows * V_cube.n_cols);
  }

  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double prior_prec = 1.0 / (prior_sd * prior_sd);
  const double var_post  = 1.0 / (sum_yt1_yt1 * sigma2_inv + prior_prec);
  const double mean_post = var_post * (sum_yt_yt1 * sigma2_inv +
                                       prior_mean * prior_prec);

  // truncated-normal independence proposal on (-0.99, 0.99) via inverse CDF
  const double lo = -0.99, hi = 0.99;
  const double sd_post = sqrt(var_post);
  double F_lo = R::pnorm(lo, mean_post, sd_post, 1, 0);
  double F_hi = R::pnorm(hi, mean_post, sd_post, 1, 0);
  if(!(F_hi > F_lo)) return rho_current;
  double u = R::runif(0.0, 1.0) * (F_hi - F_lo) + F_lo;
  u = std::max(1e-12, std::min(1.0 - 1e-12, u));
  double rho_prop = R::qnorm(u, mean_post, sd_post, 1, 0);
  rho_prop = std::max(lo + 1e-12, std::min(hi - 1e-12, rho_prop));

  // MH correction with the stationary t = 1 factor (proposal density and
  // truncation constants cancel; the ratio reduces to f1(prop) / f1(curr))
  double rc = std::max(lo, std::min(hi, rho_current));
  double log_f1_prop = 0.5 * n1 * log(1.0 - rho_prop * rho_prop) -
    0.5 * (1.0 - rho_prop * rho_prop) * S1 * sigma2_inv;
  double log_f1_curr = 0.5 * n1 * log(1.0 - rc * rc) -
    0.5 * (1.0 - rc * rc) * S1 * sigma2_inv;

  if(log(R::runif(0.0, 1.0)) < log_f1_prop - log_f1_curr) {
    return rho_prop;
  }
  return rho_current;
}

//' Sample innovation variance for dynamic latent factors
//'
//' @param U_cube 3D array of U positions (n x R x T)
//' @param V_cube 3D array of V positions (n x R x T)
//' @param rho_uv AR(1) parameter
//' @param symmetric Whether network is symmetric
//' @return Updated sigma_uv value
// [[Rcpp::export]]
double sample_sigma_uv(const arma::cube& U_cube, const arma::cube& V_cube,
                       double rho_uv, bool symmetric) {

  const int n = U_cube.n_rows;
  const int R = U_cube.n_cols;
  const int T = U_cube.n_slices;

  double ss = 0.0;

  // stationary-init contribution (t = 1 slice): coherent with the kernels'
  // stationary AR(1) initial condition u_1 ~ N(0, sigma^2 / (1 - rho^2))
  const double one_m_rho2 = 1.0 - std::min(rho_uv * rho_uv, 0.9801);
  ss += one_m_rho2 * accu(U_cube.slice(0) % U_cube.slice(0));
  int count_init = n * R;
  if(!symmetric) {
    ss += one_m_rho2 * accu(V_cube.slice(0) % V_cube.slice(0));
    count_init += (int)(V_cube.n_rows * V_cube.n_cols);
  }

  for(int t = 1; t < T; t++) {
    const arma::mat& U_t = U_cube.slice(t);
    const arma::mat& U_t1 = U_cube.slice(t-1);

    arma::mat innov_U = U_t - rho_uv * U_t1;
    ss += accu(innov_U % innov_U);

    if(!symmetric) {
      const arma::mat& V_t = V_cube.slice(t);
      const arma::mat& V_t1 = V_cube.slice(t-1);

      arma::mat innov_V = V_t - rho_uv * V_t1;
      ss += accu(innov_V % innov_V);
    }
  }

  const int count = n * R * (T - 1) * (symmetric ? 1 : 2) + count_init;

  const double shape = count / 2.0 + 1.0;
  const double scale = ss / 2.0 + 1.0;
  const double sigma2_new = 1.0 / R::rgamma(shape, 1.0/scale);

  return sqrt(sigma2_new);
}

// gaussian log marginal of the data under an isotropic prior, in the
// s2-scaled parameterization used by the dynamic uv sampler. prior precision
// is c_prior * I; post precision B = V'V + c_prior * I; eta_post / eta_prior
// are the s2-scaled mean numerators. the data normalizing constant is common
// across competing priors and is dropped, so only differences are meaningful.
static double snap_log_marginal(const arma::mat& iB, double logdet_B,
                                const arma::vec& eta_post,
                                const arma::vec& eta_prior,
                                double c_prior, double s2, int R) {
  double q_post = arma::dot(eta_post, iB * eta_post);
  double q_prior = arma::dot(eta_prior, eta_prior) / c_prior;
  return 0.5 * (R * std::log(c_prior) - logdet_B + (q_post - q_prior) / s2);
}

// draw one latent vector, choosing per actor between an ar(1) drift prior and
// a diffuse snap prior. fills update_row with the new position and returns the
// snap indicator. c_stat is the s2-scaled stationary-init prior precision
// used at t = 0 (no backward neighbour).
static int snap_draw_row(arma::rowvec& update_row, const arma::mat& VtV,
                         const arma::vec& ei_data, const arma::vec& back_mean,
                         double s2_sigma2_inv, double s2_kappa_inv,
                         double fwd_diag, const arma::vec& fwd_ei,
                         double logit_pi, double s2, int R, bool has_back,
                         double c_stat) {
  if(!has_back) {
    // t = 0: stationary AR(1) init prior, no snap choice
    double c = fwd_diag + c_stat;
    arma::mat B = VtV; B.diag() += c;
    arma::mat iB;
    if(!safe_inv_sympd_uv(B, iB)) return 0;
    arma::vec mu = iB * (ei_data + fwd_ei);
    arma::mat L = chol(s2 * iB, "lower");
    arma::vec z(R, fill::randn);
    update_row = (mu + L * z).t();
    return 0;
  }

  double c0 = s2_sigma2_inv + fwd_diag;
  double c1 = s2_kappa_inv + fwd_diag;
  arma::vec eta_prior0 = back_mean + fwd_ei;
  arma::vec eta_prior1 = fwd_ei;
  arma::vec eta_post0 = ei_data + eta_prior0;
  arma::vec eta_post1 = ei_data + eta_prior1;

  arma::mat B0 = VtV; B0.diag() += c0;
  arma::mat B1 = VtV; B1.diag() += c1;
  arma::mat iB0, iB1;
  if(!safe_inv_sympd_uv(B0, iB0) || !safe_inv_sympd_uv(B1, iB1)) return 0;

  double ld0, ld1, sgn;
  arma::log_det(ld0, sgn, B0);
  arma::log_det(ld1, sgn, B1);

  double log_ml0 = snap_log_marginal(iB0, ld0, eta_post0, eta_prior0, c0, s2, R);
  double log_ml1 = snap_log_marginal(iB1, ld1, eta_post1, eta_prior1, c1, s2, R);

  double log_odds = logit_pi + log_ml1 - log_ml0;
  log_odds = std::max(-20.0, std::min(20.0, log_odds));
  double prob_snap = 1.0 / (1.0 + std::exp(-log_odds));
  int d = (arma::randu() < prob_snap) ? 1 : 0;

  const arma::mat& iB = d ? iB1 : iB0;
  const arma::vec& eta_post = d ? eta_post1 : eta_post0;
  arma::vec mu = iB * eta_post;
  arma::mat L = chol(s2 * iB, "lower");
  arma::vec z(R, fill::randn);
  update_row = (mu + L * z).t();
  return d;
}

//' Update dynamic latent positions with snap-shift model selection
//'
//' Like rUV_dynamic_fc_cpp but, for t > 0, chooses per actor between an AR(1)
//' drift prior and a diffuse N(0, kappa^2 I) snap prior via a Gaussian
//' log-marginal-likelihood model selection, drawing a Bernoulli snap indicator
//' delta and sampling the latent position from the selected posterior.
//'
//' @param U_current Current 3D array of U positions (n x R x T)
//' @param V_current Current 3D array of V positions (n x R x T)
//' @param ET 3D array of residuals (n x n x T)
//' @param rho_uv AR(1) autoregressive parameter for the drift prior
//' @param sigma_uv Innovation standard deviation for the drift prior
//' @param s2 Dyadic variance
//' @param kappa Diffuse snap-prior standard deviation (kappa^2 >> sigma_uv^2)
//' @param pi_snap Prior snap probability
//' @param delta_u_current Current sender-side snap indicators from the previous sweep
//' @param delta_v_current Current receiver-side snap indicators from the previous sweep
//' @param shrink Whether to apply shrinkage
//' @param symmetric Whether network is symmetric
//' @return List with updated U, V arrays and delta_u, delta_v snap indicators
// [[Rcpp::export]]
List rUV_dynamic_snap_fc_cpp(arma::cube U_current, arma::cube V_current,
	                             const arma::cube& ET, double rho_uv,
	                             double sigma_uv, double s2, double kappa,
	                             double pi_snap,
	                             const arma::mat& delta_u_current,
	                             const arma::mat& delta_v_current,
	                             bool shrink, bool symmetric) {

  const int n = U_current.n_rows;
  const int R = U_current.n_cols;
  const int T = U_current.n_slices;

  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2 * sigma2_inv;
  const double s2_rho_sigma2 = s2 * rho_uv * sigma2_inv;
  const double s2_rho2_sigma2 = s2 * rho_uv * rho_uv * sigma2_inv;
  const double s2_kappa_inv = s2 / (kappa * kappa);
  // stationary AR(1) init prior at t = 0 (s2-scaled precision)
  const double s2_stat_init = s2 * (1.0 - rho_uv * rho_uv) * sigma2_inv;
  const double logit_pi = std::log(pi_snap + 1e-300) -
    std::log(1.0 - pi_snap + 1e-300);

  arma::cube& U_new = U_current;
  arma::cube& V_new = V_current;
  arma::mat delta_u(n, T, fill::zeros);
  arma::mat delta_v(n, T, fill::zeros);

  for(int t = 0; t < T; t++) {
    const arma::mat& E_t = ET.slice(t);
    arma::mat& U_t = U_new.slice(t);
    arma::mat& V_t = V_new.slice(t);

    arma::mat VtV = V_t.t() * V_t;
    arma::mat E_all_U = E_t * V_t;

    for(int i = 0; i < n; i++) {
      arma::vec ei = E_all_U.row(i).t();

      double fwd_diag = 0.0;
      arma::vec fwd_ei(R, fill::zeros);
	      if(t < T-1 && delta_u_current(i, t+1) < 0.5) {
	        fwd_diag = s2_rho2_sigma2;
	        fwd_ei = s2_rho_sigma2 * U_current.slice(t+1).row(i).t();
	      }
      arma::vec back_mean(R, fill::zeros);
      if(t > 0) {
        back_mean = s2_rho_sigma2 * U_current.slice(t-1).row(i).t();
      }

      arma::rowvec new_row = U_t.row(i);
      int d = snap_draw_row(new_row, VtV, ei, back_mean,
                            s2_sigma2_inv, s2_kappa_inv, fwd_diag, fwd_ei,
                            logit_pi, s2, R, t > 0, s2_stat_init);
      U_t.row(i) = new_row;
      delta_u(i, t) = d;
    }

    if(!symmetric) {
      arma::mat UtU = U_t.t() * U_t;
      arma::mat E_all_V = E_t.t() * U_t;

      for(int j = 0; j < n; j++) {
        arma::vec ej = E_all_V.row(j).t();

        double fwd_diag = 0.0;
        arma::vec fwd_ej(R, fill::zeros);
	        if(t < T-1 && delta_v_current(j, t+1) < 0.5) {
	          fwd_diag = s2_rho2_sigma2;
	          fwd_ej = s2_rho_sigma2 * V_current.slice(t+1).row(j).t();
	        }
        arma::vec back_mean(R, fill::zeros);
        if(t > 0) {
          back_mean = s2_rho_sigma2 * V_current.slice(t-1).row(j).t();
        }

        arma::rowvec new_row = V_t.row(j);
        int d = snap_draw_row(new_row, UtU, ej, back_mean,
                              s2_sigma2_inv, s2_kappa_inv, fwd_diag, fwd_ej,
                              logit_pi, s2, R, t > 0, s2_stat_init);
        V_t.row(j) = new_row;
        delta_v(j, t) = d;
      }
    } else {
      V_new.slice(t) = U_new.slice(t);
      delta_v.col(t) = delta_u.col(t);
    }
  }

  // gauge fixing: one common rotation across all time slices (see
  // rUV_dynamic_fc_cpp); per-period rotation and scale balancing removed.
  if(shrink) {
    common_rotation_uv(U_new, V_new);
  }

  return List::create(Named("U") = U_new, Named("V") = V_new,
                      Named("delta_u") = delta_u, Named("delta_v") = delta_v);
}

//' Update dynamic latent positions with heavy-tailed (Student-t) AR(1) innovations
//'
//' Like rUV_dynamic_fc_cpp but each AR(1) innovation is Student-t rather than
//' Gaussian, via a scale-mixture: the innovation for \code{u[t, i]} has variance
//' \code{sigma^2 / lambda[t, i]} with \code{lambda[t, i]} distributed as
//' Gamma(nu/2, nu/2). Provides a
//' continuous heavy-tailed alternative to the discrete snap-shift model.
//'
//' @param U_current Current 3D array of U positions (n x R x T)
//' @param V_current Current 3D array of V positions (n x R x T)
//' @param ET 3D array of residuals (n x n x T)
//' @param rho_uv AR(1) autoregressive parameter
//' @param sigma_uv Innovation scale
//' @param s2 Dyadic variance
//' @param nu Student-t degrees of freedom
//' @param lambda_u Current local scales for U (n x T)
//' @param lambda_v Current local scales for V (n x T)
//' @param shrink Whether to apply shrinkage
//' @param symmetric Whether network is symmetric
//' @return List with updated U, V arrays and lambda_u, lambda_v local scales
// [[Rcpp::export]]
List rUV_dynamic_t_fc_cpp(arma::cube U_current, arma::cube V_current,
                          const arma::cube& ET, double rho_uv,
                          double sigma_uv, double s2, double nu,
                          arma::mat lambda_u, arma::mat lambda_v,
                          bool shrink, bool symmetric) {

  const int n = U_current.n_rows;
  const int R = U_current.n_cols;
  const int T = U_current.n_slices;

  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2 * sigma2_inv;
  const double s2_rho_sigma2 = s2 * rho_uv * sigma2_inv;
  const double s2_rho2_sigma2 = s2 * rho_uv * rho_uv * sigma2_inv;
  // t = 0 initial condition: exactly the Gaussian stationary AR(1) marginal
  // N(0, sigma^2 / (1 - rho^2)), so the transition-based hyper updates stay
  // coherent by construction (the t mixture applies to transitions only)
  const double s2_stat_init = s2 * (1.0 - rho_uv * rho_uv) * sigma2_inv;

  arma::cube& U_new = U_current;
  arma::cube& V_new = V_current;

  for(int t = 0; t < T; t++) {
    const arma::mat& E_t = ET.slice(t);
    arma::mat& U_t = U_new.slice(t);
    arma::mat& V_t = V_new.slice(t);

    arma::mat VtV = V_t.t() * V_t;
    arma::mat E_all_U = E_t * V_t;

    for(int i = 0; i < n; i++) {
      arma::vec ei = E_all_U.row(i).t();
      arma::mat VtV_mod = VtV;

      // backward term scaled by this period's local precision
      if(t > 0) {
        double lam = lambda_u(i, t);
        VtV_mod.diag() += lam * s2_sigma2_inv;
        ei += lam * s2_rho_sigma2 * U_current.slice(t-1).row(i).t();
      } else {
        VtV_mod.diag() += s2_stat_init;
      }
      // forward term scaled by the next period's local precision
      if(t < T-1) {
        double lam_next = lambda_u(i, t+1);
        VtV_mod.diag() += lam_next * s2_rho2_sigma2;
        ei += lam_next * s2_rho_sigma2 * U_current.slice(t+1).row(i).t();
      }

      arma::mat iVtV;
      if(!safe_inv_sympd_uv(VtV_mod, iVtV)) continue;
      arma::vec mu = iVtV * ei;
      arma::mat L = chol(s2 * iVtV, "lower");
      arma::vec z(R, fill::randn);
      U_t.row(i) = (mu + L * z).t();
    }

    if(!symmetric) {
      arma::mat UtU = U_t.t() * U_t;
      arma::mat E_all_V = E_t.t() * U_t;

      for(int j = 0; j < n; j++) {
        arma::vec ej = E_all_V.row(j).t();
        arma::mat UtU_mod = UtU;

        if(t > 0) {
          double lam = lambda_v(j, t);
          UtU_mod.diag() += lam * s2_sigma2_inv;
          ej += lam * s2_rho_sigma2 * V_current.slice(t-1).row(j).t();
        } else {
          UtU_mod.diag() += s2_stat_init;
        }
        if(t < T-1) {
          double lam_next = lambda_v(j, t+1);
          UtU_mod.diag() += lam_next * s2_rho2_sigma2;
          ej += lam_next * s2_rho_sigma2 * V_current.slice(t+1).row(j).t();
        }

        arma::mat iUtU;
        if(!safe_inv_sympd_uv(UtU_mod, iUtU)) continue;
        arma::vec mu = iUtU * ej;
        arma::mat L = chol(s2 * iUtU, "lower");
        arma::vec z(R, fill::randn);
        V_t.row(j) = (mu + L * z).t();
      }
    } else {
      V_new.slice(t) = U_new.slice(t);
    }
  }

  // sample the local scales from the ar(1) residuals (t scale-mixture)
  arma::mat lambda_u_new = lambda_u;
  arma::mat lambda_v_new = lambda_v;
  double shape = (nu + R) / 2.0;
  for(int t = 1; t < T; t++) {
    for(int i = 0; i < n; i++) {
      arma::rowvec eps = U_new.slice(t).row(i) - rho_uv * U_new.slice(t-1).row(i);
      double rate = (nu + arma::dot(eps, eps) * sigma2_inv) / 2.0;
      lambda_u_new(i, t) = R::rgamma(shape, 1.0 / rate);
    }
    if(!symmetric) {
      for(int j = 0; j < n; j++) {
        arma::rowvec eps = V_new.slice(t).row(j) - rho_uv * V_new.slice(t-1).row(j);
        double rate = (nu + arma::dot(eps, eps) * sigma2_inv) / 2.0;
        lambda_v_new(j, t) = R::rgamma(shape, 1.0 / rate);
      }
    } else {
      lambda_v_new.col(t) = lambda_u_new.col(t);
    }
  }

  // gauge fixing: one common rotation across all time slices (see
  // rUV_dynamic_fc_cpp); per-period rotation and scale balancing removed.
  if(shrink) {
    common_rotation_uv(U_new, V_new);
  }

  return List::create(Named("U") = U_new, Named("V") = V_new,
                      Named("lambda_u") = lambda_u_new, Named("lambda_v") = lambda_v_new);
}
