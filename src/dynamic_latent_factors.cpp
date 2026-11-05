// [[Rcpp::depends(RcppArmadillo)]]
// state arrays are taken as const references and copied inside each kernel:
// RcppArmadillo wraps the caller's memory for an arma::cube argument but copies
// an arma::mat, so a by-value cube parameter would write into the R object and
// defeat the revert on a refused draw

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// a precision whose diagonally scaled inverse reaches 0.01 / eps leaves under
// two digits to resolve it, and whether its inversion succeeds depends on rounding
static bool resolved_inverse_uv(const arma::mat& A, const arma::mat& inverse) {
  return inverse.is_finite() &&
    arma::max(inverse.diag() % A.diag()) * 100.0 * arma::datum::eps < 1.0;
}

// Preserve the supplied Gaussian precision. A diagonally scaled solve can
// improve its conditioning; adding a ridge would instead change the target.
static bool safe_inv_sympd_uv(const arma::mat& M, arma::mat& out) {
  arma::mat A = 0.5 * (M + M.t());
  if (!A.is_finite()) return false;
  if (arma::inv_sympd(out, A) && out.is_finite()) return resolved_inverse_uv(A, out);
  arma::vec diagonal = A.diag();
  if (arma::any(diagonal <= 0.0)) return false;
  arma::vec scale = 1.0 / arma::sqrt(diagonal);
  arma::mat normalized = A;
  normalized.each_col() %= scale;
  normalized.each_row() %= scale.t();
  arma::mat inverse;
  if (!arma::inv_sympd(inverse, normalized) || !inverse.is_finite()) return false;
  inverse.each_col() %= scale;
  inverse.each_row() %= scale.t();
  out = inverse;
  return resolved_inverse_uv(A, out);
}

// Gauge-fixing rotation applied identically to every time slice. A single
// COMMON orthogonal rotation of (U_t, V_t) for all t is a symmetry of both
// the likelihood (U_t V_t' invariant) and the isotropic AR(1) prior
// (innovations u_t - rho u_{t-1} are rotated jointly, preserving their law),
// so it leaves every rotation-invariant functional untouched. Per-period
// rotations are NOT a symmetry here (they change the innovations), so no
// per-period rotation is applied.
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
//' The per-period observation model is the pair model the rest of the
//' sampler uses: the two cells of a dyad, (e_ij, e_ji), have covariance
//' s2 * matrix(c(1, rho, rho, 1), 2, 2). Conditioning one actor's
//' coordinates on every other actor's leaves an exact Gaussian full
//' conditional in which the dyadic correlation enters three ways: the
//' likelihood precision is scaled by 1 / (1 - rho^2), the linear term
//' carries -rho times the reciprocal residual, and the self-cell (i, i)
//' contributes a rank-one correction. Passing \code{rho = 0} reproduces
//' the independent-dyad kernel exactly.
//'
//' @param U_current Current 3D array of U positions (n x R x T)
//' @param V_current Current 3D array of V positions (n x R x T)
//' @param ET 3D array of residuals (n x n x T)
//' @param rho_uv AR(1) autoregressive parameter for latent positions
//' @param sigma_uv Innovation standard deviation for latent positions
//' @param s2 Dyadic variance
//' @param max_abs Optional positive coordinate bound; nonpositive or infinite
//'   values disable it. Each actor conditional respects the bound.
//' @param shrink Whether to apply shrinkage
//' @param symmetric Whether network is symmetric
//' @param rho Dyadic correlation of the (i, j) and (j, i) residual pair
//' @return List with updated U and V arrays
// [[Rcpp::export]]
List rUV_dynamic_fc_cpp(const arma::cube& U_current, const arma::cube& V_current,
                        const arma::cube& ET, double rho_uv,
                        double sigma_uv, double s2, bool shrink,
                        bool symmetric, double rho = 0.0, double max_abs = -1.0) {

  const int n = U_current.n_rows;
  const int R = U_current.n_cols;
  const int T = U_current.n_slices;

  // the dyadic correlation whitens the pair, so the likelihood variance of
  // a single cell given its reciprocal is s2 * (1 - rho^2). a symmetric fit
  // mirrors V into U and carries its own gauge, so it stays uncorrelated.
  if (!symmetric && (!std::isfinite(rho) || std::abs(rho) >= 1.0))
    Rcpp::stop("Directed dyadic correlation must lie strictly between -1 and 1.");
  const bool use_rho = (!symmetric) && rho != 0.0;
  const double s2t = use_rho ? s2 * ((1.0 - rho) * (1.0 + rho)) : s2;
  if (!std::isfinite(s2t) || s2t <= 0.0)
    Rcpp::stop("Conditional residual variance must be finite and positive.");

  // ar(1) prior terms need s2 scaling to match the s2-parameterized
  // posterior: VtV_mod = V'V + s2*(ar terms), var = s2 * inv(VtV_mod)
  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2t * sigma2_inv;
  const double s2_rho_sigma2 = s2t * rho_uv * sigma2_inv;
  const double s2_rho2_sigma2 = s2t * rho_uv * rho_uv * sigma2_inv;
  // stationary AR(1) initial condition at t = 0: u_1 ~ N(0, sigma^2/(1-rho^2)),
  // i.e. s2-scaled prior precision s2 * (1 - rho^2) / sigma^2
  const double s2_stat_init = s2t * (1.0 - rho_uv * rho_uv) * sigma2_inv;

  // rcpparmadillo wraps the caller's array instead of copying it, so the
  // kernel draws into copies of its own: a sweep the caller refuses has to
  // leave the positions it was given exactly as they were
  arma::cube U_new = U_current;
  arma::cube V_new = V_current;
  const bool bounded = std::isfinite(max_abs) && max_abs > 0.0;
  int bound_rejections = 0;

  arma::mat VtV(R, R);
  arma::mat UtU(R, R);
  arma::mat cholDecomp(R, R);

  // a symmetric fit hands over V = U L with L the signed diagonal of the
  // product U L U'; each direction's l is read back from the pair
  arma::rowvec Lsym(R, fill::ones);
  if(symmetric) {
    for(int r = 0; r < R; r++) {
      double num = 0.0, den = 0.0;
      for(int t = 0; t < T; t++) {
        num += arma::dot(V_current.slice(t).col(r), U_current.slice(t).col(r));
        den += arma::dot(U_current.slice(t).col(r), U_current.slice(t).col(r));
      }
      if(den > 1e-300 && std::isfinite(num / den)) Lsym(r) = num / den;
      else {
        // Tiny squared coordinates can underflow although V = U * lambda
        // still identifies a finite signed weight. Recover it from the
        // largest nonzero coordinate, where the ratio is best resolved.
        double largest = 0.0, numerator = 0.0, denominator = 0.0;
        for(int t = 0; t < T; t++) {
          for(int i = 0; i < n; i++) {
            const double u = U_current(i, r, t), v = V_current(i, r, t);
            if(!std::isfinite(u) || !std::isfinite(v))
              Rcpp::stop("Symmetric latent coordinates must be finite.");
            if(std::abs(u) > largest) {
              largest = std::abs(u); numerator = v; denominator = u;
            }
          }
        }
        if(largest > 0.0) {
          const double weight = numerator / denominator;
          if(!std::isfinite(weight))
            Rcpp::stop("Unable to recover a finite symmetric signed weight.");
          Lsym(r) = weight;
        }
      }
    }
  }

  for(int t = 0; t < T; t++) {
    const arma::mat& E_t = ET.slice(t);
    arma::mat& U_t = U_new.slice(t);
    arma::mat& V_t = V_new.slice(t);

    if(symmetric) {
      // gibbs sweep of the symmetric product u_t l u_t': row i is drawn from
      // the regression of the off-diagonal cells of row i on the current
      // mirrored coordinates v_j = l u_j of every other actor, rows already
      // drawn in this period at their new values. the self cell is quartic
      // in u_i and is left out; the z step redraws the diagonal from its own
      // conditional given the new positions
      arma::mat Vs = U_t.each_row() % Lsym;
      arma::mat VtVs = Vs.t() * Vs;
      arma::mat Es = E_t * Vs;
      for(int i = 0; i < n; i++) {
        const arma::rowvec vi = Vs.row(i);
        const arma::rowvec u_old = U_t.row(i);
        arma::vec ei = Es.row(i).t() - E_t(i, i) * vi.t();
        arma::mat P = VtVs - vi.t() * vi;
        if(t > 0) {
          P.diag() += s2_sigma2_inv;
          ei += s2_rho_sigma2 * U_new.slice(t-1).row(i).t();
        } else {
          P.diag() += s2_stat_init;
        }
        if(t < T-1) {
          P.diag() += s2_rho2_sigma2;
          ei += s2_rho_sigma2 * U_new.slice(t+1).row(i).t();
        }
        arma::mat iP;
        if(!safe_inv_sympd_uv(P, iP)) Rcpp::stop("Unable to resolve the latent-factor Gaussian precision.");
        arma::vec mu = iP * ei;
        cholDecomp = chol(s2t * iP, "lower");
        arma::vec randVec(R, fill::randn);
        U_t.row(i) = (mu + cholDecomp * randVec).t();
      if (bounded && arma::abs(U_t.row(i)).max() > max_abs) {
        U_t.row(i) = u_old;
        ++bound_rejections;
      }
        const arma::rowvec vnew = U_t.row(i) % Lsym;
        VtVs += vnew.t() * vnew - vi.t() * vi;
        Es += E_t.col(i) * (vnew - vi);
        Vs.row(i) = vnew;
      }
      V_new.slice(t) = U_new.slice(t);
      continue;
    }

    VtV = V_t.t() * V_t;

    // full residual and the reciprocal-cell sums the correlated linear term
    // needs. both are carried forward by rank-one updates as each actor is
    // drawn, so every draw conditions on the coordinates already moved in
    // this sweep (a batched version would condition on stale ones and would
    // not be a gibbs sweep)
    arma::mat Resid;
    arma::mat Wrec;
    if(use_rho) {
      Resid = E_t - U_t * V_t.t();
      Wrec = Resid.t() * V_t;
    }

    // batch compute all ei = V_t^T * E_t^T as a single BLAS call.
    // E_all_U = E_t * V_t gives n x R, each row i = E_t.row(i) * V_t
    arma::mat E_all_U = E_t * V_t;  // n x R (one BLAS gemm instead of n gemv calls)

    for(int i = 0; i < n; i++) {
      arma::vec ei = E_all_U.row(i).t();  // Extract row (O(R) not O(nR))
      arma::mat VtV_mod = VtV;
      arma::rowvec u_old = U_t.row(i);

      if(use_rho) {
        u_old = U_t.row(i);
        const arma::rowvec vi = V_t.row(i);
        ei -= rho * (Wrec.row(i).t() + arma::dot(u_old, vi) * vi.t());
        VtV_mod -= rho * (vi.t() * vi);
      }

      if(t > 0) {
        VtV_mod.diag() += s2_sigma2_inv;
        ei += s2_rho_sigma2 * U_new.slice(t-1).row(i).t();
      } else {
        // stationary AR(1) prior at t = 0 (the model's marginal init)
        VtV_mod.diag() += s2_stat_init;
      }
      if(t < T-1) {
        VtV_mod.diag() += s2_rho2_sigma2;
        ei += s2_rho_sigma2 * U_new.slice(t+1).row(i).t();
      }

      arma::mat iVtV;
      if(!safe_inv_sympd_uv(VtV_mod, iVtV)) Rcpp::stop("Unable to resolve the latent-factor Gaussian precision.");
      arma::vec mu = iVtV * ei;

      cholDecomp = chol(s2t * iVtV, "lower");
      arma::vec randVec(R, fill::randn);
      U_t.row(i) = (mu + cholDecomp * randVec).t();
      if (bounded && arma::abs(U_t.row(i)).max() > max_abs) {
        U_t.row(i) = u_old;
        ++bound_rejections;
      }

      if(use_rho) {
        arma::vec Vdu = V_t * (U_t.row(i) - u_old).t();
        Resid.row(i) -= Vdu.t();
        Wrec -= Vdu * V_t.row(i);
      }
    }

    if(!symmetric) {
      UtU = U_t.t() * U_t;

      // batch compute all ej = U_t^T * E_t as a single BLAS call
      arma::mat E_all_V = E_t.t() * U_t;  // n x R

      arma::mat Mrec;
      if(use_rho) Mrec = Resid * U_t;

      for(int j = 0; j < n; j++) {
        arma::vec ej = E_all_V.row(j).t();  // Extract row
        arma::mat UtU_mod = UtU;
        arma::rowvec v_old = V_t.row(j);

        if(use_rho) {
          v_old = V_t.row(j);
          const arma::rowvec uj = U_t.row(j);
          ej -= rho * (Mrec.row(j).t() + arma::dot(v_old, uj) * uj.t());
          UtU_mod -= rho * (uj.t() * uj);
        }

        if(t > 0) {
          UtU_mod.diag() += s2_sigma2_inv;
          ej += s2_rho_sigma2 * V_new.slice(t-1).row(j).t();
        } else {
          UtU_mod.diag() += s2_stat_init;
        }
        if(t < T-1) {
          UtU_mod.diag() += s2_rho2_sigma2;
          ej += s2_rho_sigma2 * V_new.slice(t+1).row(j).t();
        }

        arma::mat iUtU;
        if(!safe_inv_sympd_uv(UtU_mod, iUtU)) Rcpp::stop("Unable to resolve the latent-factor Gaussian precision.");
        arma::vec mu = iUtU * ej;

        cholDecomp = chol(s2t * iUtU, "lower");
        arma::vec randVec(R, fill::randn);
        V_t.row(j) = (mu + cholDecomp * randVec).t();
      if (bounded && arma::abs(V_t.row(j)).max() > max_abs) {
        V_t.row(j) = v_old;
        ++bound_rejections;
      }

        if(use_rho) {
          arma::vec Udv = U_t * (V_t.row(j) - v_old).t();
          Resid.col(j) -= Udv;
          Mrec -= Udv * U_t.row(j);
        }
      }
    } else {
      V_new.slice(t) = U_new.slice(t);
    }
  }

  // gauge fixing: one common rotation across all time slices (posterior
  // symmetry). per-period Procrustes and U/V column scale balancing are
  // state-dependent maps that break the AR(1) prior's law, so the kernel
  // applies neither.
  // a symmetric fit's product is U L U' with a signed diagonal L, which a
  // common rotation of the coordinates does not preserve; the directed
  // product U V' is rotation invariant
  if(shrink && !symmetric && !bounded) {
    common_rotation_uv(U_new, V_new);
  }

  if (bounded) return List::create(Named("U") = U_new, Named("V") = V_new, Named("bound_rejections") = bound_rejections);
  return List::create(Named("U") = U_new, Named("V") = V_new);
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

// log density of the next position x_next under each branch's prior with the
// current position integrated out: the drift branch gives
// N(rho^2 x_prev, sigma^2 (1 + rho^2) I), the snap branch
// N(0, (sigma^2 + rho^2 kappa^2) I). the forward transition is a factor of
// both branch integrands, and folding it into a normalised prior drops this
// ratio, so the snap log odds carries it explicitly whenever the next period
// drifts from this one
static double snap_forward_log_ratio(const arma::vec& x_prev,
                                     const arma::vec& x_next, double rho_uv,
                                     double sigma2, double kappa2, int R) {
  const double v0 = sigma2 * (1.0 + rho_uv * rho_uv);
  const double v1 = sigma2 + rho_uv * rho_uv * kappa2;
  const arma::vec r0 = x_next - (rho_uv * rho_uv) * x_prev;
  return -0.5 * R * (std::log(v1) - std::log(v0)) -
    0.5 * arma::dot(x_next, x_next) / v1 + 0.5 * arma::dot(r0, r0) / v0;
}

// draw one latent vector, choosing per actor between an ar(1) drift prior and
// a diffuse snap prior. fills update_row with the new position and returns the
// snap indicator. c_stat is the s2-scaled stationary-init prior precision
// used at t = 0 (no backward neighbour). fwd_log_ratio is
// snap_forward_log_ratio() when the next period drifts from this one, else 0.
static int snap_draw_row(arma::rowvec& update_row, const arma::mat& VtV,
                         const arma::vec& ei_data, const arma::vec& back_mean,
                         double s2_sigma2_inv, double s2_kappa_inv,
                         double fwd_diag, const arma::vec& fwd_ei,
                         double logit_pi, double s2, int R, bool has_back,
                         double c_stat, double fwd_log_ratio) {
  if(!has_back) {
    // t = 0: stationary AR(1) init prior, no snap choice
    double c = fwd_diag + c_stat;
    arma::mat B = VtV; B.diag() += c;
    arma::mat iB;
    if(!safe_inv_sympd_uv(B, iB)) Rcpp::stop("Unable to resolve the latent-factor Gaussian precision.");
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
  if(!safe_inv_sympd_uv(B0, iB0) || !safe_inv_sympd_uv(B1, iB1)) Rcpp::stop("Unable to resolve the latent-factor Gaussian precision.");

  double ld0, ld1, sgn;
  arma::log_det(ld0, sgn, B0);
  arma::log_det(ld1, sgn, B1);

  double log_ml0 = snap_log_marginal(iB0, ld0, eta_post0, eta_prior0, c0, s2, R);
  double log_ml1 = snap_log_marginal(iB1, ld1, eta_post1, eta_prior1, c1, s2, R);

  double log_odds = logit_pi + log_ml1 - log_ml0 + fwd_log_ratio;
  if (std::isnan(log_odds)) Rcpp::stop("Invalid snap posterior odds.");
  double prob_snap;
  if (log_odds >= 0.0) prob_snap = 1.0 / (1.0 + std::exp(-log_odds));
  else { const double odds = std::exp(log_odds); prob_snap = odds / (1.0 + odds); }
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
//' @param max_abs Optional positive coordinate bound; nonpositive or infinite
//'   values disable it. Each actor conditional respects the bound.
//' @param shrink Whether to apply shrinkage
//' @param symmetric Whether network is symmetric
//' @param rho Dyadic correlation of the (i, j) and (j, i) residual pair
//' @return List with updated U, V arrays and delta_u, delta_v snap indicators
// [[Rcpp::export]]
List rUV_dynamic_snap_fc_cpp(const arma::cube& U_current, const arma::cube& V_current,
	                             const arma::cube& ET, double rho_uv,
	                             double sigma_uv, double s2, double kappa,
	                             double pi_snap,
	                             const arma::mat& delta_u_current,
	                             const arma::mat& delta_v_current,
	                             bool shrink, bool symmetric,
	                             double rho = 0.0, double max_abs = -1.0) {

  const int n = U_current.n_rows;
  const int R = U_current.n_cols;
  const int T = U_current.n_slices;

  // pair whitening, as in rUV_dynamic_fc_cpp
  if (!symmetric && (!std::isfinite(rho) || std::abs(rho) >= 1.0))
    Rcpp::stop("Directed dyadic correlation must lie strictly between -1 and 1.");
  const bool use_rho = (!symmetric) && rho != 0.0;
  const double s2t = use_rho ? s2 * ((1.0 - rho) * (1.0 + rho)) : s2;
  if (!std::isfinite(s2t) || s2t <= 0.0)
    Rcpp::stop("Conditional residual variance must be finite and positive.");

  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2t * sigma2_inv;
  const double s2_rho_sigma2 = s2t * rho_uv * sigma2_inv;
  const double s2_rho2_sigma2 = s2t * rho_uv * rho_uv * sigma2_inv;
  const double s2_kappa_inv = s2t / (kappa * kappa);
  // unscaled variances for the forward-transition ratio of the snap odds
  const double sigma2 = sigma_uv * sigma_uv;
  const double kappa2 = kappa * kappa;
  // stationary AR(1) init prior at t = 0 (s2-scaled precision)
  const double s2_stat_init = s2t * (1.0 - rho_uv * rho_uv) * sigma2_inv;
  if (!std::isfinite(pi_snap) || pi_snap < 0.0 || pi_snap > 1.0)
    Rcpp::stop("Snap probability must lie between zero and one.");
  const double logit_pi = std::log(pi_snap) - std::log1p(-pi_snap);

  // rcpparmadillo wraps the caller's array instead of copying it, so the
  // kernel draws into copies of its own: a sweep the caller refuses has to
  // leave the positions it was given exactly as they were
  arma::cube U_new = U_current;
  arma::cube V_new = V_current;
  const bool bounded = std::isfinite(max_abs) && max_abs > 0.0;
  int bound_rejections = 0;
  arma::mat delta_u(n, T, fill::zeros);
  arma::mat delta_v(n, T, fill::zeros);

  for(int t = 0; t < T; t++) {
    const arma::mat& E_t = ET.slice(t);
    arma::mat& U_t = U_new.slice(t);
    arma::mat& V_t = V_new.slice(t);

    arma::mat VtV = V_t.t() * V_t;
    arma::mat E_all_U = E_t * V_t;

    arma::mat Resid;
    arma::mat Wrec;
    if(use_rho) {
      Resid = E_t - U_t * V_t.t();
      Wrec = Resid.t() * V_t;
    }

    for(int i = 0; i < n; i++) {
      arma::vec ei = E_all_U.row(i).t();
      arma::mat VtV_i = VtV;
      arma::rowvec u_old = U_t.row(i);
      if(use_rho) {
        u_old = U_t.row(i);
        const arma::rowvec vi = V_t.row(i);
        ei -= rho * (Wrec.row(i).t() + arma::dot(u_old, vi) * vi.t());
        VtV_i -= rho * (vi.t() * vi);
      }

      double fwd_diag = 0.0;
      arma::vec fwd_ei(R, fill::zeros);
      const bool fwd_on = (t < T-1) && (delta_u_current(i, t+1) < 0.5);
      if(fwd_on) {
        fwd_diag = s2_rho2_sigma2;
        fwd_ei = s2_rho_sigma2 * U_new.slice(t+1).row(i).t();
      }
      arma::vec back_mean(R, fill::zeros);
      if(t > 0) {
        back_mean = s2_rho_sigma2 * U_new.slice(t-1).row(i).t();
      }
      double fwd_lr = 0.0;
      if(t > 0 && fwd_on) {
        fwd_lr = snap_forward_log_ratio(U_new.slice(t-1).row(i).t(),
                                        U_new.slice(t+1).row(i).t(),
                                        rho_uv, sigma2, kappa2, R);
      }

      arma::rowvec new_row = U_t.row(i);
      int d = snap_draw_row(new_row, VtV_i, ei, back_mean,
                            s2_sigma2_inv, s2_kappa_inv, fwd_diag, fwd_ei,
                            logit_pi, s2t, R, t > 0, s2_stat_init, fwd_lr);
      if (bounded && arma::abs(new_row).max() > max_abs) {
        new_row = u_old;
        d = static_cast<int>(delta_u_current(i, t));
        ++bound_rejections;
      }
      U_t.row(i) = new_row;
      delta_u(i, t) = d;

      if(use_rho) {
        arma::vec Vdu = V_t * (U_t.row(i) - u_old).t();
        Resid.row(i) -= Vdu.t();
        Wrec -= Vdu * V_t.row(i);
      }
    }

    if(!symmetric) {
      arma::mat UtU = U_t.t() * U_t;
      arma::mat E_all_V = E_t.t() * U_t;

      arma::mat Mrec;
      if(use_rho) Mrec = Resid * U_t;

      for(int j = 0; j < n; j++) {
        arma::vec ej = E_all_V.row(j).t();
        arma::mat UtU_j = UtU;
        arma::rowvec v_old = V_t.row(j);
        if(use_rho) {
          v_old = V_t.row(j);
          const arma::rowvec uj = U_t.row(j);
          ej -= rho * (Mrec.row(j).t() + arma::dot(v_old, uj) * uj.t());
          UtU_j -= rho * (uj.t() * uj);
        }

        double fwd_diag = 0.0;
        arma::vec fwd_ej(R, fill::zeros);
        const bool fwd_on = (t < T-1) && (delta_v_current(j, t+1) < 0.5);
        if(fwd_on) {
          fwd_diag = s2_rho2_sigma2;
          fwd_ej = s2_rho_sigma2 * V_new.slice(t+1).row(j).t();
        }
        arma::vec back_mean(R, fill::zeros);
        if(t > 0) {
          back_mean = s2_rho_sigma2 * V_new.slice(t-1).row(j).t();
        }
        double fwd_lr = 0.0;
        if(t > 0 && fwd_on) {
          fwd_lr = snap_forward_log_ratio(V_new.slice(t-1).row(j).t(),
                                          V_new.slice(t+1).row(j).t(),
                                          rho_uv, sigma2, kappa2, R);
        }

        arma::rowvec new_row = V_t.row(j);
        int d = snap_draw_row(new_row, UtU_j, ej, back_mean,
                              s2_sigma2_inv, s2_kappa_inv, fwd_diag, fwd_ej,
                              logit_pi, s2t, R, t > 0, s2_stat_init, fwd_lr);
        if (bounded && arma::abs(new_row).max() > max_abs) {
          new_row = v_old;
          d = static_cast<int>(delta_v_current(j, t));
          ++bound_rejections;
        }
        V_t.row(j) = new_row;
        delta_v(j, t) = d;

        if(use_rho) {
          arma::vec Udv = U_t * (V_t.row(j) - v_old).t();
          Resid.col(j) -= Udv;
          Mrec -= Udv * U_t.row(j);
        }
      }
    } else {
      V_new.slice(t) = U_new.slice(t);
      delta_v.col(t) = delta_u.col(t);
    }
  }

  // gauge fixing: one common rotation across all time slices (see
  // rUV_dynamic_fc_cpp); per-period rotation and scale balancing removed.
  // a symmetric fit's product is U L U' with a signed diagonal L, which a
  // common rotation of the coordinates does not preserve; the directed
  // product U V' is rotation invariant
  if(shrink && !symmetric && !bounded) {
    common_rotation_uv(U_new, V_new);
  }

  if (bounded) return List::create(Named("U") = U_new, Named("V") = V_new,
                      Named("delta_u") = delta_u, Named("delta_v") = delta_v, Named("bound_rejections") = bound_rejections);
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
//' @param max_abs Optional positive coordinate bound; nonpositive or infinite
//'   values disable it. Each actor conditional respects the bound.
//' @param shrink Whether to apply shrinkage
//' @param symmetric Whether network is symmetric
//' @param rho Dyadic correlation of the (i, j) and (j, i) residual pair
//' @return List with updated U, V arrays and lambda_u, lambda_v local scales
// [[Rcpp::export]]
List rUV_dynamic_t_fc_cpp(const arma::cube& U_current, const arma::cube& V_current,
                          const arma::cube& ET, double rho_uv,
                          double sigma_uv, double s2, double nu,
                          const arma::mat& lambda_u, const arma::mat& lambda_v,
                          bool shrink, bool symmetric, double rho = 0.0, double max_abs = -1.0) {

  const int n = U_current.n_rows;
  const int R = U_current.n_cols;
  const int T = U_current.n_slices;

  // pair whitening, as in rUV_dynamic_fc_cpp
  if (!symmetric && (!std::isfinite(rho) || std::abs(rho) >= 1.0))
    Rcpp::stop("Directed dyadic correlation must lie strictly between -1 and 1.");
  const bool use_rho = (!symmetric) && rho != 0.0;
  const double s2t = use_rho ? s2 * ((1.0 - rho) * (1.0 + rho)) : s2;
  if (!std::isfinite(s2t) || s2t <= 0.0)
    Rcpp::stop("Conditional residual variance must be finite and positive.");

  const double sigma2_inv = 1.0 / (sigma_uv * sigma_uv);
  const double s2_sigma2_inv = s2t * sigma2_inv;
  const double s2_rho_sigma2 = s2t * rho_uv * sigma2_inv;
  const double s2_rho2_sigma2 = s2t * rho_uv * rho_uv * sigma2_inv;
  // t = 0 initial condition: exactly the Gaussian stationary AR(1) marginal
  // N(0, sigma^2 / (1 - rho^2)), so the transition-based hyper updates stay
  // coherent by construction (the t mixture applies to transitions only)
  const double s2_stat_init = s2t * (1.0 - rho_uv * rho_uv) * sigma2_inv;

  // rcpparmadillo wraps the caller's array instead of copying it, so the
  // kernel draws into copies of its own: a sweep the caller refuses has to
  // leave the positions it was given exactly as they were
  arma::cube U_new = U_current;
  arma::cube V_new = V_current;
  const bool bounded = std::isfinite(max_abs) && max_abs > 0.0;
  int bound_rejections = 0;

  for(int t = 0; t < T; t++) {
    const arma::mat& E_t = ET.slice(t);
    arma::mat& U_t = U_new.slice(t);
    arma::mat& V_t = V_new.slice(t);

    arma::mat VtV = V_t.t() * V_t;
    arma::mat E_all_U = E_t * V_t;

    arma::mat Resid;
    arma::mat Wrec;
    if(use_rho) {
      Resid = E_t - U_t * V_t.t();
      Wrec = Resid.t() * V_t;
    }

    for(int i = 0; i < n; i++) {
      arma::vec ei = E_all_U.row(i).t();
      arma::mat VtV_mod = VtV;
      arma::rowvec u_old = U_t.row(i);

      if(use_rho) {
        u_old = U_t.row(i);
        const arma::rowvec vi = V_t.row(i);
        ei -= rho * (Wrec.row(i).t() + arma::dot(u_old, vi) * vi.t());
        VtV_mod -= rho * (vi.t() * vi);
      }

      // backward term scaled by this period's local precision
      if(t > 0) {
        double lam = lambda_u(i, t);
        VtV_mod.diag() += lam * s2_sigma2_inv;
        ei += lam * s2_rho_sigma2 * U_new.slice(t-1).row(i).t();
      } else {
        VtV_mod.diag() += s2_stat_init;
      }
      // forward term scaled by the next period's local precision
      if(t < T-1) {
        double lam_next = lambda_u(i, t+1);
        VtV_mod.diag() += lam_next * s2_rho2_sigma2;
        ei += lam_next * s2_rho_sigma2 * U_new.slice(t+1).row(i).t();
      }

      arma::mat iVtV;
      if(!safe_inv_sympd_uv(VtV_mod, iVtV)) Rcpp::stop("Unable to resolve the latent-factor Gaussian precision.");
      arma::vec mu = iVtV * ei;
      arma::mat L = chol(s2t * iVtV, "lower");
      arma::vec z(R, fill::randn);
      U_t.row(i) = (mu + L * z).t();
      if (bounded && arma::abs(U_t.row(i)).max() > max_abs) {
        U_t.row(i) = u_old;
        ++bound_rejections;
      }

      if(use_rho) {
        arma::vec Vdu = V_t * (U_t.row(i) - u_old).t();
        Resid.row(i) -= Vdu.t();
        Wrec -= Vdu * V_t.row(i);
      }
    }

    if(!symmetric) {
      arma::mat UtU = U_t.t() * U_t;
      arma::mat E_all_V = E_t.t() * U_t;

      arma::mat Mrec;
      if(use_rho) Mrec = Resid * U_t;

      for(int j = 0; j < n; j++) {
        arma::vec ej = E_all_V.row(j).t();
        arma::mat UtU_mod = UtU;
        arma::rowvec v_old = V_t.row(j);

        if(use_rho) {
          v_old = V_t.row(j);
          const arma::rowvec uj = U_t.row(j);
          ej -= rho * (Mrec.row(j).t() + arma::dot(v_old, uj) * uj.t());
          UtU_mod -= rho * (uj.t() * uj);
        }

        if(t > 0) {
          double lam = lambda_v(j, t);
          UtU_mod.diag() += lam * s2_sigma2_inv;
          ej += lam * s2_rho_sigma2 * V_new.slice(t-1).row(j).t();
        } else {
          UtU_mod.diag() += s2_stat_init;
        }
        if(t < T-1) {
          double lam_next = lambda_v(j, t+1);
          UtU_mod.diag() += lam_next * s2_rho2_sigma2;
          ej += lam_next * s2_rho_sigma2 * V_new.slice(t+1).row(j).t();
        }

        arma::mat iUtU;
        if(!safe_inv_sympd_uv(UtU_mod, iUtU)) Rcpp::stop("Unable to resolve the latent-factor Gaussian precision.");
        arma::vec mu = iUtU * ej;
        arma::mat L = chol(s2t * iUtU, "lower");
        arma::vec z(R, fill::randn);
        V_t.row(j) = (mu + L * z).t();
      if (bounded && arma::abs(V_t.row(j)).max() > max_abs) {
        V_t.row(j) = v_old;
        ++bound_rejections;
      }

        if(use_rho) {
          arma::vec Udv = U_t * (V_t.row(j) - v_old).t();
          Resid.col(j) -= Udv;
          Mrec -= Udv * U_t.row(j);
        }
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
  // a symmetric fit's product is U L U' with a signed diagonal L, which a
  // common rotation of the coordinates does not preserve; the directed
  // product U V' is rotation invariant
  if(shrink && !symmetric && !bounded) {
    common_rotation_uv(U_new, V_new);
  }

  if (bounded) return List::create(Named("U") = U_new, Named("V") = V_new,
                      Named("lambda_u") = lambda_u_new, Named("lambda_v") = lambda_v_new, Named("bound_rejections") = bound_rejections);
  return List::create(Named("U") = U_new, Named("V") = V_new,
                      Named("lambda_u") = lambda_u_new, Named("lambda_v") = lambda_v_new);
}
