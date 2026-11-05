#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

namespace {
// minimum-norm solution with the same relative eigenvalue cutoff as the r solver
arma::vec als_factor_solve(const arma::mat& gram, const arma::vec& rhs) {
	arma::vec values;
	arma::mat vectors;
	if (!arma::eig_sym(values, vectors, 0.5 * (gram + gram.t()))) {
		return arma::pinv(gram) * rhs;
	}
	double threshold = std::max(values.max(), 0.0) * std::sqrt(std::numeric_limits<double>::epsilon());
	arma::uvec keep = arma::find(values > threshold);
	if (keep.is_empty()) return arma::vec(rhs.n_elem, arma::fill::zeros);
	arma::mat basis = vectors.cols(keep);
	return basis * ((basis.t() * rhs) / values.elem(keep));
}

arma::vec als_factor_trace(const arma::mat& weights, const arma::mat& factors,
	const arma::vec& variance, bool rows) {
	arma::mat ridge = arma::diagmat(1.0 / arma::clamp(variance, 1e-8, arma::datum::inf));
	arma::vec total(factors.n_cols, arma::fill::zeros);
	arma::uword count = rows ? weights.n_rows : weights.n_cols;
	for (arma::uword i = 0; i < count; ++i) {
		arma::vec w = rows ? arma::vec(weights.row(i).t()) : arma::vec(weights.col(i));
		arma::mat gram = factors.t() * (factors.each_col() % w) + ridge;
		arma::mat inverse;
		if (!arma::inv(inverse, gram)) inverse = arma::pinv(gram);
		total += inverse.diag();
	}
	return total;
}

// ---- small-rank path ------------------------------------------------------
// the row and column systems of a sweep are rank x rank with a small rank, so
// their grams are accumulated for all rows (columns) in one pass over the
// weights and solved in place, in place of an eigendecomposition with its heap
// allocations per row. the minimum-norm eigen solution equals the exact
// solution whenever no eigenvalue falls below the cutoff (sqrt(eps) times the
// largest), which the cholesky factorization below certifies through
//   lambda_min / lambda_max >= det / trace^rank ;
// a system that cannot be certified takes the eigen path above.

const int kFastRank = 6;

// cholesky factor of the rank x rank symmetric matrix g (column-major) into
// the lower triangle of L; true when it is positive definite and the ratio of
// its smallest to its largest eigenvalue is certified above 1e-7
bool chol_certified(const double* g, int r, double* L) {
	double det = 1.0, tr = 0.0;
	for (int j = 0; j < r; ++j) {
		double s = g[j + r * j];
		tr += s;
		for (int k = 0; k < j; ++k) s -= L[j + r * k] * L[j + r * k];
		if (!(s > 0.0) || !std::isfinite(s)) return false;
		const double d = std::sqrt(s);
		L[j + r * j] = d;
		det *= s;
		for (int i = j + 1; i < r; ++i) {
			double t = g[i + r * j];
			for (int k = 0; k < j; ++k) t -= L[i + r * k] * L[j + r * k];
			L[i + r * j] = t / d;
		}
	}
	double trace_power = 1.0;
	for (int k = 0; k < r; ++k) trace_power *= tr;
	return std::isfinite(det) && det > 1e-7 * trace_power;
}

// solve L L' x = rhs
void chol_solve(const double* L, int r, const double* rhs, double* x) {
	double y[kFastRank];
	for (int i = 0; i < r; ++i) {
		double s = rhs[i];
		for (int k = 0; k < i; ++k) s -= L[i + r * k] * y[k];
		y[i] = s / L[i + r * i];
	}
	for (int i = r - 1; i >= 0; --i) {
		double s = y[i];
		for (int k = i + 1; k < r; ++k) s -= L[k + r * i] * x[k];
		x[i] = s / L[i + r * i];
	}
}

// the system g x = rhs: certified cholesky solve, else the eigen solution
void solve_small(const double* g, const double* rhs, int r, double* x) {
	double L[kFastRank * kFastRank];
	if (chol_certified(g, r, L)) {
		chol_solve(L, r, rhs, x);
		return;
	}
	arma::mat gram(r, r);
	arma::vec b(r);
	for (int k = 0; k < r; ++k) {
		b[k] = rhs[k];
		for (int l = 0; l < r; ++l) gram(k, l) = g[k + r * l];
	}
	arma::vec sol = als_factor_solve(gram, b);
	for (int k = 0; k < r; ++k) x[k] = sol[k];
}

// diagonal of the inverse of g: certified cholesky, else arma::inv / pinv
void inverse_diag_small(const double* g, int r, double* out) {
	double L[kFastRank * kFastRank];
	if (chol_certified(g, r, L)) {
		double e[kFastRank], col[kFastRank];
		for (int k = 0; k < r; ++k) {
			for (int l = 0; l < r; ++l) e[l] = (l == k) ? 1.0 : 0.0;
			chol_solve(L, r, e, col);
			out[k] = col[k];
		}
		return;
	}
	arma::mat gram(r, r), inverse;
	for (int k = 0; k < r; ++k) for (int l = 0; l < r; ++l) gram(k, l) = g[k + r * l];
	if (!arma::inv(inverse, gram)) inverse = arma::pinv(gram);
	for (int k = 0; k < r; ++k) out[k] = inverse(k, k);
}

// gram sums sum_j w_ij f_j f_j' and (when y is given) right-hand sides
// sum_j w_ij y_ij f_j for every row i, f_j being row j of `f` (nc x r).
// gram is laid out by upper-triangle pair then row (pair * nr + i), rhs by
// factor column then row (k * nr + i), so each accumulation is a contiguous
// pass over a column of w (and y)
void accumulate_rows(const double* w, const double* y, const double* f,
	int nr, int nc, int r, std::vector<double>& gram, std::vector<double>& rhs) {
	std::fill(gram.begin(), gram.end(), 0.0);
	if (y) std::fill(rhs.begin(), rhs.end(), 0.0);
	for (int j = 0; j < nc; ++j) {
		const double* wj = w + (std::size_t) nr * j;
		double fj[kFastRank];
		for (int k = 0; k < r; ++k) fj[k] = f[j + (std::size_t) nc * k];
		int pair = 0;
		for (int k = 0; k < r; ++k) {
			for (int l = k; l < r; ++l, ++pair) {
				const double c = fj[k] * fj[l];
				double* gp = gram.data() + (std::size_t) pair * nr;
				for (int i = 0; i < nr; ++i) gp[i] += wj[i] * c;
			}
		}
		if (y) {
			const double* yj = y + (std::size_t) nr * j;
			for (int k = 0; k < r; ++k) {
				double* rp = rhs.data() + (std::size_t) k * nr;
				for (int i = 0; i < nr; ++i) rp[i] += wj[i] * yj[i] * fj[k];
			}
		}
	}
}

// the same sums for every column j: sum_i w_ij f_i f_i' and sum_i w_ij y_ij f_i
// with f_i row i of `f` (nr x r); gram is (pair * nc + j), rhs (k * nc + j)
void accumulate_cols(const double* w, const double* y, const double* f,
	int nr, int nc, int r, std::vector<double>& gram, std::vector<double>& rhs) {
	for (int j = 0; j < nc; ++j) {
		const double* wj = w + (std::size_t) nr * j;
		int pair = 0;
		for (int k = 0; k < r; ++k) {
			const double* fk = f + (std::size_t) nr * k;
			for (int l = k; l < r; ++l, ++pair) {
				const double* fl = f + (std::size_t) nr * l;
				double s = 0.0;
				for (int i = 0; i < nr; ++i) s += wj[i] * fk[i] * fl[i];
				gram[(std::size_t) pair * nc + j] = s;
			}
		}
		if (y) {
			const double* yj = y + (std::size_t) nr * j;
			for (int k = 0; k < r; ++k) {
				const double* fk = f + (std::size_t) nr * k;
				double s = 0.0;
				for (int i = 0; i < nr; ++i) s += wj[i] * yj[i] * fk[i];
				rhs[(std::size_t) k * nc + j] = s;
			}
		}
	}
}

// expand the accumulated pair sums of unit `idx` (of `count`) into the full
// symmetric r x r matrix, adding the ridge to the diagonal
void expand_gram(const std::vector<double>& gram, int count, int idx, int r,
	const double* ridge, double* g) {
	int pair = 0;
	for (int k = 0; k < r; ++k) {
		for (int l = k; l < r; ++l, ++pair) {
			const double v = gram[(std::size_t) pair * count + idx];
			g[k + r * l] = v;
			g[l + r * k] = v;
		}
		g[k + r * k] += ridge[k];
	}
}

inline double clamped_scale(double s) {
	return std::min(std::max(s, 1e-8), arma::datum::inf);
}

Rcpp::List als_factor_update_small(const arma::mat& response, const arma::mat& weights,
	const arma::mat& u_init, const arma::mat& v_init, const arma::vec& su_init,
	const arma::vec& sv_init, double prior_scale, double prior_df, int maxit, double tol) {
	const int nr = (int) response.n_rows, nc = (int) response.n_cols;
	const int r = (int) u_init.n_cols;
	const int npair = r * (r + 1) / 2;
	const double* w = weights.memptr();
	const double* y = response.memptr();

	arma::mat u = u_init, v = v_init;
	std::vector<double> gram_r((std::size_t) npair * nr), rhs_r((std::size_t) r * nr);
	std::vector<double> gram_c((std::size_t) npair * nc), rhs_c((std::size_t) r * nc);

	double ridge[kFastRank], g[kFastRank * kFastRank], b[kFastRank], x[kFastRank], dg[kFastRank];

	// posterior-variance trace terms from the incoming factors
	arma::vec su(r), sv(r);
	{
		for (int k = 0; k < r; ++k) ridge[k] = 1.0 / clamped_scale(su_init[k]);
		accumulate_rows(w, nullptr, v.memptr(), nr, nc, r, gram_r, rhs_r);
		std::vector<double> total(r, 0.0);
		for (int i = 0; i < nr; ++i) {
			expand_gram(gram_r, nr, i, r, ridge, g);
			inverse_diag_small(g, r, dg);
			for (int k = 0; k < r; ++k) total[k] += dg[k];
		}
		arma::vec uss = arma::sum(arma::square(u), 0).t();
		for (int k = 0; k < r; ++k) su[k] = (prior_scale + uss[k] + total[k]) / (nr + prior_df);
	}
	{
		for (int k = 0; k < r; ++k) ridge[k] = 1.0 / clamped_scale(sv_init[k]);
		accumulate_cols(w, nullptr, u.memptr(), nr, nc, r, gram_c, rhs_c);
		std::vector<double> total(r, 0.0);
		for (int j = 0; j < nc; ++j) {
			expand_gram(gram_c, nc, j, r, ridge, g);
			inverse_diag_small(g, r, dg);
			for (int k = 0; k < r; ++k) total[k] += dg[k];
		}
		arma::vec vss = arma::sum(arma::square(v), 0).t();
		for (int k = 0; k < r; ++k) sv[k] = (prior_scale + vss[k] + total[k]) / (nc + prior_df);
	}
	double lu[kFastRank], lv[kFastRank];
	for (int k = 0; k < r; ++k) {
		lu[k] = 1.0 / clamped_scale(su[k]);
		lv[k] = 1.0 / clamped_scale(sv[k]);
	}

	arma::mat product(nr, nc, arma::fill::zeros);
	double old_objective = arma::datum::inf;
	int iterations = 0;
	for (int it = 0; it < maxit; ++it) {
		Rcpp::checkUserInterrupt();
		// rows: u_i = (v' diag(w_i) v + ridge)^+ v' (w_i y_i)
		accumulate_rows(w, y, v.memptr(), nr, nc, r, gram_r, rhs_r);
		for (int i = 0; i < nr; ++i) {
			expand_gram(gram_r, nr, i, r, lu, g);
			for (int k = 0; k < r; ++k) b[k] = rhs_r[(std::size_t) k * nr + i];
			solve_small(g, b, r, x);
			for (int k = 0; k < r; ++k) u(i, k) = x[k];
		}
		// columns, from the updated rows
		accumulate_cols(w, y, u.memptr(), nr, nc, r, gram_c, rhs_c);
		for (int j = 0; j < nc; ++j) {
			expand_gram(gram_c, nc, j, r, lv, g);
			for (int k = 0; k < r; ++k) b[k] = rhs_c[(std::size_t) k * nc + j];
			solve_small(g, b, r, x);
			for (int k = 0; k < r; ++k) v(j, k) = x[k];
		}
		product = u * v.t();
		double data_term = 0.0;
		{
			const double* pp = product.memptr();
			const std::size_t n = (std::size_t) nr * nc;
			for (std::size_t q = 0; q < n; ++q) {
				const double d = y[q] - pp[q];
				data_term += w[q] * d * d;
			}
		}
		double penalty = 0.0;
		for (int k = 0; k < r; ++k) {
			double uu = 0.0, vv = 0.0;
			for (int i = 0; i < nr; ++i) uu += u(i, k) * u(i, k);
			for (int j = 0; j < nc; ++j) vv += v(j, k) * v(j, k);
			penalty += uu / clamped_scale(su[k]) + vv / clamped_scale(sv[k]);
		}
		const double objective = data_term + penalty;
		++iterations;
		bool done = std::isfinite(old_objective) &&
			std::abs(old_objective - objective) / std::max(1.0, std::abs(old_objective)) < tol;
		old_objective = objective;
		if (done) break;
	}
	return Rcpp::List::create(Rcpp::Named("U") = u, Rcpp::Named("V") = v,
		Rcpp::Named("L") = R_NilValue, Rcpp::Named("Omat") = product,
		Rcpp::Named("sig2_u") = Rcpp::NumericVector(su.begin(), su.end()),
		Rcpp::Named("sig2_v") = Rcpp::NumericVector(sv.begin(), sv.end()),
		Rcpp::Named("obj") = old_objective, Rcpp::Named("iters") = iterations);
}
}

namespace {
Rcpp::List als_factor_update_core(const arma::mat& response, const arma::mat& weights,
	const arma::mat& u_init, const arma::mat& v_init, const arma::vec& su_init,
	const arma::vec& sv_init, double prior_scale, double prior_df,
	int maxit, double tol) {
	const arma::uword nr = response.n_rows, nc = response.n_cols, rank = u_init.n_cols;
	if (weights.n_rows != nr || weights.n_cols != nc || u_init.n_rows != nr ||
		v_init.n_rows != nc || v_init.n_cols != rank || su_init.n_elem != rank ||
		sv_init.n_elem != rank || rank == 0 || maxit < 1) {
		Rcpp::stop("Incompatible dimensions for the ALS factor update.");
	}
	if (rank <= (arma::uword) kFastRank) {
		return als_factor_update_small(response, weights, u_init, v_init, su_init,
			sv_init, prior_scale, prior_df, maxit, tol);
	}
	arma::mat u = u_init, v = v_init;
	arma::vec su = (prior_scale + arma::sum(arma::square(u), 0).t() +
		als_factor_trace(weights, v, su_init, true)) / (nr + prior_df);
	arma::vec sv = (prior_scale + arma::sum(arma::square(v), 0).t() +
		als_factor_trace(weights, u, sv_init, false)) / (nc + prior_df);
	arma::mat lu = arma::diagmat(1.0 / arma::clamp(su, 1e-8, arma::datum::inf));
	arma::mat lv = arma::diagmat(1.0 / arma::clamp(sv, 1e-8, arma::datum::inf));
	arma::mat product(nr, nc, arma::fill::zeros);
	double old_objective = arma::datum::inf;
	int iterations = 0;
	for (int it = 0; it < maxit; ++it) {
		Rcpp::checkUserInterrupt();
		for (arma::uword i = 0; i < nr; ++i) {
			arma::vec w = weights.row(i).t();
			arma::mat gram = v.t() * (v.each_col() % w) + lu;
			u.row(i) = als_factor_solve(gram, v.t() * (w % response.row(i).t())).t();
		}
		for (arma::uword j = 0; j < nc; ++j) {
			arma::vec w = weights.col(j);
			arma::mat gram = u.t() * (u.each_col() % w) + lv;
			v.row(j) = als_factor_solve(gram, u.t() * (w % response.col(j))).t();
		}
		product = u * v.t();
		double objective = arma::accu(weights % arma::square(response - product)) +
			arma::dot(arma::sum(arma::square(u), 0).t(), 1.0 / arma::clamp(su, 1e-8, arma::datum::inf)) +
			arma::dot(arma::sum(arma::square(v), 0).t(), 1.0 / arma::clamp(sv, 1e-8, arma::datum::inf));
		++iterations;
		bool done = std::isfinite(old_objective) &&
			std::abs(old_objective - objective) / std::max(1.0, std::abs(old_objective)) < tol;
		old_objective = objective;
		if (done) break;
	}
	return Rcpp::List::create(Rcpp::Named("U") = u, Rcpp::Named("V") = v,
		Rcpp::Named("L") = R_NilValue, Rcpp::Named("Omat") = product,
		Rcpp::Named("sig2_u") = Rcpp::NumericVector(su.begin(), su.end()),
		Rcpp::Named("sig2_v") = Rcpp::NumericVector(sv.begin(), sv.end()),
		Rcpp::Named("obj") = old_objective, Rcpp::Named("iters") = iterations);
}
}

// penalized factor sweeps with fixed hyperparameters during each m-step
// [[Rcpp::export(rng = false)]]
Rcpp::List als_factor_update_cpp(const arma::mat& response, const arma::mat& weights,
	const arma::mat& u_init, const arma::mat& v_init, const arma::vec& su_init,
	const arma::vec& sv_init, double prior_scale, double prior_df,
	int maxit = 15, double tol = 1e-6) {
	return als_factor_update_core(response, weights, u_init, v_init, su_init, sv_init,
		prior_scale, prior_df, maxit, tol);
}

// ---- balanced factors -----------------------------------------------------
// the penalized update followed by the balancing of the two factors that
// .ae_lowrank_em_als applies in r (R/ame_als_lowrank.R): the factors are
// rotated, at a fixed product u v', so that they carry equal scale under
// their own prior variances. the arithmetic is the r code's -- a householder
// qr of each factor (linpack's dqrdc2, with its limited pivoting), the svd of
// the product of the r factors, and the assignment of its singular directions
// to the factor columns in the order of the prior variances. the r path is
// kept for every state this routine does not take (a column order that
// differs from the variance order, a failed svd): the result then carries
// balanced = false and the unbalanced update.
namespace {
// positions of the values in decreasing order; ties keep their original
// order and non-finite values come last, as order(x, decreasing = TRUE)
std::vector<int> order_decreasing(const std::vector<double>& v) {
	std::vector<int> idx(v.size());
	for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = (int) i;
	std::stable_sort(idx.begin(), idx.end(), [&v](int a, int b) {
		const bool na = std::isnan(v[a]), nb = std::isnan(v[b]);
		if (na || nb) return !na && nb;
		return v[a] > v[b];
	});
	return idx;
}

double nrm2(int n, const double* x) {
	double s = 0.0;
	for (int i = 0; i < n; ++i) s += x[i] * x[i];
	return std::sqrt(s);
}

double dot(int n, const double* x, const double* y) {
	double s = 0.0;
	for (int i = 0; i < n; ++i) s += x[i] * y[i];
	return s;
}

// linpack dqrdc2 as qr() calls it (tol = 1e-7): the reduced matrix x (n x p,
// column-major), the householder scalars, the column order after the limited
// pivoting (a column whose remaining norm has fallen below tol times its
// original norm is cycled to the end) and the rank. false only when n <= p
struct LinpackQr {
	std::vector<double> x, qraux;
	std::vector<int> pivot;   // 0-based original column of each position
	int rank;
};

bool qr_linpack(const arma::mat& a, LinpackQr& q) {
	const int n = (int) a.n_rows, p = (int) a.n_cols;
	if (n <= p) return false;
	const double tol = 1e-7;
	std::vector<double>& x = q.x;
	std::vector<double>& qraux = q.qraux;
	x.assign(a.memptr(), a.memptr() + (std::size_t) n * p);
	qraux.assign(p, 0.0);
	q.pivot.resize(p);
	for (int j = 0; j < p; ++j) q.pivot[j] = j;
	std::vector<double> work2(p);
	for (int j = 0; j < p; ++j) {
		qraux[j] = nrm2(n, &x[(std::size_t) n * j]);
		work2[j] = qraux[j];
		if (work2[j] == 0.0) work2[j] = 1.0;
	}
	int k = p + 1;               // fortran's k (1-based); columns l >= k are the negligible ones
	for (int l = 1; l <= p; ++l) {   // fortran l (1-based), lup = min(n, p) = p
		// cycle negligible columns from l to the end
		while (!(l >= k || qraux[l - 1] >= work2[l - 1] * tol)) {
			for (int i = 0; i < n; ++i) {
				const double t = x[i + (std::size_t) n * (l - 1)];
				for (int j = l; j < p; ++j) x[i + (std::size_t) n * (j - 1)] = x[i + (std::size_t) n * j];
				x[i + (std::size_t) n * (p - 1)] = t;
			}
			const int pv = q.pivot[l - 1];
			const double tq = qraux[l - 1], tw = work2[l - 1];
			for (int j = l; j < p; ++j) {
				q.pivot[j - 1] = q.pivot[j];
				qraux[j - 1] = qraux[j];
				work2[j - 1] = work2[j];
			}
			q.pivot[p - 1] = pv;
			qraux[p - 1] = tq;
			work2[p - 1] = tw;
			--k;
		}
		const int l0 = l - 1;       // 0-based
		if (l == n) continue;
		double* xl = &x[l0 + (std::size_t) n * l0];
		double nrmxl = nrm2(n - l0, xl);
		if (nrmxl == 0.0) continue;
		if (xl[0] != 0.0) nrmxl = std::copysign(nrmxl, xl[0]);
		const double scale = 1.0 / nrmxl;
		for (int i = 0; i < n - l0; ++i) xl[i] *= scale;
		xl[0] += 1.0;
		for (int j = l0 + 1; j < p; ++j) {
			double* xj = &x[l0 + (std::size_t) n * j];
			const double t = -dot(n - l0, xl, xj) / xl[0];
			for (int i = 0; i < n - l0; ++i) xj[i] += t * xl[i];
			if (qraux[j] != 0.0) {
				const double qq = std::fabs(xj[0]) / qraux[j];
				double tt = 1.0 - qq * qq;
				tt = std::max(tt, 0.0);
				if (std::fabs(tt) < 1e-6) {
					qraux[j] = nrm2(n - l0 - 1, xj + 1);
				} else {
					qraux[j] *= std::sqrt(tt);
				}
			}
		}
		qraux[l0] = xl[0];
		xl[0] = -nrmxl;
	}
	q.rank = std::min(k - 1, n);
	return true;
}

// y <- Q y for the first `rank` householder reflections of a reduced qr
// (qr.qy applies min(rank, n - 1) of them); y is n x m
void qr_apply_q(LinpackQr& q, int n, arma::mat& y) {
	const int m = (int) y.n_cols;
	const int ju = std::min(q.rank, n - 1);
	for (int l = ju - 1; l >= 0; --l) {
		if (q.qraux[l] == 0.0) continue;
		double* xl = &q.x[l + (std::size_t) n * l];
		const double saved = xl[0];
		xl[0] = q.qraux[l];
		for (int c = 0; c < m; ++c) {
			double* yc = y.colptr(c) + l;
			const double t = -dot(n - l, xl, yc) / xl[0];
			for (int i = 0; i < n - l; ++i) yc[i] += t * xl[i];
		}
		xl[0] = saved;
	}
}

// upper triangle of the reduced matrix, columns put back in their original
// order (qr.R(qr)[, order(qr$pivot)])
arma::mat qr_r_factor(const LinpackQr& q, int n, int p) {
	arma::mat r(p, p, arma::fill::zeros);
	for (int j = 0; j < p; ++j) {
		const int orig = q.pivot[j];
		for (int i = 0; i <= j; ++i) r(i, orig) = q.x[i + (std::size_t) n * j];
	}
	return r;
}

// sqrt(colSums(f^2)) as r evaluates it: squares rounded, summed in long double
std::vector<double> column_norms_sq(const arma::mat& f) {
	std::vector<double> out(f.n_cols);
	for (arma::uword k = 0; k < f.n_cols; ++k) {
		long double acc = 0.0L;
		const double* c = f.colptr(k);
		for (arma::uword i = 0; i < f.n_rows; ++i) {
			const double sq = c[i] * c[i];
			acc += sq;
		}
		out[k] = (double) acc;
	}
	return out;
}
}

// [[Rcpp::export(rng = false)]]
Rcpp::List als_factor_balanced_cpp(const arma::mat& response, const arma::mat& weights,
	const arma::mat& u_init, const arma::mat& v_init, const arma::vec& su_init,
	const arma::vec& sv_init, double prior_scale, double prior_df,
	int maxit = 15, double tol = 1e-6) {
	const arma::uword nr = response.n_rows, nc = response.n_cols, rank = u_init.n_cols;
	if (weights.n_rows != nr || weights.n_cols != nc || u_init.n_rows != nr ||
		v_init.n_rows != nc || v_init.n_cols != rank || su_init.n_elem != rank ||
		sv_init.n_elem != rank || rank == 0 || maxit < 1) {
		Rcpp::stop("Incompatible dimensions for the ALS factor update.");
	}
	Rcpp::List out = als_factor_update_core(response, weights, u_init, v_init, su_init,
		sv_init, prior_scale, prior_df, maxit, tol);
	out["balanced"] = false;
	const int r = (int) rank;
	if (rank > (arma::uword) kFastRank) return out;
	arma::mat u = Rcpp::as<arma::mat>(out["U"]), v = Rcpp::as<arma::mat>(out["V"]);
	std::vector<double> su = Rcpp::as<std::vector<double> >(out["sig2_u"]);
	std::vector<double> sv = Rcpp::as<std::vector<double> >(out["sig2_v"]);
	// the variance order and the factor-norm order must agree
	std::vector<double> key(r), fkey(r);
	std::vector<double> uu = column_norms_sq(u), vv = column_norms_sq(v);
	for (int k = 0; k < r; ++k) {
		key[k] = std::sqrt(su[k] * sv[k]);
		fkey[k] = std::sqrt(uu[k] * vv[k]);
	}
	std::vector<int> index = order_decreasing(key);
	if (index != order_decreasing(fkey)) return out;
	LinpackQr qu, qv;
	if (!qr_linpack(u, qu) || !qr_linpack(v, qv)) return out;
	arma::mat ru = qr_r_factor(qu, (int) nr, r), rv = qr_r_factor(qv, (int) nc, r);
	arma::mat product = ru * rv.t();
	arma::mat svd_u, svd_v;
	arma::vec d;
	if (!arma::svd(svd_u, d, svd_v, product, "dc")) return out;
	arma::mat yu(nr, r, arma::fill::zeros), yv(nc, r, arma::fill::zeros);
	yu.rows(0, r - 1) = svd_u;
	yv.rows(0, r - 1) = svd_v;
	qr_apply_q(qu, (int) nr, yu);
	qr_apply_q(qv, (int) nc, yv);
	arma::mat ub(nr, r, arma::fill::zeros), vb(nc, r, arma::fill::zeros);
	for (int k = 0; k < r; ++k) {
		const int col = index[k];
		const double ratio = std::pow(su[col] / sv[col], 0.25);
		const double fu = std::sqrt(d[k]) * ratio, fv = std::sqrt(d[k]) / ratio;
		for (arma::uword i = 0; i < nr; ++i) ub(i, col) = yu(i, k) * fu;
		for (arma::uword j = 0; j < nc; ++j) vb(j, col) = yv(j, k) * fv;
	}
	out["U"] = ub;
	out["V"] = vb;
	out["Omat"] = arma::mat(ub * vb.t());
	out["balanced"] = true;
	return out;
}
