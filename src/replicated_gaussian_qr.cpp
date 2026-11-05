// Joint Gaussian full conditional of regression coefficients and additive
// actor effects for replicated (longitudinal) unipartite relational data.
//
// Each period's cell grid is split into its symmetric and antisymmetric
// parts, whitened by the reciprocal variances s2 (1 + rho) and s2 (1 - rho).
// An orthonormal Helmert basis separates the actor directions: every
// Helmert direction h carries one symmetric and one antisymmetric
// observation per period on its k actor coordinates, and the overall-mean
// direction carries one symmetric observation. Covariate columns are split
// into their projections on these directions and an orthogonal residual.
// The actor coordinates are eliminated block by block with a small QR, the
// coefficients are drawn from their marginal Gaussian (column-equilibrated
// QR of the stacked square-root information), and the actor coordinates are
// then drawn given the coefficients. Every design-only quantity is computed
// once by replicated_gaussian_design_cpp and reused across sweeps.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
using namespace arma;
using namespace Rcpp;

namespace {

// H' v for the orthonormal Helmert basis H (n x (n-1)): column h has
// 1/sqrt((h+1)(h+2)) on entries 0..h and -(h+1)/sqrt((h+1)(h+2)) at h+1.
vec helmert_t(const vec& v) {
	const uword n = v.n_elem;
	vec out(n ? n - 1 : 0);
	double cumulative = 0;
	for (uword h = 0; h + 1 < n; ++h) {
		cumulative += v(h);
		out(h) = (cumulative - double(h + 1) * v(h + 1)) /
			std::sqrt(double(h + 1) * double(h + 2));
	}
	return out;
}

// H f, the length-n actor vector of Helmert coordinates f (length n-1).
vec helmert(const vec& f, uword n) {
	vec out(n, fill::zeros);
	double tail = 0;
	for (uword ii = n; ii > 0; --ii) {
		const uword i = ii - 1;
		if (i + 1 < n) tail += f(i) / std::sqrt(double(i + 1) * double(i + 2));
		out(i) = tail;
		if (i >= 1) out(i) -= double(i) * f(i - 1) / std::sqrt(double(i) * double(i + 1));
	}
	return out;
}

// Root G (2 x k) of the actor covariance, k = rank. A positive definite
// covariance uses its Cholesky factor; otherwise the positive eigen
// directions are kept, so disabled or tied margins stay exactly supported.
mat actor_root(const mat& Sab) {
	mat S = 0.5 * (Sab + Sab.t()), G;
	if (S(0, 0) > 0 && S(1, 1) > 0 && chol(G, S, "lower")) return G;
	vec ev; mat V;
	if (!eig_sym(ev, V, S)) Rcpp::stop("Actor covariance eigendecomposition failed.");
	const double tolerance = 4 * std::numeric_limits<double>::epsilon() * std::max(abs(ev).max(), 0.0);
	uvec keep = find(ev > tolerance);
	G = V.cols(keep);
	if (keep.n_elem) G.each_row() %= sqrt(ev.elem(keep)).t();
	return G;
}

// Square root C (C'C = P) of a symmetric positive semidefinite precision.
mat precision_root(const mat& P) {
	const uword p = P.n_rows;
	mat S = 0.5 * (P + P.t()), R;
	if (!p || chol(R, S)) return R;
	vec ev; mat V;
	if (!eig_sym(ev, V, S)) Rcpp::stop("Coefficient prior eigendecomposition failed.");
	ev.elem(find(ev < 0)).zeros();
	R = V.t();
	R.each_col() %= sqrt(ev);
	return R;
}

// Gaussian with square-root information rows A theta ~ y (precision A'A,
// linear term A'y): returns the mean and a root L with covariance L L'.
// Columns are equilibrated before an orthogonal factorization, so the draw
// is invariant to the units of each coefficient and the row weights never
// enter through a squared normal-equation sum. A system that is not
// numerically of full rank falls back to the eigendecomposition of the
// equilibrated precision with a relative eigenvalue floor.
void gaussian_from_rows(mat A, const vec& y, vec& mean, mat& root) {
	const uword p = A.n_cols;
	if (!p) { mean.reset(); root.reset(); return; }
	vec d(p);
	for (uword j = 0; j < p; ++j) {
		d(j) = norm(A.col(j), 2);
		if (!(d(j) > 0) || !std::isfinite(d(j))) d(j) = 1;
	}
	A.each_row() /= d.t();
	mat Q, R;
	bool full = A.n_rows >= p && qr_econ(Q, R, A);
	if (full) {
		vec rd = abs(R.diag());
		full = rd.min() > 16 * p * std::numeric_limits<double>::epsilon() * rd.max();
	}
	if (full) {
		root = solve(trimatu(R), eye(p, p));
		mean = root * (Q.t() * y);
	} else {
		mat Ps = A.t() * A;
		Ps = 0.5 * (Ps + Ps.t());
		vec ev; mat V;
		if (!eig_sym(ev, V, Ps)) Rcpp::stop("Coefficient precision eigendecomposition failed.");
		const double floor = std::max(ev.max(), 1.0) * std::sqrt(std::numeric_limits<double>::epsilon());
		ev.elem(find(ev < floor)).fill(floor);
		mean = V * ((V.t() * (A.t() * y)) / ev);
		root = V;
		root.each_row() /= sqrt(ev).t();
	}
	mean /= d;
	root.each_col() /= d;
}

// Design cache as views on the R-owned memory (no copy per sweep).
cube cube_view(const List& design, const char* name) {
	NumericVector x = design[name];
	IntegerVector dims = x.attr("dim");
	if (dims.size() != 3) Rcpp::stop("Invalid Gaussian design cache.");
	return cube(x.begin(), dims[0], dims[1], dims[2], false, true);
}

struct Design {
	uword n, p, nt;
	cube QP, QM, Rp, Rm, bs, ba;
	mat b0;
	explicit Design(const List& design) :
		n(as<uword>(design["n"])), p(as<uword>(design["p"])), nt(as<uword>(design["nt"])),
		QP(cube_view(design, "QP")), QM(cube_view(design, "QM")),
		Rp(cube_view(design, "Rp")), Rm(cube_view(design, "Rm")),
		bs(cube_view(design, "bs")), ba(cube_view(design, "ba")),
		b0(as<mat>(design["b0"])) {}
};

struct System {
	uword n, p, k, nr;
	bool center;
	mat G;              // 2 x k actor root
	vec beta_mean;      // p
	mat beta_root;      // p x p, covariance beta_root beta_root'
	mat Rc, R0;         // k x k block roots (Helmert, common)
	cube C;             // k x p x nr coefficient coupling
	mat D;              // k x nr block responses
};

System build_system(const Design& ds, const cube& Z, const mat& Sab, const vec& rho_in,
		double s2, const mat& P0, const vec& m0, bool center) {
	const uword n = ds.n, p = ds.p, nt = ds.nt;
	if (Z.n_rows != n || Z.n_cols != n || Z.n_slices != nt)
		Rcpp::stop("Gaussian response dimensions disagree with the design.");
	if (!Z.is_finite()) Rcpp::stop("Gaussian response must be finite.");
	if (Sab.n_rows != 2 || Sab.n_cols != 2 || !Sab.is_finite())
		Rcpp::stop("Actor covariance must be a finite 2 x 2 matrix.");
	if (!std::isfinite(s2) || s2 <= 0) Rcpp::stop("Residual variance must be positive and finite.");
	if (P0.n_rows != p || P0.n_cols != p || m0.n_elem != p || !P0.is_finite() || !m0.is_finite())
		Rcpp::stop("Coefficient prior must be finite and match the design.");
	vec rho = rho_in.n_elem == 1 ? vec(nt, fill::value(rho_in(0))) : rho_in;
	if (rho.n_elem != nt || !rho.is_finite() || any(abs(rho) >= 1))
		Rcpp::stop("Dyadic correlation must lie strictly between -1 and 1.");

	System sy;
	sy.n = n; sy.p = p; sy.center = center;
	sy.G = actor_root(Sab);
	const uword k = sy.G.n_cols, nh = n - 1;
	sy.k = k; sy.nr = center ? nh : n;
	vec wp = 1 / sqrt(s2 * (1 + rho)), wm = 1 / sqrt(s2 * (1 - rho));
	const double sn = std::sqrt(2.0 / n);

	// square-root information rows of the coefficient marginal
	std::vector<mat> rows;
	std::vector<vec> resp;
	if (p) {
		mat C = precision_root(P0);
		rows.push_back(C);
		resp.push_back(C * m0);
	}
	// period projections of the response
	mat ys(nh, nt), ya(nh, nt);
	vec z0(nt);
	const uword npair_p = n * (n + 1) / 2, npair_m = n * (n - 1) / 2;
	vec sp(npair_p), sm(npair_m);
	const double r2 = std::sqrt(2.0);
	for (uword t = 0; t < nt; ++t) {
		const mat& Zt = Z.slice(t);
		vec rs = sum(Zt, 1), cs = sum(Zt, 0).t();
		z0(t) = accu(rs) / n;
		ys.col(t) = sn * helmert_t(0.5 * (rs + cs));
		ya.col(t) = sn * helmert_t(0.5 * (rs - cs));
		if (p) {
			uword ip = 0, im = 0;
			for (uword j = 0; j < n; ++j) {
				for (uword i = 0; i < j; ++i) {
					sp(ip++) = (Zt(i, j) + Zt(j, i)) / r2;
					sm(im++) = (Zt(i, j) - Zt(j, i)) / r2;
				}
				sp(ip++) = Zt(j, j);
			}
			rows.push_back(wp(t) * ds.Rp.slice(t));
			resp.push_back(wp(t) * (ds.QP.slice(t).t() * sp));
			rows.push_back(wm(t) * ds.Rm.slice(t));
			resp.push_back(wm(t) * (ds.QM.slice(t).t() * sm));
		}
	}

	// Helmert actor blocks: rows (2t, 2t+1) are the symmetric and
	// antisymmetric observations of period t
	sy.C.zeros(k, p, sy.nr); sy.D.zeros(k, sy.nr);
	mat Qt;  // (2nt + k) x 2nt: transpose of the observed rows of Q
	if (k) {
		mat E(2 * nt + k, k, fill::zeros);
		const double sh = std::sqrt(n / 2.0);
		for (uword t = 0; t < nt; ++t) {
			E.row(2 * t) = sh * wp(t) * (sy.G.row(0) + sy.G.row(1));
			E.row(2 * t + 1) = sh * wm(t) * (sy.G.row(0) - sy.G.row(1));
		}
		E.rows(2 * nt, 2 * nt + k - 1) = eye(k, k);
		mat Q, R;
		if (!qr(Q, R, E)) Rcpp::stop("Actor block QR failed.");
		sy.Rc = R.rows(0, k - 1);
		Qt = Q.rows(0, 2 * nt - 1).t();
	}
	mat B(2 * nt, p);
	vec y(2 * nt);
	for (uword h = 0; h < nh; ++h) {
		for (uword t = 0; t < nt; ++t) {
			if (p) {
				B.row(2 * t) = wp(t) * ds.bs.slice(t).row(h);
				B.row(2 * t + 1) = wm(t) * ds.ba.slice(t).row(h);
			}
			y(2 * t) = wp(t) * ys(h, t);
			y(2 * t + 1) = wm(t) * ya(h, t);
		}
		if (k) {
			mat TB = Qt * B;
			vec ty = Qt * y;
			sy.C.slice(h) = TB.rows(0, k - 1);
			sy.D.col(h) = ty.head(k);
			if (p) {
				rows.push_back(TB.rows(k, TB.n_rows - 1));
				resp.push_back(ty.tail(ty.n_elem - k));
			}
		} else if (p) {
			rows.push_back(B);
			resp.push_back(y);
		}
	}
	// overall-mean direction: one symmetric observation per period
	mat B0(nt, p);
	vec y0(nt);
	for (uword t = 0; t < nt; ++t) {
		if (p) B0.row(t) = wp(t) * ds.b0.row(t);
		y0(t) = wp(t) * z0(t);
	}
	if (!center && k) {
		mat E(nt + k, k, fill::zeros);
		for (uword t = 0; t < nt; ++t)
			E.row(t) = std::sqrt(double(n)) * wp(t) * (sy.G.row(0) + sy.G.row(1));
		E.rows(nt, nt + k - 1) = eye(k, k);
		mat Q, R;
		if (!qr(Q, R, E)) Rcpp::stop("Actor mean QR failed.");
		sy.R0 = R.rows(0, k - 1);
		mat Q0t = Q.rows(0, nt - 1).t();
		mat TB = Q0t * B0;
		vec ty = Q0t * y0;
		sy.C.slice(n - 1) = TB.rows(0, k - 1);
		sy.D.col(n - 1) = ty.head(k);
		if (p) {
			rows.push_back(TB.rows(k, TB.n_rows - 1));
			resp.push_back(ty.tail(ty.n_elem - k));
		}
	} else if (p) {
		rows.push_back(B0);
		resp.push_back(y0);
	}
	if (p) {
		uword m = 0;
		for (const mat& r : rows) m += r.n_rows;
		mat A(m, p);
		vec yy(m);
		uword at = 0;
		for (uword i = 0; i < rows.size(); ++i) {
			if (!rows[i].n_rows) continue;
			A.rows(at, at + rows[i].n_rows - 1) = rows[i];
			yy.subvec(at, at + rows[i].n_rows - 1) = resp[i];
			at += rows[i].n_rows;
		}
		gaussian_from_rows(A, yy, sy.beta_mean, sy.beta_root);
	}
	if (!sy.beta_mean.is_finite() || !sy.beta_root.is_finite())
		Rcpp::stop("Coefficient conditional moments are not finite.");
	return sy;
}

// actor block roots and their inverses, Helmert blocks first, then the mean
mat block_inverse(const System& sy, uword h) {
	const mat& R = (!sy.center && h == sy.n - 1) ? sy.R0 : sy.Rc;
	return solve(trimatu(R), eye(sy.k, sy.k));
}

// physical additive effects from the block coordinates F (nr x k)
void actor_effects(const System& sy, const mat& F, vec& a, vec& b) {
	const uword n = sy.n;
	a.zeros(n); b.zeros(n);
	if (!sy.k) return;
	mat W = F * sy.G.t();  // nr x 2
	if (n > 1) {
		a = helmert(W.col(0).head(n - 1), n);
		b = helmert(W.col(1).head(n - 1), n);
	}
	if (!sy.center) {
		a += W(n - 1, 0) / std::sqrt(double(n));
		b += W(n - 1, 1) / std::sqrt(double(n));
	}
}

// thin QR of an m x p block; a block with fewer rows than columns keeps
// zero padding, so Q' y and R stay p-dimensional
void qr_factor(const mat& X, mat& Q, mat& R) {
	const uword m = X.n_rows, p = X.n_cols;
	Q.zeros(m, p); R.zeros(p, p);
	if (!m || !p) return;
	mat q, r;
	if (!qr_econ(q, r, X)) Rcpp::stop("Design residual QR failed.");
	Q.cols(0, q.n_cols - 1) = q;
	R.rows(0, r.n_rows - 1) = r;
}

}  // namespace

// Design-only quantities: residual covariate columns in pair coordinates,
// their Gram matrices, and the projections on the actor directions.
// [[Rcpp::export(rng = false)]]
List replicated_gaussian_design_cpp(const arma::cube& mX, int n) {
	if (n < 1 || mX.n_rows != uword(n) * uword(n)) Rcpp::stop("Design must have n^2 rows.");
	if (!mX.is_finite()) Rcpp::stop("Design must be finite.");
	const uword N = n, p = mX.n_cols, nt = mX.n_slices, nh = N - 1;
	const uword npair_p = N * (N + 1) / 2, npair_m = N * (N - 1) / 2;
	const double sn = std::sqrt(2.0 / N), s2n = std::sqrt(2.0 * N), r2 = std::sqrt(2.0);
	cube RP(npair_p, p, nt), RM(npair_m, p, nt);
	cube QP(npair_p, p, nt), QM(npair_m, p, nt), Rp(p, p, nt), Rm(p, p, nt);
	cube bs(nh, p, nt), ba(nh, p, nt);
	mat b0(nt, p);
	for (uword t = 0; t < nt; ++t) {
		for (uword v = 0; v < p; ++v) {
			mat x(mX.slice(t).colptr(v), N, N);
			vec rs = sum(x, 1), cs = sum(x, 0).t();
			b0(t, v) = accu(rs) / N;
			vec s = sn * helmert_t(0.5 * (rs + cs)), a = sn * helmert_t(0.5 * (rs - cs));
			bs.slice(t).col(v) = s;
			ba.slice(t).col(v) = a;
			vec u = helmert(s, N), w = helmert(a, N);
			const double c0 = b0(t, v) / N;
			uword ip = 0, im = 0;
			for (uword j = 0; j < N; ++j) {
				for (uword i = 0; i < j; ++i) {
					RP(ip++, v, t) = r2 * (0.5 * (x(i, j) + x(j, i)) - c0 - (u(i) + u(j)) / s2n);
					RM(im++, v, t) = r2 * (0.5 * (x(i, j) - x(j, i)) - (w(i) - w(j)) / s2n);
				}
				RP(ip++, v, t) = x(j, j) - c0 - 2 * u(j) / s2n;
			}
		}
		// orthogonal factors of the residual columns: rows of the square-root
		// information and the response map, both independent of rho and s2
		qr_factor(RP.slice(t), QP.slice(t), Rp.slice(t));
		qr_factor(RM.slice(t), QM.slice(t), Rm.slice(t));
	}
	return List::create(_["n"] = N, _["p"] = p, _["nt"] = nt, _["QP"] = QP, _["QM"] = QM,
		_["Rp"] = Rp, _["Rm"] = Rm, _["bs"] = bs, _["ba"] = ba, _["b0"] = b0);
}

// One draw of (beta, a, b) from the joint full conditional.
// [[Rcpp::export]]
List replicated_gaussian_draw_cpp(const List& design, const arma::cube& ZT, const arma::mat& Sab,
		const arma::vec& rho, double s2, const arma::mat& P0, const arma::vec& m0, bool center = false) {
	const Design ds(design);
	const System sy = build_system(ds, ZT, Sab, rho, s2, P0, m0, center);
	const uword p = sy.p, k = sy.k;
	vec beta = sy.beta_mean;
	if (p) beta += sy.beta_root * vec(randn<vec>(p));
	mat F(sy.nr, k);
	if (k) for (uword h = 0; h < sy.nr; ++h) {
		mat Ri = block_inverse(sy, h);
		vec e = randn<vec>(k);
		F.row(h) = (Ri * (sy.D.col(h) - sy.C.slice(h) * beta + e)).t();
	}
	vec a, b;
	actor_effects(sy, F, a, b);
	if (!beta.is_finite() || !a.is_finite() || !b.is_finite())
		Rcpp::stop("Gaussian draw is not finite.");
	return List::create(_["beta"] = beta, _["a"] = a, _["b"] = b);
}

// Exact mean and covariance of (beta, a, b) under the same conditional.
// [[Rcpp::export(rng = false)]]
List replicated_gaussian_moments_cpp(const List& design, const arma::cube& ZT, const arma::mat& Sab,
		const arma::vec& rho, double s2, const arma::mat& P0, const arma::vec& m0, bool center = false) {
	const Design ds(design);
	const System sy = build_system(ds, ZT, Sab, rho, s2, P0, m0, center);
	const uword p = sy.p, k = sy.k, n = sy.n, nr = sy.nr, dim = p + nr * k;
	// root columns: p coefficient normals, then k normals per actor block
	mat root(p + 2 * n, dim, fill::zeros);
	vec mean(p + 2 * n, fill::zeros);
	if (p) {
		root.submat(0, 0, p - 1, p - 1) = sy.beta_root;
		mean.head(p) = sy.beta_mean;
	}
	mat Fmean(nr, k);
	std::vector<mat> Fcol(dim, mat(nr, k, fill::zeros));
	if (k) for (uword h = 0; h < nr; ++h) {
		mat Ri = block_inverse(sy, h);
		Fmean.row(h) = (Ri * (sy.D.col(h) - sy.C.slice(h) * sy.beta_mean)).t();
		mat coupling = -Ri * sy.C.slice(h) * sy.beta_root;
		for (uword j = 0; j < p; ++j) Fcol[j].row(h) = coupling.col(j).t();
		for (uword j = 0; j < k; ++j) Fcol[p + h * k + j].row(h) = Ri.col(j).t();
	}
	vec a, b;
	actor_effects(sy, Fmean, a, b);
	mean.subvec(p, p + n - 1) = a;
	mean.subvec(p + n, p + 2 * n - 1) = b;
	for (uword j = 0; j < dim; ++j) {
		actor_effects(sy, Fcol[j], a, b);
		root.submat(p, j, p + n - 1, j) = a;
		root.submat(p + n, j, p + 2 * n - 1, j) = b;
	}
	return List::create(_["mean"] = mean, _["covariance"] = root * root.t(), _["rank"] = k);
}
