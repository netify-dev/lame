// Kernels for the block coordinate descent of the fast (MCMC-free) AME
// estimator (.ae_bcd_fit in R/ame_als.R).
//
// A bcd fit holds its data (working response z, weights w, dyadic design x)
// fixed across sweeps, and every block of a sweep needs those arrays only
// through sums over time (n_row x n_col matrices) and over cells
// ((p + 1)-vectors). ae_bcd_stats_cpp forms those sums once per fit; the
// other kernels carry out the dyad-level arithmetic of one sweep in a single
// pass each, in place of the many whole-matrix operations the same
// arithmetic costs at R level. Column-major storage throughout.

#include <Rcpp.h>
#include <vector>
#include "ae_kernels.h"

using namespace Rcpp;

// sums of the fixed data of a bcd fit
//   sz  : sum_t w z                              (n_row x n_col)
//   gm  : sum_t w [1, x_1, ..., x_p]             (n_row * n_col x (p + 1))
//   xtz : sum over cells of w [1, x] z           (p + 1)
//   zz  : sum over cells of w z^2
//   xtx : sum over cells of w [1, x]' [1, x]     ((p + 1) x (p + 1))
//   gr, gc, szr, szc : row and column totals of gm and sz
// warr, zf : n_row x n_col x t arrays (zf zero where unobserved);
// xarr     : n_row x n_col x p x t design array. a cell of zero weight
// contributes nothing, whatever the design holds there.
// [[Rcpp::export(rng = false)]]
List ae_bcd_stats_cpp(const NumericVector& warr, const NumericVector& zf,
                      const NumericVector& xarr, int nr, int nc, int tt, int p, int threads = 1) {
	const std::size_t n_dyad = (std::size_t) nr * (std::size_t) nc;
	const int pa = p + 1;
	if ((std::size_t) warr.size() != n_dyad * tt || (std::size_t) zf.size() != n_dyad * tt ||
	    (std::size_t) xarr.size() != n_dyad * (std::size_t) p * tt) {
		stop("Incompatible dimensions for the bcd data sums.");
	}
	NumericMatrix sz(nr, nc);
	NumericMatrix gm((int) n_dyad, pa);
	NumericVector xtz(pa);
	NumericMatrix xtx(pa, pa);
	double zz = 0.0;
	const double* w = warr.begin();
	const double* z = zf.begin();
	const double* x = xarr.begin();
	double* szp = sz.begin();
	double* gmp = gm.begin();
	std::vector<double> xtz_acc(pa, 0.0), xtx_acc(pa * pa, 0.0);
	// the columns are independent: each yields its own partial sums, which are added in the order of the
	// columns afterwards, so that the result does not depend on the number of threads
	std::vector<double> colxtz((std::size_t) nc * pa, 0.0), colxtx((std::size_t) nc * pa * pa, 0.0),
		colzz((std::size_t) nc, 0.0);
	const int nt = lame_ae::resolve_threads(threads);
	(void) nt;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (int j = 0; j < nc; ++j) {
		double* col_xtz = &colxtz[(std::size_t) j * pa];
		double* col_xtx = &colxtx[(std::size_t) j * pa * pa];
		std::vector<double> xr(pa);
		double col_zz = 0.0;
		// the periods outermost, so that a period of a column is a run of contiguous cells in the weights, the
		// response and every plane of the design; a dyad still accumulates its periods in order
		for (int t = 0; t < tt; ++t) {
			const std::size_t base = (std::size_t) nr * j;
			for (int i = 0; i < nr; ++i) {
				const std::size_t dyad = (std::size_t) i + base;
				const std::size_t cell = dyad + n_dyad * (std::size_t) t;
				const double wc = w[cell];
				if (!(wc > 0.0)) continue;
				const double zc = z[cell];
				gmp[dyad] += wc;
				szp[dyad] += wc * zc;
				col_zz += wc * zc * zc;
				xr[0] = 1.0;
				for (int k = 0; k < p; ++k) {
					const double xv = x[dyad + n_dyad * ((std::size_t) k + (std::size_t) p * t)];
					xr[k + 1] = xv;
					gmp[dyad + n_dyad * (std::size_t) (k + 1)] += wc * xv;
				}
				for (int k = 0; k < pa; ++k) {
					const double wxk = wc * xr[k];
					col_xtz[k] += wxk * zc;
					for (int l = 0; l <= k; ++l) col_xtx[k + pa * l] += wxk * xr[l];
				}
			}
		}
		colzz[j] = col_zz;
	}
	for (int j = 0; j < nc; ++j) {
		const double* col_xtz = &colxtz[(std::size_t) j * pa];
		const double* col_xtx = &colxtx[(std::size_t) j * pa * pa];
		for (int k = 0; k < pa; ++k) {
			xtz_acc[k] += col_xtz[k];
			for (int l = 0; l <= k; ++l) xtx_acc[k + pa * l] += col_xtx[k + pa * l];
		}
		zz += colzz[j];
	}
	for (int k = 0; k < pa; ++k) {
		xtz[k] = xtz_acc[k];
		for (int l = 0; l <= k; ++l) {
			xtx(k, l) = xtx_acc[k + pa * l];
			xtx(l, k) = xtx_acc[k + pa * l];
		}
	}
	// row and column totals of gm and sz: the whole of the dyad-level data a
	// sweep without a multiplicative term needs
	NumericMatrix gr(nr, pa), gc(nc, pa);
	NumericVector szr(nr), szc(nc);
	for (int j = 0; j < nc; ++j) {
		const std::size_t off = (std::size_t) nr * j;
		for (int k = 0; k < pa; ++k) {
			const double* g = gmp + n_dyad * (std::size_t) k + off;
			double col = 0.0;
			for (int i = 0; i < nr; ++i) { gr(i, k) += g[i]; col += g[i]; }
			gc(j, k) = col;
		}
		double col = 0.0;
		for (int i = 0; i < nr; ++i) { szr[i] += szp[off + i]; col += szp[off + i]; }
		szc[j] = col;
	}
	return List::create(Named("sz") = sz, Named("gm") = gm, Named("xtz") = xtz,
	                    Named("zz") = zz, Named("xtx") = xtx, Named("gr") = gr,
	                    Named("gc") = gc, Named("szr") = szr, Named("szc") = szc);
}

// sum_t w [1, x_1, ..., x_p] for every dyad (n_row * n_col x (p + 1)), the gm of ae_bcd_stats_cpp for weights of
// either sign (the derivatives of a weight, which the laplace correction sums; ae_bcd_stats_cpp counts positive
// weights only). a cell of weight zero contributes nothing, whatever the design holds there
// [[Rcpp::export(rng = false)]]
NumericMatrix ae_weighted_sums_cpp(const NumericVector& warr, const NumericVector& xarr, int nr, int nc, int tt,
                                   int p, int threads = 1) {
	const std::size_t n_dyad = (std::size_t) nr * (std::size_t) nc;
	const int pa = p + 1;
	if ((std::size_t) warr.size() != n_dyad * tt || (std::size_t) xarr.size() != n_dyad * (std::size_t) p * tt) {
		stop("Incompatible dimensions for the weighted sums.");
	}
	NumericMatrix gm((int) n_dyad, pa);
	const double* w = warr.begin();
	const double* x = xarr.begin();
	double* gmp = gm.begin();
	const int nt = lame_ae::resolve_threads(threads);
	(void) nt;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nt) schedule(static) if (nt > 1)
#endif
	for (int j = 0; j < nc; ++j) {
		for (int i = 0; i < nr; ++i) {
			const std::size_t dyad = (std::size_t) i + (std::size_t) nr * j;
			double s_w = 0.0;
			for (int t = 0; t < tt; ++t) {
				const double wc = w[dyad + n_dyad * (std::size_t) t];
				if (wc == 0.0) continue;
				s_w += wc;
				for (int k = 0; k < p; ++k) {
					gmp[dyad + n_dyad * (std::size_t) (k + 1)] += wc * x[dyad + n_dyad * ((std::size_t) k + (std::size_t) p * t)];
				}
			}
			gmp[dyad] = s_w;
		}
	}
	return gm;
}

// gm' vec(c) with c_ij = a_i + b_j + omat_ij: the dyad-level part of the
// right-hand side of the regression block's normal equations
// [[Rcpp::export(rng = false)]]
NumericVector ae_bcd_gtc_cpp(const NumericVector& a, const NumericVector& b,
                             const NumericMatrix& omat, const NumericMatrix& gm) {
	const int nr = omat.nrow(), nc = omat.ncol(), pa = gm.ncol();
	if (a.size() != nr || b.size() != nc || gm.nrow() != nr * nc) {
		stop("Incompatible dimensions for the bcd regression right-hand side.");
	}
	const std::size_t n_dyad = (std::size_t) nr * (std::size_t) nc;
	NumericVector out(pa);
	std::vector<double> acc(pa, 0.0), cbuf(nr);
	const double* om = omat.begin();
	const double* gp = gm.begin();
	for (int j = 0; j < nc; ++j) {
		const double bj = b[j];
		const double* omj = om + (std::size_t) nr * j;
		for (int i = 0; i < nr; ++i) cbuf[i] = a[i] + bj + omj[i];
		for (int k = 0; k < pa; ++k) {
			const double* g = gp + n_dyad * (std::size_t) k + (std::size_t) nr * j;
			double s = 0.0;
			for (int i = 0; i < nr; ++i) s += g[i] * cbuf[i];
			acc[k] += s;
		}
	}
	for (int k = 0; k < pa; ++k) out[k] = acc[k];
	return out;
}

// the residual sums the additive block needs, given the regression
// coefficients theta = [mu, beta]:
//   phi_ij  = sum_k gm_ijk theta_k            (= sum_t w (mu + x'beta))
//   g_ij    = sz_ij - phi_ij - wsum_ij omat_ij
//   base_rs = row sums of g, base_cs = column sums of g
// [[Rcpp::export(rng = false)]]
List ae_bcd_rs_cpp(const NumericVector& theta, const NumericMatrix& gm,
                   const NumericMatrix& sz, const NumericMatrix& wsum,
                   const NumericMatrix& omat) {
	const int nr = sz.nrow(), nc = sz.ncol(), pa = gm.ncol();
	if (theta.size() != pa || gm.nrow() != nr * nc || wsum.nrow() != nr ||
	    wsum.ncol() != nc || omat.nrow() != nr || omat.ncol() != nc) {
		stop("Incompatible dimensions for the bcd additive right-hand side.");
	}
	const std::size_t n_dyad = (std::size_t) nr * (std::size_t) nc;
	NumericMatrix phi(nr, nc);
	NumericVector rs(nr), cs(nc);
	const double* gp = gm.begin();
	const double* szp = sz.begin();
	const double* wp = wsum.begin();
	const double* om = omat.begin();
	double* ph = phi.begin();
	std::vector<double> rs_acc(nr, 0.0);
	for (int j = 0; j < nc; ++j) {
		const std::size_t off = (std::size_t) nr * j;
		for (int i = 0; i < nr; ++i) ph[off + i] = gp[off + i] * theta[0];
		for (int k = 1; k < pa; ++k) {
			const double tk = theta[k];
			const double* g = gp + n_dyad * (std::size_t) k + off;
			for (int i = 0; i < nr; ++i) ph[off + i] += g[i] * tk;
		}
		double cs_j = 0.0;
		for (int i = 0; i < nr; ++i) {
			const double v = szp[off + i] - ph[off + i] - wp[off + i] * om[off + i];
			rs_acc[i] += v;
			cs_j += v;
		}
		cs[j] = cs_j;
	}
	for (int i = 0; i < nr; ++i) rs[i] = rs_acc[i];
	return List::create(Named("base_rs") = rs, Named("base_cs") = cs, Named("phi") = phi);
}

// the weighted per-dyad mean residual that the multiplicative block fits,
//   rbar_ij = (sz_ij - phi_ij - dmu wsum_ij - wsum_ij (a_i + b_j)) / wsum_ij
// on observed dyads (wsum > 0), zero elsewhere; dmu is the change in the
// intercept since phi was formed (the grand mean folded in by centring)
// [[Rcpp::export(rng = false)]]
NumericMatrix ae_bcd_rbar_cpp(const NumericVector& a, const NumericVector& b,
                              const NumericMatrix& phi, double dmu,
                              const NumericMatrix& sz, const NumericMatrix& wsum) {
	const int nr = sz.nrow(), nc = sz.ncol();
	if (a.size() != nr || b.size() != nc || phi.nrow() != nr || phi.ncol() != nc ||
	    wsum.nrow() != nr || wsum.ncol() != nc) {
		stop("Incompatible dimensions for the bcd multiplicative target.");
	}
	NumericMatrix out(nr, nc);
	const double* ph = phi.begin();
	const double* szp = sz.begin();
	const double* wp = wsum.begin();
	double* o = out.begin();
	for (int j = 0; j < nc; ++j) {
		const std::size_t off = (std::size_t) nr * j;
		const double bj = b[j];
		for (int i = 0; i < nr; ++i) {
			const double wc = wp[off + i];
			if (wc > 0.0) {
				const double res = szp[off + i] - (ph[off + i] + dmu * wc) - wc * (a[i] + bj);
				o[off + i] = res / wc;
			}
		}
	}
	return out;
}

// the two dyad-level sums of the residual sum of squares of a whole iterate,
//   s1 = sum c_ij u_ij ,  s2 = sum wsum_ij c_ij^2 ,
// with c_ij = a_i + b_j + omat_ij and u_ij = sz_ij - phi_ij - dmu wsum_ij
// [[Rcpp::export(rng = false)]]
NumericVector ae_bcd_sse_cpp(const NumericVector& a, const NumericVector& b,
                             const NumericMatrix& omat, const NumericMatrix& phi,
                             double dmu, const NumericMatrix& sz,
                             const NumericMatrix& wsum) {
	const int nr = sz.nrow(), nc = sz.ncol();
	if (a.size() != nr || b.size() != nc || omat.nrow() != nr || omat.ncol() != nc ||
	    phi.nrow() != nr || phi.ncol() != nc || wsum.nrow() != nr || wsum.ncol() != nc) {
		stop("Incompatible dimensions for the bcd residual sum of squares.");
	}
	const double* om = omat.begin();
	const double* ph = phi.begin();
	const double* szp = sz.begin();
	const double* wp = wsum.begin();
	double s1 = 0.0, s2 = 0.0;
	for (int j = 0; j < nc; ++j) {
		const std::size_t off = (std::size_t) nr * j;
		const double bj = b[j];
		double c1 = 0.0, c2 = 0.0;
		for (int i = 0; i < nr; ++i) {
			const double c = a[i] + bj + om[off + i];
			const double wc = wp[off + i];
			c1 += c * (szp[off + i] - (ph[off + i] + dmu * wc));
			c2 += wc * c * c;
		}
		s1 += c1;
		s2 += c2;
	}
	return NumericVector::create(s1, s2);
}

// the additive block of a sweep: the sender / receiver effects that minimise
// the weighted residual given the rest, by the gauss-seidel recursion of the
// r code (up to 200 sweeps to a 1e-10 step), started from (a0, b0) -- or, when
// `direct` is supplied, from the exact solution of the normal equations
//   cnt_row a_i + sum_j wsum_ij b_j = base_rs_i ,
//   sum_i wsum_ij a_i + cnt_col b_j = base_cs_j
// computed from the cholesky factor in `direct` (see .ae_additive_factor) and
// moved along the a + t, b - t gauge to the count-weighted sum of the column
// effects that the recursion conserves. returns the effects and the largest
// residual of the normal equations with the scale it is judged against.
// symmetric models carry one effect vector (b = a) and the recursion
// a = (base_rs - wsum a) / cnt.
static void tri_solve_upper_system(const double* ch, int n, std::vector<double>& x) {
	// solve ch' ch x = rhs in place, ch upper triangular (column-major)
	for (int i = 0; i < n; ++i) {
		const double* col = ch + (std::size_t) n * i;
		double s = x[i];
		for (int k = 0; k < i; ++k) s -= col[k] * x[k];
		x[i] = s / col[i];
	}
	for (int i = n - 1; i >= 0; --i) {
		const double* col = ch + (std::size_t) n * i;
		x[i] /= col[i];
		const double xi = x[i];
		for (int k = 0; k < i; ++k) x[k] -= col[k] * xi;
	}
}

// y = wsum b (axpy form) or y = wsum' a (dot form)
static void matvec_w(const double* w, int nr, int nc, const double* b, double* y) {
	std::fill(y, y + nr, 0.0);
	for (int j = 0; j < nc; ++j) {
		const double bj = b[j];
		const double* wj = w + (std::size_t) nr * j;
		for (int i = 0; i < nr; ++i) y[i] += wj[i] * bj;
	}
}
static void matvec_wt(const double* w, int nr, int nc, const double* a, double* y) {
	for (int j = 0; j < nc; ++j) {
		const double* wj = w + (std::size_t) nr * j;
		double s = 0.0;
		for (int i = 0; i < nr; ++i) s += wj[i] * a[i];
		y[j] = s;
	}
}

// largest residual of the normal equations (zero for actors with no weight)
static double additive_resid(const double* w, int nr, int nc, const double* base_rs,
	const double* base_cs, const double* a, const double* b, const double* cnt_row,
	const double* cnt_col, std::vector<double>& tmp_r, std::vector<double>& tmp_c) {
	matvec_w(w, nr, nc, b, tmp_r.data());
	matvec_wt(w, nr, nc, a, tmp_c.data());
	double worst = 0.0;
	for (int i = 0; i < nr; ++i) {
		if (cnt_row[i] == 0.0) continue;
		const double r = std::fabs(cnt_row[i] * a[i] + tmp_r[i] - base_rs[i]);
		if (std::isnan(r)) return r;
		if (r > worst) worst = r;
	}
	for (int j = 0; j < nc; ++j) {
		if (cnt_col[j] == 0.0) continue;
		const double r = std::fabs(cnt_col[j] * b[j] + tmp_c[j] - base_cs[j]);
		if (std::isnan(r)) return r;
		if (r > worst) worst = r;
	}
	return worst;
}

// [[Rcpp::export(rng = false)]]
List ae_additive_cpp(const NumericMatrix& wsum, const NumericVector& base_rs,
                     const NumericVector& base_cs, const NumericVector& a0,
                     const NumericVector& b0, const NumericVector& cnt_row,
                     const NumericVector& cnt_col, bool symmetric, SEXP direct) {
	const int nr = wsum.nrow(), nc = wsum.ncol();
	if (base_rs.size() != nr || base_cs.size() != nc || a0.size() != nr || b0.size() != nc ||
	    cnt_row.size() != nr || cnt_col.size() != nc) {
		stop("Incompatible dimensions for the bcd additive block.");
	}
	std::vector<double> a(a0.begin(), a0.end()), b(b0.begin(), b0.end());
	std::vector<double> pos_r(nr), pos_c(nc);
	for (int i = 0; i < nr; ++i) pos_r[i] = std::max(cnt_row[i], 1.0);
	for (int j = 0; j < nc; ++j) pos_c[j] = std::max(cnt_col[j], 1.0);
	const double* w = wsum.begin();

	if (!Rf_isNull(direct)) {
		List d(direct);
		NumericMatrix ch = d["ch"];
		IntegerVector rows = d["rows"];
		const int n = ch.nrow();
		std::vector<double> x(n);
		if (symmetric) {
			for (int k = 0; k < n; ++k) x[k] = base_rs[rows[k] - 1];
			tri_solve_upper_system(ch.begin(), n, x);
			std::fill(a.begin(), a.end(), 0.0);
			for (int k = 0; k < n; ++k) a[rows[k] - 1] = x[k];
			b = a;
		} else {
			IntegerVector free_cols = d["free"];
			IntegerVector cols = d["cols"];
			NumericVector cnt_act = d["cnt_act"];
			const int nrow_act = rows.size();
			for (int k = 0; k < nrow_act; ++k) x[k] = base_rs[rows[k] - 1];
			for (int k = 0; k < free_cols.size(); ++k) x[nrow_act + k] = base_cs[free_cols[k] - 1];
			tri_solve_upper_system(ch.begin(), n, x);
			std::fill(a.begin(), a.end(), 0.0);
			std::fill(b.begin(), b.end(), 0.0);
			for (int k = 0; k < nrow_act; ++k) a[rows[k] - 1] = x[k];
			for (int k = 0; k < free_cols.size(); ++k) b[free_cols[k] - 1] = x[nrow_act + k];
			double num = 0.0, den = 0.0, start = 0.0;
			for (int k = 0; k < cols.size(); ++k) { num += cnt_act[k] * b[cols[k] - 1]; den += cnt_act[k]; }
			for (int j = 0; j < nc; ++j) start += cnt_col[j] * b0[j];
			const double shift = (num - start) / den;
			for (int k = 0; k < nrow_act; ++k) a[rows[k] - 1] += shift;
			for (int k = 0; k < cols.size(); ++k) b[cols[k] - 1] -= shift;
		}
	}

	std::vector<double> tmp_r(nr), tmp_c(nc), a_prev(nr), b_prev(nc);
	for (int sweep = 0; sweep < 200; ++sweep) {
		a_prev = a; b_prev = b;
		if (symmetric) {
			matvec_w(w, nr, nc, a.data(), tmp_r.data());
			for (int i = 0; i < nr; ++i) a[i] = (base_rs[i] - tmp_r[i]) / pos_r[i];
			for (int i = 0; i < nr; ++i) if (cnt_row[i] == 0.0) a[i] = 0.0;
			b = a;
		} else {
			matvec_w(w, nr, nc, b.data(), tmp_r.data());
			for (int i = 0; i < nr; ++i) a[i] = (base_rs[i] - tmp_r[i]) / pos_r[i];
			for (int i = 0; i < nr; ++i) if (cnt_row[i] == 0.0) a[i] = 0.0;
			matvec_wt(w, nr, nc, a.data(), tmp_c.data());
			for (int j = 0; j < nc; ++j) b[j] = (base_cs[j] - tmp_c[j]) / pos_c[j];
			for (int j = 0; j < nc; ++j) if (cnt_col[j] == 0.0) b[j] = 0.0;
		}
		double da = 0.0, db = 0.0;
		for (int i = 0; i < nr; ++i) da = std::max(da, std::fabs(a[i] - a_prev[i]));
		for (int j = 0; j < nc; ++j) db = std::max(db, std::fabs(b[j] - b_prev[j]));
		if (da < 1e-10 && db < 1e-10) break;
	}
	const double resid = additive_resid(w, nr, nc, base_rs.begin(), base_cs.begin(), a.data(),
		b.data(), cnt_row.begin(), cnt_col.begin(), tmp_r, tmp_c);
	double scale = 1.0;
	for (int i = 0; i < nr; ++i) scale = std::max(scale, std::fabs(base_rs[i]));
	for (int j = 0; j < nc; ++j) scale = std::max(scale, std::fabs(base_cs[j]));
	return List::create(Named("a") = NumericVector(a.begin(), a.end()),
	                    Named("b") = NumericVector(b.begin(), b.end()),
	                    Named("resid") = resid, Named("scale") = scale);
}

// the residual of the additive normal equations at given effects (the
// certificate that follows an exact re-solve)
// [[Rcpp::export(rng = false)]]
double ae_additive_resid_cpp(const NumericMatrix& wsum, const NumericVector& base_rs,
                             const NumericVector& base_cs, const NumericVector& a,
                             const NumericVector& b, const NumericVector& cnt_row,
                             const NumericVector& cnt_col) {
	const int nr = wsum.nrow(), nc = wsum.ncol();
	std::vector<double> tmp_r(nr), tmp_c(nc);
	return additive_resid(wsum.begin(), nr, nc, base_rs.begin(), base_cs.begin(), a.begin(),
		b.begin(), cnt_row.begin(), cnt_col.begin(), tmp_r, tmp_c);
}

// weighted residual sum of squares of a whole iterate straight from the cells,
//   sum_ijt w_ijt (z_ijt - mu - x_ijt' beta - c_ij)^2 ,
// with c the n_row x n_col matrix of effects (a_i + b_j + o_ij). the reference
// form of the sums-based expression the sweeps use; periods are summed in
// long double and added in order
// [[Rcpp::export(rng = false)]]
double ae_bcd_sse_exact_cpp(const NumericVector& warr, const NumericVector& zf,
                            const NumericVector& xarr, int nr, int nc, int tt, int p,
                            double mu, const NumericVector& beta, const NumericMatrix& cmat) {
	const std::size_t n_dyad = (std::size_t) nr * (std::size_t) nc;
	if ((std::size_t) warr.size() != n_dyad * tt || (std::size_t) zf.size() != n_dyad * tt ||
	    (std::size_t) xarr.size() != n_dyad * (std::size_t) p * tt || beta.size() != p ||
	    (std::size_t) cmat.size() != n_dyad) {
		stop("Incompatible dimensions for the bcd residual sum of squares.");
	}
	const double* w = warr.begin();
	const double* z = zf.begin();
	const double* x = xarr.begin();
	const double* c = cmat.begin();
	double total = 0.0;
	for (int t = 0; t < tt; ++t) {
		long double period = 0.0L;
		const std::size_t cell0 = n_dyad * (std::size_t) t;
		for (std::size_t dyad = 0; dyad < n_dyad; ++dyad) {
			const double wc = w[cell0 + dyad];
			double eta = mu + c[dyad];
			for (int k = 0; k < p; ++k) {
				eta += beta[k] * x[dyad + n_dyad * ((std::size_t) k + (std::size_t) p * t)];
			}
			const double r = z[cell0 + dyad] - eta;
			period += (long double) (wc * (r * r));
		}
		total += (double) period;
	}
	return total;
}
