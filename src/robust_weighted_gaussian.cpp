// Weighted Gaussian conditional kernels for persistent pair variance.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

static mat ar_precision(int T, double rho, double variance) {
  mat Q(T,T,fill::zeros);
  Q(0,0) = 1-rho*rho;
  for(int t=1;t<T;++t) {
    Q(t,t) += 1; Q(t-1,t-1) += rho*rho;
    Q(t,t-1) -= rho; Q(t-1,t) -= rho;
  }
  return Q/variance;
}
static mat actor_precision(int T, double rho, double variance, const vec& tau) {
  if (tau.n_elem != (uword)T || !tau.is_finite() || any(tau <= 0)) stop("Invalid actor precision path");
  mat Q(T,T,fill::zeros);Q(0,0)=tau[0]*(1-rho*rho);
  for(int t=1;t<T;++t){Q(t,t)+=tau[t];Q(t-1,t-1)+=tau[t]*rho*rho;Q(t,t-1)-=tau[t]*rho;Q(t-1,t)-=tau[t]*rho;}
  return Q/variance;
}
static vec normal_draw(const mat& Q, const vec& h) {
  mat R=chol(Q); vec z(h.n_elem);
  for(uword k=0;k<z.n_elem;++k) z[k]=R::rnorm(0,1);
  vec m=solve(trimatu(R),solve(trimatl(R.t()),h));
  return m+solve(trimatu(R),z);
}
// [[Rcpp::export]]
arma::vec robust_weighted_beta(const arma::cube& residual, const arma::cube& weights,
    const arma::mat& design, double variance, const arma::vec& prior_precision) {
  int p=design.n_cols; mat Q=diagmat(prior_precision); vec h(p,fill::zeros);
  vec r=vectorise(residual),w=vectorise(weights)/variance;
  for(uword k=0;k<w.n_elem;++k) if(w[k]>0) {
    vec x=design.row(k).t(); Q+=w[k]*(x*x.t());h+=w[k]*r[k]*x;
  }
  return normal_draw(Q,h);
}
// [[Rcpp::export]]
arma::vec robust_weighted_mu(const arma::cube& residual,const arma::cube& weights,
    double variance,double rho,double innovation_variance) {
  int T=residual.n_slices;mat Q=ar_precision(T,rho,innovation_variance);vec h(T,fill::zeros);
  for(int t=0;t<T;++t) for(uword k=0;k<residual.slice(t).n_elem;++k) {
    double w=weights.slice(t)[k]/variance;
    if(w>0){Q(t,t)+=w;h[t]+=w*residual.slice(t)[k];}
  }
  return normal_draw(Q,h);
}
// Joint actor trajectory and actor level, with exact zero-sum constraints.
// [[Rcpp::export]]
Rcpp::List robust_weighted_actor(const arma::cube& residual,const arma::cube& weights,
    double variance,double rho,double innovation_variance,const arma::mat& level_cov,
    const arma::vec& other_level,bool row_side,
    Rcpp::Nullable<Rcpp::NumericVector> innovation_precision = R_NilValue) {
  int n=row_side?residual.n_rows:residual.n_cols,T=residual.n_slices,q=T+1;
  mat A = innovation_precision.isNull() ? ar_precision(T,rho,innovation_variance) :
    actor_precision(T,rho,innovation_variance,Rcpp::as<arma::vec>(innovation_precision.get()));
  mat P=inv_sympd(level_cov);
  vec onesT(T,fill::ones),q1=A*onesT;mat base(q,q,fill::zeros);
  base.submat(0,0,T-1,T-1)=A;base.submat(0,T,T-1,T)=-q1;
  base.submat(T,0,T,T-1)=-q1.t();int side=row_side?0:1;
  base(T,T)=accu(q1)+P(side,side);
  std::vector<mat> S(n);mat sumS(q,q,fill::zeros),draws(n,q);
  for(int i=0;i<n;++i) {
    mat Q=base;vec h(q,fill::zeros);h[T]=-P(side,1-side)*other_level[i];
    for(int t=0;t<T;++t) {
      int m=row_side?residual.n_cols:residual.n_rows;
      for(int j=0;j<m;++j) {
        double w=(row_side?weights(i,j,t):weights(j,i,t))/variance;
        if(w>0){Q(t,t)+=w;h[t]+=w*(row_side?residual(i,j,t):residual(j,i,t));}
      }
    }
    S[i]=inv_sympd(Q);sumS+=S[i];draws.row(i)=normal_draw(Q,h).t();
  }
  vec shift=solve(sumS,sum(draws,0).t());
  for(int i=0;i<n;++i)draws.row(i)-=(S[i]*shift).t();
  return List::create(_["path"]=draws.cols(0,T-1),_["level"]=draws.col(T));
}
// Independent observation errors with heterogeneous precision. Each actor's
// entire rank-by-time path is drawn jointly; no artificial diagonal is used.
// [[Rcpp::export]]
Rcpp::List robust_weighted_uv(const arma::cube& residual,const arma::cube& weights,
    const arma::cube& U_current,const arma::cube& V_current,double variance,
    double rho,double innovation_variance) {
  cube U=U_current,V=V_current;int T=U.n_slices,R=U.n_cols;
  mat base=kron(ar_precision(T,rho,innovation_variance),eye(R,R));
  for(int side=0;side<2;++side) {
    cube& target=side==0?U:V;const cube& partner=side==0?V:U;
    for(uword i=0;i<target.n_rows;++i) {
      mat Q=base;vec h(T*R,fill::zeros);
      for(int t=0;t<T;++t) for(uword j=0;j<partner.n_rows;++j) {
        double w=(side==0?weights(i,j,t):weights(j,i,t))/variance;
        if(w>0) {
          vec x=partner.slice(t).row(j).t();double y=side==0?residual(i,j,t):residual(j,i,t);
          Q.submat(t*R,t*R,(t+1)*R-1,(t+1)*R-1)+=w*(x*x.t());
          h.subvec(t*R,(t+1)*R-1)+=w*y*x;
        }
      }
      vec z=normal_draw(Q,h);
      for(int t=0;t<T;++t)target.slice(t).row(i)=z.subvec(t*R,(t+1)*R-1).t();
    }
  }
  return List::create(_["U"]=U,_["V"]=V);
}
