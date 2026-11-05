#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
static mat upper_solve393(const mat& R,const mat& y){
 mat x;if(!solve(x,trimatu(R),y,solve_opts::fast+solve_opts::no_approx)||!x.is_finite())stop("Joint Gaussian triangular solve unresolved.");return x;
}
// Preparation keeps compact roots and the thin orthogonal response maps. It is
// performed once for a design, independently of rho, Sab and random draws.
// [[Rcpp::export(rng = false)]]
List joint_levels_context_cpp(const arma::mat& D,const arma::mat& H,const arma::mat& prior,const arma::uvec& from,
 const arma::uvec& to,const arma::uvec& partner){
 const uword n=H.n_rows,k=H.n_cols,p=prior.n_rows,q=p+2*k,m=D.n_rows;
 if(n<2||p<1||m==0||k!=n-1||prior.n_cols!=p||D.n_cols!=q||from.n_elem!=m||to.n_elem!=m||partner.n_elem!=m||!D.is_finite()||!H.is_finite()||!prior.is_finite())stop("Invalid joint Gaussian design.");
 if(any(from>=n)||any(to>=n)||any(partner>=m))stop("Invalid joint Gaussian dyad indices.");
 // The compact elimination uses the complete balanced dyad grid and an
 // orthonormal centered actor basis, so verify these context invariants once.
 const double basis_tol=128*std::numeric_limits<double>::epsilon()*n;
 if(max(max(abs(H.t()*H-eye(k,k))))>basis_tol||max(abs(sum(H,0)))>basis_tol)stop("Joint actor basis must be centered and orthonormal.");
 if(any(vectorise(prior!=prior.t())))stop("Joint coefficient prior must be symmetric.");
 if(m%(n*(n-1))!=0)stop("Joint Gaussian design must contain complete periods.");
 const uword repetitions=m/(n*(n-1));umat counts(n,n,fill::zeros);
 for(uword i=0;i<m;++i)++counts(from(i),to(i));
 for(uword i=0;i<n;++i)for(uword j=0;j<n;++j)if(counts(i,j)!=(i==j?0:repetitions))stop("Joint Gaussian design must contain balanced complete periods.");
 for(uword i=0;i<m;++i)if(partner(partner(i))!=i||from(i)!=to(partner(i))||to(i)!=from(partner(i))||from(i)==to(i))stop("Invalid joint Gaussian partner mapping.");
 for(uword i=0;i<m;++i)for(uword j=0;j<k;++j)if(D(i,p+j)!=H(from(i),j)||D(i,p+k+j)!=H(to(i),j))stop("Joint actor design disagrees with the supplied basis.");
 uvec pair=find(regspace<uvec>(0,m-1)<partner);
 const uword nt=m/(n*(n-1));if(nt*n*(n-1)!=m)stop("Joint Gaussian design must contain complete periods.");
 uvec fi=from.elem(pair),tj=to.elem(pair);mat X=D.cols(0,p-1);
 mat plus=(X.rows(pair)+X.rows(partner.elem(pair)))/std::sqrt(2.0);
 mat minus=(X.rows(pair)-X.rows(partner.elem(pair)))/std::sqrt(2.0);
 mat margins_plus(n,p,fill::zeros),margins_minus(n,p,fill::zeros);
 for(uword i=0;i<pair.n_elem;++i){margins_plus.row(fi(i))+=plus.row(i);margins_plus.row(tj(i))+=plus.row(i);margins_minus.row(fi(i))+=minus.row(i);margins_minus.row(tj(i))-=minus.row(i);}
 const double sp=std::sqrt(double(nt)*(n-2)),sm=std::sqrt(double(nt)*n);
 mat projected_plus(k,p,fill::zeros),projected_minus=H.t()*margins_minus/sm;
 if(n>2)projected_plus=H.t()*margins_plus/sp;
 mat node_plus=H*projected_plus,node_minus=H*projected_minus;
 for(uword i=0;i<pair.n_elem;++i){if(n>2)plus.row(i)-=(node_plus.row(fi(i))+node_plus.row(tj(i)))/sp;minus.row(i)-=(node_minus.row(fi(i))-node_minus.row(tj(i)))/sm;}
 mat Qp,Rp,Qm,Rm;if(!qr_econ(Qp,Rp,plus)||!qr_econ(Qm,Rm,minus))stop("Joint Gaussian residual design QR failed.");
 mat prior_beta(p,p);if(!chol(prior_beta,prior))stop("Joint coefficient prior must be positive definite.");
 return List::create(_["X"]=X,_["H"]=H,_["prior"]=prior,_["prior_beta"]=prior_beta,_["from"]=from,_["to"]=to,_["partner"]=partner,_["pair"]=pair,_["projected_plus"]=projected_plus,_["projected_minus"]=projected_minus,_["Qplus"]=Qp,_["Rplus"]=Rp,_["Qminus"]=Qm,_["Rminus"]=Rm,_["nt"]=nt);
}
// [[Rcpp::export(rng = false)]]
List joint_levels_parameters_cpp(const List& cx,const arma::vec& z,const arma::mat& Sab,double rho){
 mat H=cx["H"],X=cx["X"],Pbeta=cx["prior_beta"],Qp=cx["Qplus"],Qm=cx["Qminus"],Rp=cx["Rplus"],Rm=cx["Rminus"],Cp=cx["projected_plus"],Cm=cx["projected_minus"];
 uvec pair=cx["pair"],partner=cx["partner"],from=cx["from"],to=cx["to"];uword n=H.n_rows,p=Pbeta.n_rows,k=H.n_cols,q=p+2*k,nt=as<uword>(cx["nt"]);
 if(z.n_elem!=X.n_rows||!z.is_finite()||!std::isfinite(rho)||std::abs(rho)>=1||Sab.n_rows!=2||Sab.n_cols!=2||!Sab.is_finite())stop("Invalid joint Gaussian response or covariance.");
 mat G;if(!chol(G,mat(0.5*(Sab+Sab.t())),"lower"))stop("Joint actor covariance must be positive definite.");
 const double wp=1/std::sqrt(1+rho),wm=1/std::sqrt(1-rho),sp=std::sqrt(double(nt)*(n-2)),sm=std::sqrt(double(nt)*n);
 uvec fi=from.elem(pair),tj=to.elem(pair);
 vec zp=(z.elem(pair)+z.elem(partner.elem(pair)))/std::sqrt(2.0),zm=(z.elem(pair)-z.elem(partner.elem(pair)))/std::sqrt(2.0);
 vec mp(n,fill::zeros),mm(n,fill::zeros);
 for(uword i=0;i<pair.n_elem;++i){mp(fi(i))+=zp(i);mp(tj(i))+=zp(i);mm(fi(i))+=zm(i);mm(tj(i))-=zm(i);}
 vec dp(k,fill::zeros),dm=H.t()*mm/sm;if(n>2)dp=H.t()*mp/sp;
 vec hp=H*dp,hm=H*dm;
 for(uword i=0;i<pair.n_elem;++i){if(n>2)zp(i)-=(hp(fi(i))+hp(tj(i)))/sp;zm(i)-=(hm(fi(i))-hm(tj(i)))/sm;}
 mat E(4,2,fill::zeros);E.row(0)=(sp*wp/std::sqrt(2.0))*(G.row(0)+G.row(1));E.row(1)=(sm*wm/std::sqrt(2.0))*(G.row(0)-G.row(1));E.rows(2,3)=eye(2,2);
 mat Qa,Ra0;if(!qr(Qa,Ra0,E))stop("Joint Gaussian actor QR failed.");mat Ra=Ra0.rows(0,1);mat Ri=upper_solve393(Ra,eye(2,2));
 mat A(2*k+Rp.n_rows+Rm.n_rows+p,p,fill::zeros);vec y(A.n_rows,fill::zeros);std::vector<mat> blocks(k);std::vector<vec> responses(k);
 for(uword h=0;h<k;++h){mat C(4,p,fill::zeros);C.row(0)=wp*Cp.row(h);C.row(1)=wm*Cm.row(h);C=Qa.t()*C;vec v={wp*dp(h),wm*dm(h),0,0};v=Qa.t()*v;blocks[h]=C.rows(0,1);responses[h]=v.head(2);A.rows(2*h,2*h+1)=C.rows(2,3);y.subvec(2*h,2*h+1)=v.tail(2);}
 uword at=2*k;A.rows(at,at+Rp.n_rows-1)=wp*Rp;y.subvec(at,at+Rp.n_rows-1)=wp*(Qp.t()*zp);at+=Rp.n_rows;
 A.rows(at,at+Rm.n_rows-1)=wm*Rm;y.subvec(at,at+Rm.n_rows-1)=wm*(Qm.t()*zm);at+=Rm.n_rows;A.rows(at,at+p-1)=Pbeta;
 vec scales(p);for(uword j=0;j<p;++j){double scale=0;for(uword i=0;i<A.n_rows;++i)scale=std::hypot(scale,A(i,j));if(!(scale>0)||!std::isfinite(scale))stop("Joint Gaussian column scale unresolved.");scales(j)=scale;A.col(j)/=scale;}
 mat Q,R;if(!qr_econ(Q,R,A))stop("Joint Gaussian coefficient QR failed.");
 vec beta=upper_solve393(R,Q.t()*y);beta/=scales;mat Lb=upper_solve393(R,eye(p,p));Lb.each_col()/=scales;
 vec mean(q,fill::zeros);mean.head(p)=beta;mat noise(q,q,fill::zeros);noise.submat(0,0,p-1,p-1)=Lb;
 for(uword h=0;h<k;++h){uvec indices={p+h,p+k+h};mat conditional=G*Ri;mean.elem(indices)=conditional*(responses[h]-blocks[h]*beta);noise.submat(indices,regspace<uvec>(0,p-1))=-conditional*blocks[h]*Lb;uvec nc={p+2*h,p+2*h+1};noise.submat(indices,nc)=conditional;}

 return List::create(_["mean"]=mean,_["noise_root"]=noise);
}
