// Exact observed-probit Langevin transitions, prior refresh, sign reflections, and replica exchange.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
struct TemperedBinaryContext {
 int n,r,tt,p,base;
 arma::mat H,B,X,counts,hss,incidence_sum;
 arma::vec beta_precision;
 std::vector<arma::mat> incidence;
 arma::uvec ii,jj,time;
 arma::vec sy;
 double av,bv,ash,bs,rm,rs;
 bool hierarchical;
 bool collapsed;
 TemperedBinaryContext(List x):n(as<int>(x["n"])),r(as<int>(x["r"])),tt(as<int>(x["tt"])),p(as<int>(x["p"])),H(as<arma::mat>(x["H"])),B(as<arma::mat>(x["B"])),X(as<arma::mat>(x["X"])),ii(as<arma::uvec>(x["ii"])),jj(as<arma::uvec>(x["jj"])),time(as<arma::uvec>(x["time"])),sy(as<arma::vec>(x["sy"])),av(as<double>(x["av"])),bv(as<double>(x["bv"])),ash(as<double>(x["as"])),bs(as<double>(x["bs"])),rm(as<double>(x["rm"])),rs(as<double>(x["rs"])),hierarchical(as<bool>(x["hierarchical"])),collapsed(as<bool>(x["collapsed"])) {
 if(n<4||r<1||r>=n||tt<2||p<1||hierarchical||
    H.n_rows!=(unsigned)n||H.n_cols!=(unsigned)(n-1)||B.n_rows!=(unsigned)p||B.n_cols!=(unsigned)p||
    X.n_rows!=sy.n_elem||X.n_cols!=(unsigned)p||ii.n_elem!=sy.n_elem||jj.n_elem!=sy.n_elem||time.n_elem!=sy.n_elem||
    !H.is_finite()||!B.is_finite()||!X.is_finite()||!sy.is_finite()||
    !R_finite(av)||av<=0||!R_finite(bv)||bv<=0||!R_finite(ash)||ash<=0||!R_finite(bs)||bs<=0||
    !R_finite(rm)||std::abs(rm)>=.99||!R_finite(rs)||rs<=0)
   stop("Invalid tempered binary context.");
 for(unsigned z=0;z<sy.n_elem;z++) if(ii[z]>=jj[z]||jj[z]>=(unsigned)n||time[z]>=(unsigned)tt||std::abs(sy[z])!=1)
   stop("Invalid observed unordered dyads in tempered context.");
 double width=(double)n*r+r+p+(double)(n-1)*(tt+1);
 if(width>2147483644.0)stop("Tempered state exceeds the integer dimension limit.");
 base=(int)width;counts.zeros(n,n);hss.zeros(n-1,tt);incidence_sum.zeros(n,n);for(int t=0;t<tt;t++)incidence.push_back(arma::mat(n,n,arma::fill::zeros));for(unsigned z=0;z<sy.n_elem;z++){unsigned i=ii[z],j=jj[z],t=time[z];counts(i,j)++;counts(j,i)++;incidence[t](i,i)++;incidence[t](j,j)++;incidence[t](i,j)++;incidence[t](j,i)++;hss.col(t)+=arma::square(H.row(i)+H.row(j)).t();}for(int t=0;t<tt;t++)incidence_sum+=incidence[t];beta_precision=arma::sum(arma::square(X*B),0).t();}
};
struct TemperedBinaryEvaluation {double ll,lp;arma::vec g;bool valid;};
arma::vec tempered_binary_scale_density(double theta) {
 const double b=.5, pi=3.14159265358979323846;
 double at=std::abs(theta),lp,grad;
 if(at<=b/10){
  double x=(theta/b)*(theta/b),term=1,sum=1,weighted=0;
  for(int k=1;k<200;k++){double next=-term*(2*k+1)*x;if(std::abs(next)>std::abs(term))break;term=next;sum+=term;weighted+=k*term;if(std::abs(term)<1e-16*std::abs(sum))break;}
  lp=-std::log(b)-.5*std::log(2*pi)+std::log(sum);grad=theta==0?0:2*weighted/(theta*sum);
 }else{
  double z=b/at,rr=std::exp(.5*std::log(2*pi)+.5*z*z+R::pnorm(-z,0,1,1,1)),dd=1-z*rr;
  lp=std::log(b)-.5*std::log(2*pi)-2*std::log(at)+std::log(dd);grad=(-2-z*(z-(1+z*z)*rr)/dd)/theta;
 }
 return arma::vec({lp,grad});
}
TemperedBinaryEvaluation tempered_binary_target(const arma::vec& q,const TemperedBinaryContext& c,double temp){
 TemperedBinaryEvaluation e;e.ll=0;e.lp=R_NegInf;e.g.zeros(q.n_elem);e.valid=false;
 if(q.n_elem!=(unsigned)(c.base+3)||!q.is_finite()||!R_finite(temp)||temp<0||temp>1)return e;
 double va=std::exp(q[c.base]),sigma=std::exp(q[c.base+1]),th=std::tanh(q[c.base+2]),rho=.99*th,dr=.99*(1-th*th);
 if(!R_finite(va)||va<=0||!R_finite(sigma)||sigma<=0)return e;
 int k=0;arma::mat Uw=arma::reshape(q.subvec(k,k+c.n*c.r-1),c.n,c.r);k+=c.n*c.r;
 arma::vec L=q.subvec(k,k+c.r-1);k+=c.r;arma::vec beta=c.B*q.subvec(k,k+c.p-1);k+=c.p;
 arma::vec level=std::sqrt(va)*c.H*q.subvec(k,k+c.n-2);k+=c.n-1;
 arma::mat innov=arma::reshape(q.subvec(k,k+(c.n-1)*c.tt-1),c.n-1,c.tt),dev=innov,dd=innov;
 dev.col(0)=sigma/std::sqrt(1-rho*rho)*innov.col(0);dd.col(0)=dev.col(0)*rho/(1-rho*rho);
 for(int t=1;t<c.tt;t++){dev.col(t)=rho*dev.col(t-1)+sigma*innov.col(t);dd.col(t)=dev.col(t-1)+rho*dd.col(t-1);}
 arma::mat path=c.H*dev;path.each_col()+=level;
 arma::vec sc(c.r,arma::fill::ones);if(c.hierarchical)sc=arma::exp(q.subvec(c.base+3,c.base+2+c.r));
 arma::mat U=Uw;U.each_row()%=sc.t();arma::mat factor=(U.each_row()%L.t())*U.t();
 arma::mat D(c.n,c.n,arma::fill::zeros),ga(c.n,c.tt,arma::fill::zeros);arma::vec gb(c.p,arma::fill::zeros);
 for(unsigned z=0;z<c.sy.n_elem;z++){
  unsigned i=c.ii[z],j=c.jj[z],t=c.time[z];double eta=arma::dot(c.X.row(z),beta)+path(i,t)+path(j,t)+factor(i,j);
  double lp=R::pnorm(c.sy[z]*eta,0,1,1,1);double g=c.sy[z]*std::exp(R::dnorm(eta,0,1,1)-lp);
  if(!R_finite(lp)||!R_finite(g))return e;
  e.ll+=lp;D(i,j)+=g;D(j,i)+=g;ga(i,t)+=g;ga(j,t)+=g;gb+=g*c.X.row(z).t();
 }
 arma::mat du=D*U,gu=du;gu.each_row()%=L.t();gu.each_row()%=sc.t();
 arma::vec gl=.5*arma::sum(U%du,0).t(),glev=std::sqrt(va)*c.H.t()*arma::sum(ga,1);
 arma::mat gd=c.H.t()*ga;for(int t=c.tt-2;t>=0;t--)gd.col(t)+=rho*gd.col(t+1);gd*=sigma;gd.col(0)/=std::sqrt(1-rho*rho);
 k=0;e.g.subvec(k,k+c.n*c.r-1)=arma::vectorise(gu);k+=c.n*c.r;e.g.subvec(k,k+c.r-1)=gl;k+=c.r;
 e.g.subvec(k,k+c.p-1)=c.B.t()*gb;k+=c.p;e.g.subvec(k,k+c.n-2)=glev;k+=c.n-1;e.g.subvec(k,c.base-1)=arma::vectorise(gd);
 e.g.head(c.base)=temp*e.g.head(c.base)-q.head(c.base);
 double normal_sum=arma::dot(q.head(c.n*c.r),q.head(c.n*c.r))+arma::dot(q.subvec(c.n*c.r+c.r,c.base-1),q.subvec(c.n*c.r+c.r,c.base-1));if(!c.collapsed)normal_sum+=arma::dot(L,L);
 double prior=-.5*normal_sum-c.av*q[c.base]-c.bv/va-2*c.ash*q[c.base+1]-c.bs/(sigma*sigma)-.5*std::pow((rho-c.rm)/c.rs,2)+(std::log(.99)+2*(std::log(2.0)-std::abs(q[c.base+2])-std::log1p(std::exp(-2*std::abs(q[c.base+2])))));
 e.g[c.base]=temp*.5*arma::dot(arma::sum(ga,1),level)-c.av+c.bv/va;
 e.g[c.base+1]=temp*arma::accu(ga%(c.H*dev))-2*c.ash+2*c.bs/(sigma*sigma);
 e.g[c.base+2]=temp*arma::accu(ga%(c.H*dd))*dr-(rho-c.rm)/(c.rs*c.rs)*dr-2*th;
 if(c.hierarchical){
  for(int j=0;j<c.r;j++){double ls=q[c.base+3+j],iv=std::exp(-2*ls);prior+=-2*ls-.5*iv;e.g[c.base+3+j]=temp*L[j]*arma::dot(du.col(j),U.col(j))-2+iv;}
 }
 if(c.collapsed){for(int j=0;j<c.r;j++){int pos=c.n*c.r+j;arma::vec pg=tempered_binary_scale_density(q[pos]);prior+=pg[0];e.g[pos]=temp*gl[j]+pg[1];}}
 e.lp=temp*e.ll+prior;e.valid=R_finite(e.lp)&&e.g.is_finite();return e;
}
// [[Rcpp::export]]
List tempered_binary_evaluate(arma::vec q,List context,double temperature=1){TemperedBinaryContext c(context);TemperedBinaryEvaluation e=tempered_binary_target(q,c,temperature);return List::create(_["logdensity"]=e.lp,_["loglik"]=e.ll,_["gradient"]=e.g,_["valid"]=e.valid);}
// diagonal positive Fisher metric; likelihood curvature uses its unit upper bound.
arma::vec tempered_binary_precision(const arma::vec& q,const TemperedBinaryContext& c,double temperature){
 if(q.n_elem!=(unsigned)(c.base+3)||!q.is_finite()||!R_finite(temperature)||temperature<0||temperature>1)stop("Invalid tempered metric state.");
 arma::vec prec(q.n_elem,arma::fill::ones);
 double va=std::exp(q[c.base]),sigma=std::exp(q[c.base+1]),th=std::tanh(q[c.base+2]),rho=.99*th,dr=.99*(1-th*th);
 arma::mat U=arma::reshape(q.head(c.n*c.r),c.n,c.r);arma::vec L=q.subvec(c.n*c.r,c.n*c.r+c.r-1);
 int levpos=c.n*c.r+c.r+c.p,ip=levpos+c.n-1;
 arma::vec level=std::sqrt(va)*c.H*q.subvec(levpos,ip-1);
 arma::mat innov=arma::reshape(q.subvec(ip,c.base-1),c.n-1,c.tt),dev=innov,dd=innov;
 dev.col(0)=sigma/std::sqrt(1-rho*rho)*innov.col(0);dd.col(0)=dev.col(0)*rho/(1-rho*rho);
 for(int t=1;t<c.tt;t++){dev.col(t)=rho*dev.col(t-1)+sigma*innov.col(t);dd.col(t)=dev.col(t-1)+rho*dd.col(t-1);}
 arma::mat path=c.H*dev,dpath=c.H*dd;
 arma::mat uu=arma::square(U);
 for(int k=0;k<c.r;k++){arma::vec cu=c.counts*uu.col(k);prec.subvec(k*c.n,(k+1)*c.n-1)+=temperature*L[k]*L[k]*cu;prec[c.n*c.r+k]+=temperature*.5*arma::dot(uu.col(k),cu);}
 prec.subvec(c.n*c.r+c.r,levpos-1)+=temperature*c.beta_precision;
 prec.subvec(levpos,ip-1)+=temperature*va*arma::sum(c.hss,1);
 for(int a=0;a<c.tt;a++){arma::vec v(c.n-1,arma::fill::zeros);for(int t=a;t<c.tt;t++)v+=std::pow(rho,2*(t-a))*c.hss.col(t);double fac=temperature*sigma*sigma;if(a==0)fac/=(1-rho*rho);prec.subvec(ip+a*(c.n-1),ip+(a+1)*(c.n-1)-1)+=fac*v;}
 prec[c.base]+=temperature*.25*arma::as_scalar(level.t()*c.incidence_sum*level);
 for(int t=0;t<c.tt;t++){prec[c.base+1]+=temperature*arma::as_scalar(path.col(t).t()*c.incidence[t]*path.col(t));prec[c.base+2]+=temperature*dr*dr*arma::as_scalar(dpath.col(t).t()*c.incidence[t]*dpath.col(t));}
 if(c.collapsed)for(int k=0;k<c.r;k++)prec[c.n*c.r+k]+=24/(1+100*L[k]*L[k]);
 prec[c.base]+=c.bv/va;prec[c.base+1]+=4*c.bs/(sigma*sigma);prec[c.base+2]+=dr*dr/(c.rs*c.rs)+2*(1-th*th);
 return prec;
}
struct TemperedBinaryProposal {arma::vec q,candidate,sd,drift,new_sd,new_drift;TemperedBinaryEvaluation ev;double prob,ratio;bool accepted,valid;};
TemperedBinaryProposal tempered_binary_proposal(const arma::vec& q,const TemperedBinaryContext& c,double temp,double eps,const TemperedBinaryEvaluation& old){
 TemperedBinaryProposal z;z.q=q;z.ev=old;z.prob=0;z.ratio=R_NegInf;z.accepted=false;z.valid=false;
 z.sd=eps/arma::sqrt(tempered_binary_precision(q,c,temp));z.drift=z.sd%arma::tanh(.5*z.sd%old.g);
 arma::vec noise(q.n_elem);for(unsigned k=0;k<noise.n_elem;k++)noise[k]=R::rnorm(0,1);
 z.candidate=q+z.drift+z.sd%noise;TemperedBinaryEvaluation en=tempered_binary_target(z.candidate,c,temp);
 if(!en.valid)return z;
 z.new_sd=eps/arma::sqrt(tempered_binary_precision(z.candidate,c,temp));z.new_drift=z.new_sd%arma::tanh(.5*z.new_sd%en.g);
 if(!z.sd.is_finite()||!z.new_sd.is_finite()||arma::any(z.sd<=0)||arma::any(z.new_sd<=0))return z;
 arma::vec reverse=(q-z.candidate-z.new_drift)/z.new_sd;
 arma::vec forward_noise=(z.candidate-q-z.drift)/z.sd;
 double forward=-arma::accu(arma::log(z.sd))-.5*arma::dot(forward_noise,forward_noise),backward=-arma::accu(arma::log(z.new_sd))-.5*arma::dot(reverse,reverse);
 z.ratio=en.lp-old.lp+backward-forward;z.valid=R_finite(z.ratio);if(!z.valid)return z;
 z.prob=std::exp(std::min(0.0,z.ratio));z.accepted=R::runif(0,1)<z.prob;
 if(z.accepted){z.q=z.candidate;z.ev=en;}return z;
}
// [[Rcpp::export]]
List tempered_binary_transition(arma::vec q,List context,double temperature,double stepsize){TemperedBinaryContext c(context);TemperedBinaryEvaluation old=tempered_binary_target(q,c,temperature);if(!old.valid)stop("invalid initial state");TemperedBinaryProposal z=tempered_binary_proposal(q,c,temperature,stepsize,old);return List::create(_["q"]=z.q,_["candidate"]=z.candidate,_["sd"]=z.sd,_["drift"]=z.drift,_["new_sd"]=z.new_sd,_["new_drift"]=z.new_drift,_["log_ratio"]=z.ratio,_["accepted"]=z.accepted,_["valid"]=z.valid);}

arma::vec tempered_binary_prior(const TemperedBinaryContext& c) {
 arma::vec q(c.base+3);
 for(int k=0;k<c.base;k++)q[k]=R::rnorm(0,1);
 if(c.collapsed)for(int k=0;k<c.r;k++)q[c.n*c.r+k]/=R::rgamma(1,2);
 q[c.base]=-std::log(R::rgamma(c.av,1/c.bv));
 q[c.base+1]=-.5*std::log(R::rgamma(c.ash,1/c.bs));
 double rho;do{rho=R::rnorm(c.rm,c.rs);}while(std::abs(rho)>=.99);
 q[c.base+2]=std::atanh(rho/.99);return q;
}
// [[Rcpp::export]]
List tempered_binary_reflection(arma::vec q,List context,double temperature,int column){
 TemperedBinaryContext c(context);TemperedBinaryEvaluation old=tempered_binary_target(q,c,temperature);if(!old.valid||column<0||column>=c.r)stop("Invalid signed-weight reflection state.");arma::vec qp=q;qp[c.n*c.r+column]*=-1;TemperedBinaryEvaluation en=tempered_binary_target(qp,c,temperature);
 double ratio=temperature*(en.ll-old.ll);bool accepted=en.valid&&std::log(R::runif(0,1))<ratio;
 return List::create(_["q"]=accepted?qp:q,_["proposal"]=qp,_["accepted"]=accepted,_["log_ratio"]=ratio,_["valid"]=en.valid);
}
// [[Rcpp::export]]
List tempered_binary_block(arma::mat q,List context,arma::vec temperatures,arma::vec stepsizes,int start,int count,int burn,int thin){
 TemperedBinaryContext c(context);
 if(q.n_rows!=(unsigned)(c.base+3)||q.n_cols<1||q.n_cols>1000||!q.is_finite()||
    temperatures.n_elem!=q.n_cols||stepsizes.n_elem!=q.n_cols||!temperatures.is_finite()||!stepsizes.is_finite()||
    arma::any(temperatures<0)||arma::any(temperatures>1)||arma::any(arma::diff(temperatures)>=0)||
    arma::any(stepsizes<=0)||start<0||count<1||burn<0||thin<1||((double)start+count)>2147483647.0)
   stop("Invalid tempered block state or sampling schedule.");
 int K=q.n_cols;arma::vec acc(K,arma::fill::zeros),bad=acc,att=acc,swapacc(K-1,arma::fill::zeros),swapatt=swapacc;
 arma::vec reflection_acc(K,arma::fill::zeros),reflection_att(K,arma::fill::zeros);int prior_refreshes=0;
 std::vector<arma::vec> draws,mode_trace;arma::vec ll(K);std::vector<TemperedBinaryEvaluation> ev;
 for(int j=0;j<K;j++){ev.push_back(tempered_binary_target(q.col(j),c,temperatures[j]));if(!ev.back().valid)stop("invalid initial tempered state");ll[j]=ev.back().ll;}
 for(int s=start+1;s<=start+count;s++){
  if(s%10==0)checkUserInterrupt();
  for(int j=0;j<K;j++){
   if(temperatures[j]==0){arma::vec qp=tempered_binary_prior(c);TemperedBinaryEvaluation en=tempered_binary_target(qp,c,0);if(!en.valid)stop("Nonfinite exact prior refresh");q.col(j)=qp;ev[j]=en;ll[j]=en.ll;prior_refreshes++;continue;}
   TemperedBinaryProposal z=tempered_binary_proposal(q.col(j),c,temperatures[j],stepsizes[j],ev[j]);att[j]++;bad[j]+=!z.valid;acc[j]+=z.accepted;
   q.col(j)=z.q;ev[j]=z.ev;ll[j]=z.ev.ll;
   if(s<=burn){double logeps=std::log(stepsizes[j])+(z.prob-.57)/std::sqrt((double)s+10);stepsizes[j]=std::min(2.5,std::max(1e-5,std::exp(logeps)));}
  }
  if(s%5==0)for(int j=0;j<K-1;j++){
   int column=(int)std::floor(R::runif(0,c.r));arma::vec qp=q.col(j);qp[c.n*c.r+column]*=-1;
   TemperedBinaryEvaluation en=tempered_binary_target(qp,c,temperatures[j]);reflection_att[j]++;bad[j]+=!en.valid;
   if(en.valid&&std::log(R::runif(0,1))<temperatures[j]*(en.ll-ll[j])){q.col(j)=qp;ev[j]=en;ll[j]=en.ll;reflection_acc[j]++;}
  }
  for(int j=s%2;j<K-1;j+=2){swapatt[j]++;double lr=(temperatures[j]-temperatures[j+1])*(ll[j+1]-ll[j]);if(std::log(R::runif(0,1))<lr){q.swap_cols(j,j+1);std::swap(ll[j],ll[j+1]);ev[j]=tempered_binary_target(q.col(j),c,temperatures[j]);ev[j+1]=tempered_binary_target(q.col(j+1),c,temperatures[j+1]);swapacc[j]++;}}
  if(s>burn&&(s-burn)%thin==0){draws.push_back(q.col(0));arma::vec modes(K);for(int j=0;j<K;j++)modes[j]=arma::accu(q.col(j).subvec(c.n*c.r,c.n*c.r+c.r-1)<0);mode_trace.push_back(modes);}
 }
 arma::mat out(q.n_rows,draws.size()),mt(K,draws.size());for(unsigned i=0;i<draws.size();i++){out.col(i)=draws[i];mt.col(i)=mode_trace[i];}
 return List::create(_["q"]=q,_["stepsizes"]=stepsizes,_["draws"]=out,_["mode_trace"]=mt,_["attempts"]=att,_["accepts"]=acc,_["nonfinite"]=bad,_["swap_attempts"]=swapatt,_["swap_accepts"]=swapacc,_["loglik"]=ll,_["reflection_attempts"]=reflection_att,_["reflection_accepts"]=reflection_acc,_["prior_refreshes"]=prior_refreshes);
}
