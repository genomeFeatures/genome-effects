#include "nnmf.h"
#include "lib/logFile.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
/////////////////////////////////////////////////////////////////////////////////////////////////////////
double euclid(const dMat& a,const dMat& b) {
  double euc=0.0;
  for(int i=0;i<a.size();++i) 
    for(int j=0;j<a[0].size();++j) {
      //std::cout  << "[" << i << ":" << j << "] euc= " << euc << ", a= " << a[i][j] << ", b= " << b[i][j] << '\n';
      euc+= (a[i][j]-b[i][j])*(a[i][j]-b[i][j]);
    }
  return sqrt(euc);///static_cast<double>(a.size()*a[0].size());
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
double div(const dMat& a,const dMat& b) {
  double div=0.0;
  for(int i=0;i<a.size();++i)
    for(int j=0;j<a[0].size();++j) 
      //div += fabs(a[i][j]*log(a[i][j]+1e-20) - b[i][j]);
      div += a[i][j]*log(b[i][j]+1e-100) - b[i][j];
  return div;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
dMat matmult(const dMat& A, const dMat& B) {
  dMat C(A.size(),dVec(B[0].size(),0.0));
  int N=A.size(),M=B[0].size(),r=A[0].size();
  for(int i=0;i<N;++i)
    for(int j=0;j<M;++j)
      for(int k=0;k<r;++k) {
	//	std::cout << i << " " << j << " " << k << std::endl;
	C[i][j] += A[i][k]*B[k][j];
      }
  return C;
}  
/////////////////////////////////////////////////////////////////////////////////////////////////////////
NNMF::NNMF(const dMat& Vmat,int r, double tol,int maxIt,bool quiet,int metric):
  V_(Vmat),r_(5),tol_(0.001),metric_(metric),maxIt_(maxIt),itConv_(0) {
  if(r>0 && r<20) r_=r;
  if(tol_>1e-10 && tol_<10.0) tol_=tol;
  if(V_.empty()) return;
  int N=V_.size(), M=V_[0].size();
  for(int i=0;i<N;i++) if(V_[i].size()!=M) return;
  // okay to here- now set up the necessary arrays and matrix views
  ////////////////////
  // gsl initialization stuff
  time_t t0; time(&t0); 
  const gsl_rng_type *T;
  gsl_rng *rnd;
  gsl_rng_env_setup();
  T= gsl_rng_default;
  gsl_rng_default_seed= t0%10000;
  rnd= gsl_rng_alloc(T);
  // end gsl initialization
  int NM=N*M,Nr=N*r_,rM=r_*M,rr=r_*r_;
  //
  H_= dMat(r_,dVec(M,0.5));
  for(int i=0;i<r_;++i) for(int j=0;j<M;++j) H_[i][j]= gsl_rng_uniform(rnd);
  W_= dMat(N,dVec(r_,1.0));
  for(int i=0;i<N;++i) for(int j=0;j<r_;++j) W_[i][j]= gsl_rng_uniform(rnd);
  gsl_rng_free(rnd);
  dMat WH= matmult(W_,H_);
  // diff_= euclid(V_,WH);
  diff_= div(V_,WH);
  if(!quiet) std::cout << "\nN= " << N << ", M=" << M << ", r= " << r_ << ", diff= " << diff_ << '\n';

  double lastD=diff_,delta=1.0;
  while(itConv_++ < maxIt && delta > tol_) { // && delta > tol_) {
    dMat lastWH(WH);
    // update paired columns of W and rows of H 
    for(int a=0;a<r_;++a) {
      // column (N long) of W first
      double colsum=0.0;
      dVec newW(N,0.0);
      for(int i=0;i<N;++i) {
	double sum=0.0;
	for(int mu=0;mu<M;++mu) 
	  if(WH[i][mu]>0.0) sum += H_[a][mu]*V_[i][mu]/WH[i][mu];
	  else if (H_[a][mu]*V_[i][mu] !=0.0)
	    std::cout << "** trouble in river city (col update), a,mu,i=" << a << "," << mu << ',' << i << " **\n";
	newW[i]= W_[i][a] * sum;
	colsum+= newW[i];
      }
      // normalize the column of W
      for(int i=0;i<N;++i) W_[i][a] = newW[i]/colsum;
      // remake an updated WH
      //WH= matmult(W_,H_);
      // now the corresponding row (M long) of H
      dVec newH(M,0.0);
      for(int mu=0;mu<M;++mu) {
	double sum=0.0;
	for(int i=0;i<N;++i) 
	  if(WH[i][mu]>0.0) sum += W_[i][a]*V_[i][mu]/WH[i][mu];
	  else if (H_[a][mu]*V_[i][mu] !=0.0)
	    std::cout << "** trouble in river city (row update), a,mu,i=" << a << "," << mu << ',' << i << " **\n";
	newH[mu]= H_[a][mu] * sum;
      }
      H_[a]=newH;
      // remake the updated WH
      WH= matmult(W_,H_);
    }
    // now recalculate diff_
    WH= matmult(W_,H_);
    diff_= div(V_,WH);
    res_= euclid(V_,WH);
    delta=0.0;
    if(diff_+lastD !=0.0) delta= 2.0*fabs(diff_-lastD)/fabs(diff_+lastD);
    lastD=diff_;
    if(!quiet) 
      std::cout << "# iteration " << itConv_ << ", F= " << std::setprecision(10) << diff_ << ", delta= " 
		<< std::setprecision(5) << delta << " #" << std::endl;
  }
  return;
}

