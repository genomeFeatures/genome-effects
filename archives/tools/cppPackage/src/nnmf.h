//////////////////////////////////////////////////////////////////////////
#ifndef __nnmf_h
#define __nnmf_h
#include <vector>
#include <string>
//////////////////////////////////////////////////////////////////////////
typedef std::vector<double> dVec;
typedef std::vector<dVec> dMat;
//////////////////////////////////////////////////////////////////////////
class NNMF {
 private:
  dMat V_,W_,H_;
  int r_, metric_,maxIt_,itConv_;
  double tol_,diff_,res_;
  NNMF() {}
 public:
  //NNMF(const NNMF& o):V_(o.V_),W_(o.W_),H_(o.H_),r_(o.r_),metric_(o.metric_),tol_(o.tol_),
  //  diff_(o.diff_),maxIt_(o.maxIt_),itConv_(o.itConv_) {}
  NNMF(const dMat& V,int r, double tol=0.001,int maxIt=1000,bool quiet=false,int metric=0);
  dMat bases() const { return dMat(W_); }
  dMat weights() const { return dMat(H_); }
  int nBases() const { return r_; }
  int Iterations() const { return itConv_; }
  double tolerance() const { return tol_; }
  double Fscore() const { return diff_; }
  double Residue() const { return res_; }
};
//////////////////////////////////////////////////////////////////////////
#endif
//////////////////////////////////////////////////////////////////////////
