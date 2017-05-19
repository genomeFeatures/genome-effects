///////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <strstream>
#include "../lib/logFile.h"
#include "../lib/wordStat.h"
#include "../lib/util.h"
///////////////////////////////////////////////////////////////////////////////
iVec PrincipalPeriods::operator()(const string& w) const {
  iVec p;
  for(int i=1;i<w.length();i++) {
    bool per=true; int j=0;
    while(((j+i)<w.length())&&per) {
      per= (w[j]==w[i+j]); ++j;
    }
    if(per) {
      bool princ=true;
      for(unsigned k=0;(k<p.size())&&princ;k++) princ= (i % p[k]);
      if(princ) p.push_back(i);
    }
  }
  return iVec(p);
}
///////////////////////////////////////////////////////////////////////////////
iVec AllPeriods::operator()(const string& w) const {
  iVec p;
  for(int i=1;i<w.length();i++) {
    bool per=true; int j=0;
    while(((j+i)<w.length())&&per) {
      per= (w[j]==w[i+j]); ++j;
    }
    if(per) p.push_back(i);
  }
  return iVec(p);
}
///////////////////////////////////////////////////////////////////////////////
double SchbathA::operator()(const string& w,const dVec& Ppre) const {
  if(Ppre.size()!=w.length()) {
    WriteToLogs()(string("SchbathA: incorrect prefix probability vector size"));
    return 0.0;
  }
  iVec pi=PrincipalPeriods()(w);
  double A=0.0;
  for(unsigned i=0;i<pi.size();i++) {
    A += Ppre[pi[i]]; 
  }
  return A;
}
///////////////////////////////////////////////////////////////////////////////
CompPoiss::CompPoiss(const string& w,int n,int x,const dVec& Ppre,double Pw):
w_(w),n_(n),x_(x),Ppre_(Ppre),Pw_(Pw),Pval_() {
  iVec pi=PrincipalPeriods()(w_);
  double A=((pi.size()>0)?SchbathA()(w_,Ppre_):0.0); 
  int l=w_.length();
  int Lmin((pi.size()>0)?pi[0]*(n_-1)+l:n_*l);
  Pval_=0.0;
  double logA=log(A),logsA=log(1.0-A)-log(A),lfn=logFac()(n_),lfn1=logFac()(n_-1);
  double nlogA=n_*logA;
  for(unsigned L=Lmin;L<=x_;L++) {
    double lambda= (L-l+1)*(1.0-A)*Pw_;
    double loglam= log(lambda);
    if(pi.size()>0) {
      double sum=0.0;
      for(unsigned j=1;j<=n_;j++) 
	sum += exp(lfn1-logFac()(j-1)-logFac()(n_-j)-logFac()(j)+static_cast<double>(j)*(loglam+logsA));
      Pval_+= exp(-1.0*lambda + nlogA + log(sum));
    } else {
      Pval_+= exp(-1.0*lambda + n_*loglam - lfn);
    }
  }
  Pval_ *= Pw_/(1.0-A);
  char line[501]; std::strstream buff(line,500);
  buff << w << ": n=" << n << ", x=" << x << ", Pw=" << Pw_ << ", A=" << A << ", Lmin=" << Lmin << ", P=" << Pval_ << '\0';
  WriteToLog()(line);
}
///////////////////////////////////////////////////////////////////////////////
