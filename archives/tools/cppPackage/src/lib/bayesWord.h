//////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __bayesWord_h
#define __bayesWord_h
//////////////////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
using std::string; using std::vector;
typedef vector<double> dVec;
typedef vector<dVec>   dMat;
typedef vector<string> strVec;
typedef vector<int>    iVec;
//////////////////////////////////////////////////////////////////////////////////////////////
class BayesWordModel {
 public:
  typedef enum{DNA,RNA,Protein} SeqType;
 protected:
  dVec bFr_,bMod_;
  dMat wFr_,wMod_;
  double bPs_,wPs_;
  int N_,nHit_; SeqType seqT_;
  inline int BaseToIndex(char b) const;
 public:
  BayesWordModel():bFr_(),bMod_(),wFr_(),wMod_(),bPs_(),wPs_(),seqT_(DNA),N_(),nHit_() {}
  BayesWordModel(const BayesWordModel& o):
    bFr_(o.bFr_),bMod_(o.bMod_),wFr_(o.wFr_),wMod_(o.wMod_),bPs_(o.bPs_),wPs_(o.wPs_),
    seqT_(o.seqT_),N_(o.N_),nHit_(o.nHit_) {}
  BayesWordModel(const strVec& s,const iVec& p,int N,double pScaleM,double pScaleBG,
		 const dVec& totalBGCounts=dVec(),int omit= -1,const dMat& prior=dMat(),
		 SeqType seqType=DNA);
  BayesWordModel(const dMat& wordModel,const dVec& bgModel,int n,SeqType seqType=DNA);//:
  BayesWordModel(const dMat& wCounts,const dVec& bgCounts,int n,double pScaleM,double pScaleBG,
		 SeqType seqType=DNA);
  //
  BayesWordModel& operator=(const BayesWordModel& o) { 
    if(*this == o) return *this;
    bFr_=o.bFr_; bMod_=o.bMod_; wFr_=o.wFr_; wMod_=o.wMod_; bPs_=o.bPs_; wPs_=o.wPs_; 
    N_=o.N_; seqT_=o.seqT_; nHit_=o.nHit_;
    return *this;
  }
  bool operator==(const BayesWordModel& o) const {
    return ((bFr_==o.bFr_)&&(bMod_==o.bMod_)&&(wFr_==o.wFr_)&&(wMod_==o.wMod_)&&(bPs_==o.bPs_)&&
	    (wPs_==o.wPs_)&&(N_==o.N_)&&(nHit_==o.nHit_)&&(seqT_==o.seqT_));
  }
  //
  double logPword(const string& w)        const;
  double logPwordBG(const string& w)      const;
  double logPwordRatio(const string& w)   const;
  double logOddsWord(const string& w)     const;
  double MAP()                            const;
  dVec SeqProbVector(const string& s)     const;
  dVec SeqProbRatioVector(const string& s)const;
  dVec SeqProbVector(const string& s,double& pMiss)     const;
  dVec SeqProbRatioVector(const string& s,double& pMiss)const;
  dVec SummedProbVector(const strVec& s)  const;
  dVec TotalCounts()                      const;
  //
  dVec logBGModel()            const { return dVec(bMod_); }
  dMat logWordModel()          const { return dMat(wMod_); }
  dVec BGCounts()              const { return dVec(bFr_); }
  dMat WordCounts()            const { return dMat(wFr_); }
  bool IsProtein()             const { return (seqT_==Protein); }
  bool IsDNA()                 const { return (seqT_==DNA); }
  bool IsRNA()                 const { return (seqT_==RNA); }
  int N()                      const { return N_; }
  double bgPseudocountWeight() const { return bPs_; }
  double wPseudocountWeight()  const { return wPs_; }
};
//////////////////////////////////////////////////////////////////////////////////////////////
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
