///////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __WORDCOUNT_H
#define __WORDCOUNT_H
///////////////////////////////////////////////////////////////////////////////////////////////////
#include <functional>
#include <vector>
#include <string>
#include <map>
using std::vector; using std::string; using std::map; using std::less;
typedef vector<unsigned> uVec;
typedef vector<int> iVec;
typedef vector<iVec> iMat;
typedef vector<double> dVec;
typedef vector<dVec> dMat;
typedef vector<string> strVec;
typedef map<string,unsigned,less<string> > strUintMap;
typedef vector<strUintMap> suMapVec;
typedef map<string,uVec,less<string> > strUVecMap;
typedef vector<strUVecMap> suvMapVec;
///////////////////////////////////////////////////////////////////////////////////////////////////
const unsigned maxWordSize=20;
///////////////////////////////////////////////////////////////////////////////////////////////////
class cjWordCount {
 protected:
  unsigned mN_,nS_; // mN_ = max word length; nS_= # of seqs counted
  uVec nW_; // nW_[i] = # of words of length i+1
  suMapVec wc_; // the counts, a vector of maps from string to int
  bool ds_; // ds_ = true if we are counting double stranded
 public:
  cjWordCount():mN_(),nS_(),nW_(),wc_(),ds_() {}
  cjWordCount(const cjWordCount& o):mN_(o.mN_),nS_(o.nS_),nW_(o.nW_),wc_(o.wc_),ds_(o.ds_) {}
  cjWordCount(unsigned maxN,bool ds=false);
  cjWordCount(const strVec& seqs,unsigned maxN,bool ds=false);
  bool operator==(const cjWordCount& o) const {
    return ((mN_==o.mN_)&&(nS_==o.nS_)&&(nW_==o.nW_)&&(wc_==o.wc_)&&(ds_==o.ds_));
  }
  cjWordCount& operator=(const cjWordCount& o);
  void AddSequences(const strVec& seqs);
  virtual void AddSequence(const string& seq);
  unsigned WordCount(const string& word) const;
  double WordFrac(const string& word) const;
  double WordPerc(const string& word) const;
  double ZScore(const string& word,unsigned order) const;
  double KarlinRho(const string& w) const;
  double Entropy(unsigned N) const;
  double Information(unsigned N) const;
  unsigned MaxWordLength() const { return mN_; }
  unsigned SequenceCount() const { return nS_; }
  bool DoubleStrandCount() const { return ds_; }
  unsigned UniqueWordCount(unsigned N) const;
  unsigned TotalWordCount(unsigned N) const;
  unsigned AmbiguousCount(unsigned N) const;
  double ExpectM(const string& w,unsigned o,unsigned L) const;
  double ExpectFreq(const string& w, unsigned o) const;
  //		double ExpectP(const string& w,unsigned o,unsigned L) const;
  void WordExpCV(const string& word,unsigned order,unsigned L,
		 double& E, double& V) const;
  string GenerateSequence(unsigned N,unsigned order) const;
  strVec ListWords() const;
  string ListWordCounts() const;
  unsigned LoadCountsFromFile(const string& fn);
  double MarkovPredF(const string& w,int order) const; 
};
///////////////////////////////////////////////////////////////////////////////////////////////////
class cjPosWordCount:public cjWordCount {
 protected:
  suvMapVec wp_; 
  unsigned minSeqL_,maxSeqL_;
 public:
  cjPosWordCount():cjWordCount(),wp_(),minSeqL_(),maxSeqL_() {}
  cjPosWordCount(const cjPosWordCount& o):cjWordCount(o),wp_(o.wp_),
    minSeqL_(o.minSeqL_),maxSeqL_(o.maxSeqL_) {}
  cjPosWordCount(unsigned maxN,bool ds=false);
  cjPosWordCount(const strVec& seqs,unsigned maxN,bool ds=false);
  void AddSequence(const string& seq);
  bool operator==(const cjPosWordCount& o) const {
    return ((mN_==o.mN_)&&(nS_==o.nS_)&&(nW_==o.nW_)&&(wc_==o.wc_)&&(wp_==o.wp_));
  }
  cjPosWordCount& operator=(const cjPosWordCount& o);
  double PosMean(const string& w) const;
  double PosVar(const string& w) const;
  unsigned LongestSeq() const { return maxSeqL_; }
  unsigned ShortestSeq() const { return minSeqL_; }
  double ChiSqr(const string& w,unsigned nw,dVec& bins) const;
  double zChiSqr(const string& w,unsigned wL,int zPos,dVec& bins,double& avgBinCount) const;
  unsigned SimProfWords(const string& w,strVec& words,dVec& aBins,
			unsigned nw,double minR,double minChi,const strVec& fWords=strVec()) const;
  void GenerateRandomChi(unsigned nw,unsigned order,unsigned nTest,unsigned nSeqs,
			 unsigned seqL,dVec& chiBound) const;
  unsigned RangeWordCount(string& w,unsigned l,unsigned h) const;
  void KMeanCluster(unsigned k,unsigned n,unsigned wL,int zPos,double minAvgBinCount,
		    unsigned maxIt,dMat& bins,strVec& words,iVec& cID,dVec& Ssq,
		    dMat& clusterCenters,dMat& clusterPr) const;
};
///////////////////////////////////////////////////////////////////////////////////////////////////
//  double BestClusterP(const string& s,const string& w,int& start,int& stop) const;
///////////////////////////////////////////////////////////////////////////////////////////////////
class ClusterProfile {
 protected:
  int N_; string s_; dVec bc_; iMat wc_; iVec s1_,s2_,sz_;
 public:
  ClusterProfile():N_(),s_(),bc_(),wc_(),s1_(),s2_(),sz_() {}
  ClusterProfile(const ClusterProfile& o):N_(o.N_),s_(o.s_),bc_(o.bc_),wc_(o.wc_),s1_(o.s1_),s2_(o.s2_),sz_(o.sz_) {}
  ClusterProfile(int n,const string& s,const cjWordCount& wc);
  double pearsCorr(const ClusterProfile& o) const;// { return Pearson()(bc_,o.bc_); }
  double spearCorr(const ClusterProfile& o) const;// { return Spearman()(bc_,o.bc_); }
  dVec pVec() const { return dVec(bc_); }
  iVec startVec() const { return iVec(s1_); }
  iVec stopVec() const { return iVec(s2_); }
  iVec sizeVec() const { return iVec(sz_); }
  iMat cMat() const { return iMat(wc_); }
};
///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////
