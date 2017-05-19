///////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __WCOUNT_H
#define __WCOUNT_H
///////////////////////////////////////////////////////////////////////////////////////////////////
#include <functional>
#include <vector>
#include <string>
#include <map>
using std::vector; using std::string; using std::map; using std::less;
//vectors
typedef vector<unsigned> uVec;
typedef vector<int> iVec;
typedef vector<double> dVec;
typedef vector<string> strVec;
//matrixes
typedef vector<iVec> iMat;
typedef vector<dVec> dMat;
//maps
typedef map<string,unsigned,less<string> > strUintMap; // string to unsigned (for counting)
typedef map<string,int,less<string> > strIntMap;       // string to int
typedef map<string,uVec,less<string> > strUVecMap;     // string to uVec (for storing occurrences)
typedef map<string,iVec,less<string> > strIVecMap;     // string to iVec
//vectors of maps
typedef vector<strUintMap> suMapVec;
typedef vector<strIntMap> siMapVec;
typedef vector<strUVecMap> suvMapVec;
typedef vector<strIVecMap> sivMapVec;
///////////////////////////////////////////////////////////////////////////////////////////////////
const unsigned maxWordSize_=20;
const unsigned maxSeqLength_=1048576;
///////////////////////////////////////////////////////////////////////////////////////////
// class cjWCount- a word counting class
// variables:
//  mN_: max word length; 
//  nS_: number of seqs included in the count
//  nW_: (vector of unsigned) nW_[i] = number of words of length i+1
//  wc_: the counts, a vector of maps from string to integer different maps for each size word
//  ds_: true if we are counting double stranded
//
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
class cjWCount {
 protected:
  unsigned mN_,nS_; // mN_ = max word length; nS_= # of seqs counted
  uVec nW_;         // nW_[i] = # of words of length i+1
  siMapVec wc_;     // the counts, a vector of maps from string to int
  bool ds_;         // ds_ = true if we are counting double stranded
 public:
  cjWCount():mN_(),nS_(),nW_(),wc_(),ds_() {}
  cjWCount(const cjWCount& o):mN_(o.mN_),nS_(o.nS_),nW_(o.nW_),wc_(o.wc_),ds_(o.ds_) {}
  cjWCount(unsigned maxN,bool ds=false);
  cjWCount(const strVec& seqs,unsigned maxN,bool ds=false);
  bool operator==(const cjWCount& o) const {
    return ((mN_==o.mN_)&&(nS_==o.nS_)&&(nW_==o.nW_)&&(wc_==o.wc_)&&(ds_==o.ds_));
  }
  cjWCount& operator=(const cjWCount& o);
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
  double Expect(const string& w, unsigned o) const;
  double ExpectFreq(const string& w, unsigned o) const;
  string GenerateSequence(unsigned N,unsigned order) const;
  strVec ListWords() const;
  string ListWordCounts() const;
  unsigned LoadCountsFromFile(const string& fn);
};
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
// class cjPWCount- derived from cjWCount
//   word counting class that also maintains position information
//    new variables:
//      wp_: vector of maps from string to vector of ints
//      minSeqL_: min length of analyzed sequences
//      maxSeqL_: max length of analyzed sequences
/////////////////////////////////////////////////////////////////////////////////////////////////
class cjPWCount:public cjWCount {
 protected:
  sivMapVec wp_; 
  unsigned minSeqL_,maxSeqL_;
 public:
  cjPWCount():cjWCount(),wp_(),minSeqL_(),maxSeqL_() {}
  cjPWCount(const cjPWCount& o):cjWCount(o),wp_(o.wp_),minSeqL_(o.minSeqL_),maxSeqL_(o.maxSeqL_) {}
  cjPWCount(unsigned maxN,bool ds=false);
  cjPWCount(const strVec& seqs,unsigned maxN,bool ds=false);
  void AddSequence(const string& seq);
  bool operator==(const cjPWCount& o) const {
    return ((mN_==o.mN_)&&(nS_==o.nS_)&&(nW_==o.nW_)&&(wc_==o.wc_)&&(wp_==o.wp_));
  }
  cjPWCount& operator=(const cjPWCount& o);
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
/*
class ClusterProfile {
 protected:
  int N_; string s_; dVec bc_; iMat wc_; iVec s1_,s2_,sz_;
 public:
  ClusterProfile():N_(),s_(),bc_(),wc_(),s1_(),s2_(),sz_() {}
  ClusterProfile(const ClusterProfile& o):N_(o.N_),s_(o.s_),bc_(o.bc_),wc_(o.wc_),s1_(o.s1_),s2_(o.s2_),sz_(o.sz_) {}
  ClusterProfile(int n,const string& s,const cjWCount& wc);
  double pearsCorr(const ClusterProfile& o) const;// { return Pearson()(bc_,o.bc_); }
  double spearCorr(const ClusterProfile& o) const;// { return Spearman()(bc_,o.bc_); }
  dVec pVec() const { return dVec(bc_); }
  iVec startVec() const { return iVec(s1_); }
  iVec stopVec() const { return iVec(s2_); }
  iVec sizeVec() const { return iVec(sz_); }
  iMat cMat() const { return iMat(wc_); }
  }; */
///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////
