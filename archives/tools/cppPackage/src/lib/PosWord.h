//************************************************************************************************
#ifndef __posword_h
#define __posword_h
//************************************************************************************************
#include <vector>
#include <map>
typedef std::vector<unsigned> uVec;
typedef std::vector<std::string> strVec;
typedef std::vector<double> dVec;
typedef std::vector<dVec> dMat;
typedef std::vector<int> iVec;
typedef std::map<std::string,int,std::less<std::string> > strIntMap;
//************************************************************************************************
/// utility structure posCountVec
struct posCountVec {
  unsigned totalCount; 
  uVec pCounts;
  posCountVec():totalCount(),pCounts() {}
  posCountVec(int n):totalCount(0),pCounts(n,0) {}
  unsigned operator[](unsigned n) { return pCounts.at(n); }
  unsigned Total() const { return totalCount; }
};
//************************************************************************************************
/// useful typedefs
typedef std::map<std::string,posCountVec,std::less<std::string> > WordPosCounts;
typedef std::vector<WordPosCounts> WordPosCountVector;
//************************************************************************************************
/// class WordPositionCounts
class WordPositionCounts {
 private:
  unsigned maxK_,nSeqs_;
  int numP_;
  bool doubleStranded_;
  WordPosCountVector wpc_;
  WordPositionCounts() {}
 public:
  friend class WordWindowCount;
  WordPositionCounts(const WordPositionCounts& o):numP_(o.numP_),maxK_(o.maxK_),
    nSeqs_(o.nSeqs_),wpc_(o.wpc_) {}
  WordPositionCounts(unsigned K,int N,bool ds,const strVec& seqs=strVec());
  //
  int AddSequence(const strVec& seqs);
  int AddSequence(const std::string& seq);
  unsigned Count(const std::string& w,int m=0) const;
  unsigned Count(const std::string& w, int low, int high,int m=0) const;
  uVec WordDist(const std::string& w,int m=0) const; 
  //
  bool DSCount()               const { return doubleStranded_; }
  unsigned NPositionsCounted() const { return numP_; }
  unsigned MaxWordLength()     const { return maxK_; }
  unsigned SequencesCounted()  const { return nSeqs_; }
};
//************************************************************************************************
/// struct EVRec- a utility class designed to calculate estimates of E and V (from Schbath works)
//                given sequence w, order m, and counts of the m- and m+1-length sub words
//                (mWords are not used if m=0)
struct EVRec {
  double E,V; bool good;
  EVRec(const std::string& w,int m,const strIntMap& mp1Words,const strIntMap& mWords);
  bool Okay() const { return good; }
};
//************************************************************************************************
/// class WordWindowCount- uses information in WordPositionCounts to build a windowed version, 
//   including statistical analyses a la Schbath et al
class WordWindowCount {
 private:
  int winL_,offset_,firstP_,maxP_,k_,m_;
  WordPosCounts kCounts_,mCounts_,mp1Counts_;
  const WordPositionCounts& wpc_;
  iVec bounds_; bool good_,smooth_; double sig_;
 public:
  WordWindowCount(const WordWindowCount& o):winL_(o.winL_),offset_(o.offset_),firstP_(o.firstP_),maxP_(o.maxP_),
    k_(o.k_),m_(o.m_),kCounts_(o.kCounts_),mCounts_(o.mCounts_),mp1Counts_(o.mp1Counts_),wpc_(o.wpc_),good_(o.good_),
    smooth_(o.smooth_),sig_(o.sig_){}
  WordWindowCount(const WordPositionCounts& wpc,int k,int m,int win,int max,int first=0,int offset=0,bool smooth=false,double sig=3.0);
  //
  dVec ZVector(const std::string& w) const;
  dVec CountVector(const std::string& w) const;
  double ChiSqr(std::string& w) const;
  double SSqr(std::string& w) const;
  bool Smooth(double sig);
  //
  iVec Boundaries()     const { return iVec(bounds_); }
  int WindowWidth()     const { return winL_;         }
  int StartPosition()   const { return firstP_;       }
  int MaxProfileWidth() const { return maxP_;         }
  int WordSize()        const { return k_;            }
  int ModelOrder()      const { return m_;            }
  bool Okay()           const { return good_;         }
};
//************************************************************************************************
//************************************************************************************************
/// class WordWindowFullCount- uses information in WordPositionCounts to build a windowed version, 
//   including statistical analyses a la Schbath et al, same as above, except all words up to 
//   size k are counted
class WordWindowFullCount {
 private:
  int winL_,offset_,firstP_,maxP_,k_;
  const WordPositionCounts& wpc_;
  iVec bounds_; bool good_,smooth_; double sig_; 
 public:
  WordWindowFullCount(const WordWindowFullCount& o):winL_(o.winL_),offset_(o.offset_),firstP_(o.firstP_),maxP_(o.maxP_),
    k_(o.k_),wpc_(o.wpc_),good_(o.good_),smooth_(o.smooth_),sig_(o.sig_) {}
  WordWindowFullCount(const WordPositionCounts& wpc,int k,int win,int first=0,int offset=0,bool smooth=false,double sig=3.0);
  //
  dVec ZVector(const std::string& w,int m) const;
  dVec ZSVector(const std::string& w,int m) const;
  dVec CountVector(const std::string& w) const;
  double ChiSqr(std::string& w, int m, dVec& Z) const;
  double ZSSqr(std::string& w, int m, dVec& Z) const;
  double SSqr(std::string& w) const;
  double Count(std::string& w) const { return wpc_.Count(w); }
  bool Smooth(double sig);
  //
  iVec Boundaries()     const { return iVec(bounds_); }
  int WindowWidth()     const { return winL_;         }
  int StartPosition()   const { return firstP_;       }
  int MaxProfileWidth() const { return maxP_;         }
  int WordSize()        const { return k_;            }
  bool Okay()           const { return good_;         }
};
//************************************************************************************************
#endif
//************************************************************************************************
