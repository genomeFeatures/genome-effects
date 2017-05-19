//************************************************************************************************
#include <math.h>
#ifdef INTEL
#include <mathimf.h>
#endif
#include <iostream>
#include <fstream>
#include <algorithm>
#include "../lib/PosWord.h"
#include "../lib/logFile.h"
#include "../lib/basicSeq.h"
#include "../lib/util.h"
//************************************************************************************************
const std::string amb("ambig"); 
const std::string val("valid");
const int MaxK= 10, MaxL=100, MaxP=1000;
dMat SmoothVectors; double lastSigSmooth= -1.0;
//************************************************************************************************
///
WordPositionCounts::WordPositionCounts(unsigned K,int N,bool ds,const strVec& seqs):
maxK_(K),numP_(N),doubleStranded_(ds),wpc_() {
  WriteToLog()(std::string("WordPositionCounts::constructor"));
  WriteIntToLog()(std::string("word length= "),maxK_);
  WriteIntToLog()(std::string("prof length= "),numP_);
  WriteIntToLog()(std::string("num of seqs= "),seqs.size());
  // need to put sanity checks on variables here, but not today 17 dec 2003
  // make the vector of map-counts, first index is the length of the word
  wpc_= WordPosCountVector(maxK_,WordPosCounts());
  // now for each word length, put all words, with zeroed elements into the count map
  std::cout  << "building structure..."; std::cout.flush();
  for(int i=0;i<maxK_;++i) {
    int nWords= pow(4,i+1);
    for(int j=0;j<nWords;++j) 
      wpc_[i][std::string(IdxToSeq()(j,i+1))]= posCountVec(numP_);
    // also add an element for ambiguous words (wildcards) and one for total non-ambiguous words
    wpc_[i][amb]= posCountVec(numP_);
    wpc_[i][val]= posCountVec(numP_);
  }
  std::cout << "done. "; std::cout.flush();
  if(!seqs.empty()) AddSequence(seqs);
  return;
}
//************************************************************************************************
/// 
int WordPositionCounts::AddSequence(const strVec& seqs) {
  int update= seqs.size()/500; if(update<10) update=10;
  int updateL=update*100;
  std::cout << "Adding " << seqs.size() << " sequences, Each . equals " << update << " sequences processed"<< std::endl;
  for(int i=0;i<seqs.size();++i) {
    if(!((i+1)%update)) { 
      std::cout << ".";
      if((i+1)%updateL)  std::cout.flush();
      else               std::cout << " " << (i+1) << std::endl; 
    }
    int L=seqs[i].length();
    if(L<numP_) {
      WriteIntToLog()(std::string("short sequence, number "),i);
      WriteIntToLog()(std::string("short sequence, length "),L);
    } else if(L>numP_) {
      WriteIntToLog()(std::string("long sequence, number "),i);
      WriteIntToLog()(std::string("long sequence, length "),L);
    }
    for(int p=0;p<L && p<numP_;++p) {
      int stop1=L-p,stop2=numP_-p;
      for(int k=0;k<maxK_ && k<stop1 && k<stop2;++k) {
	std::string w(seqs[i].substr(p,k+1));
	if(ValidDNASequence()(w) || ValidRNASequence()(w)) {
	  WordPosCounts::iterator it=wpc_[k].find(w);
	  ++(it->second.pCounts[p]);	  ++(it->second.totalCount);
	  //	  ++wpc_[k][w].pCounts[p];	  ++wpc_[k][w].totalCount;
	  it=wpc_[k].find(val);
	  ++(it->second.pCounts[p]);     ++(it->second.totalCount);
	  //	  ++wpc_[k][val].pCounts[p];	  ++wpc_[k][val].totalCount;
	} else {
	  WordPosCounts::iterator it=wpc_[k].find(amb);
	  ++(it->second.pCounts[p]);     ++(it->second.totalCount);
	  //	  ++wpc_[k][amb].pCounts[p];	  ++wpc_[k][amb].totalCount;
	}
      }
    }
    ++nSeqs_;
  }
  return seqs.size();
}
//************************************************************************************************
/// 
int WordPositionCounts::AddSequence(const std::string& seq) {
  int L=seq.length();
  if(L<numP_) {
    WriteIntToLog()(std::string("short sequence, length "),L);
  } else if(L>numP_) {
    WriteIntToLog()(std::string("long sequence, length "),L);
  }
  for(int p=0;p<L && p<numP_;++p)
    for(int k=0;k<maxK_ && (p+k)<L && (p+k)<numP_;++k) {
      std::string w(seq.substr(p,k+1));
      if(ValidDNASequence()(w) || ValidRNASequence()(w)) {
	++wpc_[k][w].pCounts[p];	++wpc_[k][w].totalCount;
	++wpc_[k][val].pCounts[p];	++wpc_[k][val].totalCount;
      } else {
	++wpc_[k][amb].pCounts[p];	++wpc_[k][amb].totalCount;
      }
    }
  ++nSeqs_;
  return 1;
}
//************************************************************************************************
/// 
unsigned WordPositionCounts::Count(const std::string& w,int m) const {
  if((w.length()>maxK_ && w!=val && w!=amb) || !w.length()) {
    WriteToLogs()(std::string("incorrect length word passed to WordPositionCounts::Count ")+w);
    return 0;
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w) && w!=val && w!=amb) {
    WriteToLogs()(std::string("invalid sequence word passed to WordPositionCounts::Count ")+w);
    return 0;
  }
  WordPosCounts::const_iterator w_iter;
  if(w==val || w==amb) {
    if(m<1 || m>(maxK_+1)) {
      WriteToLogs()(std::string("invalid m passed to WordPositionCounts::Count ")+w);
      return 0;
    }
    w_iter= wpc_[m-1].find(w);
  } else
    w_iter= wpc_[w.length()-1].find(w);
  return (*w_iter).second.Total();
}
//************************************************************************************************
/// 
unsigned WordPositionCounts::Count(const std::string&w, int low, int high,int m) const {
  if((w.length()>maxK_ && w!=val && w!=amb) || !w.length()) {
    WriteToLogs()(std::string("incorrect length word passed to WordPositionCounts::Count ")+w);
    return 0;
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w) && w!=val && w!=amb) {
    WriteToLogs()(std::string("invalid sequence word passed to WordPositionCounts::Count ")+w);
    return 0;
  }
  unsigned count=0;
  int start=low; if(start<0) start=0;
  int stop=high; if(stop >numP_) stop=numP_;
  WordPosCounts::const_iterator w_iter;
  if(w==val || w==amb) {
    if(m<1 || m>(maxK_+1)) {
      WriteToLogs()(std::string("invalid m passed to WordPositionCounts::Count ")+w);
      return 0;
    }
    w_iter= wpc_[m-1].find(w);
  } else
    w_iter= wpc_[w.length()-1].find(w);
  for(int p=start;p<stop;++p)
    count += (*w_iter).second.pCounts[p];
  return count;
}
//************************************************************************************************
/// 
uVec WordPositionCounts::WordDist(const std::string& w, int m) const {
  if((w.length()>maxK_ && w!=val && w!=amb) || !w.length()) {
    WriteToLogs()(std::string("incorrect length word passed to WordPositionCounts::WordDist ")+w);
    return uVec();
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w) && w!=val && w!=amb) {
    WriteToLogs()(std::string("invalid sequence word passed to WordPositionCounts::WordDist ")+w);
    return uVec();
  }
  WordPosCounts::const_iterator w_iter;
  if(w==val || w==amb) {
    if(m<1 || m>(maxK_+1)) {
      WriteToLogs()(std::string("invalid m passed to WordPositionCounts::AmbigDist ")+w);
      return uVec();
    }
    w_iter= wpc_[m-1].find(w);
  } else
    w_iter= wpc_[w.length()-1].find(w);
  return uVec((*w_iter).second.pCounts);
}
//************************************************************************************************
EVRec::EVRec(const std::string& w,int m,const strIntMap& mp1Words,const strIntMap& mWords) {
  //  std::ofstream dumpFile("EVdump.txt",std::ios::app);
  int k= w.length(); E=V= 0.0; good=false;
  if(!k) { WriteToLogs()(std::string("empty seq passed to EVRec")); return; }
  if(m<0 || m>(k-2)) { WriteIntToLogs()(std::string("invalid order m passed to EVRec"),m); return; } 
  // now check to make sure that all the necessary words are in mWords and mp1Words
  for(int i=0;i<k-m;++i)
    if(!mp1Words.count(w.substr(i,m+1))) {
      WriteToLogs()(std::string("EVRec: mp1Words missing subword: ")+w.substr(i,m+1));
      return;
    }
  if(m) 
    for(int i=0;i<k-m+1;++i)
      if(!mWords.count(w.substr(i,m))) {
	WriteToLogs()(std::string("EVRec: mWords missing subword: ")+w.substr(i,m));
	return;
      }
  good=true;
  // okay- data checks out okay, go ahead and do the calculation
  E=1.0;
  for(int i=0;i<k-m;++i) {
    strIntMap::const_iterator mp1w=mp1Words.find(w.substr(i,m+1));
    //    dumpFile << "multiplying E (" << E << ") by word " << w.substr(i,m+1) << " -> " << mp1w->second << '\n';
    E*= static_cast<double>(mp1w->second);
    if(m && i) {
      strIntMap::const_iterator mw=mWords.find(w.substr(i,m));
      //      dumpFile << "dividing E (" << E << ") by word " << w.substr(i,m) << " -> " << mw->second << '\n';
      E/= static_cast<double>(mw->second);
    }
  }
  // special case- m=0: we need the total (valid) counts
  if(!m) {
    strIntMap::const_iterator vw=mp1Words.find(val);
    double Nval= vw->second;
    for(int i=1;i<k;++i) E /= Nval;
    //    dumpFile << "m=0, E= " << E << '\n';
    strIntMap mp1SubWords= SubWords()(w,m+1);
    double mSum=0.0;
    for(strIntMap::iterator a=mp1SubWords.begin();a!=mp1SubWords.end();++a) {
      strIntMap::const_iterator mp1w= mp1Words.find(a->first);
      mSum += static_cast<double>(a->second * a->second)/static_cast<double>(mp1w->second);
      //      dumpFile << " msum term: " << (*a).first << ", n=" << a->second << ", N=" << mp1w->second << '\n';  
    }
    V=E - E*E*mSum;
    //    dumpFile << " m=0, V= " << V << "\n\n";
    return;
  }
  //  dumpFile << "E = " << E << "\n\n";
  // next the variance calculation
  V=E;
  //  dumpFile << "variance starts equal to E-> " << V << '\n';
  // first the sum on self-intersection terms
  double sum2=0.0;
  for(int d=1;d<(k-m);++d)
    if(w.substr(0,k-d)==w.substr(d,k-d)) {
      double extE=E;
      std::string extW= w.substr(0,d) + w;
      for(int i=(k-m);i<(k+d-m);++i) {
	strIntMap::const_iterator mp1w= mp1Words.find(extW.substr(i,m+1));
	extE *= static_cast<double>(mp1w->second);
	strIntMap::const_iterator mw= mWords.find(extW.substr(i,m));
	extE /= static_cast<double>(mw->second);
      }
      sum2+= extE;
    }
  // now add to variance estimator
  V += 2.0*sum2;
  //  dumpFile << "variance second term, self interaction " << (2.0*sum2) << ", V= " << V << '\n';
  // next the sums over subwords of w
  // m+1 words first
  double sum3=0.0,sum4=0.0;
  strIntMap mp1SubWords= SubWords()(w,m+1);
  for(strIntMap::iterator a=mp1SubWords.begin();a!=mp1SubWords.end();++a) {
    strIntMap::const_iterator mp1w= mp1Words.find((*a).first);
    sum3 += static_cast<double>(a->second * a->second) / static_cast<double>(mp1w->second);
    //    dumpFile << " sum3 term: " << (*a).first << ", n=" << a->second << ", N=" << mp1w->second << '\n';  
  }
  // m words second (if m>0)
  strIntMap mSubWords= SubWords()(w,m);
  std::string last(w.substr(k-m,m));
  if(mSubWords[last]==1) mSubWords.erase(last);
  else                   --mSubWords[last];
  // this gives us counts of the subwords inside of w, not including the final m-word
  for(strIntMap::iterator a=mSubWords.begin();a!=mSubWords.end();++a) {
    strIntMap::const_iterator mw= mWords.find(a->first);
    sum4+= static_cast<double>(a->second * a->second) / static_cast<double>(mw->second);
    //    dumpFile << " sum4 term: " << a->first << ", n=" << a->second << ", N=" << mw->second << '\n';  
  }
  //  dumpFile << "variance 3rd sum= " << sum3 << "\nvariance 4th sum= " << sum4 << '\n';
  // now include the first word term
  std::string first(w.substr(0,m));
  strIntMap::const_iterator mw= mWords.find(first);
  //  dumpFile << " sum4 extra term-> n=" << mSubWords[first] << ", N=" << mw->second << '\n'; 
  sum4+= (1.0-2.0*static_cast<double>(mSubWords[first])) / static_cast<double>(mw->second);
  //  dumpFile << "variance 4th sum, with first word term-> " << sum4 << '\n';
  // now add to variance estimator
  V += E*E * (sum4 - sum3);
  //  dumpFile << "Final variance= " << V << "\n\n";
  //  dumpFile.close();
}
//************************************************************************************************
/*
class WordWindowCount {
 private:
  int winL_,offset_,firstP_,maxP_,k_,m_;
  WordPosCounts kCounts_,mCounts_,mp1Counts_;
  WordPositionCounts& wpc_;
  iVec bounds_;
 public:
  WordWindowCount(const WordWindowCount& o):winL_(o.winL_),offset_(o.offset_),firstP_(o.firstP_),maxP_(o.maxP_),
    k_(o.k_),m_(o.m_),kCounts_(o.kCounts_),mCounts_(o.mCounts_),mp1Counts_(o.mp1Counts_),wpc_(o.wpc_) {}
  WordWindowCount(const WordPositionCounts& wpc,int k,int m,int win,int max,int first=0,int offset=0);
  //
  dVec ZVector(const std::string& w) const;
  dVec CountVector(const std::string& w) const;
  double ChiSqr(std::string& w) const;
  double SSqr(std::string& w) const;
  //
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////
void GenerateSmoothVectors(dMat& sv,double smSig) {
  sv.clear();
  int halfWindow= static_cast<int>(round(4.0*smSig));
  int window= 2*halfWindow+1;
  sv= dMat(2,dVec(window,1.0));
  if(window==1) return;
  for(int i=-halfWindow;i<=halfWindow;++i) {
    int idx= i+halfWindow;
    sv[0][idx]= exp(-0.5*static_cast<double>(i*i)/(smSig*smSig));
    if(i == -halfWindow)
      sv[1][idx] = sv[0][idx];
    else
      sv[1][idx] = sv[1][idx-1]+sv[0][idx];
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void SmoothVec(const dMat& sv, dVec& v) {
  if(sv.size()!=2) return;
  int w2= sv[0].size()/2,L=v.size();
  if(w2<1) return;
  dVec nv(L,0.0);
  for(int i=0;i<v.size();++i) {
    nv[i]=0.0; int nTerms=0;
    for(int j= i-w2;j<=(i+w2) && j<L;++j) 
      if(j>=0) {
	++nTerms;
	nv[i]+= v[j]*sv[0][j-i+w2];
      }
    nv[i] /= sv[1][nTerms-1];
  }
  v=nv;
}    
//************************************************************************************************
/// class WordWindowCount constructor
WordWindowCount::WordWindowCount(const WordPositionCounts& wpc,int k,int m,int win,int max,int first,int offset,bool smooth,double sig):
wpc_(wpc),k_(k),m_(m),winL_(win),maxP_(max),firstP_(first),offset_(offset),kCounts_(),mCounts_(),mp1Counts_(),bounds_(),good_(false),
smooth_(smooth),sig_(sig) {
  // reality checks first
  if(k_<=0)     { 
    WriteIntToLogs()(std::string("negative or zero k in WordWindowCount constructor"),k_); return; }
  if(k_>MaxK)   { 
    WriteIntToLogs()(std::string("k exceeds max value in WordWindowCount constructor, resetting"),k_); k_=MaxK; }
  if(m_<0)      { 
    WriteIntToLogs()(std::string("negative m in WordWindowCount constructor"),m_); return; }
  if(m_>(k_-2)) { 
    WriteIntToLogs()(std::string("m exceeds k-2 in WordWindowCount constructor, resetting to k-2"),m_); m_=k-2; }
  if(winL_<0) {
    WriteIntToLogs()(std::string("negative window length in WordWindowcount constructor"),winL_); return; }
  if(winL_>MaxL){ 
    WriteIntToLogs()(std::string("winL exceeds max value in WordWindowCount constructor, resetting"),winL_); winL_=MaxL; }
  if(offset_>winL_) {
    WriteIntToLogs()(std::string("offset exceeds window size in WordWindowCount constructor, taking modulo"),offset_); offset_= winL_% offset_; }
  if(maxP_>MaxP) {
    WriteIntToLogs()(std::string("max profile width exceeds max value in WordWindowCount constructor, resetting"),maxP_); maxP_=MaxP; }
  // okay- if we made it to here, we're good to go
  good_ = true;
  // now find the boundaries defined by the input values, and set up the assignment array, which will be used for
  // rapid conversion of absolute position to window
  // ** important ** firstP_ is only cosmetic, such that it affects the boundary positions and values, however, the
  //                 counts in wpc are absolute to the sequences, which is to say that the first position is 0
  bounds_.reserve((maxP_-k_+1)/winL_ +1);
  iVec assign(maxP_,0);
  int currentWin=0;  bounds_.push_back(firstP_);
  WriteIntToLog()(std::string("firstP_= "),firstP_);
  WriteIntToLog()(std::string("maxP_= "),maxP_);
  WriteIntToLog()(std::string("offset_= "),offset_);
  for(int p= firstP_; p< firstP_+maxP_; ++p) {
    if(p != firstP_ && p%winL_ == offset_ && (p+k_-firstP_)<maxP_) {
      ++currentWin;
      bounds_.push_back(p);
      WriteIntToLog()(std::string("firstP_+p+k_"),p+k_-firstP_);
      WriteIntToLog()(std::string("boundary set "),p);
    }
    assign[p-firstP_]= currentWin;
  }
  bounds_.push_back(firstP_+maxP_);
  // now, use the positional counts in wpc_ to get the windowed information- remembering that for modeling purposes, the 
  // windows of the m-Words need to be extended to winL_+k_-m_+1 and for m+1-Words to winL_+k_-m_, this means that 
  // m-Words and m+1-Words may contribute to more than one window, but since we're only using them to model the 
  // k-Words based on normal approximation, it should be okay
  // kWords first:
  int nWords= pow(4,k_),nWin=bounds_.size()-1;
  for(int j=0;j<nWords+2;++j) {
    std::string word;
    if(j<nWords)       word= IdxToSeq()(j,k_);
    else if(j==nWords) word= val;
    else               word= amb;
    kCounts_[word]= posCountVec(nWin);
    uVec counts= wpc_.WordDist(word,k_);
    for(int p=0;p<counts.size();++p) {
      kCounts_[word].totalCount += counts[p];
      kCounts_[word].pCounts[assign[p]] += counts[p];
    }
  }
  // next m+1-Words
  nWords= pow(4,m_+1); int winEdge= k_-m_-1;
  for(int j=0;j<nWords+2;++j) {
    std::string word;
    if(j<nWords)       word= IdxToSeq()(j,m_+1);
    else if(j==nWords) word= val;
    else               word= amb;
    mp1Counts_[word]= posCountVec(nWin);
    uVec counts= wpc_.WordDist(word,m_+1);
    for(int p=0;p<counts.size();++p) {
      mp1Counts_[word].totalCount += counts[p];
      mp1Counts_[word].pCounts[assign[p]] += counts[p];
      // now check for overlap
      if(assign[p]!=0 && (p+firstP_ - bounds_[assign[p]])<winEdge) 
	mp1Counts_[word].pCounts[assign[p]-1] += counts[p];
    }
  }
  // finally, if m!=0, m-Words
  if(m) {
    nWords= pow(4,m_); winEdge= k_-m_;
    for(int j=0;j<nWords+2;++j) {
      std::string word;
      if(j<nWords)       word= IdxToSeq()(j,m_);
      else if(j==nWords) word= val;
      else               word= amb;
      mCounts_[word]= posCountVec(nWin);
      uVec counts= wpc_.WordDist(word,m_);
      for(int p=0;p<counts.size();++p) {
	mCounts_[word].totalCount += counts[p];
	mCounts_[word].pCounts[assign[p]] += counts[p];
	// now check for overlap
	if(assign[p]!=0 && (p+firstP_ - bounds_[assign[p]])<winEdge) 
	  mCounts_[word].pCounts[assign[p]-1] += counts[p];
      }
    }
  }
}
//************************************************************************************************
///
dVec WordWindowCount::CountVector(const std::string& w) const {
  dVec c;
  if(w.length()!=k_ && w.length()!=m_ && w.length()!=(m_+1)) {
    WriteToLogs()(std::string("incorrect length word w in WordWindowCount::CountVector")); return c; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowCount::CountVector")); return c; 
  }
  int nWin=bounds_.size()-1; c.reserve(nWin);
  if(w.length()==k_) {
    WordPosCounts::const_iterator kw= kCounts_.find(w);
    for(int i=0;i<nWin;++i)
      c.push_back(static_cast<double>(kw->second.pCounts[i]));
  } else if (w.length()==(m_+1)) {
    WordPosCounts::const_iterator pw= mp1Counts_.find(w);
    for(int i=0;i<nWin;++i)
      c.push_back(static_cast<double>(pw->second.pCounts[i]));
  } else if (w.length()) {
    WordPosCounts::const_iterator pw= mCounts_.find(w);
    for(int i=0;i<nWin;++i)
      c.push_back(static_cast<double>(pw->second.pCounts[i]));
  }    
  if(smooth_) {
    if(sig_ != lastSigSmooth) { 
      GenerateSmoothVectors(SmoothVectors,sig_);
      lastSigSmooth= sig_;
    }
    SmoothVec(SmoothVectors,c);
  }
  return c;
}
//************************************************************************************************
///
dVec WordWindowCount::ZVector(const std::string& w) const {
  //  std::ofstream dumpFile("Zdump.txt",std::ios::app);
  dVec ZVec;
  if(w.length()!=k_) {
    WriteToLogs()(std::string("incorrect length word w in WordWindowCount::ZVector")); return ZVec; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowCount::ZVector")); return ZVec; 
  }
  // get a list of the m+1- and m-words that will be needed, so the substrings only need to be
  // taken once.  Include also the necessary self-overlap words 
  strVec subMp1Words,subMWords;
  for(int p=0;p<(k_-m_+1);++p) {
    if(p<(k_-m_)) {
      std::string wmp1= w.substr(p,m_+1);
      if(std::find(subMp1Words.begin(),subMp1Words.end(),wmp1)==subMp1Words.end()) subMp1Words.push_back(wmp1);
    }
    if(m_) {
      std::string wm= w.substr(p,m_);
      if(std::find(subMWords.begin(),subMWords.end(),wm)==subMWords.end()) subMWords.push_back(wm);
    }
  }
  // now the self-overlap sub-words
  for(int d=1;d<k_-m_;++d)
    if(w.substr(0,k_-d)==w.substr(d,k_-d)) {
      std::string extWord= w.substr(0,d)+w;
      for(int p=(k_-m_);p<d;++p) 
	if(std::find(subMp1Words.begin(),subMp1Words.end(),extWord.substr(p,m_+1))==subMp1Words.end()) 
	  subMp1Words.push_back(extWord.substr(p,m_+1));
      if(m_)
	for(int p=(k_-m_+1);p<d;++p)
	  if(std::find(subMWords.begin(),subMWords.end(),extWord.substr(p,m_))==subMWords.end())
	    subMWords.push_back(extWord.substr(p,m_));
    }
  /*
  // dump the lists of words to the dump file
  dumpFile << "m+1 word list for E,V calculations:\n";
  for(int i=0;i<subMp1Words.size();++i) dumpFile << "  " << subMp1Words[i];
  dumpFile << '\n';
  if(m_) {
    dumpFile << "m word list for E,V calculations:\n";
    for(int i=0;i<subMWords.size();++i) dumpFile << "  " << subMWords[i];
    dumpFile << '\n';
  }
  */
  // loop over all windows for this word
  int nWin= bounds_.size()-1; ZVec.reserve(nWin);
  for(int i=0;i<nWin;++i) {
    strIntMap mWords,mp1Words;
    for(int j=0;j<subMp1Words.size();++j) {
      WordPosCounts::const_iterator p1w= mp1Counts_.find(subMp1Words[j]);
      mp1Words[subMp1Words[j]]= p1w->second.pCounts[i];
      //      dumpFile << " subword-m+1 " << subMp1Words[j] << ", count[" << i << "]= " << p1w->second.pCounts[i] << '\n'; 
    }
    if(m_)
      for(int j=0;j<subMWords.size();++j) {
	WordPosCounts::const_iterator pw= mCounts_.find(subMWords[j]);
	mWords[subMWords[j]]= pw->second.pCounts[i];
	//	dumpFile << " subword-m " << subMWords[j] << ", count[" << i << "]= " << pw->second.pCounts[i] << '\n';
      }
    else {
      WordPosCounts::const_iterator pw= mp1Counts_.find(val);
      mp1Words[val]= pw->second.pCounts[i];
      //      dumpFile << " valid 1-mers [" << i << "]= " << pw->second.pCounts[i] << '\n';
    } 
    EVRec er(w,m_,mp1Words,mWords);
    WordPosCounts::const_iterator kw= kCounts_.find(w);
    double z= (static_cast<double>(kw->second.pCounts[i])-er.E)/sqrt(er.V);
    if(er.V==0.0) z=0.0;
    ZVec.push_back(z);
    //    dumpFile << w << '\t' << bounds_[i] << '\t' << k_ << ',' << m_ << '\t' << (kw->second.pCounts[i]) << '\t' << er.E << '\t' << sqrt(er.V) << " (" << er.V << ")\t" << z << '\n';
  }
  //  dumpFile.close();
  if(smooth_) {
    if(sig_ != lastSigSmooth) { 
      GenerateSmoothVectors(SmoothVectors,sig_);
      lastSigSmooth= sig_;
    }
    SmoothVec(SmoothVectors,ZVec);
  }

  return dVec(ZVec);
}
//************************************************************************************************
///
double  WordWindowCount::SSqr(std::string& w) const {
  std::ofstream dumpFile("Sdump.txt",std::ios::app);
  if(w.length()!=k_) {
    WriteToLogs()(std::string("incorrect length word w in WordWindowCount::SSqr")); return -1.0; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowCount::SSqr")); return -1.0; 
  }
  double nWin= bounds_.size()-3; int n= bounds_.size()-1;
  // get the fractional sizes of the end windows
  double fracWin0= static_cast<double>(bounds_[1]-bounds_[0])/static_cast<double>(winL_);
  double fracWinN= static_cast<double>(bounds_[n]-bounds_[n-1]-k_+1)/static_cast<double>(winL_);
  nWin += fracWin0+fracWinN;
  // calculate the average per window
  WordPosCounts::const_iterator kw= kCounts_.find(w);
  double mean= static_cast<double>(kw->second.Total())/nWin;
  // scaled counts for the end windows
  double scCount0= static_cast<double>(kw->second.pCounts[0])/fracWin0;
  double scCountN= static_cast<double>(kw->second.pCounts[n])/fracWinN;
  // do the sum, starting with the ends, then the remainder
  double sSqr= (scCount0-mean)*(scCount0-mean)+(scCountN-mean)*(scCountN-mean);
  for(int i=1;i<n;++i)
    sSqr += (kw->second.pCounts[i]-mean)*(kw->second.pCounts[i]-mean);
  sSqr /= mean;
  dumpFile << w << '\t' << kw->first << '\t' << kw->second.Total() << '\t' << mean << '\t' << sSqr << '\t' << fracWin0 << '\t' << fracWinN << '\t' << kw->second.pCounts.size();
  dumpFile << '\t' << bounds_[0] << '\t' << bounds_[1] << '\t' << bounds_[n-1] << '\t' << bounds_[n] << '\n';
  dumpFile.close();
  return sSqr;
}
//************************************************************************************************
///
double WordWindowCount::ChiSqr(std::string& w) const {
  //  std::ofstream dumpFile("Zdump.txt",std::ios::app);
  if(w.length()!=k_) {
    WriteToLogs()(std::string("incorrect length word w in WordWindowCount::ChiSqr")); return -1.0; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowCount::ChiSqr")); return -1.0; 
  }
  int n= bounds_.size()-1;
  dVec Z= ZVector(w);
  double chisqr=0.0;
  for(int i=0;i<n;++i)
    chisqr += Z[i]*Z[i];
  //  dumpFile << w << ", X-sqr= " << chisqr << "\n\n";
  //  dumpFile.close();
  return chisqr;
}
//************************************************************************************************
/// 
bool WordWindowCount::Smooth(double sig) {
  if(sig<0.0 || sig>50.0) return false;
  sig_=sig; smooth_=true;
  return true;
}
//************************************************************************************************
/// class WordWindowCount constructor
WordWindowFullCount::WordWindowFullCount(const WordPositionCounts& wpc,int k,int win,int first,int offset,bool smooth,double sig):
wpc_(wpc),k_(k),winL_(win),maxP_(),firstP_(first),offset_(offset),bounds_(),good_(false),smooth_(smooth),sig_(sig) {
  // reality checks first
  if(k_<=0)     { 
    WriteIntToLogs()(std::string("negative or zero k in WordWindowCount constructor"),k_); return; }
  if(k_>MaxK)   { 
    WriteIntToLogs()(std::string("k exceeds max value in WordWindowCount constructor, resetting"),k_); k_=MaxK; }
  if(winL_<0) {
    WriteIntToLogs()(std::string("negative window length in WordWindowcount constructor"),winL_); return; }
  if(winL_>MaxL){ 
    WriteIntToLogs()(std::string("winL exceeds max value in WordWindowCount constructor, resetting"),winL_); winL_=MaxL; }
  if(offset_>winL_) {
    WriteIntToLogs()(std::string("offset exceeds window size in WordWindowCount constructor, taking modulo"),offset_); offset_= winL_% offset_; }
  // okay- if we made it to here, we're good to go
  good_ = true;
  maxP_= wpc.NPositionsCounted();
  // now find the boundaries defined by the input values, and set up the assignment array, which will be used for
  // rapid conversion of absolute position to window
  // ** important ** firstP_ is only cosmetic, such that it affects the boundary positions and values, however, the
  //                 counts in wpc are absolute to the sequences, which is to say that the first position is 0
  bounds_.reserve((maxP_-k_+1)/winL_ +1);
  iVec assign(maxP_,0);
  int currentWin=0;  bounds_.push_back(firstP_);
  //  WriteIntToLog()(std::string("WordWindowFullCount constructor: firstP_= "),firstP_);
  //  WriteIntToLog()(std::string("WordWindowFullCount constructor: maxP_= "),maxP_);
  //  WriteIntToLog()(std::string("WordWindowFullCount constructor: offset_= "),offset_);
  for(int p= firstP_; p< firstP_+maxP_; ++p) {
    if(p != firstP_ && p%winL_ == offset_ && (p+k_-firstP_-1)<maxP_) {
      ++currentWin;
      bounds_.push_back(p);
      //      WriteIntToLog()(std::string("WordWindowFullCount constructor: boundary set "),p);
    }
    assign[p-firstP_]= currentWin;
  }
  bounds_.push_back(firstP_+maxP_);
  //  WriteIntToLog()(std::string("WordWindowFullCount constructor: boundary set "),firstP_+maxP_);
  //  WriteIntToLog()(std::string("WordWindowFullCount constructor: effective last boundary "),firstP_+maxP_-k_+1);  
  /*
  // now, use the positional counts in wpc_ to get the windowed information- remembering that for modeling purposes, the 
  // windows of the m-Words need to be extended to winL_+k_-m_+1 and for m+1-Words to winL_+k_-m_, this means that 
  // m-Words and m+1-Words may contribute to more than one window, but since we're only using them to model the 
  // k-Words based on normal approximation, it should be okay
  // kWords first:
  int nWords= pow(4,k_),nWin=bounds_.size()-1;
  for(int j=0;j<nWords+2;++j) {
    std::string word;
    if(j<nWords)       word= IdxToSeq()(j,k_);
    else if(j==nWords) word= val;
    else               word= amb;
    kCounts_[word]= posCountVec(nWin);
    uVec counts= wpc_.WordDist(word,k_);
    for(int p=0;p<counts.size();++p) {
      kCounts_[word].totalCount += counts[p];
      kCounts_[word].pCounts[assign[p]] += counts[p];
    }
  }
  // next m+1-Words
  nWords= pow(4,m_+1); int winEdge= k_-m_-1;
  for(int j=0;j<nWords+2;++j) {
    std::string word;
    if(j<nWords)       word= IdxToSeq()(j,m_+1);
    else if(j==nWords) word= val;
    else               word= amb;
    mp1Counts_[word]= posCountVec(nWin);
    uVec counts= wpc_.WordDist(word,m_+1);
    for(int p=0;p<counts.size();++p) {
      mp1Counts_[word].totalCount += counts[p];
      mp1Counts_[word].pCounts[assign[p]] += counts[p];
      // now check for overlap
      if(assign[p]!=0 && (p+firstP_ - bounds_[assign[p]])<winEdge) 
	mp1Counts_[word].pCounts[assign[p]-1] += counts[p];
    }
  }
  // finally, if m!=0, m-Words
  if(m) {
    nWords= pow(4,m_); winEdge= k_-m_;
    for(int j=0;j<nWords+2;++j) {
      std::string word;
      if(j<nWords)       word= IdxToSeq()(j,m_);
      else if(j==nWords) word= val;
      else               word= amb;
      mCounts_[word]= posCountVec(nWin);
      uVec counts= wpc_.WordDist(word,m_);
      for(int p=0;p<counts.size();++p) {
	mCounts_[word].totalCount += counts[p];
	mCounts_[word].pCounts[assign[p]] += counts[p];
	// now check for overlap
	if(assign[p]!=0 && (p+firstP_ - bounds_[assign[p]])<winEdge) 
	  mCounts_[word].pCounts[assign[p]-1] += counts[p];
      }
    }
  }
  */
}
//************************************************************************************************
///
dVec WordWindowFullCount::CountVector(const std::string& w) const {
  dVec c;
  if(w.length()>k_ || !w.length()) {
    WriteIntToLogs()(std::string("incorrect length word w in WordWindowFullCount::CountVector-> ")+w,w.length()); return c; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowFullCount::CountVector -> ")+w); return c; 
  }
  int nWin=bounds_.size()-1; c.reserve(nWin);
  for(int i=0;i<nWin;++i) 
    c.push_back(static_cast<double>(wpc_.Count(w,bounds_[i]-firstP_,bounds_[i+1]-firstP_)));
  //
  if(smooth_) {
    if(sig_ != lastSigSmooth) { 
      GenerateSmoothVectors(SmoothVectors,sig_);
      lastSigSmooth= sig_;
    }
    SmoothVec(SmoothVectors,c);
  }
  return c;
}
//************************************************************************************************
///
double WordWindowFullCount::SSqr(std::string& w) const {
  if(w.length()>k_ || !w.length()) {
    WriteIntToLogs()(std::string("incorrect length word w in WordWindowFullCount::SSqr-> ")+w,w.length()); return -1.0; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowFullCount::SSqr -> ")+w); return -1.0; 
  }
  //  std::ofstream dumpFile("Sdump.txt",std::ios::app);
  if(!wpc_.Count(w)) return 0.0;
  double nWin= bounds_.size()-3; int n= bounds_.size()-1;
  // get the fractional sizes of the end windows
  double fracWin0= static_cast<double>(bounds_[1]-bounds_[0])/static_cast<double>(winL_);
  double fracWinN= static_cast<double>(bounds_[n]-bounds_[n-1]-k_+1)/static_cast<double>(winL_);
  nWin += fracWin0+fracWinN;
  // calculate the average per window
  dVec counts= CountVector(w);
  double mean= 0.0; for(int i=0;i<counts.size();++i) mean +=counts[i];
  mean /= nWin;
  // scaled counts for the end windows
  double scCount0= static_cast<double>(counts[0])/fracWin0;
  double scCountN= static_cast<double>(counts[n])/fracWinN;
  // do the sum, starting with the ends, then the remainder
  double sSqr= (scCount0-mean)*(scCount0-mean)+(scCountN-mean)*(scCountN-mean);
  for(int i=1;i<n;++i)
    sSqr += (counts[i]-mean)*(counts[i]-mean);
  sSqr /= mean;
  //  dumpFile << w << '\t' << kw->first << '\t' << kw->second.Total() << '\t' << mean << '\t' << sSqr << '\t' << fracWin0 << '\t' << fracWinN << '\t' << kw->second.pCounts.size();
  //  dumpFile << '\t' << bounds_[0] << '\t' << bounds_[1] << '\t' << bounds_[n-1] << '\t' << bounds_[n] << '\n';
  //  dumpFile.close();
  return sSqr;
}
//************************************************************************************************
///
dVec WordWindowFullCount::ZVector(const std::string& w,int m) const {
  //  std::ofstream dumpFile("Zdump.txt",std::ios::app);
  dVec ZVec;
  if(w.length()>k_ || w.length()<2) {
    WriteToLogs()(std::string("incorrect length word w in WordWindowFullCount::ZVector-> ")+w); return ZVec; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowFullCount::ZVector-> ")+w); return ZVec; 
  }
  if(m<0 || m>(w.length()-2)) {
    WriteIntToLogs()(std::string("invalid order m in WordWindowFullCount::ZVector-> "),m); return ZVec;
  }
  // get a list of the m+1- and m-words that will be needed, so the substrings only need to be
  // taken once.  Include also the necessary self-overlap words 
  strVec subMp1Words,subMWords;
  for(int p=0;p<(k_-m+1);++p) {
    if(p<(k_-m)) {
      std::string wmp1= w.substr(p,m+1);
      if(std::find(subMp1Words.begin(),subMp1Words.end(),wmp1)==subMp1Words.end()) subMp1Words.push_back(wmp1);
    }
    if(m) {
      std::string wm= w.substr(p,m);
      if(std::find(subMWords.begin(),subMWords.end(),wm)==subMWords.end()) subMWords.push_back(wm);
    }
  }
  // now the self-overlap sub-words
  for(int d=1;d<k_-m;++d)
    if(w.substr(0,k_-d)==w.substr(d,k_-d)) {
      std::string extWord= w.substr(0,d)+w;
      for(int p=(k_-m);p<d;++p) 
	if(std::find(subMp1Words.begin(),subMp1Words.end(),extWord.substr(p,m+1))==subMp1Words.end()) 
	  subMp1Words.push_back(extWord.substr(p,m+1));
      if(m)
	for(int p=(k_-m+1);p<d;++p)
	  if(std::find(subMWords.begin(),subMWords.end(),extWord.substr(p,m))==subMWords.end())
	    subMWords.push_back(extWord.substr(p,m));
    }
  // loop over all windows for this word;  the only tricky thing is expanding the windows rightward for the
  //  m- and m+1-subwords. since the window is being defined for the word w.  For the Z calculation to make sense,
  //  we need to include all valid subwords that would be contained within the last word starting in each window
  // For instance, if we are counting hexamers, then for m= k-2=4, we need an expansion of the window by one for
  //  the hexamers, and by two for the tetramers.
  int nWin= bounds_.size()-1; ZVec.reserve(nWin);
  for(int i=0;i<nWin;++i) {
    int leftBound= bounds_[i],rightBound= bounds_[i+1];
    int mRightBound= rightBound+k_-m, mp1RightBound= rightBound+k_-m-1;
    strIntMap mWords,mp1Words;
    for(int j=0;j<subMp1Words.size();++j) {
      mp1Words[subMp1Words[j]]= wpc_.Count(subMp1Words[j],leftBound-firstP_,mp1RightBound-firstP_);
      //      dumpFile << " subword-m+1 " << subMp1Words[j] << ", count[" << i << "]= " << p1w->second.pCounts[i] << '\n'; 
    }
    if(m)
      for(int j=0;j<subMWords.size();++j) {
	mWords[subMWords[j]]= wpc_.Count(subMWords[j],leftBound-firstP_,mRightBound-firstP_);
	//	dumpFile << " subword-m " << subMWords[j] << ", count[" << i << "]= " << pw->second.pCounts[i] << '\n';
      }
    else {
      mp1Words[val]= wpc_.Count(val,leftBound-firstP_,mp1RightBound-firstP_,m+1);
      //      dumpFile << " valid 1-mers [" << i << "]= " << pw->second.pCounts[i] << '\n';
    } 
    EVRec er(w,m,mp1Words,mWords);
    double z= (static_cast<double>(wpc_.Count(w,leftBound-firstP_,rightBound-firstP_))-er.E)/sqrt(er.V);
    if(er.V==0.0) z=0.0;
    ZVec.push_back(z);
    //    dumpFile << w << '\t' << bounds_[i] << '\t' << k_ << ',' << m << '\t' << (kw->second.pCounts[i]) << '\t' << er.E << '\t' << sqrt(er.V) << " (" << er.V << ")\t" << z << '\n';
  }
  //  dumpFile.close();
  if(smooth_) {
    if(sig_ != lastSigSmooth) { 
      GenerateSmoothVectors(SmoothVectors,sig_);
      lastSigSmooth= sig_;
    }
    SmoothVec(SmoothVectors,ZVec);
  }
  return dVec(ZVec);
}
//************************************************************************************************
///
dVec WordWindowFullCount::ZSVector(const std::string& w,int m) const {
  //  std::ofstream dumpFile("Zdump.txt",std::ios::app);
  dVec ZVec;
  if(w.length()>k_ || w.length()<2) {
    WriteToLogs()(std::string("incorrect length word w in WordWindowFullCount::ZSVector-> ")+w); return ZVec; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowFullCount::ZSVector-> ")+w); return ZVec; 
  }
  if(m<0 || m>(w.length()-2)) {
    WriteIntToLogs()(std::string("invalid order m in WordWindowFullCount::ZSVector-> "),m); return ZVec;
  }
  // get a list of the m+1- and m-words that will be needed, so the substrings only need to be
  // taken once.  Include also the necessary self-overlap words 
  strVec subMp1Words,subMWords;
  for(int p=0;p<(k_-m+1);++p) {
    if(p<(k_-m)) {
      std::string wmp1= w.substr(p,m+1);
      if(std::find(subMp1Words.begin(),subMp1Words.end(),wmp1)==subMp1Words.end()) subMp1Words.push_back(wmp1);
    }
    if(m) {
      std::string wm= w.substr(p,m);
      if(std::find(subMWords.begin(),subMWords.end(),wm)==subMWords.end()) subMWords.push_back(wm);
    }
  }
  // now the self-overlap sub-words
  for(int d=1;d<k_-m;++d)
    if(w.substr(0,k_-d)==w.substr(d,k_-d)) {
      std::string extWord= w.substr(0,d)+w;
      for(int p=(k_-m);p<d;++p) 
	if(std::find(subMp1Words.begin(),subMp1Words.end(),extWord.substr(p,m+1))==subMp1Words.end()) 
	  subMp1Words.push_back(extWord.substr(p,m+1));
      if(m)
	for(int p=(k_-m+1);p<d;++p)
	  if(std::find(subMWords.begin(),subMWords.end(),extWord.substr(p,m))==subMWords.end())
	    subMWords.push_back(extWord.substr(p,m));
    }
  // loop over all windows for this word;  the only tricky thing is expanding the windows rightward for the
  //  m- and m+1-subwords. since the window is being defined for the word w.  For the Z calculation to make sense,
  //  we need to include all valid subwords that would be contained within the last word starting in each window
  // For instance, if we are counting hexamers, then for m= k-2=4, we need an expansion of the window by one for
  //  the hexamers, and by two for the tetramers.
  int nWin= bounds_.size()-1; ZVec.reserve(nWin);
  for(int i=0;i<nWin;++i) {
    int leftBound= bounds_[i],rightBound= bounds_[i+1];
    int mRightBound= rightBound+k_-m, mp1RightBound= rightBound+k_-m-1;
    strIntMap mWords,mp1Words;
    for(int j=0;j<subMp1Words.size();++j) {
      mp1Words[subMp1Words[j]]= wpc_.Count(subMp1Words[j],leftBound-firstP_,mp1RightBound-firstP_);
      //      dumpFile << " subword-m+1 " << subMp1Words[j] << ", count[" << i << "]= " << p1w->second.pCounts[i] << '\n'; 
    }
    if(m)
      for(int j=0;j<subMWords.size();++j) {
	mWords[subMWords[j]]= wpc_.Count(subMWords[j],leftBound-firstP_,mRightBound-firstP_);
	//	dumpFile << " subword-m " << subMWords[j] << ", count[" << i << "]= " << pw->second.pCounts[i] << '\n';
      }
    else {
      mp1Words[val]= wpc_.Count(val,leftBound-firstP_,mp1RightBound-firstP_,m+1);
      //      dumpFile << " valid 1-mers [" << i << "]= " << pw->second.pCounts[i] << '\n';
    } 
    EVRec er(w,m,mp1Words,mWords);
    double wCount=static_cast<double>(wpc_.Count(w,leftBound-firstP_,rightBound-firstP_)); 
    double zs= (wCount-er.E)*(wCount-er.E);
    if(er.E==0.0 && wCount==0.0) zs= 0.0;
    else if(er.E==0) zs= 1.0e30;
    else          zs/= er.E;
    ZVec.push_back(zs);
  }
  //  dumpFile.close();
  if(smooth_) {
    if(sig_ != lastSigSmooth) { 
      GenerateSmoothVectors(SmoothVectors,sig_);
      lastSigSmooth= sig_;
    }
    SmoothVec(SmoothVectors,ZVec);
  }
  return dVec(ZVec);
}
//************************************************************************************************
///
double WordWindowFullCount::ChiSqr(std::string& w,int m, dVec& Z) const {
  //  std::ofstream dumpFile("Zdump.txt",std::ios::app);
  if(w.length()!=k_) {
    WriteToLogs()(std::string("incorrect length word w in WordWindowFullCount::ChiSqr-> ")+w); return -1.0; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowCount::ChiSqr-> ")+w); return -1.0; 
  }
  if(m<0 || m>(w.length()-2)) {
    WriteIntToLogs()(std::string("invalid order m in WordWindowFullCount::ChiSqr-> "),m); return -1.0;
  }
  if(!wpc_.Count(w,m)) return 0.0;
  int n= bounds_.size()-1;
  Z.clear();
  Z= ZVector(w,m);
  double chisqr=0.0;
  for(int i=0;i<n;++i)
    chisqr += Z[i]*Z[i];
  //  dumpFile << w << ", X-sqr= " << chisqr << "\n\n";
  //  dumpFile.close();
  return chisqr;
}
//************************************************************************************************
///
double WordWindowFullCount::ZSSqr(std::string& w,int m, dVec& Z) const {
  //  std::ofstream dumpFile("Zdump.txt",std::ios::app);
  if(w.length()!=k_) {
    WriteToLogs()(std::string("incorrect length word w in WordWindowFullCount::ZSSqr-> ")+w); return -1.0; 
  }
  if(!ValidDNASequence()(w) && !ValidRNASequence()(w)) {
    WriteToLogs()(std::string("invalid word w in WordWindowCount::ZSSqr-> ")+w); return -1.0; 
  }
  if(m<0 || m>(w.length()-2)) {
    WriteIntToLogs()(std::string("invalid order m in WordWindowFullCount::ZSSqr-> "),m); return -1.0;
  }
  if(!wpc_.Count(w,m)) return 0.0;
  int n= bounds_.size()-1;
  Z.clear();
  Z= ZSVector(w,m);
  double zssqr=0.0;
  for(int i=0;i<n;++i)
    zssqr += Z[i];
  //  dumpFile << w << ", X-sqr= " << chisqr << "\n\n";
  //  dumpFile.close();
  return zssqr;
}
//************************************************************************************************
/// 
bool WordWindowFullCount::Smooth(double sig) {
  if(sig<0.0 || sig>50.0) return false;
  sig_=sig; smooth_=true;
  return true;
}
//************************************************************************************************
