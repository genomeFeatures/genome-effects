////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __align_h
#define __align_h
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <string>
using std::vector;
using std::string;
typedef vector<int> iVec;
typedef vector<iVec> iMat;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
class TJalign {
 protected:
  int gInit_,gExt_,match_,miss_,score_;
  string s1_,s2_,a1_,a2_; iMat S_;
  virtual bool BuildMatrixes_()=0;
  TJalign():gInit_(),gExt_(),match_(),miss_(),score_(),s1_(),s2_(),a1_(),a2_(),S_() {}
 public:
  TJalign(const TJalign& o):gInit_(o.gInit_),gExt_(o.gExt_),match_(o.match_),miss_(o.miss_),\
    score_(o.score_),s1_(o.s1_),s2_(o.s2_),a1_(o.a1_),a2_(o.a2_),S_(o.S_) {}
  TJalign(const string& s1,const string& s2,int match,int mis,int gInit,int gExt);
  TJalign(const string& s1,const string& s2,const iMat& S,int gInit,int gExt);
  //
  int Alignment(string& a1,string& a2) const { a1=a1_; a2=a2_; return score_; }
  int operator()() const { return score_; }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
class SmithWaterman:public TJalign {
 protected:
  virtual bool BuildMatrixes_();
  SmithWaterman():TJalign() {}
  int start1_,stop1_,start2_,stop2_;
 public:
  SmithWaterman(const string& s1,const string& s2,int match,int mis,int gInit,int gExt);
  SmithWaterman(const string& s1,const string& s2,const iMat& S,int gInit,int gExt);
  int Start1() { return start1_; }
  int Stop1() { return stop1_; }
  int Start2() { return start2_; }
  int Stop2() { return stop2_; }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
class NeedlemanWunsch:public TJalign {
 protected:
  bool pEnd1_,pEnd2_;
  virtual bool BuildMatrixes_();
  NeedlemanWunsch():TJalign() {}
 public:
  NeedlemanWunsch(const NeedlemanWunsch& o):TJalign(o),pEnd1_(o.pEnd1_),pEnd2_(o.pEnd2_) {}
  NeedlemanWunsch(const string& s1,const string& s2,bool e1,bool e2,int match,int mis,int gInit,int gExt);
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
#endif
