///////////////////////////////////////////////////////////////////////////////////
#ifndef __vectorF_h
#define __vectorF_h
#include "util.h"
#include <iostream>
#include <vector>
#include <string>
#include <functional>
using std::vector; using std::binary_function; using std::unary_function; using std::string;
typedef vector<double> dVec;
typedef vector<int> iVec;
///////////////////////////////////////////////////////////////////////////////////
// hadamard product: results in a vector for which each element is the product
//                   of the corresponding elements of the two input vectors
class HadProduct:public binary_function<dVec,dVec,dVec>  {
 public:
  dVec operator()(const dVec& v1,const dVec& v2) const;
};
///////////////////////////////////////////////////////////////////////////////////
class InnerProduct:public binary_function<dVec,dVec,double>  {
 public:
  double operator()(const dVec& v1,const dVec& v2) const;
};
///////////////////////////////////////////////////////////////////////////////////
class VectorSum:public binary_function<dVec,bool,double>  {
 public:
  double operator()(dVec& v1,bool normalize=false) const;
};
///////////////////////////////////////////////////////////////////////////////////
class VectorSumSq:public binary_function<dVec,bool,double>  {
 public:
  double operator()(dVec& v1,bool normalize=false) const;
};
///////////////////////////////////////////////////////////////////////////////////
class VectorChiSq:public binary_function<dVec,bool,double>  {
 public:
  double operator()(dVec& v1,bool normalize=false) const;
};
///////////////////////////////////////////////////////////////////////////////////
class InnerProductN:public binary_function<dVec,dVec,double>  {
 public:
  double operator()(const dVec& v1,const dVec& v2) const;
};
////////////////////////////////////////////////////////////////////////////////////
class EuclideanDistance:public binary_function<dVec,dVec,double>  {
 public:
  double operator()(const dVec& v1,const dVec& v2) const;
};
////////////////////////////////////////////////////////////////////////////////////
class Pearson:public binary_function<dVec,dVec,double>  {
 public:
  double operator()(const dVec& v1,const dVec& v2) const;
};
////////////////////////////////////////////////////////////////////////////////////
class Spearman:public binary_function<dVec,dVec,double>  {
 public:
  double operator()(const dVec& v1,const dVec& v2) const;
};
////////////////////////////////////////////////////////////////////////////////////
// returns a rank vector corresponding to the input vector
class RankVector:public unary_function<dVec,dVec> {
 public:
  dVec operator()(const dVec& v) const;
};
////////////////////////////////////////////////////////////////////////////////////
class SortedIndexVector:public unary_function<dVec,iVec> {
 public:
  iVec operator()(const dVec& v) const;
};
////////////////////////////////////////////////////////////////////////////////////
// generates an N-long vector of ones
class Ones:public unary_function<int,dVec> {
 public:
  dVec operator()(int N) const;
};
////////////////////////////////////////////////////////////////////////////////////
// generates a zero vector of length N
class Zeroes:public unary_function<int,dVec> {
 public:
  dVec operator()(int N) const;
};
////////////////////////////////////////////////////////////////////////////////////
// generates a random vector of length N
class RandomVector:public unary_function<int,dVec> {
 public:
  dVec operator()(int N) const;
};
////////////////////////////////////////////////////////////////////////////////////
// writes the vector to a stream, tab-delimited, and includes an end of line character
std::ostream& operator << (std::ostream& os, const dVec& v); 
////////////////////////////////////////////////////////////////////////////////////
class TextBarGraph:public binary_function<dVec,int,string> {
 public:
  string operator()(const dVec& d,int width) const;
};
////////////////////////////////////////////////////////////////////////////////////
class TextBarGraphComp:public binary_function<dVec,int,string> {
 public:
  string operator()(const dVec& d,int width) const;
};
////////////////////////////////////////////////////////////////////////////////////
class HorizontalTextBarGraph:public binary_function<dVec,int,string> {
 public:
  string operator()(const dVec& d,int width) const;
};
////////////////////////////////////////////////////////////////////////////////////
class Histogram {
 protected:
  static const int MaxBins=1002;
  iVec bins_;
  dVec bounds_;
  double minBound_,binSize_;
  int maxCount_;
  Histogram();
 public:
  Histogram(const Histogram& h):bins_(h.bins_),bounds_(h.bounds_),minBound_(h.minBound_),binSize_(h.binSize_),maxCount_(h.maxCount_) {}
  Histogram(int nBins,double minBound,double binSize);
  Histogram(const dVec& bounds);
  //
  double MinBound() const { return minBound_; }
  double BinSize() const { return binSize_; }
  int MaxCount() const { return maxCount_; }
  int NBins() const { return bins_.size(); }
  //
  int Add(double v);
  int Add(const dVec& v);
  int Collapse(int max);
  std::string TableList(bool transpose,bool header,int p=-1,int w=-1) const;
  std::string HorizontalBarGraph(int width) const;
  std::string VerticalBarGraph(int width) const;
  std::string HorizontalBarGraph(int M,int N) const;
  std::string VerticalBarGraph(int M,int N) const;
};
////////////////////////////////////////////////////////////////////////////////////
class KS_pairtest {
 private: 
  dVec data1_,data2_;
  double d_,p_;
 public:
  KS_pairtest(const dVec& d1, const dVec& d2);
  int N() const { return data1_.size(); }
  double D() const { return d_; }
  double P() const { return p_; }
};
////////////////////////////////////////////////////////////////////////////////////
class ChiSq_pairtest {
 private: 
  dVec data1_,data2_; bool cons_;
  double x2_,p_; int df_;
 public:
  ChiSq_pairtest(const dVec& d1, const dVec& d2, bool c=false);
  int N() const { return data1_.size(); }
  double X2() const { return x2_; }
  double P() const { return p_; }
  int DF() const { return df_; }
};
////////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////////
