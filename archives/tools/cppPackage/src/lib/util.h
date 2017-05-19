//////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __util_h
#define __util_h
//////////////////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
#include <functional>
#include <utility>
#include <limits>
#include <gsl/gsl_rng.h>
#include <math.h>
using std::string; using std::vector; using std::binary_function; using std::unary_function;
typedef vector<string> strVec;
typedef vector<double> dVec;
typedef vector<int> iVec;
//static const double NEG_INF= -1.0* pow(10.0,std::numeric_limits<double>().max_exponent10);
static const double NEG_INF= -1.0e6;//* pow(10.0,std::numeric_limits<double>().max_exponent10);
static const double MIN_LOG= log(std::numeric_limits<double>().min_exponent10);
static const double LOG2_NEGATIVE_INFINITY= -1.0e6;
//////////////////////////////////////////////////////////////////////////////////////////////
class LoadList:public unary_function<string,strVec> {
 public:
  strVec operator()(const string& fn) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class split:public binary_function<string,char,strVec> {
 public:
  strVec operator()(const string& s,char splitChar='\0') const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class join:public binary_function<strVec,char,string> {
 public:
  string operator()(const strVec& sv,char joinChar='\t') const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class rand3:public unary_function<long*,double> {
 public:
  double operator()(long *seed) const;
}; 
//////////////////////////////////////////////////////////////////////////////////////////////
class gasdev:public unary_function<long*,double> {
 public:
  double operator()(long *seed) const;
}; 
//////////////////////////////////////////////////////////////////////////////////////////////
class gammln:public unary_function<double,double> {
 public:
  double operator()(double xx) const;
}; 
//////////////////////////////////////////////////////////////////////////////////////////////
class erfcc:public unary_function<double,double> {
 public:
  double operator()(double xx) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
//class round:public unary_function<int,double> {
// public:
//  int operator()(double xx) const;
//};
//////////////////////////////////////////////////////////////////////////////////////////////
class SafeLog2:public unary_function<double,double> {
 public:
  double operator()(double x) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class SafePow2: public unary_function<double,double> {
 public:
  double operator()(double x) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class SafeLog:public unary_function<double,double> {
 public:
  double operator()(double x) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class SafeExp:public unary_function<double,double> {
 public:
  double operator()(double x) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class betai {
 protected:
  double res_;
  betai() {}
 public:
  betai(const betai& o): res_(o.res_) {}
  betai(double a,double b,double x);
  double operator()() const { return res_; }
};
//////////////////////////////////////////////////////////////////////////////////////////////
class logFac:public unary_function<int,double> {
 public:
  double operator()(int x) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class DrawFromProbVector:public binary_function<dVec,long,int> { 
 public:
  int operator()(const dVec& v,long& seed) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class DrawAgainstSingleProb:public binary_function<double,long,bool> { 
 public:
  bool operator()(double v,long& seed) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class InvertProbVector:public binary_function<dVec,bool,dVec> { 
 public:
  dVec operator()(const dVec& v,bool normalize) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class HeatProbVector:public binary_function<dVec,double,dVec> { 
 public:
  dVec operator()(const dVec& v,double T) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class ScaledInvChiSqVector:public unary_function<dVec,dVec> {
 public:
  dVec operator()(const dVec& v) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
typedef std::pair<int,int> intRange;
inline intRange make_ordered_IntRange(int i0,int i1) {
  if(i0<=i1) return std::make_pair(i0,i1);
  return std::make_pair(i1,i0);
}
inline intRange make_IntRange(int i0,int i1) {
  return std::make_pair(i0,i1); }
inline bool overlap (const intRange& r1,const intRange& r2) {
  return !(r2.second<=r1.first || r1.second<=r2.first); }
inline bool contains (const intRange& r1,const intRange& r2) {
  return (r1.first<=r2.first && r1.second>=r2.second); }
inline bool contains (const intRange& r,int v) { 
  return (r.first<=v && v<r.second); }
inline bool operator< (int v,const intRange& r) {
  return v<r.first; }
inline bool operator> (int v,const intRange& r) {
  return v>=r.second; }
inline int length (const intRange& r) {
  return r.second-r.first;
}
bool merge(intRange& r1,const intRange& r2);
//////////////////////////////////////////////////////////////////////////////////////////////
typedef vector<intRange> irVec;
//////////////////////////////////////////////////////////////////////////////////////////////
class GetAccession:public unary_function<std::string,std::string> {
 public:
  std::string operator()(std::string line) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class IntToString:public std::unary_function<int,std::string> {
 public:
  std::string operator()(int i) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class DoubleToString:public std::binary_function<double,int,std::string> {
 public:
  std::string operator()(double d,int dd= -1) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class DoubleToFixedString:public std::binary_function<double,int,std::string> {
 public:
  std::string operator()(double d,int dd= -1) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class MakeFilenameSafe:public std::unary_function<std::string,std::string> {
 public:
  std::string operator()(const std::string& fn) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class CumulativePoisson {
 private:
  double cp_;
 public:
  CumulativePoisson(int x,int n,double p);
  double operator()() const { return cp_; }
  double P() const { return cp_; }
};
//////////////////////////////////////////////////////////////////////////////////////////////
class CumulativeBinary {
 private:
  double cp_;
 public:
  CumulativeBinary(int x,int n,double p);
  double operator()() const { return cp_; }
  double P() const { return cp_; }
};
//////////////////////////////////////////////////////////////////////////////////////////////
class GSLrng {
 private:
  const gsl_rng_type *T_;
  gsl_rng *r_;
  unsigned seed_;
 public:
  GSLrng();
  GSLrng(unsigned seed,const gsl_rng_type *T= 0);
  ~GSLrng(){ gsl_rng_free(r_); }
  //
  unsigned get()                      { return gsl_rng_get(r_); }
  double uniform()                    { return gsl_rng_uniform(r_); }
  double uniform_pos()                { return gsl_rng_uniform_pos(r_); }
  unsigned uniform_int(unsigned range){ return gsl_rng_uniform_int(r_,range); }
  //
  std::string name() const { return std::string(gsl_rng_name(r_)); }
  unsigned max()     const { return gsl_rng_max(r_); }
  unsigned min()     const { return gsl_rng_min(r_); }
  unsigned seed()    const { return seed_; }
  void reseed(unsigned seed) { seed_=seed; gsl_rng_set(r_,seed_); }
};
//////////////////////////////////////////////////////////////////////////////////////////////
class AddLog2:public binary_function <double,double,double> {
 public:
  double operator()(double log2val1,double log2val2) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class AddLog:public std::binary_function<double,double,double> {
 public:
  double operator()(double logp,double logq) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
#endif
//////////////////////////////////////////////////////////////////////////////////////////////
