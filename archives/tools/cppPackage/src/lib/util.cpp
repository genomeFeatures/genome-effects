//////////////////////////////////////////////////////////////////////////////////////////////
#define _GNU_SOURCE 1
#include <math.h>
#ifdef INTEL
#include <mathimf.h>
#endif
#include <time.h>
#include <ctype.h>
#include "../lib/util.h"
#include "../lib/logFile.h"
#include <iostream>
#include <strstream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_randist.h>
const int MaxFilenameLength=128;
//////////////////////////////////////////////////////////////////////////////////////////////
strVec LoadList::operator()(const string& fn) const {
  std::ifstream inFile(fn.c_str(),std::ios::in);
  if(inFile.fail()) {
    WriteToLogs()(string("Error opening list file: ")+fn);
    return strVec();
  }
  strVec ret; ret.reserve(1000);
  char *lineIn= new char[5001]; 
  inFile.getline(lineIn,5000);
  while(!inFile.eof()) {
    ret.push_back(string(lineIn));
    inFile.getline(lineIn,5000);
  }
  inFile.close();
  delete[] lineIn;
  return ret;
}
//////////////////////////////////////////////////////////////////////////////////////////////
strVec split::operator()(const string& s,char splitChar) const {
  strVec fields; string ss;
  if(splitChar==' ') {
    int i=0; while(s[i]==' ') ++i;
    while(i<s.length()) {
      if(s[i]!=' ') ss += s[i];
      else if(!i) ss+= s[i];
      else if(s[i-1]!=' ') ss+=s[i];
      ++i;
    }
  } else if(splitChar=='\0') {
    int i=0; while(isspace(s[i])) ++i;
    while(i<s.length()) {
      if(!isspace(s[i])) ss += s[i];
      else if(!i) ss+= s[i];
      else if(!isspace(s[i-1])) ss+=s[i];
      ++i;
    }
  } else
    ss=s;
  int sCount=0;
  for(int i=0;i<ss.length();++i) if(ss[i]==splitChar)++sCount;
  fields.reserve(sCount+1);
  //  if(splitChar=='+') WriteToLog()(std::string("splitting: ")+std::string(ss));
  bool inWord=false; string word;
  for(int i=0;i<ss.length();++i) {
    bool inWordi= !(ss[i]==splitChar || (splitChar=='\0' && isspace(ss[i])));
    if(inWord && !inWordi) {
      fields.push_back(word);
      //      if(splitChar=='+') WriteToLog()(std::string("word: ")+word);
      word.clear();
    } else if (!inWord && !inWordi && splitChar!='\0') {
      fields.push_back(word);
    } else if (inWordi) 
      word += ss[i];
    inWord=inWordi;
  }
  if(inWord) fields.push_back(word);
  else if(splitChar!='\0') fields.push_back(word);
  return fields;
}
//////////////////////////////////////////////////////////////////////////////////////////////
string join::operator()(const strVec& sv,char joinChar) const {
  if(sv.empty()) return string();
  string s;
  s+= sv[0];
  for(unsigned i=1;i<sv.size();i++) { s+= joinChar; s+= sv[i]; }
  return string(s);
}
//////////////////////////////////////////////////////////////////////////////////////////////
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
double rand3::operator()(long *idum) const {
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  //
  if((*idum<0)||(iff==0)) {
    iff=1;
    mj=labs(MSEED-labs(*idum));
    mj %= MBIG;
    ma[55]= mj;
    mk=1;
    for(i=1;i<=54;i++) {
      ii=(21*i)%55;
      ma[ii]=mk;
      mk=mj-mk;
      if(mk<MZ) mk += MBIG;
      mj=ma[ii];
    }
    for(k=1;k<=4;k++)
      for(i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if(ma[i]<MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if(++inext == 56) inext=1;
  if(++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if(mj<MZ) mj += MBIG;
  ma[inext]=mj;
  return (mj*FAC);
}
//////////////////////////////////////////////////////////////////////////////////////////////
double gasdev::operator()(long *idum) const {
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  if(*idum < 0) iset=0;
  if(iset == 0) {
    do {
      v1= 2.0*rand3()(idum)-1.0;
      v2= 2.0*rand3()(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0);
    fac= sqrt(-2.0*log(rsq)/rsq);
    gset= v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////
double SafeLog2::operator()(double x) const { 
  if(x > 0.0) return log2(x); 
  if(x < 0.0) WriteFloatToErr()(std::string("negative value in SafeLog2: "),x);
  return LOG2_NEGATIVE_INFINITY; 
}
//////////////////////////////////////////////////////////////////////////////////////////////
double SafePow2::operator()(double x) const {
  if(x <= LOG2_NEGATIVE_INFINITY) return 0.0;
  return pow(2.0,x);
}
//////////////////////////////////////////////////////////////////////////////////////////////
double SafeLog::operator()(double x) const {
  if(x > 0.0) return log(x);
  if(x < 0.0) WriteFloatToErr()(std::string("negative value in SafeLog: "),x);
  return NEG_INF;
}
//////////////////////////////////////////////////////////////////////////////////////////////
double SafeExp::operator()(double x) const {
  if(x <= NEG_INF) return 0.0;
  return exp(x);
}
//////////////////////////////////////////////////////////////////////////////////////////////
//int round::operator()(double xx) const { 
//  if((xx-floor(xx))<0.5) return static_cast<int>(floor(xx)); 
//  return static_cast<int>(floor(xx)+1.0);
//}
//////////////////////////////////////////////////////////////////////////////////////////////
double gammln::operator()(double xx) const {
  static double cof[6]= {76.18009172947146,-86.50532032941677,24.01409824083091,
			 -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  double x,y;
  y=x=xx;
  double tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  double ser=1.000000000190015;
  for(int j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}
//////////////////////////////////////////////////////////////////////////////////////////////
double erfcc::operator()(double x) const {
  double t,z,ans;
  
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+
  				      t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
  //return gsl_sf_erfc(x);
}
//////////////////////////////////////////////////////////////////////////////////////////////
#define MAXIT 10000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
double betacf(double a,double b,double x) {
  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;
  //
  WriteFloatToLog()(string("betacf in, x="),x);
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if(fabs(d)<FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for(m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if(fabs(d)<FPMIN)d=FPMIN;
    c=1.0+aa/c;
    if(fabs(c)<FPMIN)c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if(fabs(d)<FPMIN)d=FPMIN;
    c=1.0+aa/c;
    if(fabs(c)<FPMIN)c=FPMIN;
    d=1.0/d;
    del=d*c;
    h*=del;
    if(fabs(del-1.0)<EPS)break;
  }
  if(m>MAXIT) {
    WriteToLogs()(string("a or b too big, or MAXIT too small in betacf"));
    WriteFloatToErr()(string("a:"),a);
    WriteFloatToErr()(string("b:"),b);
    WriteFloatToErr()(string("x:"),x);
    WriteFloatToErr()(string("del:"),del);
    WriteIntToErr()(string("m:"),m);
    exit(-1);
  }
  return h;
}
//////////////////////////////////////////////////////////////////////////////////////////////
betai::betai(double a,double b,double x) {
  if((x<0.0)||(x>1.0)) {
    WriteToLogs()(string("bad x value in betai"));
    exit(-1);
  }
  double bt;
  if((x==0.0)||(x==1.0)) 
    bt=0.0;
  else 
    bt=exp(gammln()(a+b)-gammln()(a)-gammln()(b)+a*log(x)+b*log(1.0-x));
  if(x<(a+1.0)/(a+b+2.0))
    res_= bt*betacf(a,b,x)/a;
  else
    res_= 1.0-bt*betacf(b,a,1.0-x)/b; 
}
static bool firstLF(true);
typedef vector<double> dVec;
static dVec lf;
//////////////////////////////////////////////////////////////////////////////////////////////
double logFac::operator()(int x) const {
  if(x>2000) return gammln()(static_cast<double>(x));
  if(firstLF) {
    lf.reserve(2001); double curr=0.0; lf.push_back(0.0); lf.push_back(0.0);
    for(unsigned i=2;i<=2000;i++) {
      curr+= log(static_cast<double>(i));
      lf.push_back(curr);
    }
    firstLF=false;
  }
  return lf[x];
}
//////////////////////////////////////////////////////////////////////////////////////////////
int DrawFromProbVector::operator()(const dVec& v,long& seed) const {
  if(v.empty()) {
    WriteToErr()(string("zero size vector v in DrawFromProbVector"));
    return -1;
  }
  double sum=0.0; 
  for(unsigned i=0;i<v.size();i++) 
    if(v[i]>=0.0) sum+=v[i];
    else { WriteToErr()(string("negative element in vector v in DrawFromProbVector")); return -1; }
  double r= rand3()(&seed);
  int p=0; double cdf=v[0]/sum;
  while ((r>cdf)&&(p<(v.size()-1))) {
    cdf += v[++p]/sum;
    if(cdf>(1.0+1e-3)) {
      WriteFloatToErr()(string("cdf exceeds 1.0 in drawRandPos: "),cdf);
      WriteIntToErr()(string("v.size()= "),v.size());
      WriteIntToErr()(string("p = "),p);
      WriteFloatToErr()(string("sum= "),sum);
      char line[1001]; std::strstream buff(line,1000);
      buff << "values: "; for(unsigned i=0;i<v.size();i++) buff << " " << v[i];
      buff << "\nscaled: "; for(unsigned i=0;i<v.size();i++) buff << " " << (v[i]/sum);
      buff << '\0';
      WriteToErr()(string(line));
      return p;
    }
  }
  return p;
} 
//////////////////////////////////////////////////////////////////////////////////////////////
bool DrawAgainstSingleProb::operator()(double v,long& seed) const {
  if((v<0.0)||(v>1.0)) {
    WriteFloatToLog()(string("v out of range in DrawAgainstSingleProb"),v);
    WriteFloatToErr()(string("v out of range in DrawAgainstSingleProb"),v);
    return false;
  }
  double r= rand3()(&seed);
  return r<v;
} 
//////////////////////////////////////////////////////////////////////////////////////////////
dVec InvertProbVector::operator()(const dVec& v,bool normalize) const {
  if(v.empty()) { WriteToErr()(string("zero size vector v in InvertProbVector")); return dVec();  }
  // check for zero or negative values
  bool zeroFound=false;
  for(unsigned i=0;i<v.size();++i)
    if(v[i]==0.0) zeroFound=true;
    else if (v[i]<0.0) { WriteToErr()(string("negative element in vector v in InvertProbVector")); return dVec(); }
  // invert the array
  dVec inv; double sum=0.0;
  for(unsigned i=0;i<v.size();++i) {
    if(!zeroFound)     inv.push_back(1.0/v[i]);
    else if(v[i]!=0.0) inv.push_back(0.0);
    else               inv.push_back(1.0);
    if(normalize) sum+= inv[i];
  }
  if(normalize) for(unsigned i=0;i<inv.size();i++) inv[i] /= sum;
  return dVec(inv);
}
//////////////////////////////////////////////////////////////////////////////////////////////
dVec HeatProbVector::operator()(const dVec& v,double T) const {
  if(v.empty()) { WriteToErr()(string("zero size vector v in InvertProbVector")); return dVec(); }
  // 
  dVec hProb(v.size(),double(0.0)); 
  for(unsigned i=0;i<v.size();++i) {
    if(v[i]<0.0) { WriteToLogs()(string("negative element in vector v in InvertProbVector")); return dVec(); }
    if(v[i]>0.0) hProb[i]= exp(log(v[i])/T);
  }
  return hProb;
}
//////////////////////////////////////////////////////////////////////////////////////////////
dVec ScaledInvChiSqVector::operator()(const dVec& vec) const {
  if(vec.size()!=5) {
    WriteToLogs()(string("ScaledInvChiSqVector: vec.size()!=5")); 
    return dVec();
  }
  double X=vec[0],X1=vec[1],dX=vec[2],v=vec[3],s2=vec[4];
  dVec ret; ret.reserve(static_cast<int>((X1-X)/dX));
  double s=sqrt(s2),vo2=v/2.0;
  while(X<=X1) {
    if(X<=0.0) ret.push_back(0.0);
    else       ret.push_back(exp(vo2*log(vo2)-gammln()(vo2)+v*log(s)-(vo2+1)*log(X)-vo2*s2/X));
    X += dX;
  }
  return ret;
}
//////////////////////////////////////////////////////////////////////////////////////////////
bool merge(intRange& r1,const intRange& r2) {
  if(!overlap(r1,r2)) return false;
  if(r2.first<r1.first) r1.first=r2.first;
  if(r2.second>r1.second) r1.second=r2.second;
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////
std::string GetAccession::operator() (std::string line) const {
  if(line.empty()) return std::string();
  int start= (line[0]=='>')?1:0;
  int p= line.find(' ');
  if(p>=0) line=line.substr(start,p-start);
  strVec f=split()(line,'|');
  return *(f.rbegin());
}
//////////////////////////////////////////////////////////////////////////////////////////////
std::string IntToString::operator() (int i) const {
  char line[25]; std::strstream buff(line,24,std::ios::out);
  buff << i << '\0';
  return std::string(line);
}
//////////////////////////////////////////////////////////////////////////////////////////////
std::string DoubleToString::operator() (double d,int dd) const {
  char line[25]; std::strstream buff(line,24,std::ios::out);
  if(dd<=0)  buff << d << '\0';
  else       buff << std::setprecision(dd) << d << '\0';
  return std::string(line);
}
//////////////////////////////////////////////////////////////////////////////////////////////
std::string DoubleToFixedString::operator() (double d,int dd) const {
  char line[25]; std::strstream buff(line,24,std::ios::out);
  if(dd<0 || dd>10)     buff << d << '\0';
  else if(!dd) buff << round(d) << '\0';
  else {
    double factor= pow(10,dd);
    buff << (round(d*factor)/factor) << '\0';
  }
  return std::string(line);
}
//////////////////////////////////////////////////////////////////////////////////////////////
std::string MakeFilenameSafe::operator()(const std::string& fn) const {
  std::string safeFN; bool IllegalChars=false;
  for(int i=0;i<fn.length() && i<MaxFilenameLength;++i)
    if( ((fn[i]>='a')&&(fn[i]<='z')) || ((fn[i]>='A')&&(fn[i]<='Z')) || ((fn[i]>='0')&&(fn[i]<='9')) 
	|| (fn[i]=='.') || (fn[i]==':') || (fn[i]=='-') || (fn[i]=='+') || (fn[i]=='_') || (fn[i]=='/'))
      safeFN += fn[i];
    else
      IllegalChars=true;
  if(IllegalChars) WriteToLogs()(std::string("Proposed file name \"")+fn+std::string(" has illegal characters, reduced to: ")+safeFN);
  return safeFN;
}
//////////////////////////////////////////////////////////////////////////////////////////////
CumulativePoisson::CumulativePoisson(int x,int n,double p) {
  double mu=p* static_cast<double>(n);
  cp_=gsl_ran_poisson_pdf(0.0,mu);
  for(int i=1;i<=x;++i) 
    cp_+=gsl_ran_poisson_pdf(i,mu);
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
CumulativeBinary::CumulativeBinary(int x,int n,double p) {
  cp_=gsl_ran_binomial_pdf(0,p,n);
  for(int i=1;i<=x;++i) 
    cp_+=gsl_ran_binomial_pdf(i,p,n);
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////
bool GSLrngEnv=false;
//////////////////////////////////////////////////////////////////////////////////////////////
GSLrng::GSLrng():seed_(0) {
  if(!GSLrngEnv) {
    gsl_rng_env_setup();
    GSLrngEnv=true;
  }
  T_= gsl_rng_default;
  r_= gsl_rng_alloc(T_);
  gsl_rng_set(r_,seed_);
}
//////////////////////////////////////////////////////////////////////////////////////////////
GSLrng::GSLrng(unsigned seed,const gsl_rng_type *T):seed_(seed) {
  if(!GSLrngEnv) {
    gsl_rng_env_setup();
    GSLrngEnv=true;
  }
  if(T) T_=T;
  else  T_= gsl_rng_default;
  r_=gsl_rng_alloc(T_);
  gsl_rng_set(r_,seed_);
}
//////////////////////////////////////////////////////////////////////////////////////////////
double AddLog2::operator()(double log2val1,double log2val2) const {
  double max, diff;
  if(log2val2<=LOG2_NEGATIVE_INFINITY && log2val1<=LOG2_NEGATIVE_INFINITY) 
    return LOG2_NEGATIVE_INFINITY;
  if(log2val2<=LOG2_NEGATIVE_INFINITY) return log2val1;
  if(log2val1<=LOG2_NEGATIVE_INFINITY) return log2val2;
  if (log2val1 > log2val2) { max = log2val1; diff = log2val2 - log2val1; } 
  else                     { max = log2val2; diff = log2val1 - log2val2; }
  if(diff <= MIN_LOG) return max;
  // Now diff <= 0 so Math.exp(diff) will not overflow
  return max + log2(1 + pow(2,diff));
}
//////////////////////////////////////////////////////////////////////////////////////////////
double AddLog::operator()(double logp,double logq) const {
  if(logp <= NEG_INF && logq <= NEG_INF) return NEG_INF;
  if(logp <= NEG_INF) return logq;
  if(logq <= NEG_INF) return logp;
  double max,diff;
  if(logp>logq) { max=logp; diff= logq-logp; }
  else          { max=logq; diff= logp-logq; }
  if(diff <= MIN_LOG) return max;
  //eventually replace the calls below to log and exp with a table lookup for ~35x speedup
  return max + log(1 + exp(diff));
}
//////////////////////////////////////////////////////////////////////////////////////////////




