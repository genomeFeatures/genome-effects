#include <time.h>
#define _GNU_SOURCE 1
#include <math.h>
#include <iomanip>
#include <strstream>
#include <algorithm>
#include <gsl/gsl_sort_double.h>
#include "../lib/vectorF.h"
#include "../lib/util.h"
#include "../lib/logFile.h"
typedef vector<int> iVec;
///////////////////////////////////////////////////////////////////
double VectorSum::operator()(dVec& v1,bool normalize) const {
  double s=0.0;
  for(unsigned i=0;i<v1.size();++i) s+= v1[i];
  if(normalize) for(unsigned i=0;i<v1.size();++i) v1[i] /= s;
  return double(s);
}
///////////////////////////////////////////////////////////////////
double VectorSumSq::operator()(dVec& v1,bool normalize) const {
  double s=0.0;
  for(unsigned i=0;i<v1.size();++i) s+= v1[i]*v1[i];
  if(normalize) for(unsigned i=0;i<v1.size();++i) v1[i] /= s;
  return double(s);
}
///////////////////////////////////////////////////////////////////
double VectorChiSq::operator()(dVec& v1,bool normalize) const {
  double s=0.0,n=static_cast<double>(v1.size());
  for(unsigned i=0;i<v1.size();++i) s+= v1[i];
  double avg = s/n,chi=0.0;
  for(unsigned i=0;i<v1.size();++i) chi+= (v1[i]-avg)*(v1[i]-avg);
  chi /= avg;
  if(normalize) for(unsigned i=0;i<v1.size();++i) v1[i] /= s;
  return double(chi);
}
///////////////////////////////////////////////////////////////////
dVec HadProduct::operator()(const dVec& x,const dVec& y) const {
  if(x.size()!=y.size()) return dVec(); // only defined for equal size vectors
  dVec p; p.reserve(x.size());
  for(unsigned i=0;i<x.size();++i) p.push_back(x[i]*y[i]);
  return dVec(p);
}
///////////////////////////////////////////////////////////////////
double InnerProduct::operator()(const dVec& x,const dVec& y) const {
  if(x.size()!=y.size()) return -1.0e10; // only defined for equal size vectors
  double xy=0.0;
  for(unsigned i=0;i<x.size();++i) xy += x[i]*y[i];
  return xy;
}
///////////////////////////////////////////////////////////////////
double InnerProductN::operator()(const dVec& x,const dVec& y) const {
  if(x.size()!=y.size()) return -1.0e10; // only defined for equal size vectors
  double ssX=0.0,ssY=0.0,xy=0.0;
  for(unsigned i=0;i<x.size();++i) {
    ssX += x[i]*x[i]; ssY += y[i]*y[i]; xy += x[i]*y[i];
  }
  if((ssX==0.0)||(ssY==0.0)) return 0.0;
  return xy/sqrt(ssX*ssY);
}
///////////////////////////////////////////////////////////////////
double EuclideanDistance::operator()(const dVec& x,const dVec& y) const {
  if(x.size()!=y.size()) return -1.0e10; // only defined for equal size vectors
  double d=0.0;
  for(unsigned i=0;i<x.size();++i)
    d += (x[i]-y[i])*(x[i]-y[i]);
  return sqrt(d);
}
///////////////////////////////////////////////////////////////////
double Pearson::operator()(const dVec& x,const dVec& y) const { 
  if(x.size()!=y.size()) return -1.0e10; // only defined for equal size vectors
  double xBar=0.0,yBar=0.0;
  for(unsigned i=0;i<x.size();++i) {
    xBar += x[i]; yBar += y[i];
  }
  yBar /= (double)x.size(); xBar /= (double)x.size();
  double sum1=0.0,sum2=0.0,sum3=0.0;
  for(unsigned i2=0;i2<x.size();i2++) {
    sum1 += (x[i2]-xBar)*(y[i2]-yBar); 
    sum2 += (x[i2]-xBar)*(x[i2]-xBar);
    sum3 += (y[i2]-yBar)*(y[i2]-yBar);
  }
  double r;
  if((sum2!=0.0)&&(sum3!=0.0))r = sum1/sqrt(sum2*sum3);
  else if(sum1==0.0) r= 0.0;
  else WriteToLogs()(string("sum problem in Pearson"));
  return r;
}
///////////////////////////////////////////////////////////////////
class rec { 
public:
  int rank; double val;
  rec(int r,double v):rank(r),val(v) {}
  rec(const rec& o):rank(o.rank),val(o.val) {}
  rec():rank(),val() {}
  rec& operator=(const rec& o) { rank=o.rank; val=o.val; return *this; }
  bool operator<(const rec& o) const { return val<o.val; }
};
typedef std::vector<rec> recVec;
//
dVec RankVector::operator()(const dVec& v) const {
  recVec rv; rv.reserve(v.size());
  for(unsigned i=0;i<v.size();++i) rv.push_back(rec(i,v[i]));
  std::sort(rv.begin(),rv.end());
  // next copy to the double vector, accounting for ties
  dVec R(v.size(),0.0);
  unsigned i=0;
  while(i<(v.size()-1)) 
    if(v[rv[i].rank]!=v[rv[i+1].rank]) {
      R[rv[i].rank]=i+1;
      ++i;
    } else {
      unsigned it;
      for(it=i+1;(it<v.size())&&(v[rv[it].rank]==v[rv[i].rank]);it++);
      double rank=0.5*(i+it+1);
      for(unsigned ii=i;ii<=(it-1);++ii)R[rv[ii].rank]=rank;
      i=it;
    }
  if(i==(v.size()-1)) R[rv[i].rank]=rv[i].rank+1;
  return dVec(R);
}
///////////////////////////////////////////////////////////////////
iVec SortedIndexVector::operator()(const dVec& v) const {
  iVec r; for(unsigned i=0;i<v.size();++i) r.push_back(i);
  // first do a bubble sort on index
  for(unsigned i=0;i<(v.size()-1);++i)
    for(unsigned j=i+1;j<v.size();++j) 
      if(v[r[i]]<v[r[j]]) { unsigned t=r[i]; r[i]=r[j]; r[j]=t; }
  return iVec(r);
}
///////////////////////////////////////////////////////////////////
double Spearman::operator()(const dVec& x,const dVec& y) const {
  if(x.size()!=y.size()) return -1.0e10; // only defined for equal size vectors
  dVec R=RankVector()(x), S=RankVector()(y);
  return Pearson()(R,S);
}
///////////////////////////////////////////////////////////////////
unsigned rowSize=10;
std::ostream& operator << (std::ostream& os,const dVec& v) {
  for(unsigned i=1;i<=v.size();++i) 
    if((i%rowSize)&&(i%rowSize == 1)) os << v[i-1];
    else if (i%rowSize)               os << "\t" << v[i-1];
    else                              os << "\t" << v[i-1] << std::endl;
  if(v.size()%rowSize) os << std::endl;
  return os;
}
///////////////////////////////////////////////////////////////////
dVec Ones::operator()(int N) const {
  if(N<=0) {
    WriteIntToErr()(string("Invalid N in Ones"),N);
    return dVec();
  }
  dVec o;
  for(unsigned i=0;i<N;++i) o.push_back(1.0);
  return dVec(o);
}
///////////////////////////////////////////////////////////////////
dVec Zeroes::operator()(int N) const {
  if(N<=0) {
    WriteIntToErr()(string("Invalid N in Zeroes"),N);
    return dVec();
  }
  dVec z;
  for(unsigned i=0;i<N;++i) z.push_back(0.0);
  return dVec(z);
}
///////////////////////////////////////////////////////////////////
dVec RandomVector::operator()(int N) const {
  if(N<=0) {
    WriteIntToErr()(string("Invalid N in RandomVector"),N);
    return dVec();
  }
  long seed; time_t t0; time(&t0); seed= -1*t0;
  dVec r;
  for(unsigned i=0;i<N;++i) r.push_back(rand3()(&seed));
  return dVec(r);
}
///////////////////////////////////////////////////////////////////
string TextBarGraph::operator()(const dVec& d,int width) const {
  double maxVal= 0.0;
  for(unsigned i=0;i<d.size();++i) if(d[i]>maxVal)maxVal=d[i];
  double vSteps=20.0;
  string retString; retString= string();
  for(unsigned i=0;i<d.size();i+=width) {
    for(unsigned row=0;row<vSteps;row++) {
      char line2[width+3]; std::strstream buff2(line2,width+2,std::ios::out);
      for(unsigned j=0;(j<width)&&((i+j)<d.size());++j) 
	if((vSteps*d[i+j]/maxVal)>(vSteps-row-1)) buff2 << '*';
	else                                       buff2 << ' ';
      buff2 << std::endl << '\0'; retString += string(line2);
    }
    char line1[width+2]; std::strstream buff1(line1,width+1,std::ios::out);
    for(unsigned p=0;(p<width)&&((i+p)<d.size());p+=10) buff1 << std::setw(10) << (p+10+i);
    buff1 << "\n\n" << '\0' ; retString += string(line1);
  }  
  return string(retString);
}
///////////////////////////////////////////////////////////////////
string HorizontalTextBarGraph::operator()(const dVec& d,int width) const {
  double maxVal= 0.0;
  for(unsigned i=0;i<d.size();++i) if(d[i]>maxVal)maxVal=d[i];
  int prec= static_cast<int>(log10(maxVal))+1;
  string retString; retString= string();
  for(unsigned i=0;i<d.size();++i) {
    int nMarks= static_cast<int>(d[i]/maxVal*static_cast<double>(width));
    char line[501]; std::strstream buff(line,500);
    if(!(i%5)) buff << std::setw(3) << i << ' ';
    else    buff << "    ";
    buff << std::setprecision(prec) << std::setw(6) << d[i] << ' ';
    for(unsigned j=0;j<nMarks;++j) buff << '*';
    buff << std::endl << '\0';
    retString += string(line);
  }
  return string(retString);
}
///////////////////////////////////////////////////////////////////
string TextBarGraphComp::operator()(const dVec& d,int width) const {
  double maxVal= 0.0; int N=d.size();
  int window=N/width; if(N%width) ++window; 
  WriteIntToLog()(string("TextBarGraphComp, size= "),N);
  WriteIntToLog()(string("TextBarGraphComp, width= "),width);
  WriteIntToLog()(string("TextBarGraphComp, window= "),window);
  dVec d2; d2.reserve(width);
  for(unsigned i=0;i<N;i+=window) {
    d2.push_back(0.0);
    for(unsigned j=0;(j<window)&&((i+j)<N);++j) if(d[i+j]>=0.0)d2[i/window]+=d[i+j];
    if(d2[i/window]>maxVal)maxVal=d2[i/window];
  }
  int log10mv= log10(maxVal);
  if(maxVal<10)        maxVal=10;
  else if(maxVal<100)  maxVal= (static_cast<int>(maxVal)/10+1)*10; 
  else if(maxVal<1000) maxVal= (static_cast<int>(maxVal)/100+1)*100;
  else             
  WriteIntToLog()(string("maxVal in TextBarGraphComp:"),maxVal);
  double vSteps=10.0; 
  std::string rString; 
  //
  for(unsigned row=0;row<vSteps;++row) {
    if(!(row%2)) { 
      char l[16]; std::strstream b(l,15,std::ios::out); 
      b << std::setw(4) << static_cast<int>(maxVal*(vSteps-row)/vSteps) << "-" << '\0';
      rString += string(l);
    } else
      rString += string("    -");
    for(unsigned i=0;(i<width)&&(i<d2.size());++i) 
      if((vSteps*d2[i]/maxVal)>(vSteps-row-1)) rString += '*';
      else                                     rString += ' ';
    rString += "\n";
  }
  rString += "     ";
  for(unsigned i=0;(i<width)&&(i<d2.size());i+=10) {
    for(unsigned j=0;j<4;++j) rString += '-'; rString += '|';
    for(unsigned j=0;j<4;++j) rString += '-'; rString += '|';
  }
  rString += "\n     ";
  for(unsigned i=0;(i<width)&&(i<d2.size());i+=10) {
    char line[12]; std::strstream buff(line,11,std::ios::out);    
    buff << std::setw(10) << ((i+10)*window) << '\0';
    rString += string(line);
  }
  rString += "\n\n";
  //
  return string(rString);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
  iVec bins_;
  dVec bounds_;
  double minBound_,binSize_;
  int maxCount_;
  Histogram();
 */
//////////////////////////////////////////////////////////////////////////////////////////////////////////
Histogram::Histogram(int nBins,double minBound,double binSize):bins_(),bounds_(),minBound_(minBound),binSize_(binSize),maxCount_(0) {
  if(nBins<2)       { WriteIntToErr()(std::string("nBins < 2 in Histogram constructor"),nBins); return; }
  if(nBins>MaxBins) { WriteIntToErr()(std::string("nBins > maxBins in Histogram constructor"),nBins); return; }
  if(binSize_<=0.0) { WriteFloatToErr()(std::string("binSize_ <= 0.0 in Histogram constructor"),binSize_); return; }
  //
  bins_= iVec(nBins,int(0));
  bounds_= dVec(nBins-1,minBound_);
  for(int i=1;i<nBins-1;++i) bounds_[i]= bounds_[i-1]+binSize_;
  WriteIntToLog()(std::string("Histogram constructor, bins_.size()= "),bins_.size());
  WriteIntToLog()(std::string("Histogram constructor, bounds_.size()= "),bounds_.size());
  WriteFloatToLog()(std::string("Histogram constructor, minBound= "),minBound_);
  WriteIntToLog()(std::string("Histogram constructor, minBin_= "),binSize_);
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
Histogram::Histogram(const dVec& bounds):bins_(),bounds_(bounds),minBound_(),binSize_(-1.0),maxCount_(0) {
  if(bounds_.empty()) { WriteToErr()(std::string("Empty bounds vector in Histogram constructor")); return; }
  bool monotonic=true;
  for(int i=1;i<bounds_.size() && monotonic;++i) monotonic= (bounds_[i-1]<=bounds_[i]); // equal is meaningless, but not an error
  if(!monotonic) { WriteToErr()(std::string("bounds vector is not monotonic increasing in Histogram constructor")); return; }
  // okay, proceed
  bins_= iVec(bounds_.size()+1,int(0));
  minBound_= bounds_[0];
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram::Add(double v) {
  int i=0;
  while(i<bounds_.size() && bounds_[i]<=v) ++i;
  ++bins_[i];
  if(bins_[i]>maxCount_) maxCount_=bins_[i];
  return i;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram::Add(const dVec& v) {
  for(int i=0;i<v.size();++i) Add(v[i]);
  return v.size();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
int Histogram::Collapse(int max) {
  dVec newBounds; iVec newBins;
  int m=(bins_.size()-1)/(max-1);
  if((bins_.size()-1)%(max-1)) ++m;
  binSize_*=m;
  newBins.reserve(bins_.size()/static_cast<double>(m));
  newBounds.reserve(bins_.size()/static_cast<double>(m));
  newBins.push_back(bins_[0]); newBounds.push_back(bounds_[0]); maxCount_=newBins[0];
  int j=1;
  while(j<bins_.size()) {
    int sum=0;
    for(int k=0;k<m && j<bins_.size();++k) sum+=bins_[j++];
    if(sum>maxCount_) maxCount_=sum;
    newBins.push_back(sum);
    if(j<bounds_.size()) newBounds.push_back(bounds_[j-1]);
  }
  bins_=newBins; bounds_=newBounds;
  return bins_.size();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string Histogram::TableList(bool transpose,bool header,int p,int w) const {
  std::string tbl;
  if(!transpose) {
    if(header) tbl= std::string(">=\t<\tN\n");
    for(int i=0;i<bins_.size();++i) {
      if(!i) tbl += std::string("-inf\t");
      else   tbl += DoubleToFixedString()(bounds_[i-1],p) + std::string("\t");
      if(i == bounds_.size()) tbl += std::string("inf\t");
      else                    tbl += DoubleToFixedString()(bounds_[i],p) + std::string("\t");
      tbl += IntToString()(bins_[i]) + std::string("\n");
    }
  } else { // transposed table
    // min row first
    if(header) tbl += std::string(">=\t");
    for(int i=0;i<bins_.size();++i)
      if(!i) tbl += std::string("-inf");
      else   tbl += std::string("\t") + DoubleToFixedString()(bounds_[i-1],p);
    tbl += std::string("\n");
    // max row next
    if(header) tbl += std::string("<\t");
    for(int i=0;i<bins_.size();++i)
      if(!i)                       tbl += DoubleToFixedString()(bounds_[i],p);
      else if(i == bounds_.size()) tbl += std::string("\tinf");
      else                         tbl += std::string("\t")+ DoubleToFixedString()(bounds_[i],p);
    tbl += std::string("\n");
    // count row last
    if(header) tbl += std::string("N\t");
    for(int i=0;i<bins_.size();++i) {
      if(i) tbl += std::string("\t");
      tbl += IntToString()(bins_[i]);
    }
  }
  return tbl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string Histogram::HorizontalBarGraph(int width) const {
  dVec v; v.reserve(bins_.size());
  for(int i=0;i<bins_.size();++i) v.push_back(bins_[i]);
  return HorizontalTextBarGraph()(v,width);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string Histogram::VerticalBarGraph(int width) const {
  dVec v; v.reserve(bins_.size());
  for(int i=0;i<bins_.size();++i) v.push_back(bins_[i]);
  return TextBarGraph()(v,width);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string Histogram::HorizontalBarGraph(int M,int N) const {
  dVec v; v.reserve(bins_.size());   double maxVal= -1.0e10;
  for(int i=0;i<bins_.size();++i) {
    v.push_back(bins_[i]);
    if(v[i]>maxVal) maxVal=v[i];
  }
  int prec= static_cast<int>(log10(maxVal))+1;
  string retString; retString= string();
  for(unsigned i=0;i<v.size();++i) {
    int nMarks= static_cast<int>(v[i]/maxVal*static_cast<double>(M));
    char line[501]; std::strstream buff(line,500);
    if(!(i%5)) buff << std::setw(3) << i << ' ';
    else    buff << "    ";
    buff << std::setprecision(prec) << std::setw(6) << v[i] << ' ';
    for(unsigned j=0;j<nMarks;++j) buff << '*';
    buff << std::endl << '\0';
    retString += string(line);
  }
  return string(retString);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string Histogram::VerticalBarGraph(int M,int N) const {
  dVec v; v.reserve(bins_.size());  double maxVal= -1.0e10;
  for(int i=0;i<bins_.size();++i) {
    v.push_back(bins_[i]);
    if(v[i]>maxVal)maxVal=v[i];
  }
  double vSteps=N;
  string retString; retString= string();
  for(unsigned i=0;i<v.size();i+=M) {
    for(unsigned row=0;row<vSteps;row++) {
      char line2[M+3]; std::strstream buff2(line2,M+2,std::ios::out);
      for(unsigned j=0;(j<M)&&((i+j)<v.size());++j) 
	if((vSteps*v[i+j]/maxVal)>(vSteps-row-1)) buff2 << '*';
	else                                      buff2 << ' ';
      buff2 << std::endl << '\0'; retString += string(line2);
    }
    char line1[M+2]; std::strstream buff1(line1,M+1,std::ios::out);
    for(unsigned p=0;(p<M)&&((i+p)<v.size());p+=10) buff1 << std::setw(10) << (p+10+i);
    buff1 << "\n\n" << '\0' ; retString += string(line1);
  }  
  return string(retString);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
const double EPS1=0.001;
const double EPS2=1.0e-16;
double probks(double alam) {
  int j;
  double a2,fac=2.0,sum=0.0,term,termbf=0.0;
  //
  a2= -2.0*alam*alam;
  for(j=1;j<=100;++j) {
    term=fac*exp(a2*j*j);
    sum += term;
    if(abs(term) <= EPS1*termbf || abs(term) <= EPS2*sum) return sum;
    fac= -fac;
    termbf=abs(term);
  }
  return 1.0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
KS_pairtest::KS_pairtest(const dVec& dd1,const dVec& dd2):data1_(dd1),data2_(dd2),d_(-1.0),p_(-1.0) {
  //  if(data1_.size()!=data2_.size()) {
  //    WriteToLogs()(std::string("vector sizes don't agree in KS_pairtest"));
  //    return;
  //  }
  std::sort(data1_.begin(),data1_.end());
  std::sort(data2_.begin(),data2_.end());
  double fn1=0.0,fn2=0.0,en1=data1_.size(),en2=data2_.size(),en,d1,d2,dt;
  d_=0.0;
  int j1=0,j2=0;
  while(j1<data1_.size() && j2<data2_.size()) {
    if((d1=data1_[j1]) <= (d2=data2_[j2])) fn1= j1++/en1;
    if(d2<=d1) fn2=j2++/en2;
    if((dt=fabs(fn2-fn1)) > d_) d_=dt;
  }
  en=sqrt(en1*en2/(en1+en2));
  p_=probks((en+0.12+0.11/en)*d_);
  return;
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////
const int ITMAX=500;
const double EPS=3.0e-10;
void gser(double *gamser,double a,double x, double *gln) {
  int n;
  double sum,del,ap;
  //
  *gln=gammln()(a);
  if(x<=0.0) {
    if(x < 0.0) { WriteToLogs()(std::string("x less than 0 in routine gser")); exit(-1); }
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for(n=1;n<=ITMAX;++n) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    WriteToLogs()(std::string("a too large, ITMAX too small in routine gser")); 
    WriteFloatToErr()(std::string("a= "),a);
    WriteFloatToErr()(std::string("x= "),x);
    exit(-1);
    return;
  }
}
//
const double FPMIN=1.0e-30;
void gcf( double *gammcf, double a, double x, double *gln ) {
  int i;
  double an,b,c,d,del,h;
  //   
  *gln = gammln()( a );
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  //
  for( i = 1; i <= ITMAX; ++i ) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
       
    if ( fabs( d ) < FPMIN ) d = FPMIN;
    c = b + an / c;
    if ( fabs( c ) < FPMIN ) c = FPMIN;
    d   = 1.0 / d;
    del = d * c;
    h  *= del;
    if ( fabs( del - 1.0 ) < EPS) break;
  }

  if (i > ITMAX) { 
    WriteToLogs()(std::string("a too large, ITMAX too small in gcf")); 
    WriteFloatToErr()(std::string("a= "),a);
    WriteFloatToErr()(std::string("x= "),x);
    exit(-1); }
  *gammcf = exp( -x + a * log( x ) -( *gln ) ) * h;
}
//
double gammq(double a,double x) {
  double gamser,gammcf,gln;

  if(x<0.0 || a<=0.0) { WriteToLogs()(std::string("Invalid arguments in routine gammq")); exit(-1); }
  if(x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}
/*
class ChiSq_pairtest {
 private: 
  dVec data1_,data2_; bool cons_;
  double x2_,p_; int df_;
 public:
  ChiSq_pairtest(const dVec& d1, const dVec& d2, bool c=false);
 */
/////////////////////////////////////////////////////////////////////////////////////////////////////
ChiSq_pairtest::ChiSq_pairtest(const dVec& dd1,const dVec& dd2,bool con):data1_(dd1),data2_(dd2),cons_(con),df_(0),x2_(0.0),p_(1.0) {
  if(data1_.size()!=data2_.size()) {
    WriteToLogs()(std::string("vector sizes don't agree in ChiSq_pairtest")); return;
  }
  double R=VectorSum()(data1_),S=VectorSum()(data2_);
  int nBins= data1_.size();
  if(!cons_) df_=nBins;
  else       df_=nBins-1;
  //
  double rootSR= sqrt(S/R),rootRS= sqrt(R/S);
  for(int i=0;i<nBins;++i) 
    if(data1_[i]==0.0 && data2_[i]==0.0) 
      --df_;
    else {
      double temp= rootSR*data1_[i] - rootRS*data2_[i];
      x2_ += temp*temp/(data1_[i]+data2_[i]);
    }

  //  WriteFloatToErr()(std::string("R= "),R);
  //  WriteFloatToErr()(std::string("S= "),S);
  //  WriteFloatToErr()(std::string("chi sqr= "),x2_);
  //  WriteFloatToErr()(std::string("deg fr = "),df_);
  p_= gammq(0.5*df_,0.5*x2_);
}
