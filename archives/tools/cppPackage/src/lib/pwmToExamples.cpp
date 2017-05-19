// fa2Dust.cpp : Defines the entry point for the console application.
//
#include <unistd.h>
#include <time.h>
#include <math.h>
#ifdef INTEL
#include <mathimf.h>
#endif
#include <cstdlib>
#include <iostream>
#include <strstream>
#include <iomanip>
#include <fstream>
#include <map>
#include <functional>
#include "../lib/fastaSeq.h"
#include "../lib/basicSeq.h"
#include "../lib/logFile.h"
#include "../lib/arg.h"
#include "../lib/util.h"
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>
typedef std::vector<double> dVec;
typedef std::vector<dVec> dMat;
typedef std::vector<dMat> dTen;
typedef vector<std::string> strVec;
typedef std::map<int,double,std::less<int> > idMap;
typedef std::vector<idMap> idMapVec;
//
bool oFile=false,empModel=false,verbose=true,sepFiles=false;
//int inMu=100; double inSig=5, 
long seed1; int iseed1=0;
std::string ofPrefix("pwm"); int nBases=4,nEx=1000; 
//
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  //  if(a.GetFloatValue(std::string("w"),std::string("insert sigma"),inSig,-1,100,true,true)) ++pCount;
  if(a.GetString(std::string("o"),std::string("output prefix"),ofPrefix)) { oFile=true; ++pCount; }
  if(a.KeyPresent(std::string("quiet"))) { ++pCount; verbose=false; }
  if(a.KeyPresent(std::string("sep"))) { ++pCount; sepFiles=true; }
  if(a.GetIntValue(std::string("n"),std::string("nEx"),nEx,1,1000000)) ++pCount; 
  if(a.GetIntValue(std::string("s"),std::string("seed1"),iseed1,-1000000,1000000)) { ++pCount; seed1=iseed1; }
  else { time_t t0; time(&t0); seed1= -1*t0; }
  if(pCount!=a.numOptArgs()) {
    WriteToErr()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteToLog()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLog()(std::string("processed= "),pCount);
    WriteIntToLog()(std::string("entered= "),a.numOptArgs());
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ProbVec {
private:
  dVec p_,c_;
public:
  ProbVec::ProbVec(const dVec& p) {
    int N=p.size(); p_.reserve(N); c_.reserve(N);
    if(!N) return;
    if(p[0]<0.0) return;
    p_.push_back(p[0]); c_.push_back(p[0]);
    for(int i=1;i<N;++i) 
      if(p[i]<0.0) return;
      else { p_.push_back(p[i]); c_.push_back(c_[i-1]+p[i]); }
    double sum=c_.back();
    for(int i=0;i<N;++i) { p_[i] /= sum; c_[i] /=sum; }
    return;
  }
  ProbVec::ProbVec() {}
  ProbVec::ProbVec(const ProbVec& o):p_(o.p_),c_(o.c_) {}
  ProbVec& operator=(const ProbVec& o) { p_=o.p_; c_=o.c_; return *this; }
  dVec pVector() const { return p_; }
  double operator[](unsigned i) const { if(i<p_.size()) return p_[i]; else return -1.0; }
  int rndIdx(double p) const{
    int r=0; while(c_[r]<p && r<c_.size()) ++r; 
    return r;
  }
};
typedef std::vector<ProbVec> ProbMat;
typedef std::vector<ProbMat> ProbTen;
///////////////////////////////////////////////////////////////
bool LoadSeqModel(const std::string& f,dTen& m,strVec& label) {
  std::ifstream mFile(f.c_str(),std::ios::in);
  if(mFile.fail()) return false;
  m.clear(); label.clear();
  char lineIn[501];
  // read to the first def line
  while(!mFile.eof()&&lineIn[0]!='>') mFile.getline(lineIn,500);
  while(!mFile.eof()) {
    dMat tempMat; 
    std::string line(lineIn),tempLabel;
    if(line.length()) tempLabel= line.substr(1,line.length()-1);
    mFile.getline(lineIn,500); line= std::string(lineIn);
    while(line[0]!='>' && !mFile.eof()) {
      strVec fields= split()(line);
      if(fields.size()>=nBases && line[0] != '#') {
	dVec temp(nBases,0.0);
	double sum=0.0;
	for(int i=0;i<4;++i) {
	  temp[i]= atof(fields[i].c_str());
	  sum += temp[i];
	}
	for(int i=0;i<4;++i) temp[i]/=sum;
	tempMat.push_back(temp);
      }
      mFile.getline(lineIn,500);
      line= std::string(lineIn);
    }
    if(tempMat.size()) {
      label.push_back(tempLabel);
      m.push_back(tempMat);
    }
    /*
    if(line.empty()) continue;
    if(line[0]=='>') { label = line.substr(1,line.length()-1); continue; }
    if(line[0]=='#') continue;
    strVec fields= split()(lineIn,' '); 
    if(fields.size()!=nBases) continue;
    dVec temp; double sum=0.0;
    for(int i=0;i<fields.size();++i) { temp.push_back(std::atoi(fields[i].c_str())); sum += temp[i]; }
    // normalize
    for(int i=0;i<temp.size();++i) temp[i]/=sum; 
    m.push_back(temp);
    */
  } 
  mFile.close();
  return true;
}
///////////////////////////////////////////////////////////////
std::string MakeSig(const dMat& m,long& seed) {
  if(m.empty()) return std::string();
  std::string rs;
  for(dMat::const_iterator it=m.begin();it!=m.end();++it) {
    int idx= DrawFromProbVector()(*it,seed);
    rs += DNAbases[idx];
  }
  return rs;
}
///////////////////////////////////////////////////////////////
std::string MakeInsert(const ProbMat& m,gsl_rng *rnd) {
  std::string ins;
  for(ProbMat::const_iterator it=m.begin();it!=m.end();++it) {
    int idx= it->rndIdx(gsl_rng_uniform(rnd));
    ins += DNAbases[idx];
  }
  return ins;
}
///////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
  if(!LogSetup()(argc,argv)) { std::cout << "error opening logs" << std::endl; exit(-1);}
  //
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<2)||(arg1==std::string("-h"))||(arg1==std::string("-help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]" << std::endl;
    std::cout << "       " << argv[0] << " modfilename [-o string] [-n int] [-s int] [-sep] [-quiet]" << std::endl; 
    WriteToLog()(std::string("usage message displayed"));
    WriteToErr()(std::string("usage message displayed"));
    exit(0); 
  }
  //
  cjArgs arguments(argc,argv,1);
  ProcessArgs(arguments);
  std::string mfname= arguments.RequiredArgs()[0];
  //// Load the sequence content models
  std::cout << "Loading sequence content models from file " << mfname << "..."; std::cout.flush();
  dTen seqModels; strVec seqLabels;
  if(!LoadSeqModel(mfname,seqModels,seqLabels)) {
    WriteToLogs()(std::string("error sequence model from ")+mfname);
    std::cout << "\n** Error loading sequence models from file " << mfname << ". exiting. **\n";
    exit(-1);
  }
  int nSeqModels= seqModels.size();
  std::cout << "done. " << nSeqModels << " models loaded.\n";
  //// final sanity check
  if(!nSeqModels) {
    WriteToLogs()(std::string("no sequence models"));
    std::cout << "** No sequence models found! exiting. **\n";
    exit(-1);
  }
  //// now convert the raw data to probability vectors (vectors and matrixes of type ProbVec)
  std::cout << "Converting models..."; std::cout.flush();
  int nModels=nSeqModels;//,L=seqs.LongestSeq();
  ProbTen seqModelMat; seqModelMat.reserve(nModels);
  for(int m=0;m<nModels;++m) {
    ProbMat tempMat; 
    tempMat.reserve(seqModels[m].size());
    for(int i=0;i<seqModels[m].size();++i)
      tempMat.push_back(ProbVec(seqModels[m][i]));
    seqModelMat.push_back(tempMat);
  }
  std::cout << "done.\n";
  // dump the models to test for correct read
  if(verbose) {
    std::cout << "\nSequence Models:\n";
    for(int m=0;m<nModels;++m) {
      std::cout << "model " << m << "; label= " << seqLabels[m] << "\n\tC\tT\tA\tG\n";
      for(int p=0;p<seqModelMat[m].size();++p) {
	std::cout << (p+1);
	for(int j=0;j<nBases;++j) std::cout << '\t' << seqModelMat[m][p][j];
	std::cout << '\n';
      }
    }
  }
  //////////// start gsl initialization ////////////
  time_t t0; time(&t0);
  const gsl_rng_type *T;
  gsl_rng *rnd;
  gsl_rng_env_setup();
  T= gsl_rng_default;
  gsl_rng_default_seed= t0%10000;
  rnd= gsl_rng_alloc(T);
  ///////////// end gsl initialization /////////////
  std::string ofname(ofPrefix);
  std::string sepFileLabel(".A");
  if(sepFiles) ofname += sepFileLabel;
  ofname += std::string(".logoEx.txt");
  std::ofstream ofile(ofname.c_str(),std::ios::out);
  for(int m=0;m<nModels;++m) {
    ofile << "#model " << m << ':' << seqLabels[m] << '\n';
    for(int x=0;x<nEx;++x) {
      std::string signal= MakeInsert(seqModelMat[m],rnd);
      ofile << signal << '\n';
    }
    ofile << '\n';
    if(sepFiles && (m+1)<nModels) {
      ofile.close();
      ++sepFileLabel[1];
      ofname= ofPrefix + sepFileLabel + std::string(".logoEx.txt");
      ofile.open(ofname.c_str(),std::ios::out);
    }
  }
  ofile.close();
  //
  WriteToLogTS()(std::string("Normal termination."));
  WriteToErrTS()(std::string("Normal termination."));
  return 0;
}
