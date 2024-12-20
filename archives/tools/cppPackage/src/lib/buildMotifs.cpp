/////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <limits>
#include <gsl/gsl_rng.h>
#include "../lib/arg.h"
#include "../lib/logFile.h"
#include "../lib/util.h"
#include "../lib/basicSeq.h"
#include "../lib/vectorF.h"
/////////////////////////////////////////////////////////////////////////////////
typedef std::vector<std::string> strVec;
typedef std::vector<double> dVec;
typedef std::vector<dVec> dMat;
//AddLog2 log2plus();
/////////////////////////////////////////////////////////////////////////////////
const int A=0,C=1,G=2,T=3,U=3;
/////////////////////////////////////////////////////////////////////////////////
inline int baseIndex(char c) {
  if(c=='A' || c=='a') return A;
  if(c=='C' || c=='c') return C;
  if(c=='G' || c=='g') return G;
  if(c=='T' || c=='U' || c=='t' || c=='u') return T;
  return -1;
}
///////////////////////////////////////////////////////////////////////
bool noLogs=false,quiet=false,forceOvlp=false,useChi2=false; 
int minOvlp=2,maxIt=1000,platIt=25,nStarts=20,minMotifL=4,maxMotifL=20,simCounts=1000;
std::string ofPrefix("motifs.build"),bgFilename;
double initWeight=0.1,minWeightFrac=0.001,wgtFudge=1.0;;
///////////////////////////////////////////////////////////////////////
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  if(a.KeyPresent(std::string("nologs"))) { noLogs=true; LogsOff(); ++pCount; }
  if(a.KeyPresent(std::string("quiet")) || a.KeyPresent(std::string("Q"))) { quiet=true; ++pCount; }
  if(a.KeyPresent(std::string("ovlp"))) { forceOvlp=true; ++pCount; }
  if(a.KeyPresent(std::string("X2"))) { useChi2=true; ++pCount; }
  if(a.KeyPresent(std::string("multinomial"))) { useChi2=false; ++pCount; }
  if(a.GetIntValue(std::string("I"),std::string("max It"),maxIt,10,100000,true,true)) ++pCount;
  if(a.GetIntValue(std::string("P"),std::string("plateau Iterations"),platIt,10,100000,true,true)) ++pCount;
  if(a.GetIntValue(std::string("S"),std::string("n starts"),nStarts,0,10000,true,true)) ++pCount;
  if(a.GetIntValue(std::string("M"),std::string("max motif L"),maxMotifL,4,100,true,true)) ++pCount;
  if(a.GetIntValue(std::string("N"),std::string("simulated counts"),simCounts,100,100000,true,true)) ++pCount;
  if(a.GetIntValue(std::string("m"),std::string("min motif L"),minMotifL,3,maxMotifL,true,true)) ++pCount;
  if(a.GetIntValue(std::string("x"),std::string("min word-motif overlap"),minOvlp,1,100,true,true)) ++pCount;
  if(a.GetFloatValue(std::string("w"),std::string("bg weight"),initWeight,0.0,1.0,true,true)) ++pCount;
  if(a.GetFloatValue(std::string("f"),std::string("weight ff"),wgtFudge,0.0,3.0,true,true)) ++pCount;
  if(a.GetString(std::string("o"),std::string("output prefix"),ofPrefix)) ++pCount;
  if(a.GetString(std::string("B"),std::string("bgFilename"),bgFilename)) ++pCount;
  if(pCount!=a.numOptArgs()) {
    WriteToLogs()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLogs()(std::string("entered= "),a.numOptArgs());
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// routine that opens a statistics file (as generated by fa2Stats) and generates a 0-order (nucleotide) base-2 log vector of frequencies
dVec logBGVector(const std::string& fname) {
  dVec bg;
  std::ifstream inFile(fname.c_str(),std::ios::in);
  char line[101];
  inFile.getline(line,100);
  // read past
  while(!inFile.eof() && line[0]=='#') inFile.getline(line,100);
  for(int i=0;i<4;++i) {
    strVec f= split()(std::string(line),'\t');
    if(f.size()<3) return dVec();
    bg.push_back(atof(f[2].c_str()));
    inFile.getline(line,100);
  }
  inFile.close();
  dVec logBG;
  for(int i=0;i<4;++i) logBG.push_back(log2(bg[i]+1.0e-10));
  return logBG;
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
  int rndIdx(double p) {
    int r=0; while(c_[r]<p && r<c_.size()) ++r; 
    return r;
  }
};
/////////////////////////////////////////////////////////////////////////////////
// July 2007- trying an experiment- instead of averaging, use the maximum value instead...
// didn't work;  reverting
double totalLogP(const std::string& w,const dMat& logM, const dVec& logBG) {
  double totalP=0.0,maxP= -1.0e100;;
  int K=w.length(),N=logM.size();
  for(int p=(minOvlp-K);p<=N-minOvlp;++p) {
    double logP= 0.0;
    for(int j=0;j<K;++j) 
      if((p+j)<0 || (p+j)>=N) logP += logBG[baseIndex(w[j])];
      else                    logP += logM[p+j][baseIndex(w[j])];
    totalP += pow(2.0,logP);
    if(logP>maxP) maxP=logP;
  }
  //  return log2(totalP/static_cast<double>(N-K+1));
  return log2(totalP/static_cast<double>(N-2*minOvlp+K+1));
  //return maxP;
}
/////////////////////////////////////////////////////////////////////////////////
int bestPosP(const std::string& w,const dMat& logM,const dVec& logBG) { 
  int K=w.length(),N=logM.size();
  double bestP=-1000; int bestPos=0;
  for(int p=minOvlp-K;p<=N-minOvlp;++p) {
    double logP= 0.0;
    for(int j=0;j<K;++j) 
      //if((p+j)>=0 && (p+j)<N) logP += logM[p+j][baseIndex(w[j])]-logBG[baseIndex(w[j])];
      if((p+j)<0 || (p+j)>=N) logP += logBG[baseIndex(w[j])];
      else                    logP += logM[p+j][baseIndex(w[j])];
    if(logP > bestP) { bestP=logP; bestPos=p; }
    //else if(logP==bestP && abs(p-K+1)<abs(bestPos-K+1)) { bestP=logP; bestPos=p; }
  }
  return bestPos;
}
/////////////////////////////////////////////////////////////////////////////////
dVec posP(const std::string& w,const dMat& logM,const dVec& logBG) {
  int K=w.length(),N=logM.size();
  dVec pp(N+K-2*minOvlp+1,0.0);
  for(int p=minOvlp-K;p<=N-minOvlp;++p) {
    double logP= 0.0;
    for(int j=0;j<K;++j) 
      //  if((p+j)>=0 && (p+j)<N) logP += logM[p+j][baseIndex(w[j])]-logBG[baseIndex(w[j])];
      if((p+j)<0 || (p+j)>=N) logP += logBG[baseIndex(w[j])];
      else                    logP += logM[p+j][baseIndex(w[j])];
    pp[p+K-minOvlp]=pow(2.0,logP);
  }
  return pp;
}
/////////////////////////////////////////////////////////////////////////////////
void updateModel(const std::string& w,int pos,double val,dMat& logM,bool improve) {
  for(int k=0;k<w.length();++k) 
    if((pos+k)>=0 && (pos+k)<logM.size()) { // seed in the word
      // pull the existing nucleotide vector first
      dVec temp(4,0.0); for(int kk=0;kk<4;++kk) temp[kk]= pow(2.0,logM[pos+k][kk]);
      // add (or subtract) val (set externally) to the indicated base
      double fac= improve?1.0:-1.0;
      temp[baseIndex(w[k])] += val*fac;
      double sum=1.0+fac*val;
      if(temp[baseIndex(w[k])]<=0.0) {
	temp[baseIndex(w[k])]=1.0e-10;
	sum=0.0;
	for(int i=0;i<4;++i) sum+=temp[i];
      }
      for(int kk=0;kk<4;++kk) logM[pos+k][kk] = log2(temp[kk]/sum);
    }
  return;
}
/////////////////////////////////////////////////////////////////////////////////
typedef struct wgtWord{
  std::string word; double weight,log2weight;
  wgtWord(std::string w1,double w2):word(w1),weight(w2) {log2weight=log2(weight/100.0);}
  bool operator< (const wgtWord& o) const { return weight>o.weight; }
};
typedef std::vector<wgtWord> wgtWordVec;
wgtWordVec MakeWeightedWordVector(const strVec& words,const dVec& weights) {
  if(words.size()!=weights.size()) return wgtWordVec();
  wgtWordVec wwv; wwv.reserve(words.size());
  for(int i=0;i<words.size();++i) wwv.push_back(wgtWord(words[i],weights[i]));
  std::sort(wwv.begin(),wwv.end());
  wgtWordVec temp(wwv);
  wwv.clear();
  double minWeight= temp[0].weight*minWeightFrac;
  for(int i=0;temp[i].weight>=minWeight && i<temp.size();++i) wwv.push_back(temp[i]);
  return wwv;
}
typedef std::vector<wgtWordVec> wgtWordMat;
/////////////////////////////////////////////////////////////////////////////////
int main(int argc,char *argv[]) {
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<2)||(arg1==std::string("-h"))||(arg1==std::string("-help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]\n";
    std::cout << "       " << argv[0] << " filename [-o filestring] [-r int] [-l int] [-n int] [-d] [-quiet] [-nologs]\n"; 
    exit(0); 
  }
  if(!LogSetup()(argc,argv)) { std::cout << "error creating logs. exiting.\n"; exit(-1); }
  //////////// start gsl initialization ////////////
  time_t t0; time(&t0);
  const gsl_rng_type *T;
  gsl_rng *rnd;
  gsl_rng_env_setup();
  T= gsl_rng_default;
  gsl_rng_default_seed= t0%10000;
  rnd= gsl_rng_alloc(T);
  ///////////// end gsl initialization /////////////
  cjArgs arguments(argc,argv,1);
  ProcessArgs(arguments);
  std::string fname(arguments.RequiredArgs()[0]);
  //DumpParameters(std::cout,argc,argv);
  //std::ofstream test("testDump.txt",std::ios::out);
  //DumpParameters(test,argc,argv);
  //test.close();
  // Load the file
  std::ifstream inFile(fname.c_str(),std::ios::in);
  char lineIn[10001];
  inFile.getline(lineIn,10000);
  strVec f= split()(std::string(lineIn));
  if(f.empty()) {
    std::cout << "invalid entry in " << fname << "; exiting\n"; exit(-1);
  }
  // the first line will have a bunch of information we need
  int K= f[0].length(),r=f.size()-1; // K= word size, r= # of bases
  if(minOvlp>K) minOvlp=K;
  if(minOvlp<1) minOvlp=1;
  //  if(minMotifL<K) minMotifL=K;
  //if(maxMotifL<K) maxMotifL=K;
  int M= pow(4,K);                   // M= max number of words
  // set initial matrixes
  strVec words; words.reserve(M);   // word vector
  dMat baseV(r,dVec());             // baseVector has base vector weights
  for(int i=0;i<r;++i) baseV[i].reserve(M); // up to M entries for each vector
  dMat rawBG(r,dVec(5,0.0));        // rawBG will contain the raw BG nucleotide counts 
  // now loop to end of file
  // 01 August 2007 a fudge fix-  the vanilla NMF is producing a "flattened" motif, as the 
  // NMF weight matrix is not accurately reproducing the counts from a test set.  the best fix
  // is likely the non-smooth NMF, but in the meantime, introduce a fudge factor, a power by which
  // all values in a column (all weights for a given matrix) will be exponentiated- in effect
  // making the distribution more extreme
  int m=0;
  dVec wgtSum(r,0.0);
  dMat tempWgtMatrix(r,dVec());
  while(!inFile.eof() && ++m <= M) {
    if(f.size()==(r+1) && f[0].length()==K) {
      words.push_back(f[0]); 
      for(int j=0;j<r;++j) {
	double weight= atof(f[j+1].c_str())/100.0;
	if(wgtFudge != 1.0) {
	  weight= pow(weight,wgtFudge);
	  wgtSum[j] += weight;
	}
	tempWgtMatrix[j].push_back(weight);
      }
    }
    inFile.getline(lineIn,10000);
    f= split()(std::string(lineIn));
  }
  inFile.close();
  // now, continuing the fudge fix, if the wgtFudge factor is not 1.0, we need to renormalize each vector	
  if(wgtFudge != 1.0) 
    for(int j=0;j<r;++j)
      for(int m=0;m<tempWgtMatrix[j].size();++m)
	tempWgtMatrix[j][m] /= wgtSum[j];
  // end of fudge fix; now to create the other necessary data, baseV and rawBG
  baseV= dMat(r,dVec(tempWgtMatrix[0].size(),0.0));
  for(int j=0;j<r;++j) {
    for(int m=0;m<tempWgtMatrix[j].size();++m) {
      baseV[j][m]= static_cast<double>(simCounts)*tempWgtMatrix[j][m];
      for(int w=0;w<K;++w) {
	double d = tempWgtMatrix[j][m]*static_cast<double>(simCounts);
	rawBG[j][baseIndex(words[m][w])] += d;
	rawBG[j][4] +=  d;
      }
    }
  }
  //while(!inFile.eof() && ++m <= M) {
  //  if(f.size()==(r+1) && f[0].length()==K) {
  //    words.push_back(f[0]); 
  //    for(int j=0;j<r;++j) {
  //	double weight= static_cast<double>(simCounts)/100.0*atof(f[j+1].c_str());
  //	if(wgtFudge != 1.0) {  
  //	  weight= pow(weight,wgtFudge);
  //	  wgtSum[j] += weight;
  //	}
  //	baseV[j].push_back(weight);
  //	for(int w=0;w<K;++w) { 
  //	  rawBG[j][baseIndex(f[0][w])] += weight;
  //	  rawBG[j][4]+=weight;
  //	}
  //      }
  //    }
  //    inFile.getline(lineIn,10000);
  //    f= split()(std::string(lineIn));
  //  }
  //  inFile.close();
  ////
  //if(wgtFudge != 1.0) {
  //  for(int j=0;j<r;++j)
  //    for(int m=0;m<baseV[j].size();++m)
  //	baseV[j][m] /= wgtSum[j];
  //}
  //////////// dump loading to output for a check /////////////////
  std::cout << m << " lines read, " << words.size() << " records generated\n\nbackground models:\n";
  std::cout << "           A         C         G         T\n";
  // background is either from an input file (generated by fa2Stats) or empiric from the word lists
  dMat logBG(r,dVec(4,0.0));
  dVec fLogBG;
  if(bgFilename.length()) fLogBG= logBGVector(bgFilename);
  for(int i=0;i<r;++i) {
    std::cout << i << ':';
    if(fLogBG.size()) 
      logBG[i]= fLogBG;
    else
      for(int j=0;j<4;++j) {
	logBG[i][j]= log2(rawBG[i][j]/rawBG[i][4]);
      }
    for(int j=0;j<4;++j) std::cout << std::setw(10) << std::setprecision(4) << logBG[i][j];
    std::cout << "\n  ";
    for(int j=0;j<4;++j) std::cout << std::setw(10) << rawBG[i][j];
    std::cout << "\n";
  }
  std::cout << "\ntop 10 words of each base vec\n";
  wgtWordMat wwm(r,wgtWordVec());
  // create the matrix that gives the weighting of each word for each word in each motif- this is 
  // the maximum likelihood probability of the probability for that word in this model
  dMat weightMat(r,dVec());
  for(int i=0;i<r;++i) {
    wwm[i]= MakeWeightedWordVector(words,baseV[i]);
    std::cout << i << ":"; 
    double weightSum=0.0; int wCount=wwm[i].size();
    for(int j=0;j<wCount;++j) {
      if(j<10) std::cout << wwm[i][j].word << "[" << std::setprecision(5) << wwm[i][j].weight << "] ";
      weightSum += wwm[i][j].weight;
    }
    weightMat[i].reserve(wCount);
    for(int j=0;j<wCount;++j) weightMat[i].push_back(log2(wwm[i][j].weight/weightSum));
    std::cout << "\n";
  }
  std::cout.flush();
  /////////////////// end of output dump ///////////////////      
  // make the probability vector for each base
  ProbVec *baseProbs= new ProbVec[r];  
  for(int i=0;i<r;++i) {
    dVec scaledBaseVec(baseV[i]);
    for(int ii=0;ii<scaledBaseVec.size();++ii) scaledBaseVec[ii]= pow(baseV[i][ii],2.0);
    baseProbs[i]= ProbVec(scaledBaseVec);
  }
  // now, go through each base vector to try to generate a motif
  std::string ofname(ofPrefix+std::string(".motifs"));
  std::ofstream mFile(ofname.c_str(),std::ios::out);
  dVec bestMotifScore(r,-1.0e10);
  // motif lengths are sampled uniformly between minMotifL and maxMotifL (inclusive)
  dVec LengthProbVec(maxMotifL-minMotifL+1,1.0);
  ProbVec lengthProbs(LengthProbVec);
  // check file
  std::string cfname(ofPrefix+std::string(".check.txt"));
  std::ofstream cFile(cfname.c_str(),std::ios::out);
  // loop over all r motif models in the input file
  for(int i=0;i<r;++i) {
    std::cout << "# ** building motif " << (i+1) << " of " << r << " **\n";  
    // start with a matrix of bg-vectors that are a weighted average of 1/4 each and the empiric background
    dVec initBG(4,0.25 * (1.0-initWeight)); 
    for(int k=0;k<4;++k) initBG[k]= log2(initBG[k]+ initWeight*pow(2.0,logBG[i][k]));
    // first generate the bg prob vector, which is the probability of observing word wwm[i][j].word under
    // the bg nucleotide frequencies
    dVec bgProb(wwm[i].size(),0.0); double totalBGprob=NEG_INF;
    for(int j=0;j<wwm[i].size();++j) {
      for(int k=0;k<wwm[i][j].word.length();++k)
	//bgProb[j] += logBG[i][baseIndex(wwm[i][j].word[k])];
	bgProb[j] += initBG[baseIndex(wwm[i][j].word[k])];
      totalBGprob = AddLog2()(totalBGprob,bgProb[j]);
      //      totalBGprob += pow(2.0,bgProb[j]);
    }
    //normalize to using only the words in the list
    //totalBGprob = log2(totalBGprob);
    for(int j=0;j<bgProb.size();++j) bgProb[j] -= totalBGprob;
    dMat maxLogModel;//(motifSize,initBG); 
    double maxMultiProb=-1.0e10,maxTotalModProb,bestTotalModProb;
    double minChi2=1e100;
    dVec bestModProb,maxModProb;
    for(int s=0;s<nStarts;++s) {
      int motifSize= minMotifL + lengthProbs.rndIdx(gsl_rng_uniform(rnd));
      dMat logModel(motifSize,initBG);
      // July 2007-  try guiding the initial guess based on the best word in the set
      // if the motif is longer than the word, choose a random place to seed the word
      int seedStart=0;
      if(motifSize > K) 
	seedStart= static_cast<int>(gsl_rng_uniform(rnd) * static_cast<double>(motifSize-K));	
      for(int j=0;j<motifSize && j<K;++j) {
	dVec newStart(4,-5.05889368905);
	newStart[baseIndex(wwm[i][0].word[j])] = -0.13606154958;
	logModel[j+seedStart]= newStart; 
      }
      // end of motif guess
      // get the total probability for each word in the current model
      dVec modProb(wwm[i].size(),0.0);
      for(int j=0;j<wwm[i].size();++j) modProb[j]= totalLogP(wwm[i][j].word,logModel,initBG);
      // calculate the total distribution probability (logarithm of a multinomial without the constant)
      double totalModProb=NEG_INF;
      for(int j=0;j<modProb.size();++j) //totalModProb+= pow(2.0,modProb[j]); 
	totalModProb = AddLog2()(totalModProb,modProb[j]);
	//totalModProb= log2(totalModProb);
      double logDistP=0.0,chi2=0.0;
      double E=0.0,M=0.0;
      for(int j=0;j<wwm[i].size();++j) {
	logDistP += wwm[i][j].weight * (modProb[j]-totalModProb-bgProb[j]);
	E = pow(2.0,weightMat[i][j]); M = pow(2.0,modProb[j]-totalModProb);
	chi2 += (E-M)*(E-M)*1.0/E;//(E+M);
	//std::cout << "E= " << E << ",M= " << M << ",E-M= " << (E-M) << ",chi2= " << chi2 << '\n' ;
      }
      //chi2=1e6;
      //std::cout << "initial logDistP= " << logDistP << ", initial chi2= " << chi2 << ", totalModProb= " << totalModProb << '\n';
      // here's the guts of it.  Select a word to improve based on the observed probabilities
      dMat oldLogModel(logModel),bestLogModel(logModel); 
      double oldLogDistP=logDistP,bestLogDistP=logDistP,oldChi2=chi2,bestChi2=chi2;
      int plateau=0,iter=0;
      while(iter<maxIt && plateau<platIt) {
	if(iter) {
	  totalModProb=NEG_INF;
	  for(int j=0;j<modProb.size();++j) //totalModProb+= pow(2.0,modProb[j]); 
	    totalModProb= AddLog2()(totalModProb,modProb[j]);
	  //	  totalModProb= log2(totalModProb);
	}
	// new criterion for selecting word to update- weighted based on absolute difference in weight fraction and probability fraction
	//      double f1= pow(2.0,weightMat[i][j]),f2=pow(2.0,maxModProb[j]-maxTotalModProb);
	dVec pDiff(wwm[i].size(),0.0);
	E=0.0; M=0.0;
	for(int j=0;j<pDiff.size();++j) {
	  E = pow(2.0,weightMat[i][j]); M = pow(2.0,modProb[j]-totalModProb);
	  //if(iter % 2) 
	  pDiff[j]= (E-M) * (E-M);
	  //else         
	  //  pDiff[j]= (E-M)*(E-M)/E;
	  //pDiff[j]= fabs(pow(2.0,weightMat[i][j])-pow(2.0,modProb[j]-totalModProb));
	}
	ProbVec selector(pDiff);
	// select the word to use in the update, based on the baseProb vector for this model (weighted towards higher ranked words)
	int pos,which;
	//which= baseProbs[i].rndIdx(gsl_rng_uniform(rnd));
	which= selector.rndIdx(gsl_rng_uniform(rnd));
	// select the position in the model at which the word will be updated, sampled according to the probability of the word occuring there
	dVec pp= posP(wwm[i][which].word,logModel,initBG);
	for(int kk=0;kk<pp.size();++kk) pp[kk]= pow(pp[kk],2.0);
	ProbVec P(pp);//posP(wwm[i][which].word,logModel));
	pos= P.rndIdx(gsl_rng_uniform(rnd)) + minOvlp - K;
	updateModel(wwm[i][which].word,pos,wwm[i][which].weight/static_cast<double>(simCounts),logModel,(modProb[which]-totalModProb) < weightMat[i][which]);
	//
	for(int j=0;j<wwm[i].size();++j) {
	  modProb[j]= totalLogP(wwm[i][j].word,logModel,initBG);
	  //if(j<10) std::cout << wwm[i][j].word << '\t' << wwm[i][j].log2weight << '\t' << bgProb[j] << '\t' << 
	  //	     modProb[j] << '\t' << (wwm[i][j].log2weight-modProb[j]) << '\n';
	}
	oldLogDistP=logDistP; oldLogModel=logModel; oldChi2=chi2;
	logDistP=0.0; chi2=0.0;
	E=0.0; M=0.0;
	for(int j=0;j<wwm[i].size();++j) {
	  logDistP += wwm[i][j].weight * (modProb[j]-totalModProb-bgProb[j]);
	  E = pow(2.0,weightMat[i][j]); M = pow(2.0,modProb[j]-totalModProb);
	  chi2 += 1.0*(E-M)*(E-M)/E;//(E+M);
	}
	//std::cout << "chi2= " << chi2 << ", oldChi2= " << oldChi2 << ", bestChi2= " << bestChi2 << '\n';
	if((useChi2 && chi2<bestChi2) || ( !useChi2 && (logDistP>bestLogDistP))) { 
	  plateau=0;
	  bestLogDistP=logDistP; bestLogModel=logModel; bestModProb=modProb; bestTotalModProb=totalModProb; bestChi2=chi2;
	} else
	  ++plateau;
	if((useChi2 && chi2>oldChi2) || ( !useChi2 && (logDistP<oldLogDistP))) {
	  double rr=gsl_rng_uniform(rnd);
	  double pp=useChi2?pow(2.0,10*(oldChi2-chi2)):pow(2.0,1.0*(logDistP-oldLogDistP));
	  // mcmc rejection condition.
	  // temporary print out to make sure that it's running correctly
	  //std::cout << rr << '\t' << pp << '\t' << (logDistP-oldLogDistP)  << '\t' << logDistP << '\t' << oldLogDistP;
	  if(rr>pp) {
	    //std::cout << "\treject\n";
	    logDistP=oldLogDistP; 
	    logModel=oldLogModel;
	    chi2=oldChi2;
	  } //else
	    //std::cout << "\taccept\n";
	}
	++iter;
      }
      std::cout << "# Start " << s << ", bestLogDistP= " << std::setprecision(4) << bestLogDistP << ", bestX2= "
		<< bestChi2 << ", it= " << iter << ", motifSize= " << bestLogModel.size();
      if(!s || (useChi2 && (bestChi2<minChi2)) || (!useChi2 && (bestLogDistP>maxMultiProb))) {
	minChi2= bestChi2;
	bestMotifScore[i]= bestLogDistP;
	maxMultiProb=bestLogDistP;
	maxLogModel=bestLogModel;
	maxModProb= bestModProb;
	maxTotalModProb= bestTotalModProb;
	std::cout << " ** top **\n";
	for(int jj=0;jj<75;++jj) std::cout << "#";
	std::cout << " working Model ";
	for(int jj=0;jj<75;++jj) std::cout << "#";
	std::cout << '\n';
	for(int jj=0;jj<maxLogModel.size();++jj) {
	  std::cout << "# " << std::setw(3) << jj << ":";
	  dVec pp(4,0.0);
	  for(int kk=0;kk<4;++kk) pp[kk]= pow(2.0,maxLogModel[jj][kk]);
	  double temp=pp[0]; pp[0]=pp[1]; pp[1]=pp[3]; pp[3]=pp[2]; pp[2]=temp;
	  for(int kk=0;kk<4;++kk) std::cout << std::setw(8) << std::setprecision(2) << pp[kk]; 
	  std::cout << '\t' << SeqLogoString()(pp,100) << '\n';
	}
	for(int jj=0;jj<165;++jj) std::cout << "#"; 
      }
      std::cout << std::endl;
    } // end of nStarts loop
    for(int jj=0;jj<165;++jj) std::cout << "#";
    std::cout << "\n";
    for(int jj=0;jj<76;++jj) std::cout << "#";
    std::cout << " Best  Model ";
    for(int jj=0;jj<76;++jj) std::cout << "#";
    std::cout << '\n';
    std::cout << "# best model, logOdds= " << std::setprecision(4) << maxMultiProb;
    std::cout << ", Chi-Squared= " << minChi2 << "\n#             A         C         G         T   IC + Logo \n";
    mFile << "model " << (i+1) << ", multinomial log-odds= " << maxMultiProb << ", chi2= " << minChi2 << '\n';
    for(int j=0;j<maxLogModel.size();++j) {
      std::cout << "# " << std::setw(2) << j << ":";
      mFile << j;
      double IC=2.0; dVec pp(4,0.0);
      for(int k=0;k<4;++k) {
	double p= pow(2.0,maxLogModel[j][k]); 
	std::cout << std::setw(10) << p;
	mFile << '\t' << std::setprecision(9) << p;
	IC += p*maxLogModel[j][k]; pp[k]=p;
      }
      //      std::cout << std::setw(10) << IC << '\n';
      dVec logoP(4,0.0); logoP[0]=pp[1]; logoP[1]=pp[3]; logoP[2]=pp[0]; logoP[3]=pp[2];
      std::string logo= SeqLogoString()(logoP,100);
      std::cout << "   " << logo << '\n';
      mFile << '\t' << IC;
      //for(int k=0;k<4;++k) mFile << '\t' << (pp[k]*IC);
      mFile << '\t' << logo << '\n';
    }
    for(int jj=0;jj<165;++jj) std::cout << "#";
    std::cout << "\n";
    for(int jj=0;jj<165;++jj) std::cout << "#";
    std::cout << "\n\n";
    // dump probs and weights to a file for checking
    cFile << "Model " << i << ", " << wwm[i].size() << " words included\n";
    for(int j=0;j<wwm[i].size();++j) {
      double f1= pow(2.0,weightMat[i][j]),f2=pow(2.0,maxModProb[j]-maxTotalModProb);
      cFile << wwm[i][j].word << '\t' << f1 << '\t' << f2 << '\t' << wwm[i][j].weight;
      cFile << '\t' << pow(2.0,bgProb[j]) << '\t' << pow(2.0,maxModProb[j]) << '\t' << f2;
      cFile << '\t' << fabs(f1-f2) << '\t' << ((f1-f2)/f1) << '\n';
      //      cFile << '\t' << fabs(pow(2.0,weightMat[i][j])-pow(2.0,maxModProb[j]-maxTotalModProb)/pow(2.0,weightMat[i][j])) << '\n';
    }
  }
  cFile.close();
  mFile << "Model scores:\n";
  std::cout << "Model scores:\n";
  for(int i=0;i<r;++i) {
    mFile << std::setprecision(4) << bestMotifScore[i] << '\n';
    std::cout << std::setprecision(4) << " " << bestMotifScore[i];
  }
  mFile.close();
  std::cout << '\n';
  return 0;
}
