///////////////////////////////////////////////////////////////////////
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib/arg.h"
#include "lib/logFile.h"
#include "lib/util.h"
#include "lib/vectorF.h"
#include "nnmf.h"
#include <string>
#include <vector>
///////////////////////////////////////////////////////////////////////
typedef std::vector<double> dVec;
typedef std::vector<dVec> dMat;
typedef std::vector<std::string> strVec;
///////////////////////////////////////////////////////////////////////
bool quiet=false,noLogs=false,hnorm=false,vnorm=false,transpose=true;
string ofPrefix("nnmf"); double tol=0.001;
int nBases=5,maxIt=1000,nStarts=20;//NRand=0,LRand=0,ORand=0,LineL=80,nBins=51;
///////////////////////////////////////////////////////////////////////
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  if(a.KeyPresent(std::string("nologs"))) { noLogs=true; LogsOff(); ++pCount; }
  if(a.KeyPresent(std::string("hnorm"))) { hnorm=true; ++pCount; }
  if(a.KeyPresent(std::string("vnorm"))) { vnorm=true; ++pCount; }
  if(a.KeyPresent(std::string("trans"))) { transpose=false; ++pCount; }
  if(a.KeyPresent(std::string("quiet")) || a.KeyPresent(std::string("Q"))) { quiet=true; ++pCount; }
  if(a.GetIntValue(std::string("r"),std::string("n bases"),nBases,2,25,true,true)) ++pCount;
  if(a.GetIntValue(std::string("I"),std::string("max It"),maxIt,10,100000,true,true)) ++pCount;
  if(a.GetIntValue(std::string("S"),std::string("n starts"),nStarts,0,10000,true,true)) ++pCount;
  if(a.GetFloatValue(std::string("t"),std::string("tolerance"),tol,0.0,100.0,true,true)) ++pCount;
  if(a.GetString(std::string("o"),std::string("output prefix"),ofPrefix)) ++pCount;
  if(pCount!=a.numOptArgs()) {
    WriteToLogs()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLogs()(std::string("entered= "),a.numOptArgs());
  }
  return;
}
///////////////////////////////////////////////////////////////////
typedef struct wordWeightRec {
  std::string word;  double weight;
  wordWeightRec(const std::string& wo,double we):word(wo),weight(we) {}
};
bool operator < (const wordWeightRec& r1, const wordWeightRec& r2) { return r1.weight > r2.weight; }
typedef std::vector<wordWeightRec> wwVec;
///////////////////////////////////////////////////////////////////
dMat MakeMotif(const wwVec& wv,double base,std::ostream& os) {
  if(wv.empty() || base<0.0) return dMat();
  int k=wv[0].word.length(); 
  int km1=k-1;
  dMat motifScore(3*k-2,dVec(5,base));
  for(int i=0;i<motifScore.size();++i) motifScore[i][4]*=4.0;
  //
  for(int w=0;w<wv.size();++w) {
    double maxscore=0.0; int maxo=3;
    if(w) 
      for(int o= -km1;o<=km1;++o) {
	double score=0.0;      
	for(int p=0;p<k;++p) {
	  if(wv[w].word[p]=='A') score += motifScore[km1+o+p][0];
	  if(wv[w].word[p]=='C') score += motifScore[km1+o+p][1];
	  if(wv[w].word[p]=='G') score += motifScore[km1+o+p][2];
	  if(wv[w].word[p]=='T') score += motifScore[km1+o+p][3];
	}
	//os << "  update: " << wv[w].word << "\t" << o << '\t' << score << '\n';
	if(score>maxscore || (score==maxscore && abs(o)<abs(maxo))) {
	  maxscore=score;
	  maxo=o;
	}
      }
    else
      maxo=0;
    os << w << ":" << wv[w].word << "\t" << maxo << '\t' << maxscore << '\n';
    if(maxscore>(4.0*base) || !w) {
      for(int p=0;p<k;++p) {
	double factor= pow(4.0,-1.0*static_cast<double>(abs(maxo)));
	if(maxo<0 && (p+maxo)<0) factor *= pow(4.0,static_cast<double>(-p-maxo));
	else if(maxo>0 && (p+maxo-km1)>0) factor *= pow(4.0,static_cast<double>(p+maxo-km1));
	motifScore[km1+maxo+p][4] += wv[w].weight * factor;
	if(wv[w].word[p]=='A') motifScore[km1+maxo+p][0] += wv[w].weight * factor;
	else if(wv[w].word[p]=='C') motifScore[km1+maxo+p][1] += wv[w].weight * factor;
	else if(wv[w].word[p]=='G') motifScore[km1+maxo+p][2] += wv[w].weight * factor;
	else if(wv[w].word[p]=='T') motifScore[km1+maxo+p][3] += wv[w].weight * factor;
      }
    }
  }
  return motifScore;
}
///////////////////////////////////////////////////////////////////
void DumpParameters(std::ostream& os,int argc,char *argv[]) {
  time_t t0; time(&t0);
  os << "############################ " << argv[0] << "Execution parameters ############################\n";
  os << "# execution time: " << ctime(&t0);
  os << "# cmd line:"; for(int i=0;i<argc;++i) os << ' ' << argv[i]; os << "\n"; 
  os << "# input file                     = " << argv[1] << "\n";
  os << "# output prefix                  = " << ofPrefix << "\n";
  os << "# n Bases    [-r]                = " << nBases << '\n';
  os << "# tolerance  [-t]                = " << tol << '\n';
  os << "# max Iter   [-I]                = " << maxIt << '\n';
  os << "# n Starts   [-S]                = " << nStarts << '\n';
  os << "# quiet                          = " << (quiet?"on\n":"off\n");
  os << "# nologs                         = " << (noLogs?"on\n":"off\n");
  os << "# normalization [-hnorm|-vnorm]  = ";
  if(hnorm)       os << "horizontal\n";
  else if (vnorm) os << "vertical\n";
  else            os << "none\n";
  os << "# transpose [-trans]             = " << (transpose?"on\n":"off\n");
  os << "#####################################################################################\n";
  return;
}
///////////////////////////////////////////////////////////////////
void DumpMatrixes(const dMat& b,const dMat& w,const strVec& windows,const strVec& words) {
  std::cout << "bases\n";
  std::string bfilename(ofPrefix+std::string(".bases.txt"));
  std::ofstream bfile(bfilename.c_str(),std::ios::out);
  for(int i=0;i<b.size();++i) {
    if(!transpose) std::cout << windows[i];
    else           std::cout << words[i];
    if(!transpose) bfile << windows[i];
    else           bfile << words[i];
    for(int j=0;j<b[i].size();++j) {
      std::cout << "\t" << std::setprecision(6) << 100.0*b[i][j];
      bfile << "\t" << std::setprecision(6) << 100.0*b[i][j];
    }
    std::cout << '\n'; bfile << '\n';
  }
  bfile.close();
  std::cout << "weights(transposed)\n";
  std::string wfilename(ofPrefix+std::string(".weights.txt"));
  std::ofstream wfile(wfilename.c_str(),std::ios::out);
  for(int i=0;i<w[0].size();++i) {
    if(!transpose) std::cout << words[i];
    else           std::cout << windows[i];
    if(!transpose) wfile << words[i];
    else           wfile << windows[i];
    for(int j=0;j<w.size();++j) {
      wfile << "\t" << std::setprecision(6) << w[j][i];
      std::cout << "\t" << std::setprecision(6) << w[j][i];
    }
    std::cout << '\n'; wfile << '\n';
  }
  wfile.close();
  return;
}
///////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<2)||(arg1==std::string("-h"))||(arg1==std::string("-help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]\n";
    std::cout << "       " << argv[0] << " filename [-o filestring] [-r int] [-l int] [-n int] [-d] [-quiet] [-nologs]\n"; 
    std::cout <<"\nWhere: filename --> Input file matrix\n";
    std::cout << "       " <<"-o       --> Output prefix\n"; 
    std::cout << "       " <<"-r       -->  n Bases\n"; 
    std::cout << "       " <<"-l       -->\n"; 
    std::cout << "       " <<"-n       -->\n"; 
    std::cout << "       " <<"-d       -->\n"; 
    std::cout << "       " <<"-quiet   --> Not Verbose\n"; 
    std::cout << "       " <<"-nologs  --> Do not generate log files\n\n"; 
    exit(0); 
  }
  if(!LogSetup()(argc,argv)) { std::cout << "error creating logs. exiting.\n"; exit(-1); }
//
//
  cjArgs arguments(argc,argv,1);
  ProcessArgs(arguments);
  std::string fname(arguments.RequiredArgs()[0]);
  DumpParameters(std::cout,argc,argv);
  // Load the file
  std::ifstream inFile(fname.c_str(),std::ios::in);
  char lineIn[10001];
  inFile.getline(lineIn,10000);
  strVec f= split()(std::string(lineIn));
  int M= f.size()-1;
  strVec windows; windows.reserve(M); 
  for(int i=1;i<=M;++i) windows.push_back(f[i]);
  inFile.getline(lineIn,10000);
  f= split()(std::string(lineIn));
  int K= f[0].length();
  int N=int(pow(4.0,K));
  strVec words; words.reserve(N);
  dMat counts; counts.reserve(N); 
  dVec totalCounts; totalCounts.reserve(N);
  int i=0;
  std::cout << "K= " << K << "; M= " << M << "; N= " << N;
  while(!inFile.eof() && i<N) {
    if(i) {
      inFile.getline(lineIn,10000);
      f= split()(std::string(lineIn));
    }
    words.push_back(f[0]);
    dVec temp; temp.reserve(500);
    double total=0.0;
    for(int j=1;j<f.size();++j) {
      temp.push_back(atof(f[j].c_str()));
      total += temp.back();
    }
    totalCounts.push_back(total);
    if(temp.size()==M) counts.push_back(temp);
    ++i;
  }
  std::cout << ", count= " << counts.size() << std::endl; 
  inFile.close();
  // vertical normalization: normalize to the total counts in the column (position)
  if(vnorm) {
    double maxsum=0.0;
    for(int j=0;j<counts[0].size();++j) {
      double sum=0.0;
      for(int i=0;i<counts.size();++i) sum += counts[i][j];
      if(sum>maxsum) maxsum=sum;
    }
    for(int j=0;j<counts[0].size();++j)
      for(int i=0;i<counts.size();++i) 
	counts[i][j] /= maxsum;
  }

  // horizontal normalization: normalize all words to constant counts
  else if(hnorm)
    for(int i=0;i<counts.size();++i) {
      double sum=0.0;
      for(int j=0;j<counts[i].size();++j) sum += counts[i][j];
      for(int j=0;j<counts[i].size();++j) counts[i][j] /= sum;
    }
  // by default we transpose
  dMat cc(counts);
  if(!transpose) {
    cc= dMat(M,dVec(counts.size(),0.0));
    for(int i=0;i<counts.size();++i)
      for(int j=0;j<counts[i].size();++j)
	cc[j][i]=counts[i][j];
  } 
  //
  std::cout << "####### Starting NNMF build number 1 ######"; 
  time_t lastLoop; time(&lastLoop);
  NNMF decomp(cc,nBases,tol,maxIt,quiet);
  dVec fScores; fScores.reserve(nStarts); fScores.push_back(decomp.Fscore());
  double score=decomp.Fscore(),res=decomp.Residue(),totalTime=0;
  double bestScore=score, bestRes=res;
  int totalIt= decomp.Iterations(),iter= decomp.Iterations();
  std::string pfname(ofPrefix+std::string(".prog.txt"));
  std::ofstream pfile(pfname.c_str(),std::ios::out);
  DumpParameters(pfile,argc,argv);
  for(int i=0;i<nStarts;++i) {
    bool best= i==0;
    if(i) {
      std::cout << "\n####### Starting NNMF build number " << (i+1) << " ######" << std::endl; 
      NNMF tdecomp(cc,nBases,tol,maxIt,quiet);
      totalIt+= (iter=tdecomp.Iterations());
      score= tdecomp.Fscore();
      res= tdecomp.Residue();
      fScores.push_back(score);
      if(score>bestScore) {
	decomp=tdecomp;
	bestScore= score;
	bestRes= res;
	best=true;
      }
    } 
    if(best) {
      dMat b= decomp.bases(),w=decomp.weights();
      DumpMatrixes(b,w,windows,words);
    }
    pfile << "########  F = " << std::setprecision(8) << score << ", res= " << res << ", it= " << iter;
    std::cout << "########  F = " << std::setprecision(8) << score << ", res= " <<  res << ", it= " << iter;
    if(score==bestScore) { 
      std::cout << " * top * ";
      pfile << " *top * ";
    }
    std::cout << "  ########\n";
    pfile << "  ########\n";
    time_t loopEnd; time(&loopEnd); 
    int delta= loopEnd-lastLoop; totalTime += delta;
    double avgTime= totalTime/static_cast<double>(i+1);
    double avgIter= static_cast<double>(totalIt)/static_cast<double>(i+1);
    pfile << "########  Iterations= " << iter << ", avg iterations= " << std::setprecision(6) 
	  << avgIter << "  ########\n########  Loop time= " << delta << " s, avg loop time= " 
	  << std::setprecision(6) << avgTime << " s ########\n########  Projected remaining time = " 
	  << std::setprecision(6) << static_cast<double>(nStarts-i-1)*avgTime << " s  ########" << std::endl;
    if(!quiet) {
      std::cout << "########  Iterations= " << iter << ", avg iterations= " << std::setprecision(6) 
		<< avgIter << "  ########\n########  Loop time= " << delta << " s, avg loop time= " 
		<< std::setprecision(6) << avgTime << " s ########\n########  Projected remaining time = " 
		<< std::setprecision(6) << static_cast<double>(nStarts-i-1)*avgTime << " s  ########\n";;
    }
    lastLoop=loopEnd;
  }
  pfile << "\n################# Sampled F-scores (" << nStarts << " seeds.) ###############\n";
  std::sort(fScores.begin(),fScores.end());
  for(int i=0;i<nStarts;++i) {
    pfile << std::setw(12) << std::setprecision(8) << fScores[i];
    if ((i+1)==nStarts || !((i+1)%8)) pfile << "\n";
  }
  std::cout << "################# Final F-score = " << std::setprecision(10) << decomp.Fscore() << " ################\n";
  pfile << "################# Final F-score = " << std::setprecision(10) << decomp.Fscore() << " ################\n";
  std::cout << "################# Final Residue = " << std::setprecision(10) << decomp.Residue() << " ################\n";
  pfile << "################# Final Residue = " << std::setprecision(10) << decomp.Residue() << " ################\n";
  pfile.close();
  WriteToLogTS()(std::string("Normal termination."));
  WriteToErrTS()(std::string("Normal termination."));
  return 0;
}
