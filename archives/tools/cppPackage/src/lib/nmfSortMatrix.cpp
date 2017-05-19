//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include "../lib/logFile.h"
#include "../lib/util.h"
#include "../lib/arg.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <functional>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef std::vector<double> dVec;
typedef std::vector<int> iVec;
typedef std::vector<dVec> dMat;
typedef std::vector<std::string> strVec;
typedef std::map<std::string,dVec,std::less<std::string> > strDvecMap;
typedef std::map<std::string,double,std::less<std::string> > strDblMap;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool noLogs=false,quiet=false,mask=false;//,absolute=false; 
//double threshold=0.5; //int zPos=0;
std::string ofPrefix("nmfSorted");
double maskTh=2.0,wordTh=0.0;
dVec threshold; 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  if(a.KeyPresent(std::string("nologs"))) { noLogs=true; LogsOff(); ++pCount; }
  //  if(a.KeyPresent(std::string("abs"))) { absolute=true; ++pCount; }
  //if(a.KeyPresent(std::string("mask"))) { mask=true; ++pCount; }
  if(a.KeyPresent(std::string("quiet")) || a.KeyPresent(std::string("Q"))) { quiet=true; ++pCount; }
  //if(a.GetFloatValue(std::string("m"),std::string("mask threshold"),maskTh,-4.0,4.0,true,true)) ++pCount;
  //if(a.GetFloatValue(std::string("w"),std::string("word threshold"),wordTh,-4.0,4.0,true,true)) ++pCount;
  //if(a.GetFloatVec(std::string("t"),std::string("score threshold vector"),threshold)) ++pCount;
  if(a.GetString(std::string("o"),std::string("output prefix"),ofPrefix)) ++pCount;
  //if(a.GetString(std::string("B"),std::string("bgFilename"),bgFilename)) ++pCount;
  if(pCount!=a.numOptArgs()) {
    WriteToLogs()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLogs()(std::string("entered= "),a.numOptArgs());
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int LoadLabeledMatrix(const std::string& fname,dMat& mat,strVec& labels) {
  std::ifstream infile(fname.c_str(),std::ios::in);
  char line[1001];
  infile.getline(line,1000);
  strVec f= split()(std::string(line),'\t');
  int nCols= f.size()-1;
  //  labels.clear(); mat.clear();
  if(nCols<1) return 0;
  labels.push_back(f[0]);
  dVec temp(nCols,0.0);
  for(int i=0;i<nCols;++i) temp[i]=atof(f[i+1].c_str());
  dVec sum=temp;
  mat.push_back(temp);
  while(!infile.eof()) {
    infile.getline(line,1000);
    std::string lineIn(line);
    f= split()(lineIn);
    if(f.size()==nCols+1) {
      labels.push_back(f[0]);
      for(int i=0;i<nCols;++i) {
	temp[i]=atof(f[i+1].c_str());
	sum[i] += temp[i];
      }
      mat.push_back(temp);
    }
  }
  infile.close();
  //for(int i=0;i<nCols;++i)
  //  for(int j=0;j<mat.size();++j)
  //    mat[j][i] = log2(mat[j][i]/sum[i]);
  return nCols;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct idxMode{
  int idx,mode;
  idxMode(int i,int m):idx(i),mode(m) {}
  bool operator < (const idxMode& o) const { return mode < o.mode; }
};
typedef std::vector<idxMode> imVec;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct rnkWord{
  double score;
  std::string word;
  rnkWord(std::string w,double s):word(w),score(s) {}
  bool operator < (const rnkWord& o) const { return score > o.score; }
};
typedef std::vector<rnkWord> rwVec;
typedef std::vector<rwVec>   rwMat;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc,char *argv[]) {
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<2)||(arg1==std::string("-h"))||(arg1==std::string("-help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]\n";
    std::cout << "       " << argv[0] << " fileprefix [-o filestring] [-s float] [-quiet] [-nologs]\n"; 
    exit(0); 
  }
  if(!LogSetup()(argc,argv)) { std::cout << "error creating logs. exiting.\n"; exit(-1); }
  cjArgs arguments(argc,argv,1);
  ProcessArgs(arguments);
  std::string bfname(arguments.RequiredArgs()[0]+std::string(".bases.txt"));
  std::string wfname(arguments.RequiredArgs()[0]+std::string(".weights.txt"));
  //
  std::cout << "Loading base file: " << bfname << "..."; std::cout.flush();
  dMat bases; bases.reserve(10); strVec words; words.reserve(500);
  int nBases= LoadLabeledMatrix(bfname,bases,words);
  //
  std::cout << "done.\nLoading weight file: " << wfname << "..."; std::cout.flush();
  dMat weights; weights.reserve(100); strVec posLabels; iVec pos; posLabels.reserve(500); pos.reserve(500);
  int nModels= LoadLabeledMatrix(wfname,weights,posLabels);
  //
  std::cout << "done.\n\nWriting summary of models..."; std::cout.flush();
  std::cout << "# nModels= " << nModels << "\n# nBases=" << nBases << '\n';
  std::cout << "# WeightLength= " << weights.size() << "\n# BaseLength= " << bases.size() << '\n';
  std::cout << "# word count= " << words.size() << "\n# pos count= " << posLabels.size() << '\n';
  int W= words[0].length(),step= atof(posLabels[2].c_str())-atof(posLabels[1].c_str());
  std::cout << "# Weight step= " << step << ", Word length= " << W << '\n';
  std::cout << "Searching for distribution modes..."; std::cout.flush();
  iVec modePos(nModels,0); imVec modes;
  for(int i=0;i<nModels;++i) {
    for(int j=1;j<weights.size();++j)
      if(weights[j][i]>weights[modePos[i]][i]) 
	modePos[i]= j;
    std::cout << std::setw(4) << modePos[i];
    modes.push_back(idxMode(i,modePos[i]));
  }
  std::sort(modes.begin(),modes.end());
  //
  std::cout << "\nRe-sorted elements:\n";
  for(int i=0;i<modes.size();++i)
    std::cout << (modes[i].idx+1) << " => " << (i+1) << " [mode= " << modes[i].mode << "]\n";
  //
  std::string wfname2(ofPrefix+std::string(".weights.txt"));
  std::cout << "\n\nWriting sorted weight file " << wfname2 << "..."; std::cout.flush();
  std::ofstream wfile2(wfname2.c_str(),std::ios::out);
  //std::cout << "\n\n#weights\n";
  dVec colSum(nModels,0.0);
  for(int i=0;i<weights.size();++i) {
    wfile2 << posLabels[i];
    for(int j=0;j<nModels;++j) {
      wfile2 << '\t' << weights[i][modes[j].idx];
      colSum[j] += weights[i][modes[j].idx];
    }
    wfile2 << '\n';
  }
  wfile2.close();
  // now write the normalized version
  std::string wfname3(ofPrefix+std::string(".nweights.txt"));
  std::cout << "done\nWriting sorted normalized weight file " << wfname3 << "..."; std::cout.flush();
  std::ofstream wfile3(wfname3.c_str(),std::ios::out);
  //std::cout << "\n\n#weights\n";
  for(int i=0;i<weights.size();++i) {
    wfile3 << posLabels[i];
    for(int j=0;j<nModels;++j) 
      wfile3 << '\t' << weights[i][modes[j].idx]/colSum[j];
    wfile3 << '\n';
  }
  wfile3.close();
  //
  std::string bfname2(ofPrefix+std::string(".bases.txt"));
  std::cout << "done\nWriting sorted base file " << bfname2 << "..."; std::cout.flush();
  std::ofstream bfile2(bfname2.c_str(),std::ios::out);
  for(int i=0;i<bases.size();++i) {
    bfile2 << words[i];
    for(int j=0;j<nModels;++j)
      bfile2 << '\t' << bases[i][modes[j].idx];
    bfile2 << '\n';
  }
  bfile2.close();
  std::cout << "done.\nGenerating Scored Word lists..."; std::cout.flush();
  rwMat rankedWords(nModels,rwVec());
  for(int i=0;i<nModels;++i) {
    rankedWords[i].reserve(words.size());
    for(int j=0;j<words.size();++j)
      rankedWords[i].push_back(rnkWord(words[j],bases[j][modes[i].idx]));
    std::sort(rankedWords[i].begin(),rankedWords[i].end());
  }
  std::string rfname(ofPrefix+std::string(".rwords.txt"));
  std::ofstream rfile(rfname.c_str(),std::ios::out);
  for(int i=0;i<words.size();++i) 
    for(int j=0;j<nModels;++j)
      if(j != (nModels-1)) rfile << rankedWords[j][i].word << ": " << std::setprecision(3) << rankedWords[j][i].score << '\t'; 
      else                 rfile << rankedWords[j][i].word << ": " << std::setprecision(3) << rankedWords[j][i].score << '\n'; 
  rfile.close();
  std::cout << "done.\n";
  //
  //std::string ofname(ofPrefix+std::string(".sites.txt"));
  //std::ofstream ofile(ofname.c_str(),std::ios::out);
  //ofile.close();
  //dfile.close();
  //
  return 0;
}
