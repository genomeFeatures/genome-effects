///////////////////////////////////////////////////////////////////////////////////////////////////////
#include <math.h>
#ifdef INTEL
#include <mathimf.h>
#endif
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "../lib/basicSeq.h"
#include "../lib/arg.h"
#include "../lib/PosWord.h"
#include "../lib/fastaSeq.h"
#include "../lib/logFile.h"
#include <vector>
#include <algorithm>
///////////////////////////////////////////////////////////////////////////////////////////////////////
typedef std::vector<strVec> strMat;
typedef std::vector<double> dVec;
typedef std::vector<dVec> dMat;
///////////////////////////////////////////////////////////////////////////////////////////////////////
bool ds=false,useZS=false,smooth=false,S2only=false;
std::string ofName("zStats"),testSeq("AATAAA");
int Order=0,WindL=10,WordL=6,SeqL=300,zPnt=0,minWordC=0;
double pcWeight=0.01,smSig=5.0;//minSSqr=50.0,smSig=5.0;
///////////////////////////////////////////////////////////////////////////////////////////////////////
struct WordX2 { 
  std::string word; double X2; dVec V; 
  WordX2(std::string w, double x, dVec v):word(w),X2(x),V(v) {}
  WordX2(const WordX2& o):word(o.word),X2(o.X2),V(o.V) {}
  WordX2():word(),X2(),V() {}
  bool operator < (const WordX2& o) const { return (X2>o.X2 || (X2==o.X2 && word<o.word)); }
};
typedef std::vector<WordX2> WordX2Vec;
///////////////////////////////////////////////////////////////////////////////////////////////////////
struct WordRankVec { 
  std::string word; iVec rank; dVec s;
  WordRankVec(std::string w, int n):word(w),rank(n,0),s(n,0.0) {}
  WordRankVec():word(),rank(),s() {}
  int rankSum() const { 
    int s=0; for(iVec::const_iterator i=rank.begin();i!=rank.end();++i) s+= (*i); return s; }
  bool operator < (const WordRankVec& o) const { 
    int s=rankSum(), os=o.rankSum(); return (s<os || (s==os && word<o.word)); }
};
typedef std::vector<WordRankVec> WordRvVec;
///////////////////////////////////////////////////////////////////////////////////////////////////////
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  if(a.KeyPresent(std::string("d"))) { ds=true; ++pCount; }
  if(a.KeyPresent(std::string("zs"))) { useZS=true; ++pCount; }
  if(a.KeyPresent(std::string("sm"))) { smooth=true; ++pCount; }
  if(a.KeyPresent(std::string("s2"))) { S2only=true; ++pCount; }
  if(a.GetIntValue(string("k"),string("word length"),WordL,0,10,true,true)) ++pCount;
  if(a.GetIntValue(string("r"),string("order"),Order,-1,6,true,true)) ++pCount;
  if(a.GetIntValue(string("L"),string("sequence profile length"),SeqL,WordL,5000,true,true)) ++pCount;
  if(a.GetIntValue(string("z"),string("zero point"),zPnt,-10000,10000)) ++pCount;
  if(a.GetIntValue(string("w"),string("window length"),WindL,0,100)) ++pCount;
  if(a.GetIntValue(string("M"),string("min word count"),minWordC,0,10000)) ++pCount;
  if(a.GetFloatValue(string("p"),string("pseudo count weight"),pcWeight,0.0,1.0)) ++pCount;
  //  if(a.GetFloatValue(string("S"),string("min s-Squared"),minSSqr,0.0,1000.0)) ++pCount;
  if(a.GetFloatValue(string("sig"),string("smooth sigma"),smSig,0.0,50.0)) ++pCount;
  if(a.GetString(std::string("o"),std::string("outfile prefix"),ofName)) ++pCount; 
  if(a.GetString(std::string("ts"),std::string("test sequence"),testSeq)) ++pCount; 
  if(pCount!=a.numOptArgs()) {
    WriteToLogs()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLog()(std::string("processed= "),pCount);
    WriteIntToLog()(std::string("entered= "),a.numOptArgs());
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void PrintParams(std::ostream& os,int argc, char* argv[]) {
  time_t t0; time(&t0);
  os << "################################################################################\n";
  os << "#                Window Count: Nucleic Acid Profiled Counts                    #\n"; 
  os << "################################################################################\n";
  os << "# " << ctime(&t0);
  os << "# cmdLine: " << argv[0]; for(int i=1;i<argc;++i) os << " " << argv[i]; os << '\n';
  os << "################################################################################\n";
  os << "# Input Fasta File                        = " << argv[1] << '\n';
  os << "############################# Parameter settings: ##############################\n";
  os << "# Word Length                         [k] = " << WordL << '\n';
  os << "# Max Markov Chain Order              [r] = " << Order << '\n';
  os << "# Sequence Profile Length             [L] = " << SeqL << '\n';
  os << "# Sequence Zero Length                [z] = " << zPnt << '\n';
  os << "# Counting Window Size                [w] = " << WindL << '\n';
  os << "# Pseudo-count weight                 [p] = " << pcWeight << '\n';
  os << "# Output File Prefix                  [o] = " << ofName << '\n';
  os << "# Double Strand Counting              [d] = ";
  if(ds) os << "On\n"; else os << "Off\n";
  os << "# Window Smoothing                   [sm] = ";
  if(smooth) 
    os << "On\n# Smoothing sigma (gaussian window) [sig] = " << smSig << '\n';
  else
    os << "Off\n";
  os << "################################################################################\n\n";
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
typedef std::vector<strVec> strMat;
///////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<2)||(arg1==std::string("-h"))||(arg1==std::string("-help"))||(arg1==std::string("--help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]\n";
    std::cout << "       " << argv[0] << " filename [-d] [-o outFile] [-r order] [-w window size]\n";
    std::cout << "                         [-z zero Point] [-k # of clusters] [-L profile length]\n"; 
    exit(0); 
  }
//
  if(!LogSetup()(argc,argv)) { std::cout << "error opening logs" << std::endl; exit(-1);}
  cjArgs arguments(argc,argv,1);
  ProcessArgs(arguments);
  std::string fname= arguments.RequiredArgs()[0];
  PrintParams(std::cout,argc,argv);
//
  TJfastaSequences seqs;
  std::cout << "Loading sequences from: " << fname.c_str() << "..."; std::cout.flush();
  seqs.LoadSequences(fname);
  std::cout << "done.\nAnalyzing " << seqs.size() << " sequences for word counts..."; std::cout.flush();
  strVec rawSeqs; rawSeqs.reserve(seqs.size());
  for(int i=0;i<seqs.size();++i) rawSeqs.push_back(seqs[i].Sequence());
  WordPositionCounts wpc(WordL,SeqL,ds,rawSeqs); 
  std::cout << "done.\nBuilding windowed count..."; std::cout.flush();
  WordWindowFullCount wwc(wpc,WordL,WindL,-1*zPnt,0,smooth,smSig);
  std::cout << "done.\n"; std::cout.flush();
  int nWords= pow(4,WordL); 
  WordRvVec wordRanks; wordRanks.reserve(nWords);
  if(!S2only)
      for(int i=0;i<nWords;++i) wordRanks.push_back(WordRankVec(IdxToSeq()(i,WordL),Order+2));
  else
      for(int i=0;i<nWords;++i) wordRanks.push_back(WordRankVec(IdxToSeq()(i,WordL),1));
  iVec bounds= wwc.Boundaries();
  //
  int update= (nWords<50)?1:nWords/50; dMat perc;
  int loopStop= S2only?1:Order+2;
  for(int o=0;o<loopStop;++o) {
    if(!o) std::cout << "Calculating S-squared values";  
    else   std::cout << "Calculating X-squared values, order= " << (o-1);
    std::cout.flush();
    WordX2Vec words; words.reserve(nWords); 
    for(int i=0;i<nWords;++i) {
      if(!((i+1)%update)) { std::cout << "."; std::cout.flush(); }
      std::string w= IdxToSeq()(i,WordL);
      //    double x2= wwc.SSqr(w);
      double x2; dVec v;
      if(o) { 
	if(!useZS) x2= wwc.ChiSqr(w,o-1,v);
	else       x2= wwc.ZSSqr(w,o-1,v);
      } else  {
	x2= wwc.SSqr(w);
	v= wwc.CountVector(w);
      }
      if(wwc.Count(w)>=minWordC) words.push_back(WordX2(w,x2,v));
    }
    std::cout << "done. Sorting..."; std::cout.flush();
    std::sort(words.begin(),words.end());
    // collect the percentile
    perc.push_back(dVec());
    perc[o].push_back(words[nWords/1000].X2); perc[o].push_back(words[nWords/100].X2);
    perc[o].push_back(words[nWords/10].X2);   perc[o].push_back(words[nWords/4].X2);
    perc[o].push_back(words[nWords/2].X2);    perc[o].push_back(words[nWords*3/4].X2);
    perc[o].push_back(words[nWords*9/10].X2); perc[o].push_back(words[nWords*99/100].X2);
    perc[o].push_back(words[nWords*999/1000].X2);
    //
    std::cout << "done.\nTop 10 words based on x2" << std::endl;
    for(int i=0;i<10 && i<words.size();++i)
      std::cout << words[i].word << '\t' << words[i].X2 << '\n';
    std::cout << "\nCollecting ranks and writing output file..."; std::cout.flush();
    std::string cfname(ofName),sfname(ofName);
    if(!o) {
      cfname += std::string(".s2.counts.txt");
      sfname += std::string(".s2.txt");
    } else {
      char cOrder= '0'+static_cast<char>(o-1); 
      std::string co; co += cOrder;
      cfname += std::string(".x2_") + co + std::string(".counts.txt");
      sfname += std::string(".x2_") + co + std::string(".txt");
    }
    std::ofstream cFile(cfname.c_str(),std::ios::out);
    std::ofstream sFile(sfname.c_str(),std::ios::out);
    cFile << "Word"; 
    if(o) sFile << "Word\tX^2"; else sFile << "Word\tS^2\n";
    for(int i=0;i<bounds.size()-1;++i) cFile << '\t' << bounds[i]; 
    cFile << '\n';
    for(int i=0;i<words.size();++i) {
      wordRanks[SeqToIdx()(words[i].word)].rank[o]= i;
      wordRanks[SeqToIdx()(words[i].word)].s[o]= words[i].X2;
      sFile << words[i].word << '\t' << words[i].X2 << '\n';
      cFile << words[i].word;
      for(int j=0;j<words[i].V.size();++j) cFile << '\t' << words[i].V[j];
      cFile << '\n';
    }
    cFile.close();
    std::cout << "done." << std::endl;  
  }
  std::cout << "Sorting based on RankSum..."; std::cout.flush();
  std::sort (wordRanks.begin(),wordRanks.end());
  std::cout << "done.\nTop 20:\nSeq\tsum\tS2";
  if(!S2only)  for(int j=1;j<wordRanks[0].rank.size();++j) std::cout << "\tX2(" << j-1 << ")";
  std::cout << "\n";
  std::string rfname(ofName+std::string(".ranks.txt")),chfname(ofName+std::string(".chi.txt"));
  std::ofstream RankFile(rfname.c_str(),std::ios::out),ChiFile(chfname.c_str(),std::ios::out);
  RankFile << "Rank\tWord\tsum\tS2"; ChiFile << "Word\tS2";
  for(int j=1;j<wordRanks[0].rank.size();++j) {
    RankFile << "\tX2(" << (j-1) << ")";
    ChiFile << "\tX2(" << (j-1) << ")";
  }
  RankFile << "\n"; ChiFile << "\n";
  for(int i=0;i<wordRanks.size();++i) {
    if(i<20) {
      std::cout << wordRanks[i].word << '\t' << wordRanks[i].rankSum();
      for(int j=0;j<wordRanks[i].rank.size();++j) std::cout << '\t' << wordRanks[i].rank[j];
      std::cout << '\n';
    }
    RankFile << i << '\t' << wordRanks[i].word << '\t' << wordRanks[i].rankSum();
    for(int j=0;j<wordRanks[i].rank.size();++j)
      RankFile << '\t' << wordRanks[i].rank[j];
    RankFile << '\n';
    ChiFile << wordRanks[i].word;
    for(int j=0;j<wordRanks[i].s.size();++j)
      ChiFile << '\t' << wordRanks[i].s[j];
    ChiFile << '\n';
  }
  RankFile.close(); ChiFile.close();
  std::cout << "percentile comparisons\n\t.001\t.01\t.1\t.25\t.5\t.75\t.9\t.99\t.999\n";
  for(int i=0;i<perc.size();++i) {
    std::cout << i;
    for(int j=0;j<perc[i].size();++j) std::cout << '\t' << perc[i][j];
    std::cout << '\n';
  }
  WriteToLogTS()(std::string("Normal termination."));
  WriteToErrTS()(std::string("Normal termination."));
  return 0;
}
