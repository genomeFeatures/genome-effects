///////////////////////////////////////////////////////////////////////
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "../lib/fastaSeq.h"
#include "../lib/logFile.h"
#include "../lib/basicSeq.h"
#include "../lib/wordCount.h"
#include "../lib/arg.h"
#include <vector>
#include <string>
///////////////////////////////////////////////////////////////////////
using std::string; using std::vector;
bool oFile=false,ds=false,listAll=false,noLogs=false,quiet=false;
std::string ofPrefix;
int ORand=0,nBins=51;
///////////////////////////////////////////////////////////////////////
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  if(a.KeyPresent(std::string("nologs"))) { noLogs=true; LogsOff(); ++pCount; }
  if(a.KeyPresent(std::string("d"))) { ds=true; ++pCount; }
  if(a.KeyPresent(std::string("quiet")) || a.KeyPresent(std::string("Q"))) { quiet=true; ++pCount; }
  if(a.GetIntValue(std::string("r"),std::string("order"),ORand,0,7,true,true)) ++pCount;
  if(a.GetIntValue(std::string("b"),std::string("bin count"),nBins,10,500,true,true)) ++pCount;
  if(a.GetString(std::string("o"),std::string("output prefix"),ofPrefix)) ++pCount;
  if(a.KeyPresent(string("a"))) { ++pCount; listAll=true; }
  if(pCount!=a.numOptArgs()) {
    WriteToErr()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteToLog()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLog()(std::string("processed= "),pCount);
    WriteIntToLog()(std::string("entered= "),a.numOptArgs());
  }
  if(ds) WriteToLog()(std::string("double-stranded statistics"));
  else   WriteToLog()(std::string("single-stranded statistics"));
  return;
}
///////////////////////////////////////////////////////////////////
void DumpParameters(std::ostream& os,int argc,char *argv[]) {
  time_t t0; time(&t0);
  os << "############################ " << argv[0] << "Execution parameters ############################\n";
  os << "# execution time: " << ctime(&t0);
  os << "# cmd line:"; for(int i=0;i<argc;++i) os << ' ' << argv[i]; os << "\n"; 
  os << "# input file        = " << argv[1] << "\n";
  os << "# output prefix     = " << ofPrefix << "\n";
  os << "# order      [-r]   = " << ORand << '\n';
  os << "# bin count  [-b]   = " << nBins << '\n';
  os << "# ds                = " << (ds?"on\n":"off\n");
  os << "# quiet             = " << (quiet?"on\n":"off\n");
  os << "# nologs            = " << (noLogs?"on\n":"off\n");
  os << "#####################################################################################\n";
  return;
}
///////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<2)||(arg1==std::string("-h"))||(arg1==std::string("-help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]\n";
    std::cout << "       " << argv[0] << " filename [-o filestring] [-r int] [-d] [-quiet] [-nologs]\n"; 
    exit(0); 
  }
  if(!LogSetup()(argc,argv)) { std::cout << "error creating logs. exiting.\n"; exit(-1); }
//
  cjArgs arguments(argc,argv,1);
  if(arguments.numOptArgs()) ProcessArgs(arguments);
  std::string fname= arguments.RequiredArgs()[0];
//
  if(!(nBins%2)) ++nBins;
  if(!quiet) DumpParameters(std::cout,argc,argv);
  std::ofstream outFile;
  if(ofPrefix.length()) { 
      std::string ofname(ofPrefix+std::string(".stat"));
      outFile.open(ofname.c_str(),std::ios::out);
      DumpParameters(outFile,argc,argv);
  }
//
  if(!quiet) { std::cout << "Loading sequences from " << fname << "..."; std::cout.flush(); }
  TJfastaSequences seqs;
  seqs.LoadSequences(fname);
  if(seqs.empty()) {
    WriteToLogs()(std::string("no sequences read in"));
    exit(-1);
  }
  unsigned nSeqs= seqs.size(),totalChar=seqs.TotalCharacters();
  unsigned longSeq= seqs.LongestSeq(),shortSeq= seqs.ShortestSeq();
  double avgSeq= seqs.AverageLength(), medSeq= seqs.MedianLength();
  unsigned binW= static_cast<unsigned>(medSeq)/((nBins-1)/2);
  if(binW<5) binW=5;
  uVec hist= seqs.LengthHist(nBins,binW);
  //
  if(!quiet) { 
    std::cout << "done. " << seqs.size() << " loaded.\n";
    std::cout << "\n# Sequence count  = " << nSeqs 
	      << "\n# Character count = " << totalChar
	      << "\n# Shortest seq    = " << shortSeq
	      << "\n# Longest seq     = " << longSeq 
	      << "\n# Average length  = " << avgSeq 
	      << "\n# Median  length  = " << medSeq 
	      << "\n# size histogram. bin size " << binW << " up to " << (nBins*binW) << "\n";
    for(unsigned k=0;k<nBins;++k)   std::cout << ((k+1)*binW) << '\t' << hist[k] << '\n';

    std::cout << "\nCollecting stats..." << std::endl; 
  }
  //
  cjWordCount stats(ORand+1);
  for(unsigned j=0;j<seqs.size();j++) {
    if(!((j+1)%25)) {
      std::cout << ".";
      if((j+1)%2500) std::cout.flush();
      else           std::cout << (j+1) << " / " << seqs.size() << std::endl;
    }
    stats.AddSequence(seqs[j].Sequence());
    if(ds) stats.AddSequence(ReverseComplement()(seqs[j].Sequence()));
  }
//
  WriteToLog()(std::string("FASTA File: ")+fname);
  WriteIntToLog()(std::string("sequence count: "),nSeqs);
  WriteIntToLog()(std::string("character count: "),totalChar);
  WriteIntToLog()(std::string("shortest seq: "),shortSeq);
  WriteIntToLog()(std::string("longest seq: "),longSeq);
  WriteFloatToLog()(std::string("average length: "),avgSeq);
  WriteFloatToLog()(std::string("median  length: "),medSeq);
  WriteIntToLog()(std::string("Size histogram, bin size "),binW);
  for(unsigned k=0;k<nBins;k++) 
    WriteIntToLog()(std::string(),hist[k]);
  //
  if(nSeqs) {
	std::ostream& outRef(ofPrefix.length()?outFile:std::cout);
    if(!listAll) {
      strVec words= stats.ListWords();
      for(unsigned ww=0;ww<words.size();ww++) {
	outRef << words[ww].c_str() << '\t' << stats.WordCount(words[ww]);
	outRef << '\t' << stats.WordFrac(words[ww]);
	if(words[ww].length()>2)
	  outRef << '\t' << stats.ZScore(words[ww],words[ww].length()-2);
	outRef << std::endl;
      }
    } else 
      for(unsigned n=0;n<=ORand;++n) {
	int nWords= pow(4,n+1);
	for(int i=0;i<nWords;++i) {
	  string w= IdxToSeq()(i,n+1);
	  outRef << w << '\t' << stats.WordCount(w) << '\t' << stats.WordFrac(w) << std::endl;
	}
      }
    
    for(unsigned ox=0;ox<=ORand;++ox) 
      outRef << "Ambig-" << ox << '\t' << stats.AmbiguousCount(ox+1) << std::endl;
    if(oFile) outFile.close();
  }
  //
  WriteToLogTS()(std::string("Normal termination."));
  WriteToErrTS()(std::string("Normal termination."));
  return 0;
}
