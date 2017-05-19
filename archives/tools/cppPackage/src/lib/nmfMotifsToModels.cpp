#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../lib/util.h"
#include "../lib/logFile.h"
#include "../lib/arg.h"
typedef std::vector<double> dVec;
typedef std::vector<dVec> dMat;
//////////////////////////////////////////////////////////////////////////////////////
double icTh= 0.1; int minW=4;
std::string ofPrefix ("MtoM");
//////////////////////////////////////////////////////////////////////////////////////
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  //  if(a.GetFloatValue(std::string("w"),std::string("insert sigma"),inSig,-1,100,true,true)) ++pCount;
  if(a.GetString(std::string("o"),std::string("output prefix"),ofPrefix)) ++pCount; 
  //  if(a.KeyPresent(std::string("gcf"))) { ++pCount; makeGCF=true; }
  if(a.GetIntValue(std::string("m"),std::string("minW"),minW,1,20)) ++pCount; 
  if(a.GetFloatValue(std::string("t"),std::string("IC threshold"),icTh,0.0,2.0)) ++pCount;
  if(pCount!=a.numOptArgs()) {
    WriteToErr()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteToLog()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLog()(std::string("processed= "),pCount);
    WriteIntToLog()(std::string("entered= "),a.numOptArgs());
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[]) {
  if(!LogSetup()(argc,argv)) { std::cout << "error opening logs" << std::endl; exit(-1);}
  //
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<2)||(arg1==std::string("-h"))||(arg1==std::string("-help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]" << std::endl;
    std::cout << "       " << argv[0] << " prefix [-o string] [-t float] [-m int]" << std::endl; 
    WriteToLog()(std::string("usage message displayed"));
    WriteToErr()(std::string("usage message displayed"));
    exit(0); 
  }
  //
  cjArgs arguments(argc,argv,1);
  ProcessArgs(arguments);
  // get the input file; the input argument is the prefix, needs to have ".motifs" appended
  std::string ifPrefix= arguments.RequiredArgs()[0];
  std::string ifname = ifPrefix + std::string(".motifs");
  std::ifstream inFile(ifname.c_str(),std::ios::in);
  // set up the output file: ofPrefix + ".models.txt"
  std::string ofname(ofPrefix + std::string(".models.txt"));
  std::ofstream ofile(ofname.c_str(),std::ios::out);
  std::cout << "Loading and processing data from file: " << ifname << ", writing output to: " << ofname << std::endl;
  // process the file
  char lineIn[10001]; std::string line(lineIn); int mCount=0;
  while(!inFile.eof()) {
    while(line.substr(0,5) != std::string("model")) {
      inFile.getline(lineIn,10001);
      line= std::string(lineIn);
    }
    std::string header(line);
    //  std::string header(lineIn);
    inFile.getline(lineIn,10001); line=std::string(lineIn);
    int state= 0; bool noHead=true;
    dMat motif; dVec icVec;
    while(line.substr(0,5) != std::string("model") && !inFile.eof()) {
      strVec f= split()(std::string(lineIn));
      //      std::cout << "line size= " << f.size() << '\n';
      if(f.size()>=7) {
	double ic= atof(f[5].c_str());
	icVec.push_back(ic);
	dVec v; 
	for(int i=1;i<f.size()-2;++i) v.push_back(atof(f[i].c_str()));
	if(v.size()>3) motif.push_back(v);
      }
      inFile.getline(lineIn,10001);
      line=std::string(lineIn);
    }
    if(icVec.size() != motif.size()) {
      std::cout << "ic and motif not same size.  exiting\n";
      exit(-1);
    }
    std::cout << "Motif " << mCount;
    ++mCount;
    int start=0,stop=icVec.size();
    while(icVec[start]<icTh && start<icVec.size()) ++start;
    while(stop>start && icVec[stop-1]<icTh) --stop;
    if((stop-start)<icVec.size()) {
      std::cout << " clipped: ";
      if(start != 0) std::cout << start << " positions at start;";
      if(stop!=icVec.size()) std::cout << (icVec.size()-stop) << " positions at end;";
    }
    if((stop-start)>=minW) {
      ofile << ">" << header << '\n';
      for(int i=start;i<stop;++i)
	ofile << motif[i][1] << '\t' << motif[i][3] << '\t' << motif[i][0] << '\t' << motif[i][2] << "\t // " << "ic= " << icVec[i] << '\n';
      std::cout << " okay, width= " << (stop-start) << '\n';
    } else
      std::cout << " dropped; width too small: " << (stop-start) << '\n';
  }
  inFile.close();
  ofile.close();
  std::cout << "done.\n";
  return 0;
}
