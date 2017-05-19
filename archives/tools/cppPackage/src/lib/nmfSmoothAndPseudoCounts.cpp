#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../lib/util.h"
#include "../lib/arg.h"
#include "../lib/logFile.h"
//////////////////////////////////////////////////////////////////////////////////////
typedef std::vector<std::string> strVec;
typedef std::vector<double> dVec;
typedef std::vector<dVec> dMat;
//////////////////////////////////////////////////////////////////////////////////////
int smoothOrder=3,maxRows=1000000;
double pseudoCounts=1; bool makeGCT=false;
std::string ofPrefix ("smoothedS2");
//////////////////////////////////////////////////////////////////////////////////////
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  //  if(a.GetFloatValue(std::string("w"),std::string("insert sigma"),inSig,-1,100,true,true)) ++pCount;
  if(a.GetString(std::string("o"),std::string("output prefix"),ofPrefix)) ++pCount; 
  if(a.KeyPresent(std::string("gct"))) { ++pCount; makeGCT=true; }
  if(a.GetIntValue(std::string("r"),std::string("smoothOrder"),smoothOrder,1,20)) ++pCount; 
  if(a.GetIntValue(std::string("N"),std::string("maxRows"),maxRows,1,100000000)) ++pCount; 
  if(a.GetFloatValue(std::string("p"),std::string("pseudocounts"),pseudoCounts,0.0,100.0)) ++pCount;
  if(pCount!=a.numOptArgs()) {
    WriteToErr()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteToLog()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLog()(std::string("processed= "),pCount);
    WriteIntToLog()(std::string("entered= "),a.numOptArgs());
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ** smoothVector **: routine that returns a smoothing of an input vector of doubles.
//  Pascal's triangle is used for combining coefficients
//   order: [integer] the number of values on each side to include in the smoothing.
//   norm:  [boolean] switch to normalize the vector to unit sum
//   v:     [vector of doubles] input data
dVec smoothVector(const dVec& v,int order,bool norm=false) {
  if(order < 1) return v;
  if(v.size() < (2*order+1)) return v;
  if(order > 20) order=20;
  // generate pascal's triangle coefficients
  dVec coeff(2*order+1,1.0);
  for(int i=1;i<2*order;++i) 
    for(int m=i;m>0;--m)
      coeff[m] += coeff[m-1];
  // smoothed vector spaceholder, sum to normalize the vector to sum to 1 (probability)
  dVec sv(v.size(),0.0); double sum=0.0;
  for(int i=0;i<v.size();++i) {
    double cSum=0.0; // cSum is the sum of coefficients used, to deal with edge effects
    for(int j= -1*order;j<=order;++j) 
      if((i+j)>=0 && (i+j)<v.size()) {
        cSum += coeff[j+order];
        sv[i] += coeff[j+order] * v[i+j];
      }
    sv[i] /= cSum;
    if(norm) sum += sv[i];
  }
  if(norm && sum!=1.0) for(int i=0;i<sv.size();++i) sv[i] /= sum;
  return sv;
}		
//////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc,char *argv[]) {
  if(!LogSetup()(argc,argv)) { std::cout << "error opening logs" << std::endl; exit(-1);}
  //
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<2)||(arg1==std::string("-h"))||(arg1==std::string("-help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]" << std::endl;
    std::cout << "       " << argv[0] << " prefix [-o string] [-r int] [-p float] [-N int] [-gct]" << std::endl; 
    WriteToLog()(std::string("usage message displayed"));
    WriteToErr()(std::string("usage message displayed"));
    exit(0); 
  }
  //
  cjArgs arguments(argc,argv,1);
  ProcessArgs(arguments);
  // vectors to hold the results
  strVec rowLabels,colLabels;
  dMat processedData;
  // reserve at least 1000 rows worth of data
  rowLabels.reserve(1000); processedData.reserve(1000);
  // get the input file; the input argument is the prefix, needs to have ".s2.counts.txt" appended
  std::string ifPrefix= arguments.RequiredArgs()[0];
  std::string mfname = ifPrefix + std::string(".s2.counts.txt");
  std::cout << "Loading and processing data from file: " << mfname << "...";  std::cout.flush();
  std::ifstream inFile(mfname.c_str(),std::ios::in);
  // get the headers first (will be 1 entry too long- ignore the 0th, over the rowLabels)
  char lineIn[10001];
  inFile.getline(lineIn,10001);
  colLabels= split()(std::string(lineIn));
  int nColsP1= colLabels.size();
  // loop through the file until either all rows or maxRows are processed
  int rowCount=0;
  while(!inFile.eof() && rowCount < maxRows) {
    inFile.getline(lineIn,10001);
    strVec f= split()(std::string(lineIn));
    if(f.size()>=nColsP1) {
      rowLabels.push_back(f[0]);
      dVec v; v.reserve(nColsP1);
      for(int i=1;i<f.size();++i) v.push_back(pseudoCounts+atof(f[i].c_str()));
      dVec sm= smoothVector(v,smoothOrder);
      processedData.push_back(sm);
      ++rowCount;
    }
  }
  inFile.close();
  // set up the mandatory output file
  std::string ofname(ofPrefix+std::string(".counts.txt"));
  std::cout << "done. Processed " << processedData.size() << " records.\nWriting smoothed data file: " << ofname << "..."; 
  std::cout.flush();
  std::ofstream ofile(ofname.c_str(),std::ios::out);
  // write the headers first (skipping the first)
  ofile << colLabels[0];
  for(int i=1;i<nColsP1;++i) ofile << '\t' << colLabels[i];
  ofile << '\n';
  // now loop through the rest of the data
  for(int i=0;i<processedData.size();++i) {
    ofile << rowLabels[i];
    for(int j=0;j<processedData[i].size();++j) ofile << '\t' << processedData[i][j];
    ofile << '\n';
  }
  ofile.close();
  std::cout << "done.\n";
  // finally, create the gct format file if desired
  if(makeGCT) {
    //first load the s2 vector, which will be used for the description field in the gct format
    strVec s2Vector; s2Vector.reserve(processedData.size());
    std::string s2fname= ifPrefix + std::string(".s2.txt");
    std::cout << "GCF file creation: Loading S2 file: " << s2fname << "...";  std::cout.flush();
    std::ifstream s2file(s2fname.c_str(),std::ios::in);
    s2file.getline(lineIn,10000); //skip the first line (headers)
    int s2Count=0;
    while(!s2file.eof() && s2Count != processedData.size()) {
      s2file.getline(lineIn,10000);
      strVec f= split()(std::string(lineIn));
      if(f.size()>=2) s2Vector.push_back(f[1]);
    }
    // now create the gcf file
    std::string gfname= ofPrefix + std::string(".gct");
    std::cout << "done.\nWriting GCT file: " << gfname << "..."; std::cout.flush();
    std::ofstream gfile(gfname.c_str(),std::ios::out);
    gfile << "#1.2\n" << processedData.size() << '\t' << (colLabels.size()-1) << "\nNAME\tDescription";
    for(int i=1;i<colLabels.size();++i) gfile << '\t' << colLabels[i];
    gfile << '\n';
    for(int i=0;i<processedData.size();++i) {
      gfile << rowLabels[i] << '\t' << s2Vector[i];
      for(int j=0;j<processedData[i].size();++j) gfile << '\t' << processedData[i][j];
      gfile << '\n';
    }
    gfile.close();
    std::cout << "done.\n";
  }
  // all done
  return 0;
}
