//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <time.h>
#include "../lib/logFile.h"
#include "../lib/util.h"
#include "../lib/arg.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef std::vector<double> dVec;
typedef std::vector<int> iVec;
typedef std::vector<dVec> dMat;
typedef std::vector<std::string> strVec;
typedef std::vector<strVec> strMat;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool noLogs=false,quiet=false; 
//double threshold=0.5; //int zPos=0;
std::string ofPrefix("default"),inPrefix,inPrefix2,inPrefix3,inPrefix4;
int nWords=10,nModels=8;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ProcessArgs(const cjArgs& a) {
  unsigned pCount=0;
  if(a.KeyPresent(std::string("nologs"))) { noLogs=true; LogsOff(); ++pCount; }
  if(a.KeyPresent(std::string("quiet")) || a.KeyPresent(std::string("Q"))) { quiet=true; ++pCount; }
  if(a.GetString(std::string("o"),std::string("output prefix"),ofPrefix)) ++pCount;
  if(a.GetIntValue(std::string("r"),std::string("number of models"),nModels,1,20)) ++pCount;
  if(a.GetIntValue(std::string("nw"),std::string("number of words"),nWords,5,100)) ++pCount;
  //if(a.GetFloatValue(std::string("m"),std::string("mask threshold"),maskTh,-4.0,4.0,true,true)) ++pCount;
  //if(a.GetFloatValue(std::string("w"),std::string("word threshold"),wordTh,-4.0,4.0,true,true)) ++pCount;
  //if(a.GetFloatVec(std::string("t"),std::string("score threshold vector"),threshold)) ++pCount;
  if(a.GetString(std::string("i"),std::string("inPrefix4"),inPrefix4)) ++pCount;
  if(pCount!=a.numOptArgs()) {
    WriteToLogs()(std::string("arguments processed and entered don't agree- check for invalid args"));
    WriteIntToLogs()(std::string("entered= "),a.numOptArgs());
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc,char *argv[]) {
  std::string arg1;
  if(argc>1) arg1= std::string(argv[1]);
  if((argc<4)||(arg1==std::string("-h"))||(arg1==std::string("-help"))) { 
    std::cout << "Usage: " << argv[0] << " [-h|-help]\n";
    std::cout << "       " << argv[0] << " fileprefix1 fileprefix2 fileprefix3 [-o filestring] [-nw int] [-quiet] [-nologs]\n"; 
    exit(0); 
  }
  if(!LogSetup()(argc,argv)) { std::cout << "error creating logs. exiting.\n"; exit(-1); }
  cjArgs arguments(argc,argv,3);
  ProcessArgs(arguments);
  inPrefix= arguments.RequiredArgs()[0]; 
  inPrefix2= arguments.RequiredArgs()[1];
  inPrefix3= arguments.RequiredArgs()[2];
  std::cout << "inPrefixes: " << inPrefix << ", " << inPrefix2 << ", " << inPrefix3 << std::endl;
  if(ofPrefix==std::string("default")) ofPrefix= inPrefix;
  // set up the filenames
  std::string htmlFile   = ofPrefix + std::string(".html");
  std::string baseFile   = inPrefix + std::string(".bases.txt");
  std::string weightFile = inPrefix + std::string(".weights.txt");
  std::string nWeightFile= inPrefix + std::string(".nweights.txt");
  std::string rWordFile  = inPrefix + std::string(".rwords.txt");
  std::string motifFile  = inPrefix + std::string(".motifs");
  std::string modelFile  = inPrefix + std::string(".models.txt");
  std::string sCountFile = inPrefix + std::string(".counts.txt");
  std::string countFile  = inPrefix2 + std::string(".s2.counts.txt");
  std::string s2File     = inPrefix2 + std::string(".s2.txt");
  std::string inputFile  = (inPrefix4==std::string())?inPrefix3 + std::string(".fa"):inPrefix4;
  std::string statFile   = inPrefix3 + std::string(".stat");
  std::string logFile1   = inPrefix + std::string(".logs.txt");
  std::string errFile1   = inPrefix + std::string(".errs.txt");
  std::string logFile2   = inPrefix2 + std::string(".logs.txt");
  std::string errFile2   = inPrefix2 + std::string(".errs.txt");
  std::string logoBase   = inPrefix + std::string(".A.png");

  int nCols= 2 * (nModels/2 + nModels%2);
  double cellWidth= 100.0/static_cast<double>(nCols);
  std::ofstream hfile(htmlFile.c_str(),std::ios::out);
  hfile << "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">\n\n";
  hfile << "<html>\n<head>\n<title>NMF analysis: " << inPrefix << "</title>\n</head>\n\n";
  //
  hfile << "<body>\n<table border=\"0\" cellpadding=\"2\" bgcolor=\"#CCCCCC\">\n";
  hfile << "<tr bgcolor=\"#FFFFFF\">\n";
  hfile << "<th colspan=\"" << nCols << "\"><h1 align=\"center\">NMF motif analysis for sequence file: <a href=\"";
  hfile << inputFile << "\">" << inputFile << "</a></h1></th>\n";
  hfile << "</tr><tr bgcolor=\"#FFFFFF\">";
  hfile << "<th colspan=\"" << nCols << "\"><font size=\"+2\">Positioning Distribution</font></th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << (nCols/2) << "\">Weighted by relative contribution</th>\n";
  hfile << "<th colspan=\"" << (nCols/2) << "\">Normalized to uniform overall probability</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<td colspan=\"" << (nCols/2) << "\" align=\"center\"><img src=\"" << inPrefix << ".w.png\" width=\"500\" height=\"300\" /></td>\n";
  hfile << "<td colspan=\"" << (nCols/2) << "\" align=\"center\"><img src=\"" << inPrefix << ".n.png\" width=\"500\" height=\"300\" /></td>\n";
  hfile << "</tr><tr bgcolor=\"#FFFFFF\">\n<th colspan=\"" << nCols << "\"><font size=\"+2\">Sequence content</font></th>\n</tr><tr>\n";
  for(int i=0;i<nModels;++i)
    hfile << "<th width=\"" << cellWidth << "\">Element " << (i+1) << "</th>\n";
  hfile << "</tr><tr>\n"; hfile.flush();
  char label= 'A';
  for(int i=0;i<nModels;++i) 
    hfile << "<td align=\"center\"><img src=\"" << inPrefix << "." << label++ << ".png\" width=\"110\" height=\"55\" /></td>\n";
  hfile << "</tr><tr>\n";
  // read in the rword file
  std::ifstream rfile(rWordFile.c_str(),std::ios::in);
  strMat rankedWords(nModels,strVec());
  for(int i=0;i<nWords && !rfile.eof();++i) {
    char lineIn[5001];
    rfile.getline(lineIn,5000);
    strVec rw= split()(std::string(lineIn),'\t');
    std::cout << lineIn << ':' << rw.size() << std::endl; 
    for(int j=0;j<rw.size() && j< rankedWords.size();++j) 
      rankedWords[j].push_back(rw[j]);
  }
  // now print out the ranked words in the cells below the logos
  for(int i=0;i<nModels;++i) {
    hfile << "<td align=\"center\"><pre>\n";
    for(int j=0;j<nWords && j<rankedWords[0].size();++j)
      hfile << rankedWords[i][j] << '\n';
    hfile << "</pre></td>\n";
  }
  // finish the table
  hfile << "</tr><tr bgcolor=\"#FFFFFF\">\n";
  hfile << "<th colspan=\"" << nCols << "\"><font size=\"+2\">Raw output files</font></th>\n";
  hfile << "</tr><tr>\n";
  //
  int colSpan1=7,colSpan2=13;
  if(nCols==4)      { colSpan1=1; colSpan2=3; }
  else if(nCols==6) { colSpan1=2; colSpan2=4; }
  else if(nCols==8) { colSpan1=3; colSpan2=5; }
  else if(nCols==10){ colSpan1=3; colSpan2=7; }
  else if(nCols==12){ colSpan1=4; colSpan2=8; }
  else if(nCols==14){ colSpan1=5; colSpan2=9; }
  else if(nCols==16){ colSpan1=5; colSpan2=11; }
  else if(nCols==18){ colSpan1=6; colSpan2=12; }
  // 
  hfile << "<th colspan=\"" << colSpan1 << "\" align=\"right\" bgcolor=\"#CCCCCC\">File</th>\n";
  hfile << "<th align=\"left\" colspan=\"" << colSpan2 << "\" bgcolor=\"#EEEEEE\">Description</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#EEEEEE\" align=\"right\"><a href=\"" << weightFile << "\">" << weightFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#CCCCCC\" align=\"left\">Position distribution file (unnormalized)</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#CCCCCC\" align=\"right\"><a href=\"" << nWeightFile << "\">" << nWeightFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#EEEEEE\" align=\"left\">Position distribution file (normalized)</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#EEEEEE\" align=\"right\"><a href=\"" << baseFile << "\">" << baseFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#CCCCCC\" align=\"left\">Sequence content file (words by motif))</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#CCCCCC\" align=\"right\"><a href=\"" << rWordFile << "\">" << rWordFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#EEEEEE\" align=\"left\">Sequence content file (ranked words in each motif)</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#EEEEEE\" align=\"right\"><a href=\"" << motifFile << "\">" << motifFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#CCCCCC\" align=\"left\">Projected motifs (raw nmfBuildMotifs output)</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#CCCCCC\" align=\"right\"><a href=\"" << modelFile << "\">" << modelFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#EEEEEE\" align=\"left\">Processed motifs (trimmed of low IC edges)</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#EEEEEE\" align=\"right\"><a href=\"" << sCountFile << "\">" << sCountFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#CCCCCC\" align=\"left\">Smoothed Windowed Word Count file</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#CCCCCC\" align=\"right\"><a href=\"" << countFile << "\">" << countFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#EEEEEE\" align=\"left\">Raw Windowed Word Count file</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#EEEEEE\" align=\"right\"><a href=\"" << s2File << "\">" << s2File << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#CCCCCC\" align=\"left\">S2 (chi-squared statistic) file</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#CCCCCC\" align=\"right\"><a href=\"" << statFile << "\">" << statFile << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#EEEEEE\" align=\"left\">Nucleotide/dinucleotide statistics file</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#EEEEEE\" align=\"right\"><a href=\"" << logFile1 << "\">" << logFile1 << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#CCCCCC\" align=\"left\">Concatenated log files</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#CCCCCC\" align=\"right\"><a href=\"" << errFile1 << "\">" << errFile1 << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#EEEEEE\" align=\"left\">Concatenated err files</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#EEEEEE\" align=\"right\"><a href=\"" << logFile2 << "\">" << logFile2 << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#CCCCCC\" align=\"left\">Concatenated log files (pre-nmf, independent of r)</th>\n";
  hfile << "</tr><tr>\n";
  hfile << "<th colspan=\"" << colSpan1 << "\" bgcolor= \"#CCCCCC\" align=\"right\"><a href=\"" << errFile2 << "\">" << errFile2 << "</a></th>\n";
  hfile << "<th colspan=\"" << colSpan2 << "\" bgcolor= \"#EEEEEE\" align=\"left\">Concatenated err files (pre-nmf, independent of r)</th>\n";
  hfile << "</tr>\n</table>\n</body>\n</html>\n";
  //
  return 0;
}

