/**************************************************************************
  Purpose:	Develop C++ libraries with commonelly used functions
                for Bioinformatics and Computational Biology applications

  Author: 	Lucie Hutchins, Scientific Software Engineer
  Department:	Department of Research, Bioinformatics
  Institution:	The Jackson Laboratory
  Date:		September , 2011
  modified:     November , 2013
   contact: lucie.hutchins@jax.org
            lucie.hutchins@yahoo.com
            lucie.hutchins@genomeeffects.org
http://www.cplusplus.com/doc/tutorial/files/
*******************************************************************/
#include "ggenome.h"
#include <ctype.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>
////////////////
#include <sys/types.h>
#include <sys/param.h>
////////////////
#include <errno.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <cstring>
using namespace std;

char pathname[MAXPATHLEN];
void createDir(string dir);
void displaySequence(ostream& outf,string line,int len,string seq);
void downloadRemoteFile(string localFile,string url){string cmd="curl -s '"+url+"' -o "+localFile; system(cmd.c_str());}

bool checkLocalGenome(string organismVersion,string strain,string genomeBase){
 DIR* dp;string genomePath=genomeBase+"/"+organismVersion;
 if(!strain.empty())genomePath+="/"+strain+"_genome";dp = opendir(genomePath.c_str());
 bool genomeExists=(dp== NULL)?false:true;if(dp!= NULL){try{closedir(dp);}catch(exception&e){}}
}
int getNearestTx(featureIndexMap& transcriptIds,string chrom, long pos);
long getRegionEnd(featureIndexMap& transcriptIds,string chrom, long pos);
//null constructor default to current working directory
ggenome::ggenome(){
 DIR* dp; string toolBaseDir="";getwd(pathname);toolBaseDir=string(pathname);
  struct stat sb; 
  //local server setting 
  if(toolBaseDir.substr((toolBaseDir.length()-1),1)=="/")
     toolBaseDir =toolBaseDir.substr((0,toolBaseDir.length()-2));
  string tool_configfile,genomeUrl;tool_configfile=toolBaseDir+"/tool_config.txt";
  ///load user settings 
  loadConfigData(tool_configfile);
  genomeData.toolBase=genomeData.toolBase.empty()?toolBaseDir:genomeData.toolBase; 
  if(genomeData.toolBase.substr((genomeData.toolBase.length()-1),1)=="/")
     genomeData.toolBase =genomeData.toolBase.substr((0,genomeData.toolBase.length()-2));
  //stores path to tool data directory
  genomeData.dataDirName = genomeData.toolBase+"/tools_data";
  if((dp = opendir(genomeData.dataDirName.c_str())) == NULL)createDir(genomeData.dataDirName);
  //stores paths to local copy of genome -- I don't need this local config any more
  genomeData.genomeBase=genomeData.dataDirName+"/genomes";
  if(genomeData.localGenomeConfigFile.empty())genomeData.localGenomeConfigFile=genomeData.genomeBase+"/local_organisms_w_chromosome_data.txt";
  genomeData.organismsWithGenome=genomeData.genomeBase+"/service_organisms_w_chromosome_data.txt";
  genomeData.annotationsXmlFile=genomeData.dataDirName+"/annotations.xml";
  genomeData.annotationsTexFile=genomeData.dataDirName+"/annotations.txt";
  bool load=false; //update the local copy of annotions list every 14 days 
  if (stat(genomeData.annotationsTexFile.c_str(), &sb) == -1)load=true;
  else if((lastModifiedDayCount(genomeData.annotationsTexFile))>14||lastModifiedDayCount(genomeData.annotationsTexFile)<0)load=true;
  if(load){setAnnotatedOrg();}loadAnnotatedOrganisms();
  //remote server setting 
  genomeData.server=genomeData.server.empty()?"http://demon.jax.org":genomeData.server;
  genomeData.transcriptService="/transcriptdb/web-services/";
  genomeData.remoteGenomeConfigFile=genomeData.transcriptService+"service_organisms_w_chromosome_data.txt";
  genomeData.webserviceUrl=genomeData.server+genomeData.transcriptService;
  genomeUrl=genomeData.server+genomeData.remoteGenomeConfigFile;
  //check if remote config file exists if not, download a copy
  if(stat(genomeData.organismsWithGenome.c_str(), &sb) == -1)
    downloadRemoteFile(genomeData.organismsWithGenome,genomeUrl);
  loadOrganismAllMap();setaminoAcid();
 
}
ggenome::ggenome(string toolBaseDir,string host)
{
 DIR* dp;// string currentWorkingDir="";getwd(pathname);currentWorkingDir=string(pathname);
  struct stat sb; 
  //local server setting 
  if(toolBaseDir.empty()){toolBaseDir=string(getwd(pathname));}
  if(toolBaseDir.substr((toolBaseDir.length()-1),1)=="/")
     toolBaseDir =toolBaseDir.substr((0,toolBaseDir.length()-2));
  string tool_configfile,genomeUrl;tool_configfile=toolBaseDir+"/tool_config.txt";
  ///load user settings 
  loadConfigData(tool_configfile);genomeData.toolBase=genomeData.toolBase.empty()?toolBaseDir:genomeData.toolBase; 
  if(genomeData.toolBase.substr((genomeData.toolBase.length()-1),1)=="/")
     genomeData.toolBase =genomeData.toolBase.substr((0,genomeData.toolBase.length()-2));
  if(!host.empty())genomeData.host=host;
  genomeData.dataDirName = genomeData.toolBase+"/tools_data";genomeData.genomeBase=genomeData.dataDirName+"/genomes";
  if(genomeData.localGenomeConfigFile.empty())
     genomeData.localGenomeConfigFile=genomeData.genomeBase+"/local_organisms_w_chromosome_data.txt";
  if(!host.empty() &&(toLowerCase(host)=="ucsc"))
     genomeData.organismsWithGenome=genomeData.genomeBase+"/"+host+"_organisms_w_chromosome_data.txt";
  else genomeData.organismsWithGenome=genomeData.genomeBase+"/service_organisms_w_chromosome_data.txt";
  //stores path to tool data directory
  if((dp = opendir(genomeData.dataDirName.c_str())) == NULL)createDir(genomeData.dataDirName);
  genomeData.annotationsXmlFile=genomeData.dataDirName+"/annotations.xml";
  genomeData.annotationsTexFile=genomeData.dataDirName+"/annotations.txt";
  if(!host.empty()){
    genomeData.annotationsXmlFile=genomeData.dataDirName+"/"+host+"_annotations.xml";
    genomeData.annotationsTexFile=genomeData.dataDirName+"/"+host+"_annotations.txt";
  }
  bool load=false; //update the local copy of annotions list every 14 days 
  //remote server setting 
  genomeData.server=genomeData.server.empty()?"http://demon.jax.org":genomeData.server;
  if(genomeData.transcriptService.empty())genomeData.transcriptService="/transcriptdb/web-services/";
  if(!host.empty() &&(toLowerCase(host)=="ucsc"))
      genomeData.remoteGenomeConfigFile=genomeData.transcriptService+host+"_organisms_w_chromosome_data.txt";
  else genomeData.remoteGenomeConfigFile=genomeData.transcriptService+"service_organisms_w_chromosome_data.txt";
  genomeData.webserviceUrl=genomeData.server+genomeData.transcriptService;
  genomeUrl=genomeData.server+genomeData.remoteGenomeConfigFile;
  if (stat(genomeData.annotationsTexFile.c_str(), &sb) == -1)load=true;
  else if((lastModifiedDayCount(genomeData.annotationsTexFile))>14||lastModifiedDayCount(genomeData.annotationsTexFile)<0)load=true;
  if(load){setAnnotatedOrg();}loadAnnotatedOrganisms();
  //check if remote config file exists if not, download a copy
  if(stat(genomeData.organismsWithGenome.c_str(), &sb) == -1){
    downloadRemoteFile(genomeData.organismsWithGenome,genomeUrl);
    //cout<<"consturctor with host and dir downloaded:"<<genomeData.organismsWithGenome<<";"<<genomeUrl<<endl;
  }loadOrganismAllMap();setaminoAcid();
}
ggenome::ggenome(string host)
{
 DIR* dp; string currentWorkingDir="";getwd(pathname);currentWorkingDir=string(pathname);
  struct stat sb; 
  //local server setting 
  if(currentWorkingDir.substr((currentWorkingDir.length()-1),1)=="/")
     currentWorkingDir =currentWorkingDir.substr((0,currentWorkingDir.length()-2));
  string tool_configfile,genomeUrl;tool_configfile=currentWorkingDir+"/tool_config.txt";
  ///load user settings 
  loadConfigData(tool_configfile);
  genomeData.toolBase=genomeData.toolBase.empty()?currentWorkingDir:genomeData.toolBase; 
  if(genomeData.toolBase.substr((genomeData.toolBase.length()-1),1)=="/")
     genomeData.toolBase =genomeData.toolBase.substr((0,genomeData.toolBase.length()-2));
  if(!host.empty())genomeData.host=host;
  //stores path to tool data directory
  genomeData.dataDirName = genomeData.toolBase+"/tools_data";
  if((dp = opendir(genomeData.dataDirName.c_str())) == NULL)createDir(genomeData.dataDirName);
  genomeData.genomeBase=genomeData.dataDirName+"/genomes";
  if(!host.empty() &&(toLowerCase(host)=="ucsc"))
     genomeData.organismsWithGenome=genomeData.genomeBase+"/"+host+"_organisms_w_chromosome_data.txt";
  else genomeData.organismsWithGenome=genomeData.genomeBase+"/service_organisms_w_chromosome_data.txt";
  genomeData.annotationsXmlFile=genomeData.dataDirName+"/annotations.xml";
  genomeData.annotationsTexFile=genomeData.dataDirName+"/annotations.txt";
  if(!host.empty()){
    genomeData.annotationsXmlFile=genomeData.dataDirName+"/"+host+"_annotations.xml";
    genomeData.annotationsTexFile=genomeData.dataDirName+"/"+host+"_annotations.txt";
  }
  bool load=false; //update the local copy of annotions list every 14 days 
  //remote server setting 
  genomeData.server=genomeData.server.empty()?"http://demon.jax.org":genomeData.server;
  if(genomeData.transcriptService.empty())genomeData.transcriptService="/transcriptdb/web-services/";
  if(!host.empty() &&(toLowerCase(host)=="ucsc"))
        genomeData.remoteGenomeConfigFile=genomeData.transcriptService+host+"_organisms_w_chromosome_data.txt";
  else genomeData.remoteGenomeConfigFile=genomeData.transcriptService+"service_organisms_w_chromosome_data.txt";
  genomeData.webserviceUrl=genomeData.server+genomeData.transcriptService;
  genomeUrl=genomeData.server+genomeData.remoteGenomeConfigFile;
  if (stat(genomeData.annotationsTexFile.c_str(), &sb) == -1)load=true;
  else if((lastModifiedDayCount(genomeData.annotationsTexFile))>14||lastModifiedDayCount(genomeData.annotationsTexFile)<0)load=true;
  if(load){setAnnotatedOrg();}loadAnnotatedOrganisms();
  //check if remote config file exists if not, download a copy
  if(stat(genomeData.organismsWithGenome.c_str(), &sb) == -1){
    downloadRemoteFile(genomeData.organismsWithGenome,genomeUrl);
    //cout<<"consturctor with host downloaded:"<<genomeData.organismsWithGenome<<";"<<genomeUrl<<endl;
  }
  loadOrganismAllMap();setaminoAcid();
}
/***********************************************************************************
  loadConfigData: sets class member variable with data provided in the config file
  The default chromosome file prefix and sufix is empty string
*********************************************************************************/
void ggenome::loadConfigData(string config_file){ string line,token;
  if(!config_file.empty()){ ifstream inf(config_file.c_str()); 
     if(inf){
         while(getline(inf,line)){
            if(line.empty())continue;
            if(line.find("WEBSERVER=")!=string::npos){
              token="WEBSERVER=";genomeData.server=line.substr(line.find(token)+token.length());
            }else if(line.find("TRANSCRIPT_SERVICE=")!=string::npos){
              token="TRANSCRIPT_SERVICE=";
              genomeData.transcriptService=line.substr(line.find(token)+token.length());
            }else if(line.find("REMOTE_GENOMECONFIGFILE_PATH=")!=string::npos){
              token="REMOTE_GENOMECONFIGFILE_PATH=";
              genomeData.remoteGenomeConfigFile=line.substr(line.find(token)+token.length());
            }else if(line.find("TOOL_BASE=")!=string::npos){
              token="TOOL_BASE=";genomeData.toolBase=line.substr(line.find(token)+token.length());
            }else if(line.find("GENOME_BASE=")!=string::npos){
              token="GENOME_BASE=";genomeData.genomeBase=line.substr(line.find(token)+token.length());
            }else if(line.find("LOCAL_GENOMECONFIGFILE_PATH=")!=string::npos){
               token="LOCAL_GENOMECONFIGFILE_PATH=";
              genomeData.localGenomeConfigFile=line.substr(line.find(token)+token.length());
            }else if(line.find("CHROM_PREFIX=")!=string::npos){
               token="CHROM_PREFIX=";
              genomeData.chromosomeFilePrefix=line.substr(line.find(token)+token.length());
            }else if(line.find("CHROM_SUFFIX=")!=string::npos){ token="CHROM_SUFFIX=";
              genomeData.chromosomeFileSuffix=line.substr(line.find(token)+token.length());}
         }
      }
   }
 // now load the organismMap
}
/***********************************************************************************************
 trimStr: removes leading and trailing blanks, and cariage return character from a given string
************************************************************************************************/
std::string ggenome::trimStr(string line)
{
  string token(line);
  if(!token.empty()){
      string::size_type pos =token.find_last_not_of(' ');if(pos!= string::npos){token.erase(pos + 1);}
      pos = token.find_first_not_of(' ');if(pos != string::npos) token.erase(0, pos);
      string ch="\r",replacement="";
      if(token.find(ch)!=string::npos){
         while(token.find(ch)!=string::npos){token.replace(token.find(ch),ch.length(),replacement);}
      }
      ch="\n";
      if(token.find(ch)!=string::npos){
         while(token.find(ch)!=string::npos){token.replace(token.find(ch),ch.length(),replacement);}
      }
  }
 return token;
}
void createDir(string dir){
   long lSize=double(dir.length()+1)+200;
    // allocate memory to contain the portion of data desired:
   char *temp = (char*) malloc (sizeof(char)*(lSize+1));   
    sprintf(temp,"mkdir %s",dir.c_str());
    system(temp);free(temp);
}
// displaySequence(outf,bG,line,len,seq);
void displaySequence(ostream& outf,string header,int len,string seq){
 outf<<">"<<header<<endl;if(len==0){len=seq.length();}
 long current_start=0; string tempstr;
  while(current_start< seq.length()){
        if((current_start+len)>=seq.length())tempstr=seq.substr(current_start);
        else tempstr=seq.substr(current_start,len);
        outf<<tempstr<<endl;current_start+=len;
  }
}
std::string ggenome::filterChrom(string chr){
  chr=toUpperCase(chr);size_t pos =chr.find("CHR");
  if(pos!=string::npos)chr=chr.replace(pos,3,"");pos =chr.find(">"); 
  if(pos!=string::npos)chr=chr.replace(pos,1,"");
  string token="'",replacement="";
  if(chr.find(token)!=string::npos){
     while(chr.find(token)!=string::npos){chr.replace(chr.find(token),token.length(),replacement);}
  }
  token='"';
  if(chr.find(token)!=string::npos){
      while(chr.find(token)!=string::npos){chr.replace(chr.find(token),token.length(),replacement);}
  }
  return chr;
}
/************************************************************************
loads the list of organisms with genome data into memory(all versions)
Both organisms and organismAll maps are loaded.
Note: The organism files, both the local and remote copies,are expected to have the following fields:
  1. organismGroup                                                 -> optional
  2. organismName -- organism common name (mouse, rat,human, ...)  -> optional
  3. currentBuild -- mm9 for mouse, hg19 for human                 -> required
  4. genomeDataPath  -- absolute path to genome data directory     -> required
     (directory containing the *.dat files)
Data structure organismAll stores the paths to the chromosome data both
     local and remote.
*************************************************************************/
void ggenome::loadOrganismAllMap(){
 ifstream infile(genomeData.organismsWithGenome.c_str());const char* delimiter=",\t";
 string line;string organismGroup,organism,build="",remotechromosomeDataDir="";DIR* dp;
 if(infile){ getline(infile,line);if(!organismAll.empty()) organismAll.clear();if(!organisms.empty()) organisms.clear();
    while(getline(infile,line)){if(!line.empty())
     { int i=0;char * pch;pch= strtok((char*)line.c_str(),delimiter);
        while(pch != NULL){//parse line fields into variables
           if(i==0)organismGroup=string(pch);else if(i==1)organism=string(pch);
           else if(i==2)build=string(pch);else if(i==3)remotechromosomeDataDir=string(pch);
           pch = strtok (NULL, delimiter);++i;
         } //this is where I check the local chrom index directory for this organism version
         //localchromosomeDataDir
         organismAll[build].localChromosomeDataDir=genomeData.genomeBase+"/"+build;
         organisms[build].localChromosomeDataDir=genomeData.genomeBase+"/"+build;
         dp=opendir(organisms[build].localChromosomeDataDir.c_str());  
         if(dp==NULL)createDir(organisms[build].localChromosomeDataDir);
         organismAll[build].organismGroup=organismGroup;organisms[build].organismGroup=organismGroup;
         organismAll[build].chromosomeDataDir=remotechromosomeDataDir;organismAll[build].organism=organism;
         organisms[build].chromosomeDataDir=remotechromosomeDataDir;
         size_t c_pos=0, t_pos=0;string currentPrefix,thisPrefix; intVec posList,tposList;
         if(organisms.find(organism)!=organisms.end()){  //update organism/current version map
            getDigitPositions(organisms[organism].currentBuild, posList); getDigitPositions(build,tposList);
            if(!posList.empty()){c_pos=posList.front()-1; currentPrefix=organisms[organism].currentBuild.substr(0,c_pos);}
            if(!tposList.empty()){t_pos=tposList.front()-1;thisPrefix=build.substr(0,t_pos);}
            if((c_pos>0 && t_pos>0)&&(currentPrefix.compare(thisPrefix)==0)){ 
                 int currentVersion=getNumericVersion(organisms[organism].currentBuild),thisVesrion=getNumericVersion(build);
                 if((currentVersion-thisVesrion)<0){
                     organisms[organism].currentBuild=build; organisms[build].currentBuild=build;
                     organisms[organism].chromosomeDataDir=remotechromosomeDataDir;
            }}
         }else{ organisms[organism].currentBuild=build;organisms[build].currentBuild=build;
                organisms[organism].chromosomeDataDir=remotechromosomeDataDir;
         }
         organisms[organism].organismGroup=organismGroup; 
       }
      }
 }
 infile.close();
}

/********************************************************************
 Returns the current genome db build for a given organism (mm9 - mouse)
*********************************************************************/
std::string ggenome::getCurrentBuild(string organism){
   string currentBuild="";organism=toLowerCase(organism);
   if(organisms.empty())loadOrganismAllMap();
   if(organisms.find(organism)!=organisms.end())currentBuild=organisms[organism].currentBuild;
  return currentBuild;
}

/********************************************************************
 Returns the group of a given organism (mammal , ...)
*********************************************************************/
std::string ggenome::getOrganismGroup(string organism){
   organism=toLowerCase(organism);string organismGroup;
   if(organisms.empty())loadOrganismAllMap();
   if(organisms.find(organism)!=organisms.end())organismGroup=organisms[organism].organismGroup;
  return organismGroup;
}

/********************************************************************
 Returns the genome data path of a given organism (/data/seq/mammal/mouse/...)
*********************************************************************/
std::string ggenome::getRemoteGenomeDir(string build){
  if(organismAll.empty())loadOrganismAllMap();
  if(organismAll.find(build)!=organismAll.end())return organismAll[build].chromosomeDataDir;
  else  return "";
}

/**********************************************************************************************
 Returns the path to local copy of the gemome of a given organism strain and version(AJ ,mm9)
*************************************************************************************************/
std::string ggenome::getLocalGenomeDir(string build){
  if(organismAll.empty())loadOrganismAllMap();
  if(organismAll.find(build)!=organismAll.end())return organismAll[build].localChromosomeDataDir;
  else return "";
}

/**************************************************************************
Returns the organism name associated with the specified build (mm9 - mouse)
**************************************************************************/
std::string ggenome::getOrganismName(string organismVersion){
  string organism="";
 if(organisms.empty()||organismAll.empty())loadOrganismAllMap();
  if(organismAll.find(organismVersion)!=organismAll.end())
    organism=organismAll[organismVersion].organism;
  return organism;
}

/**********************************************************************
 Returns true if organism exists in our organisms with genome database
**********************************************************************/
bool ggenome::organismExists(string organism){
   bool orgExists=false; organism=toLowerCase(organism);
   if(organisms.empty())loadOrganismAllMap();
   if(organisms.find(organism)!=organisms.end())orgExists=true;
   else if(organismAll.find(organism)!=organismAll.end())orgExists=true;
  return orgExists;
}

/********************************************************************
list the current version of each organism
*********************************************************************/
void ggenome::listCurrentGenomes(){
  organismMap::iterator iter;
   string header="List of organism current assembly version ";
  string db=(genomeData.host.empty())?"Graber Transcript":"UCSC";
   header+= "- Using "+db+" database";
  cout<<"\n******************************************\n"<<header<<"\n*******\n";
  cout<<"organismGroup\torganism ( current Version )\n";
  if(organisms.empty())loadOrganismAllMap();
  for(iter=organisms.begin();iter!=organisms.end();++iter){
      cout<<iter->second.organismGroup<<"\t"<<iter->first<<" ("<<iter->second.currentBuild<<")"<<endl;
   }
}

/********************************************************************
list the all the versions of each organism
*********************************************************************/
void ggenome::listAllGenomes(){
  organismBMap::iterator iter;
  string header="List of organism with available chromosome assembly data";
  string db=(genomeData.host.empty())?"Graber Transcript":"UCSC";
   header+= "- Using "+db+" database";
  cout<<"\n******************************************\n"<<header<<"\n*******\n";
  cout<<"organismGroup\torganism ( Assembly Version )\n";
  if(organisms.empty())loadOrganismAllMap();
  for(iter= organismAll.begin();iter!= organismAll.end();++iter){
      cout<<iter->second.organismGroup<<"\t"<<iter->second.organism<<" ("<<iter->first<<")"<<endl;
   }
}

/********************************************************************
getChromosomeFileName formats and return the chromosome 
absolute file name path
*********************************************************************/
string ggenome::getChromosomeFileName(string chromosome,string strain,string build){
    string chromosome_f=filterChrom(chromosome); string chromPath="";
    chromPath=getLocalGenomeDir(build);if(!strain.empty())chromPath+="/"+strain+"_genome";
    string chrom=chromPath+"/dat/"+trimStr(genomeData.chromosomeFilePrefix)+trimStr(chromosome_f)+
                trimStr(genomeData.chromosomeFileSuffix);
  return chrom;   
}
/**********************************************************************
 lastModifiedDayCount returns the number of days since a given file was
  modified till today.
**********************************************************************/
int ggenome::lastModifiedDayCount(string filename){
    struct stat stat_base; 
    if(stat(filename.c_str(), &stat_base)<0)return -20;
    else{
       struct tm  *ts; ts = localtime(&stat_base.st_mtime);
       int modifiedDay=(*ts).tm_mday;int modifiedMonth=(*ts).tm_mon+1;
       time_t rawtime;struct tm *timeinfo; time ( &rawtime );timeinfo = localtime ( &rawtime ); 
       int today=(*timeinfo).tm_mday;int diff=-20;int thismonth=(*timeinfo).tm_mon+1;
       if(modifiedMonth==thismonth)diff=today-modifiedDay;
      return diff;
    }
}

void ggenome::lastModifiedSince(lastModified &fileLMD, string filename){
    struct stat stat_base; stat(filename.c_str(), &stat_base);
    struct tm  *ts; ts = localtime(&stat_base.st_mtime);
    fileLMD.month=(*ts).tm_mon+1; fileLMD.year=(*ts).tm_year+1900;
    fileLMD.day=(*ts).tm_mday;
}
/*************************************************************************************************
Returns the number of lines a given file contains.
Handy for some QA
***************************************************************************************************/
int ggenome::getFileLineCount(string filename){
  int linecount=0; string line;long lSize=double(filename.length()+1)+double(filename.length()+1)+100; 
  ifstream inf; char *temp = (char*) malloc (sizeof(char)*(lSize+1)); string line_count="temp_list.txt";
 if (temp != NULL){ 
     sprintf(temp,"wc -l %s >%s",filename.c_str(),line_count.c_str());
     try{system(temp);}catch(exception& e){}
     try{
          inf.open(line_count.c_str());
          while(getline(inf,line)){
                if(line.empty())continue; char *pch,*temp2;
                lSize=double(line.length()+1); temp2 = (char*) malloc (sizeof(char)*(lSize+1));
                if(temp2 == NULL) continue;
                strcpy(temp2,line.c_str());pch = strtok(temp2,"\t"); if(linecount<=0){linecount=atoi(pch);}
                free (temp2);
           }inf.close();
      }catch(exception &e){}
  }free (temp); 
 return linecount;
}
 
/********************************************************************************************************
#  specifies the mutation type(1->transition or  2->transversion)
#  Note: Purine: A,G ; Pyrimidine: C, T
#  Transition: In molecular biology,transition is a mutation changing a purine to another purine nucleotide 
#              (A <-> G) or a pyrimidine to another pyrimidine nucleotide (C <-> T). 
#              Approximately two out of every three (~67%)single nucleotide polymorphisms (SNPs) are transitions.
#
#  Transversion : In molecular biology, transversion refers to the substitution of a purine for a pyrimidine 
                  or vice versa. 
#                 It can only be reverted by a spontaneous reversion. Because this type of mutation changes 
#                 the chemical structure dramatically, the consequences of this change tend to be more severe 
#                 and less common than that of transitions.Transversions can 
#                 be caused by ionizing radiation and alkylating agents.
**************************************************************************************************************/
int ggenome::getMutationType(string ref_allele,string other_allele){
   int mutation_type=0; //Transversion
    ref_allele=toUpperCase(ref_allele);other_allele=toUpperCase(other_allele);
   if(((ref_allele == "A") &&(other_allele== "G"))or((ref_allele=="G") &&(other_allele=="A"))
          or((ref_allele=="C") &&(other_allele=="T"))
          or((ref_allele=="T") &&(other_allele=="C"))){mutation_type=1;  //transition
    }else mutation_type=2;
   return mutation_type;
}

/*************************************************************************************************
Returns the amino acid sequence equivalent of the CDS sequence input.
Most protein sequences are derived from translations of CoDing Sequence (CDS) derived from gene predictions. 
A CoDing Sequence (CDS) is a region of DNA or RNA whose sequence determines the sequence of amino acids in a protein. 
It should not be mixed up with an Open Reading Frame (ORF), which is a series of DNA codons that does not contain any STOP codons. 
All CDS are ORFs, but not all ORFs are CDS...
/*keep the stats of CDS sequences that have :
   1. both start and stop codons are first and last codons of the coding region
   2. start codon is first,stop codon is not last
   3. start codon is not first,stop codon is last 
   4. both start and stop codons are not first and last codons of the coding region
***************************************************************************************************/
std::string ggenome::cds2aaSeq(string seq){ string protein="";
  if(!seq.empty()){ if(codon2aa.empty())setaminoAcid();
     if(seq.length()>=3){seq=toUpperCase(seq);int remainder=seq.length()%3,i=0;
        string startcodon=seq.substr(0,3),tempstr=seq,starttoken="",endtoken="";
        if(remainder!=0){
            if(startCodon.find(startcodon)!=startCodon.end()){ //the first codon is ATG/AUG
                endtoken=(remainder==1)?"X":"XX"; ++remainder;tempstr=seq.substr(0,seq.length()-remainder);
            }else{ //the first codon is incomplete
              tempstr=seq.substr(remainder,seq.length()-1);starttoken=(remainder==1)?"X":"XX";}
         }while(i<tempstr.length()){string codon=codon2aa[tempstr.substr(i,3)];i+=3;protein+=codon;}
         if(endtoken!="")protein+=endtoken;if(starttoken!="")protein=starttoken+protein;
  }}
  return protein;
}
/******************************************************************************************************
  Class constructor with config file/Directory name:
*******************************************************************************************************/
genomeSequence::genomeSequence(string toolBaseDir,string host):ggenome(toolBaseDir,host){}
genomeFeature::genomeFeature(string toolBaseDir,string host):ggenome(toolBaseDir,host){}

genomeSequence::genomeSequence(string host):ggenome(host){}
genomeFeature::genomeFeature(string host):ggenome(host){}

genomeFile::genomeFile(){} //use the base class null contstructor
genomeFile::genomeFile(string toolBaseDir,string host):ggenome(toolBaseDir,host){}

genomeFeature::genomeFeature():ggenome(){}
genomeSequence::genomeSequence():ggenome(){}
/********************************************************************
 displays the content of the config file 
**************************************************************************/
void genomeSequence::displayConfigData()
{
   cout<<"Genome base directory is : "<<genomeData.genomeBase<<endl;
   cout<<"Genome list file: "<<genomeData.organismsWithGenome<<endl;
   cout<<"Chromosome file prefix: "<<genomeData.chromosomeFilePrefix<<endl;
   cout<<"Chromosome file sufix: "<<genomeData.chromosomeFileSuffix<<endl;
}
/*********************************************************************************************************
 getGenomeSeq: given a chromosome data file , getGenomeSeq will return a subsequence 
               between loc1 and loc2. The function returns an empty string if wrong input
 Input : chromosome_file_name, loc1,loc2  where loc1<loc2
 Note: loc2 is in 1-base cordinates standard,loc1 coordinate base is specified using start_base
***********************************************************************************************************/
string genomeSequence::getGenomeSeq(string ChromosomeSeqFileName,long loc1, long loc2,int start_base)
{
       string buffer; if(loc1<0||loc2<0)return "";
       if(loc1>loc2){ long temp=loc1; loc1=loc2; loc2=temp;} 
       if(start_base==1)--loc1;  // the first bp of the file follows array indexing starting at 0
	FILE* fp;long lSize=0;char * buffer2;size_t k;int n;lSize=long(fabs(double(loc2-loc1)));
	fp = fopen64(ChromosomeSeqFileName.c_str(),"rb"); //for large files, use fopen64
	if(fp==NULL)return "";                           // allocate memory to contain the portion of data desired:
        buffer2 = (char*) malloc (sizeof(char)*(lSize+1)); if (buffer2 == NULL)return ""; // copy the read region into the buffer:
        fseeko64(fp,loc1,sizeof(char));k = fread (buffer2,1,lSize,fp);
        if (k != lSize)return ""; buffer2[lSize]='\0'; buffer = string(buffer2);
        fclose(fp);free (buffer2); return buffer;
}
/****************************************************************************************************
 displayTranscript: will display the formatted genomic sequence of the specified transcript cordinates.
 The sequence can be one of the following:
 a.) a dna sequence of a concatenation of all the exons of that transcript. exon1exon2.....,exonk
     where k is the last exon of the transcript (ex)
 b.) a dna sequence of all the coding exons of a transcript- CDS sequence (cds)
 c.) an aminoacid sequence of all the coding exons of a transcript - aa sequence  (aa)
 d.) a dna sequence of every exon of that transcript.- one exon per line  (exon)
 e.) a dna sequence of every intron of that transcript.- one intron per line (intron)
 f.) a dna sequence of a concatenation of all the exons and introns of that transcript.Introns and exons are separated by a commas 
     format: exon1,intron1,exon2,intron2,.....,intronk-1,exonk where k is the last exon of the transcript (tex)
 g.) a genomic sequence delimitted by chromStart and chromEnd - default
 h.) a genomic sequence delimitted by txStart-offset and txStart+offset -> transcript start site (tss)
 i.) a genomic sequence delimitted by txeEnd-offset and txEnd+offset -> Terminal transcript site (tts)
 j.) exon-intron junction - exonEnd+-offset and -> (ei)
 k.) intron-exon junction - intronEnd+-offset and -> (ie)
 l.) a dna sequence of all exons in non coding region of a transcript- 5' and 3' utrs -> (utr)

Promoters are regions of DNA that promote transcription and, in eukaryotes, are found at -30, -75, and -90 base pairs upstream from the transcription start site (abbreviated to TSS)
 Input : bG is a structure holding the transcript dada 
 Note: genomic cordinates are in 1-base cordinates standard
*****************************************************************************************************/
// void displayTranscript(ostream& outf,struct tabularData& bG,int len,string seqType,int featuresStart_base);
void genomeSequence::displayTranscript(ostream& outf,struct tabularData& bG,int len,string seqType,int start_base,string strain,int offset){
   string seq="",temp,strandlabel="forward strand"; bool reverse_flag=false;if((bG.strand=="-")){reverse_flag=true;strandlabel="reverse strand";}
   char chromStart[15], chromEnd[15]; sprintf(chromStart,"%d",long(bG.chromStart));sprintf(chromEnd,"%d",long(bG.chromEnd));
   char exonRank[50],exonCount[8]; offset=(offset>0)?offset:100; //offset is used for tss,tts,ie,and ei queries
   if((seqType=="cds")||(seqType=="aa")){
       sprintf(chromStart,"%d",long(bG.cdsStarts[0]));sprintf(chromEnd,"%d",long(bG.cdsEnds[0]));
   }bG.chrom=filterChrom(bG.chrom); string header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel;
   string chrom_file=getChromosomeFileName(bG.chrom,strain,bG.organismVersion); 
   exonRankMap exonMap;getExonRanksMap(bG,exonMap);exonRankMap::iterator iter,initer;
   if(seqType=="cds"){//cds query: now compute the sequence of coding exons . 
      if(bG.cdsStarts[0]!=bG.cdsEnds[0]){
       seq="";getCodingSeq(seq,bG,chrom_file,start_base);
       sprintf(exonRank,":featureLength %d:feature %s",seq.size(),"CDS");
       header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel+exonRank;
       displaySequence(outf,header,len,seq);
      }
   }else if(seqType=="utr"){seq="";getUTRSeq(seq,bG,chrom_file,start_base);
       if(seq !=""){ std::size_t commasPos = seq.find(",");
         if(commasPos!=std::string::npos){
            string fivePrimSeq=seq.substr(0,commasPos),threePrimSeq=(seq.substr(commasPos).length()>1)?seq.substr(commasPos+1):"";
            if(!fivePrimSeq.empty()){
                if(reverse_flag){sprintf(chromStart,"%d",long(bG.cdsEnds[0]));sprintf(chromEnd,"%d",long(bG.chromEnd));}
                else{sprintf(chromStart,"%d",long(bG.chromStart));sprintf(chromEnd,"%d",long(bG.cdsStarts[0]));}
                sprintf(exonRank,":featureLength %d:feature %s",fivePrimSeq.size(),"5'UTR");
                header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel+exonRank;
                displaySequence(outf,header,len,fivePrimSeq);
            }
            if(!threePrimSeq.empty()){
                if(reverse_flag){sprintf(chromStart,"%d",long(bG.chromStart));sprintf(chromEnd,"%d",long(bG.cdsStarts[0]));}
                else{sprintf(chromStart,"%d",long(bG.cdsEnds[0]));sprintf(chromEnd,"%d",long(bG.chromEnd));}
                sprintf(exonRank,":featureLength %d:feature %s",threePrimSeq.size(),"3'UTR");
                header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel+exonRank;
                displaySequence(outf,header,len,threePrimSeq);
            }
          }
       }
   }else if(seqType=="aa"){//amino acid query: now compute the amino acid sequence of coding exons . 
       if(bG.cdsStarts[0]!=bG.cdsEnds[0]){
         seq="";getCodingSeq(seq,bG,chrom_file,start_base);seq=cds2aaSeq(seq);displaySequence(outf,header,len,seq);
       }
   }else if(seqType=="ex"){ //concatenate exons of each transcript. 
     if(!exonMap.empty()){
       for(iter=exonMap.begin();iter!=exonMap.end();++iter){
         if(reverse_flag){
              temp=getGenomeSeq(chrom_file,iter->second.start,iter->second.end,start_base);
              temp=reverseSeq(complementSeq(temp));seq+=temp;
          }else {seq+=getGenomeSeq(chrom_file,iter->second.start,iter->second.end,start_base);}
       }}
      sprintf(exonRank,":featureCount %d:featureLength %d:feature %s",long(exonMap.size()),seq.length(),"All Exons");
      header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel+exonRank;
      displaySequence(outf,header,len,seq);
   }else if(seqType=="exon"){ //do not concatenate exons
      if(!exonMap.empty()){
          for(iter=exonMap.begin();iter!=exonMap.end();++iter){
              long exonStart=iter->second.start,exonEnd=iter->second.end;int exonLen=fabs(exonEnd-exonStart);
              sprintf(chromStart,"%d",exonStart);sprintf(chromEnd,"%d",exonEnd);seq="";
            sprintf(exonRank,":rank %d:featureCount %d:featureLength %d:feature %s",long(iter->first),long(exonMap.size()),exonLen,seqType.c_str());
              header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel+exonRank;
              seq=getGenomeSeq(chrom_file,exonStart,exonEnd,start_base);
              if(reverse_flag){seq=reverseSeq(complementSeq(seq));} displaySequence(outf,header,len,seq);
      }}
   }else if(seqType=="intron"){ //do not concatenate introns
      if(!exonMap.empty()){
          if(exonMap.size()>1){ //single exon transcripts do not have exons
            for(iter=exonMap.begin();iter!=exonMap.end();++iter){
              if(iter->first==exonMap.size())continue;initer=iter; ++initer; 
              long intronStart=iter->second.end,intronEnd=initer->second.start;seq="";
              sprintf(chromStart,"%d",intronStart);sprintf(chromEnd,"%d",intronEnd);
              if(reverse_flag){
                 intronStart=initer->second.end;intronEnd=iter->second.start;
                 sprintf(chromStart,"%d",intronStart);sprintf(chromEnd,"%d",intronEnd);
              } int intronLen=fabs(intronEnd-intronStart);
       sprintf(exonRank,":rank %d:featureCount %d:featureLength %d:feature %s",long(iter->first),long(exonMap.size())-1,intronLen,seqType.c_str());
              header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel+exonRank;
              seq=getGenomeSeq(chrom_file,intronStart,intronEnd,start_base);
              if(reverse_flag){seq=reverseSeq(complementSeq(seq));}displaySequence(outf,header,len,seq);
      }}}
   }else if(seqType=="tex"){ //concatenate exon,intron,exon,.. 
     if(!exonMap.empty()){ long intronStart=0,intronEnd=0;string exons="";seq="";
       for(iter=exonMap.begin();iter!=exonMap.end();++iter){ 
         if(reverse_flag){
              temp=getGenomeSeq(chrom_file,iter->second.start,iter->second.end,start_base);
              temp=reverseSeq(complementSeq(temp));exons+=temp;
              if(iter->first==1)seq+=temp;else seq+=","+temp; //now get the intron sequence if applicable
              if(iter->first!=exonMap.size()){initer=iter; ++initer;
                 intronStart=initer->second.end;intronEnd=iter->second.start;
                 temp=getGenomeSeq(chrom_file,intronStart,intronEnd,start_base);
                 temp=reverseSeq(complementSeq(temp));seq+=","+temp;
              } 
          }else { temp=getGenomeSeq(chrom_file,iter->second.start,iter->second.end,start_base);exons+=temp;
              if(iter->first==1)seq+=temp;else seq+=","+temp;
              if(iter->first!=exonMap.size()){initer=iter; ++initer;
                 intronStart=initer->second.end;intronEnd=iter->second.start;
                 temp=getGenomeSeq(chrom_file,intronStart,intronEnd,start_base);seq+=","+temp;
          }}
       }sprintf(exonRank,":exonCount %d:exonsLength %d:feature %s",long(exonMap.size()),exons.length(),"Transcript");header+=exonRank;
      displaySequence(outf,header,len,seq);
  }}else if((seqType=="ie")||(seqType=="ei")){ //compute feature junction
         if(!exonMap.empty()){ long regionStart=0,regionEnd=0;seq="";
            for(iter=exonMap.begin();iter!=exonMap.end();++iter){
                if(iter->first==1){   //first exon has only exon-intron junction
                   if(seqType=="ei"){
                      regionStart=(reverse_flag)?(iter->second.start-offset):(iter->second.end-offset);
                      regionEnd=(reverse_flag)?(iter->second.start+offset):(iter->second.end+offset); 
                }}else if(iter->first==exonMap.size()){ //last exon has only intron-exon junction
                   if(seqType=="ie"){
                      regionStart=(reverse_flag)?(iter->second.end-offset):(iter->second.start-offset);
                      regionEnd=(reverse_flag)?(iter->second.end+offset):(iter->second.start+offset); 
                }}else{
                   if(seqType=="ie"){
                      regionStart=(reverse_flag)?(iter->second.end-offset):(iter->second.start-offset);
                      regionEnd=(reverse_flag)?(iter->second.end+offset):(iter->second.start+offset); 
                   }else{
                      regionStart=(reverse_flag)?(iter->second.start-offset):(iter->second.end-offset);
                      regionEnd=(reverse_flag)?(iter->second.start+offset):(iter->second.end+offset); 
                }}
                if(regionStart==0)continue;
                seq=getGenomeSeq(chrom_file,regionStart,regionEnd,start_base);
                seq=(reverse_flag)?reverseSeq(complementSeq(seq)):seq;
                sprintf(exonRank,":rank %d:exonCount %d:feature %s",long(iter->first),long(exonMap.size()),seqType.c_str());
                sprintf(chromStart,"%d",regionStart);sprintf(chromEnd,"%d",regionEnd);
                header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel+exonRank;
                displaySequence(outf,header,len,seq); regionStart=0;regionEnd=0;
        }}
  }else if((seqType=="tss")||(seqType=="tts")){ 
     //transcript start site +- offset - default offset  size 10bp
      long regionStart=(reverse_flag)?bG.chromEnd-offset:bG.chromStart-offset;
      long regionEnd=(reverse_flag)?bG.chromEnd+offset:bG.chromStart+offset;
      if(seqType=="tts"){ //terminal transcript site +- offset - default offset  size 110bp
         regionStart=(reverse_flag)?bG.chromStart-offset:bG.chromEnd-offset;
         regionEnd=(reverse_flag)?bG.chromStart+offset:bG.chromEnd+offset;
      }seq=getGenomeSeq(chrom_file,regionStart,regionEnd,start_base);
      seq=(reverse_flag)?reverseSeq(complementSeq(seq)):seq;sprintf(exonRank,":exonCount %d:feature %s",long(exonMap.size()),seqType.c_str());
      sprintf(chromStart,"%d",regionStart);sprintf(chromEnd,"%d",regionEnd);
      header=bG.name+"\t"+bG.chrom+":"+chromStart+"-"+chromEnd+":"+strandlabel+exonRank;
      displaySequence(outf,header,len,seq); 
   }else{ //genomic region
      seq=getGenomeSeq(chrom_file,long(bG.chromStart),long(bG.chromEnd),start_base);
      if(reverse_flag){seq=reverseSeq(complementSeq(seq));}displaySequence(outf,header,len,seq); 
   }
  //I need to implement junctions, introns, and terminal exons sequence extraction
}
/***********************************************************************************
 COMPUTE the coding sequence of a coding transcript. 
 Input: transcript structure, chromosome filename, feature starts coordinate base flag-
  0 means feature start follow zero-base standard, 1 mean one-base standard
************************************************************************************/
void genomeSequence::getCodingSeq(string &seq,struct tabularData & transcript,string chrom_file,int start_base){
 exonRankMap::iterator iter_lower,iter_upper;string temp;seq="";exonRankMap exonMap;
 bool reverse_flag=false;if(transcript.strand!="+")reverse_flag=true;
 exonCord target_exon; getOverlapExon(transcript,transcript.cdsStarts[0],target_exon);
 if(target_exon.rank>0){
    getExonRanksMap(transcript,exonMap);
    iter_lower= exonMap.lower_bound(target_exon.firstCoding_rank);iter_upper=exonMap.upper_bound(target_exon.lastCoding_rank);
    if(reverse_flag){
       iter_lower= exonMap.lower_bound(target_exon.lastCoding_rank);iter_upper=exonMap.upper_bound(target_exon.firstCoding_rank);}
    while(iter_lower!=iter_upper){
          if(reverse_flag){ //reverse strand
             if(iter_lower->first==target_exon.lastCoding_rank){ //this is the first coding exon
                  temp=getGenomeSeq(chrom_file,iter_lower->second.start,transcript.cdsEnds[0],start_base);
                  temp=reverseSeq(complementSeq(temp));seq+=temp;
             }else if(iter_lower->first==target_exon.firstCoding_rank){//this is the last coding exon
                  temp=getGenomeSeq(chrom_file,transcript.cdsStarts[0],iter_lower->second.end,start_base);
                  temp=reverseSeq(complementSeq(temp));seq+=temp;
             }else{temp=getGenomeSeq(chrom_file,iter_lower->second.start,iter_lower->second.end,start_base);
                  temp=reverseSeq(complementSeq(temp));seq+=temp;}
          }else{ //forward strand
                 if(iter_lower->first==target_exon.firstCoding_rank){//this is the first coding exon
                    seq+=getGenomeSeq(chrom_file,transcript.cdsStarts[0],iter_lower->second.end,start_base);
                 }else if(iter_lower->first==target_exon.lastCoding_rank){//this is the last coding exon
                    seq+=getGenomeSeq(chrom_file,iter_lower->second.start,transcript.cdsEnds[0],start_base);
                 }else{//fully coding exon
                 seq+=getGenomeSeq(chrom_file,iter_lower->second.start,iter_lower->second.end,start_base);}
         }++iter_lower;
   }}
}
/***********************************************************************************
 COMPUTE the sequence of the untranslated region of a coding transcript. 
 returns the 5' and the 3' utr sequences as one string separated with a commas -
 
 Input: transcript structure, chromosome filename, feature starts coordinate base flag-
  0 means feature start follow zero-base standard, 1 mean one-base standard

************************************************************************************/
void genomeSequence::getUTRSeq(string &seq,struct tabularData & transcript,string chrom_file,int start_base){
  exonRankMap::iterator iter_lower,iter_upper,iter;string temp,str;seq="";exonRankMap exonMap;
  bool reverse_flag=false;if(transcript.strand!="+")reverse_flag=true;
  exonCord target_exon;getOverlapExon(transcript,transcript.cdsStarts[0],target_exon);
  if(target_exon.rank>0){
     getExonRanksMap(transcript,exonMap); //returns exons of a transcript sorted by exon rank 
     //Returns an iterator pointing to the first exon in the container 
     //whose rank is not considered to go before firstCodingExon_rank (i.e., either it is equivalent or goes after).
     iter_lower= exonMap.lower_bound(target_exon.firstCoding_rank);
    //Returns an iterator pointing to the first exon in the container whose rank is considered to go after lastCodingExon_rank.
     iter_upper=exonMap.upper_bound(target_exon.lastCoding_rank);
     //get the 5' utr 
     if(reverse_flag){ --iter_upper;//adjust pointer to cds start exon 
         temp=getGenomeSeq(chrom_file,transcript.cdsEnds[0],iter_upper->second.end,start_base);
         temp=reverseSeq(complementSeq(temp));
         for(iter=exonMap.begin();iter!=iter_upper;++iter){//exon fully on 5'utr
             str=getGenomeSeq(chrom_file,iter->second.start,iter->second.end,start_base);
             str=reverseSeq(complementSeq(str));seq+=str;
         }seq+=temp;
     }else{ temp=getGenomeSeq(chrom_file,iter_lower->second.start,transcript.cdsStarts[0],start_base);
          for(iter=exonMap.begin();iter!=iter_lower;++iter){//exon fully on 5'utr
             seq+=getGenomeSeq(chrom_file,iter->second.start,iter->second.end,start_base);
          }seq+=temp;
      }seq+=",";  //get the 3' utr 
      if(reverse_flag){ 
         temp=getGenomeSeq(chrom_file,iter_lower->second.start,transcript.cdsStarts[0],start_base);//last coding exon
         temp=reverseSeq(complementSeq(temp));seq+=temp;++iter_lower;
         while(iter_lower!=exonMap.end()){ //exon fully on 3'utr
              str=getGenomeSeq(chrom_file,iter_lower->second.start,iter_lower->second.end,start_base);
              str=reverseSeq(complementSeq(str));seq+=str;++iter_lower;
          }
      }else{--iter_upper;//adjust pointer
          seq+=getGenomeSeq(chrom_file,transcript.cdsEnds[0],iter_upper->second.end,start_base); ++iter_upper;
          while(iter_upper!=exonMap.end()){//exon fully on 3'utr
              seq+=getGenomeSeq(chrom_file,iter_upper->second.start,iter_upper->second.end,start_base);++iter_upper;
      }}
   }
  
}
/***********************************************************************************
 returns the exon structure of the exon that ovelap the specified base pair position 
 or sets the exon structure members to their default value if specified base pair position 
    hits on intron or intergenic region
 Input: transcript structure, base pair position, exon structure -
  (start,end,rank,firstCoding_exon rank, last coding exon rank)
 
************************************************************************************/
void ggenome::getOverlapExon(struct tabularData & transcript,long bpPos,struct exonCord& exon){
   //returns the exon of this transcript that contain the feature, exons are stored sorted by exon_end
  if(transcript.cdsStarts[0]==transcript.cdsEnds[0]){
     exon.rank=0;exon.start=0;exon.end=0;
  }else{
     exonsMap exonMap; getExonMap(transcript,exonMap);exonsMap::iterator exonLowerIter;exonLowerIter=exonMap.lower_bound(bpPos);
     if(exonLowerIter!=exonMap.end()){
        if(exonLowerIter->second.first<=bpPos && bpPos<=exonLowerIter->first){
           exon.rank=exonLowerIter->second.second;exon.start=exonLowerIter->second.first;exon.end=exonLowerIter->first;}
        else{exon.rank=0;exon.start=0;exon.end=0;} 
     }else{exon.rank=0;exon.start=0;exon.end=0;}  
     if(!transcript.cdsStarts.empty()){//get first and last coding exons ranks
        exonLowerIter=exonMap.lower_bound(transcript.cdsStarts[0]);exon.firstCoding_rank=exonLowerIter->second.second;
        exonLowerIter=exonMap.lower_bound(transcript.cdsEnds[0]);exon.lastCoding_rank=exonLowerIter->second.second;
  }}
}

/***********************************************************************************
 returns the exon structure of the exon that is closest to the specified base pair position 
 Input: transcript structure, base pair position, exon structure -
  (start,end,rank,firstCoding_exon rank, last coding exon rank)

************************************************************************************/
void ggenome::getNearestExon(struct tabularData & transcript,long bpPos,struct exonCord& exon){
   //returns the exon of this transcript that is the closest to the feature, exons are stored sorted by exon_end
  exonsMap exonMap; getExonMap(transcript,exonMap); 
  exonsMap::iterator exonUpperIter;exonUpperIter=exonMap.upper_bound(bpPos);
  if(exonUpperIter!=exonMap.end()){
    exon.rank=exonUpperIter->second.second;exon.start=exonUpperIter->second.first;
    exon.end=exonUpperIter->first; --exonUpperIter;
    if((exon.start-bpPos)>(bpPos-exonUpperIter->first)){
      exon.rank=exonUpperIter->second.second;exon.start=exonUpperIter->second.first;
      exon.end=exonUpperIter->first;
    }
  }else{exon.rank=0;exon.start=0;exon.end=0;}//get first and last coding exons rank
  if(!transcript.cdsStarts.empty()){
     exonUpperIter=exonMap.lower_bound(transcript.cdsStarts[0]);exon.firstCoding_rank=exonUpperIter->second.second;
     exonUpperIter=exonMap.lower_bound(transcript.cdsEnds[0]);exon.lastCoding_rank=exonUpperIter->second.second;
  }
 
}
/***********************************************************************
returns the SNP distance from CDS start or zero if SNP on intron or UTR
************************************************************************/
int ggenome::getPosInCds(struct tabularData & transcript,long bpPos){
 exonRankMap::iterator iter_lower,iter_upper;int snpLocInCds=0;exonRankMap exonRanks;
  bool reverse_flag=false;if(transcript.strand!="+")reverse_flag=true;
 exonCord target_exon; 
 getOverlapExon(transcript,bpPos,target_exon);getExonRanksMap(transcript,exonRanks);
 if((target_exon.rank==0)||transcript.cdsStarts.empty())return 0;//SNP on intron
 if(!reverse_flag){
   if((target_exon.rank<target_exon.firstCoding_rank)||(target_exon.rank>target_exon.lastCoding_rank))return 0;//SNP in UTR
 }else{
   if((target_exon.rank<target_exon.lastCoding_rank)||(target_exon.rank>target_exon.firstCoding_rank))return 0;//SNP in UTR
 }
 try{ if(bpPos<transcript.cdsStarts[0]||bpPos>transcript.cdsEnds[0])return 0;
 }catch(exception& e){return 0;}
  iter_lower= exonRanks.lower_bound(target_exon.firstCoding_rank);iter_upper=exonRanks.upper_bound(target_exon.rank);
 if(reverse_flag){
  iter_lower= exonRanks.lower_bound(target_exon.lastCoding_rank);iter_upper=exonRanks.upper_bound(target_exon.rank);}
 while(iter_lower!=iter_upper){
     if(reverse_flag){ //reverse strand
         if(iter_lower->first==target_exon.rank){
            if(target_exon.rank==target_exon.lastCoding_rank)snpLocInCds+=transcript.cdsEnds[0]-bpPos;
            else snpLocInCds+=iter_lower->second.end-bpPos;
         }else if(iter_lower->first==target_exon.lastCoding_rank){ //this is the first coding exon
             snpLocInCds=transcript.cdsEnds[0]-iter_lower->second.start;
         }else{snpLocInCds+=iter_lower->second.end-iter_lower->second.start;}//this entire length of this exon is to be added
     }else{ //forward strand
           if(iter_lower->first==target_exon.rank){
               if(target_exon.rank==target_exon.firstCoding_rank)snpLocInCds+=bpPos-transcript.cdsStarts[0];
               else snpLocInCds+=bpPos-iter_lower->second.start;
           }else if(iter_lower->first==target_exon.firstCoding_rank){//this is the first coding exon
               snpLocInCds=iter_lower->second.end-transcript.cdsStarts[0];
           }else{snpLocInCds+=iter_lower->second.end-iter_lower->second.start;}//this entire length of this exon is to be added
     }++iter_lower;
 }
 if(reverse_flag)++ snpLocInCds;
 return snpLocInCds;
}
/***********************************************************************************************************
 public member: For a given transcript, returns the a map of exon ranks with the associated cordinates
 Assumption: exon coordinates in gene Prediction file are sorted by ascending order

Note: I may have to enforce this assumption by first storing exons into a map using exon starts as keys
************************************************************************************************************/
void ggenome::getExonRanksMap(struct tabularData & transcript,exonRankMap& exonMap){
  int i=0; int rank=1;if(!exonMap.empty())exonMap.clear();
  int exoncount=(int)transcript.exonStarts.size(), exonends=transcript.exonEnds.size();
  if(exoncount>0 &&(exoncount==exonends)){
     if((transcript.strand=="-")){ i=exonends;--i;
         while(i>=0){
              exonMap[rank].start=long(transcript.exonStarts[i]);
              exonMap[rank].end=long(transcript.exonEnds[i]);++rank;--i;
          }
      }else{i=0;
        while(i<exonends){
           exonMap[rank].start=long(transcript.exonStarts[i]);
           exonMap[rank].end=long(transcript.exonEnds[i]);++rank;++i;
          }
      }
   }
}
/***********************************************************************************************************
 public member: For a given transcript, returns a map of exon cordinates with the associated ranks
  exonEnds are keys
 Assumption: exon coordinates in gene Prediction file are sorted by ascending order
************************************************************************************************************/
void ggenome::getExonMap(struct tabularData & transcript,exonsMap& exonMap){
  int i=0; int rank=1;if(!exonMap.empty())exonMap.clear();
  int exoncount=(int)transcript.exonStarts.size(), exonends=transcript.exonEnds.size();
  if(exoncount>0 &&(exoncount==exonends)){
     if((transcript.strand=="-")){ i=exonends;--i;
         while(i>=0){
              exonMap[long(transcript.exonEnds[i])]=std::make_pair(long(transcript.exonStarts[i]),rank);
              ++rank;--i;
          }
      }else{i=0;
        while(i<exonends){
            exonMap[long(transcript.exonEnds[i])]=std::make_pair(long(transcript.exonStarts[i]),rank);
            ++rank;++i;
          }
   }}
}


/************************** reverseSeq ******************************************************************
 reverseSeq: given a DNA sequence , reverseSeq will return the reversed
  version of the original sequence without modifying the original sequence
 Input : DNA sequence
***********************************************************************************************************/
string genomeSequence::reverseSeq(string sequence){
       string tempstring;
       if(!sequence.empty()){
         for(int j= sequence.length()-1;j>=0;--j)tempstring +=sequence.at(j);
       }
  return tempstring;
}

/************************** ComplementSeq ***************************************************************
 ComplementSeq: given a DNA sequence , ComplementSeq will return the  complemented
                version of the original sequence without modifying the original sequence
 Input : DNA sequence
***********************************************************************************************************/
string genomeSequence::complementSeq(string sequence){
  string tempstring(sequence);int i =0;
  if(!sequence.empty()){
       for(i=0;i<sequence.length(); ++i ) {
               if (toupper(tempstring.at(i))=='T')tempstring[i]='A';
               else if (toupper(tempstring.at(i))=='A')tempstring[i]='T';
               else if (toupper(tempstring.at(i))=='C')tempstring[i]='G';
               else if (toupper(tempstring.at(i))=='G')tempstring[i]='C';
       }
   }
  return tempstring;
}

/************************** isCpGSite ***************************************************************
Assumption: one-base coordinate system

 isCpGSite: given a chromosome_file.dat file , and a bp_position , isCpGSite
            will return true if the SNP is located on a CG island
 Input : chromosome file, bp_position
Note:
    CpG sites or CG sites are regions of DNA where a cytosine nucleotide occurs next to a guanine nucleotide in the linear sequence of bases  along its length. "CpG" is shorthand for "—C—phosphate—G—", that is, cytosine and guanine separated by only one phosphate; phosphate links any two nucleosides together in DNA. The "CpG" notation is used to distinguish this linear sequence from the CG base-pairing of cytosine and guanine.
CpG dinucleotides have long been observed to occur with a much lower frequency in the sequence of vertebrate genomes than would be expected due to random chance.
***********************************************************************************************************/
bool genomeSequence::isCpGSite(string chromosomeSeqFileName, long bpPos, const char* strand,int start_base){
   string sequence = getGenomeSeq(chromosomeSeqFileName,bpPos,bpPos+1,start_base);
   if(sequence.length()==2){ 
     if(strcmp(strand,"+")!=0)sequence=complementSeq(reverseSeq(sequence));
     if(strcmp(toUpperCase(sequence).c_str(),"CG")==0)return true;
     else return false;
   }
 return false;
}

/****************************************************************************************
 returns the number of fields given the file line and a delimiter
Note: strtok modifies the original string, so be careful and make sure you use a
      copy of the original string instead
*****************************************************************************************/
int genomeFile::getFieldsCount(char* line, const char* delimiter){
 int count=0;long lSize=0;char *buffer2;lSize=long(fabs(double(strlen(line))));
  buffer2 =(char*) malloc (sizeof(char)*(lSize+1));
  strcpy(buffer2,line);char *pch;pch= strtok(buffer2,delimiter);
  while(pch!=NULL){
        ++count; pch = strtok (NULL, delimiter); 
   }free (buffer2); 
  return count;
}

/****************************************************************************************
 returns a map of fields and their indice given the file header and the delimiter
 Note: Do I really need the extra step to make a copy of the input string?
       I could get rid of the buffer step and save time and space 
     Store field names in lower cases 
Note: strtok modifies the original string, so be careful and make sure you use a
      copy of the original string instead
*****************************************************************************************/
void genomeFile::getFieldsIndex(char* headerLine, const char* delimiter,fieldIndexMap& fieldsMap){
  int current_index=0;
  long lSize=0;char * buffer2;lSize=long(fabs(double(strlen(headerLine))))+1;
  buffer2 =(char*) malloc (sizeof(char)*(lSize+1));strcpy(buffer2,headerLine);
  char * pch;pch= strtok(buffer2,delimiter);
  while (pch != NULL){ 
        fieldsMap[toLowerCase(string(pch))]=current_index; ++current_index; 
        pch = strtok (NULL, delimiter);
   }
  free (buffer2); 
}

/****************************************************************************************
parses fields from a given line given a delimiter into a vector of string
 Can be used for returning a vector of exon starts or exon ends
Note: strtok modifies the original string, so be careful and make sure you use a
      copy of the original string instead
*****************************************************************************************/
void  genomeFile::getFieldsData(char* line, const char* delimiter,strVec &fieldsData){
  if(!fieldsData.empty())fieldsData.clear();
 long lSize=0;char *buffer2;lSize=long(fabs(double(strlen(line))));
  buffer2 =(char*)malloc(sizeof(char)*(lSize+1)); strcpy(buffer2,line);
  char *pch;pch= strtok(buffer2,delimiter); 
  while (pch != NULL){fieldsData.push_back(string(pch));pch = strtok (NULL, delimiter);}
  free (buffer2); 
  //return fieldsData;
}
/****************************************************************************************
parses fields from a given line given a delimiter into a vector integers
Note: strtok modifies the original string, so be careful and make sure you use a
      copy of the original string instead
*****************************************************************************************/
void genomeFile::getIntFieldsData(char* line, const char* delimiter,intVec &fieldsData){
   string token(line);
  //I will get rid of tempstr when I modity isAllDigit to take cstring instead of string
  if(!fieldsData.empty())fieldsData.clear();
  if(line != NULL){
     if(isAllDigit(token))fieldsData.push_back(atoi(line));
     else{ 
        long lSize=0;char *buffer2;lSize=long(fabs(double(strlen(line))));
        buffer2 =(char*)malloc(sizeof(char)*(lSize+1));
        strcpy(buffer2,line);char *pch;pch= strtok(buffer2,delimiter);
        while (pch != NULL){ fieldsData.push_back(atoi(pch));pch = strtok (NULL, delimiter);}
        free (buffer2); 
     }
  }
   //return fieldsData;
}

/****************************************************************************************
returns true if a given SNP file does not have the appropriate format
*****************************************************************************************/
bool genomeFile::isBadSNPHeader(char* headerLine,const char* delimiter){
   fieldIndexMap fieldsMap;bool isBadHeader=false;
   getFieldsIndex(headerLine,delimiter,fieldsMap);
   if((fieldsMap.find("chromosome") == fieldsMap.end())&&
       (fieldsMap.find("chrom") == fieldsMap.end())&&
       (fieldsMap.find("affyid") == fieldsMap.end())&&
       (fieldsMap.find("chr") == fieldsMap.end()))isBadHeader=true;
  if((fieldsMap.find("position") == fieldsMap.end())&&
      (fieldsMap.find("pos") == fieldsMap.end())&&
      (fieldsMap.find("start") == fieldsMap.end())&&
      (fieldsMap.find("chromstart") == fieldsMap.end())&&
      (fieldsMap.find("bppos") == fieldsMap.end())&&
      (fieldsMap.find("location") == fieldsMap.end())&&
      (fieldsMap.find("bpposition") == fieldsMap.end()))isBadHeader=true;
  return isBadHeader;
}

/****************************************************************************************
returns true if a given tabular file does not have the appropriate format
*****************************************************************************************/
bool genomeFile::isBadTabularHeader(char* headerLine,const char* delimiter){
   fieldIndexMap fieldsMap;bool isBadHeader=false;
   getFieldsIndex(headerLine,delimiter,fieldsMap);
   if((fieldsMap.find("chromosome") == fieldsMap.end())&&
      (fieldsMap.find("chrom") == fieldsMap.end())&&
      (fieldsMap.find("chr") == fieldsMap.end()))isBadHeader=true;
   if((fieldsMap.find("chromstart") == fieldsMap.end())&&
       (fieldsMap.find("chrom_start") == fieldsMap.end())&&
       (fieldsMap.find("chrstart") == fieldsMap.end())&&
       (fieldsMap.find("chr_start") == fieldsMap.end())&&
       (fieldsMap.find("locationstart") == fieldsMap.end())&&
       (fieldsMap.find("location_start") == fieldsMap.end())&&
       (fieldsMap.find("locstart") == fieldsMap.end())&&
       (fieldsMap.find("loc_start") == fieldsMap.end())&&
       (fieldsMap.find("start") == fieldsMap.end()))isBadHeader=true;
   if((fieldsMap.find("chromend") == fieldsMap.end())&&
       (fieldsMap.find("chrom_end") == fieldsMap.end())&&
       (fieldsMap.find("chrend") == fieldsMap.end())&&
       (fieldsMap.find("chr_end") == fieldsMap.end())&&
       (fieldsMap.find("locationend") == fieldsMap.end())&&
       (fieldsMap.find("location_end") == fieldsMap.end())&&
       (fieldsMap.find("locend") == fieldsMap.end())&&
       (fieldsMap.find("loc_end") == fieldsMap.end())&&
       (fieldsMap.find("end") == fieldsMap.end()))isBadHeader=true;
  return isBadHeader;
}

/****************************************************************************************
Stores line fields into a tabular structure defined in ggenome class
set default index for each field chrom->0; chromStart->1; chromEnd->2; strand->3; name->4
*****************************************************************************************/
void genomeFile::getTabularFields(char* headerLine, char* line,const char* delimiter,struct tabularData& tabularFields){
    int index=0;const char* commas=",";
    fieldIndexMap fieldsMap; strVec fieldsData; getFieldsData(line,delimiter,fieldsData);
    getFieldsIndex(headerLine,delimiter,fieldsMap); 
    if(fieldsData.size()!=fieldsMap.size())return;
    bool has_chrom=true,has_chromStart=true,has_chromEnd=true;
  /* if((fieldsMap.find("chromosome") == fieldsMap.end())&&
      (fieldsMap.find("chrom") == fieldsMap.end())&&
      (fieldsMap.find("chr") == fieldsMap.end())){fieldsMap["Chromosome"]=0;has_chrom=false;}
   if(fieldsData.size()>1){
   if((fieldsMap.find("chromstart") == fieldsMap.end())&&
       (fieldsMap.find("chrom_start") == fieldsMap.end())&&
       (fieldsMap.find("chrstart") == fieldsMap.end())&&
       (fieldsMap.find("chr_start") == fieldsMap.end())&&
       (fieldsMap.find("locationstart") == fieldsMap.end())&&
       (fieldsMap.find("location_start") == fieldsMap.end())&&
       (fieldsMap.find("txstart") == fieldsMap.end())&&
       (fieldsMap.find("locstart") == fieldsMap.end())&&
       (fieldsMap.find("loc_start") == fieldsMap.end())&&
       (fieldsMap.find("start") == fieldsMap.end())){fieldsMap["chromStart"]=1;has_chromStart=false;}
   }
   if(fieldsData.size()>2){
    if((fieldsMap.find("chromend") == fieldsMap.end())&&
       (fieldsMap.find("chrom_end") == fieldsMap.end())&&
       (fieldsMap.find("chrend") == fieldsMap.end())&&
       (fieldsMap.find("chr_end") == fieldsMap.end())&&
       (fieldsMap.find("locationend") == fieldsMap.end())&&
       (fieldsMap.find("location_end") == fieldsMap.end())&&
       (fieldsMap.find("locend") == fieldsMap.end())&&
       (fieldsMap.find("txend") == fieldsMap.end())&&
       (fieldsMap.find("loc_end") == fieldsMap.end())&&
       (fieldsMap.find("end") == fieldsMap.end())){fieldsMap["chromEnd"]=2;has_chromEnd=false;}
    }
  if(!has_chromEnd && !has_chrom && !has_chromStart){
    if(fieldsData.size()>3){  //I'm not sure about this decision, what if only the strand info was missing
                             //and the rest of the fields were provided with the labels
      if(fieldsMap.find("strand")== fieldsMap.end())fieldsMap["Strand"]=3;
    }
    if(fieldsData.size()>4){  //'m not sure about this decision, same concern as for the strand
       if(fieldsMap.find("name")== fieldsMap.end())fieldsMap["Name"]=4;
     }
  }
   double  int score;
   string name, organism, organismVersion,gene,source;
       string strand; 
        //transcript features (exons, introns, cds,...)
  */
  if(fieldsMap.find("chromosome") != fieldsMap.end()){
        index=fieldsMap["chromosome"];
        if(index<fieldsData.size()) tabularFields.chrom=fieldsData[index];
        else tabularFields.chrom="-1";
    }else if(fieldsMap.find("chrom") != fieldsMap.end()){
        index=fieldsMap["chrom"];
        if(index<fieldsData.size()) tabularFields.chrom=fieldsData[index];
        else  tabularFields.chrom="-1";
    }else if(fieldsMap.find("chr") != fieldsMap.end()){
        index=fieldsMap["chr"];
        if(index<fieldsData.size()) tabularFields.chrom=fieldsData[index];
        else  tabularFields.chrom="-1";
    }
   if(fieldsMap.find("chromstart") != fieldsMap.end()){
        index=fieldsMap["chromstart"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1;  
     }else if(fieldsMap.find("chrom_start") != fieldsMap.end()){
        index=fieldsMap["chrom_start"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1;  
     }else if(fieldsMap.find("start") != fieldsMap.end()){
        index=fieldsMap["start"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1;  
     }else if(fieldsMap.find("chrstart") != fieldsMap.end()){
        index=fieldsMap["chrstart"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1; 
     }else if(fieldsMap.find("txstart") != fieldsMap.end()){
        index=fieldsMap["txstart"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1; 
     }
    else if(fieldsMap.find("chr_start") != fieldsMap.end()){
        index=fieldsMap["chr_start"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1;  
     }else if(fieldsMap.find("loc_start") != fieldsMap.end()){
        index=fieldsMap["loc_start"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1; 
     }else if(fieldsMap.find("locstart") != fieldsMap.end()){
        index=fieldsMap["locstart"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1;   
     }else if(fieldsMap.find("locationstart") != fieldsMap.end()){
        index=fieldsMap["locationstart"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1; 
     }else if(fieldsMap.find("location_start") != fieldsMap.end()){
        index=fieldsMap["location_start"];
        if(index<fieldsData.size())tabularFields.chromStart=atol(fieldsData[index].c_str());
        else tabularFields.chromStart=-1;  
     }
     if(fieldsMap.find("chromend") != fieldsMap.end()){
        index=fieldsMap["chromend"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
        
     }else if(fieldsMap.find("txend") != fieldsMap.end()){
        index=fieldsMap["txend"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
     }
     else if(fieldsMap.find("chrom_end") != fieldsMap.end()){
        index=fieldsMap["chrom_end"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
     }else if(fieldsMap.find("chrend") != fieldsMap.end()){
        index=fieldsMap["chrend"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
     }else if(fieldsMap.find("chr_end") != fieldsMap.end()){
        index=fieldsMap["chr_end"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
     }else if(fieldsMap.find("locationend") != fieldsMap.end()){
        index=fieldsMap["locationend"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
     }else if(fieldsMap.find("location_end") != fieldsMap.end()){
        index=fieldsMap["location_end"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
     }else if(fieldsMap.find("locend") != fieldsMap.end()){
        index=fieldsMap["locend"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
     }else if(fieldsMap.find("loc_end") != fieldsMap.end()){
        index=fieldsMap["loc_end"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;
    }else if(fieldsMap.find("end") != fieldsMap.end()){
        index=fieldsMap["end"];
        if(index<fieldsData.size())tabularFields.chromEnd=atol(fieldsData[index].c_str());
        else tabularFields.chromEnd=-1;   
     }
    if(fieldsMap.find("exonstarts")!= fieldsMap.end()){ index=fieldsMap["exonstarts"];
      getIntFieldsData((char*)fieldsData[index].c_str(), commas,tabularFields.exonStarts);
   }
   if(fieldsMap.find("exonends")!= fieldsMap.end()){ 
     index=fieldsMap["exonends"]; getIntFieldsData((char*)fieldsData[index].c_str(), commas,tabularFields.exonEnds);
   }
   if(fieldsMap.find("cdsstart")!= fieldsMap.end()){ index=fieldsMap["cdsstart"];
      getIntFieldsData((char*)fieldsData[index].c_str(), commas,tabularFields.cdsStarts);
   }else if(fieldsMap.find("cdsstarts")!= fieldsMap.end()){ index=fieldsMap["cdsstarts"];
      getIntFieldsData((char*)fieldsData[index].c_str(), commas,tabularFields.cdsStarts);
   }
   if(fieldsMap.find("cdsend")!= fieldsMap.end()){ index=fieldsMap["cdsend"];
     getIntFieldsData((char*)fieldsData[index].c_str(), commas,tabularFields.cdsEnds);
   }else  if(fieldsMap.find("cdsends")!= fieldsMap.end()){ index=fieldsMap["cdsends"];
     getIntFieldsData((char*)fieldsData[index].c_str(), commas,tabularFields.cdsEnds);
   }
   if(fieldsMap.find("score")!= fieldsMap.end())
        tabularFields.score=atoi(fieldsData[fieldsMap["score"]].c_str());
   if(fieldsMap.find("name")!= fieldsMap.end())tabularFields.name=fieldsData[fieldsMap["name"]];
   if(fieldsMap.find("source")!= fieldsMap.end())tabularFields.source=fieldsData[fieldsMap["source"]];
   if(fieldsMap.find("gene")!= fieldsMap.end())tabularFields.gene=fieldsData[fieldsMap["gene"]];
   else if(fieldsMap.find("name2")!= fieldsMap.end()){tabularFields.gene=fieldsData[fieldsMap["name2"]];} 
    if(fieldsMap.find("organism")!= fieldsMap.end())tabularFields.organism=fieldsData[fieldsMap["organism"]];
    if(fieldsMap.find("organismversion")!= fieldsMap.end())
       tabularFields.organismVersion=fieldsData[fieldsMap["organismversion"]];
    if(fieldsMap.find("strand") != fieldsMap.end()){
       index=fieldsMap["strand"];
       if(index<fieldsData.size())tabularFields.strand=fieldsData[index];
       else  tabularFields.strand="+";
    }else  tabularFields.strand="+";
   
}
/****************************************************************************************
Returns the SNP structure of this line
*****************************************************************************************/
void genomeFile::getSNPFields(char* headerLine,char* line,const char* delimiter,struct snpData& snpFields){
    int index=0; fieldIndexMap fieldsMap; strVec fieldsData;getFieldsData(line,delimiter,fieldsData);
    getFieldsIndex(headerLine,delimiter,fieldsMap);
   if(fieldsMap.find("snp_localid") != fieldsMap.end()){     //SNP local id (an integer value from cgdsnpdb)
         snpFields.snp_localID=atol(fieldsData[fieldsMap["snp_localid"]].c_str());
    }else if(fieldsMap.find("snplocalid") != fieldsMap.end()){
         snpFields.snp_localID=atol(fieldsData[fieldsMap["snplocalid"]].c_str());
    }else{snpFields.snp_localID=0;}
    if(fieldsMap.find("id") != fieldsMap.end()){ //SNP accession id
         snpFields.submitterID=fieldsData[fieldsMap["id"]];
    }else if(fieldsMap.find("snpid") != fieldsMap.end()){
          snpFields.submitterID=fieldsData[fieldsMap["snpid"]];
    }else if(fieldsMap.find("snpaccession") != fieldsMap.end()){
          snpFields.submitterID=fieldsData[fieldsMap["snpaccession"]];
    }else{snpFields.submitterID="";}
    if(fieldsMap.find("chromosome") != fieldsMap.end()){
        index=fieldsMap["chromosome"];
        if(index<fieldsData.size())snpFields.chrom=fieldsData[index];
        else snpFields.chrom="-1";
    }else if(fieldsMap.find("chrom") != fieldsMap.end()){
        index=fieldsMap["chrom"];
        if(index<fieldsData.size())snpFields.chrom=fieldsData[index];
        else snpFields.chrom="-1";
    }else if(fieldsMap.find("chr") != fieldsMap.end()){
        index=fieldsMap["chr"];
        if(index<fieldsData.size())snpFields.chrom=fieldsData[index];
        else snpFields.chrom="-1";
    }else snpFields.chrom="-1";
    if(fieldsMap.find("position")!= fieldsMap.end()){
      index=fieldsMap["position"];
      if(index<fieldsData.size())snpFields.pos=atol(fieldsData[index].c_str());
       else snpFields.pos=-1;
    }else if(fieldsMap.find("pos") != fieldsMap.end()){
        index=fieldsMap["pos"];
       if(index<fieldsData.size())snpFields.pos=atol(fieldsData[index].c_str());
       else snpFields.pos=-1;
    }
    else if(fieldsMap.find("bppos") != fieldsMap.end()){
       index=fieldsMap["bppos"];
       if(index<fieldsData.size())snpFields.pos=atol(fieldsData[index].c_str());
       else snpFields.pos=-1;
    }
    else if(fieldsMap.find("bpposition") != fieldsMap.end()){
        index=fieldsMap["bpposition"];
       if(index<fieldsData.size())snpFields.pos=atol(fieldsData[index].c_str());
       else snpFields.pos=-1;
    }
    else if(fieldsMap.find("location") != fieldsMap.end()){
        index=fieldsMap["location"];
       if(index<fieldsData.size())snpFields.pos=atol(fieldsData[index].c_str());
       else snpFields.pos=-1;
    }else if(fieldsMap.find("loc") != fieldsMap.end()){
        index=fieldsMap["loc"];
       if(index<fieldsData.size())snpFields.pos=atol(fieldsData[index].c_str());
       else snpFields.pos=-1;
    }
   else snpFields.pos=-1;
    if(fieldsMap.find("strand") != fieldsMap.end()){
       index=fieldsMap["strand"];
       if(index<fieldsData.size())snpFields.strand=fieldsData[index];
       else  snpFields.strand="+";
    }else  snpFields.strand="+";
   if(fieldsMap.find("reference_base") != fieldsMap.end()){
      index=fieldsMap["reference_base"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
   }else if(fieldsMap.find("referencebase") != fieldsMap.end()){
     index=fieldsMap["referencebase"];
     if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
     else  snpFields.referenceBase="";
    }
    else if(fieldsMap.find("ref_base") != fieldsMap.end()){
      index=fieldsMap["ref_base"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
    }
    else if(fieldsMap.find("refbase") != fieldsMap.end()){
      index=fieldsMap["refbase"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
   }
    else if(fieldsMap.find("ref_allele") != fieldsMap.end()){
      index=fieldsMap["ref_allele"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
   }
    else if(fieldsMap.find("allele") != fieldsMap.end()){
      index=fieldsMap["allele"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
   }
   else if(fieldsMap.find("refallele") != fieldsMap.end()){
      index=fieldsMap["refallele"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
   }
   else if(fieldsMap.find("reference_allele") != fieldsMap.end()){
      index=fieldsMap["reference_allele"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
   }
  else if(fieldsMap.find("referenceallele") != fieldsMap.end()){
      index=fieldsMap["referenceallele"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
   }
   else if(fieldsMap.find("ref") != fieldsMap.end()){
      index=fieldsMap["ref"];
      if(index<fieldsData.size())snpFields.referenceBase=fieldsData[index];
      else  snpFields.referenceBase="";
   }
   else if(fieldsMap.find("alleles") != fieldsMap.end()){
     index=fieldsMap["alleles"];
     if(index<fieldsData.size()){size_t pos =fieldsData[index].find_first_of("/");
        if(pos!=string::npos){
           snpFields.referenceBase=fieldsData[index].substr(0,pos);
           snpFields.consensusBase=fieldsData[index].substr(pos+1);
        }else{snpFields.referenceBase="";snpFields.consensusBase="";}
     }
   }
   if(fieldsMap.find("consensus_base") != fieldsMap.end()){
     index=fieldsMap["consensus_base"];
     if(index<fieldsData.size())snpFields.consensusBase=fieldsData[index];
     else  snpFields.consensusBase="";
   }
   else if(fieldsMap.find("consensusbase") != fieldsMap.end()){
     index=fieldsMap["consensusbase"];
     if(index<fieldsData.size())snpFields.consensusBase=fieldsData[index];
     else  snpFields.consensusBase="";
   }
   else if(fieldsMap.find("snpallele") != fieldsMap.end()){
     index=fieldsMap["snpallele"];
     if(index<fieldsData.size())snpFields.consensusBase=fieldsData[index];
     else  snpFields.consensusBase="";
   }
   else if(fieldsMap.find("snp_allele") != fieldsMap.end()){
     index=fieldsMap["snp_allele"];
     if(index<fieldsData.size())snpFields.consensusBase=fieldsData[index];
     else  snpFields.consensusBase="";
   }
   else if(fieldsMap.find("other_allele") != fieldsMap.end()){
     index=fieldsMap["other_allele"];
     if(index<fieldsData.size())snpFields.consensusBase=fieldsData[index];
     else  snpFields.consensusBase="";
   }
   else if(fieldsMap.find("alt") != fieldsMap.end()){
     index=fieldsMap["alt"];
     if(index<fieldsData.size())snpFields.consensusBase=fieldsData[index];
     else  snpFields.consensusBase="";
   }
   else if(fieldsMap.find("otherallele") != fieldsMap.end()){
     index=fieldsMap["otherallele"];
     if(index<fieldsData.size())snpFields.consensusBase=fieldsData[index];
     else  snpFields.consensusBase="";
   }
}

/****************************************************************************************
Returns the psl structure of this line
*****************************************************************************************/
void genomeFile::getPslFields(char* headerLine,char* line,const char* delimiter,struct pslData & psl){
    fieldIndexMap fieldsMap; strVec fieldsData;getFieldsData(line,delimiter,fieldsData);
    getFieldsIndex(headerLine,delimiter,fieldsMap);
    psl.matches=atoi(fieldsData[fieldsMap["matches"]].c_str());
    psl.misMatches=atoi(fieldsData[fieldsMap["misMatches"]].c_str());
    psl.repMatches=atoi(fieldsData[fieldsMap["repMatches"]].c_str());
    psl.nCount=atoi(fieldsData[fieldsMap["nCount"]].c_str());
    psl.qNumberInsert=atoi(fieldsData[fieldsMap["qNumberInsert"]].c_str());
    psl.qBaseInsert=atoi(fieldsData[fieldsMap["qBaseInsert"]].c_str());
    psl.tNumberInsert=atoi(fieldsData[fieldsMap["tNumberInsert"]].c_str());
    psl.tBaseInsert=atoi(fieldsData[fieldsMap["tBaseInsert"]].c_str());
    psl.strand=fieldsData[fieldsMap["strand"]];
    psl.qName=fieldsData[fieldsMap["qName"]];
    psl.qSize=atoi(fieldsData[fieldsMap["qSize"]].c_str());
    psl.qStart=atoi(fieldsData[fieldsMap["qStart"]].c_str());
    psl.qEnd=atoi(fieldsData[fieldsMap["qEnd"]].c_str());
    psl.tName=fieldsData[fieldsMap["tName"]];
    psl.tSize=atol(fieldsData[fieldsMap["tSize"]].c_str());
    psl.tStart=atol(fieldsData[fieldsMap["tStart"]].c_str());
    psl.tEnd=atol(fieldsData[fieldsMap["tEnd"]].c_str());
    psl.blockCount=atoi(fieldsData[fieldsMap["nblockCount"]].c_str());
    psl.blockSizes=fieldsData[fieldsMap["blockSizes"]];
    psl.qStarts=fieldsData[fieldsMap["qStarts"]];
    psl.tStarts=fieldsData[fieldsMap["tStarts"]]; 
}
///convert a gtf file (feature specific) into a genePred format file
// read the gtf by transcript then sort exons to generate exon starts and
// exon ends - gtf files do not have a header but follow a standard
// Each line is a feature and contains 9 columns
//feature [CDS,"CDS", "start_codon", "stop_codon",exon,..
/// attributes:gene_id "ENSMUSG00000078517";transcript_id "ENSMUST00000144335"; exon_count "2"; exon_rank "1"; 
/*
4	mm9-ensGene	CDS	138908535	138908632	.	+	-1	gene_id "ENSMUSG00000078517";transcript_id "ENSMUST00000042096"; exon_count "23"; exon_rank "1";
4	mm9-ensGene	exon	138910076	138910200	.	+	2	gene_id "ENSMUSG00000078517";transcript_id "ENSMUST00000042096"; exon_count "23"; exon_rank "2";
4	mm9-ensGene	CDS	138910076	138910200	.	+	2	gene_id "ENSMUSG00000078517";transcript_id "ENSMUST00000042096"; exon_count "23"; exon_rank "2";
4	mm9-ensGene	exon	138910727	138910792	.	+	1	gene_id "ENSMUSG00000078517";transcript_id "ENSMUST00000042096"; exon_count "23"; exon_rank "3";
AB000381 Twinscan  CDS          700   707   .   +   2  gene_id "001"; transcript_id "001.1";
AB000381 Twinscan  start_codon  380   382   .   +   0  gene_id "001"; transcript_id "001.1";
AB000381 Twinscan  stop_codon   708   710   .   +   0  gene_id "001"; transcript_id "001.1";
start and end are in 1-base standard

Plus strand: start is always <= end
CDS STOP= max(stop_codon.start.start_codon.end)
CDS START= min(stop_codon.start,stop_codon.end)

MInus strand: make sure start is always <=end and that start_codon<= stop_codon
CDS STOP= min(stop_codon.start,stop_codon.end)
CDS START= max(start_codon.start,start_codon.end)

Fields must be tab-separated. Also, all but the final field in 
each feature line must contain a value; "empty" columns should be denoted with a '.'

seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
source - name of the program that generated this feature, or the data source (database or project name)
feature - feature type name, e.g. Gene, Variation, Similarity
start - Start position of the feature, with sequence numbering starting at 1.
end - End position of the feature, with sequence numbering starting at 1.
score - A floating point value.
strand - defined as + (forward) or - (reverse).
frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
*/
void genomeFile::gtf2genePred(istream& gtff,ostream& genPredf){
  const char* delimiter="\t";strVec fieldsData;string line;
  while(getNextBlock(gtff,transcript)){  //read  by transcript block
        displayTranscript(transcript,genPredf);
        //getFieldsData((char*)line.c_str(),delimiter,fieldsData);
 }}
//read transcript features into transcript structure from a gtf file
bool genomeFile::getNextBlock(istream& gtff,const char* delimiter,struct tabularData& transcript){
   bool hasMore=true; //flag to signal more features to read for a transcript
   while(){
   }
}
//displays a transcript in gtf format or genePrediction format
void genomeFile::displayTranscript(struct tabularData& transcript,ostream& genPredf,bool gtf){
}
/*********************************************************
Notes: binFile takes a tab-delimitted (fileType=1) or 
       a commas separated file(fileType=2) input file with two
       columns (x_axis=length/position and  y_axis=count/frequency) then bin
       data according to the user specified step. Assuming no duplicate rows. Data is sorted by
       x_axis column  
***************************************************************/
void genomeFile::binFile(istream& inf,int xIndex,int yIndex,int stepSize,int fileType,ostream& out){
     const char* delimiter="\t"; if(fileType==2)delimiter=",";
     intIntMap binMap, xyMap;strVec fieldsData; //here we are loading data in memory for fast processing
                                                //assuming data not too big
     string line;
     while(getline(inf,line)){ //load data into a map to ensure sorted data by x_axis
         getFieldsData((char*)line.c_str(),delimiter,fieldsData);
         if((xIndex<fieldsData.size())&&(yIndex<fieldsData.size())){
            if(isAllDigit(fieldsData[xIndex])&&isAllDigit(fieldsData[yIndex])){
               xyMap[atoi(fieldsData[xIndex].c_str())]=atoi(fieldsData[yIndex].c_str());
             }
         }
     }
    intIntMap::iterator iter;int index=0;// now bin data file with 10 steps
   for(iter=xyMap.begin();iter!=xyMap.end();++iter){
       if(index==0)index=iter->first;
       if(iter->first >=(index+stepSize)) index+=stepSize;
       binMap[index]+= iter->second;
    }
   for(iter=binMap.begin();iter!=binMap.end();++iter){
       out<<iter->first<<delimiter<<iter->second<<"\n";
    }
   
}
/*****************************************************************************
Returns the header line of this file
*****************************************************************************/
string genomeFile::getTabularHeaderLine(const char* delimiter,istream& inf){
  string line,temp; string chrom[2],pos[2]; chrom[0]="chromosome";chrom[1]="chrom";
  pos[0]="position";pos[1]="pos";
  while(getline(inf,line)){ //load data into a map to ensure sorted data by x_axis
     int i=0;bool found=false;temp=line;line=toLowerCase(line);
     for(int i=0;i<2;++i){
         if(line.find(chrom[i])!=std::string::npos){found=true;}}
     if(found){ found=false; for(int i=0;i<2;++i){if(line.find(pos[i])!=std::string::npos){found=true;}}}
     if(!found){temp="";}else break;
   }
 return temp;
}
/*********************************************************************************************
 Class Null constructor 
**********************************************************************************************/
/** Loads codon to amino-acid mapping into memory
*/
void ggenome::setaminoAcid(){   //codon to amino-acid mapping  
        stopCodons["UGA"]="*";
        stopCodons["TGA"]="*";
        stopCodons["UAA"]="*";
        stopCodons["TAA"]="*";
        stopCodons["UGA"]="*";
        stopCodons["TGA"]="*";
        startCodon["AUG"]="M";
        startCodon["ATG"]="M";
        // first letter is U/T and second letter is A
       codon2aa["UAG"]="*";
       codon2aa["UAA"]="*";
       codon2aa["UAU"]="Y";
       codon2aa["UAC"]="Y";
        codon2aa["TAT"]="Y";
        codon2aa["TAC"]="Y";     
        codon2aa["TAG"]="*";
        codon2aa["TAA"]="*";
        startCodon;
        // first letter is G and second letter is C
        codon2aa["GCT"]="A";
        codon2aa["GCC"]="A";
        codon2aa["GCA"]="A";
        codon2aa["GCG"]="A";
       codon2aa["GCU"]="A";
       codon2aa["GCC"]="A";
       codon2aa["GCA"]="A";
       codon2aa["GCG"]="A";
        // first letter is G and second letter is A
        codon2aa["GAA"]="E";
        codon2aa["GAG"]="E";
        codon2aa["GAT"]="D";
        codon2aa["GAC"]="D";
       codon2aa["GAA"]="E";
       codon2aa["GAG"]="E";
       codon2aa["GAU"]="D";
       codon2aa["GAC"]="D";
        // first letter is G and second letter is G
        codon2aa["GGT"]="G";
        codon2aa["GGC"]="G";
        codon2aa["GGA"]="G";
        codon2aa["GGG"]="G";
       codon2aa["GGU"]="G";
       codon2aa["GGC"]="G";
       codon2aa["GGA"]="G";
       codon2aa["GGG"]="G";
       // first letter is G and second letter is T/U
        codon2aa["GTT"]="V";
        codon2aa["GTC"]="V";
        codon2aa["GTA"]="V";
        codon2aa["GTG"]="V";
       codon2aa["GUU"]="V";
       codon2aa["GUC"]="V";
       codon2aa["GUA"]="V";
       codon2aa["GUG"]="V";
        // first letter is C and second letter is G
        codon2aa["CGT"]="R";
        codon2aa["CGC"]="R";
        codon2aa["CGA"]="R";
        codon2aa["CGG"]="R";
       codon2aa["CGU"]="R";
       codon2aa["CGC"]="R";
       codon2aa["CGA"]="R";
       codon2aa["CGG"]="R";
        // first letter is C and second letter is C
       codon2aa["CCU"]="P";
       codon2aa["CCC"]="P";
       codon2aa["CCA"]="P";
       codon2aa["CCG"]="P";
        codon2aa["CCT"]="P";
        codon2aa["CCC"]="P";
        codon2aa["CCA"]="P";
        codon2aa["CCG"]="P";
        // first letter is C and second letter is U
       codon2aa["CUU"]="L";
       codon2aa["CUC"]="L";
       codon2aa["CUA"]="L";
       codon2aa["CUG"]="L";
        codon2aa["CTT"]="L";
        codon2aa["CTC"]="L";
        codon2aa["CTA"]="L";
        codon2aa["CTG"]="L";
         // first letter is C and second letter is A
       codon2aa["CAA"]="Q";
       codon2aa["CAG"]="Q";
       codon2aa["CAU"]="H";
       codon2aa["CAC"]="H";
        codon2aa["CAA"]="Q";
        codon2aa["CAG"]="Q";
        codon2aa["CAT"]="H";
        codon2aa["CAC"]="H";
        // first letter is A and second letter is A
       codon2aa["AAU"]="N";
       codon2aa["AAC"]="N";
       codon2aa["AAA"]="K";
       codon2aa["AAG"]="K";
        codon2aa["AAA"]="K";
        codon2aa["AAG"]="K";
        codon2aa["AAT"]="N";
        codon2aa["AAC"]="N";
       // first letter is A and second letter is G
       codon2aa["AGU"]="S";
       codon2aa["AGC"]="S";
        codon2aa["AGT"]="S";
        codon2aa["AGC"]="S";
       codon2aa["AGA"]="R";
       codon2aa["AGG"]="R";
        codon2aa["AGA"]="R";
        codon2aa["AGG"]="R";
        // first letter is A and second letter is U/T
       codon2aa["AUU"]="I";
       codon2aa["AUC"]="I";
       codon2aa["AUA"]="I";
       codon2aa["AUG"]="M";
        codon2aa["ATT"]="I";
        codon2aa["ATC"]="I";
        codon2aa["ATA"]="I";
        codon2aa["ATG"]="M";
        // first letter is A and second letter is C
        codon2aa["ACT"]="T";
        codon2aa["ACC"]="T";
        codon2aa["ACA"]="T";
        codon2aa["ACG"]="T";
       codon2aa["ACU"]="T";
       codon2aa["ACC"]="T";
       codon2aa["ACA"]="T";
       codon2aa["ACG"]="T";
        // first letter is U/T and second letter is G
       codon2aa["UGU"]="C";
       codon2aa["UGC"]="C";
       codon2aa["UGG"]="W";
       codon2aa["UGA"]="*";
        codon2aa["TGG"]="W";
        codon2aa["TGT"]="C";
        codon2aa["TGC"]="C";
        codon2aa["TGA"]="*";
        // first letter is U/T and second letter is A
       codon2aa["UAG"]="*";
       codon2aa["UAA"]="*";
       codon2aa["UAU"]="Y";
       codon2aa["UAC"]="Y";
        codon2aa["TAT"]="Y";
        codon2aa["TAC"]="Y";     
        codon2aa["TAG"]="*";
        codon2aa["TAA"]="*";
         // first letter is U/T and second letter is U/T
       codon2aa["UUA"]="L";
       codon2aa["UUG"]="L";
       codon2aa["UUU"]="F";
       codon2aa["UUC"]="F";
        codon2aa["TTA"]="L";
        codon2aa["TTG"]="L";
        codon2aa["TTT"]="F";
        codon2aa["TTC"]="F";
        // first letter is U/T and second letter is C
       codon2aa["UCU"]="S";
       codon2aa["UCC"]="S";
       codon2aa["UCA"]="S";
       codon2aa["UCG"]="S";
        codon2aa["TCT"]="S";
        codon2aa["TCC"]="S";
        codon2aa["TCA"]="S";
        codon2aa["TCG"]="S";
        oneChar2Symbol['A']="Ala";oneChar2Symbol['C']="Cys";oneChar2Symbol['D']="Asp";oneChar2Symbol['E']="Glu";
        oneChar2Symbol['F']="Phe";oneChar2Symbol['G']="Gly";oneChar2Symbol['H']="His";oneChar2Symbol['I']="Ile";
        oneChar2Symbol['K']="Lys";oneChar2Symbol['L']="Leu";oneChar2Symbol['M']="Met";oneChar2Symbol['N']="Asn";
        oneChar2Symbol['P']="Pro";oneChar2Symbol['Q']="Gln";oneChar2Symbol['R']="Arg";oneChar2Symbol['S']="Ser";
        oneChar2Symbol['T']="Thr";oneChar2Symbol['V']="Val";oneChar2Symbol['W']="Trp";oneChar2Symbol['Y']="Tyr";
        oneChar2Symbol['*']="Stop";

        oneLetter2Name['A']="Alanine";oneLetter2Name['C']="Cysteine";oneLetter2Name['D']="Aspartic acid";
        oneLetter2Name['E']="Glutamic acid";oneLetter2Name['F']="Phenylalanine";oneLetter2Name['G']="Glycine";
        oneLetter2Name['H']="Histidine";oneLetter2Name['I']="Isoleucine";oneLetter2Name['K']="Lysine";
        oneLetter2Name['L']="Leucine";oneLetter2Name['M']="Methionine";
        oneLetter2Name['N']="Asparagine";oneLetter2Name['P']="Proline";oneLetter2Name['Q']="Glutamine";
        oneLetter2Name['R']="Arginine";oneLetter2Name['S']="Serine";oneLetter2Name['T']="Threonine";
        oneLetter2Name['V']="Valine";oneLetter2Name['W']="Tryptophan";oneLetter2Name['Y']="Tyrosine";
        oneLetter2Name['*']="Stop";
        symbol2OneChar["ala"]="A";symbol2OneChar["cys"]="C";symbol2OneChar["asp"]="D";symbol2OneChar["glu"]="E";
        symbol2OneChar["phe"]="F";symbol2OneChar["gly"]="G";symbol2OneChar["his"]="H";symbol2OneChar["ile"]="I";
        symbol2OneChar["lys"]="K";symbol2OneChar["leu"]="L";symbol2OneChar["met"]="M";symbol2OneChar["asn"]="N";
        symbol2OneChar["pro"]="P";symbol2OneChar["gln"]="Q";symbol2OneChar["arg"]="R";symbol2OneChar["ser"]="S";
        symbol2OneChar["thr"]="T";symbol2OneChar["val"]="V";symbol2OneChar["trp"]="W";symbol2OneChar["tyr"]="Y";
         symbol2OneChar["stop"]="*";
        name2OneChar["alanine"]="A";name2OneChar["cysteine"]="C";name2OneChar["aspartic acid"]="D";
        name2OneChar["glutamic acid"]="E";name2OneChar["phenylalanine"]="F";name2OneChar["glycine"]="G";
        name2OneChar["histidine"]="H";name2OneChar["isoleucine"]="I";name2OneChar["lysine"]="K";
        name2OneChar["leucine"]="L";name2OneChar["methionine"]="M";name2OneChar["asparagine"]="N";
        name2OneChar["proline"]="P";name2OneChar["glutamine"]="Q";name2OneChar["arginine"]="R";
        name2OneChar["serine"]="S";name2OneChar["threonine"]="T";name2OneChar["valine"]="V";
        name2OneChar["tryptophan"]="W";name2OneChar["tyrosine"]="Y";
        name2OneChar["stop"]="stop";
        //D, and E are electrically charged negavively. K,R, and H are positiveley charged positively
        acidity['A']="Neutral";acidity['C']="Neutral";acidity['D']="Acidic";acidity['E']="Acidic";
        acidity['F']="Neutral";acidity['G']="Neutral";acidity['H']="Weakly basic";acidity['I']="Neutral";
        acidity['K']="Basic";acidity['L']="Neutral";acidity['M']="Neutral";acidity['N']="Neutral";
        acidity['P']="Neutral";acidity['Q']="Neutral";acidity['R']="Strongly basic";acidity['S']="Neutral";
        acidity['T']="Neutral";acidity['V']="Neutral";acidity['W']="Neutral";acidity['Y']="Neutral";

        polarity['A']="Non-polar";polarity['C']="Slightly-polar";polarity['D']="Polar";polarity['E']="Polar";
        polarity['F']="Non-polar";polarity['G']="Non-polar";polarity['H']="Polar";polarity['I']="Non-polar";
        polarity['K']="Polar";polarity['L']="Non-polar";polarity['M']="Non-polar";polarity['N']="Polar";
        polarity['P']="Non-polar";polarity['Q']="Polar";polarity['R']="Polar";polarity['S']="Polar";
        polarity['T']="Polar";polarity['V']="Non-polar";polarity['W']="Slightly-polar";polarity['Y']="Polar";
       /**
        Hydropathy (hydrophobicity vs. hydrophilicity or lipophobicity vs. lipophilicity)
        is usually characterized by numbers (hydrophobic moments, HM) from -7.5 (Arg) to 3.1 (Ile), 
       whereas hydrophobicity is a measure of how strongly the side chains are pushed out of water. 
       The more positive a number, the more the amino acid residue will tend not to be in an aqueous environment. 
       Negative numbers indicate hydrophilic side chains, with more negative numbers indicating greater affinity for water.

        Molecules with similar hydropathy have affinity to each other, they are compatible; 
       molecules with different hydropathy repel each other, and they are not compatible. 
       To express this numerically, we use the hydropathy compatibility index (HCI) and collect 
       these indices (20 × 20) in the matrix. HCIs were calculated using the formula 
       HCI = 20 - abs( [HM(A) - HM(B)] × 19/10.6])
       where HM(A) and HM(B) are the hydrophobic moments of the amino acids A and B 
       and HM(Arg)-HM (Ile) = 10.6. This formula gives the maximal index (20) for 
      identical amino acids (closest hydrophobicity) and the minimal value (1) 
       for the two hydrophobically most distant amino acids (Arg and Ile). The "|" indicate absolute values
       
       **/
        //hydrophobic
        hydropathyIndex['A']=1.8;hydropathyIndex['C']=2.5;hydropathyIndex['M']=1.9;
        hydropathyIndex['F']=2.8;hydropathyIndex['I']=4.5;hydropathyIndex['L']=3.8;hydropathyIndex['V']=4.2; 
        hydropathyIndex['G']=-0.4;
        //hydrophilic
        hydropathyIndex['H']=-3.2;hydropathyIndex['R']=-4.5;hydropathyIndex['K']=-3.9;//charged positively
        hydropathyIndex['E']=-3.5;hydropathyIndex['D']=-3.5;  //charged negatively
        hydropathyIndex['N']=-3.5;hydropathyIndex['P']=-1.6;hydropathyIndex['Q']=-3.5;
        hydropathyIndex['S']=-0.8;hydropathyIndex['T']=-0.7;hydropathyIndex['W']=-0.9;hydropathyIndex['Y']=-1.3;
       
       /// amino acid Molecular Weight (g/mol) 
       ///Approximate Molecular Weight of a Protein: M.W. of protein = # amino acids x 110 Da 
        mass['A']=89.1;mass['C']=121.2;mass['D']=133.1;mass['E']=147.1;
        mass['F']=165.2;mass['G']=75.1;mass['H']=155.2;mass['I']=131.2;
        mass['K']=146.2;mass['L']=131.2;mass['M']=149.2;mass['N']=132.1;
        mass['P']=115.1;mass['Q']=146.2;mass['R']=174.2;mass['S']=105.1;
        mass['T']=119.1;mass['V']=117.1;mass['W']=204.2;mass['Y']=181.1;
        
        ///volume of the amino acid in Ångstrom* (ang^3)
        volume['A']=88.6;volume['C']=108.5;volume['D']=111.1;volume['E']=138.4;
        volume['F']=189.9;volume['G']=60.1;volume['H']=153.2;volume['I']=166.7;

        volume['K']=168.6;volume['L']=166.7;volume['M']=162.9;volume['N']=114.1;
        volume['P']=112.7;volume['Q']=143.8;volume['R']=173.4;volume['S']=89.0;
        volume['T']=116.1;volume['V']=140.0;volume['W']=227.8;volume['Y']=193.6;
}
/** isValidCodon returns true if the current codon is valid  or false otherwise *********/
bool aminoAcid::isValidCodon(string codon){
   codon= toUpperCase(trimStr(codon));
   if(codon2aa.find(codon)!=codon2aa.end())return true;
   else return false;
}
/*********************** displays the aminoacid table *****************************/
void aminoAcid::displayTable(){
   strMap::iterator iter = symbol2OneChar.begin();string codon="Asn"; 
  cout<<"Name\tthreeLetterSymbol\toneLetterSymbol\tAcidity\tPolarity\thydroPathyIndex\tMass (g/mole)\tVolume (ang^3)\tCodons\n";
   while(iter != symbol2OneChar.end()){
      string codonlist=""; strVec::iterator viter; strVec codonVec=getCodonList(iter->first,0);
      if(!codonVec.empty()){
        for(viter=codonVec.begin();viter!=codonVec.end();++viter)
        {codonlist+=(*viter)+",";}
      }
      cout<<getAaName(iter->first)<<"\t"<<getThreeLetterSymbol(iter->first)<<"\t"<<iter->second<<"\t";
      cout<<getAaAcidity(iter->first)<<"\t"<<getAaPolarity(iter->first)<<"\t";
      cout<<getAaHydropathyIndex(iter->first)<<"\t"<<getAaMass(iter->first)<<"\t";
      cout<<getAaVolume(iter->first)<<"\t"<<codonlist<<endl;++iter;
    }
}
/**************************************************************************************
returns a list of three-character codon for the specified aa three-letter symbol
   or aminoacid comon name. if isDNA I true, return a DNA 
sequence otherwise returns an RNA sequence
*************************************************************************************/
strVec aminoAcid::getCodonList(string aa, bool isDNA){
 strVec codons;  strMap::iterator iterStart=codon2aa.begin(), iterEnd=codon2aa.end();
 string token=getOneLetterSymbol(aa);
 if(!token.empty()){
    while(iterStart!=iterEnd){
        if(iterStart->second ==token)codons.push_back(iterStart->first);++iterStart;
    }
 }return codons;
}
/********************************************************************
this function returns a one character 
aminoacid code for the specified codon/commonName
**********************************************************************/
string aminoAcid::getOneLetterSymbol(string codon){ string token="";
  if(isValidCodon(codon)){
     if(codon2aa.find(toUpperCase(trimStr(codon)))!=codon2aa.end())
        token= codon2aa[toUpperCase(trimStr(codon))];
   }
   else{
      if(symbol2OneChar.find(toLowerCase(trimStr(codon)))!=symbol2OneChar.end())
         token= symbol2OneChar[toLowerCase(trimStr(codon))];
      else if(name2OneChar.find(toLowerCase(trimStr(codon)))!=name2OneChar.end())
         token=name2OneChar[toLowerCase(trimStr(codon))];
   }
 return token;
}

/***********************************************************************************************************
   private member:
   setAnnotatedOrg generates the current list of all the organisms (common name,version)
   with annotations using the specified host, default Graber Transcript database
************************************************************************************************************/
void ggenome::setAnnotatedOrg(){    
string url= genomeData.webserviceUrl+"?s=1"; if(!genomeData.host.empty())url+="&host="+genomeData.host;
string query="use XML::Simple;use Data::Dumper;use LWP 5.64;my $browser = LWP::UserAgent->new;\
              $browser->timeout(10);$browser->env_proxy;$browser->agent('Mozilla/5.0');\
              $xml = new XML::Simple (KeyAttr=>[]);$url=\""+url+"\";\
              my $xml_file=\""+genomeData.annotationsXmlFile+"\"; $out_file=\""+genomeData.annotationsTexFile+"\";\
              my $response = $browser->get($url, \":content_file\" => $xml_file);\
              open(OUT,\"> $out_file\");if(OUT){ print OUT \"Organism\\tVersion\\tOrganism_id\\tVersion_id\\n\";\
              $data = $xml->XMLin(\"$xml_file\");\
              if(ref($data->{organism}) ne 'ARRAY'){ eval{\
                 $organism=$data->{organism}; $organism_name=$organism->{name};\
                 $organism_version=$organism->{version};$id=$organism->{organism_id};$version_id=$organism->{version_id};\
                 print OUT \"$organism_name\\t$organism_version\\t$id\\t$version_id\\n\"; };\
               }else{\
                 foreach $organism(@{$data->{organism}}){ eval{\
                    $organism_name=$organism->{name};$organism_version=$organism->{version};\
                    $id=$organism->{organism_id};$version_id=$organism->{version_id};\
                    print OUT \"$organism_name\\t$organism_version\\t$id\\t$version_id\\n\";};\
                 }}close(OUT);if(-f \"$xml_file\"){my $unlink= unlink \"$xml_file\";}}";//
      long lSize=double(query.length()+1)+100;char *temp = (char*) malloc (sizeof(char)*(lSize+1));
      if (temp != NULL){sprintf(temp,"perl -e '%s'",query.c_str());system(temp);}
      free (temp);
}
/***********************************************************************************************************
   private member:
   setPredictionOrg generates the current list of all the gene predictions
   for a given organism version 
************************************************************************************************************/
void genomeFeature::setOrgPredictions(string organismVersion){
string predictions_xml_file=genomeData.dataDirName+"/"+organismVersion+"-genePredictions.xml";
string predictions_text_file=genomeData.dataDirName+"/"+organismVersion+"-genePredictions.txt";
if(!genomeData.host.empty()){
  predictions_xml_file=genomeData.dataDirName+"/"+genomeData.host+"_"+organismVersion+"-genePredictions.xml";
  predictions_text_file=genomeData.dataDirName+"/"+genomeData.host+"_"+organismVersion+"-genePredictions.txt";
}
struct stat sb;
if(stat(predictions_text_file.c_str(),&sb)==-1){
  string url= genomeData.webserviceUrl+"?l=1&v="+organismVersion;
  if(!genomeData.host.empty())url+="&host="+genomeData.host;
  string query="use XML::Simple;use Data::Dumper;use LWP 5.64;my $browser = LWP::UserAgent->new;\
              $browser->timeout(10);$browser->env_proxy;$browser->agent('Mozilla/5.0');\
              $xml = new XML::Simple (KeyAttr=>[]);$url=\""+url+"\";\
              my $xml_file=\""+predictions_xml_file+"\"; $out_file=\""+predictions_text_file+"\";\
              my $response = $browser->get($url, \":content_file\" => $xml_file);open(OUT,\"> $out_file\");\
              if(OUT){ print OUT \"prediction_name\\tprediction_id\\n\";\
                   $data = $xml->XMLin(\"$xml_file\");\
                   if(ref($data->{prediction}) ne 'ARRAY'){eval{\
                     $prediction=$data->{prediction}; $prediction_name=$prediction->{name};chomp($prediction_name);\
                     $id=$prediction->{id};print OUT \"$prediction_name\\t$id\\n\";\
                   };}else{foreach $prediction(@{$data->{prediction}}){eval{\
                      $prediction_name=$prediction->{name};chomp($prediction_name);$id=$prediction->{id};\
                      print OUT \"$prediction_name\\t$id\\n\";};\
                    }}\
              }close(OUT);if(-f \"$xml_file\"){my $unlink= unlink \"$xml_file\";}";
      long lSize=double(query.length()+1)+100;
      char *temp = (char*) malloc (sizeof(char)*(lSize+1));
      if (temp != NULL){ 
             sprintf(temp,"perl -e '%s'",query.c_str());system(temp);
       }
    free (temp);
 }//download feature config
 predictions_xml_file=genomeData.dataDirName+"/ucsc_"+organismVersion+"-features.xml";
 predictions_text_file=genomeData.dataDirName+"/ucsc_"+organismVersion+"-features.txt";
 if(stat(predictions_text_file.c_str(),&sb)==-1){
    string url= genomeData.webserviceUrl+"?featureList=1&host=ucsc&v="+organismVersion;
    string query="use XML::Simple;use Data::Dumper;use LWP 5.64;my $browser = LWP::UserAgent->new;\
              $browser->timeout(10);$browser->env_proxy;$browser->agent('Mozilla/5.0');\
              $xml = new XML::Simple (KeyAttr=>[]);$url=\""+url+"\";\
              my $xml_file=\""+predictions_xml_file+"\"; $out_file=\""+predictions_text_file+"\";\
              my $response = $browser->get($url, \":content_file\" => $xml_file);open(OUT,\"> $out_file\");\
              if(OUT){ print OUT \"prediction_name\\tprediction_id\\n\";\
                   $data = $xml->XMLin(\"$xml_file\");\
                   if(ref($data->{prediction}) ne 'ARRAY'){eval{\
                     $prediction=$data->{prediction}; $prediction_name=$prediction->{name};\
                     $id=$prediction->{id};chomp($prediction_name);print OUT \"$prediction_name\\t$id\\n\";\
                   };}else{foreach $prediction(@{$data->{prediction}}){eval{\
                      $prediction_name=$prediction->{name};chomp($prediction_name);$id=$prediction->{id};\
                      print OUT \"$prediction_name\\t$id\\n\";};\
                    }}\
              }close(OUT);if(-f \"$xml_file\"){my $unlink= unlink \"$xml_file\";}";
       long lSize=double(query.length()+1)+100;
       char *temp = (char*) malloc (sizeof(char)*(lSize+1));
       if (temp != NULL){sprintf(temp,"perl -e '%s'",query.c_str());system(temp);}
       free (temp);
   }
}
/***********************************************************************************************************
   private member:
   setOrgChromosomes generates the list of all the chromosomes
   for a given organism version 
************************************************************************************************************/
void genomeFeature::setOrgChromosomes(string organismVersion){
 string chromosomes_xml_file=genomeData.dataDirName+"/"+organismVersion+"-chromosomes.xml";
 string chromosomes_text_file=genomeData.dataDirName+"/"+organismVersion+"-chromosomes.txt";
 string url= genomeData.webserviceUrl+"?v="+organismVersion+"&chrl=1";
 string query="use XML::Simple;use Data::Dumper;use LWP 5.64;my $browser = LWP::UserAgent->new;\
              $browser->timeout(10);$browser->env_proxy;$browser->agent('Mozilla/5.0');\
              $xml = new XML::Simple (KeyAttr=>[]);$url=\""+url+"\";\
              my $xml_file=\""+chromosomes_xml_file+"\"; $out_file=\""+chromosomes_text_file+"\";\
              my $response = $browser->get($url, \":content_file\" => $xml_file);\
              open(OUT,\"> $out_file\");if(OUT){ print OUT \"Chromosome\\tchromID\\tchromSize\\n\";\
              $data = $xml->XMLin(\"$xml_file\");\
              if(ref($data->{chromosome}) ne 'ARRAY'){ eval{\
                 $chrom=$data->{chromosome}; $name=$chrom->{name};\
                 $size=$chrom->{size};$id=$chrom->{id};print OUT \"$name\\t$id\\t$size\\n\"; };\
               }else{\
                 foreach $chrom(@{$data->{chromosome}}){ eval{\
                    $name=$chrom->{name};\
                    $size=$chrom->{size};$id=$chrom->{id};print OUT \"$name\\t$id\\t$size\\n\"; };\
                 }}close(OUT);if(-f \"$xml_file\"){my $unlink= unlink \"$xml_file\";}}";
      long lSize=double(query.length()+1)+100;char *temp = (char*) malloc (sizeof(char)*(lSize+1));
      if (temp != NULL){ sprintf(temp,"perl -e '%s'",query.c_str());system(temp);}
      free (temp);
}
/********************************************************* in used
loads the list of organisms that have annotations 
in the database into a dictionary
Rule of thumb: make all the look up keys to be lower cases
**********************************************************/
void ggenome::loadAnnotatedOrganisms(){
 ifstream infile(genomeData.annotationsTexFile.c_str());const char* delimiter=",\t";
 if(!infile){setAnnotatedOrg();}if(infile)infile.close();
 infile.open(genomeData.annotationsTexFile.c_str());
 if(infile){
     string line;string organism,build;
     getline(infile,line);
     if(!organismInDB.empty()) organismInDB.clear();
     if(!organismIdMap.empty())organismIdMap.clear();
     while(getline(infile,line))
      {
          if(!line.empty())
          {
               long lSize=0;char *buffer2;lSize=long(fabs(double(line.length())));
               buffer2 =(char*)malloc(sizeof(char)*(lSize+1));int i=0,id=0,version_id=0;
               strcpy(buffer2,line.c_str());char * pch;pch= strtok(buffer2,delimiter);
               while (pch != NULL){
                      if(i==0)organism=string(pch);
                      else if(i==1)build=string(pch);
                      else if(i==2)id=atoi(pch); else if(i==3)version_id=atoi(pch);
                      pch = strtok (NULL, delimiter);++i;
               }
               organismInDB[toLowerCase(build)].organism=organism; 
               organismInDB[toLowerCase(build)].organism_id=id;
               organismInDB[toLowerCase(build)].organismVersion_id=version_id;
               organismIdMap[toLowerCase(organism)]=id;
         }
     }
 }
}
/******************************************************* in used
downloads annotations from database server 
Note: this will be cached on the server (web-service server) then wget the file
     because I don't want to rely on the client perl modules configuration
The  cached annotations files are updated every 30 days
A call to the server: service?host=ucsc&v=mm9&a=ensGene&annot=1
***********************************************************************/
void genomeFeature::loadAnnotationTranscripts(string organismVersion,string annotSource, string& annotations_text_file){
  annotations_text_file=genomeData.dataDirName+"/"+organismVersion+"-"+annotSource+"-transcripts.txt";
  string annot_zipfile=genomeData.dataDirName+"/"+organismVersion+"-"+annotSource+"-transcripts.txt.gz";
  string url= genomeData.webserviceUrl+"?v="+organismVersion+"&annot=1&host="+genomeData.host+"&a="+annotSource;
  if((annotSource.find("miRNA")!=std::string::npos)||(annotSource.find("snp")!=std::string::npos)||
      (annotSource.find("microsat")!=std::string::npos)||(annotSource.find("simpleRepeat")!=std::string::npos)){
      annotations_text_file=genomeData.dataDirName+"/ucsc_"+organismVersion+"-"+annotSource+".txt";
      annot_zipfile=genomeData.dataDirName+"/ucsc_"+organismVersion+"-"+annotSource+".txt.gz";
      url= genomeData.webserviceUrl+"?v="+organismVersion+"&annot=1&host=ucsc&a="+annotSource;
  }
 
  if(!genomeData.host.empty()){
      annotations_text_file=genomeData.dataDirName+"/"+genomeData.host+"_"+organismVersion+"-"+annotSource+"-transcripts.txt";
      annot_zipfile=genomeData.dataDirName+"/"+genomeData.host+"_"+organismVersion+"-"+annotSource+"-transcripts.txt.gz";
     if((annotSource.find("miRNA")!=std::string::npos)||(annotSource.find("snp")!=std::string::npos)||
        (annotSource.find("microsat")!=std::string::npos)||(annotSource.find("simpleRepeat")!=std::string::npos)){
         annotations_text_file=genomeData.dataDirName+"/ucsc_"+organismVersion+"-"+annotSource+".txt";
         annot_zipfile=genomeData.dataDirName+"/ucsc_"+organismVersion+"-"+annotSource+".txt.gz";
         url= genomeData.webserviceUrl+"?v="+organismVersion+"&annot=1&host=ucsc&a="+annotSource;
     }
  }
  struct stat sb;
  if((stat(annotations_text_file.c_str(), &sb)==-1)||sb.st_size<=200||(lastModifiedDayCount(annotations_text_file))>30 ){
   string nurl="temp.txt"; string query="";
   if(stat(annot_zipfile.c_str(),&sb)!=-1){ //remove existing zip file before downloading a new one
      query="rm "+annot_zipfile;
      long lSize=double(query.length()+1)+100;char *temp = (char*) malloc (sizeof(char)*(lSize+1));
      sprintf(temp,"%s",query.c_str());system(temp);free (temp);
    }
   query="use XML::Simple;use LWP 5.64;my $browser = LWP::UserAgent->new;\
              $browser->timeout(10);$browser->env_proxy;$browser->agent('Mozilla/5.0');\
              $xml = new XML::Simple (KeyAttr=>[]);$url=\""+url+"\";\
              my $xml_file=temp.xml; $out_file=\"temp.txt\";\
              my $response = $browser->get($url, \":content_file\" => $xml_file);\
              open(OUT,\"> $out_file\");if(OUT){\
              $data = $xml->XMLin(\"$xml_file\");\
              if(ref($data->{annotfile}) ne 'ARRAY'){ eval{\
                 $chrom=$data->{annotfile}; $name=$chrom->{name};\
                 print OUT \"$name\\n\"; };\
               }else{\
                 foreach $chrom(@{$data->{annotfile}}){ eval{\
                    $name=$chrom->{name};\
                    print OUT \"$name\\n\"; };\
                 }}close(OUT);if(-f \"$xml_file\"){my $unlink= unlink \"$xml_file\";}}";
    //I need to replace this query with wget -O 'url' 
    long lSize=double(query.length()+1)+100;char *temp = (char*) malloc (sizeof(char)*(lSize+1));
    if (temp != NULL){ sprintf(temp,"perl -e '%s'",query.c_str());system(temp);}free (temp);
    if(stat(nurl.c_str(), &sb)!=-1){ //get file path from temp.txt
       ifstream inf(nurl.c_str()); string annoturl; getline(inf,annoturl);//cout<<"The url is :"<<annoturl<<endl;
       if(stat(annotations_text_file.c_str(), &sb)!=-1){
          query="mv "+annotations_text_file+" old_"+annotations_text_file;
          lSize=double(query.length()+1)+100;char *temp2 = (char*) malloc (sizeof(char)*(lSize+1));
          sprintf(temp2,"%s",query.c_str());system(temp2);free (temp2);
        }//now download the annoattions
        query="wget -O "+annot_zipfile+" '"+annoturl+"' -q";
        lSize=double(query.length()+1)+100;char *temp2 = (char*) malloc (sizeof(char)*(lSize+1));
        sprintf(temp2,"%s",query.c_str());system(temp2);free (temp2);
        if(stat(annot_zipfile.c_str(), &sb)!=-1){
           query="gunzip -q "+annot_zipfile; 
          lSize=double(query.length()+1)+100;char *temp3 = (char*) malloc (sizeof(char)*(lSize+1));
          sprintf(temp3,"%s",query.c_str());system(temp3);free (temp3);
        }
   }}
}
         

 /********************************************************************************* in used
  Index gene annotations into memory using genomic region uniq id.
  Generates an id for every uniq genomic region [chr,start,end]and stores into a map
  called transcriptIds [pair(chrom,start)][end]=id. 
  These ids will be used to build transcript networks called transcripts <id,listTx>.
  transcriptIds is used to get overlapping transcripts for a given SNP.
  The transcriptIds data structures speed up key lookup using the map lower and upper bound functions on pair(chrom,start)
  .Transcripts is use to get the list of transcripts for a given txid

Note I need to use the same structure I used for connected components
exons[std::make_pair(std::make_pair(fields[0],atoi(fields[1])),"start")].push_back(std::make_pair(atoi(fields[2]),count));
exons[std::make_pair(std::make_pair(fields[0],atoi(fields[2])),"end")].push_back(std::make_pair(atoi(fields[1]),count));

transcriptIds[std::make_pair(std::make_pair(bG.chrom,bG.chromEnd),"end")][std::make_pair(bG.strand,bG.chromStart)]=id;
transcriptIds[std::make_pair(std::make_pair(bG.chrom,bG.chromStart),"start")][std::make_pair(bG.strand,bG.chromEnd)]=id;


typedef std::pair<string,int> strIntMap;
typedef std::pair<strIntMap,string>charIntcharPair;
typedef std::vector<intPair> nodeIndex;          //stores <exonEnd,index> pair
typedef std::map<charIntcharPair,nodeIndex> exonIndex2Map; //stores chr,v1_end,type,v2_start,strand
**********************************************************************************************************************/
void genomeFeature::IndexAnnotationsByRegion(featureIndexMap& transcriptIds,featureList& transcripts,
                ifstream& inf,const char* delimiter,genomeFile &fileObject){ 
    string headerLine,line; int id=0; tabularData bG;
    if(!transcriptIds.empty())transcriptIds.clear();if(!transcripts.empty())transcripts.clear(); 
    getline(inf,headerLine);featureIndexMap::iterator featureIter;longIntMaps::iterator listIter;
    while(getline(inf,line)){ 
          fileObject.getTabularFields((char*)headerLine.c_str(),(char*)line.c_str(),delimiter,bG);
          if((bG.chrom=="-1")||(bG.chromStart<=0)||(bG.chromEnd<=0)) continue; ++id;
          bG.chrom=filterChrom(bG.chrom);
          featureIter=transcriptIds.find(std::make_pair(bG.chrom,bG.chromStart));
          if(featureIter==transcriptIds.end()){ //new genomic region=new id
             transcriptIds[std::make_pair(bG.chrom,bG.chromStart)][bG.chromEnd]=id;transcripts[id].push_back(bG);
          }else{ //existing genomic region, get this id 
             listIter=transcriptIds[std::make_pair(bG.chrom,bG.chromStart)].find(bG.chromEnd);
             if(listIter==transcriptIds[std::make_pair(bG.chrom,bG.chromStart)].end()){
                transcriptIds[std::make_pair(bG.chrom,bG.chromStart)][bG.chromEnd]=id;transcripts[id].push_back(bG);
             }else{ //transcript may share the same genomic region id but could have different set of exons
                int thisId=transcriptIds[std::make_pair(bG.chrom,bG.chromStart)][bG.chromEnd];
                transcripts[thisId].push_back(bG);
             }
          }
     }
}

/***************************************************************** in used
 Returns the transcript with nearest exon from specified bp position
 transcriptIds[std::make_pair(bG.chrom,bG.chromStart)][bG.chromEnd]=id
******************************************************************/
int genomeFeature::getNearestTranscript(string chrom, long pos,featureIndexMap& transcriptIds){
   featureIndexMap::iterator txUpperIter,txLowerIter;longIntMaps::iterator subiter;int leftregionid=0,rightregionid=0;
 long leftEnds=0,rightTxStart=0;string strand,thischrom;
 txUpperIter=transcriptIds.upper_bound(std::make_pair(chrom,pos));txLowerIter=transcriptIds.begin();
 //get coordinate and id of first tx right of pos
 if(txUpperIter!=transcriptIds.end()){rightTxStart=txUpperIter->first.second;rightregionid=txUpperIter->second.begin()->second;
   thischrom=txUpperIter->first.first;--txUpperIter;} bool more=true;
 while(more){ 
      for(subiter=txUpperIter->second.begin();subiter!=txUpperIter->second.end();++subiter){
          if(subiter->first<pos){//find tx on 5' end of SNP whose txend is the closest to pos
             if(leftEnds==0){leftEnds=subiter->first;leftregionid=subiter->second;}
             else if(leftEnds<subiter->first){leftEnds=subiter->first;leftregionid=subiter->second;}
       }} --txUpperIter;
     try{ if(txUpperIter->first.first==chrom)more=true;else{more=false;}}
     catch(exception& e){more=false;} 
  }
  if((pos-leftEnds)<(rightTxStart-pos)){return leftregionid;}else{return rightregionid;}
}
/***************************************************************** ******* in used
 Returns a vector of transcriptIds that overlap the specified bp position 
*****************************************************************************/
void genomeFeature::getOverlapTranscripts(string chrom, long pos,featureIndexMap& transcriptIds,intVec& txList){
 featureIndexMap::iterator txUpperIter,txLowerIter;longIntMaps::iterator subiter;
 //get a pointer to the first location where chromStart>= snp.pos ; transcripts are stored sorted by chromStart
 txUpperIter=transcriptIds.upper_bound(std::make_pair(chrom,pos));txLowerIter=transcriptIds.begin();
 if(!txList.empty())txList.clear(); if(txUpperIter!=transcriptIds.end())txUpperIter++; bool more=true;
 while(more){ 
      for(subiter=txUpperIter->second.begin();subiter!=txUpperIter->second.end();++subiter){
          if((txUpperIter->first.second<=pos)&&(subiter->first>=pos)){txList.push_back(subiter->second);}}
    --txUpperIter;
     try{ if(txUpperIter->first.first==chrom)more=true;else{more=false;}}
     catch(exception& e){more=false;} 
  }
}
/********************************************************************************** 
  Index gene annotations into memory using transcript accession id and gene name.
  This speed up term lookup
***********************************************************************************/
void genomeFeature::IndexAnnotationsByAcc(accessionFeatureList& genes,accessionFeatureList& transcripts,
           ifstream& inf,const char* delimiter,genomeFile &fileObject){
   string headerLine,line; int id=0; tabularData bG;
   if(!genes.empty())genes.clear();if(!transcripts.empty())transcripts.clear(); 
   getline(inf,headerLine);//name,gene
   while(getline(inf,line)){
       fileObject.getTabularFields((char*)headerLine.c_str(),(char*)line.c_str(),delimiter,bG);
       if((bG.chrom=="-1")||(bG.chromStart<=0)||(bG.chromEnd<=0)) continue; ++id;
       if(!bG.name.empty())transcripts[bG.name].push_back(bG);
       if(!bG.gene.empty())genes[bG.gene].push_back(bG);          
   }  
}

/************ returns a local id associated to a given organism version) *******
int genomeFeature::getOrganismVersionId(string organismVersion){
   if(organismInDB.find(toLowerCase(organismVersion))!=organismInDB.end())
     return organismInDB[toLowerCase(organismVersion)].organismVersion_id;
   else return 0;
}

/********************************************************* in used
Displays the list of organisms that have annotations 
in for the selected database
**********************************************************/
void genomeFeature::displayOrganismInDB(){
 if(organismInDB.empty())loadAnnotatedOrganisms();
 organismInDbMap::iterator iter;string db="Graber Transcript database";
 if(!genomeData.host.empty())db=toLowerCase(genomeData.host)+" database";
 cout<<"\n*******\nOrganims With Annotation Data in "<<db<<"\n*******\nOrganism\tVersion\n";
 for(iter=organismInDB.begin();iter!=organismInDB.end();++iter){
    cout<<iter->second.organism<<"\t"<<iter->first<<endl;
  }
}
/********************************************************* in used
Displays the current list of all the gene predictions
for a given organism version 
**********************************************************/
void genomeFeature::displayOrgAnnotList(string organismVersion,const char* delimiter,genomeFile &fileObject){
  string predictions_text_file=genomeData.dataDirName+"/"+organismVersion+"-genePredictions.txt";
  string db="Graber Transcript database";strVec fields;fieldIndexMap fieldsMap;
  if(!genomeData.host.empty()){ db=toLowerCase(genomeData.host)+" database";
     predictions_text_file=genomeData.dataDirName+"/"+genomeData.host+"_"+organismVersion+"-genePredictions.txt";
  }
  struct stat sb; if(stat(predictions_text_file.c_str(), &sb)==-1)setOrgPredictions(organismVersion);
  ifstream inf(predictions_text_file.c_str());string line; //
  cout<<"\n*******\nAvailable Gene Prediction sources for "<<organismVersion<<" in "<<db<<"\n*******\n";
  string header; getline(inf,header); fileObject.getFieldsIndex((char*)header.c_str(),delimiter,fieldsMap);
  fieldIndexMap::iterator iter=fieldsMap.find("prediction_name"); int index=(iter!=fieldsMap.end())?iter->second:-1; 
  if(index>=0)cout<<"Prediction Source Name"<<endl;
  while(getline(inf,line)){ fileObject.getFieldsData((char*)line.c_str(),delimiter,fields);
    if((index>=0)&&(index<fields.size()))cout<<fields[index]<<endl;
  }
}
void genomeFeature::generateFeatureList(string organismVersion,const char* delimiter,genomeFile &fileObject,strMap& featMap){
  setOrgPredictions(organismVersion);//generate feature annotation configs for this organism version
  string predictions_text_file=genomeData.dataDirName+"/ucsc_"+organismVersion+"-features.txt";
  struct stat sb; 
  if(stat(predictions_text_file.c_str(), &sb)!=-1){ ifstream inf(predictions_text_file.c_str());string line;
     fieldIndexMap fieldsMap;strVec fields;
     string header; getline(inf,header); fileObject.getFieldsIndex((char*)header.c_str(),delimiter,fieldsMap);
     fieldIndexMap::iterator iter=fieldsMap.find("prediction_name"); int index=(iter!=fieldsMap.end())?iter->second:-1; 
     //if(index>=0)cout<<"Features List"<<endl;
     while(getline(inf,line)){ fileObject.getFieldsData((char*)line.c_str(),delimiter,fields);
        if((index>=0)&&(index<fields.size())){featMap[fields[index]]=fields[index];/*cout<<fields[index]<<endl;*/}
     }
  }
}


