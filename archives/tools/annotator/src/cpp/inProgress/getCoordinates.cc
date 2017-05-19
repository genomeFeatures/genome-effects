/*************************************************************************
  Project:	Graber utils standard code effort
  Purpose:	generates the coordinates  of transcripts/features for 
                the specified organism version and annotation source.

Note: The output of this program can be used by other programs, for example
      snpAnnotator, getTxFastaSeq 

  Author: 	Lucie Hutchins, Scientific Software Engineer
  Department:	Department of Research, Bioinformatics
  Institution:	The Jackson Laboratory
  Date:		March , 2013

Usage: ./getCoordinates -v organismVersion -t featureType [-T itemTerm][-o resultFilename] [-a annotationSource] 

Note: 
1. If the option -T is not specified, the program will generate all the feature coordinates
   for the specified organism version
2. if the -a option is not specified, features from all the annotation sources are generated
   the features are displayed by annotation source
Supported result formats:
a. -f genePrediction -- gene prediction format (transcript base format) --> default
b. -f gtf  -- gtf format (feature base format)
***************************************************************************/
#include "ggenome.h"
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/param.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

const int MBP=1000000;
int main(int argc, char** argv)
{
  if(argc <=1){
     printf("run ./getCoordinates --help for help\n");
     printf("run ./getCoordinates --orgVersions to list all the organism versions currently stored in our database\n");
     printf("run ./getCoordinates -v organismVersion --annotations  to list all the annotation sources for the specified organism version\n");
     exit(1);
  }
  int c;static int version_flag, annot_flag,help_flag;
  string organismVersion="",dirName="",itemTerm="",annotation="",type="";
   while (1)
   {   static struct option long_options[] =
        {
           {"help", no_argument,       &help_flag, 1},
           {"orgVersions", no_argument,    &version_flag, 1},
           {"annotations", no_argument,    &annot_flag, 1},
           {"outfile", required_argument,0,'f'},
           {"itemType",    required_argument, 0, 't'},
           {"organismVersion", required_argument,0,'v'},
           {"annotSource", required_argument,0,'a'},  
           {0, 0, 0, 0}
         };
         int option_index = 0;  char *chr= NULL;
         c = getopt_long (argc, argv, "f:v:t:a:",long_options, &option_index);
         if(c == -1){break;}
         switch (c)
          {
             case 0: break;
             case 'v':if (optarg)organismVersion=optarg;break;
             case 'f':if (optarg)outFileName=optarg;break; 
             case 't':if(optarg)itemType=optarg;break;
             case 'a':if(optarg)annotation=optarg;break;
             case '?':break;
             default:
                printf("run ./getCoordinates --help to display usage\n");
                abort();
           }
     }
    tabularData bG;string seq;genomeSequence mysSeq(dirName);genomeFile gFile;
    if(!dirName.empty()){
        if(dirName.substr((dirName.length()-1),1)=="/")dirName =dirName.substr(0,dirName.length()-1);
    }
    if(help_flag){string desc="";
      string message="This program generates the coordinates of genomic items ";
      message+=" for the specified organism version and annotation source.Genomic items are either transcripts, or transcript features\n";
      printf("*********************************************************************************************\n");
      printf("\nNote: %s \n",message.c_str());
      printf("\nProgram usage: ./prog_name -v organismVersion -f itemType -d resultFilename [-a annotationSource] \n"); 
      printf("**********************************************************************************************\nUse:\n");
      printf("-f ouput file name, default stdout. if -t transcripts => a genePrediction format file is generated\n");
      printf("       if -t features => a gtf format file is generated \n");
      printf("-t the type of the genomic items (-t transcripts -> genePred annotations (default),-t features -> gtf annotations)\n");
      printf("--help to display this help \n --singleH to specify that each genomic region has one header line\n");
      printf("--orgVersions to list all the organism versions currently stored in our database\n");
      printf("--annotations  to list all the annotation sources for the specified organism version\n");
      printf("-v Organism version to use.Example -v mm9, -v hg19 (mm9 is the default organism version)\n");
      printf("-a annotation source.Example -a mm9_ensGene, -v hg19_ensGene (mm9_ensGene is the default annotation source)\n");
      printf("\nExamples:\n ./getCoordinates --orgVersions to list all the organism versions currently stored in our database\n");
      printf("./getCoordinates -v mm9 --annotations  to list all the annotation sources for mm9\n");
      printf("\n ./getCoordinates -v mm9 [-d dirName] -f transcripts --> generates ensembl transcript annotations for mm9\n");
      printf("\n ./getCoordinates -v mm9 [-d dirName] -f features --> generates ensembl transcript features annotations for mm9\n");
      printf("\n ./getCoordinates -v mm9 [-d dirName] -f transcripts -a mm9_mgiGene --> generates MGI annotations for mm9\n");
      printf("\n ./getCoordinates -v mm9 [-d dirName] -f features -a mm9_mgiGene --> generates MGI features annotations for mm9\n");
      printf("**********************************************************\n");
      mysSeq.listAllGenomes(); exit(0);
     }
     if(organismVersion.empty()){if(organism=="")organism="mouse";}
     else{ 
         if(organism=="")organism=mysSeq.getOrganismName(organismVersion);
         if(mysSeq.toUpperCase(organism)!=mysSeq.toUpperCase(mysSeq.getOrganismName(organismVersion))){
            cout<<"Your organism "<<organism<<" does not have version "<<organismVersion<<endl;
            cout<<"Make sure the config file organisms_w_chromosome_data.txt is accessible from the working directory.\n\
                  Run the program with the --help option for help.\n";
            return 0;
          }
     }
    if(!mysSeq.organismExists(organism)){
        cout<<"We do not have the genome data for "<<organism<<" --"<<organismVersion<<endl;
        cout<<"Make sure the config file organisms_w_chromosome_data.txt is accessible from the working directory.\n\
               Run the program with the --help option for a complete organisms list\n";
        return 0;
   }
   if(organismVersion.empty())organismVersion=mysSeq.getCurrentBuild(organism);
   ofstream output_f;string headerLine,line; bool stdoutt=true,stdinh=true;
   if(!outFileName.empty()){output_f.open(outFileName.c_str());if(output_f)stdoutt=false;}
   ifstream inf; genomeFile fileObject;
   if(!inputFile.empty()){ inf.open(inputFile.c_str());
       if(inf)getline(inf,headerLine); //get the header line and store the index of each field
       else{fprintf (stderr, ": Couldn't open file %s; %s\n", inputFile.c_str(), strerror (errno));
            exit (EXIT_FAILURE);
       }stdinh=false;
    }else{getline(cin,headerLine);}
    if(stdinh){ // reading from stdin
       while(getline(cin,line)){
             fileObject.getTabularFields(headerLine,line,delimiter,bG);
             if((bG.chrom=="-1")||(bG.chromStart<=0)||(bG.chromEnd<=0))continue;
             if(bG.organism.empty())bG.organism=organism;
             if(bG.organismVersion.empty())bG.organismVersion=organismVersion;
             if(unit){bG.chromStart=bG.chromStart*MBP;bG.chromEnd=bG.chromEnd*MBP;} 
             if(!stdoutt)mysSeq.displayTranscript(outputFormat,output_f,bG,sh_flag,line,len,seqType,codonType);
             else mysSeq.displayTranscript(outputFormat,cout,bG,sh_flag,line,len,seqType,codonType);
       }
    }
    else{
       while(getline(inf,line)){
             fileObject.getTabularFields(headerLine,line,delimiter,bG);
             if((bG.chrom=="-1")||(bG.chromStart<=0)||(bG.chromEnd<=0))continue;
             if(bG.organism.empty())bG.organism=organism;
             if(bG.organismVersion.empty())bG.organismVersion=organismVersion;
             if(unit){bG.chromStart=bG.chromStart*MBP;bG.chromEnd=bG.chromEnd*MBP;}
             if(!stdoutt)mysSeq.displayTranscript(outputFormat,output_f,bG,sh_flag,line,len,seqType,codonType);
             else mysSeq.displayTranscript(outputFormat,cout,bG,sh_flag,line,len,seqType,codonType);   
       }
       inf.close();
    }
  if(output_f)output_f.close();
  string message="Program complete\n"; 
  return 0;
 }
