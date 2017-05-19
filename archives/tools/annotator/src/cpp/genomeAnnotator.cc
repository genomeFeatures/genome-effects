/********************************************************************************
  Purpose:      statistical analysis of genomic features -
  Goal:         creating multitasked bioinformatics tools
 
  Author: 	Lucie Hutchins, Scientific Software Engineer
  Department:	Department of Research, Bioinformatics
  Institution:	The Jackson Laboratory
  Date:		September , 2011
  modified:     November , 2013
  contact: lucie.hutchins@jax.org
            lucie.hutchins@yahoo.com
            lucie.hutchins@genomeeffects.org 
 Summary: 
   This tool acts as a BioInformatic tools multiplexer. It selects and run the appropriate 
   tool for the specific task. It can be expanded by adding more options

Note: I need to have pfam to gene converter (pfam <-> Ensembl for example)
********************************************************************************/
#include "ggenome.h"
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/param.h>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
char pathname2[MAXPATHLEN]; const int MBP=1000000;

void genome_validator(genomeFile& FileObject,genomeSequence& SeqObject,
               string headerLine,istream& inf,ostream& output_f,const char* delimiter,
               string strain,string organismVersion,int snpStart_base,strMap& featMap,genomeFeature& feature);
void loadFeature(featureMap& featureM,string filename,genomeFile& FileObject,const char* delimiter,string feature);
void genome_annotate(genomeSequence& SeqObject,genomeFeature& feature,string headerLine,istream& inf,ostream& output_f,
                     featureIndexMap& transcriptIds,featureList& transcripts,bool coding_flag,
                     const char* delimiter,string organism,string organismVersion,
                     int snpStart_base,int cordStart_base,string annotationSource);

void getFlanking(genomeFile& FileObject,genomeSequence& SeqObject,string headerLine,istream& inf,
                 ostream& output_f,const char* delimiter,string strain,
                 string organismVersion,int snpStart_base,int offset);

void getRegionSequence(featureList& transcripts,ostream& output_f,string organism,string organismVersion,
                       int featuresStart_base,string seqType,int len,string strain,genomeSequence& SeqObject,int offset);

//tool usage functions - every tool has a usage function
void displayUsage(string type);void displayDefaultUsage();void displaySeqUsage();
void displayFlankingUsage();void displayValidatorUsage();void displayAnnotatorUsage();
//returns the header line of the file
void getHeaderLine(string& headerLine,istream& inf,const char* delimiter,string& preHeader,genomeFile& FileObject);
//generates the tool_config file if not exists -- I get rid of the config file
void generate_config(ostream& outf); 
int main(int argc, char** argv)
{  
 if(argc <=1){
     printf("type ./genomeAnnotator --help for help\n");
     exit(1);
  } 
  string infileName="",pos="", chromosome="",dirName="",organismVersion="",organism="";
  string ref_allele="",consensus_allele="",strand="+",output_filename,annotationFile,annotationSource,localGenomeDir; 
  int snpStart_base=1,end_base=1,cordStart_base=0;long position=0;//default feature start coordinate to 0-base gene Pred format
  static int help_flag,validate_flag,coding_flag,annot_flag,org_list_flag,genome_list_flag,annot_list_flag,oneBaseStart_flag,snpZeroBaseStart_flag;
  static int flanking_flag,sequence_flag;int offset=0;string seqType;int len; string host="";string preHeader="",strain="";
  int c;
  while (1)
   {   static struct option long_options[] =
        {  {"help", no_argument,  &help_flag, 1},
            {"orgList", no_argument,  &org_list_flag, 1},
            {"genomeList", no_argument,  &genome_list_flag, 1},
            {"annotList", no_argument,  &annot_list_flag, 1},

            {"validate",no_argument,&validate_flag,1},
            {"annotate",no_argument,&annot_flag,1},
            {"coding", no_argument,  &coding_flag, 1},
            {"flanking", no_argument,  &flanking_flag, 1},
            {"sequence", no_argument,  &sequence_flag, 1},

            {"oneBaseStart", no_argument,  &oneBaseStart_flag, 1},
            {"snpZeroBaseStart",  no_argument,  &snpZeroBaseStart_flag, 1},

            {"resultFile", required_argument,0,'F'},
            {"snpFile",  required_argument,0, 'f'},
            {"annotFile",    required_argument, 0, 'A'},

           
            {"localGenome",  required_argument,0, 'g'},
            {"strain",  required_argument,0, 'S'},
            {"dir",    required_argument, 0, 'd'},

            {"len",  required_argument,0, 'l'},
            {"seqType",required_argument,0, 'T'},
            {"offset",  required_argument,0, 'o'},
            
            {"organismVersion", required_argument,0,'v'},
            {"annotSource", required_argument,0,'a'},
            {"host", required_argument,0,'H'},
            
            {"bpPos", required_argument,0,'p'},
            {"chrom",  required_argument,0, 'c'},
            {"strand",    required_argument, 0, 't'},
            {"refAllele", required_argument,0,'r'},
            {"snpAllele",  required_argument,0, 'C'},
            {0, 0, 0, 0}
         };
         /* getopt_long stores the option index here. */
         int option_index = 0;  char *chr= NULL;c = getopt_long (argc, argv, "F:f:d:v:a:A:c:o:p:s:r:C:t:l:T:H:g:",long_options, &option_index);
         if(c == -1){break;}
         switch (c)
          {
             case 0: break;
             case 'f':infileName=optarg;break;
             case 'F':if (optarg)output_filename=optarg;break;
             case 'd':dirName= optarg;break;
             case 'v':if (optarg)organismVersion=optarg;break;
             case 'a':annotationSource=optarg;break;
             case 'A':annotationFile= optarg;break;
             case 'H':if (optarg)host=optarg;break;
             case 'g':if (optarg)localGenomeDir=optarg;break;
             case 'c':chromosome= optarg;break;
             case 'p':if (optarg)pos= optarg;break;
             case 't':if (optarg)strand=optarg;break; 
             case 'o':if (optarg)offset=atoi(optarg);break;
             case 'l':if (optarg)len=atoi(optarg);break;
             case 'T':seqType= optarg;break;
             case 'r':if (optarg)ref_allele=optarg;break;
             case 'C':if (optarg)consensus_allele=optarg;break;
             case '?':break;
             default:
                printf("enter --help to display usage\n");
                abort();
           }
   }
  string headerLine,line; char* delimiter="\t"; bool stdoutt=true,stdint=true;tabularData bG;//string seq; 
  if(!dirName.empty()){ if(dirName.substr((dirName.length()-1),1)=="/")dirName =dirName.substr(0,dirName.length()-1);}
  DIR* dp; string currentWorkingDir=""; struct stat sb; 
  if(!dirName.empty()){ currentWorkingDir=dirName;}else{getwd(pathname2);currentWorkingDir=string(pathname2);}
  if(currentWorkingDir.substr((currentWorkingDir.length()-1),1)=="/")currentWorkingDir =currentWorkingDir.substr(0,currentWorkingDir.length()-1);
  string confi_file=currentWorkingDir+"/tool_config.txt";
  string seqtype="";
  if(validate_flag)seqtype="validate";else if(flanking_flag)seqtype="flanking";else if(sequence_flag)seqtype="sequence";
  else if(annot_flag)seqtype="annotate";
  if (stat(confi_file.c_str(), &sb) == -1){ofstream outf(confi_file.c_str());generate_config(outf);}
  if(help_flag){displayUsage(seqtype);exit(1);}
  genomeSequence mysSeq(dirName,host);if(genome_list_flag){ mysSeq.listAllGenomes();exit(1);}
  organism=mysSeq.getOrganismName(organismVersion);genomeFile gFile(dirName,host); ifstream inf; 
  genomeFeature features(host); if(org_list_flag){features.displayOrganismInDB();exit(1);}
 if(organism.empty()){ //this is where we will download the genome if needed
   cout<<"You need to index the genome for  "<<organism<<" --"<<organismVersion<<endl;
   cout<<"run the program with as ./genomeAnnotator -v "<<organismVersion<<" --host ucsc --download\n ";
   exit(1);
 }if(annot_list_flag){features.displayOrgAnnotList(organismVersion,delimiter,gFile);exit(1);}
 time_t rawtime2; time (&rawtime2);printf ("Program started at: %s", ctime (&rawtime2));
  // set both SNPs and gene annotations cordinates standards
 if(snpZeroBaseStart_flag)snpStart_base=0;if(oneBaseStart_flag)cordStart_base=1;ofstream output_f;bool stdout=true;
 if(!output_filename.empty()){output_f.open(output_filename.c_str());if(output_f)stdout=false;}
 if(chromosome.empty()){ bool stdin=true; 
     featureIndexMap  transcriptIds;featureList transcripts;//genomeFeature features;
     if(!validate_flag &&!flanking_flag){ //only index annotations if running the annotator or the getSequence programs
        if(!annotationFile.empty()){ inf.open(annotationFile.c_str());
            features.IndexAnnotationsByRegion(transcriptIds,transcripts,inf,delimiter,gFile);
            inf.close();
         }
        else{ //if the annotations file for this annotation source does not exist,
              // or the local copy is more than two weeks old, then download from server
            features.loadAnnotationTranscripts(organismVersion,annotationSource,annotationFile);
            inf.open(annotationFile.c_str());
            features.IndexAnnotationsByRegion(transcriptIds,transcripts,inf,delimiter,gFile);inf.close();
         }
     }
     //now download and index feature coordinates [miRNA,snp,repeat,...]  ---------- I'm here  
     //I'll move this under snp validator
    strMap featMap;features.generateFeatureList(organismVersion,delimiter,gFile,featMap);
    if(sequence_flag){ //adjust the header line - either the SNP file header or the annotations file header
      if(!annotationFile.empty()){inf.open(annotationFile.c_str()); 
          if(inf){getHeaderLine(headerLine,inf,delimiter,preHeader,gFile); stdin=false;}
      }else{getHeaderLine(headerLine,cin,delimiter,preHeader,gFile);}
    }else{
      if(!infileName.empty()){inf.open(infileName.c_str()); getHeaderLine(headerLine,inf,delimiter,preHeader,gFile); stdin=false;
      }else{getHeaderLine(headerLine,cin,delimiter,preHeader,gFile);}
    }
    if(validate_flag){ //run the validator
      if(!stdout)output_f<<preHeader;else cout<<preHeader;
      if(stdin){
          if(!stdout)genome_validator(gFile,mysSeq,headerLine,cin,output_f,delimiter,strain,organismVersion,snpStart_base,featMap,features);
          else genome_validator(gFile,mysSeq,headerLine,cin,cout,delimiter,strain,organismVersion,snpStart_base,featMap,features);
       }else{ 
          if(!stdout)genome_validator(gFile,mysSeq,headerLine,inf,output_f,delimiter,strain,organismVersion,snpStart_base,featMap,features);
          else genome_validator(gFile,mysSeq,headerLine,inf,cout,delimiter,strain,organismVersion,snpStart_base,featMap,features);
        }
    }else if(sequence_flag){ //run getSequence 
          if(!stdout)getRegionSequence(transcripts,output_f,organism,organismVersion,cordStart_base,seqType,len,strain,mysSeq,offset);
         else  getRegionSequence(transcripts,cout,organism,organismVersion,cordStart_base,seqType,len,strain,mysSeq,offset);
    }else if(flanking_flag){  //run getFlanking
         if(stdin){
          if(!stdout)getFlanking(gFile,mysSeq,headerLine,cin,output_f,delimiter,strain,organismVersion,snpStart_base,offset);
          else getFlanking(gFile,mysSeq,headerLine,cin,cout,delimiter,strain,organismVersion,snpStart_base,offset);
       }else{ 
          if(!stdout)getFlanking(gFile,mysSeq,headerLine,inf,output_f,delimiter,strain,organismVersion,snpStart_base,offset);
          else getFlanking(gFile,mysSeq,headerLine,inf,cout,delimiter,strain,organismVersion,snpStart_base,offset);
        }
    }else{ //run the annotator
         if(!stdout)output_f<<preHeader;else cout<<preHeader;
         if(stdin){
             if(!stdout)genome_annotate(mysSeq,features,headerLine,cin,output_f,transcriptIds,transcripts,coding_flag, delimiter, strain,
                     organismVersion,snpStart_base,cordStart_base,annotationSource);
             else genome_annotate(mysSeq,features,headerLine,cin,cout,transcriptIds,transcripts,coding_flag, delimiter, strain,
                     organismVersion,snpStart_base,cordStart_base,annotationSource);
         }else{
             if(!stdout)genome_annotate(mysSeq,features,headerLine,inf,output_f,transcriptIds,transcripts,coding_flag, delimiter, strain,
                     organismVersion,snpStart_base,cordStart_base,annotationSource);
             else genome_annotate(mysSeq,features,headerLine,inf,cout,transcriptIds,transcripts,coding_flag, delimiter, strain,
                     organismVersion,snpStart_base,cordStart_base,annotationSource);
         }
     }
  time_t rawtime; time (&rawtime); printf ("Program ended at: %s", ctime (&rawtime));
 }
 else{ // single genomic position
    if(!mysSeq.isAllDigit(pos)){
       printf("The SNP base pair position should be a numeric value\nrun ./genomeSNPValidator --help for help \n");return 0;
     } 
    position =atol(pos.c_str());chromosome=mysSeq.filterChrom(chromosome);
    string genomeFilename=mysSeq.getChromosomeFileName(chromosome,strain,organismVersion);
    string varch,mutation="N/A";int is_conflict=0;
    string local_ref_allele=mysSeq.toUpperCase(mysSeq.getGenomeSeq(genomeFilename,position,position,snpStart_base));
    if(!ref_allele.empty()){
     if(local_ref_allele!=ref_allele){ //reference alleles mismatch
        if(local_ref_allele==mysSeq.complementSeq(ref_allele))is_conflict=2;//possible reverse strand
        else if(local_ref_allele==consensus_allele)is_conflict=5; //allele swap
        else if(local_ref_allele==mysSeq.complementSeq(consensus_allele)){
                 varch=ref_allele;ref_allele=consensus_allele;
                 consensus_allele=varch; varch=""; is_conflict=2; //possible reverse strand
        }else is_conflict=3;  //our reference allele does not match either the complement of the specified
                             // source reference allele or the other allele
     }if(is_conflict==5){   //allele swap
        varch=ref_allele; ref_allele=consensus_allele;consensus_allele=varch;
     }else if(is_conflict==2){//reverse strand
          ref_allele=mysSeq.complementSeq(ref_allele);consensus_allele=mysSeq.complementSeq(consensus_allele);
     }}string iscpg="No";
    cout<<"organism (organismVersion)\tchromosome:position\tlocal_ref_allele/consensus_allele\tstrand\tis_conflict\tmutationType\tiscpgSite"<<endl;
    if(mysSeq.isCpGSite(genomeFilename,position,strand.c_str(),snpStart_base))iscpg="Yes";
     if(!consensus_allele.empty()){
        int mutation_type=mysSeq.getMutationType(local_ref_allele,consensus_allele);
         if(mutation_type==2)mutation="Transversion";else mutation="Transition";
     }
    cout<<organism<<"("<<organismVersion<<")\t"<<chromosome<<":"<<position<<
         "\t"<<local_ref_allele<<"/"<<consensus_allele<<"\t"<<strand<<"\t"<<is_conflict<<"\t"<<mutation<<"\t"<<iscpg<<endl;
   }
  
 return 0;
}
void getHeaderLine(string& headerLine,istream& inf,const char* delimiter,string& preHeader,genomeFile& FileObject){ 
 preHeader="";size_t pos; string temp="",prefix="";
 while(getline(inf,headerLine)){if(headerLine.empty())continue;
    pos=headerLine.find("#");temp=headerLine;
    if(pos!=std::string::npos){temp=headerLine.substr(pos+1);prefix=headerLine.substr(0,1);}
    if((!FileObject.isBadSNPHeader((char*)temp.c_str(),delimiter))||
       (!FileObject.isBadTabularHeader((char*)temp.c_str(),delimiter))){headerLine=temp;break;}//
    else  preHeader+=headerLine+"\n";
 }preHeader+=prefix;
}
//consider cases the genome data for this organism does not exists - no chromosme data but we do have annotations
//Note: some SNPs lines have the allele fields (Ref, Alt)-- A,      G,T --  as a commas separated list
// In that case, get thecartesian product of RefxAlt and process each heterogenous combination
/****
b.) The result of the SNP coding version will have the following fields appended to it:
           Exon,ExonID(s),Intron,Intergenic,codingClass,TranscriptID(s),GeneSymbol,geneStrand,ensemblDbSchema
      Where:
         Exon: Yes/No -- Yes if SNP is locate on exon
         ExonID(s): a commas separated list of all the exons of this gene containing the SNP
                where each exon has additional information including: snpFrame,CDS_len,snpPosInCDS,codonPosInProtein    
         Intron: Yes/No -- Yes if SNP is locate on intron
         Intergenic: Yes/No -- Yes if SNP is locate on intergenic region
         Function_class: a commas separated list of all the functional implications of this SNP on this gene
         TranscriptID(s): a commas separated list of all the transcripts of this gene containing the SNP
         GeneSymbol: gene symbol 
So pointer and reference occupies same amount of memory - as they are both addresses. As a general rule,
    Use references in function parameters and return types to define attractive interfaces.
    Use pointers to implement algorithms and data structures.
******/
void genome_annotate(genomeSequence& SeqObject,genomeFeature& feature,string headerLine,istream& inf,ostream& output_f,
                     featureIndexMap& transcriptIds,featureList& transcripts,bool coding_flag,
                    const char* delimiter, string strain,string organismVersion,
                    int snpStart_base,int cordStart_base,string annotationSource){
  fieldIndexMap fieldsMap; strVec fieldsData,headerFields;snpData snp;string line;genomeFile FileObject; //genomeFeature feature;
  aminoAcid codonData;FileObject.getFieldsIndex((char*)headerLine.c_str(),delimiter,fieldsMap);int miscount=0,is_conflict=0;
  intVec txIdsList;exonsMap exonMap;string ref_aa,consensus_aa;
  string mutation="Transition",genomeFilename,local_ref_allele="",iscpg="";
  string gene="---";int transcrID=0;string codon=".",c_codon=".",r_aa=".",c_aa=".",snpClass=".",exon_ID=".",tx_id=".",gene_strand=".";
  string is_exon="No",is_intron="No",is_utr="No",is_intergenic="No",is_cds="No";
  FileObject.getFieldsData((char*)headerLine.c_str(),delimiter, headerFields);
  output_f<<headerLine<<delimiter<<"is_Exon"<<delimiter<<"is_Intron"<<delimiter<<"is_UTR"<<delimiter<<"is_Coding"<<delimiter<<"is_Intergenic"<<
         delimiter<<"ExonID(s)"<<delimiter<<"TranscriptID(s)"<<delimiter<<"Gene_Symbol"<< delimiter<<"gene_strand"<<
         delimiter<<"functionClass"<<delimiter<<"aa_characteristics"<<
         delimiter<<"Is_CpG"<<delimiter<<"MutationType"<<delimiter<<"QA_flag"<<
         delimiter<<"annotationSource"<<delimiter<<"OrganismAssemblyVersion"<<endl;struct stat sb;
  while(getline(inf,line)){ if(line.empty())continue; 
       FileObject.getSNPFields((char*)headerLine.c_str(),(char*)line.c_str(),delimiter,snp);is_conflict=0;miscount=0;
       FileObject.getFieldsData((char*)line.c_str(),delimiter,fieldsData); string local_ref_allele;
       if(headerFields.size()!=fieldsData.size()){is_conflict=17;}
       if(miscount>0){ int count=fieldsData.size();while(count<headerFields.size()){line+="\t.";++count;}}
       if(snp.pos>0){
          snp.chrom=SeqObject.filterChrom(snp.chrom); iscpg="NA";
          string genomeFilename=SeqObject.getChromosomeFileName(snp.chrom,strain,organismVersion);
          gene=".";transcrID=0;codon=".";c_codon=".";r_aa=".";c_aa=".";snpClass=".";exon_ID=".";tx_id=".";gene_strand=".";
          int exonRank=0,exonStart=0,exonEnd=0,intronRank=0;is_exon="No";is_intron="No";is_utr="No";is_intergenic="No";is_cds="No";
          string annotation="";long cds_start,cds_end;string varch;is_conflict=0; string aa_char=".";
          /////////////// do some QA on this SNP ///////////////////////////////////////////
          if(stat(genomeFilename.c_str(), &sb) != -1){ 
             local_ref_allele=SeqObject.toUpperCase(SeqObject.getGenomeSeq(genomeFilename,snp.pos,snp.pos,snpStart_base));
             if(local_ref_allele!=snp.referenceBase){   //reference alleles mismatch
                if(local_ref_allele==SeqObject.complementSeq(snp.referenceBase))is_conflict=2;//possible reverse strand
                else if(local_ref_allele==snp.consensusBase)is_conflict=5; //allele swap
                else if(local_ref_allele==SeqObject.complementSeq(snp.consensusBase)){
                 varch=snp.referenceBase;snp.referenceBase=snp.consensusBase;
                 snp.consensusBase=varch; varch=""; is_conflict=2; //possible reverse strand
                }else is_conflict=3;  //our reference allele does not match either the complement of the specified
                                   // source reference allele or the other allele
             }if(is_conflict==5){//allele swap
               varch=snp.referenceBase; snp.referenceBase=snp.consensusBase;snp.consensusBase=varch;
             }else if(is_conflict==2){//reverse strand
              snp.referenceBase=SeqObject.complementSeq(snp.referenceBase);
              snp.consensusBase=SeqObject.complementSeq(snp.consensusBase);
            }iscpg="No";
            if(SeqObject.isCpGSite(genomeFilename,snp.pos,snp.strand.c_str(),snpStart_base))iscpg="Yes";   
          } 
         mutation="Transition";if(SeqObject.getMutationType(local_ref_allele,snp.consensusBase)==2)mutation="Transversion";
         feature.getOverlapTranscripts(snp.chrom,snp.pos,transcriptIds,txIdsList);
         //////////////////////////////////////////////////////////
          if(txIdsList.empty()){ //intergenic SNP  //get neareast exon from this location
               snpClass="Intergenic";is_intergenic="Yes"; featureVec::iterator txIter;
               int txid=feature.getNearestTranscript(snp.chrom, snp.pos,transcriptIds);
               if(transcripts.find(txid)!=transcripts.end()){ //for each transcript that contain this SNP annotate
                      for(txIter=transcripts[txid].begin();txIter!=transcripts[txid].end();++txIter){
                        tx_id+= (*txIter).name;
                }}
              
               string annotation=line+delimiter+is_exon+delimiter+is_intron+delimiter+is_utr+delimiter+is_cds+delimiter+is_intergenic;
                      annotation+=delimiter+exon_ID+delimiter+tx_id+delimiter+gene+delimiter+gene_strand+delimiter+snpClass;
                      annotation+=delimiter+aa_char+delimiter+iscpg+delimiter+mutation;
              output_f<<annotation<<delimiter<<is_conflict<<delimiter<<annotationSource<<delimiter<<organismVersion<<endl;
          }else{ intVec::iterator vecIter; featureVec::iterator txIter;
            //foreach region id that contains this SNP, get list of ASSOCIATED transcripts and process
            //I need to store this info by gene, then display it later
               for(vecIter=txIdsList.begin();vecIter!=txIdsList.end();++vecIter){
                   if(transcripts.find(*vecIter)!=transcripts.end()){ //for each transcript that contain this SNP annotate
                      for(txIter=transcripts[*vecIter].begin();txIter!=transcripts[*vecIter].end();++txIter){
                         if(line.find((*txIter).name)==std::string::npos)continue;//I will remove this later
                          is_exon="No";is_intron="No";is_utr="No";is_intergenic="No";is_cds="No";
                          gene=".";transcrID=0;snpClass=".";exon_ID="."; tx_id=".";gene_strand=".";
                          exonRank=0;exonStart=0;exonEnd=0;intronRank=0;
                          cds_start=((*txIter).cdsStarts[0]>0)?(*txIter).cdsStarts[0]:0;
                          cds_end=((*txIter).cdsEnds[0]>0)?(*txIter).cdsEnds[0]:0; 
                          cds_start=(cds_end==cds_start)?0:cds_start;
                          cds_end=(cds_end==cds_start)?0:cds_end;
                          exonCord target_exon;feature.getOverlapExon(*txIter,snp.pos,target_exon);
                          if(target_exon.rank==0){snpClass="Intron";is_intron="Yes";
                            feature.getNearestExon(*txIter,snp.pos,target_exon);
                          }else{ is_exon="Yes";
                             if(cds_start==0)snpClass="Exon - non coding transcript";
                             else{
                                 if(snp.pos<cds_start){is_utr="Yes";
                                    if((*txIter).strand=="+"){snpClass="5'UTR - Exon";}else{snpClass="3'UTR - Exon";}
                                 }else if(snp.pos>cds_end){is_utr="Yes";
                                    if((*txIter).strand=="+"){snpClass="3'UTR - Exon";}else{snpClass="5'UTR - Exon";}
                                 }else{snpClass="Exon-Coding ";is_cds="Yes";}
                         }}
                        long snpLocInCds=0,cds_len=0;int postRelProt=0;string seq=""; 
                        string r_codon=".",r_aa=".",s_codon=".",s_aa="."; int posInProtein=0,frame=-1;
                        bool reverse_flag=false;if((*txIter).strand!="+")reverse_flag=true;
                       if((coding_flag)){ //check if we need to include the codon info - only if chromosome file exists
                          if((is_cds=="Yes")){snpLocInCds=SeqObject.getPosInCds(*txIter,snp.pos);
                             if(snpLocInCds>0){
                               SeqObject.getCodingSeq(seq,*txIter,genomeFilename,cordStart_base);int offset=0;
                               if((seq.length()%3)!=0){
                                  if((seq.length()%3)==1){offset=2;seq="NN"+seq;} //start codon has one leading base pair nna-aaa-aaa-aaa
                                  else{offset=1; seq="N"+seq;}//start codon has two leading base pair naa-aaa-aaa-aaa
                               }posInProtein=(snpLocInCds+offset)/3 ;frame=(snpLocInCds+offset)%3;
                               if((snpLocInCds+offset)%3==0){frame=3;}else ++posInProtein; 
                               try{ r_codon=seq.substr((snpLocInCds+offset)-frame,3);}catch(exception& e){ r_codon=".";} 
                               ref_aa=codonData.getThreeLetterSymbol(r_codon);
                               if(r_codon.length()==3){ s_codon=r_codon;
                                   if(snp.strand!=(*txIter).strand)s_codon.replace(frame-1,1,SeqObject.complementSeq(snp.consensusBase));
                                   else s_codon.replace(frame-1,1,snp.consensusBase);
                               }s_aa=codonData.getThreeLetterSymbol(s_codon);
                             snpClass=(s_aa==ref_aa)?"Exon-Coding Synonymous":"Exon-Coding NonSynonymous";
                       }}}
                      annotation=line+delimiter+is_exon+delimiter+exon_ID+delimiter+is_intron+delimiter+is_intergenic;
                      annotation+=delimiter+is_utr+delimiter+is_cds+delimiter+snpClass;
                      annotation+=delimiter+(*txIter).name+delimiter+(*txIter).gene;
                      annotation+=delimiter+(*txIter).strand+delimiter+annotationSource;
                      annotation+=delimiter+organismVersion+delimiter+mutation+delimiter+iscpg;
                      output_f<<annotation<<delimiter<<"rank ";
                      if(is_cds=="No")output_f<<target_exon.rank<<endl;
                      else{output_f<<target_exon.rank<<delimiter<<frame<<","<<
                              snpLocInCds<<","<<posInProtein<<" Codon:"<<ref_aa<<","<<r_codon<<","<<s_aa<<","<<s_codon<<endl;
                      }
                     
           }}}}
         
          
}}}
/*************************************************************************************************************
 SNPsValidator: runs QA on user's input SNPs file using the specified reference genome and sets flags accordingly.
 If UCSC has annotations for some genomic features (dbSNP SNPs, simple Repeat,miRNA,and Microsatellite )
 for the specified organism assembly version, the validator will also check
 whether or not your SNPs call overlap these annotations and sets annotation flags accordingly.
If the local index of reference genome is not found, SNPs calls will be checked agains 
 
class= enum('unknown', 'single', 'in-del', 'het', 'microsatellite', 'named', 'mnp', 'insertion', 'deletion')
***************************************************************************************************************/
void genome_validator(genomeFile& FileObject,genomeSequence& SeqObject,string headerLine,istream& inf,ostream& output_f,
                      const char* delimiter,string strain,string organismVersion,int start_base,strMap& featMap,genomeFeature& features){
  fieldIndexMap fieldsMap; strVec fieldsData,headerFields;snpData snp;string line;
  FileObject.getFieldsIndex((char*)headerLine.c_str(),delimiter,fieldsMap);int miscount=0,is_conflict=0;
  strMap featureFiles;
  string mutation="Transition",genomeFilename,local_ref_allele="",iscpg="";
  FileObject.getFieldsData((char*)headerLine.c_str(),delimiter, headerFields);struct stat sb;
  output_f<<headerLine<<delimiter<<delimiter<<"OrganismAssemblyVersion"<<delimiter<<"MutationType"<<
          delimiter<<"QA_flag"<<delimiter<<"prog_refAllele/prog_concensusAllele"<<delimiter<<"Is_CpG";
  if(!featMap.empty()){ strMap::iterator iter; string feature_file="";//load dbSNP,simpleRepeat,miRNA if any
       for(iter=featMap.begin();iter!=featMap.end();++iter){
           features.loadAnnotationTranscripts(organismVersion,iter->first,feature_file);
           if(stat(feature_file.c_str(),&sb)!=-1){ //index 
              featureFiles[iter->first]=feature_file;
              if(iter->first.find("snp")!=std::string::npos)
                 output_f<<delimiter<<"db"<<iter->first<<"-rs#:class:strand";
              else continue; // output_f<<delimiter<<iter->first;
  }}}output_f<<endl; 
  //load features
 featureMap snps,microsat,repeat,mirna;//features are mapped by genomic coordinates <pair(pair(chrom,start),strand),<end, list(hits)>
 strMap::iterator iter;featureMap::iterator fiter;
 if(!featureFiles.empty()){
    for(iter=featureFiles.begin();iter!=featureFiles.end();++iter){
         if(iter->first.find("snp")!=std::string::npos)loadFeature(snps,iter->second,FileObject,delimiter,iter->first);
         else continue;
        // else if(iter->first.find("microsat")!=std::string::npos)loadFeature(microsat,iter->second,FileObject,delimiter,iter->first);
        // else if(iter->first.find("simpleRepeat")!=std::string::npos)loadFeature(repeat,iter->second,FileObject,delimiter,iter->first);
        // else if(iter->first.find("miRNA")!=std::string::npos)loadFeature(mirna,iter->second,FileObject,delimiter,iter->first);
         
 }}
 while(getline(inf,line)){ if(line.empty())continue;
       FileObject.getSNPFields((char*)headerLine.c_str(),(char*)line.c_str(),delimiter,snp);is_conflict=0;miscount=0;
       FileObject.getFieldsData((char*)line.c_str(),delimiter,fieldsData); string local_ref_allele;
       if(headerFields.size()!=fieldsData.size()){is_conflict=17;}
       if(snp.pos>0){
          snp.chrom=SeqObject.filterChrom(snp.chrom);string varch;is_conflict=0;
          string genomeFilename=SeqObject.getChromosomeFileName(snp.chrom,strain,organismVersion);
          if(stat(genomeFilename.c_str(),&sb)!=-1){
             local_ref_allele=SeqObject.toUpperCase(SeqObject.getGenomeSeq(genomeFilename,snp.pos,snp.pos,start_base));
             if(local_ref_allele!=snp.referenceBase){   //reference alleles mismatch
               if(local_ref_allele==SeqObject.complementSeq(snp.referenceBase))is_conflict=2;//possible reverse strand
               else if(local_ref_allele==snp.consensusBase)is_conflict=5; //allele swap
               else if(local_ref_allele==SeqObject.complementSeq(snp.consensusBase)){
                    varch=snp.referenceBase;snp.referenceBase=snp.consensusBase;
                    snp.consensusBase=varch; varch=""; is_conflict=2; //possible reverse strand
                }else is_conflict=3;  //our reference allele does not match either the complement of the specified
                                      // source reference allele or the other allele
            }//allele swap
            if(is_conflict==5){varch=snp.referenceBase; snp.referenceBase=snp.consensusBase;snp.consensusBase=varch;}
            else if(is_conflict==2){//reverse strand
                 snp.referenceBase=SeqObject.complementSeq(snp.referenceBase);
                 snp.consensusBase=SeqObject.complementSeq(snp.consensusBase);
            }iscpg="No";if(SeqObject.isCpGSite(genomeFilename,snp.pos,snp.strand.c_str(),start_base))iscpg="Yes";
          }else is_conflict=-1;
       } mutation="Transition";longstrVec::iterator liter;
       if(SeqObject.getMutationType(local_ref_allele,snp.consensusBase)==2)mutation="Transversion";
       output_f<<line<<delimiter<<organismVersion<<delimiter<<mutation<<delimiter<<is_conflict<<
               delimiter<<local_ref_allele<<"/"<<snp.consensusBase<<delimiter<<iscpg;
       string name=""; bool more=true; 
       for(iter=featMap.begin();iter!=featMap.end();++iter){name="";
           if(iter->first.find("snp")!=std::string::npos){
              fiter=snps.upper_bound(std::make_pair(snp.chrom,snp.pos));if(fiter!=snps.end())fiter++;
            }else continue;
            /*else if(iter->first.find("miRNA")!=std::string::npos){
              fiter=mirna.upper_bound(std::make_pair(snp.chrom,snp.pos));if(fiter!=mirna.end())fiter++;
            }else if(iter->first.find("microsat")!=std::string::npos){
              fiter=microsat.upper_bound(std::make_pair(snp.chrom,snp.pos));if(fiter!=microsat.end())fiter++;
            }else if(iter->first.find("simpleRepeat")!=std::string::npos){
              fiter=repeat.upper_bound(std::make_pair(snp.chrom,snp.pos));if(fiter!=repeat.end())fiter++;
            }*/more=true; 
            while(more){ 
                  for(liter=fiter->second.begin();liter!=fiter->second.end();++liter){
                     if((fiter->first.second<=snp.pos)&&(liter->first>=snp.pos)){
                         for (unsigned i=0; i<liter->second.size(); ++i)name+=liter->second[i];
                  }}--fiter;
                  try{ if(fiter->first.first==snp.chrom)more=true;else{more=false;}}
                  catch(exception& e){more=false;} 
            }name=(name=="")?".":name;output_f<<delimiter<<name;
       }output_f<<endl;
  }
}
void loadFeature(featureMap& featureM,string filename,genomeFile& FileObject,const char* delimiter,string feature){
 string line;ifstream feat(filename.c_str()); strVec fieldsData;getline(feat,line);//remove header lines
 fieldIndexMap fieldsMap;int chromindex=-1,strandindex=-1,startindex=-1,endindex=-1,nameindex=-1,classindex=-1;
 FileObject.getFieldsIndex((char*)line.c_str(),delimiter,fieldsMap);
 if(fieldsMap.find("chrom") != fieldsMap.end())chromindex=fieldsMap["chrom"];
 if(fieldsMap.find("chromstart") != fieldsMap.end())startindex=fieldsMap["chromstart"];
 if(fieldsMap.find("chromend") != fieldsMap.end())endindex=fieldsMap["chromend"];
 if(fieldsMap.find("name") != fieldsMap.end())nameindex=fieldsMap["name"];
 if(fieldsMap.find("strand") != fieldsMap.end())strandindex=fieldsMap["strand"];
 if(fieldsMap.find("class") != fieldsMap.end())classindex=fieldsMap["class"];
 while(getline(feat,line)){
     FileObject.getFieldsData((char*)line.c_str(),delimiter,fieldsData);
     if(fieldsData.size()==fieldsMap.size()){
          try{
          string strand=(strandindex>=0)? fieldsData[strandindex] :"all";string name=fieldsData[nameindex];
          long chromStart=(startindex>=0)?atol((char*)fieldsData[startindex].c_str()):0;
          long chomEnd=(endindex>=0)?atol((char*)fieldsData[endindex].c_str()):0;
          string chrom=(chromindex>=0)?FileObject.filterChrom(fieldsData[chromindex]):"";
          if(feature.find("snp")!=std::string::npos)name+=":"+fieldsData[classindex];
          name+=":"+strand;
          if(chrom!="")featureM[std::make_pair(chrom,chromStart)][chomEnd].push_back(name);
         }catch(exception& e){}
 }}}
/**** *********************************************************************************************
Gets the  upstream and downstream flanking sequences from the specified position
and returns a sequence s such that s= upstream_seq_offsetLen -pos-downstream_seq_offsetLlen
The program appends the following three fields to the result:
a. Flanking Coordinates -> chromStart-ChromEnd
b. Flanking Sequence  -> fasta sequence of the flanking
c. OrganismAssemblyVersion-Strain -> the reference genome and strain used to extract the sequence
Note: The result is in the same format as the input SNPs file
***************************************************************************************************/
void getFlanking(genomeFile& FileObject,genomeSequence& SeqObject,string headerLine,istream& inf,
                 ostream& output_f,const char* delimiter,string strain,string organismVersion,int snpStart_base,int offset){
  fieldIndexMap fieldsMap; strVec fieldsData,headerFields;snpData snp;string line,flanking="";
  FileObject.getFieldsIndex((char*)headerLine.c_str(),delimiter,fieldsMap);
  string genomeFilename;FileObject.getFieldsData((char*)headerLine.c_str(),delimiter, headerFields);long flanking_start=0,flanking_end=0;
  output_f<<headerLine<<delimiter<<"Flanking Coordinates"<<delimiter<<"Flanking Sequence"<<delimiter<<"OrganismAssemblyVersion-Strain"<<endl;
  offset=(offset<=0)?25:offset;
  while(getline(inf,line)){ if(line.empty())continue;
       FileObject.getSNPFields((char*)headerLine.c_str(),(char*)line.c_str(),delimiter,snp);
       FileObject.getFieldsData((char*)line.c_str(),delimiter,fieldsData); 
       if(snp.pos>0){
          snp.chrom=SeqObject.filterChrom(snp.chrom);
          genomeFilename=SeqObject.getChromosomeFileName(snp.chrom,strain,organismVersion);
          flanking_start=snp.pos-offset;flanking_end=snp.pos+offset;
          flanking=SeqObject.getGenomeSeq(genomeFilename,flanking_start,flanking_end,snpStart_base);
          if(snp.strand!="+") flanking=SeqObject.reverseSeq(SeqObject.complementSeq(flanking));
          output_f<<line<<delimiter<<flanking_start<<"-"<<flanking_end<<delimiter<<flanking<<delimiter<<organismVersion<<
                  "-"<<strain<<endl;
  }}
}
/**** *********************************************************************************************
generates and display the formatted genomic sequence for the specified sequence type (-T seqType ).
By default -if sequence type is not specified - returns a genomic sequence delimitted by chromStart and chromEnd.
seqType can be one of the following:
 ex     - a dna sequence of a concatenation of all the exons of that transcript. exon1exon2.....,exonk
          where k is the last exon of the transcript
 cds    - a dna sequence of all the coding exons of a transcript-> CDS sequence 
 aa     - an aminoacid sequence of all the coding exons of a transcript - aa sequence  
 exon   - a dna sequence of every exon of that transcript.- one exon per line 
 intron - a dna sequence of every intron of that transcript.- one intron per line 
 tex    - a dna sequence of a concatenation of all the exons and introns of that transcript.Introns and exons are separated by a commas 
            format: exon1,intron1,exon2,intron2,.....,intronk-1,exonk where k is the last exon of the transcript 
 tss    - a genomic sequence delimitted by txStart-offset and txStart+offset -> transcript start site
 tts    - a genomic sequence delimitted by txEnd-offset and txEnd+offset -> Terminal transcript site
 ei     - an exon-intron junction - exonEnd-offset and exonEnd+offset
 ie     - an  intron-exon junction - intronEnd-offset and intronEnd+offset
 utr    - a dna sequence of all exons in the 5'UTR of a coding transcript (if any)
            and
          a dna sequence of all exons in the 3'UTR of a coding transcript (if any)
Input:
a. transcripts -> annotations index 
b. featuresStart_base -> specifies whether feature start coordinates are 0-base or 1-base 
c. --strain -> the reference strain for the specified organimsVersion - needed for local genome indexing - default none
d. -T -> the type of the sequence to extract - see above
e. --len -> result sequence is displayed in chuncks size - default extracted sequence length
f. -v -> organismVersion ,organism assembly version to use for sequence extraction
g. --offset -> offset size to use for sequence types (tss,tts,ei,and ie)

Output: The program returns a fasta file and every sequence in the file contains a header
The header line is a tab-delimited line with two fields: the transcriptID field and the Attributes field
The Attributes field - Every sequence type has the following core attributes separated with a ":" character.
a.chrom; b. chromStart - chromEnd; c.strand; other attributes (key-value pairs) include: featurerank,featureCount,featureLength,and feature 

Promoters are regions of DNA that promote transcription and, in eukaryotes, are found at -30, -75, and -90 base pairs upstream from the
*******************************************************************************************************/
void getRegionSequence(featureList& transcripts,ostream& output_f,string organism,string organismVersion,
                       int featuresStart_base,string seqType,int len,string strain,genomeSequence& SeqObject,int offset){
  int count=0;featureVec::iterator txIter;featureList::iterator idIter; 
 for(idIter=transcripts.begin();idIter!=transcripts.end();++idIter){
      for(txIter=idIter->second.begin();txIter!=idIter->second.end();++txIter){
          if((*txIter).organism.empty())(*txIter).organism=organism; 
          if((*txIter).organismVersion.empty())(*txIter).organismVersion=organismVersion;
          if((seqType=="cds")||(seqType=="aa")){
             if((*txIter).cdsStarts[0]!=(*txIter).cdsEnds[0]){ ++count;
                 SeqObject.displayTranscript(output_f,*txIter,len,seqType,featuresStart_base,strain,offset);
          }}
          else SeqObject.displayTranscript(output_f,*txIter,len,seqType,featuresStart_base,strain,offset);++count;
          if(count==20)break;
      }
    if(count==20)break;
  }
}
void displayUsage(string type){
  if(type=="sequence")displaySeqUsage();
  else if(type=="flanking")displayFlankingUsage();
  else if(type=="validate")displayValidatorUsage();
  else if(type=="annotate")displayAnnotatorUsage();
  else displayDefaultUsage();
}
void displayDefaultUsage(){
 string note="Genome Annotator v1.0";
 string contact="lucie.hutchins@jax.org OR lucie.hutchins@genomeeffects.org"; 
 string message="Genome Annotator is a BioInformatic tools multiplexer that selects and run the appropriate tool for the specified task at hand.";
 message+="\nThe tool can be used as: \
           \nsnpAnnotator         -> annotates your SNPs file using the specified genes and genome feature annotations and the reference genome \
           \nsnpValidator         -> validates your SNP calls against the specified reference genome, known dbSNPs, and other genome features\
          \nsnpFlanking extractor -> extracts the associated SNP flanking sequence for each SNP in your file using the specified reference genome\
      \ngenomeSequence extractor  -> extracts the genome sequence of the specified genome features using the specified reference genome\
      \nfeatureCoordinates extractor  -> extracts genome coordinates of the specified features using the specified annotation source and organismVersion\n";
      printf("*********************************************************************************************\n");
      printf("                            %s ",note.c_str());
      printf("\n\n*********************************************************************************************\n");
      printf("Note: %s ",message.c_str());
      printf("\n***** Program usage  *****\n***Browsing the server");
      printf("\n./genomeAnnotator --orgList --host ucsc          #displays list of organisms that have gene annotations in UCSC "); 
      printf("\n./genomeAnnotator --orgList                     #displays list of organisms that have gene annotations on our server"); 
      printf("\n./genomeAnnotator -v mm9 --annotList            #displays all gene annotations we have for mm9 in our integrated database"); 
      printf("\n./genomeAnnotator -v mm9 --annotList --host ucsc   #displays all gene annotations UCSC has for mm9 "); 
      printf("\n./genomeAnnotator --genomeList --host ucsc        #displays list of organisms that have chrom assembly data in UCSC "); 
      printf("\n./genomeAnnotator --genomeList                   #displays list of organisms that have chrom assembly data on our server\n"); 
      printf("\n***Running the Tool");
      printf("\n./genomeAnnotator --help  --annotate      #to display the snpAnnotator usage"); 
      printf("\n./genomeAnnotator --help  --validate      #to display the snpValidator usage");
      printf("\n./genomeAnnotator --help  --flanking      #to display the snpFlanking extractor usage");
      printf("\n./genomeAnnotator --help  --sequence      #to display the genomeSequence extractor usage\n"); 
      printf("\n./genomeAnnotator --help  --coordinate    #to display the featureCoordinates extractor usage\n"); 
      printf("\nFor issues and suggestions, please contact the tool developper at: %s\n",contact.c_str()); 
      printf("*********************************************************************************************\n");
}
void displayAnnotatorUsage(){
string usage="The program annotates your SNPs file using the specified genes and genome feature annotations and the reference genome.\
\nIf a local index of the reference genome is not found,the codon and amino acid information will not be included.\
\n\nThe program appends the following fields to the result:\
\na. is_Exon           # Yes/No -- Yes if SNP is locate on exon\
\nb. is_Intron         # Yes/No -- Yes if SNP is locate on intron\
\nc. is_UTR            # Yes/No -- Yes if SNP is locate on 5'UTR or 3'UTR\
\ne. is_Coding         # Yes/No -- Yes if SNP is locate on coding region of the exon\
\nf. is_Intergenic     # Yes/No -- Yes if SNP is locate on intergenic region\
\ng. ExonID(s)         #a commas separated list of all the exons of this gene containing the SNP\
where each exon has additional information including: snpFrame,CDS_len,snpPosInCDS,codonPosInProtein\
\nh. TranscriptID(s)   #a commas separated list of all the transcripts of this gene containing the SNP\
\ni. Gene_Symbol       #gene symbol\
\ng. gene_strand       #gene strand\
\nk. functionClass     #a commas separated list of all the functional implications of this SNP on this gene\
\nl. aa_characteristics #Ref-Polarity:acidity:mass:volume-hydropathyIndex/Alt-Polarity:acidity:mass:volume-hydropathyIndex\
\nm. is_CpG                          #Yes,No, or NA if you do not have a local genome index for the specified organism version\
\nn. MutationType                    #Transversion,Transition, or NA if your input SNP line does not provide the SNP allele\
\no. QA_flag                         # 0 -> no issue found,2 -> SNP detected on - strand,5 -> ref and alt alleles were swapped,\
\np. annotationSource \
\nq. OrganismAssemblyVersion-Strain  #the reference genome and strain used to validate the ref allele call\
\nNote: The result is in the same format as the input SNPs file.\
\n*** Program Usage:\
\n./genomeAnnotator --coding -v organismAssemblyVersion [--strain genomeStrainName] -f snpsFile [--snpZeroBaseStart]\
-a annotationsSourceName [--oneBaseStart][--strain genomeStrainName][--host annotationHost][-F resultFile] \nOR\
\n./genomeAnnotator --coding -v organismAssemblyVersion [--strain genomeStrainName] -f snpsFile [--snpZeroBaseStart]\
-A annotationsFileName [--oneBaseStart][--strain genomeStrainName][--host annotationHost][-F resultFile]\
\nOR\
\n./genomeAnnotator -v organismAssemblyVersion [--strain genomeStrainName] -f snpsFile [--snpZeroBaseStart]\
-A annotationsFileName [--oneBaseStart][--strain genomeStrainName][--host annotationHost][-F resultFile]\
\n*** WHERE:\
\n-v        #organismVersion ,organism assembly version used to validate the ref allele call, example: -v mm9 \
\n-f        #SNPs input file in tab-delimited format. for example: -f mySNPsFile.vcf\
\n-F        #result file name -default standard out\
\n--coding  #specifies that the annotator should include the codon/aminoacid information - default none\
\n-A        #gene annotations input file in bed/genePred format. for example: -A myDownloadedGenePredictionFile.txt\
\n--oneBaseStart      #specifies that feature start coordinates are 1-base , the default is 0-base starts and 1-base ends\
\n--snpZeroBaseStart  #specifies that SNP coordinates are 0-base , the default is 1-base starts and 1-base ends\
\n--strain  #the reference strain for the specified organimsVersion - needed for local genome indexing - default none,example: -v mm9 --strain AJ\
\n\n";

printf("*********************************************************************************************\n");
 printf("%s ",usage.c_str());
}
void displayValidatorUsage(){
 string usage="The program runs QA on user's input SNPs file using the specified reference genome and sets flags accordingly.\
\nIf UCSC has feature annotations(dbSNP SNPs, simple Repeat,miRNA,and Microsatellite ) \
 for the specified organism assembly version,\nthe validator will also checks\
 whether or not your SNPs call overlap these annotations and sets annotation flags accordingly.\
\nIf the local index of reference genome is not found, SNPs calls will be checked against\
\nUCSC available features(dbSNP SNPs, simple Repeat,miRNA,and Microsatellite )\
\n\nThe program appends the following fields to the result:\
\na. OrganismAssemblyVersion-Strain #the reference genome and strain used to validate the ref allele call\
\nb. MutationType    #Transversion,Transition, or NA if your input SNP line does not provide the SNP allele\
\nc. QA_flag         # 0 -> no issue found,2 -> SNP detected on - strand,5 -> ref and alt alleles were swapped,\
\n                     3 -> ambigous ref and alt alleles call, -1 -> the local genome index was not found\
\nd. prog_refAllele/prog_concensusAllele   #What the ref and alt alleles look like relative to the + strand after our QA\
\ne. Is_CpG          #Yes,No, or NA if you do not have a local genome index for the specified organism version\
\nf. miRNA           #this column exists only if UCSC has miRNA data for the specified organism version.\
\ng. microsat        #this column exists only if UCSC has microsatellite data for the specified organism version.\
\nh. simpleRepeat    #this column exists only if UCSC has simpleRepeat data for the specified organism version.\
\ni. dbsnpXXX-rs#    #this column exists only if UCSC has SNPs data for the specified organism version.XXX is dbSNP version\
\nk. dbsnpXXX-class  #this column exists only if UCSC has SNPs for the specified organism version.\
\n                   class= enum('unknown', 'single', 'in-del', 'het', 'microsatellite', 'named', 'mnp', 'insertion', 'deletion')\
\nNote: The result is in the same format as the input SNPs file.\
\n*** Program Usage:\
./genomeAnnotator --validate -v organismAssemblyVersion -f snpsFile [--snpZeroBaseStart][--strain genomeStrainName][-F resultFile]\
\n*** WHERE:\
\n-v        #organismVersion ,organism assembly version used to validate the ref allele call, example: -v mm9 \
\n-f        #SNPs input file in tab-delimited format. for example: -f mySNPsFile.vcf\
\n-F        #result file name -default standard out\
\n--snpZeroBaseStart  #specifies that SNP coordinates are 0-base , the default is 1-base starts and 1-base ends\
\n--strain  #the reference strain for the specified organimsVersion - needed for local genome indexing - default none,example: -v mm9 --strain AJ\
\n\n";

printf("*********************************************************************************************\n");
 printf("%s ",usage.c_str());
}
/******* Displays get SNP flanking extractor tool usage *************************************************************/
void displayFlankingUsage(){
string usage="The program gets the  upstream and downstream flanking sequences from the specified bp position \
and returns a sequence s such that s is beteewn [pos-offset and pos+offset]. \
The program appends the following three fields to the result:\
\na. Flanking Coordinates -> chromStart-ChromEnd\
\nb. Flanking Sequence  -> fasta sequence of the flanking\
\nc. OrganismAssemblyVersion-Strain -> the reference genome and strain used to extract the sequence\
\nNote: The result is in the same format as the input SNPs file.";
 usage+="\n\n*** Program Usage:\n";
 usage+="./genomeAnnotator --flanking -v organismAssemblyVersion -f snpsFile [--offset offsetLen]";
 usage+="[--snpZeroBaseStart][--strain genomeStrainName][-F resultFile]\n";
usage+="*** WHERE:\
\n-v        #organismVersion ,organism assembly version to use for sequence extraction, example: -v mm9 \
\n-f        #SNPs input file in tab-delimited format. for example: -f mySNPsFile.vcf\
\n-F        #result file name -default standard out\
\n--snpZeroBaseStart  #specifies that SNP coordinates are 0-base , the default is 1-base starts and 1-base ends\
\n--strain  #the reference strain for the specified organimsVersion - needed for local genome indexing - default none,example: -v mm9 --strain AJ\
\n--offset  #size of the offset to use(default size: 25) - example: --offset 50\n";
usage+="\n*** Example output:\
\n4621385	10	10793539	A	C	2	0	Intronic	.	10793514-10793564\
ATAAAACATATTTCACATAAGTTAAAATTTAAGGTGTGGACATTTATGAAG	mm9-\
\n4621576	10	122874541	T	C	1	0	Intronic	.	122874516-122874566\
AACAATGCTTACCTCTGTGTTTCCTTATGTAGATCTAATAACACAAGGTTT	mm9-\n";

printf("*********************************************************************************************\n");
 printf("%s ",usage.c_str());
}
/******* Displays sequence extractor tool usage *************************************************************/
void displaySeqUsage(){
 string usage="Extracts and displays the formatted genomic sequence of the specified sequence type (-T seqType ) ";
        usage+="using the reference genome and the specified annotation coordinates file or annotations source name.";
        usage+="By default - if no seqType specified- returns genomic sequences delimitted by chromStart and chromEnd.";
  usage+="The program returns a fasta file and every sequence in the file contains a header\
\nThe header line is a tab-delimited line with two fields: the transcriptID field and the Attributes field.\
\nThe Attributes field - Every sequence type has the following core attributes separated with a ':' character:\
\na.chrom; b. chromStart - chromEnd; c.strand; other attributes (key-value pairs) include:featurerank/rank,featureCount/exonCount,featureLength,and feature.\n";
 usage+="*** Example output:\n>ENSMUST00000115538	1:4775720-4775920:reverse strand:exonCount 5:feature tss\
\nCTGGAAGCTGGAAGACATAGCCTACCGACGCGGATGCGGTCGCTCTGCAGCACGACAGGGGCGGGGCAGTACGGCCGCTGCAGCGCGACAGGGGCCGGGCGG\
AGCCGTAAAGCAACGCGGAGTCTACGCCGCTTCCTGGGCGCCCTTGAAGCACCGCGGGTCATGGCTGGCACGGCGCGCGGCTGCGGGACCAGCCTGGA\n";
 usage+=">ENSMUST00000116652	1:4483487-4486444:reverse strand:featureLength 249:feature 5'UTR\
\nCTGAAGTGCGGTTGGCCCCAACACTCCTCCCAAAGTATCTATCAAGAGAATGGTCAGCAGAAGTTAGATCTAGTGAGCAGCACCTCCAGACATCTGA\
ATTTCAGCCTTCCTATTTCCCCAAGAGGTCTTGGCGCCAGCGCCCGGCTCCAGCCAGTTTTCCCCAAGGCTAGCTTCCGATCCCTGCCTCAGGGTCGG\
GGGAAGCGGCGTGTCCCGTGGCCATAGCAGAGCTCGGGGTCGGTCTGGAGAGCC\
\n>ENSMUST00000116652	1:4481793-4481796:reverse strand:featureLength 3:feature 3'UTR\
\nCGG\n";
 usage+=">ENSMUST00000130201	1:4774086-4774286:reverse strand:rank 2:exonCount 5:feature ie\
\nCTGGTGCAATGGATAAACAGCATTAATATGGAATATGCTATGAGATGAGACAGTTTTCAGAAATGTATCTGAAGATCTTAAAACTCACTTTTT\
CTTGTAGGAAAGACGTCCAAGAGATCGGAGAAGAGGTAGGAAGTGTGGCAGAGGCCATAAAGGAGAAAGGCAGAGAGGAACCCGGCCAAGGCTGGGCTTTGAGGGAG\n";
usage+=">ENSMUST00000132625	1:4775553-4775753:reverse strand:rank 1:exonCount 2:feature ei\
\nTGGCACGGCGCGCGGCTGCGGGACCAGCCTGGACCTGCTGCGGTCCTTGCCGAGAGTGAGCCTGGCCAACCTGAAGCCCAGTCCTAACTCCAGA\
AAACGGGTAAACGCTTGGCGGCGGTGCGGGGAGTGCAGGGGCCGCGTCACTCTGGCGATGGACGCGTGCGGGCGCCTGGCTCCCTTCCGCGTGCGCCCAGTGTTAA\n";
  usage+="*** Program Usage:\n";
  usage+="./genomeAnnotator --sequence -v organismAssemblyVersion -A annotationsCoordFile [-T seqType][--len seqChunksLen][--offset offsetLen]";
  usage+="[--oneBaseStart][--strain genomeStrainName][--host annotationHost][-F resultFile]\n      OR\n";
  usage+="./genomeAnnotator --sequence -v organismAssemblyVersion -a annotationsSourceName [-T seqType]";
  usage+="[--len seqChunksLen][--offset offsetLen][--oneBaseStart][--strain genomeStrainName][--host annotationHost][-F resultFile]\n";
 usage+="*** WHERE:\
\n-v        #organismVersion ,organism assembly version to use for sequence extraction, example: -v mm9 \
\n-a        #gene annotations source name for example : -a refGene \
\n-A        #gene annotations input file in bed/genePred format. for example: -A myDownloadedGenePredictionFile.txt\
\n-T        #the seqType of the sequence to extract - see section on seqType\
\n-F        #result file name -default standard out\
\n--oneBaseStart  #specifies that feature start coordinates are 1-base , the default is 0-base starts and 1-base ends\
\n--len     #specifies the sequence chuncks size in the result - default extracted sequence length, example: --len 50 \
\n--strain  #the reference strain for the specified organimsVersion - needed for local genome indexing - default none,example: -v mm9 --strain AJ\
\n--offset  #size of the offset to use for sequence types tss,tts,ei,and ie only. default size: 100, example: --offset 50\
\n--host    #remote gene prediction server name, example: --host ucsc\n";

 usage+="\n*** seqType can be one of the following:\
 \nex     - a concatenation of a dna sequence of all the exons of a transcript. exon1exon2.....,exonk\
 \n         where k is the last exon of the transcript\
 \ncds    - a concatenation of dna sequence of all the coding exons of a transcript-> CDS sequence\
 \naa     - a concatenation of aminoacid sequence of all the coding exons of a transcript - aa sequence";
  usage+="\nexon   - a dna sequence of every exon of a transcript, one exon per line\
 \nintron - a dna sequence of every intron of that transcript.- one intron per line";
  usage+="\ntex    - a dna sequence of a concatenation of all the exons and introns of a transcript.\
 \n         Introns and exons are separated by a commas.format: exon1,intron1,exon2,intron2,.....,intronk-1,exonk \
 \n         where k is the last exon of the transcript\
 \ntss    - a genomic sequence delimitted by txStart-offset and txStart+offset -> transcript start site\
 \ntts    - a genomic sequence delimitted by txEnd-offset and txEnd+offset -> Terminal transcript site\
 \nei     - an exon-intron junction delimitted by exonEnd-offset and exonEnd+offset\
 \nie     - an intron-exon junction delimitted by intronEnd-offset and intronEnd+offset\
 \nutr    - a dna sequence of all exons in the 5'UTR of a coding transcript (if any)\
 \n         and a dna sequence of all exons in the 3'UTR of a coding transcript (if any)\n\n";
 printf("*********************************************************************************************\n");
 printf("%s ",usage.c_str());
}
//generates the tool_configuration file
void generate_config(ostream& outf){
 cout<<"generating config file-- Not I don't really need the config file- given that the annotator has a defined dir structure\n";
string config="#Remote server setting\n\
#### I should actually have this on the web server################\
##  Web server url\n\
WEBSERVER=http://demon.jax.org\n\
## transcript web service path -- relative to WEBSERVER\n\
TRANSCRIPT_SERVICE=/transcriptdb/web-services/\n\
#MGI_ANNOT_SERVER=ftp://ftp.informatics.jax.org/pub/mgigff/\n\
#MGI_ANNOT_FILE=MGI.gff3.gz\n\
*********************************\
# Local server setting\n\
## Set the absolute path to the tool root directory, where the tool will be writing and reading\n\
##   default to current working directory\n\
TOOL_BASE=\n\
##############################################################\
## Set the path to the base directory of genome indexes.\n\
## The base directory is called 'genomes' and is created under\n\
## TOOL_BASE. Genomes indexes are stored here by organism assembly version. When you download the tool,it comes with Mouse and Human genomes\
## (hg19,hg18,mm10,mm9).If you want to use another organism assembly version from UCSC, run the annotator with the --download option.\n\
## If you want to use you local genome, run the annotator first with the --indexLocalGenome path2localGenome -v mm9 --strain localgenomeName \
## for example: mm9/AJ_genome  will store indexes of AJ genome\
## The tool will create this directory if this does not exists under the tool base directory\n\
GENOME_BASE=genomes\n\
## Set the path to where to/you store a copy of genome with chromosome data config file localy -- relative to TOOL_BASE\n\
## Make sure you have write permissions on TOOL_BASE\n\
## each chromosome is a fasta file with the filename format *.fa\n\
## the chromsome file name format is specified in the tool config file\n\
##set up chromosome file name format . Example chr2.fa , the prefix will be chr, the suffix will be .fa\n \
## Example 2.fa, the prefix is empty\n\
CHROM_PREFIX=\n\
CHROM_SUFFIX=.fa";

 outf<<config<<endl;
}

