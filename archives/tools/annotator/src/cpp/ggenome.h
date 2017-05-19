/******************************************************************
  Purpose:	Develop C++ libraries with commonelly used functions
                for Bioinformatics and Computational Biology applications

  Author: 	Lucie Hutchins, Scientific Software Engineer
  Department:	Department of Research, Bioinformatics
  Institution:	The Jackson Laboratory
  Created:      September , 2011
  modified:     November   , 2013
  contact: lucie.hutchins@jax.org
           lucie.hutchins@yahoo.com
           lucie.hutchins@genomeeffects.org
*******************************************************************/

#ifndef _GGENOME__H
#define _GGENOME__H

#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<map>
#include<list>
#include<vector>
using namespace std;

//DEFINE conversion to bytes
const int GB=1073741824;
const int MB=1048576;
const int KB=1024; 

///define data structures
///
/*char* chrom_memblock;  //to temporary hold entire chromosome file in memory
  char* file_memblock;   //to temporary hold chunks of a very large file in memory
*/

typedef std::list<int> intList; 
typedef std::list<std::string> strList;
typedef std::vector<int> intVec; 
typedef std::vector<string> strVec;
typedef std::vector<struct tabularData> featureVec;



///handy for features indexing
typedef std::pair<string,long> strIntPair;//pair(chrom,pos)
typedef std::pair<strIntPair,string> strIntPairStr; //pair(pair(chrom,pos),"start") or pair(pair(chrom,pos),"end")
typedef std::pair<long,long> longPair;//pair(chrom,pos)
typedef std::pair<long, int> longIntPair; //can be handy for cordinates start and rank mapping

 ///handy to map location to count/frequency 
typedef std::map<int, int> intIntMap;
typedef std::map<long, long> longLongMap; //can be handy for cordinates start and end mapping
typedef std::map<strIntPair,long> strIntPairs; //map<pair(strand,pos),id>
typedef std::map<long,strIntPairs> posMap; //map<pos1,<pair(strand,pos2),id>>
typedef std::map<longPair,int> longIntmap; //map<pair(start,end),id>
typedef std::map<long, int> longIntMaps; //can be handy for cordinates start and rank mapping
typedef std::map<strIntPair,longIntMaps> featureIndexMap; //[pair(chrom,start)][end]=id
///stores exon by enxonEnd with associated exon start and rank
typedef std::map<long,longIntPair>exonsMap;
typedef std::map<long,strVec>longstrVec; //stores list

typedef std::map<strIntPair,longstrVec> featureMap; //<pair(chrom,start),<end, list(hits)>
typedef std::map<int,featureVec> featureList; //stores genomic ids with transcript list
typedef std::map<string,featureVec> accessionFeatureList; //stores accession ids (gene or transcript) with transcript list
///Note: I need to make sure I'm not using char* as key on maps
typedef std::map<int,intList> intListMap;   
typedef std::map<char,std::string> charStrMap;
typedef std::map<char,double> charDouble;

typedef std::map<std::string,std::string> strMap; 
typedef std::map<int,struct chromData> chromosomeMap;
///stores exon rank with associated exon cordinates
///handy to load all exons of a transcript into a map
typedef std::map<int,struct exonCord> exonRankMap;
/// map to store field name and corresponding index
/// handy for locating the field column index in a file
typedef std::map<std::string,int> fieldIndexMap;
///maps organism to their current assembly build 
///organism name is the key
typedef std::map<std::string,struct organismData> organismMap; 
///maps organism to their assembly builds, 
///organism version is the key                                                                 
typedef std::map<std::string,struct organismBuildData> organismBMap; 

///maps organism to their assembly builds, 
///db organism versionid, organismid                                                                 
typedef std::map<std::string,struct organismInDbStructure> organismInDbMap; 
///we have cases in the database where one transcript has 
///more than one translations. This structure will come handy in those
///cases 
//typedef std::map<int,translationVec> transcriptTranslationMap;

struct chromData{
      int id;
      string chromosome;
      long  size;
  };
struct exonCord{
    long start,end;int rank,firstCoding_rank,lastCoding_rank;
};
/// this structure stores all the supporting data for a given aminoacid
struct aminoacid {
    strVec codonList;
    std::string symbol,acidity,polarity;
    char oneLetterSymbol;
    double hydropathyIndex,mass,volume;
};
struct organismData{
        std::string currentBuild;
        std::string chromosomeDataDir;
        std::string localChromosomeDataDir; //path to local genome data
        std::string organismGroup;
    };

///data structure to hold data associated with a given organism assembly version
struct organismBuildData{
        std::string organism,strain;
        std::string chromosomeDataDir;       //remote path to genome data - for download pp
        std::string localChromosomeDataDir; //path to local genome data
        std::string organismGroup;
    };
///data structure to hold organism annotation data
struct organismInDbStructure{
        std::string organism;
        int organism_id,organismVersion_id;
    };
///pslData is a data structure to hold a psl record
struct pslData{
           int matches,misMatches,repMatches,nCount,qNumberInsert,qBaseInsert;
           int tNumberInsert,tBaseInsert,blockCount,qSize,qStart,qEnd;
           double tSize,tStart,tEnd;
           string strand,qName,blockSizes,qStarts,tStarts,tName;
        };

///genome data path configuration data structure 
struct genomeConfigFile{
        /// the format of the genome file name (chr1.fa)-> prefix is chr
        std::string chromosomeFilePrefix;
        std::string chromosomeFileSuffix;
        ///base directory for genome chromosome data 
        std::string genomeBase,toolBase;     //genome base and tool base local
        std::string dataDirName;            //a base directory  for tools generated data
        std::string organismsWithGenome,strain;    //a local copy of the remote genome data config file -- used for download pp
        std::string annotationsTexFile,annotationsXmlFile;///list of organism versions that have annotations in our db
        std::string webserviceUrl,server,transcriptService,remoteGenomeConfigFile,localGenomeConfigFile,host;
        
   };

///tabularData is a structure to hold fields of tab-delimited genenome annotations
/// files. These could be transcripts with the associated exons,
struct tabularData{
       string chrom;
       double chromStart, chromEnd;
       int score;
       string name, organism, organismVersion,gene,source;
       string strand; 
       intVec exonStarts, exonEnds,cdsStarts,cdsEnds; //transcript features (exons, introns, cds,...)
       
};

///snpData is a structure to hold fields of tab-delimited or commas separated
/// file types. These are SNP files 
struct snpData{
        string submitterID;
        string chrom;
        long pos,chromStart,chromEnd,snp_localID;
        string strand;
        string referenceBase,consensusBase;
};

///struct exonData is a data structure to store exon info
struct exonData{ 
       int organism_version_id;
       int chromosome_id;       // internal id for this chromosome
       int exon_start;         // start location (coord) of this exon
       int exon_end;          // end location(coord) of this transcript
       int strand;           // genome strand of this transcript (1->plus strand; 0->minus strand)
       intVec transcriptIdList;
};
///struct translationData is a data structure to store translation info
struct translationData{ 
       int cdsStart;   //the id of the exon where the translation start
       int cdsEnd;     // the id of the exon where the translation ends
       string transcript_name;int transcript_id,gene_prediction_id;
};
/// lastModified stored the file last modification information
struct lastModified{int month,year,day;};

/// class ggenome will provide global(shared) utilities to genome data manipulation .
/// This is the base class to be used by derived classes like:
/// 1. genomeFile
/// 2. genomeSequence
/// 3. aminoAcid
/// 4. genomeFeature
///This class has all the global members that deal with organisms information and some generic tasks.

class ggenome {

      protected:

         genomeConfigFile genomeData;      // data structure that hold program data paths
         organismMap organisms;           // data structure mapping organisms to their
                                          // current assembly builds
         organismBMap organismAll;       // data structure mapping organisms to all 
                                        //  the assembly builds we have
         organismInDbMap organismInDB;  // mapping organisms to all the assembly 
                                       //  builds we have
         fieldIndexMap organismIdMap;  //maps organism-organism_id
         strMap codon2aa; //codon to one letter amino-acid mapping for DNA/RNA sequence
         strMap stopCodons,startCodon;
         /// mapping between onletter and three letter annotations
         strMap symbol2OneChar; strMap name2OneChar;
         charStrMap oneChar2Symbol,oneLetter2Name,polarity,acidity;
         charDouble hydropathyIndex,mass,volume;	
        
   
      public:
         
         ggenome();
          /// Class constructor with specified tool annotations host.
         ggenome(string host);
         /// Class constructor with specified tool base directory.
         ggenome(string toolBaseDir,string host);
         ///removes leading and trailing blanks, and cariage return character from a given string
         std::string trimStr(string str);
        
         void setaminoAcid();
        ///  toUpperCase: given a DNA sequence , toUpperCase will return the  upper case
        ///  version of the sequence without modifying the original sequence
        ///  Input : string
  
         std::string toUpperCase(string str){string tempstring(str);int i =0;
        if(!str.empty()){while (tempstring[i]){tempstring[i]=toupper(tempstring[i]);i++;}}
              return tempstring;}
      
         ///  toLowerCase: given a sequence/string , this function will return the  lower case
        ///   version of the sequence without modifying the original sequence
        ///   Input : string
         std::string toLowerCase(string str){ string tempstring(str);int i =0;
              if(!str.empty()){while (tempstring[i]){tempstring[i]=tolower(tempstring[i]);i++;}}
              return tempstring;}

         ///   getReplacement:given a sequence , getReplacement will return a new
         ///   version of the original sequence with all the ocurence of token replaced/masked
         ///   by replacement
         ///  Input : a string
         std::string getReplacement(string originalString,const char replacement,const char token){
              string tempstring(originalString);int i =0;
              if(!originalString.empty()){
                for(i=0;i<originalString.length(); ++i ) {
                   if (toupper(tempstring.at(i))==toupper(token))tempstring[i]=replacement;}
              }return tempstring;}
        
         /// returns a list of positions relative to the length of sequence where token is found
         /// The positions returned follow zero-base standard [0,1,2,3,...]       
          void getPositions(string sequence,const char token,intVec& posList){int i =0;
              if(!sequence.empty()){
                 for(i=0;i<sequence.length(); ++i ){if (toupper(sequence.at(i))==toupper(token))posList.push_back(i);}
              }}
        
         /// returns a list of positions relative to the length of sequence where a digit is found
         ///The positions returned follow zero-base standard [0,1,2,3,...]       
         void getDigitPositions(string sequence,intVec& posList){;int i =0;
              if(!sequence.empty()){for(i=0;i<sequence.length(); ++i ){if(isdigit(sequence.at(i)))posList.push_back(i);}
              }}
       
         /// returns true if a given sequence is all digit (does not contain nondigit character)   
         bool isAllDigit(string sequence){bool allDigit=true;int i =0;
              if(!sequence.empty()){for(i=0;i<sequence.length(); ++i ){if(!isdigit(sequence.at(i)))allDigit=false;}
              }return allDigit;}
       
        ///  returns the numeric version of a given organism version
        ///  example: mm9 -> 9 ; mm10 -> 10. This function is handy for computing the latest version
        ///  of an organism
         int getNumericVersion(string organismVersion){int thisVersion=0;
            intVec posList;getDigitPositions(organismVersion,posList);
            if(!posList.empty()){
               size_t digitStart= posList.front();
               size_t digitEnd = (posList.back()-digitStart)+1;
               thisVersion=atoi(organismVersion.substr(digitStart,digitEnd).c_str());
            }
            return thisVersion; 
        }
        ///returns a map of all the exons of this transcript ranked
        /// this may be a function of genomeFeature I'm not sure yet
        void getExonRanksMap(struct tabularData & transcript,exonRankMap& exonMap);
        ///returns a map of all the exons of this transcript. exons cordinates are keys
        void getExonMap(struct tabularData & transcript,exonsMap& exonMap);
        ///returns a list of exons for the specified feature type: cds,3'utr,5'utr
        void getFeatureCoordinates(struct tabularData & transcript,exonRankMap& exonMap,int type);
        int getPosInCds(struct tabularData & transcript,long bpPos);
        void getOverlapExon(struct tabularData & transcript,long bpPos,struct exonCord& exon);
        void getNearestExon(struct tabularData & transcript,long bpPos,struct exonCord& exon);
        std::string filterChrom(string chr);
        ///specifies the mutation type(transition or transversion)
        ///Note: Purine: A,G ; Pyrimidine: C, T
        int getMutationType(string ref_allele,string other_allele);
        ///convert the CDS sequence (DNA/RNA) to an aminoacid sequence
        std::string cds2aaSeq(string seq);
        ///loads the configuration file into a structure -- all these functions are good candidates for virtual function
        /// because classes gsequence,genomefeatures,...could have different version implementations based on the derived class
        ///virtual void loadConfigData(void){}
        ///returns the config file structure --argument is a struct gConfigFile & config
        ///virtual void getConfigData(void);
         ///displays the config file content -- argument is a struct gConfigFile & config
         virtual void displayConfigData(void){}
       
        /********** members related to organism information manipulation starts here *****/
        void loadConfigData(string gseqConfigFilename);
           ///returns the config file structure --argument is a struct gConfigFile & config
           //void getConfigData(gseqConfigFile & config);
           ///displays the config file content -- argument is a struct gConfigFile & config
        ///loads the list of organisms with genome data into memory(all versions)
        ///both the remote list and the local list
         void loadOrganismAllMap();  

         ///displays the list of organisms with genome data with associated build
         void listCurrentGenomes(); void listAllGenomes(); 

         ///Returns the current genome db build for a given organism (mm9 - mouse)
         std::string getCurrentBuild(string organism);

         ///Returns the organism name associated with the specified build (mm9 - mouse)
        std::string getOrganismName(string organismVersion); 

        ///Returns the group of a given organism (mammal , ...)
        std::string getOrganismGroup(string organism); 
       
        ///returns true if we have the spacified organism in our database
        bool organismExists(string organism); 

         /// Returns the remote genome data path of a given organism version
         std::string getRemoteGenomeDir(string organismVersion);

         /// Returns the base genome data path of a given organism version(AJ ,mm9)
         std::string getLocalGenomeDir(string build);
         ///loads organisms with annotations into data dictionaries
         void loadAnnotatedOrganisms();
         void setAnnotatedOrg();
         ///return the genome file name of a given chromosome name, organism common name, and version
         std::string getChromosomeFileName(string chromosome,string organism,string organismBuild);
        /// reruns the number of days since filename was last modified
         int  lastModifiedDayCount(string filename);
        /// returns a structure containing the month, year, and day filename was modified
        void lastModifiedSince(lastModified &fileLMD, string filename);
        ///Returns the number of lines a given file contains.
        int getFileLineCount(string filename);
        ///Null destructor;
	~ggenome(){ organismAll.clear();organisms.clear();}
};         

/// class gsequence will provide utilities to sequence data manipulation .
/// members of this class will allow you to perform the following tasks on a given genomic region:
/// 1. Extract the corresponding fasta sequence from the chromosome data file using the sequence start and end coordinates
/// 2. complement a given sequence
/// 3. reverse a given sequence
/// 4. find CpG sites
/// 5. specify mutation type (transition, transversion)
/// 6. and more
///  Note: The chromosome data file is expected to be a one-line file where all the fasta sequences of a given chromosome
///   are concatenated into one (see data under /data/seq/organism_group_name/organismName/version/dat/)

class genomeSequence: public ggenome {
    public:
             genomeSequence(); //will use base class null constructor
            /// Class constructor with specified tool annotations host.
             genomeSequence(string host);
            /// Class constructor with specified genome directory.
            genomeSequence(string genomeBaseDir,string host);
            ///return a subsequence between loc1 and loc2
            ///The function returns an empty string if wrong input
            std::string getGenomeSeq(string ChromosomeSeqFileName,long loc1, long loc2,int start_base);
            ///reverseSeq: given a DNA sequence , reverseSeq will return the reversed
            ///version of the original sequence without modifying the original sequence
            std::string reverseSeq(string sequence);

            /// ComplementSeq: given a DNA sequence , ComplementSeq will return the  complemented
            /// version of the original sequence without modifying the original sequence
            std::string complementSeq(string sequence);

           /********** members related to SNPs  annotation starts here ************************/
           ///returns true if found 
           bool isCpGSite(string ChromosomeSeqFileName,long bpPos, const char* strand,int start_base);
           void displayConfigData();
           ///snpLocInCds stores the position of the SNP relative to coding sequence
           void getCodingSeqOld(string &sequence,exonRankMap& exonMap,long cds_start,long cds_end,
                            bool reverse_flag, string chrom_file ,long bpPos,long &snpLocInCds,long &cds_len,int start_base);
           void getCodingSeq(string &seq,struct tabularData & transcript,string chrom_file,int start_base);
           ///returns the 5' and the 3' utr sequence separated with a commas in a string
           void getUTRSeq(string &seq,struct tabularData & transcript,string chrom_file,int start_base);
           ///display the fasta sequence associated with a given genomic region(s)
          // void displayTranscript(int outputFormat,ostream& outf,struct tabularData& bG,
          //                     int sh_flag,char* line,int len,string seqType,string codonType,int start_base);
           void displayTranscript(ostream& outf,struct tabularData& bG,int len,string seqType,int featuresStart_base,string strain,int offset);
           ///Null destructor;
           ~genomeSequence(){} 
   private:
     genomeConfigFile gseqConfigFile;
};

/// class gfile will provide utilities to files manipulation .
/// members of this class will allow the following tasks on a given file:
/// 1. To parse different file formats and read fields into easy retrieval structures
/// 2. Generate a bin file given and input file and the bin step size
/// 3. and more

class genomeFile: public ggenome {
    public:
             
                /// Class constructor.
                genomeFile();
                genomeFile(string genomeBaseDir,string host);
                ///returns the number of fields given the file line and a delimiter
                int getFieldsCount(char* line, const char* delimiter);

                ///returns the mapping between the column headers and their index
                void getFieldsIndex(char* headerLine, const char* delimiter,fieldIndexMap& fieldsMap);
                
                ///parses fields from a given line given a delimiter into a vector of fields
                void getFieldsData(char* line, const char* delimiter,strVec &fields);
                ///parses fields from a given line given a delimiter into a vector of fields
                void getIntFieldsData(char* line, const char* delimiter,intVec &fields);
                ///returns true if a given SNP file does not have the appropriate format
                bool isBadSNPHeader(char* line,const char* delimiter);
                ///returns true if a given TEXT file does not have the appropriate format
                bool isBadTabularHeader(char* line,const char* delimiter);
                /// returns a snpData line with fields of interest into a structure 
                void getSNPFields(char* headerLine,char* line,const char* delimiter,struct snpData& snpFields);
                /// returns a tabularData line with fields of interest into a structure 
                void getTabularFields(char* headerLine,char* line,const char* delimiter,
                                      struct tabularData& tabularFields);
                /// returns a pslData line into a structure 
                void getPslFields(char* headerLine,char* line,const char* delimiter,struct pslData& psl);               
                ///Notes: binFile takes a tab-delimitted (fileType=1) or 
                /// a commas separated file(fileType=2) input file with two
                ///columns (x_axis=length/position and  y_axis=count/frequency) then bin
                ///data according to the user specified step. Data is sorted by x_axis column  
                void binFile(istream& inf,int xIndex,int yIndex,int stepSize,int fileType,ostream& out);
                ///convert a gtf file (feature specific) into a genePred format file
                void gtf2genePred(istream& gtff,ostream& genPredf);
                ///convert a genePred format file into a gtf file (feature specific) 
                void genePred2gtf(istream& gtff,ostream& genPredf);
                ///convert a genome fasta format file into a data file (index file) 
                void fasta2dat(istream& gtff,ostream& genPredf);
                ///reads gtf file in block by transcript
                bool getNextBlock(istream& gtff,const char* delimiter,struct tabularData& transcript);
                ///displays a transcript in gtf format or genePrediction format
                void displayTranscript(struct tabularData& transcript,ostream& genPredf,bool gtf); 
                string getTabularHeaderLine(const char* delimiter,istream& inf);
                ///Null destructor;
		~genomeFile(){}
     private:
};

///Class aminoacid
// members of this class allow you to perform the following tasks :
/// 1. Display the aminoacid main table :
///    a.  Name , b. threeLetterSymbol, c.oneLetterSymbol,
///    d. Acidity, e.Polarity , f.hydroPathyIndex , g.Mass (g/mole) , h.Volume (ang^3), i.Codon(s)
/// 2. get a list of all the codons for a given aa
/// 3. map amino acid to codon and vise versa:  codon->name->threeLetterSymbol->oneLetterSymbol --
///4. map amino acid to its characteristics (Acidity, Polarity , hydroPathyIndex , Mass ,Volume )
/// .......
///

class aminoAcid: public ggenome {
       public:
              ///Null constructor;
              aminoAcid(){ setaminoAcid();}
            
              void displayTable();
              ///Display the aminoacid main table
              bool isValidCodon(string codon);
              /// getCodon returns a list of three-character codon for the specified aa three-letter symbol
              /// or aminoacid comon name. if isDNA I true, return a DNA sequence otherwise returns an RNA sequence
              strVec getCodonList(string aa, bool isDNA);
              void getAminoacidData(string codon,struct aminoacid& aa);
              ///this function returns a one character aminoacid code for the specified codon/commonName 
              string getOneLetterSymbol(string codon);
              ///this function returns a three character aminoacid code for the specified codon/commonName
              string getThreeLetterSymbol(string codon){
                  string token= getOneLetterSymbol(codon);
                  if(!token.empty()){
                     if(oneChar2Symbol.find(token[0])!=oneChar2Symbol.end())token= oneChar2Symbol[token[0]];
                     else token="";
                  }
                  return token;
               }
               ///this function returns the comon name of the  aminoacid  for the specified codon/symbol
               string getAaName(string codon){
                  string token= getOneLetterSymbol(codon);
                  if(!token.empty()){
                      if(oneLetter2Name.find(token[0])!=oneLetter2Name.end())token= oneLetter2Name[token[0]];
                      else token="";
                   }return token;
               }
               ///this function returns the polirity of the  aminoacid  for the specified codon/symbol/commonName
              string getAaPolarity(string codon){
                  string token= getOneLetterSymbol(codon);
                  if(!token.empty()){
                     if(polarity.find(token[0])!=polarity.end())token= polarity[token[0]];
                     else token="--";
                  }
                  return token;
               }
             ///this function returns the acidity of the  aminoacid  for the specified codon/symbol/commonName
             string getAaAcidity(string codon){string token= getOneLetterSymbol(codon);
                  if(!token.empty()){
                     if(acidity.find(token[0])!=acidity.end())token= acidity[token[0]];
                     else token="--";
                  }
                  return token;
               }
            ///this function returns the hydropath Index of the  aminoacid  for the specified codon/symbol/commonName
            double getAaHydropathyIndex(string codon){
                  string token= getOneLetterSymbol(codon);double index=0;
                  if(!token.empty()){
                      if(hydropathyIndex.find(token[0])!=hydropathyIndex.end())index= hydropathyIndex[token[0]];
                   }
                  return index;
               } 
            ///this function returns the mass in g/mole of the  aminoacid  for the specified codon/symbol/commonName
            double getAaMass(string codon){
                  string token= getOneLetterSymbol(codon);double index=0;
                  if(!token.empty()){
                     if(mass.find(token[0])!=mass.end())index= mass[token[0]];
                   }
                  return index;
               }
            ///this function returns the volume in ang^3 of the  aminoacid  for the specified codon/symbol/commonName
            double getAaVolume(string codon){
                  string token= getOneLetterSymbol(codon);double index=0;
                  if(!token.empty()){
                      if(volume.find(token[0])!=volume.end())index= volume[token[0]];
                  }return index;
               } 
           // void setaminoAcid();
              ~aminoAcid(){}
};

/// Class genomeFeature provides utilities to genome annotations  manipulation .
/// members of this class allow you to perform the following tasks :
/// 1. generate genomic cordinates of different genomic features (transcripts, exons, introns, codonS, feature junctions,miRNA,...)
/// 2. extract features either from the database (very large chunks of data) or from the webservice (small chunks)
/// 3. .....

class genomeFeature: public ggenome{
       public:
        ///Null constructor. 
        genomeFeature();
         /// Class constructor with specified tool annotations host.
        genomeFeature(string host);
        genomeFeature(string toolBase,string host);
        ///lists annotated organisms and versions
        void displayOrganismInDB();
        void displayOrgAnnotList(string organismVersion,const char* delimiter,genomeFile &fileObject); 
        void generateFeatureList(string organismVersion,const char* delimiter,genomeFile &fileObject,strMap& featMap);
        ///loads organisms with their local id into a dictionary
        int getOrganismId(string organism);
        ///loads organisms versions with their local id into a dictionary
        int getOrganismVersionId(string organismVersion);
        /// returns the chromosome id given a chromosome name and organism version
        int getChromosomeId(string chromosome,string organismVersion);
        void getChromosomeData(string chromosome,string organismVersion, struct chromData& chrom);
        
        ///Index gene annotations for easy lookup
        void IndexAnnotationsByRegion(featureIndexMap& transcriptIds,featureList& transcripts,
           ifstream& inf,const char* delimiter,genomeFile &fileObject);
        void IndexAnnotationsByAcc(accessionFeatureList& genes,accessionFeatureList& transcripts,
           ifstream& inf,const char* delimiter,genomeFile &fileObject);
        void loadAnnotationTranscripts(string organismVersion,string annotationSource,string& annotations_text_file);
        //void getNearestTranscript(string chrom, long pos,featureIndexMap& transcriptIds,featureList& transcripts,tabularData& tx);
        int getNearestTranscript(string chrom, long pos,featureIndexMap& transcriptIds);
        void getOverlapTranscripts(string chrom, long pos,featureIndexMap& transcriptIds,intVec& txList);
      
       private:
        
         fieldIndexMap chromosomeIdMap;  //maps chromosome-chromosome_id
         chromosomeMap chrMap;           // maps chr_id- chr data
         ///loads chromosomes with their local id into a dictionary
         void loadChromosome(string organismVersion);
         void setOrgChromosomes(string organismVersion);
         /// generates the current list of all the gene predictions
         /// for a given organism version 
         void setOrgPredictions(string organismVersion);
};

#endif

