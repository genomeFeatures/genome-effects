/*
Note: be careful when using char* as a key of a map or pair
      you can have strange behavior when displaying data
      For example it was outputing a $ sign instead of the key
      use c++ string instead
*/
#include <getopt.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;
using namespace boost;

const int MB=1000000; double step=10.0;
typedef std::pair<int, int> intPair;
typedef std::pair<string,intPair> strInPair;
typedef std::pair<string,int> strIntMap;
typedef std::pair<strIntMap,string>charIntcharPair;
typedef std::pair<string,string> strPair;
typedef std::vector<char*> Fields; 
typedef std::vector<intPair> nodeIndex;          //stores <exonEnd,index> pair
typedef std::vector<int> Component;

typedef std::map<int,struct exon>junctionNetwork;  //structure to store junctions. Every junction is a node where
                                                       //the junction id is the junction id (or node id)
typedef std::map<int,Component> exonJunctionList;  //keeps track of junction ids associted to an exon id

typedef std::map<charIntcharPair,nodeIndex> exonIndex2Map; //stores chr,v1_end,type,v2_start,strand
typedef std::map<strInPair,int> exonIndexMap;

typedef std::map<int,int> edges;
typedef std::map<int,bool> updated;
typedef std::map<intPair,int> Intmap;
typedef std::map<strPair,Intmap> junctions_sorted ;

typedef adjacency_list <vecS, vecS, undirectedS> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertices_size_type VertexIndex;
typedef graph_traits<Graph>::vertices_size_type VertexIndex;
void displayBoost(Graph& G,ofstream& networks,junctionNetwork& junctions);
//junctionmap[std::make_pair(fields[0],fields[5])][std::pair(v1_end,v2_start)]=count;

//
// Graph node. a structure that stores node information
// junctions, exons, and reads, nodes --
/*
Steps:
A. input files:RUM_Unique,junctions_high-quality.bed,novel_inferred_internal_exons_quantifications_NAME,inferred_internal_exons.bed
     
   1. novel_inferred_internal_exons_quantifications_NAME - 1-base start coordinate standard
       This file gives the counts and normalized counts for all features (transcripts, exons, introns)
      we will be using exons count of exons found in RUM_Unique. Example file
   number of reads used for normalization: 11180902
      Type	Location           	min	max	Ucount	NUcount	Length
  exon 	1:8437963-8438126	13.6338	13.6338	25	0	164
  exon 1	1:9566026-9566304	8.3347	8.3347	26	0	279
  exon 2	1:9901630-9901893	0.6775	0.6775	2	0	264
  exon 3	1:9908842-9908874	2.7102	2.7102	1	0	33
  exon 4	1:9932080-9932148	58.3292	58.3292	45	0	69
  exon 5	1:9932367-9932475	31.1802	31.1802	38	0	109
  exon 6	1:10134473-10134559	45.2331	45.2331	44	0	87
  exon 7	1:11174069-11174220	8.8261	8.8261	15	0	152

  2. inferred_internal_exons.bed - 0-base start coordinate standard
  track	name="Inferred Internal Exons"  description="Inferred Internal Exons"   visibility=3    itemRgb="On"
1	3411731	3411931	.	0	.	3411731	3411931	16,78,139
1	4767528	4767652	.	0	.	4767528	4767652	16,78,139
1	4772570	4772736	.	0	.	4772570	4772736	16,78,139
1	4773953	4774108	.	0	.	4773953	4774108	16,78,139
1	4798456	4798488	.	0	.	4798456	4798488	16,78,139
1	4818581	4818647	.	0	.	4818581	4818647	16,78,139
1	4820265	4820313	.	0	.	4820265	4820313	16,78,139
1	4822308	4822379	.	0	.	4822308	4822379	16,78,139
1	4826998	4827072	.	0	.	4826998	4827072	16,78,139


   2. RUM_Unique - the unique mappers - gives one alignment per line -reads 
     format:
     i) the sequence number
    ii) chromosome
    iii) spans of the alignment in genome coordinates
    iv) the sequence of the alignment
     - all sequence is plus strand
     - sequence has a colon ":" where there is a junction
     - sequence has a +XXX+ if XXX is an insertion, e.g. +AG+ means AG
       inserted in the sequenced genome w.r.t. the reference
   Example: 
   seq.11560931a	1	3048080-3048163, 3048193-3048197	-	TTTTTTATATCCATTGTAATTCTTGTCAATTAACTAATTGTGTAGATCCTAATTTAGATAAGGAATCTACTGTTATGTTTTTGT:TAGAG
seq.541237a	1	3124999-3125005, 3125052-3125144	-	AGAGAGA:GAGAGAGAGAGAGATGTTTGTACAGAAATTTCTGTCTGAATAACATATACTTATTGCATTATTAAACTGGTGTCATAGGTAATAAAATAATCA
seq.4898027a	1	3190432-3190493, 3190496-3190533	+	CTGAGAACCTAGTGGTCTGTTCTCTTCTGGTTTGTCTTCATCTACCAGCTATTACTCCTCTT:GCTCAGAGATAAATTTAAAATACAAACATTGGAAAGAG
seq.5038720a	1	3190497-3190573, 3190575-3190597	+	CTCAGAGATAAATTTAAAATACAAACATTGGAAAGAGGAAGGAGAAAAGCACAAAATATTTAAGGCAAATTAAGCAA:TATGTACAGTCTCTTTGTGCCAG
seq.6522522a	1	3190588-3190641, 3190643-3190679	+	TTTGTGCCAGGAGGATGTGACTTTCTTTTTAAATCCCAAACCATGCTTTAAAAA:TTCAAATAAAGCAAGTCTGGTGTGAGTTGGGTTGCTT
 3. junctions_high-quality.bed : This file gives information on junctions
  The 'high quality' bed file has junctions that have known splice signals and are crossed by at least one uniquely mapping read with at least 8 bases on each side of the junction. The known junctions (those in your transcript database) are colored blue, the others are colored green. Junctions with non-canonical splice signals are colored slightly lighter. The score is the number of uniquely mapping reads crossing the junction with at least 8 bases on each side.
  So for the junction file keep also the score info - left exon score and right exon score for junctionStart and junctionEnd respectively
  
Relations: Provide an opinion on the strength of the evidence that a particular Relation is meaningful. 
We define the relation betwen three features: exon,junction, and reads.
The strength of a ralation is computed using  exon and junctions scores

New Steps:
1. index exons - load both exons files (novel_inferred_internal_exons_quantifications_NAME, and feature_quantifications_NAME)
   [chr][start][end].score=score,[chr][start][end].id=count
2. index reads where there is a junctions (the sequence contains ":") - load reads file RUM_Unique
   [pair(chr,strand)][pair(end,start)]=read_id
    [pair(chr,strand)][pair(start,end)]=read_id
  some reads have more than two pairs of coordinates, for we will have to loop
  through all the pairs in the read 
  read- exon => readEnd_i= exonEnd
3. index junctions - load junction file junctions_high-quality.bed
  read- junction=> j_start+50 = read_pair.endi, and (j_end-50)+1 = read_pair.startk
  only for reads that overlap +-8 bp with junction coordinates

4. index network (read-junctions-exons)

Note:
-- Features starts are in zero-base standard in:
   a. inferred_internal_exons.bed
   b. junctions_high-quality.bed

-- Features starts are in 1-base standard in:
   a. novel_inferred_internal_exons_quantifications_NAME
   b. feature_quantifications_NAME
  Some internal exons are only in the feature_quantifications_NAME and 
  some are only in novel_inferred_internal_exons_quantifications_NAME

1. generate exon ids
2. generate junction ids <id,junction>
   How: 
   A. get list of all exons whose exon_end = this junction_start
   B. get list of all exons whose exon_start= this junction_end
   Question: should all these sets share the same junction id?
  Answer: maybe not
  [
   generate as many junctions as we have AxB (cartesian product)
   set of all ordered pairs (a, b) where a ∈ A and b ∈ B.
   Note : the problem with the cross product approach is it will create more false positive. 
        So I suggest whe store v1s and v2s as lists. 
        this way we still identify correctly isolated junctions (not connected to another junction)
  ]
  foreach exon on A and B lists, 
    create a structure that stores exonid->listof junctions ids - this comes
   handy when building edges network
 

3. create edges <leftJunctionID,rightJunctionID>
   assumptions:
   a. two junctions of the an edge should be on the same chromosome and strand
   b. two junctions a and b share the same edge if and only if: a.right_exon=b.left_exon 
  Process:
  foreach junction
     A. get list of junctions associated with leftexon that are on the same strand
     B. get list of junctions associated with rightexon that are on the same strand
     create edges cross product of AxB 

read: from RUM_Unique
seq.101332a	1	3456693-3456737, 3495631-3495685	+	GGTGGAGATATGCCCGTGTGCAGTGAGATCTGTTGGATTCTAACG:GATTTCTCTGGAATATTGGCACAGGCTTGGCTTTCCCTCTACCCCAATCTCTTTT

Exon: in feature_quantification file
 exon 1	1:3456636-3456737	0.8768	0.8768	1	0	102

Junction: (search for junctions where e_end-50 = j_start, or e_start+offset=j_end
1	3456687	3495680	1	1	+	3456687	3495680	0,205,102	2	50,50	0,38943

Read:
seq.5784526a	1	9906358-9906383, 9908842-9908874, 9910165-9910205	-	TCAATCACTTGGTATTTTCCTAAGAG:AAGCAAATAATATCACTCCCCACACAGCAAAAG:TTCCTTTAGGAATTAGTGACAACTTTAGGCGTCTCCTCTCC
exon from novel
 exon 3	1:9908842-9908874	2.7102	2.7102	1	0	33

Read:
seq.5998543a	1	9901829-9901893, 9905203-9905237	+	CACAGGTCTGCTTTGTTGCTGTTAAGACTCTCTCATTAATTGGACAGTTACAGACTGAAACTCAA:ATCAATATTGTGCTGAAGTTAACACATTTACCAGC
exon from novel
 exon 2	1:9901630-9901893	0.6775	0.6775	2	0	264


*/
//
struct exon{ //exons are the Vs
   string strand,chrom;int left_exonEnd,right_exonStart; //junction id
   Component left_exon,right_exon,left_junction,right_junction; //list of potential exons and junctions linked to this junction
};


//a network is built using junctions and exons to create edges

/******** get fields into a vector of char ***************/
void getFields(Fields & lineFields,string line,const char* delimiter){
   if(!lineFields.empty())lineFields.clear();size_t lineSize=line.length(); 
   char* pch=strtok((char*)line.c_str(),delimiter);
   while(pch!=NULL){lineFields.push_back(pch);pch=strtok(NULL,delimiter);}
}
/// let try
//get list of chromosome from file -- only used when input file is very large
void getChromList(exonIndex2Map & exons,Fields &chromVector){
  if(!chromVector.empty())chromVector.clear();
  if(exons.empty())return;exonIndex2Map::iterator iter;std::map<string,bool> mapkeys;
  for(iter=exons.begin();iter!=exons.end();++iter){
     if(mapkeys.find(iter->first.first.first)!=mapkeys.end())continue;
     chromVector.push_back((char*)iter->first.first.first.c_str());mapkeys[iter->first.first.first]=true;
  }mapkeys.clear();
}
/***************************************************************************************************************
 Finds all the connected sub-graphs (i.e. "transcripts") 
 given a network of junctions-exons of a genomic region, global coordinates start and end of that region (global min and global max).
 The network is a map of junctions with associated list of children junctions

I need to re implement this using the new logic of global min and global max of a uniq genomic region
**************************************************************************************************************/
void FindGraphComponents(network regionMap,int start,int end,componentList bins) {
    if(junctions[node_id].right_junction.empty()){//components.push_back(node_id);
       string chrom,strand;int v1_end,v2_start;Component exonStarts,exonEnds;
       while(components.size()>0){
           int junctionid=components.back();components.pop_back();if(nodeList[junctionid].empty())skipList[junctionid]=true;
           if(chrom.empty())chrom=junctions[junctionid].chrom;if(strand.empty())strand=junctions[junctionid].strand;
           v1_end=junctions[junctionid].left_exonEnd;v2_start=junctions[junctionid].right_exonStart;
           exonEnds.push_back(v1_end);exonStarts.push_back(v2_start);
       }
       if((exonEnds.size()>0)&&(exonStarts.size()>0)){
         std::sort(exonEnds.begin(),exonEnds.end()); std::sort(exonStarts.begin(),exonStarts.end());
         out<<chrom<<"\t"<<strand<<"\t"<<exonEnds[0]<<",";
         if(exonStarts.size()>1){
            for(int i=0;i<exonStarts.size()-1;++i){out<<exonStarts[i]<<"-"<<exonEnds[i+1]<<",";}
         } out<<exonStarts[exonStarts.size()-1];out<<endl;
       }
    }
    else{
       components.push_back(node_id);
         while(nodeList[node_id].size()>0) {
            int child_nodeid=nodeList[node_id].back();nodeList[node_id].pop_back(); 
            // components.push_back(child_nodeid);
            if(skipList.find(child_nodeid)==skipList.end())
             FindGraphComponents(nodeList,junctions,child_nodeid,components,countP,out,skipList);
    } }
}
bool exonExists(int exonid,Component& exonList){
 bool exists=false;
 if(!exonList.empty()){ Component::iterator exon;
    for(exon=exonList.begin();exon!=exonList.end();++exon){if((*exon)==exonid) exists=true;}
  }return exists;
}
//filter junctions that have zero internal exons
// next filter  junctions that have exactely one internal exon
// next generate all the junctions associated to multiple exons then filter out junctions
// whose exons are not found in other junctions

/****************************************************************************************************
Junctions are nodes:
Generate an id for every exon and every junction- using the internal exons file
 and the junctions file
The junction ids will be used to build edge networks - an edge connect two junctions

Note: I store both v_start and v_end as keys to index an exon. 
      this speed exonID lookup when processing the exonJunction file (a junction (v1_end,v2_start))
     This uses twice as much memory as the number of exons but we save a lot in v id lookup
     An alternative would be to do a key scan to return all exons that have the same
     exon_end as a given v1_end (O(n))

Assumptions:
1. The internal exons file is a bed file where field[0]=chromosome, field[1]= chromStart, field[2]=chromEnd
2. The junction file is a bed file where field[0]=chromosome, field[1]= v1_end, field[2]=v2_start, field[5]=strand

*****************************************************************************************************/
void loadExonAndJuctionsMap(junctionNetwork& junctions,ifstream& exonf,ifstream& junctionf,const char* delimiter,int offset,Graph& G){ 
  exonJunctionList nodeMap,junctionChain;
  ofstream out("junction_exid.txt"); string line; int count=0; char* chromosome;
  exonIndex2Map exons; if(!exons.empty())exons.clear();Fields fields;junctionNetwork exonCordmap;
  //generate exons id
  while(getline(exonf,line)){ if(line.substr(0,5)=="track")continue;
        getFields(fields,line,delimiter); ++count;
        if(fields.size()>=3){++count;
           exons[std::make_pair(std::make_pair(fields[0],atoi(fields[1])),"start")].push_back(std::make_pair(atoi(fields[2]),count));
           exons[std::make_pair(std::make_pair(fields[0],atoi(fields[2])),"end")].push_back(std::make_pair(atoi(fields[1]),count));
           exonCordmap[count].chrom;exonCordmap[count].left_exonEnd=atoi(fields[1]);
           exonCordmap[count].right_exonStart=atoi(fields[2]);
  }}
  count=0;nodeIndex::iterator v1Iter,v2Iter; nodeIndex v1_vertexList,v2_vertexList;
  exonIndex2Map::iterator v1_vertexListIter,v2_vertexListIter;
  junctions_sorted junctionmap; int v1_end=0,v2_start=0,exonid=0,twos=0;
  ofstream out2("two-exon-junctions.txt");  string chrom,strand;
  //generate junction ids - Store junction sorted by chrom-strand-start-end
  while(getline(junctionf,line)){ 
      if(line.substr(0,5)=="track")continue; getFields(fields,line,delimiter);
      ++count;if(count%100000==0)cout<<"Processing line "<<count<<" :: "<<line<<endl;
      if(fields.size()>=6){chrom=fields[0]; strand=fields[5]; 
         v1_end=atoi(fields[1])+offset;v2_start=atoi(fields[2])-offset;
         if((v1_end)>0 &&(v2_start>0)){
            v1_vertexListIter=exons.find(std::make_pair(std::make_pair(chrom,v1_end),"end"));//get list of exons with same v1_end
            v2_vertexListIter=exons.find(std::make_pair(std::make_pair(chrom,v2_start),"start"));//get list of exons with same v2_start
            //this junction belongs to a two-exon transcript so is excluded from network
            if((v1_vertexListIter==exons.end())&&(v2_vertexListIter==exons.end())){
              out2<<chrom<<"\t"<<strand<<"\t"<<v1_end<<","<<v2_start<<endl;++twos;continue;}
            junctionmap[std::make_pair(chrom,strand)][std::make_pair(v1_end,v2_start)]+=1;
  }}}
}
void displayBoost(Graph& G,ofstream& networks,junctionNetwork& junctions){
 std::vector<int> components(num_vertices(G));int num = connected_components(G, &components[0]); 
 cout<<"Boost found a total of:"<<num<<" connected components\n";
  std::map<int,Component> componentsList;std::vector<int>::size_type vertex_id;
  for (vertex_id = 0; vertex_id != components.size(); ++vertex_id)
       componentsList[components[vertex_id]].push_back(vertex_id);
  std::map<int,Component>::iterator iter;string chrom,strand;
  int v1_end,v2_start;
  for(iter=componentsList.begin();iter!=componentsList.end();++iter){
     int node_id=iter->first;string chrom=junctions[node_id].chrom,strand=junctions[node_id].strand;
     Component exonStarts,exonEnds;
     for(int i=0; i<componentsList[iter->first].size();++i){
        node_id=componentsList[iter->first][i]; 
        int v1_end=junctions[node_id].left_exonEnd,v2_start=junctions[node_id].right_exonStart;
        if(v1_end>0){
        exonEnds.push_back(v1_end);exonStarts.push_back(v2_start);}
     }
     if((exonEnds.size()>0)&&(exonStarts.size()>0)){
        std::sort(exonEnds.begin(),exonEnds.end()); std::sort(exonStarts.begin(),exonStarts.end());
        networks<<chrom<<"\t"<<strand<<"\t"<<exonEnds[0]<<",";
        if(exonStarts.size()>1){
            for(int i=0;i<exonStarts.size()-1;++i){networks<<exonStarts[i]<<"-"<<exonEnds[i+1]<<",";}
        } networks<<exonStarts[exonStarts.size()-1];networks<<endl;
     }
  }
}
/*****************************************************************************************
******************************************************************************************/
int main(int argc, char** argv){
 string outfile,exonfile,junctionfile;int offset=0; int c; static int help_flag;
 Graph G;  
  while (1)
   {   static struct option long_options[] =
        { {"help", no_argument,       &help_flag, 1},{"outfile", required_argument,0,'f'},
          {"exonf",  required_argument,0, 'e'},{"junctionf",required_argument,0, 'j'},
          {"offset",required_argument,0, 'o'},{0, 0, 0, 0}
        };
        int option_index = 0;  char *chr= NULL; c = getopt_long (argc, argv, "o:f:e:j:",long_options, &option_index);
        if(c == -1){break;}
        switch (c)
          {case 0: break;case 'o':offset= atoi(optarg);break;
           case 'f':if (optarg) outfile=optarg;break;
           case 'e':if (optarg) exonfile= optarg;break;
           case 'j':if (optarg) junctionfile=optarg;break;case '?':break;
           default:
                printf("enter --help to display usage\n");abort();
           }
   }if(exonfile.empty()||junctionfile.empty()||help_flag){
     std::cout<<"Usage: program -e exon_filename -j junction_filename [-o offset] [-f output_filename]\n"; exit(1);
  }// get the size of the file using the struct stat
  struct stat sb; ifstream exonf(exonfile.c_str()),junctionf(junctionfile.c_str());ofstream networks(outfile.c_str());
  if (stat(exonfile.c_str(), &sb) == -1){perror("stat"); exit(EXIT_FAILURE);}  
  size_t bufsize=ceil((double)sb.st_size/(double)MB);int block_count=ceil((double)bufsize/step);
  printf("File size: %lld MB block count= %d \n",(long long) bufsize,block_count);
  //Set up exons network (Exons are nodes)
  time_t rawtime2; time (&rawtime2);const char* delimiter="\t";junctionNetwork junctions;
  loadExonAndJuctionsMap(junctions,exonf,junctionf,delimiter,offset,G);
  printf ("Program started at: %s", ctime (&rawtime2));
  displayBoost(G,networks,junctions);
 ////////////////////////
  
  time_t rawtime; time (&rawtime); printf ("Program ended at: %s", ctime (&rawtime)); 
 // cout<<"We have : "<<exons.size()<<" exons and :"<<exonNet.size()<<" junctions\n";
  
 return 0;

}

