Get Fasta Seqeuence:

Summary: 
       Given the genomic cordinates,this program extracts the fasta sequence data 
       from the specified  reference genome assembly build. The default organism is mouse and the default
       organism version is the most current version of a given organism.

Input: the input  can be either a file or a single genomic region
If input is a file:
   a. Must be a tab-delimitted file
   b. Must have the following  fields:
     1. Chromosome      
     2. ChromStart      
     3. ChromEnd  
     4. strand  (optional, default is + strand)  
 
The output format spicification:
Users can specify the output type using the -o option (tab=1,commas=2,fasta=0 default)

Program usage:
        ./prog_name [-f snpfilename][-F outputFileName][-d genomeBaseDir][-O organism][-v organismVersion]
  OR :  ./prog_name chromosome -s chromStart -e chromEnd [-d genomeBaseDir][-O organism][-v organismVersion]

 Where:
      -F ouput file name , default stdout if called without the -F option
      -f bedGraph/tabular input filename containing multiple genomic regions query,
             default stdin (reading lines from a pipe)
      -o ouput format (tab=1,commas=2,fasta=0 default)
      --help to display this help \n --singleH to specify that each genomic region has one header line
      -O Organism common name to use (example: -O mouse or -O human, ..)
      -v Organism version to use: -v mm9, -v hg19
      -c chromosome -s chromStart -e chromEnd for a single genomic range query 
      -t genomic strand (only used for a single genomic region query, -c option), default is +
      -d Directory where a config file with chromosomes data path info for every organism is stored.
         Default current working director. on rockhopper: -d /hpcdata/shared/graber_lab/genomes,on Graber Lab: -d /data/seq 
      -l Sequence chunks length\n -u Cordinate units (1=Mbp, 0=bp default)
     
Note:
 for the -d argument, we maintain in Graber Lab a text file "organisms_w_chromosome_data.txt" that contains 
 a listing of all the organisms we have chromosome data for.
 The file is a commas separated file with the following fields:
 1. Organism group,
 2. organsim name,
 3. version, and
 4. path2genomedata
For example, the entries for human look like:
mammal,human,hg19,/data/seq/mammal/human/hg19/dat
mammal,human,hg18,/data/seq/mammal/human/hg18/dat
The dat/ directory stores each chromosome data into a one-liner file

If your file does not have the required header line,we will run the program assuming that:\n";
the first four fields of your input data follow this order.
1.Chromosome,2.chromStart,3.chromEnd,4.strand.

**********************************************************
Some examples:
1. Example send result to a file: ./getFastaSeq -f inputfilename -F outputfileName -o 1 -d /data/seq
2. Example send result to stdout : ./getFastaSeq -f inputfilename -d /data/seq 
3. Example input is stdin: head -300 inputfilename |./getFastaSeq -v mm9 -F myoutput.txt -d /data/seq
4. Example query single genomic region: ./getFastaSeq -c chromosome -s chromStart -e chromEnd -o 1 -d /data/seq

