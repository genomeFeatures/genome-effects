Get Transcript Sequence:

Summary: 
       Given a bedGraph file (UCSC format) containing transcript data for a given gene prediction source,
       this program extracts the fasta sequence data from the specified reference genome assembly build. 
       The default organism is mouse and the default organism version is the most current version of a given organism.
       Users can specify the type of sequence they want using the -t option:
        -t cds -> CDS sequence,-t tx -> exons+introns, -t ex -> exons only  sequence default
       If the user wants the CDS sequence, they can use the -c option to further specify the codon type
          -c dna -> DNA sequence default, -c rna -> RNA sequence, -c aa -> one letter amino acids sequence
             the -c aa option is not implemented yet

Input: the input  is a file:
   a. Must be a tab-delimitted file
   b. Must have a header line with the following  fields (the order is not important):
     1. Chromosome      
     2. ChromStart      
     3. ChromEnd  
     4. strand 
     5. exonStarts
     6. exonEnds
     7. cdsStart
     8. cdsEnd

If your file does not have the required header line,we will run the program assuming that 
      the first four fields of your input data follow this order:
      1.Chromosome,
      2.chromStart,
      3.chromEnd,
      4.strand,
      5.exonStarts,
      6.exonEnds,
      7.cdsStart
      8.cdsEnd

The output format spicification:
Users can specify the output type using the -o option (tab=1,commas=2,fasta=0 default)

Program usage:
        ./getTxFastaSeq -f trancriptFile.txt -F trancriptFile.fa -v mm9 -t cds -c rna -d genomeBaseDir

 Where:
      -F ouput file name , default stdout if called without the -F option
      -f bedGraph/tabular input filename containing multiple genomic regions query,default stdin (reading lines from a pipe)
      -o ouput format (tab=1,commas=2,fasta=0 default)
      -t Sequence type (-t cds -> CDS sequence,-t tx -> exons+introns, -t ex -> exons only sequence default)
      -c codon type (-c dna -> DNA sequence default,-c rna -> RNA sequence, -c aa -> one letter amino acids sequence)
              the -c option is used only when -t cds is selected (coding sequence)
      --help to display this help  --singleH to specify that each genomic region has one header line
      -v Organism version to use: -v mm9, -v hg19
      -d Directory where a config file with chromosomes data path info for every organism is stored.
         Default current working director. 
         on rockhopper: -d /hpcdata/shared/graber_lab/genomes,
         on Graber Lab: -d /data/seq  
      -l Sequence chunks length  -u Cordinate units (1=Mbp, 0=bp default)

Examples:  
1. ./getTxFastaSeq -f testcds.txt -F testcdsrna.fa -v mm9 -t cds -c rna -d /data/seq

  The above command uses the transcript annotations provided in testcds.txt file to
  generate the coding sequence of of each transcript using mm9. The sequence will be in RNA format (T->U)

2. head -300 transcript.txt |./getTxFastaSeq -v mm9 -d /data/seq -t ex
   Display the the DNA sequence (exons only) to standard out of the first 299 entries of transcript.txt
   The result will be in fasta format

2. head -300 transcript.txt |./getTxFastaSeq -v mm9 -d /data/seq -t ex -o 1
   Display the the DNA sequence (exons only) to standard out of the first 299 entries of transcript.txt
   The result will be in tab-delimited format

Note:
 for the -d argument, we maintain in Graber Lab a text file "organisms_w_chromosome_data.txt" that contains 
 a listing of all the organisms we have chromosome data for.This file is kept under the genome base directory 
 (/data/seq) and is updated when we download new genome data.

 The file is a commas separated file with the following fields:
 1. Organism group,
 2. organsim name,
 3. version, and
 4. path2genomedata
For example, the entries for human look like:
mammal,human,hg19,/data/seq/mammal/human/hg19/dat
mammal,human,hg18,/data/seq/mammal/human/hg18/dat

The dat/ directory stores each chromosome data into a one-liner file
These one-liner files are generated using two perl scripts formatfastaHeaders.pl and fa2dat.pl

