Genome SNPs Validator:

Summary: 
       This program will validate all single position SNPs (novel and existing) from an input stream
       against the specified  reference genome assembly build. The default organism is mouse and the default
       organism version is the most current version of a given organism.

       The program will also check for reference allele  and snp allele consistency and sets the 
       flags accordingly. Also the program tags CpG Sites, and Mutation types (transition or transversion).

The input format spicification, must have a header line :
     1. snpid [optional]
     2. Chromosome [required]  
     3. Position  [required]  
     4. ref_allele  [required]
     5. other_allele [optional] 
     6. strand [optional]
 
The output format spicification:
     1. oginal fields list
       ....
     2. QA_tag
        0-> normal, no action taken; 
        2-> ref_allele call matches reverse strand
        3-> Bad ref_allele call, does not match; 
        4-> ref_allele matches (N)
        5-> ref_allele matches the other_allele (allele swap)
     3. CpG_flag
     4. mutation type (present only if snp_allele was provided)
     5. Fields count missmatch between header and data

Program usage on the command line:
         ./genomeSNPValidator -f <snp_file_name> [-d genome_dir][-O organism][-v organismVersion] [-o output_filename] 
    Or : ./genomeSNPValidator -f <snp_file_name> [-d genome_dir] [-t filetype][-O organism][-v organismVersion] [-o output_filename]
    Or : ./genomeSNPValidator -c <chrom> -p <position> -a <ref_allele[/other_allele]> [-s strand] [-d genome_dir]

  Where:
     -h  displays this help message
     -f  input filename (for bash queries
     -c  chromosome -p bp_position -s strand -a ref_allele/other_allele  -- for single position query
     -O  Organism common name to use (example: -O mouse or -O human, ..) 
     -v  organism version (mm9,mm8,hg19,..) default mm9 if called without the -v option
     -t  input and output format (commas=2,tab=1 default)
     -o  ouput file name , default stdout if called without the -o option
     -d  Directory where a config file with chromosomes data path info for every organism is stored. Default current working directory
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

**********************************************************
Some examples:
1. Example send result to a file: ./genomeSNPValidator -f inputfilename -o outputfileName -t 1 -d /data/seq
2. Example send result to stdout : ./genomeSNPValidator -f inputfilename -d /data/seq 
3. Example input is stdin: head -300 inputfilename |./genomeSNPValidator -v mm9 -o myoutput.txt -d /data/seq
4. Example single genomic region: ./genomeSNPValidator -c chromosome -p bp_position -a A/C -d /data/seq

