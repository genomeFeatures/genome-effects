Author: Lucie Hutchins
Date: July 2013
For generic genome features manipulation tasks wouldn't it be nice to have a single Bioinformatics tool acting as a multiplexer
that run the appropriate program based on the task at hand - instead of having to install different tools for every single task?
Think of the idea of having a single tool that, without too much hassel,can be used for
   - genome variation annotation and validation,
   - genome features sequence extraction (exons,cds,aa,ie,ei,tss,tts,introns,ex-utrs,...), 
   - file type conversion,
   - genome downloads, 
   - gene annotation download and feature coordinates extraction,
   - ...
A tool that updates its own data (gene annotations) without you having to do it manually!

Note: All the coordinate files generated from ucsc database follow the following standard:
1. all the start -> 0-base
2. all the end -> 1-base
SNPs vcf and others file types that only contain a single base pair coordinates are in 1-base

All gene annotations files should be in Bed/genePrediction format (http://genome.ucsc.edu/FAQ/FAQformat.html#format9)
- By default the annotator assumes SNP location coordinates in the SNPs file to follow 1-base standard
- By default the annotator assumes annotation coordinates in the genePrediction file to follow 0-base starts and 1-base ends standard
- By default if no input SNP file is provided, the annotator will be reading SNPs lines from standard input
- By default if no output file is provided, the annotator will be writing results lines to standard output
     if you want the result to be stored in an output file, add the option -F outputfilename
- By default if no organism version is provided, the annotator will use the current version of mouse (we can have user set default organism)
- If coordinate starts in the genePrediction file follow 1-base standard, add the --oneBaseStart flag to the program command line
- If SNP location coordinates in the SNPs file to follow 0-base standard, add the --snpZeroBaseStart flag to the program command line

You can run the progran as:
  i.  SNP Annotator 
  ii. SNP Validator 
  ii. getSNPFlanking sequence
  ii. get genome features Sequence
  ii. get genome features coordinates

      
-annotator talks to web-services via perl calls,I should probably change this to use wget instead for portability
 

