Scripts under: /hpcdata/shared/graber_lab/bin/pl 
Here you will find the scripts you'll need to run to load  and integrate novel transcripts/isoforms in the database and
to generate the transcriptome gtf file which includes both novel and known transcripts.  

You run parseAndLoad_NovelTranscripts.pl after novel transcripts/isoforms have been generated, formatted, and stored in one directory.
You run getCordinates.pl after you've run parseAndLoad_NovelTranscripts.pl  
to generate the final transcriptome gtf file. This file will have both known and novel transcripts/isoforms.

**************** step 1. run parseAndLoad_NovelTranscripts.pl ***********************************
This script porcesses novel transcripts data to assign/generate rsd [rna seq discovery] numbers 
to novel transcripts and novel genes. The sript also loads/integrate novel transcripts supporting data into the database.
 
a. Script requirements: 
    i.  The novel transcripts files are in genePred format (0-base starts, 1-base ends)
    ii. Each file name have the format: organismVersion.annotationSource.sampleId.(.+).txt  example: mm9.mm9-WSBGene.WSB_11304.liver.txt
        where (.+) is the sample brief summary if exists or "unknown"
    iii. The load_annotations.pl script should exist in the same working directory as the script parseAndLoad_NovelTranscripts.pl
    iv. Make sure you have writing privileges on the working directory where these two scripts are ran from
        otherwise create symbolic links 

 b. Script usage: ./parseAndLoad_NovelTranscripts.pl -d filesDir
    where filesDir is a directory where novel transcripts files are stored.
    example:  ./parseAndLoad_NovelTranscripts.pl -d CC_to_database
So in brief, cd to /hpcdata/shared/graber_lab/bin/pl, and run the script from there OR
   create a bin directory in your working space, cd to bin then create symbolic links 
   /hpcdata/shared/graber_lab/bin/pl/parseAndLoad_NovelTranscripts.pl, and /hpcdata/shared/graber_lab/bin/pl/load_annotations.pl
   under bin, then run parseAndLoad_NovelTranscripts.pl script.

Assumptions as decided by Nazira and I:
 1. Every filename has the format: organismVersion.annotationSource.sampleId.(.+).txt  example: mm9.mm9-WSBGene.WSB_11304.liver.txt
    where (.+)? is the sample brief summary if exists or "unknown"
 2. annotationSource name should be consistent with what's in our gene_prediction table in the database 
 3. organismVersion name should be consistent with what's in our organism_version table in the database
 4. sample tissue if exists or "unknown"
 5. You have writing permissions where you're running the scripts from
 6. Both load_annotations.pl script and parseAndLoad_NovelTranscripts.pl are in the same working directory

**************** step 2. run getCordinates.pl ***********************************
This script generates the coordinates  of transcripts/features for the specified organism version and annotation source.
The features are stored into the database using the zero-based cordinates standard for feature starts
and 1-based cordinate standard for feature ends

## Features coordinates in gtf file are generated in 1-based <start> standard and 1-base <end>
## Features coordinates in genePrediction file are generated in 0-based <start> standard and 1-base <end>

Example: perl getCordinates.pl -v mm10 -a GRCm38-ensGene  -f mm10-ensemblTransctipts.gff -F 
The above will use Ensembl annotations for mm10 to generate all the transcripts
in gtf format and store the file (mm10-ensemblTransctipts.gtf). 

Example: perl getCordinates.pl -v mm10 -a GRCm38-ensGene  -f mm10-ensemblTransctipts.txt 
The above will use Ensembl annotations for mm10 to generate all the transcripts
in gene Prediction format and store the file (mm10-ensemblTransctipts.txt). 
