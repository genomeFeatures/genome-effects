###################################################################################################
## generateGtf.pl
#  This script converts a gene prediction file to gtf format
#    
#

q{
<feature>
The following feature types are required: "CDS", "start_codon", "stop_codon". The features "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are optional. All other features will be ignored. The types must have the correct capitalization shown here.

CDS represents the coding sequence starting with the first translated codon and proceeding to the last translated codon. Unlike Genbank annotation, the stop codon is not included in the CDS for the terminal exon. The optional feature "5UTR" represents regions from the transcription start site or beginning of the known 5' UTR to the base before the start codon of the transcript. If this region is interrupted by introns then each exon or partial exon is annotated as a separate 5UTR feature. Similarly, "3UTR" represents regions after the stop codon and before the polyadenylation site or end of the known 3' untranslated region. Note that the UTR features can only be used to annotate portions of mRNA genes, not non-coding RNA genes.

The feature "exon" more generically describes any transcribed exon. Therefore, exon boundaries will be the transcription start site, splice donor, splice acceptor and poly-adenylation site. The start or stop codon will not necessarily lie on an exon boundary.

The "start_codon" feature is up to 3bp long in total and is included in the coordinates for the "CDS" features. The "stop_codon" feature similarly is up to 3bp long and is excluded from the coordinates for the "3UTR" features, if used.

The "start_codon" and "stop_codon" features are not required to be atomic; they may be interrupted by valid splice sites. A split start or stop codon appears as two distinct features. All "start_codon" and "stop_codon" features must have a 0,1,2 in the <frame> field indicating which part of the codon is represented by this feature. Contiguous start and stop codons will always have frame 0.

The "inter" feature describes an intergenic region, one which is by almost all accounts not transcribed. The "inter_CNS" feature describes an intergenic conserved noncoding sequence region. All of these should have an empty transcript_id attribute, since they are not transcribed and do not belong to any transcript. The "intron_CNS" feature describes a conserved noncoding sequence region within an intron of a transcript, and should have a transcript_id associated with it.
};
# 1. a gtf file has the following fields
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
# <seqname> is the chromosome
# <source> is the gene prediction source name
# <feature> is feature type (exon, CDS)
# [attributes] a list of attributes (gene_id;transcript_id;gene_name;feature_id; feature_rank;...)
#http://mblab.wustl.edu/GTF22.html
#
# 2. a genePred file has the following fields
# name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tid\tname2\texonFrames\n
#strain	sample	name	chrom	strand	txStart	txEnd	exCount	exStarts	exEnds
use vars qw ($opt_h $opt_f);
use Getopt::Std; getopts('hf:');
if($opt_h||!$opt_f){
 print <<HELP;
  This script converts a gene prediction file to gtf format.

HELP
exit;
}
 chomp($opt_f);

open(IN,"$opt_f");
@strains=("mm9-129S1Gene","mm9-AJGene","mm9-ensGene","mm9-CASTGene",
          "mm9-NODGene","mm9-NZOGene","mm9-PWKGene","mm9-WSBGene");
if(IN){
 #strain	sample	name	chrom	strand	txStart	txEnd	exCount	exStarts	exEnds
 $header=<IN>;
 @fileContent=<IN>;
 foreach $strain(@strains){
    $outfile="$strain-novel.gtf";open(OUT,">$outfile");
    @lines=grep(/$strain/,@fileContent);
   while(@lines){ 
      $line=shift(@lines); chomp($line);
     ($strain,$sample,$name,$chrom,$strand,$txStart,$txEnd,$exCount,$exStarts,$exEnds,$gene)=split(/\t/,$line);
      # next if(!($strain=~/mm9-C57BLGene/));
      @exonStarts=split(",",$exStarts); @exonEnds=split(",",$exEnds);$i=0;
      while($i<@exonStarts){
         print OUT "$chrom\t$sample\texon\t$exonStarts[$i]\t$exonEnds[$i]\t";
         print OUT ".\t$strand\t.\tgene_id \"$gene\"; transcript_id \"$name \";\n";++$i;
         
      }
  }
  close(OUT);
 }
 print "Program complete\n";
}
exit;
