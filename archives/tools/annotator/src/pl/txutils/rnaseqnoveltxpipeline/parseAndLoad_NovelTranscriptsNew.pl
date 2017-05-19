#!/usr/bin/perl

#############################################################################
# This script porcesses novel transcripts data to
# 1. assign/generate rsd [rna seq discovery] numbers 
# 2. generate genePrediction file and load novel transcripts in db
# 3. generate rnaseq_transcript file and load transcript-strain-sample mapping in db
#
# Note: 
#  Transcript rsd_ numbers are organism version-strain specific
#  Every novel transcript is given a rsd_ accession id (rsd_[number]) 
#  where number is a positive integer value generated when a novel transcript is detected.
# The process:
#  The following steps are taken when processing novel transcripts:
#  1. Get the database current max of rsd_[number] and current max of rsd_g[number] 
#     (the number part of the accession id) into a global variable(maxNumber)
#  2. Create and load a data structure (dbrsdmap)to store database current novel transcripts by strain
#     and the associated rsd numbers --> this is a dictionary
#    loadCurrentTxRsd($organism_version_id,$gene_pred_id) -> 
#         ${$rsd_hash_ref}{"$organism_version_id-$gene_pred_id"}{"$key"}{"$exonStarts:$exonEnds"}=$tx_name;
#         where  $key="$chrom-$strand-$txstart-$txend";
#     --> this is empty if there is no novel transcripts in db
#  3. Strore sample files by strain : $strainFiles{"$strain"}{"$sample"}=$file
#  4. Foreach $strain:
#    i. Index the new rnaSeq file into data structures to facilitate detection of identical transcripts   
#    ii. foreach novel transcript (txNew) from generated data structure, do:
#       a. assign the rsd number
#          a.1 if(Exists(dbrsdmap{txNew})) -> thisAccession=dbrsdmap{txNew}
#          a.2 else{++maxNumber;  thisAccession=_rsd_[maxNumber];dbrsdmap{txNew}=thisAccession
#       b. generate organismversion-strain-sampleID-thisAccession record 
#       c. generate genePrediction record 
#    iii load generated files into database
#
## Assumptions as decided by Nazira and I:
## 1. Every filename has the format: organismVersion.annotationSource.sampleId.(.+)?.txt  
#      example: mm9.mm9-WSBGene.WSB_11304.liver.txt where (.+)? is the sample brief summary if exists or "unknown"
## 2. annotationSource name should be consistent with what's in our gene_prediction table in the database
##    it is strain based
## 3. organismVersion name should be consistent with what's in our organism_version table in the database
## 4. sample tissue if exists or "unknown"
## Example: mm10.GRCm38-ensGene68.KO10123.testes.txt
#  
##############################################################################
use strict;use warnings;
use DBI;use vars qw ($opt_h $opt_d); use Getopt::Std; use Time::localtime;use Cwd;
$dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
getopts('hd:');
if(($opt_h)|| !($opt_d)){
    print <<HELP;
    **************************************

    Script usage: ./parseAndLoad_NovelTranscripts.pl -d filesDir
           example:  ./parseAndLoad_NovelTranscripts.pl -d CC_to_database
    where filesDir is a directory where novel transcripts files are stored.
    These files are in genePred format and file names follow the specified naming convention:
    organismVersion.annotationSource.sampleName.tissue.txt -
    example: mm10.GRCm38-ensGene68.KO10123.testes.txt
 
  

HELP
exit;
}
my $cwdir = getcwd; chomp($cwdir);#Current working directory
################################################ set up db queries 
#get the current max rsd_number for novel transcripts
$query="select MAX(cast(substring(transcript_name from locate('rsd_',transcript_name)+4) as unsigned)) ";
$query.=" from transcript_by_annotation WHERE transcript_name like 'rsd_%'";
$qh_getCurrentMaxRsdNumber= $dbh->prepare($query);
#get the current max rsd_number for novel genes
$query="select MAX(cast(substring(novel_gene_name from locate('rsd_g',novel_gene_name)+5) as unsigned)) ";
$query.=" from novel_genes WHERE novel_gene_name like 'rsd_g%'";
$qh_getCurrentGeneMaxRsdNumber= $dbh->prepare($query);
$qh_getnovelTxid=$dbh->prepare("select distinct transcript_name from transcript_by_annotation where transcript_id=? and gene_prediction_id=?");
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$qh_getOrgV = $dbh->prepare("select organism_version_id from organism_version where ucsc_db=?");
$qh_getSampleid=$dbh->prepare("select sample_id from rnaseq_sample where sample_name=?");
$qh_insertSampleid=$dbh->prepare("insert into rnaseq_sample(sample_name,sample_desc) values(?,?)");
#GET NOVEL transcripts accession ids this section is used in loadCurrentTxRsd()
$query="select distinct transcript_id,exon_start,exon_end";
$query.=" from transcript_exon t, exon e where transcript_name =? ";
$query.=" and t.organism_version_id=?  and gene_prediction_id=? ";
$query.=" and t.exon_id=e.exon_id order by exon_start ";
$qh_getexons=$dbh->prepare($query);
$query="select distinct transcript_name from transcript_by_annotation
        where organism_version_id=? and transcript_name like 'rsd_%' and gene_prediction_id=? ";
$qh_getTxRsds=$dbh->prepare($query);
$query="select tx_start,tx_end,chromosome_name, strand from transcript t, chromosome c where transcript_id=? ";
$query.=" and organism_version_id=? and c.chromosome_id=t.chromosome_id ";
$qh_getTxCord=$dbh->prepare($query);
#load rnaseq_transcript table
$qh_loadRnaSeqTx=$dbh->prepare("load data local infile ? into table rnaseq_transcript ignore 1 lines");
$query="select count(*) from rnaseq_transcript where organism_version_id=? and gene_prediction_id=?";
$qh_getRowCount=$dbh->prepare($query);

##########################################################################################
# loadCurrentRsd: index current novel transcripts for the specified
# organism version and strain  into a data structure
# that facilitates transcript -rsd accession id lookup
#
sub loadCurrentTxRsd{
  my($rsd_hash_ref,$organism_version_id,$gene_pred_id)=@_;
  $qh_getTxRsds->execute($organism_version_id,$gene_pred_id);
  if($qh_getTxRsds->rows>0){
    while(($tx_name)=$qh_getTxRsds->fetchrow_array()){
       $qh_getexons->execute($tx_name,$organism_version_id,$gene_pred_id);
       $exonStarts=""; $exonEnds="";$tx_id=0;
       next if($qh_getexons->rows<=0);
       while(($transcript_id,$exon_start,$exon_end)=$qh_getexons->fetchrow_array()){
            if($exonStarts eq ""){
               $exonStarts="$exon_start"; $exonEnds="$exon_end";$tx_id=$transcript_id;}
            else{$exonStarts.=",$exon_start"; $exonEnds.=",$exon_end";}
       }
       if($tx_id>0){ $qh_getTxCord->execute($tx_id,$organism_version_id);
         if($qh_getTxCord->rows>0){
           ($txstart,$txend,$chrom,$strand)=$qh_getTxCord->fetchrow_array();$key="$chrom-$strand-$txstart-$txend";
           ${$rsd_hash_ref}{"$organism_version_id-$gene_pred_id"}{"$key"}{"$exonStarts:$exonEnds"}=$tx_name;}
   }}}
}
################## Process starts ##################################################
$qh_getCurrentMaxRsdNumber->execute(); my $currentMaxRsdNumber=0; ##Set global rsd_ numbers counter 
if($qh_getCurrentMaxRsdNumber->rows){$currentMaxRsdNumber=$qh_getCurrentMaxRsdNumber->fetchrow_array();}
$currentMaxRsdNumber=$currentMaxRsdNumber<=0?0:$currentMaxRsdNumber; $tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
open(LOG,">novel_transcript.log");
print LOG "Max tx rsd number=$currentMaxRsdNumber Program Started $mon $mday, $yday @ $hour:$min:$sec\n";
@files=`ls $opt_d`; 
if($?!=0){
   print LOG "Failed to execute the command 'ls $opt_d': $!\n"; close(LOG);exit;
}
exit;
%novelMap=(); $tx_accession_prefix="rsd_"; %strainFiles=();$local_counter=$currentMaxRsdNumber;
## index transcript files by strain and sample
foreach $file(@files){chomp($file); 
   ($orgv,$source,$sample,$desc,$suffix)=split(/\./,$file);
   if(!($suffix)){$suffix=$desc;$desc="";}next if(!($source=~/Gene/));
   $strainFiles{"$source"}{"$sample"}=$file;$orgversion=$orgv;
 }
#Process data by strain
while(($prediction,$samplefiles)=each(%strainFiles)){
    $rnaseqTxfile="$cwdir/tempfile-$prediction-rnaseq_transcript.txt";$annotFile="$cwdir/tempfile-$prediction.txt";
    open(MAIN,">$annotFile");open(db,">$rnaseqTxfile");
    if(!MAIN || !db){print LOG "Could not open temp files for $prediction -- $!\n"; next;}
    $orgversion="";$prediction_id=0;$orgv_id=0;
    print MAIN "name\tchrom\tstrand\ttxStart\ttxEnd\texonCount\texonStarts\texonEnds\tname2\n";
    print db "organism_version_id\tgene_prediction_id\ttranscript_name\tsample_id\n";
    $qh_getPredId->execute($prediction);($prediction_id)=$qh_getPredId->fetchrow_array();
    if($prediction_id<=0){print LOG "$file  -- Bad gene annotation source name -- $prediction was not found in db\n"; next;}
    while(($sample,$file)=each(%{$samplefiles})){ 
         ($orgversion,$prediction,$sample,$desc,$suffix)=split(/\./,$file);
         if(!$suffix){ print LOG "$file -- Bad filename format -\n"; next;}#enforcing the file name format
         $qh_getSampleid->execute($sample);
         if($qh_getSampleid->rows<=0){$qh_insertSampleid->execute($sample,$desc);$qh_getSampleid->execute($sample);}
         ($sample_id)=$qh_getSampleid->fetchrow_array();$qh_getOrgV->execute($orgversion);($orgv_id)=$qh_getOrgV->fetchrow_array(); 
         if($orgv_id<=0){print LOG "$file -- Bad organismVersion name - $orgversion not found in db\n"; next;}
         if($currentMaxRsdNumber>0){  #check if you need to load the rsd map for this strain
            if(!exists($novelMap{"$orgv_id-$prediction_id"})){loadCurrentTxRsd(\%novelMap,$orgv_id,$prediction_id);}}
         open(IN,"$opt_d/$file"); if(!IN){print LOG "$opt_d/$file -- Not Loaded -- $!\n"; next;}
         #index file header field names
         $header=<IN>;chomp($header);$line="";@fields=split("\t",$header);$index=0;$headerfieldCount=scalar(@fields);
         while(@fields>0){$field=shift(@fields);chomp($field);
              if($line eq ""){$line="$field>$index";}else{$line.=",$field>$index";} ++$index;}
         $chromindex=-1;$strandindex=-1;$txStartindex=-1;$txEndindex=-1;$geneIndex=-1;$exonCountindex=-1;
         $exonStartsindex=-1;$exonEndsindex=-1;
         if($line=~/,?chrom>(\d+)/i){$chromindex=$1;}if($line=~/,strand>(\d+)/i){$strandindex=$1;}
         if($line=~/,txStart>(\d+)/i){$txStartindex=$1;}if($line=~/,txEnd>(\d+)/i){$txEndindex=$1;}
         if($line=~/,exonCount>(\d+)/i){$exonCountindex=$1;}if($line=~/,name2>(\d+)/i){$geneIndex=$1;}
         if($line=~/,exonStarts>(\d+)/i){$exonStartsindex=$1;}if($line=~/,exonEnds>(\d+)/i){$exonEndsindex=$1;}
         if($chromindex==-1||$strandindex==-1||$txStartindex==-1||$txEndindex==-1||
             $geneIndex==-1||$exonStartsindex==-1||$exonEndsindex==-1){
            print LOG "$file -- Bad header line format for a genePrediction file format-\n"; next;
         }#Index transcriptome
         #Indexes are created to detect identical transcripts across organismVersion strain/samples --
         @rowIndexmap=(); %txStartmap=();
         while(<IN>){chomp($_);split /\t/;$txstart=$_[$txStartindex];$gene="";
              $gene=$_[$geneIndex] if(scalar(@_)==$headerfieldCount;
              $rowIndexmap[$.]=join ":",($_[$exonStartsindex],$_[$exonEndsindex]);
              if(!exists($txStartmap{"$txstart"})){
                 $txStartmap{"$txstart"}=join ":",($_[$txEndindex],$_[$exonCountindex],$_[$chromindex],$_[$strandindex],$.,$gene);}
              else{
                 $txStartmap{"$txstart"}.=";".join ":",($_[$txEndindex],$_[$exonCountindex],$_[$chromindex],$_[$strandindex],$.,$gene);}
         } #end of index creation
         while(($txstart,$line)=each(%txStartmap)){ split /;/,$line;
               while(@_){ $newline=shift(@_);
                    ($txend,$excount,$chrom,$strand,$linenumber,$gene)=split(":", $newline);
                     $chrom=~s/chr//i;$chrom=~s/ch//i;$chrom=uc($chrom);$chrom="M" if($chrom eq "MT");
                     $exons=$rowIndexmap[$linenumber];$key="$chrom-$strand-$txstart-$txend";$accession="";
                     #identical transcripts across samples of same strain share the same accession number
                     if(!exists($novelMap{"$orgv_id-$prediction_id"}{"$key"}{"$exons"})){ 
                        ++$local_counter;$accession="$tx_accession_prefix$local_counter";
                        $novelMap{"$orgv_id-$prediction_id"}{"$key"}{"$exons"}="$accession";
                      }else{$accession=$novelMap{"$orgv_id-$prediction_id"}{"$key"}{"$exons"};}
                      print db "$orgv_id\t$prediction_id\t$accession\t$sample_id\n";($exstarts,$exends)=split(":",$exons);
                      print MAIN "$accession\t$chrom\t$strand\t$txstart\t$txend\t$excount\t$exstarts\t$exends\t$gene\n";   
       }} #end of while(($txstart,$line)=each(%txStartmap)){  
   } #while(($sample,$file)=each(%{$samplefiles}))
   close(db); close(MAIN); #now load rnaseq_transcript- rsd_ number to sample id mapping
   $linecount = `wc -l $rnaseqTxfile`;chomp($linecount); 
   if($linecount=~/^(\d+)\s$rnaseqTxfile/){$linecount=$1;
      if($linecount>0){--$linecount;
         $qh_loadRnaSeqTx->execute($rnaseqTxfile);$qh_getRowCount->execute($orgv_id,$prediction_id);
         ($rowCount)=$qh_getRowCount->fetchrow_array(); print LOG " $rnaseqTxfile\t$linecount\t$rowCount\n";   
       } 
   } #load annotations
   $cmd="perl $cwdir/load_annotations.pl -d $cwdir -f $annotFile -a $prediction -v $orgversion -s 0 -e 1";
   system($cmd);system("rm tempfile-*");
  #process next strain
 }# end of 
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "Program ended Total= $total -- $mon $mday, $yday @ $hour:$min:$sec\n";
close(LOG);



