#!/usr/bin/perl

#############################################################################
# This script porcesses novel transcripts data to
# 1. assign/generate rsd [rna seq discovery] numbers 
# 2. generate genePrediction file and load novel transcripts in db
# 3. generate rnaseq_transcript file and load transcript-strain-sample mapping in db
# 4. generate aggregate transcripts file load novel_genes regions mapping in db 

#
# Note: 
# Transcript rsd_ numbers are organism version-strain specific
# Aggregate transcript regions rsd_g numbers are organism version-strain specific
#
# Aggregate transcript regions are generated after all samples data of a given strain
#  have been loaded into the database 
# 
# The process:
#  Every novel transcript is given a rsd_ accession id (rsd_[number]) and every novel gene is given an rsd_g accession id (rsd_g[number])
#  where number is a positive integer value generated when a novel transcript is detected.
#  The following steps are taken when processing novel transcripts:
#  1. Get the database current max of rsd_[number] and current max of rsd_g[number] 
#     (the number part of the accession id) into a global variable(maxNumber)
#  2. Create and load a data structure (dbrsdmap)to store database current novel transcripts by strain
#     and the associated rsd numbers --> this is a dictionary
#    loadCurrentTxRsd($organism_version_id,$gene_pred_id) -> 
#         ${$rsd_hash_ref}{"$organism_version_id-$gene_pred_id"}{"$key"}{"$exonStarts:$exonEnds"}=$tx_name;
#         where  $key="$chrom-$strand-$txstart-$txend";
#     --> this is empty if there is no novel transcripts in db
#  3. 
#  4. Strore sample files by strain : $strainFiles{"$strain"}{"$sample"}=$file
#  5. Foreach $strain:
#    i. Index the new rnaSeq file into data structures to facilitate detection of identical transcripts   
#    ii. foreach novel transcript (txNew) from generated data structure, do:
#       a. assign the rsd number
#          a.1 if(Exists(dbrsdmap{txNew})) -> thisAccession=dbrsdmap{txNew}
#          a.2 else{++maxNumber;  thisAccession=_rsd_[maxNumber];dbrsdmap{txNew}=thisAccession
#       b. generate organismversion-strain-sampleID-thisAccession record 
#       c. generate genePrediction record 
# Novel genes are as organismVersion and strain specific generated as follow:
# 1. sort transcript by chr, strand,start,end
# 2. generate agregate transcript with corresponding transcript ids list
# 3. create/asign a rsd_g number to each uniq aggregate transcript
# 4. generate a db record
#

#  
##############################################################################
use DBI;use vars qw ($opt_h $opt_d); use Getopt::Std;
use Time::localtime; 

use Cwd;
$dbname="graber_transcriptdb"; $host="harlequin"; $user="lnh";$pass="lucie98";
$dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
getopts('hd:');
if(($opt_h)|| !($opt_d)){
    print <<HELP;
    **************************************

    Script usage: ./parseAndLoad_NovelTranscripts.pl -d filesDir

    where filesDir is a directory where novel transcripts files are stored.
    These files are in genePred format and file names follow the specified naming convention:
    organismVersion.annotationSource.sampleName.tissue.txt -
    example: mm10.GRCm38-ensGene68.KO10123.testes.txt
 
    example:  ./parseAndLoad_NovelTranscripts.pl -d CC_to_database

HELP
exit;
}
my $dir = getcwd; chomp($dir);#Current working directory
#get the current max rsd_number for novel transcripts
$query="select MAX(cast(substring(transcript_name from locate('rsd_',transcript_name)+4) as unsigned)) ";
$query.=" from transcript_by_annotation WHERE transcript_name like 'rsd_%'";
$qh_getCurrentMaxRsdNumber= $dbh->prepare($query);

#get the current max rsd_number for novel genes
$query="select MAX(cast(substring(novel_gene_name from locate('rsd_g',novel_gene_name)+5) as unsigned)) ";
$query.=" from novel_genes WHERE novel_gene_name like 'rsd_g%'";
$qh_getCurrentGeneMaxRsdNumber= $dbh->prepare($query);
$qh_getnovelTxid=$dbh->prepare("select transcript_name from transcript_by_annotation where transcript_id=? and gene_prediction_id=?");
my $qh_insertNewdb=$dbh->prepare("Insert into gene_prediction(gene_prediction_name) values(?)");
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
my $qh_getOrgV = $dbh->prepare("select organism_version_id from organism_version where ucsc_db=?");
my $qh_getSampleid=$dbh->prepare("select sample_id from rnaseq_sample where sample_name=?");
my $qh_insertSampleid=$dbh->prepare("insert into rnaseq_sample(sample_name,sample_desc) values(?,?)");
################################################
#GET NOVEL transcripts accession ids this section is used in loadCurrentTxRsd()
$query="select distinct transcript_id,exon_start,exon_end";
$query.=" from transcript_exon t, exon e where transcript_name =? and t.organism_version_id=?  and gene_prediction_id=?";
$query.=" and t.organism_version_id=e.organism_version_id ";
$query.=" and t.exon_id=e.exon_id order by exon_start ";
$qh_getexons=$dbh->prepare($query);
$query="select distinct transcript_name from transcript_by_annotation
        where organism_version_id=? and transcript_name like 'rsd_%' and gene_prediction_id=? ";
$qh_getTxRsds=$dbh->prepare($query);
$query="select tx_start,tx_end,chromosome_name, strand from transcript t, chromosome c where transcript_id=? ";
$query.=" and organism_version_id=? and c.chromosome_id=t.chromosome_id ";
$qh_getTxCord=$dbh->prepare($query);
###########################################
#load rnaseq_transcript table
$qh_deleternaseq=$dbh->prepare("delete from rnaseq_transcript where organism_version_id=? and gene_prediction_id=?");
$qh_loadRnaSeqTx=$dbh->prepare("load data local infile ? into table rnaseq_transcript ignore 1 lines");
$query="select count(*) from rnaseq_transcript where organism_version_id=? and gene_prediction_id=?";
$qh_getRowCount=$dbh->prepare($query);
#load novel_genes table
$qh_getCurrentNovelGene=$dbh->prepare("select distinct novel_gene_name from novel_genes where organism_version_id=? and gene_prediction_id=?");
$query="select min(tx_start),max(tx_end) from transcript t, transcript_by_annotation ta, novel_genes g 
        where novel_gene_name=? and g.transcript_name=ta.transcript_name and ta.transcript_id=t.transcript_id";
$qh_getminAndMaxCoord=$dbh->prepare($query);
$query="select ta.transcript_id from transcript_by_annotation ta, novel_genes g ";
$query.=" where novel_gene_name=? and g.transcript_name=ta.transcript_name";
$qh_novelTxlist=$dbh->prepare($query);
$query="select distinct chromosome_id, strand from transcript t,transcript_by_annotation ta, novel_genes g ";
$query.="where novel_gene_name=? and g.transcript_name=ta.transcript_name and ta.transcript_id=t.transcript_id";
$qh_getGeneChrom=$dbh->prepare($query);
## novel gene load
$qh_deleteNovel=$dbh->prepare("delete from novel_genes where organism_version_id=? and gene_prediction_id=?");
$qh_loadNovelGenes=$dbh->prepare("load data local infile ? into table novel_genes ignore 1 lines");
$query="select count(*) from novel_genes where organism_version_id=? and gene_prediction_id=?";
$qh_getGeneRowCount=$dbh->prepare($query);
$qh_deleteNovelSyn=$dbh->prepare("delete from novel_genes_synonyms where organism_version_id=? and gene_prediction_id=?");
$qh_loadNovelGeneSyn=$dbh->prepare("load data local infile ? into table novel_genes_synonyms ignore 1 lines");
$query="select count(*) from novel_genes_synonyms where organism_version_id=? and gene_prediction_id=?";
$qh_getGeneRowCountSyn=$dbh->prepare($query);
$query="delete from gene_by_annotation where organism_version_id=? and gene_prediction_id=? and gene_name ='--'"; 
$qh_deleteGenebyannot=$dbh->prepare($query);
$qh_drop=$dbh->prepare("drop temporary table if exists gene_by_annotation_temp");

$query="create temporary table gene_by_annotation_temp(
         organism_version_id tinyint default 0,gene_prediction_id tinyint default 0,
        transcript_name varchar(255),transcript_id int default 0,gene_name varchar(255),
        index(organism_version_id),index(gene_prediction_id),index(transcript_name))";
$qh_createtemp=$dbh->prepare($query);
$qh_loadGenebyannot=$dbh->prepare("load data local infile ? into table gene_by_annotation_temp ignore 1 lines");
$qh_getGenebyannotRowCount=$dbh->prepare("select count(*) from gene_by_annotation_temp where organism_version_id=? and gene_prediction_id=?"); 
$query="update gene_by_annotation_temp a, transcript_by_annotation t set a.transcript_id=t.transcript_id 
        where a.transcript_name=t.transcript_name and a.gene_prediction_id=t.gene_prediction_id";
$qh_update=$dbh->prepare($query);
$qh_insert=$dbh->prepare("insert ignore into gene_by_annotation select organism_version_id,gene_name,transcript_id,gene_prediction_id from gene_by_annotation_temp");
$query="select distinct chromosome_id,strand,tx_start,tx_end,ga.transcript_id from transcript t, transcript_by_annotation ga";
$query.=" where gene_prediction_id=? and ga.organism_version_id=? and ga.transcript_name like 'rsd_%' ";
$query.=" and ga.organism_version_id=t.organism_version_id ";
$query.=" and ga.transcript_id=t.transcript_id  order by chromosome_id,strand,tx_start ";
$qh_getRcoord=$dbh->prepare($query);
$query="select gene_name from gene_by_annotation where transcript_id=? and gene_name not like 'rsd_g%' and gene_prediction_id=?";
$qh_getThisgene=$dbh->prepare($query);
## Binary search implementation. on a sorted sequencial array, getting the min and max are O(1) constant 
## I can use this for sorted array, sorted hash keys (cpp map)
sub binarySearch{ ($array_ref,$token,$index)=@_;
 $min=0;$max=@{$array_ref}-1; $value=-1; ##implementing binary search
 while($min<=$max){  $mid=$min+ceil(($max-$min)/2);
  if($keys[$min]==$token){$value=@{$array_ref}[$min];$$index=$min; last;}
  elsif($keys[$max]==$token){$value=@{$array_ref}[$max];$$index=$max;last;}
  elsif($keys[$mid]==$token){$value=@{$array_ref}[$mid];$$index=$mid;last;}
  else{ if($keys[$mid]<$token){$min=$mid+1;}else{$max=$mid-1;}}#compute the new range to search
 } return $value
}

#######################################################################
# loadCurrentRsd: index current novel transcripts into a data structure
# in memory to facilitate tx-rsd mapping - adding the gene_prediction id
# to only return strain sepecific data
#
sub loadCurrentTxRsd{
  my($rsd_hash_ref,$organism_version_id,$gene_pred_id)=@_; $qh_getTxRsds->execute($organism_version_id,$gene_pred_id);
  if($qh_getTxRsds->rows>0){
    while(($tx_name)=$qh_getTxRsds->fetchrow_array()){
       $qh_getexons->execute($tx_name,$organism_version_id,$gene_pred_id);$exonStarts=""; $exonEnds="";$tx_id=0;
       while(($transcript_id,$exon_start,$exon_end)=$qh_getexons->fetchrow_array()){
            if($exonStarts eq ""){
               $exonStarts="$exon_start"; $exonEnds="$exon_end";$tx_id=$transcript_id;}
            else{$exonStarts.=",$exon_start"; $exonEnds.=",$exon_end";}
       }
       if($tx_id>0){
         $qh_getTxCord->execute($tx_id,$organism_version_id);
         ($txstart,$txend,$chrom,$strand)=$qh_getTxCord->fetchrow_array();
          $key="$chrom-$strand-$txstart-$txend";
          ${$rsd_hash_ref}{"$organism_version_id-$gene_pred_id"}{"$key"}{"$exonStarts:$exonEnds"}=$tx_name;
       }
    }
  }
}
#######################################################################
# generateCurrentTxAggregates: Given an organism version and gene prediction,
# this function gets all current aggregate transcript regions sorted by chrom-strand-start,
# aggregate transcript regions  are indexed by "$chromosome_id:$strand:$min-$max"}="$gene_name-$txids";
#
sub generateCurrentTxAggregates{
  my($rsd_genehash_ref,$organism_version_id,$gene_pred_id)=@_; 
  $qh_getCurrentNovelGene->execute($organism_version_id,$gene_pred_id)or die mysql_error();
  if($qh_getCurrentNovelGene->rows>0){
     while(($gene_name)=$qh_getCurrentNovelGene->fetchrow_array()){
         $min=0;$max=0;$txids="";$strand="";$chom="";
         $qh_getminAndMaxCoord->execute($gene_name);($min,$max)=$qh_getminAndMaxCoord->fetchrow_array();
         $qh_getGeneChrom->execute($gene_name);($chrom,$strand)=$qh_getGeneChrom->fetchrow_array();
         $qh_novelTxlist->execute($gene_name); 
         while(($tx_id)= $qh_novelTxlist->fetchrow_array()){$txids=($txids eq "")?"$tx_id":"$txids,$tx_id";}
         ${$rsd_genehash_ref}{"$chrom:$strand:$min-$max"}="$gene_name-$txids";
     }
  }
}
##################################################################################################
# generates aggregate transcript regions for the specified organism version and strain
# Each region has coordinates and list of transcript ids that overlap the region
# new aggregate regions are indexed by $chromosome_id:$strand:$current_start
# coordinates are sorted by chrom-strand-start -> the start of the current
#  item on the list is always>= the start of the previous item
#
# Process regions be chrom-strand
#
sub loadNovelG{
 ($region_map_ref,$org_vid,$prediction_id)=@_;
 $qh_getRcoord->execute($prediction_id,$org_vid);%chromMap=();$count=0;
 while(($chromosome_id,$strand,$tx_start,$tx_end,$transcript_id)=$qh_getRcoord->fetchrow_array()){
   ##I need to get rid of those tx that are on known genes
   $qh_getThisgene->execute($transcript_id,$prediction_id);
   if($qh_getThisgene->rows<=0){ #this 
     $key="$chromosome_id:$strand"; $chromMap{"$key"}{$tx_start}{$tx_end}=$transcript_id;++$count;
   }else{ $found=0; 
       while(($gene_name)=$qh_getThisgene->fetchrow_array()){$found=1 if($gene_name=~/\w+|\d+/);}
       if(!$found){
          $key="$chromosome_id:$strand"; $chromMap{"$key"}{$tx_start}{$tx_end}=$transcript_id;++$count;
       }
   }
 } #now process unique regions
 print "A total of $count novel transcrips on novel genes for $org_vid,$prediction_id\n";
 for my $chromStrand(keys(%chromMap)){#regions are sorted by tx_start
      $current_start=0;$current_end=0;$tx_list="";$more=0;
      for my $tx_start(sort keys(%{$chromMap{$chromStrand}})){
          $tx_id="";@ends=sort keys(%{$chromMap{"$chromStrand"}{$tx_start}});
          $end_index=@ends;--$end_index;$tx_end=$ends[$end_index];
          foreach $tend(@ends){
               $txid=$chromMap{"$chromStrand"}{$tx_start}{$tend}; 
               $qh_getnovelTxid->execute($txid,$prediction_id) or die mysql_error();
               ($accession)=$qh_getnovelTxid->fetchrow_array();
               $tx_id=($tx_id eq "")?$accession:"$tx_id,$accession";
           }
          ##
         if($current_start==0){$current_start=$tx_start;$current_end=$tx_end;$tx_list="$tx_id";}
         else{
            if(($tx_start>$current_start)&&($tx_start<$current_end)){
               if($current_end<$tx_end){$current_end=$tx_end;}  #adjust end coordinate when needed
               $tx_list.=",$tx_id";$more=1;
            }else{  #this is a new genomic region, there is no overlap with previous transcript region
                 ${$region_map_ref}{"$chromStrand:$current_start-$current_end"}{"cord"}=$current_end;
                 ${$region_map_ref}{"$chromStrand:$current_start-$current_end"}{"tx"}=$tx_list;
                 $current_start=$tx_start;$current_end=$tx_end;$tx_list="$tx_id";$more=1;
             }
         }
      }
      ##store last
      if($more){
        ${$region_map_ref}{"$chromStrand:$current_start-$current_end"}{"cord"}=$current_end;
        ${$region_map_ref}{"$chromStrand:$current_start-$current_end"}{"tx"}=$tx_list;
     }   
  }
}
###################################################################################
# returns list of genes(rsd_g) that overlap the current region of interest
# Called when the new aggregate region not found in previousely defined aggregate regions. 
#This happens in the following cases:
# 1. new transcripts,whith coordinates lower or higher than a current novel gene, 
#    overlap a current novel gene. In this case, the region will be assigned the rsd_g number
#    of the current gene. History will store the previous definition of this gene 
#    in the novel_genes_synonyms table
#
# 2. The novel aggregate region span two or more current aggregate transcripts
# In this case, the region will be assigned the rsd_g number of one of the current aggregate transcripts.
#    History will store the previous definitions of this region in the novel_genes_synonyms table
#
# Gene discovery should be run each time we load a new gene annotations set
# 
###################################################################################
sub getOverlapGenes{
  ($chromosome_id,$strand,$start,$end,$gene_map_ref,$result_map_ref)=@_;
  while(($key,$data)=each(%{$gene_map_ref})){
        ($chrom_id,$tstrand,$range)= split(":",$key);
        next if(($chromosome_id!=$chrom_id)&&($tstrand ne $strand));
        ($min,$max)=split("-",$range);
        ($gene_name,$txids)=split("-",${$gene_map_ref}{"$key"});
        if((($min<=$start)&&($max>=$start))||(($min<=$end)&&($max>=$end))){
               ${$result_map_ref}{"$gene_name"}=$txids;
        }
  }
}
##set the ids counter start for both novel genes and novel transcripts
$qh_getCurrentMaxRsdNumber->execute();$qh_getCurrentGeneMaxRsdNumber->execute();
my $currentMaxRsdNumber=0;my $GeneMaxRsdNumber=0;
if($qh_getCurrentMaxRsdNumber->rows){$currentMaxRsdNumber=$qh_getCurrentMaxRsdNumber->fetchrow_array();}
if($qh_getCurrentGeneMaxRsdNumber->rows){$GeneMaxRsdNumber=$qh_getCurrentGeneMaxRsdNumber->fetchrow_array();}
$GeneMaxRsdNumber=$GeneMaxRsdNumber<=0?0:$GeneMaxRsdNumber;
$currentMaxRsdNumber=$currentMaxRsdNumber<=0?0:$currentMaxRsdNumber;
###
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
$token=50000;$index=-1;
#$value=binarySearch(\@keys,$token,\$index); print "Found $token - $value at position $index\n ";
@files=`ls $opt_d`; open(LOG,">novel_transcript.log");
print LOG "Max tx rsd number=$currentMaxRsdNumber Program Started $mon $mday, $yday @ $hour:$min:$sec\n";
print  "Max tx rsd_ number=$currentMaxRsdNumber - Max gene rsd_g $GeneMaxRsdNumber - Program Started $mon $mday, $yday @ $hour:$min:$sec\n";
%novelMap=(); $tx_accession_prefix="rsd_"; %strainFiles=();$local_counter=$currentMaxRsdNumber;$local_gcounter=$GeneMaxRsdNumber;
#index files first by strains then by samples
## Assumptions as decided by Nazira and I:
## 1. Every filename has the format: organismVersion.annotationSource.sampleId.(.+)?.txt  example: mm9.mm9-WSBGene.WSB_11304.liver.txt
##    where (.+)? is the sample brief summary if exists or "unknown"
## 2. annotationSource name should be consistent with what's in our gene_prediction table in the database
##    it is strain based
## 3. organismVersion name should be consistent with what's in our organism_version table in the database
## 4. sample tissue if exists or "unknown"
## Example: mm10.GRCm38-ensGene68.KO10123.testes.txt
#%strainFiles=("37"=>"61","26"=>"61","28"=>"61","40"=>"61","39"=>"61","38"=>"61","33"=>"61","19"=>"61","41"=>"93");
foreach $file(@files){chomp($file); 
   ($orgv,$source,$sample,$desc,$suffix)=split(/\./,$file);
  # print "Processing $source \n"; 
   if(!($suffix)){$suffix=$desc;$desc="";}
   next if(!($source=~/Gene/));
   $strainFiles{"$source"}{"$sample"}=$file;$orgversion=$orgv;
   #print "$orgv,$source,$sample,$desc,$suffix\n";
 }
print "A total of ".keys(%strainFiles)."\n";
#exit;
while(($prediction,$samplefiles)=each(%strainFiles)){ #we process data by strain
  #print "Processing $prediction \n"; 
    open(MAIN,">tempfile-$prediction.txt");$orgversion="";%novel_gene_map=();
    print MAIN "name\tchrom\tstrand\ttxStart\ttxEnd\texonCount\texonStarts\texonEnds\tname2\n";
    open(db,">tempfile-$prediction-rnaseq_transcript.txt");
    print db "organism_version_id\tgene_prediction_id\ttranscript_name\tsample_id\n";
    $prediction=~s/B6Gene/ensGene/;$prediction_id=0;$orgv_id=0; 
    $qh_getPredId->execute($prediction);($prediction_id)=$qh_getPredId->fetchrow_array();
    if($prediction_id<=0){print LOG "Bad annotation source naming for $file  -- $prediction was not found in db\n"; next;}
    %lines=();
    # print "Processing $prediction - $prediction_id\n";
    while(($sample,$file)=each(%{$samplefiles})){ 
         ($orgversion,$prediction,$sample,$desc,$suffix)=split(/\./,$file);
         if(!($suffix)){$suffix=$desc;$desc="";}
        if($suffix){
           #get/generate database local ids for organism version,annotation source, and sample name 
           $qh_getSampleid->execute($sample);if($qh_getSampleid->rows<=0){
           $qh_insertSampleid->execute($sample,$desc);$qh_getSampleid->execute($sample);}
           ($sample_id)=$qh_getSampleid->fetchrow_array();
           $qh_getOrgV->execute($orgversion);($orgv_id)=$qh_getOrgV->fetchrow_array(); 
           if($orgv_id<=0){print LOG "Bad organism version naming for $file \n"; next;}
           #check if you need to load the rsd for this organism
           if($currentMaxRsdNumber>0){
              if(!exists($novelMap{"$orgv_id-$prediction_id"})){ loadCurrentTxRsd(\%novelMap,$orgv_id,$prediction_id);}}
          open(IN,"$opt_d/$file") or die "Bad file $!\n";
          if(IN){  
            $header=<IN>;chomp($header);$line="";@fields=split("\t",$header);$index=0;$fieldCount=scalar(@fields);
            while(@fields>0){$field=shift(@fields);chomp($field);
                if($line eq ""){$line="$field>$index";}else{$line.=",$field>$index";} ++$index;}
            $chromindex=-1;$strandindex=-1;$txStartindex=-1;$txEndindex=-1;$geneIndex=-1;
            $exonCountindex=-1;$exonStartsindex=-1;$exonEndsindex=-1;
            if($line=~/,?chrom>(\d+)/i){$chromindex=$1;}if($line=~/,strand>(\d+)/i){$strandindex=$1;}
            if($line=~/,txStart>(\d+)/i){$txStartindex=$1;}if($line=~/,txEnd>(\d+)/i){$txEndindex=$1;}
            if($line=~/,exonCount>(\d+)/i){$exonCountindex=$1;}
            if($line=~/,exonStarts>(\d+)/i){$exonStartsindex=$1;}if($line=~/,exonEnds>(\d+)/i){$exonEndsindex=$1;}
            if($line=~/,name2>(\d+)/i){$geneIndex=$1;}
            @rowIndexmap=(); %txStartmap=();
            #create indexes- generate unique transcript regions;
            #Indexes are created to detect identical transcripts across organismVersion strain/samples --
            while(<IN>){ chomp($_);split /\t/;$txstart=$_[$txStartindex];$gene="";
                $gene=$_[$geneIndex] if(scalar(@_)==$fieldCount);
                $rowIndexmap[$.]=join ":",($_[$exonStartsindex],$_[$exonEndsindex]);
                if(!exists($txStartmap{"$txstart"})){
                   $txStartmap{"$txstart"}=join ":",($_[$txEndindex],$_[$exonCountindex],$_[$chromindex],$_[$strandindex],$.,$gene);}
                else{
                   $txStartmap{"$txstart"}.=";".join ":",($_[$txEndindex],$_[$exonCountindex],$_[$chromindex],$_[$strandindex],$.,$gene);}
            } #end of index creation
            $total=0;$overlap_total=0;$i=0;#$txCoordHit=0;
            while(($txstart,$line)=each(%txStartmap)){ split /;/,$line;
                  while(@_){ $newline=shift(@_); ++$total; ++$overlap_total if($line=~/;/);
                       ($txend,$excount,$chrom,$strand,$linenumber,$gene)=split(":", $newline);
                       $chrom=~s/chr//i;$chrom=~s/ch//i;$chrom=uc($chrom);$chrom="M" if($chrom eq "MT");
                       $exons=$rowIndexmap[$linenumber];$key="$chrom-$strand-$txstart-$txend";$accession="";
                       #identical transcripts across samples of same strain share the same accession number
                       if(!exists($novelMap{"$orgv_id-$prediction_id"}{"$key"}{"$exons"})){  #create a new accession id
                          ++$local_counter;$accession="$tx_accession_prefix$local_counter";
                          $novelMap{"$orgv_id-$prediction_id"}{"$key"}{"$exons"}="$accession";
                        }else{$accession=$novelMap{"$orgv_id-$prediction_id"}{"$key"}{"$exons"};}
                       print db "$orgv_id\t$prediction_id\t$accession\t$sample_id\n";
                      ($exstarts,$exends)=split(":",$exons);
                      if(!defined($gene)||$gene eq ""){ $gene="--";
                         if(!exists($novel_gene_map{"$chrom:$strand"}{"$txstart"}{"$txend"})){
                            $novel_gene_map{"$chrom:$strand"}{"$txstart"}{"$txend"}=$accession ;
                         }else{
                           $novel_gene_map{"$chrom:$strand"}{"$txstart"}{"$txend"}.=",$accession"; 
                         }
                      }
                       print MAIN "$accession\t$chrom\t$strand\t$txstart\t$txend\t$excount\t$exstarts\t$exends\t$gene\n";
                      
                  }
            } #end of while(($txstart,$line)=each(%txStartmap)){  
         }#end of if(IN)
      } #end of if ($file=~exp)
     } #while(($sample,$file)=each(%{$samplefiles}))
    close(db); close(MAIN);
   #load rnaseq_transcript with these rsd_s
   $rnaseqTxfile="$dir/tempfile-$prediction-rnaseq_transcript.txt";$annotFile="$dir/tempfile-$prediction.txt";
   #print "$annotFile generated\n";
   $linecont = `wc -l $rnaseqTxfile`;chomp($linecont); #load ensembl exon ids
   if($linecont=~/^(\d+)\s$rnaseqTxfile/){$linecont=$1;
      if($linecont>0){--$linecont;
         $qh_loadRnaSeqTx->execute($rnaseqTxfile);
         $qh_getRowCount->execute($orgv_id,$prediction_id);($rowCount)=$qh_getRowCount->fetchrow_array();
         print LOG " $rnaseqTxfile\t$linecont\t$rowCount\n";   
       } 
   }
    #load annotations
   $cmd="perl load_annotations.pl -d $dir -f $annotFile -a $prediction -v $orgversion -s 0 -e 1";
   system($cmd); print " Annotations for $prediction loaded - now generating novel genes\n";

 ####################################################################
  ##now generate novel genes for this strain 
  %rsd_genehash=(); %uniq_region_map=();
  open(MAIN,">tempfile-$orgv_id-$prediction_id.novelGenes.txt");
  print MAIN "organism_version\ttranscript_id\tgene_name\tgene_prediction\n";
  open(db,">tempfile-$orgv_id-$prediction_id-gene_synonyms.txt"); 
  print db "organism_version\tgene_prediction\tnovel_gene_name\tsynonym_gene_name\ttranscript_list\n";
  open(NOV,">tempfile-$orgv_id-$prediction_id.novelGenesBYannot.txt"); print NOV "orgv_id\tprediction_id\ttx_name\ttxid\tnovel_gene\n";
  #load current aggregate regions for this strains
  generateCurrentTxAggregates(\%rsd_genehash,$orgv_id,$prediction_id);
 # print "Total novel genes currently found:".keys(%rsd_genehash)."\n";
  #Generate new aggregate transcript regions
 loadNovelG(\%uniq_region_map,$orgv_id,$prediction_id);
  #  now foreach novel gene in region map; get all the overlap novel genes from current;
  #  if none, assign novel gene accession id, else, get novel gene list then write record 
  #  in both novel gene table and gene synonym table 
  $nof=0;$ovelap=0;$exists=0; #print "We have ".keys(%uniq_region_map)." novel genes\n";
  while(($key,$data)=each(%uniq_region_map)){
    ($chromosome_id,$strand,$range)=split(":",$key);($start,$end)=split("-",$range);
     @tx_list=split(",",$uniq_region_map{"$key"}{"tx"});
     #next if(@tx_list<=0);
     if(exists($rsd_genehash{"$key"})){ #region match existing aggregate tx
        ($novel_gene,$oldTxList)=split("-",$rsd_genehash{"$key"});++$exists;
        while(@tx_list>0){$tx_id=shift(@tx_list);print MAIN "$orgv_id\t$tx_id\t$novel_gene\t$prediction_id\n";
            print NOV "$orgv_id\t$prediction_id\t$tx_id\t0\t$novel_gene\n";
        }
        print db "$orgv_id\t$prediction_id\t$novel_gene\t$novel_gene\t$oldTxList\n";#row in novel_genes_synonyms
     }else{ #not an exact overlap, check if region  overlaps existing novel genes before assigning rsd_g number
        %result_map=(); getOverlapGenes($chromosome_id,$strand,$start,$end,\%rsd_genehash, \%result_map);
        if(keys(%result_map)==0){ #new region, create an id
           ++$nof;
            ++$local_gcounter;$novel_gene="$tx_accession_prefix"."g"."$local_gcounter";
            while(@tx_list>0){$tx_id=shift(@tx_list);print MAIN "$orgv_id\t$tx_id\t$novel_gene\t$prediction_id\n";
               print NOV "$orgv_id\t$prediction_id\t$tx_id\t0\t$novel_gene\n";
              }
        }else{ $novel_gene="";$max=0; %temp_map=();++$ovelap;
             while(($gene,$txlist)=each(%result_map)){##SET the region id to the id of the region with max number of transcripts 
                  @txnlist=split(",",$txlist);
                  if($max==0||($max<@txnlist)){$novel_gene=$gene;$max=@txnlist;} 
                  $temp_map{"$gene"}=$txlist;  
             }
            while(@tx_list>0){$tx_id=shift(@tx_list);
                  print MAIN "$orgv_id\t$tx_id\t$novel_gene\t$prediction_id\n";
                  print NOV "$orgv_id\t$prediction_id\t$tx_id\t0\t$novel_gene\n";}
            while(($gene,$txlist)=each(%temp_map)){ 
               print db "$orgv_id\t$prediction_id\t$novel_gene\t$gene\t$txlist\n";}
        }
     }
   }# end of while(($key,$data)
   close(db); close(MAIN);close(NOV);
   #print " Novel genes for  $prediction loaded\n";
   #now load novel_genes and novel_genes_synonyms tables for this strain
   $genfile="tempfile-$orgv_id-$prediction_id.novelGenes.txt";$geneSynFile="tempfile-$orgv_id-$prediction_id-gene_synonyms.txt";
   $geneByannot="tempfile-$orgv_id-$prediction_id.novelGenesBYannot.txt";
   $linecont = `wc -l $genfile`;chomp($linecont); 
   if($linecont=~/^(\d+)\s$genfile/){$linecont=$1;
      if($linecont>0){--$linecont;
         #$qh_deleteNovel->execute($orgv_id,$prediction_id);
         $qh_loadNovelGenes->execute($genfile);
         $qh_getGeneRowCount->execute($orgv_id,$prediction_id);($rowCount)=$qh_getGeneRowCount->fetchrow_array();
         print LOG " $genfile\t$linecont\t$rowCount\n";   
       }#else{print "There are no new genes to load\n";}
    }
   $linecont = `wc -l $geneSynFile`;chomp($linecont); #load ensembl exon ids
   if($linecont=~/^(\d+)\s$geneSynFile/){$linecont=$1;
      if($linecont>0){--$linecont;
         #$qh_deleteNovelSyn->execute($orgv_id,$prediction_id);
         $qh_loadNovelGeneSyn->execute($geneSynFile);
         $qh_getGeneRowCountSyn->execute($orgv_id,$prediction_id);($rowCount)=$qh_getGeneRowCountSyn->fetchrow_array();
         print LOG " $geneSynFile\t$linecont\t$rowCount\n";  
      }#else{print "There are no new genes to load\n";}
   }
  #load gene_by_annotation table
  $linecont = `wc -l $geneByannot`;chomp($linecont); 
  if($linecont=~/^(\d+)\s$geneByannot/){$linecont=$1;
      if($linecont>0){--$linecont;
         $qh_drop->execute();$qh_createtemp->execute() or die mysql_error();
         $qh_loadGenebyannot->execute($geneByannot);
         $qh_getGenebyannotRowCount->execute($orgv_id,$prediction_id);
         ($rowCount)=$qh_getGenebyannotRowCount->fetchrow_array();
          if($linecont==$rowCount){
             $qh_update->execute();$qh_deleteGenebyannot->execute($orgv_id,$prediction_id);
             $qh_insert->execute();
          }
         print LOG " $geneByannot\t$linecont\t$rowCount\n";   
      }#else{print "There are no new genes to load\n";}
   }
   print "Novel genes loaded\n";
   system("rm tempfile-*");
  #process next strain
 }# end of 
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "Program ended Total= $total -- $mon $mday, $yday @ $hour:$min:$sec\n";



