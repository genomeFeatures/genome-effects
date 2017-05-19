#!/usr/bin/perl

###################################################################################################
## getUniqGenomeFeatureCordinates.pl
#  This script generates unique genome feature cordinates of the specified features
#   for a given organism version and gene prediction source
#  The default gene predictions are knownGene and ensGene
#  The features are stored into the database using 
#  the zero-based cordinates standard for feature starts
#  and 1-based cordinate standard for feature ends
#    
#  
# Output: a text file with the following fields
# 1. Chromosome,2. chromStart,3. chromEnd,4. strand,5. geneName, 
#    6. transcriptId, 7. featureType
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : December 2012

#   Usage: perl getUniqFeatureCordinates.pl -f <output format> -t <feature type> 
#                  -o <output file name> -v <organism version> -a <gene prediction list> -d db_host
#
#   Where: -f is output file format 1->tab-delimitted, 2->commas separated (default tab-delimitted)
#          -o  is the output file name (default standard out)
#          -t  feature type: 
#             ex    -> Exons
#             int   -> Introns  
#             tx    -> Transcripts                     
#          -v  is the ucsc organism version (default mm9)
#          -a  a commas separated list of gene prediction names (default knowGene,ensGene)
#          -l  display the list of current gene predictions for the specified organism version
#
###################################################################################################
use DBI;
use POSIX;

use vars qw ($opt_h $opt_s $opt_o $opt_t $opt_v $opt_a $opt_l $opt_f $opt_d $opt_u);
use Getopt::Std;

getopts('hlso:t:v:a:f:d:u:');
if($opt_h||(!$opt_o&&!$opt_s&&!$opt_l)) {
    print <<HELP;

This script generates unique genome feature cordinates of the specified features for a given organism version and gene prediction source.The default gene predictions are knownGene and ensGene.

Usage: perl getUniqFeatureCordinates.pl -f <output format> -t <feature type> 
                -o <outputDirectory> -v <organism version> -a <gene prediction list>
Arguments:
   -h  displays this help message
   -o  is the output directory
   -t  feature type: 
       ex    -> Exons
       int   -> Introns  
       f_ex  -> First Exons
       l_ex  -> last exons
       f_cex -> First coding exon
       l_cex -> Last coding exon      

  -f is output file format 1->tab-delimitted, 2->commas separated (default tab-delimitted)
  -v  is the ucsc organism version (default mm9)
  -a  a commas separated list of gene prediction names (default knowGene,ensGene)
  -l  display the list of current gene predictions for the specified organism version
  -d  specifies the database server (default harlequin.jax.org)

Example: perl getUniqFeatureCordinates.pl -v mm9 -l
The above will display all the gene predictions for mm9

Example: perl getUniqFeatureCordinates.pl -t tx -o features -v mm9 -a ensGene,knownGene -f 1
The above will use Ensembl  and knownGene annotations to generate all the transcripts
for mm9 and store the result in a tab-delimited file under features. 

Assumptions:
1. all the feature starts (exon,transcript,cds) are stored using UCSC standards  - zero-base
2. all the feature ends (exon,transcript,cds) are stored using UCSC standards    - 1-base

HELP
exit;
}

#set defaults
my $user ='pup';my $pwd  ='puppass';$dbname ='graber_transcriptdb';$host ='harlequin.jax.org';
$flanking=200; $organism_version="mm9";
if($opt_d){$host=$opt_d;}
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
#Steps:
# 1. Get the list of chromosome of the organism version
# 2. for each chrom
#    a. get all the provided transcript definitions on both strands sorted by exon_start then,
#       loop through every strand to filter overlapping transcript(aggregate view) and generate a new list of transcript,
#    b. For query other than tx, for each unique transcript, generate unique aggregate exons,
#       sort exons then compute the cordinates of the requested features.
#       get the list of all the transcripts that contain this feature
#    c. Check the genome strand
#       b.1. positive strand then start counting from the exon with min exon_start 
#       b.2. negative strand, start counting from the exon with max exon_end 
#
######################################################
$query="select distinct chromosome_id,strand,tx_start,tx_end,t.transcript_id from transcript t, transcript_by_annotation ta";
$query.=" where t.organism_version_id=? and t.organism_version_id=ta.organism_version_id ";
$query.=" and t.transcript_id =ta.transcript_id  and gene_prediction_id in (?) order by  chromosome_id,strand,t.tx_start";
$qh_getChromTx= $dbh->prepare($query);

$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$query="select distinct transcript_id from transcript_by_annotation 
       where organism_version_id=? and gene_prediction_id in (?)";
$qh_getTxByAnnot=$dbh->prepare($query);

my $getOverlapE="select exon_start,exon_end from exon where organism_version_id=? ";
   $getOverlapE.=" and chromosome_id=? and strand=? and exon_end >= ? and exon_start<=?  order by exon_start";
my $qh_getOverlapE = $dbh->prepare($getOverlapE)or die "Couldn't prepare statement: ".$dbh->errstr;


my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $qh_getOrgV = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getOrgUcsc="select ucsc_db  from organism_version where organism_version_id=?";
my $qh_getOrgUcsc = $dbh->prepare($getOrgUcsc)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getAllOrgV="select distinct o.organism_name,ucsc_db as version from organism o,organism_version ov
                 where o.organism_id=ov.organism_id order by organism_name";
my $qh_getAllOrgV = $dbh->prepare($getAllOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;
my $getGenePred="select gene_prediction_name from gene_prediction_by_organism op, gene_prediction p
                   where op.organism_version_id =? and op.gene_prediction_id=p.gene_prediction_id 
                   and gene_prediction_name not in ('all_est','all_mrna','estOrientInfo')";
my $qh_getGenePred = $dbh->prepare($getGenePred)or die "Couldn't prepare statement: ".$dbh->errstr;
my $getTxGene="select distinct gene_name from gene_by_annotation 
                where organism_version_id=? and transcript_id in(?) and gene_prediction_id in (?)";
my $qh_getTxGene = $dbh->prepare($getTxGene)or die "Couldn't prepare statement: ".$dbh->errstr;
my $getTxAccession="select distinct transcript_name from transcript_by_annotation 
                    where organism_version_id=? and transcript_id in(?) and gene_prediction_id in (?)";
my $qh_getTxAccession = $dbh->prepare($getTxAccession);

$query="select cdsStart,cdsEnd from transcript_translation ";
$query.=" where transcript_id in (?) and organism_version_id=? and gene_prediction_id in (?) order by cdsStart asc";
my $qh_getCDS=$dbh->prepare($query);


$organism_version=$opt_v if($opt_v);$org_vid=0;$organism_version=~s/\s+//g;
$qh_getOrgV->execute("$organism_version");
if($qh_getOrgV->rows>0){($org_vid)=$qh_getOrgV->fetchrow_array();}
if($opt_s){# distinct o.organism_name,ov.organism_version_id ,ucsc_db as version
  if($opt_l){
    print "$organism_version gene prediction sources\n";
    $qh_getGenePred->execute($org_vid);
    if($qh_getGenePred->rows>0){
       while(($gene_prediction_name)=$qh_getGenePred->fetchrow_array()){print "$gene_prediction_name\n";}}
  }
  else{ 
    $qh_getAllOrgV->execute();
    if($qh_getAllOrgV->rows>0){
      print "Organism\tVersion\n";
      while(($organism,$version)=$qh_getAllOrgV->fetchrow_array()){print "$organism\t$version\n";}}
  }
}
elsif($opt_l){
   print "$organism_version gene prediction sources\n";
   $qh_getGenePred->execute($org_vid);
   if($qh_getGenePred->rows>0){
      while(($gene_prediction_name)=$qh_getGenePred->fetchrow_array()){print "$gene_prediction_name\n";}}
}
else{
 $gene_prediction="";$pred_id=0;%predictions=();
 $gene_prediction=$opt_a;
 @predictions=split(",",$gene_prediction);$gen_prediction="";$prediction_ids="";
 foreach my $pred(@predictions){ #generate gene prediction ids for selected gene predisctions
    $qh_getPredId->execute($pred);
    if($qh_getPredId->rows>0){
     ($pred_id)=$qh_getPredId->fetchrow_array();
      if($prediction_ids eq ""){$prediction_ids="$pred_id";}
      else{$prediction_ids.=",$pred_id";}
    }
 }
 if($prediction_ids eq ""){
  print "*********\nNote: We do not have $organism_version data from the specified annotation ($gene_prediction)\n";
  print " Please run the script with the -l option to display the current annotations list.\n";
  exit(0);
 }
 sub getOverlapExons{
  ($chrom_id,$strand,$chrom_start,$chrom_end,$org_vid)= @_;
   $qh_getOverlapE->execute($org_vid,$chrom_id,$strand,$chrom_start,$chrom_end);
   if($qh_getOverlapE->rows>0){
      $current_start=0;$current_end=0;$ex_start=0;$ex_end=0;$ex_starts="";$ex_ends="";
      while(($ex_start,$ex_end)=$qh_getOverlapE->fetchrow_array()){
             $ex_start+=1;
             if($current_start==0){$current_start=$ex_start;$current_end=$ex_end;}
             else{
                 if(($ex_start>=$current_start)&&($ex_start<=$current_end)){
                     if($current_end<$ex_end){$current_end=$ex_end;}
                 }
                 else{ #this is a new genomic region, there is no overlap with previous exon
                       #$region_map{$current_start}=$current_end;
                     if($ex_starts eq ""){$ex_starts="$current_start";$ex_ends="$current_end";}
                     else{$ex_starts.=",$current_start";$ex_ends.=",$current_end";}
                     $current_start=$ex_start;$current_end=$ex_end;
                 }
             }#now add the last exon
       }
      if($ex_starts eq ""){
         if($current_start>0){$ex_starts="$current_start";$ex_ends="$current_end";}
      }
      else{
        if($current_start>0){$ex_starts.=",$current_start";$ex_ends.=",$current_end";}
      }
   }
   return "$ex_starts:$ex_ends";
 }
 ### get the list of all the transcripts of this organism version 
 $qh_getChromTx->execute($org_vid,$prediction_ids);
 if($qh_getChromTx->rows<=0){
  print "*********\nNote: nothing found for  $organism_version data from the specified annotation ($gene_prediction)\n";
  print " Please run the script with the -s and -l option to display the current annotations list. for $organism_version \n";
  exit(0);
 }
 $file_name="";$opt_o=~s/\/\s*$//;$opt_t=~s/\/\s*$//g;
 if($opt_t=~/ex/){$file_name="$opt_o/$organism_version-$gene_prediction-exons.txt";}
 else{$file_name="$opt_o/$organism_version-$gene_prediction-introns.txt";} $count=0;
 open (TS,">$file_name") or die "$!\n"; $count=0;%region_map=();%tx_map=();$tx_list="";%txmap=();
 ### set up structures to filter duplicates genomic regions
 print TS "chrom\tchromStart\tchromEnd\tstrand\tcdsStart\tcdsEnd\tfeatureStarts\tfeatureEnds\ttranscriptList\tgeneList\n";
 $current_start=0;$current_end=0;$ex_start=0;$ex_end=0;
 while(($chromosome_id,$strand,$tx_start,$tx_end,$tx_id)=$qh_getChromTx->fetchrow_array()){
        $tx_start+=1;
         if($current_start==0){$current_start=$tx_start;$current_end=$tx_end;$tx_list="$tx_id";$txmap{$tx_id}=1;}
         else{
             if(($tx_start>=$current_start)&&($tx_start<=$current_end)){
               if($current_end<$tx_end){$current_end=$tx_end;}
               if(!exists($txmap{$tx_id})){$txmap{$tx_id}=1;$tx_list.=",$tx_id";}
             }
             else{ #this is a new genomic region, there is no overlap with previous transcript region
                $region_map{$chromosome_id}{$strand}{$current_start}{"cord"}=$current_end;
                $region_map{$chromosome_id}{$strand}{$current_start}{"tx"}=$tx_list;++$count;
                $current_start=$tx_start;$current_end=$tx_end;$txmap{$tx_id}=1;$tx_list="$tx_id";
               # last if($count>20);
             }
         }
       
  }
  $region_map{$chromosome_id}{$strand}{$current_start}{"cord"}=$current_end;
  $region_map{$chromosome_id}{$strand}{$current_start}{"tx"}=$tx_list;
  #now loop through the unique transcript list to get the unique exons
 for my $chromosome(sort keys(%region_map)){
     for my $strand (sort keys (%{$region_map{$chromosome}})){
          for my $txstart(sort keys(%{$region_map{$chromosome}{$strand}})){
              $txend=$region_map{$chromosome}{$strand}{$txstart}{"cord"};
              $tx_idlist=$region_map{$chromosome}{$strand}{$txstart}{"tx"};$transcriptList="";
              $qh_getTxAccession->execute($org_vid,$tx_idlist,$prediction_ids);$geneList="";
              $qh_getCDS->execute($tx_idlist,$org_vid,$prediction_ids);$cdsStarts="";$cdsEnds="";
              while(($transcript)=$qh_getTxAccession->fetchrow_array()){
                 if($transcriptList eq ""){$transcriptList="$transcript";}
                 else{$transcriptList.=",$transcript";}}
              while(($cdsStart,$cdsEnd)=$qh_getCDS->fetchrow_array()){
                 $cdsStart+=1; #adjust tx_start to 1-base
                 if($cdsStarts eq ""){$cdsStarts=$cdsStart;$cdsEnds=$cdsEnd;
                 }else{$cdsStarts.=",$cdsStart";$cdsEnds.=",$cdsEnd";}
              }
              $qh_getTxGene->execute($org_vid,$tx_idlist,$prediction_ids);
              while(($gen)=$qh_getTxGene->fetchrow_array()){
                if($geneList eq ""){$geneList=$gen;}else{$geneList.=",$gen";}}
                # else{ #get unique exons from this region
              $ex_starts=""; $ex_ends="";
              ($ex_starts,$ex_ends)= split(":",getOverlapExons($chromosome,$strand,$txstart,$txend,$org_vid));
              if($opt_t=~/ex/){
                 if($ex_starts ne ""){
                    print TS "$chromosome\t$txstart\t$txend\t$strand\t$cdsStarts\t$cdsEnds\t$ex_starts\t$ex_ends\t$transcriptList\t$geneList\n";
                  }
               }
               else{
                   @startsexon=split(",",$ex_starts);@endsexon= split(",",$ex_ends);
                   $i=1;$intron_starts="";$intron_ends="";
                   while($i<@startsexon){
                         $int_start=$endsexon[$i-1]+1;$int_end=$startsexon[$i]-1;
                         if($i==1){$intron_starts="$int_start";$intron_ends="$int_end";}
                         else{$intron_starts.=",$int_start";$intron_ends.=",$int_end";} ++$i;
                   }
                   if($intron_ends ne ""){
            print TS "$chromosome\t$txstart\t$txend\t$strand\t$cdsStarts\t$cdsEnds\t$intron_starts\t$intron_ends\t$transcriptList\t$geneList\n";
                   }
               }
              #}   
          }
      }
   }
 print "$file_name generated\n";
}#end of query by -t option
print"program completed\n";
exit(0);

