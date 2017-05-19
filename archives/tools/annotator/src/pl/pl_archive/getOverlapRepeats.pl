#!/usr/bin/perl

###################################################################################################
## getCordinates.pl
#  This script generates all the exons that overlap repeat region
#   for a given organism version and gene prediction source
#  The default gene prediction is ensGene
#  The features are stored into the database using 
#  the zero-based cordinates standard for feature starts
#  and 1-based cordinate standard for feature ends
#    
#  
# Output: a text file with the following fields
# 1. Transcript,2. repeatid,3. repeattype,4. chromosome,5.strand, 
#    6. repeatStart, 7.repeatEnd,8.exonStart,9.exonEnd
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : February 2013

#   Usage: perl getOverlapRepeat.pl -t <feature type> 
#                  -o <output file name> -v <organism version> -a <gene prediction list>
#
#   Where: 
#          -o  is the output file name (default standard out)
#          -t  feature type: 
#                         1. repeat -> repeat
#                         2. mirna   -> miRNA  
#          -v  is the ucsc organism version (default mm9)
#          -a  a commas separated list of gene prediction names (default knowGene,ensGene)
#          -l  display the list of current gene predictions for the specified organism version
#
###################################################################################################
use DBI;
use POSIX;

use vars qw ($opt_h $opt_o $opt_t $opt_v $opt_a $opt_l);
use Getopt::Std;

getopts('hlo:t:v:a:');
if($opt_h||(!$opt_o&&!$opt_l)) {
    print <<HELP;

This script generates exons that overlap repeat/miRNA regions for
 a given organism version and gene prediction source
The default gene prediction ensGene

#   Usage: perl getOverlapRepeat.pl -t <feature type> 
#                  -o <output file name> -v <organism version> -a <gene prediction list>
#
#   Where: 
#          -o  is the output file name (default standard out)
#          -t  feature type: 
#                         1. tx -> repeat-transcript
#                         2. ex -> repeat-exon
#          -v  is the ucsc organism version (default mm9)
#          -a  a commas separated list of gene prediction names (default knowGene,ensGene)
#          -l  display the list of current gene predictions for the specified organism version

Example: perl getCordinates.pl -v mm9 -l
The above will display all the gene predictions for mm9

Example: perl getOverlapRepeat.pl -t ex -o features -v mm10 -a GRCm38-ensGene
The above will use Ensembl  gene annotations to generate all the exons
for mm10  that overlap repeat regions and store the result in a tab-delimited file under features/. 

Assumptions:
1. all the feature starts (exon,transcript,cds) are stored using UCSC standards  - zero-base
2. all the feature ends (exon,transcript,cds) are stored using UCSC standards    - 1-base

HELP
exit;
}

#set defaults
my $user ='pup';my $pwd  ='puppass';$dbname ='graber_transcriptdb';$host ='harlequin.jax.org';
 $organism_version="mm9";
if($opt_d){$host=$opt_d;}
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
######################################################
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
my $getGenePred="select gene_prediction_name from gene_prediction_by_organism op, gene_prediction p
                   where op.organism_version_id =? and op.gene_prediction_id=p.gene_prediction_id 
                   and gene_prediction_name not in ('all_est','all_mrna','estOrientInfo')";
my $qh_getGenePred = $dbh->prepare($getGenePred)or die "Couldn't prepare statement: ".$dbh->errstr;

$query="select chromosome_id,strand,repeat_feature_id,repeat_name_id,repeat_type_id,repeatChrStart,repeatChrEnd ";
$query.=" from repeats where organism_version_id=? and  chromosome_id=?";
$qh_getRepeat=$dbh->prepare($query);
my $getRepeatType="select repeat_type ,repeat_type_id from repeat_type ";
my $qh_getRepeatType = $dbh->prepare($getRepeatType)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getRepeatName="select repeat_name,repeat_name_id from repeat_name";
my $qh_getRepeatName = $dbh->prepare($getRepeatName)or die "Couldn't prepare statement: ".$dbh->errstr;

  $query=" select distinct a.transcript_name,t.exon_start,t.exon_end from exon t,transcript_exon a ";
  $query.=" where t.organism_version_id=? and t.chromosome_id=? and strand=? and t.exon_end >= ? and t.exon_start<=? ";
  $query.=" and  t.organism_version_id=a.organism_version_id  ";
  $query.=" and  a.gene_prediction_id in(?) and a.exon_id=t.exon_id order by exon_start";
my $qh_getRepExons=$dbh->prepare($query);
  $query=" select distinct a.transcript_name,tx_start,tx_end from transcript t,transcript_exon a";
  $query.=" where organism_version_id=? and chromosome_id=? and strand=? and tx_end >= ? and tx_start<=?  ";
  $query.=" and  t.organism_version_id=a.organism_version_id  ";
  $query.=" and  a.gene_prediction_id in(?) and a.transcript_id=t.transcript_id order by tx_start";
my $qh_getRepTx=$dbh->prepare($query);

my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $qh_getOrgV = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;
my $getChr="select distinct t.chromosome_id, c.chromosome_name from transcript t, chromosome c 
             where organism_version_id=? and t.chromosome_id=c.chromosome_id";
my $qh_getChr = $dbh->prepare($getChr)or die "Couldn't prepare statement: ".$dbh->errstr;


$organism_version=$opt_v if($opt_v);$org_vid=0;$organism_version=~s/\s+//g;
$qh_getOrgV->execute("$organism_version");
if($qh_getOrgV->rows>0){($org_vid)=$qh_getOrgV->fetchrow_array();}
%chrom_map=();$qh_getChr->execute($org_vid);%repeatType=();%repeatName=();
if($qh_getChr->rows>0){
  while(($chr_id,$chr_name)=$qh_getChr->fetchrow_array()){$chrom_map{$chr_id}=$chr_name;}
}
$qh_getRepeatType->execute();$qh_getRepeatName->execute();
if($qh_getRepeatType->rows>0){
  while(($rep_type,$type_id)=$qh_getRepeatType->fetchrow_array()){$repeatType{$type_id}=$rep_type;}
}
if($qh_getRepeatName->rows>0){
  while(($rep_name,$name_id)=$qh_getRepeatName->fetchrow_array()){$repeatName{$name_id}=$rep_name;}
}
if($opt_l){
   print "$organism_version gene prediction sources\n";
   $qh_getGenePred->execute($org_vid);
   if($qh_getGenePred->rows>0){
      while(($gene_prediction_name)=$qh_getGenePred->fetchrow_array()){print "$gene_prediction_name\n";}}
}
else{
 $gene_prediction="";$pred_id=0;%predictions=();$prediction_ids="";
 if($opt_a){$gene_prediction=$opt_a;
    @predictions=split(",",$gene_prediction);
   foreach my $pred(@predictions){ #generate gene prediction ids for selected gene predisctions
     $qh_getPredId->execute($pred);
     if($qh_getPredId->rows>0){($pred_id)=$qh_getPredId->fetchrow_array();
       if($prediction_ids eq ""){$prediction_ids="$pred_id";}else{$prediction_ids.=",$pred_id";}
     }
   }
 } # default all annotations
 

 $file_name="";$opt_o=~s/\/\s*$//;$opt_t=~s/\/\s*$//g;$feature_type="Transcript";
 if($opt_t=~/tx/){$file_name="$opt_o/$organism_version-$gene_prediction-repeat-transcripts.txt";}
 elsif($opt_t eq "ex"){$file_name="$opt_o/$organism_version-$gene_prediction-repeat-exons.txt";$feature_type="Exon";}
 open (TS,">$file_name") or die "$!\n"; $count=0;
 ### get the list of all the known repeats of this organism version 
 while(($chr_id,$chr_name)=each(%chrom_map)){
  $qh_getRepeat->execute($org_vid,$chr_id);
  if($qh_getRepeat->rows>0){
     $total=$qh_getRepeat->rows;
     print TS "transcript\trepeat_id\trepeat_name\trepeat_type\tchromosome\tstrand\trepeatStart\trepeatEnd\t";
     print TS "featureStart\tfeatureEnd\tfeaturetype\tsource\n";
     while(($chr_id,$strand,$repeat_id,$repeat_name_id,$repeat_type_id,$repeatChrStart,$repeatChrEnd)=$qh_getRepeat->fetchrow_array()){
         $repeat_name=$repeatName{$repeat_name_id};$repeat_type=$repeatType{$repeat_type_id};
         $repeatChrStart+=1; ++$count; print "$count of $total repeats processed\n" if($count%100000==0);
         #get overlap features
         $qh_handler;
         if($opt_t eq "ex"){
            $qh_getRepExons->execute($org_vid,$chr_id,$strand,$repeatChrStart,$repeatChrEnd,$prediction_ids);
            $qh_handler=$qh_getRepExons;
         }
         else{$qh_getRepTx->execute($org_vid,$chr_id,$strand,$repeatChrStart,$repeatChrEnd,$prediction_ids);
              $qh_handler=$qh_getRepTx;
         }
        if($qh_handler->rows>0){ # last if($count>10);
          while(($tx_name,$feature_start,$feature_end)=$qh_handler->fetchrow_array()){
            $feature_start+=1;  
            print  TS  "$tx_name\t$repeat_id\t$repeat_name\t$repeat_type\t$chr_name\t$strand\t";
            print  TS "$repeatChrStart\t$repeatChrEnd\t$feature_start\t$feature_end\t$feature_type\t$organism_version:$gene_prediction\n";
          }
        }
    }
  }
 print "Done with chromosome $chr_name\n";
 }
 print "$file_name generated\n";
}#end of query by -t option
print"program completed\n";
exit(0);

