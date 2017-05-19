#!/usr/bin/perl

use DBI;
use Time::localtime;

#***********************************************************************************************
# Given a transcript name,this script returns a list of all ovrelaping transcripts
# with the associated exon data
#
# Input : a config file containing a list of transcript names (one name per line)
#
# Output: a tab-delimetted file with the following fields:
#         a. chromosome
#         b. strand (+-)
#         c. tx_name
#         d. tx_start
#         e. tx_end
#         f. exon_start
#         g. exon_end
#
#Note : the script assumes that the input transcript list is from Jesse annotations
#       (gene_prediction_id=14)
#     
# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : Aug 2010
#
#Usage : getOverlapTranscript.pl transcript_file
#      
#********************************************************************************************
#

my ($tx_configFile)= @ARGV;
my $gene_prediction_id=14; #Jesse annotations id for the input data
my $organism_version_id=35; #default organism is mouse build 
if(!$tx_configFile){
    print "Usage: getOverlapTranscript.pl transcript_file\n";
}
# Set the path to configuration files
my $data_path="/scratch/data/graber_transcriptdb/overlapping_transcripts";
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
my $processlog="$data_path/overlapping_transcripts"."_"."$mon-$mday-$yday".".txt";
chomp($tx_configFile);

open(IN,"$tx_configFile");
if(IN){
   open(LOG,">$processlog");
   if(!LOG){                                        #if the config file has issues report and quit
      print"Can not open $processlog due to $!\n";
      close(LOG); exit(1);
   }
  
   $dbname ="graber_transcriptdb"; $host="harlequin"; $user="lnh";$pwd="lucie";
   my $dbh = DBI->connect("DBI:mysql:$dbname:$host",$user,$pwd)
                            or die "Can't connect to database!\n";
  
  my $getChr ="select chromosome_name from chromosome where chromosome_id=? ";         # get this chromosome name
  my $getTxId="select transcript_id from transcript_by_annotation ";
     $getTxId.=" where transcript_name=? and gene_prediction_id=?";                    # get this transcript id
  my $getTxCoord="select chromosome_id,strand,tx_start,tx_end from transcript where transcript_id=?";   #get this tx coordinates
     $getTxCoord.=" and organism_version_id=35";
  
  my $getTxName=" select distinct transcript_name from transcript_by_annotation where transcript_id=?";
  my $getOverlap ="select t.transcript_id,t.tx_start,t.tx_end from transcript t";
     #$getOverlap .=" ,transcript_name from transcript t, transcript_by_annotation ta ";
     $getOverlap .=" where t.chromosome_id=? and t.strand=? and t.tx_end >= ? and t.tx_start<=?";
     $getOverlap .="  and t.organism_version_id=$organism_version_id order by tx_start";
     
  my $getTxExons=" select distinct exon_id from transcript_exon where transcript_id=? "; #get exons list for this transcript
  my $getExonCoord="select exon_start,exon_end from exon where exon_id=?";            #get this exon coordinates
  
  my $qh_chrname      = $dbh->prepare($getChr)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_txid         = $dbh->prepare($getTxId)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_TxCoord      = $dbh->prepare($getTxCoord)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_Overlaplist  = $dbh->prepare($getOverlap)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_TxExons      = $dbh->prepare($getTxExons)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_ExonCoord    = $dbh->prepare($getExonCoord)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_Txname       = $dbh->prepare($getTxName)or die "Couldn't prepare statement: " . $dbh->errstr;
  print LOG "Chromosome\tStrand\ttx_name\ttx_start\ttx_end\texon_start\texon_end\n";
  while($tx_name=<IN>){  
      chomp($tx_name);
      my($chr,$chr_id,$tx_id,$transc_name,$strand,$tx_start,$tx_end,@overlap_list,@exon_list,$exon_coord);
      #get this transcript id
      $qh_txid->execute($tx_name,$gene_prediction_id) or die "Can't execute query: " . $dbh->errstr . "\n";
      if(@row=$qh_txid->fetchrow_array()){
         $tx_id=$row[0];
       }
      #get transcript coordinates
      $qh_TxCoord->execute($tx_id) or die "Can't execute query: " . $dbh->errstr . "\n";
      if(@row=$qh_TxCoord->fetchrow_array()){
         $chr_id=$row[0]; $strand=$row[1];
         $tx_start=$row[2]; $tx_end=$row[3];
       }
      # get chromosome name
      $qh_chrname->execute($chr_id) or die "Can't execute query: " . $dbh->errstr . "\n";
      if(@row=$qh_chrname->fetchrow_array()){
         $chr=$row[0];$chr=~s/^\s+//;$chr=~s/\s+$//;
       }
       #get overlap transcripts
      $qh_Overlaplist->execute($chr_id,$strand,$tx_start,$tx_end) or die "Can't execute query: " . $dbh->errstr . "\n";
      while(@row=$qh_Overlaplist->fetchrow_array()){
           $transcript_id=$row[0]; $tx_start=$row[1];$tx_end=$row[2];$transc_name=""; #$transc_name=$row[3];
           #get all the accessessions associated with this tx_id
           $qh_Txname->execute($transcript_id)or die "Can't execute query: " . $dbh->errstr . "\n";
           $j=0;
           while(@row4=$qh_Txname->fetchrow_array()){
              if($j==0){$transc_name.="$row4[0]";++$j;}
              else{$transc_name.=",$row4[0]";}
           }
           #get exon list of this transcript
           $qh_TxExons->execute($transcript_id) or die "Can't execute query: " . $dbh->errstr . "\n";
           while(@row2=$qh_TxExons->fetchrow_array()){
              $exon_id=$row2[0];
              #get exon coordinates
              $qh_ExonCoord->execute($exon_id) or die "Can't execute query: " . $dbh->errstr . "\n";
              if(@row3=$qh_ExonCoord->fetchrow_array()){
                 $ex_start=$row3[0]; $ex_end=$row3[1];
                 print LOG "$chr\t$strand\t$transc_name\t$tx_start\t$tx_end\t$ex_start\t$ex_end\n";
                }
           }
       }
   }
}
close(LOG);
print  LOG "Program complete\n";


