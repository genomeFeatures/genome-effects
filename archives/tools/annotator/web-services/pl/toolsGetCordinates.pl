#!/usr/bin/perl

###################################################################################################
## toolsGetCoordinates.pl -- 
# Fetches coordinates of specified features from specified database
# This will be ran on the server when ever the annotator client request
# non existing or expired features annotations
#
#
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   The Jackson Laboratory
#   Implimentation date : September 2013
###################################################################################################
use DBI;
use POSIX;
use XML::Simple;

# create object
$xml = XML::Simple->new (ForceArray => 1);

use vars qw ($opt_h $opt_o $opt_v $opt_a $opt_d $opt_f $opt_p);
use Getopt::Std;

getopts('ho:v:a:f:d:p:');
if($opt_h||(!$opt_v)) {
    print <<HELP;

This script generates the coordinates  of transcripts for the specified organism version and annotation source.
The features are stored into the database using the zero-based cordinates standard for feature starts
and 1-based cordinate standard for feature ends

## Features coordinates in genePrediction file are generated in 0-based <start> standard and 1-base <end>
    
Example: perl toolsGetCordinates.pl -o dirName -v mm10 -a GRCm38-ensGene 
The above will use current Ensembl annotations for mm10 to generate all the transcripts
in gene prediction format and store the file (mm10-GRCm38-ensGene-transcripts.txt) under dirName/

Arguments:
  -h  displays this help message
  -v  is the ucsc organism version (default mm9)
  -a  optional, gene prediction source to use
  -o  optional, is the output directory, default current working directory
  -f  optional, output file name (relative to the specified output directory)
 
HELP
exit;
}
$gene_prediction=($opt_a)?$opt_a:"refGene";
#set defaults
#Note:
#I need to instantiate 3 db connections: ucsc, ensembl, and graber_lab for annotation downloads
#
# add two more options: -p ensmbl or -p ucsc (to specify the db host)
#
$ucsc_host="genome-mysql.cse.ucsc.edu"; $ucsc_user="genome";
my $connection="mysql -h$ucsc_host -u$ucsc_user -A";
my $graber_user ='pup';my $graber_pwd  ='puppass';$graber_dbname ='graber_transcriptdb';$graber_host ='harlequin.jax.org';
if($opt_d){$host=$opt_d;}$hostdev="demon.jax.org";
$webservice="http://$hostdev/transcriptdb/web-services/index.php";
my $dbh = DBI->connect("DBI:mysql:database=$graber_dbname;host=$graber_host",$graber_user, $graber_pwd);
my $ucsc_dbh= DBI->connect("DBI:mysql:database=$opt_v;host=$ucsc_host",$ucsc_user,"");
#if(!$ucsc_dbh){print  "Could not connect to database :$!\n".$ucsc_dbh->errstr; exit;} 
##################################################
# get hits for the search term
##################################################
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$qh_getPredName=$dbh->prepare("select gene_prediction_name from gene_prediction where gene_prediction_id=?");

$query="select distinct transcript_id,transcript_name,gene_prediction_id from transcript_exon
       where organism_version_id=? and gene_prediction_id =?";
$qh_getTxByAnnot=$dbh->prepare($query);

my $getTxcord="select chromosome_name,strand,tx_start,tx_end from transcript e, chromosome c";
   $getTxcord.=" where transcript_id=? and organism_version_id=? and e.chromosome_id=c.chromosome_id";
   $qh_getTxcord = $dbh->prepare($getTxcord)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTxExons="select distinct chromosome_name,exon_start,exon_end,strand,exon_frame 
                from transcript_exon t,exon e,chromosome c ";
   $getTxExons.=" where transcript_id=? and transcript_name=? and t.organism_version_id=? and gene_prediction_id in (?)";
   $getTxExons.=" and t.organism_version_id=e.organism_version_id and
                  t.exon_id=e.exon_id and e.chromosome_id=c.chromosome_id";
   $getTxExons.=" order by exon_start asc";
my $qh_getTxExons = $dbh->prepare($getTxExons)or die "Couldn't prepare statement: ".$dbh->errstr;

my $qh_getCCgene=$dbh->prepare("select distinct gene_prediction_id, gene_id, gene_name from cc_founders_genes where gene_prediction_id in (?)");

$query="select cdsStart,cdsEnd from transcript_translation ";
$query.=" where transcript_id in (?) and transcript_name in(?) and organism_version_id=? and gene_prediction_id in (?) order by cdsStart asc";
my $qh_getCDS=$dbh->prepare($query);

my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $qh_getOrgV = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTxGene="select distinct gene_name from gene_by_annotation g, transcript_by_annotation t
                where g.organism_version_id=? and t.transcript_name=? and t.gene_prediction_id in (?)
                and t.transcript_id=g.transcript_id and t.gene_prediction_id=t.gene_prediction_id ";
my $qh_getTxGene = $dbh->prepare($getTxGene)or die "Couldn't prepare statement: ".$dbh->errstr;

$query="select distinct sample_name from rnaseq_sample s, rnaseq_transcript t ";
$query.=" where t.transcript_name=? and gene_prediction_id=? and organism_version_id=? ";
$query.=" and t.sample_id=s.sample_id";
$qh_getSampleName=$dbh->prepare($query); %ccgenemap=(); %txDup=();
sub getTranscript{
 ($prediction_ids,$tx_qh,$org_vid,$chrom_prefix)=@_;
 $count=0;
 while(($tx_id,$tx,$prediction_id)=$tx_qh->fetchrow_array()){
       if($prediction_ids ne ""){
            if($prediction_ids=~/^\s*(\d+)\s*$/){next if($1!=$prediction_id);}
            next if(!($prediction_ids=~/$prediction_id/));
        }
        if($prediction_id>24&&$prediction_id<=40){
           if(!exists($ccgenemap{"$prediction_id"})){$qh_getCCgene->execute($prediction_id);
               while(($gene_prediction_id, $gene_id, $gene_name )=$qh_getCCgene->fetchrow_array()){ 
                      $ccgenemap{"$gene_prediction_id"}{"$gene_id"}=$gene_name;}
            }
        }
   $qh_getTxcord->execute($tx_id,$org_vid);#get the coordinates of this transcript
   ($chrom_name,$txstrand,$tx_start,$tx_end)=$qh_getTxcord->fetchrow_array();
   $gene="";$exonStarts="";$exonEnds="";$cdsStarts="0";$cdsEnds="0";$annot="";
   $qh_getTxExons->execute($tx_id,$tx,$org_vid,$prediction_id);$qh_getPredName->execute($prediction_id);
   ($annot)=$qh_getPredName->fetchrow_array();$sample_id="";
   if($tx=~/^rsd_/){ $qh_getSampleName->execute($tx,$prediction_id,$org_vid);
       while(($sample)=$qh_getSampleName->fetchrow_array()){
          if($sample_id eq ""){$sample_id=$sample;}else{$sample_id.=",$sample";}}
   }$annot=($sample_id eq "")?$annot:$sample_id;$rows=""; $exoncount=$qh_getTxExons->rows; 
   while(($chromosome_name,$exstart,$exend,$strand)=$qh_getTxExons->fetchrow_array()){
        if($exonStarts eq ""){$exonStarts=$exstart;$exonEnds=$exend;
        }else{$exonStarts.=",$exstart";$exonEnds.=",$exend";}
   }$cdsStarts=0;$cdsEnds=0;
   $qh_getCDS->execute($tx_id,$tx,$org_vid,$prediction_id);($cdsStart,$cdsEnd)=$qh_getCDS->fetchrow_array();
   #$cdsStarts+=1; #adjust cdsStarts to 1-base
   $cdsStart=($cdsStart>0)?$cdsStart:$cdsStarts;$cdsEnd=($cdsEnd>0)?$cdsEnd:$cdsEnds;
   $qh_getTxGene->execute($org_vid,$tx,$prediction_id);$gene_name="";$gene_id="";
   $row="";
   while(($gene_name)=$qh_getTxGene->fetchrow_array()){
          $gene_id=$gene_name if(($gene_name=~/^ENS/)||($gene_name=~/^rsd_g/));
   }$gene_id=($gene_id eq "")?$gene_name:$gene_id;
    if($prediction_id>24&&$prediction_id<=40){$gene_name=$ccgenemap{"$prediction_id"}{"$gene_id"} if($gene_name eq "");}
   $gene_name=$gene_id if($gene_name eq "");
   ++$count; print "$count processed\n" if($count%10000==0);
   $chrom_name=($chrom_prefix ne "")?"$chrom_prefix$chrom_name":$chrom_name;
   $txDup{"$chrom_name\t$txstrand\t$tx_start\t$tx_end"}{"$exoncount\t$exonStarts\t$exonEnds"}.="$tx=$gene_name=$annot=$cdsStart=$cdsEnd;";
 }

}
############################
$organism_version="";$organism_version=$opt_v if($opt_v);$org_vid=0;$organism_version=~s/\s+//g;
$file_name="";$header="name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tname2\tsource\n";
if($opt_d){
  $file_name=($opt_f)?$opt_f:"$host"."_"."$organism_version-$gene_prediction-transcripts.txt";
  if($opt_o){$opt_o=~s/\/\s*$//;$file_name="$opt_o/$file_name";}
  $query="select name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,name2 from $gene_prediction";
  if($gene_prediction=~/miRNA/){
     $query="select name,chrom,strand,chromStart,chromEnd from $gene_prediction";
     $header="name\tchrom\tstrand\tchromStart\tchromEnd\tsource\n";
  }elsif(($gene_prediction=~/microsat/)||($gene_prediction=~/simpleRepeat/)){
     $query="select name,chrom,chromStart,chromEnd from $gene_prediction";
     $header="name\tchrom\tchromStart\tchromEnd\tsource\n";
  }elsif($gene_prediction=~/\s*snp\d+\s*/){
     $query="select name,chrom,strand,chromStart,chromEnd,class from $gene_prediction";
     $header="name\tchrom\tstrand\tchromStart\tchromEnd\tclass\tsource\n";
  }
  if($ucsc_dbh){$ucsc_dbh->{PrintError} = 0; # disable
     $qh_ucscgetTxByAnnot=$ucsc_dbh->prepare($query);my $rv = $qh_ucscgetTxByAnnot->execute(); 
     if($rv){open(TS,">$file_name");print TS "$header";
       if($gene_prediction=~/miRNA/){
          while(($name,$chrom,$strand,$txStart,$txEnd)=$qh_ucscgetTxByAnnot->fetchrow_array()){
                print TS "$name\t$chrom\t$strand\t$txStart\t$txEnd\t$host-$opt_v-$gene_prediction\n";}
       }elsif(($gene_prediction=~/microsat/)||($gene_prediction=~/simpleRepeat/)){
          while(($name,$chrom,$txStart,$txEnd)=$qh_ucscgetTxByAnnot->fetchrow_array()){
                print TS "$name\t$chrom\t$txStart\t$txEnd\t$host-$opt_v-$gene_prediction\n";}
       }elsif($gene_prediction=~/\s*snp\d+\s*/){
          while(($name,$chrom,$strand,$txStart,$txEnd,$class)=$qh_ucscgetTxByAnnot->fetchrow_array()){
                print TS "$name\t$chrom\t$strand\t$txStart\t$txEnd\t$class\t$host-$opt_v-$gene_prediction\n";}
       }else{
           while(($name,$chrom,$strand,$txStart,$txEnd,$cdsStart,$cdsEnd,$exonCount,
                  $exonStarts,$exonEnds,$name2)=$qh_ucscgetTxByAnnot->fetchrow_array()){
             print TS "$name\t$chrom\t$strand\t$txStart\t$txEnd\t$cdsStart\t$cdsEnd\t";$exonStarts=~s/,\s*$//;$exonEnds=~s/,\s*$//;
           print TS "$exonCount\t$exonStarts\t$exonEnds\t$name2\t$host-$opt_v-$gene_prediction\n";
           }
       }
 }}
}
else{
$qh_getOrgV->execute("$organism_version");if($qh_getOrgV->rows>0){($org_vid)=$qh_getOrgV->fetchrow_array();}
$prediction_ids="";$qh_getPredId->execute($gene_prediction);
if($qh_getPredId->rows>0){($prediction_ids)=$qh_getPredId->fetchrow_array();}
 #set the result path
$file_name="";$type=~s/\/\s*$//g;$feature_type="Transcript"; $chrom_prefix=($opt_p)?$opt_p:"";
if($opt_f){$file_name=$opt_f;}
else{$file_name="$organism_version-$gene_prediction-transcripts.txt";}
if($opt_o){$opt_o=~s/\/\s*$//;$file_name="$opt_o/$file_name";}$count=0; $qh_getTxByAnnot->execute($org_vid,$prediction_ids);
open(TS,">$file_name"); 
print TS "name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tname2\tsource\n";
if($qh_getTxByAnnot->rows>0){$count=0;
    getTranscript($prediction_ids,$qh_getTxByAnnot,$org_vid,$chrom_prefix);
    while(($key,$value)=each(%txDup)){ #this gets rid of duplicate rows
         while(($cord,$feature)=each(%{$value})){ 
             @features=split(";",$feature);%tranc=();%genes=();%annots=();%cdsStarts=();%cdsEnds=();
             foreach $item (@features){
                  ($tx,$gene,$annot,$cdsStart,$cdsEnd)=split("=",$item);
                   $tranc{"$tx"}=1;$genes{"$gene"}=1 if($gene ne "");$annots{"$annot"}=1;
                  if($cdsStart>0){$cdsStarts{"$cdsStart"}=1;$cdsEnds{"$cdsEnd"}=1;}
             }$tx=join(",",keys %tranc);$gene=join(",",keys %genes);$annot=join(",",keys %annots);
             $cdsStart=join(",",keys %cdsStarts);$cdsEnd=join(",",keys %cdsEnds);
             if(keys(%cdsStarts)==0){$cdsStart=0;$cdsEnd=0;}
             print TS "$tx\t$key\t$cdsStart\t$cdsEnd\t$cord\t$gene\t$annot\n";
         }
      }
 }
}
close(TS);
 $zipfile="$file_name"."."."gz"; $command="gzip -c $file_name > $zipfile"; system($command); 
 if(-f $file_name){system("rm $file_name");}
exit(0);

