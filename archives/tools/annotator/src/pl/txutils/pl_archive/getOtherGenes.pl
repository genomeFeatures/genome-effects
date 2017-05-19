#!/usr/bin/perl
#############################################################################
#This script assigns genes to transcript that do not ave known Ensembl gene.
#The script detects lines in the inputGenePredictionfile where the gene field is empty. The
#script uses the specified annotation source to check for overlap genes.
#  Usage: ./program -h   --for help
#  usage: ./program -f inputGenePredictionfile -v organismVersion -a genePredictionSource [-o outputFileName]
#
#
# Note: 
# The input file is a gene prediction bed format file
# This script is ran just before parseAndLoad_NovelTranscripts.pl script
#
# Author: Lucie N. Hutchins
#
use DBI;

use vars qw ($opt_h $opt_v $opt_a $opt_l $opt_f $opt_d $opt_o );
use Getopt::Std;getopts('hlv:a:f:d:o:');
if($opt_h||!$opt_v||!$opt_f||!$opt_a) {
    print <<HELP;

 This script assigns genes to transcript that do not ave known Ensembl gene.
 The script detects lines in the inputGenePredictionfile where the gene field is empty.
 The script uses the specified annotation source to check for overlap genes.

 Usage: ./program -h   --for help
 usage: ./program -f inputGenePredictionfile -v organismVersion -a genePredictionSource [-o outputFileName]
 Where inputGenePredictionfile is the name of the transcriptome input file - in genePrediction format.
 If the -o option is not specified, the outpufile name will be: [inputGenePredictionfile]-update.txt

 Note: This script is ran just before parseAndLoad_NovelTranscripts.pl script.

HELP
exit;
}
my $user ='pup';my $pwd  ='puppass';$dbname ='graber_transcriptdb';$host ='harlequin.jax.org';
if($opt_d){$host=$opt_d;}$hostdev="demon.jax.org";
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $qh_getOrgV = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$query="select distinct c.chromosome_id,chromosome_name 
       from chromosome c, transcript t where t.organism_version_id=? and t.chromosome_id=c.chromosome_id";
$getChrom=$dbh->prepare($query);#and g.gene_prediction_id=?   and  g.organism_version_id=t.organism_version_id
$query="select transcript_id from  transcript where organism_version_id=? and chromosome_id =? 
           and strand=? and ((tx_start<=? and ?<= tx_end)||(tx_start>=? and tx_end<=?)||(tx_start<=? and ?<= tx_end))";
$qh_getHits=$dbh->prepare($query);
$getGene=$dbh->prepare("select distinct gene_name from gene_by_annotation where transcript_id=? and gene_prediction_id=?");

$qh_getPredId->execute($opt_a);$qh_getOrgV->execute($opt_v);
($pred_id)=$qh_getPredId->fetchrow_array(); ($orgv_id)=$qh_getOrgV->fetchrow_array();
$out_file=($opt_o)?$opt_o:"$opt_f-update.txt";open(IN,"$opt_f");open(OUT,">$out_file");
if(!OUT || !IN){ print "File eror: $! --\n"; exit;}
#get header line: chrom\tstrand\ttxStart\ttxEnd\texonCount\texonStarts\texonEnds\tName2
$header=<IN>; chomp($header); @fields=split(/\t/,$header);%headerFields=();
for my $i (0 .. $#fields){$fieldname=lc($fields[$i]);$headerFields{"$fieldname"}=$i;}

$getChrom->execute($orgv_id);%chroms=();
while(($chrom_id,$chrom)=$getChrom->fetchrow_array()){$chroms{"$chrom"}=$chrom_id;}
while(<IN>){ chomp($_); split /\t/;
 if(scalar(@_)!=scalar(@fields)){ #missing gene 
    $chrom=$_[$headerFields{"chrom"}];$chrom=~s/chr//i;$strand=$_[$headerFields{"strand"}];
    $chrom_id=$chroms{"$chrom"}; $txStart=$_[$headerFields{"txstart"}];$txEnd=$_[$headerFields{"txend"}];
    $gene="";%genes=(); $qh_getHits->execute($orgv_id,$chrom_id,$strand,$txStart,$txStart,$txStart,$txEnd,$txEnd,$txEnd) or die mysql_error();
    while(($tx_id)= $qh_getHits->fetchrow_array()){ $getGene->execute($tx_id,$pred_id);
        while(($name)=$getGene->fetchrow_array()){$genes{"$name"}=1;}
    }
    if(keys(%genes)>0){
       @genes=keys(%genes); $gene=join(",",@genes);push(@_,$gene);print OUT join("\t",@_)."\n";
    }else{push(@_,$gene);print OUT join("\t",@_)."\n";} 
 }else{print OUT "$_\n";}
}
print "Program complete\n";

