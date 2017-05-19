#!/usr/bin/perl

###################################################################################################
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implementation date : May 2013
#
###################################################################################################
use DBI;
use POSIX;

use vars qw ($opt_h $opt_s $opt_o $opt_v $opt_l $opt_d );
use Getopt::Std;

getopts('hlso:v:d:');
if($opt_h||(!$opt_v&&!$opt_s&&!$opt_l)) {
    print <<HELP;

This script generates unique genes along the genome for a given organism version.

Usage: perl generateUniqGenes.pl -v <organism version> [-o outputDir]
Arguments:
   -h  displays this help message
   -o  is the output directory
   -v  is the ucsc organism version (default mm9)
   -l  display the list of current gene predictions for the specified organism version
   -d  specifies the database server (default harlequin.jax.org)

Example: perl generateUniqGenes.pl -v mm9 -l
The above will display all the gene predictions for mm9

Example: perl generateUniqGenes.pl -o features -v mm9 
The above will generate uniq gene regions for mm9

Assumptions:
1.Similar gene regions across different versions of the same organism share different gene rsd numbers
  because these rsd numbers are organism version specific
2.foreach version of the same organism
   a. get transcript-gene mapping 
HELP
exit;
}

#set defaults
my $user ='pup';my $pwd  ='puppass';$dbname ='graber_transcriptdb';$host ='harlequin.jax.org';
 $organism_version="mm9";
if($opt_d){$host=$opt_d;}
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
#Steps:
# 1. Get the list of chromosome of the organism version
# 2. for each chrom
#    a. get all the provided transcript definitions on both strands sorted by exon_start then,
#       loop through every strand to filter overlapping transcript(aggregate view) and generate unique genomic regions
#       with the associated transcriptList and geneList,
#    b. for each uniq genomic region,
#       b.1 assign a gene rsd number (check first from the database)
#       b.2 generate records for the new two tables
#    c. load data into the database
#     
######################################################
$query="select distinct chromosome_id,strand,tx_start,tx_end,t.transcript_id 
         from transcript t, transcript_by_annotation ta";
$query.=" where t.organism_version_id=?  and gene_prediction_id not in(3,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)
          and t.organism_version_id=ta.organism_version_id ";
$query.=" and t.transcript_id =ta.transcript_id order by  chromosome_id,strand,t.tx_start";
$qh_getChromTx= $dbh->prepare($query);


$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");

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

my $getTxGene="select distinct gene_name from gene_by_annotation where organism_version_id=? and transcript_id =? ";
my $qh_getTxGene = $dbh->prepare($getTxGene)or die "Couldn't prepare statement: ".$dbh->errstr;
$qh_getChrName=$dbh->prepare("select chromosome_name from chromosome where chromosome_id=?");

$query="select novel_gene_name from novel_genes where organism_version_id=? ";
$qh_getTxRsds=$dbh->prepare($query);
$qh_getTxList=$dbh->prepare("select distinct transcript_id from novel_genes where novel_gene_name=? ");
$qh_getTxCord=$dbh->prepare("select chromosome_id, strand,tx_start,tx_end from transcript where transcript_id =?");
$qh_getTxAcc=$dbh->prepare("select distinct transcript_name from transcript_by_annotation where transcript_id=?");
$qh_getGeneAcc=$dbh->prepare("select distinct gene_name from gene_by_annotation where transcript_id=?");
#get the current max of rsd_gnumber
$query="select MAX(cast(substring(novel_gene_name from locate('rsd_g',novel_gene_name)+5) as unsigned)) from novel_genes ";
$qh_getCurrentMaxRsdNumber= $dbh->prepare($query);

##### return the start and end cordinates of a gene #####################
sub getGeneCordinates{
 ($start_ref,$end_ref,$chrom_id_ref,$strand_ref,$gene)=@_; 
 $qh_getTxList->execute($gene);
 while(($txid)=$qh_getTxList->fetchrow_array()){ $qh_getTxCord->execute($txid);
    if($qh_getTxCord->rows>0){
       ($$chrom_id_ref,$$strand_ref,$start,$end)=$qh_getTxCord->fetchrow_array();}
    $$start_ref=(($$start_ref==0)||($$start_ref>$start))?$start:$$start_ref;
    $$end_ref=(($$end_ref==0)||($$end_ref<$end))?$end:$$end_ref;
  }
}

#######################################################################
# loadCurrentRsd: index current novel genes into a data structure
# in memory to facilitate tx-rsd mapping
sub loadCurrentRsd{
  my($rsd_hash_ref,$organism_version_id)=@_; $qh_getTxRsds->execute($organism_version_id);
  if($qh_getTxRsds->rows>0){
    while(($gene)=$qh_getTxRsds->fetchrow_array()){
       $gstart=0;$gend=0;$gene="";$chrom=0;$strand="";
       getGeneCordinates(\$gstart,\$gend,\$chrom,\$strand,$gene);
       $key="$chrom-$strand-$gstart-$gend";
       ${$rsd_hash_ref}{"$organism_version_id"}{"$key"}=$gene;
    }
  }
}
sub getTxAccList{
 ($tx_list_ref,$txid_list)=@_; split ",",$txid_list;
 foreach $txid(@_){ $qh_getTxAcc->execute($txid);
    while(($transcript)=$qh_getTxAcc->fetchrow_array()){${$tx_list_ref}{"$transcript"}=1;}
  }
}
sub getGeneAccList{
 ($gene_list_ref,$txid_list)=@_; split ",",$txid_list;
 foreach $txid(@_){ $qh_getGeneAcc->execute($txid);
     while(($gen)=$qh_getGeneAcc->fetchrow_array()){${$gene_list_ref}{"$gen"}=1;}
   }
}
#######################################################################
%chrMap=();
$qh_getCurrentMaxRsdNumber->execute();my $currentMaxRsdNumber=0;
if($qh_getCurrentMaxRsdNumber->rows){$currentMaxRsdNumber=$qh_getCurrentMaxRsdNumber->fetchrow_array();}
$currentMaxRsdNumber=$currentMaxRsdNumber<=0?0:$currentMaxRsdNumber;
%novelMap=();$local_counter=$currentMaxRsdNumber; %currentGenes=();
$gene_accession_prefix="rsd_g"; $organism_version=$opt_v if($opt_v);$org_vid=0;$organism_version=~s/\s+//g;
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
 ### get the list of all the transcripts of this organism version 
 $qh_getChromTx->execute($org_vid);
 if($qh_getChromTx->rows<=0){
  print "*********\nNote: nothing found for  $organism_version \n";
  print " Please run the script with the --help for usage \n";
  exit(0);
 }
 print "$organism_version gene prediction sources\n";
 $opt_o=~s/\/\s*$//;$opt_t=~s/\/\s*$//g; $file_name="$opt_o/$organism_version-uniqRegions.txt";}
 open (TS,">$file_name") or die "$!\n"; print TS "chrom\tgeneStart\tgeneEnd\ttranscriptsList\tgenesList\tlocal_geneID\n";
 $current_start=0;$current_end=0;$ex_start=0;$ex_end=0;$count=0;%region_map=();%rsdmap=();$tx_list="";%txmap=();
 #generate uniq genomic regions
 while(($chromosome_id,$strand,$tx_start,$tx_end,$tx_id)=$qh_getChromTx->fetchrow_array()){
        $tx_start+=1;
         if($current_start==0){$current_start=$tx_start;$current_end=$tx_end;$tx_list="$tx_id";$txmap{$tx_id}=1;}
         else{
             if(($tx_start>=$current_start)&&($tx_start<=$current_end)){
               if($current_end<$tx_end){$current_end=$tx_end;} #adjust the current end 
               if(!exists($txmap{$tx_id})){$txmap{$tx_id}=1;$tx_list.=",$tx_id";}
             }else{ #this is a new genomic region, there is no overlap with previous transcript region
                $region_map{$chromosome_id}{$strand}{$current_start}{"cord"}=$current_end;
                $region_map{$chromosome_id}{$strand}{$current_start}{"tx"}=$tx_list;++$count;
                $current_start=$tx_start;$current_end=$tx_end;$txmap{$tx_id}=1;$tx_list="$tx_id";
             }
         }   
  }$region_map{$chromosome_id}{$strand}{$current_start}{"cord"}=$current_end;
  $region_map{$chromosome_id}{$strand}{$current_start}{"tx"}=$tx_list;
  #check if you need to load the rsd for this organism
  if($currentMaxRsdNumber>0){if(!exists($rsdmap{"$orgv_id"})){
      loadCurrentRsd(\%rsdmap,$org_vid);}}
  #now loop through the uniq gene list to get the unique exons
 print "Processing the region map\n";
 for my $chromosome(sort keys(%region_map)){$chrom="";
     if(exists($chrMap{"$chromosome"})){$chrom=$chrMap{"$chromosome"};}
     else{$qh_getChrName->execute($chromosome);($chrom)=$qh_getChrName->fetchrow_array();$chrMap{"$chromosome"}=$chrom;}
      print "Processing $chrom\n";
     for my $strand (sort keys (%{$region_map{$chromosome}})){
         for my $genestart(sort keys(%{$region_map{$chromosome}{$strand}})){
              $geneend=$region_map{$chromosome}{$strand}{$genestart}{"cord"};
              $tx_idlist=$region_map{$chromosome}{$strand}{$genestart}{"tx"};
              %txAcList=();getTxAccList(\%txAcList,$tx_idlist);$transcriptList= join ",",sort keys(%txAcList);
              %geneList=();getGeneAccList(\%geneList,$tx_idlist);
              $geneList=join ",",sort keys( %geneList);$gene_rsdid="";#identical gene share the same accession id
              $key="$chromosome-$strand-$genestart-$geneend";$accession="";
              if(!exists($rsdMap{"$orgv_id"}{"$key"})){++$local_counter;
                 $accession="$gene_accession_prefix$local_counter"; $rsdMap{"$orgv_id"}{"$key"}="$accession";
              }else{$accession=$rsdMap{"$orgv_id"}{"$key"};}
              print TS "$chrom\t$genestart\t$geneend\t$strand\t$transcriptList\t$geneList\t$accession\n";  
          }
      }
   }
 print "$file_name generated\n";
#}#end of query by -t option
print"program completed\n";
exit(0);

