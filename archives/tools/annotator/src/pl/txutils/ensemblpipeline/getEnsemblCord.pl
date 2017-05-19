#!/usr/bin/perl

use DBI;
use vars qw($opt_h,$opt_v,$opt_a,$opt_f);
use Getopt::Std;
getopts('hv:a:f:');

$host="harlequin"; $usr="pup"; $pwd="puppass";
$dbname="graber_transcriptdb";
$dbh=DBI->connect("DBI:mysql:database=$dbname;host=$host",$usr,$pwd);
##############################
$query="select chromosome_name,strand,min(tx_start) as geneStart,max(tx_end) as geneEnd from transcript t,gene_by_annotation g,";
$query.=" chromosome s where g.gene_name=? and g.organism_version_id=? and g.gene_prediction_id=? ";
$query.=" and g.organism_version_id=t.organism_version_id and g.transcript_id=t.transcript_id and t.chromosome_id=s.chromosome_id";
$qh_getGeneCord=$dbh->prepare($query);
$query="select distinct gene_name,mgi_id from gene_by_annotation g,mgi_genes m,";
$query.="(select distinct transcript_id from gene_by_annotation where gene_name=? and organism_version_id=? and gene_prediction_id=?)as t ";
$query.=" where g.transcript_id=t.transcript_id and g.gene_prediction_id=? and g.gene_name=m.mgi_symbol";
$qh_getmgiGene=$dbh->prepare($query);
$qh_getGeneDesc=$dbh->prepare("select biotype,description from ensembl_genes where gene_name=? and organism_version_id=?");
$qh_getOrgv=$dbh->prepare("select organism_version_id from organism_version where ucsc_db=?");
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name=?");
if($opt_h||(!$opt_v && !$opt_a&&!$opt_f)){
   print "Usage: ./program -v organism_version -a annotation_source -f geneList_file\n";
  exit;}
chomp($opt_f);open(IN,"$opt_f");
if(IN){
 $output="$opt_a-juliet.txt"; open(OUT,">$output");$orgv_id=0;$pred_id=0;
 $qh_getOrgv->execute($opt_v);($orgv_id)= $qh_getOrgv->fetchrow_array();
 $qh_getPredId->execute($opt_a);($pred_id)=$qh_getPredId->fetchrow_array();
 print "$output file generated $orgv_id:$pred_id\n";
 %predmap=("mm9-ensGene"=>"mm9-mgiGene","GRCm38-ensGene"=>"GRCm38-mgiGene");
 $qh_getPredId->execute($predmap{"$opt_a"});($mgi_predid)=$qh_getPredId->fetchrow_array();
 if($pred_id>0 && $orgv_id>0){
    print OUT "mgi_id\tmgi_gene\tensembl_geneid\tchromosome_name\tstrand\tgeneStart\tgeneEnd\ttype\tdescription\n";
    while(<IN>){ chomp($_); $gene=$_; $gene=~s/^\s*//; $gene=~s/\s*$//;
       $qh_getGeneCord->execute($gene,$orgv_id,$pred_id) or die "bad query".mysql_error();
       ($chromosome_name,$strand,$geneStart,$geneEnd)=$qh_getGeneCord->fetchrow_array();
       $qh_getmgiGene->execute($gene,$orgv_id,$pred_id,$mgi_predid);($mgi_gene,$mgi_id)=$qh_getmgiGene->fetchrow_array();
       $qh_getGeneDesc->execute($gene,$orgv_id);($biotype,$description)=$qh_getGeneDesc->fetchrow_array();
       $geneStart=($geneStart>0)?$geneStart+1:$geneStart;
       print OUT "$mgi_id\t$mgi_gene\t$gene\t$chromosome_name\t$strand\t$geneStart\t$geneEnd\t$biotype\t$description\n";
    }
 }
}
print "Program complete\n";

exit(0);



