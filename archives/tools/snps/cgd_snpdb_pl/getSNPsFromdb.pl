#!/usr/bin/perl

######################################################################
## This script generates SNPs file by chromosome
## from the snp database
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : April 2011
#   Usage: ./getSNPsFromdb.pl 
#
#  Output:
#      string submitID;
#      string chr;
#      long pos;
#      string b6_allele
#      string snp_allele;
#      string strand;
#####################################################################
## set default variables
$data_dir=`pwd`;                    #default data directory
$local_db="cgd_snpdb37";               #default cgd snp database name
$localhost="lookingglass.jax.org";           #deafult db server
$luser="pup";$lpwd="puppass";  
   
use DBI;
use vars qw ( $opt_h $opt_o $opt_d $opt_l $opt_f);
use Getopt::Std;

getopts('ho:d:l:f:');
if($opt_h) {
    print <<HELP;

 This script generates SNP files to be used by
 the SNP annotator program. Data is generated from the specified
 SNP database

Usage:

   perl getSNPsFromdb.pl [-o <data_path>][-d <db_name>][-l <db_server>][-f <chrlist>]

Arguments:
   -h  displays this help message
   -o  path to data directory (optional)
   -f  file name containing list of chromosomes
   -d  cgd SNP database name (cgd_snpdb37, or cgdsnpdb) (optional)
   -l  cgd snp database server (cgd.jax.org default) (optional)

HELP
exit;
}

if($opt_o){$data_dir=$opt_o;} 
# check if data directory is provided by the user
$data_dir=~s/\s+$//;$data_dir=~s/\/$//;
$chr_log= "$data_dir/cgd_chr_snplist.txt";   #file to store the list of snps per chromosome   

if($opt_d){$local_db=$opt_d;}
if($opt_l){$localhost=$opt_l;}
if($opt_f){$chr_log=$opt_f;}

my $dbh = DBI->connect("DBI:mysql:$local_db:$localhost",$luser,$lpwd);
if(!$dbh){
   print  "Could not connect to the database :$!\n"; 
   exit(1);
}
$get_chr="select chromosome_id,chromosome_name from snp_chromosome where chromosome_id<=22";
$get_snps=" select snpid,bp_position,black6_allele,snp_allele ";
$get_snps.=" from snp_main where is_public=1 and chromosome_id=?";

my $qh_get_chr= $dbh->prepare($get_chr)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_get_snps= $dbh->prepare($get_snps)or die "Couldn't prepare statement: " . $dbh->errstr;


$qh_get_chr->execute(); ## get the list of chromosome
open(LOG,">$chr_log");
while(($chromosome_id,$chromosome_name)=$qh_get_chr->fetchrow_array()){
    print LOG "$data_dir/chr$chromosome_name-snps.txt\n";  
    next if(-f "$data_dir/chr$chromosome_name-snps.txt");
    open(ANT,">$data_dir/chr$chromosome_name-snps.txt");
    if(ANT){
       print ANT "snpid\tChromosome\tPosition\treference_base\tconsensus_base\tstrand\n";
       $qh_get_snps->execute($chromosome_id);
       while(($snpid,$bp_position,$black6_allele,$snp_allele)=$qh_get_snps->fetchrow_array()){
               $alleles="$black6_allele/$snp_allele";
               print ANT "$snpid\t$chromosome_name\t$bp_position\t$black6_allele\t$snp_allele\t+\n";
       }
    } 
   close(ANT);
   print " Chromosome $chromosome_name data generated\n";   
 }
print "Program complete\n";
close(LOG);
exit(0);



