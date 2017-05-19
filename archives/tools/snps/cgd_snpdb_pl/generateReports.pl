#!/usr/bin/perl

#####################################################################################
## generateReports.pl
#  This script generates SNP reports
# 
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : June 2012
#   Usage: perl load_annotator.pl  -d <directory> 
#
#   where:
#    -h  displays this help message
#    -d  A directory where to store the reports

#
###################################################################################
use DBI;
use Time::localtime;

use vars qw ($opt_h $opt_d);
use Getopt::Std;
                                                                                                             
getopts('hd:');
if($opt_h||!$opt_d) { #||!$opt_c
    print <<HELP;

 This script generates SNP reports

Usage:

   perl generateReports.pl -d <directory> 

Arguments:
   -h  displays this help message
   -d  directory

Example: perl generateReports.pl -d /scratch/data/snps
HELP
exit;
}
my $user ='lnh';
my $pwd  ='lucie98';

$dbname ='cgd_snpdb';
$host ='cgd-dev.jax.org';

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
$query="select source_name,count(distinct snpid) from snp_source s,snp_by_source sb ";
$query.=" where s.source_id=sb.source_id group by source_name";
$qh_getCountBySourec=$dbh->prepare($query);
$query="select i.snpid,i.strain_id,i.genotype_allele,g.strain_id,g.genotype_allele
         from snp_imputed i, snp_genotype g where i.snpid=g.snpid 
         and i.strain_id=g.strain_id and i.genotype_allele in ('A','T','C','G') 
         AND g.genotype_allele in ('A','T','C','G') AND i.genotype_allele !=g.genotype_allele";
$qh_getImputed_genConflict=$dbh->prepare($query);

$dir=$opt_d;$geno_conflict="$opt_d/imputed_genotype_conflict.txt";
open(OUT,">$geno_conflict");
print OUT "Source_id\tisnpid\tistrain\tigeno\tgstrain\tggeno\n";
$qh_getImputed_genConflict->execute();$count=0;
while(($isnpid,$istrain,$igeno,$gstrain,$ggeno)=$qh_getImputed_genConflict->fetchrow_array()){
   print OUT "21\t$isnpid,$istrain,$igeno,$gstrain,$ggeno\n";++$count;
   print "$count lines processed\n" if($count%100000==0);
}
print "Program Complete\n";
exit(0);






