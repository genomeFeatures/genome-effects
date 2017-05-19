#!/usr/bin/perl

##########################################################################################################
## load_genotype_conflict.pl
#
#  Tables updated: snp_genotype_allele_conflict, snp_error_log
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : June 2012
#
########################################################################################################
use DBI;
use Time::localtime;

use vars qw ($opt_h $opt_f);
use Getopt::Std;

getopts('hf:');
if($opt_h||!$opt_f) { #||!$opt_c
    print <<HELP;

 This script loads the genotype conflict data

Usage:

   perl load_genotype_conflict.pl -f <SNP conflict file>

Arguments:
   -h  displays this help message
   -f  SNP conflict file

Example: perl load_genotype_conflict.pl -f /scratch/data/snp/geno_conflict_file_name
HELP
exit;
}
my $user ='lnh';
my $pwd  ='lucie98';

$dbname ='cgd_snpdb';
$host ='cgd-dev.jax.org';

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
############## snp_imputed
############## 
####################### snp_genotype
my $qh_drop_geno_temp=$dbh->prepare("drop temporary table if exists snp_genotype_temp");
$query="create temporary table snp_genotype_temp(
          snpid int unsigned default 0,
          source_id smallint default 0,strain_id smallint default 0,
          genotype_allele char(1),is_found tinyint default 0,
          index(snpid),index(strain_id),index(genotype_allele),index(source_id))";
my $qh_create_geno_temp = $dbh->prepare($query);

$query="load data local infile ? into table snp_genotype_temp";
my $qh_load_geno = $dbh->prepare($query);

   $query="update snp_genotype_temp t, snp_genotype_conflict p set t.is_found=1 ";
   $query.=" where t.snpid=p.snpid and t.strain_id=p.strain_id  and t.source_id=p.source_id
           and t.genotype_allele in('A','T','C','G') and p.genotype_allele in('A','T','C','G') 
           and t.genotype_allele=p.genotype_allele";
my $qh_update_geno_found = $dbh->prepare($query);

$query="insert into snp_genotype_conflict select distinct snpid,source_id,strain_id,";
$query.=" genotype_allele from snp_genotype_temp where is_found=0";
my $qh_insert_geno_conf=$dbh->prepare($query);
$query="insert ignore into snp_error_log select snpid,source_id,7 ";
$query.=" from snp_genotype_temp";
my $qh_inser_geno_error=$dbh->prepare($query);
my $qh_geno_rowcount = $dbh->prepare("select count(*) as rowcount from snp_genotype_temp");

#######################################################################################
#
#used for cases the snp_allele was not specified
#

# get command line arguments 
# get command line arguments 
$snp_file=$opt_f;
#now validate the config file and get global variables
if(!(-f $snp_file)){
  print STDERR "Bad SNP file name :$snp_file\n";
  exit(1);
}
open(CONF,"$snp_file");
if(!CONF){
  print STDERR "Could not open the configuration file $config_file:$!\n";
  exit(1);
}
open(LOG,">snps/geno_conflict-insert_log.log");
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Starting load process :  $mon/$mday/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
open(SNP,"$snp_file");$header=<SNP>;
@snp_main=();
$snp_main_file="/scratch/data/snps/geno_temp.txt";
while(<SNP>){
    chomp($_);push(@snp_main,"$_");++$count;
    if($count%500000==0){
       open(MAIN,">$snp_main_file");$main_rows=0;$is_found=0;
       if(MAIN){ 
          foreach my $line(@snp_main){
             ($ids,$istrain,$igeno,$ostrain,$ogeno)=split(",",$line);
             ($source_id,$snpid)=split("\t",$ids);
            print MAIN "$snpid\t$source_id\t$istrain\t$igeno\t$is_found\n";$main_rows+=1;}
       }close(MAIN);
       #print "Progam complete\n"; exit(0);
       if(-f $snp_main_file){
         $qh_drop_geno_temp->execute();$qh_create_geno_temp->execute();
         $qh_load_geno->execute($snp_main_file);$qh_geno_rowcount->execute();
         if($qh_geno_rowcount->rows>0){
           ($load_rows)=$qh_geno_rowcount->fetchrow_array();
           if($load_rows==$main_rows){ #load was success
              $qh_update_geno_found->execute();
              $qh_insert_geno_conf->execute();$qh_inser_geno_error->execute();
              print LOG "File snp_genotype: $load_rows of $main_rows loaded\n";
           }
         }
      }
      print"$count processed \n";@snp_main=();$main_rows=0; 
    } #($count%1000000==0)
} #end of while(<SNP>)
 #load last segment 
 if(@snp_main>0){
   $qh_drop_geno_temp->execute();$qh_create_geno_temp->execute();
   $qh_load_geno->execute($snp_main_file);$qh_geno_rowcount->execute();
   if($qh_geno_rowcount->rows>0){
     ($load_rows)=$qh_geno_rowcount->fetchrow_array();
     if($load_rows==$main_rows){ #load was success
        $qh_update_geno_found->execute();
        $qh_insert_geno_conf->execute();$qh_inser_geno_error->execute();
        print LOG "File snp_genotype: $load_rows of $main_rows loaded\n";
      }
   }
 }
close(SNP);
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
close(LOG); 
print"program completed\n";
 
exit(0);

