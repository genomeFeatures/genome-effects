#!/usr/bin/perl

##########################################################################################################
## load_genotype.pl
#  This script load SNP genotype into the database.
#  The database table can be either snp_imputed if the SNP genotype has the associeted confidence score
#  or snp_genotype otherwise
#
#  Tables updated: snp_genotype/snp_imputed, snp_genotype_allele_conflict, snp_error_log
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : June 2012
#   Usage: perl load_imputed.pl  -f <SNP config file> -c <chr_list_file>
#
#   where:
#    -h  displays this help message
#    -c  A file containing a list of chromosome names (one chromosome per line)  
#    -f  SNP config file
#
########################################################################################################
use DBI;
use Time::localtime;

use vars qw ($opt_h $opt_f $opt_c);
use Getopt::Std;

getopts('hc:f:');
if($opt_h||!$opt_f) { #||!$opt_c
    print <<HELP;

 This script loads the genotype data

Usage:

   perl load_genotype.pl -f <SNP config file> -c <chromosome list>

Arguments:
   -h  displays this help message
   -f  SNP config file
   -c  file containing chromosome list

Example: perl load_genotype.pl -c chromosome_list_file.txt -f snp_source.config
HELP
exit;
}
my $user ='lnh';
my $pwd  ='sql4lucie';

$dbname ='cgdsnpdb';
$host ='cgd.jax.org';

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
$qh_getSourceid=$dbh->prepare("select source_id from snp_source where source_name=?");
############## snp_imputed
############## 
my $qh_drop_imputed_temp=$dbh->prepare("drop temporary table if exists snp_imputed_temp");
$query="create temporary table snp_imputed_temp(snpid int unsigned default 0,
            strain_id smallint default 0,ref_strain_id smallint default 0,
            genotype_allele char(1),confidence tinyint default 0,
            source_id smallint default 0,is_found tinyint default 0,has_allele_conflict tinyint default 0,
            index(has_allele_conflict),
            index(snpid),index(strain_id),index(genotype_allele),index(source_id))";
my $qh_create_imputed_temp = $dbh->prepare($query);

$query="load data local infile ? into table snp_imputed_temp";
my $qh_load_imp = $dbh->prepare($query);

   $query="update snp_imputed_temp t, snp_imputed p set t.is_found=1 ";
   $query.=" where t.snpid=p.snpid and t.strain_id=p.strain_id and t.source_id=p.source_id 
             and p.genotype_allele in('A','T','C','G') and t.genotype_allele in('A','T','C','G') 
             and t.genotype_allele=p.genotype_allele";
my $qh_update_imp_found = $dbh->prepare($query);
   $query="update snp_imputed_temp t, snp_genotype p set t.has_allele_conflict=1 ";
   $query.=" where t.snpid=p.snpid and t.strain_id=p.strain_id and t.genotype_allele in('A','T','C','G')";
   $query.=" and p.genotype_allele in('A','T','C','G') and t.genotype_allele!=p.genotype_allele ";
my $qh_update_imp_conf_geno = $dbh->prepare($query);
   $query="update snp_imputed_temp t, snp_imputed p set t.has_allele_conflict=1 ";
   $query.=" where t.snpid=p.snpid and t.strain_id=p.strain_id and t.genotype_allele in('A','T','C','G')";
   $query.=" and p.genotype_allele in('A','T','C','G') and t.genotype_allele!=p.genotype_allele ";
my $qh_update_imp_conf_imp = $dbh->prepare($query);

   $query="insert into snp_imputed(snpid,strain_id,genotype_allele,confidence,source_id) ";
   $query.="select distinct snpid,strain_id,genotype_allele,confidence,source_id
             from snp_imputed_temp where is_found=0 ";
my $qh_insert_imp = $dbh->prepare($query);
$query="insert ignore into snp_genotype_conflict select snpid,source_id,strain_id,";
$query.=" genotype_allele from snp_imputed_temp where has_allele_conflict=1";
my $qh_insert_imp_conf=$dbh->prepare($query);
$query="insert ignore into snp_error_log select snpid,source_id,7 ";
$query.=" from snp_imputed_temp where has_allele_conflict=1";
my $qh_inser_imp_error=$dbh->prepare($query);

my $qh_imp_rowcount = $dbh->prepare("select count(*) as rowcount from snp_imputed_temp");
my $analyze_genotype="analyze table snp_imputed,snp_genotype,snp_genotype_conflict, snp_error_log";
my $qh_analyze_genotype = $dbh->prepare($analyze_genotype);

####################### snp_genotype
my $qh_drop_geno_temp=$dbh->prepare("drop temporary table if exists snp_genotype_temp");
$query="create temporary table snp_genotype_temp(snpid int unsigned default 0,
            source_id smallint default 0,strain_id smallint default 0,
            genotype_allele char(1),ref_strain_id smallint default 0,
            is_found tinyint default 0,has_allele_conflict tinyint default 0,
            index(has_allele_conflict),
            index(snpid),index(strain_id),index(genotype_allele),index(source_id))";
my $qh_create_geno_temp = $dbh->prepare($query);

$query="load data local infile ? into table snp_genotype_temp";
my $qh_load_geno = $dbh->prepare($query);

   $query="update snp_genotype_temp t, snp_genotype p set t.is_found=1 ";
   $query.=" where t.snpid=p.snpid and t.strain_id=p.strain_id  and t.source_id=p.source_id
           and t.genotype_allele in('A','T','C','G') and p.genotype_allele in('A','T','C','G') 
           and t.genotype_allele=p.genotype_allele";
my $qh_update_geno_found = $dbh->prepare($query);
   $query="update snp_genotype_temp t, snp_genotype p set t.has_allele_conflict=1 ";
   $query.=" where t.snpid=p.snpid and t.strain_id=p.strain_id and t.genotype_allele in('A','T','C','G')";
   $query.=" and p.genotype_allele in('A','T','C','G') and t.genotype_allele!=p.genotype_allele ";
my $qh_update_geno_conf_geno = $dbh->prepare($query);
   $query="update snp_genotype_temp t, snp_imputed p set t.has_allele_conflict=1 ";
   $query.=" where t.snpid=p.snpid and t.strain_id=p.strain_id and t.genotype_allele in('A','T','C','G')";
   $query.=" and p.genotype_allele in('A','T','C','G') and t.genotype_allele!=p.genotype_allele ";
my $qh_update_imp_conf_geno = $dbh->prepare($query);

   $query="insert ignore into snp_genotype(snpid,source_id,strain_id,genotype_allele) ";
   $query.="select distinct snpid,source_id,strain_id,genotype_allele
             from snp_genotype_temp where is_found=0 ";
my $qh_insert_geno = $dbh->prepare($query);
$query="insert ignore into snp_genotype_conflict select distinct snpid,source_id,strain_id,";
$query.=" genotype_allele from snp_genotype_temp where has_allele_conflict=1";
my $qh_insert_geno_conf=$dbh->prepare($query);
$query="insert ignore into snp_error_log select snpid,source_id,7 ";
$query.=" from snp_genotype_temp where has_allele_conflict=1";
my $qh_inser_geno_error=$dbh->prepare($query);
my $qh_update_snpMain=$dbh->prepare("update snp_main m, snp_genotype_conflict c set m.is_conflict=1 where m.snpid=c.snpid");
my $qh_geno_rowcount = $dbh->prepare("select count(*) as rowcount from snp_genotype_temp");

#######################################################################################
#
#used for cases the snp_allele was not specified
#

# get command line arguments 
# get command line arguments 
$chr_list_file=$opt_c;$config_file=$opt_f;
#now validate the config file and get global variables
if(!(-f $config_file)){
  print STDERR "Bad configuration file name :$config_file\n";
  exit(1);
}
open(CONF,"$config_file");
if(!CONF){
  print STDERR "Could not open the configuration file $config_file:$!\n";
  exit(1);
}
@config=<CONF>;@variable=grep(/=/,@config);
my ($snp_dir,$source,$output_dir,$source_id);my $source_id=0;
$all_in_one=0;$snp_file_name="";$geno_sufix="";$snp_main_sufix="";
$cpg_sufix="";
foreach my $line(@variable){# set global variables from the config file
  chomp($line); 
  if($line=~/SNP_SOURCE=(.+)$/){$source=$1;}
  if($line=~/PIPELINE_DIR=(.+)$/){$output_dir=$1;}

  if($line=~/ALL_IN_ONE=(.+)$/){$all_in_one=$1;}
  if($line=~/SNP_FILE_NAME=(.+)$/){$snp_file_name=$1;}
  if($line=~/GENO=(.+)$/){$geno_sufix=$1;} 
  #if($line=~/TLOAD=(.+)$/){$snp_main_sufix=$1;} 
  #if($line=~/CPG=(.+)$/){$cpg_sufix=$1;} 
}
if(!(-d $output_dir)){
  print STDERR "One of the following is not a directory:$snp_dir or $output_dir.Check the config file settings\n"; 
  exit(1);
}
if($output_dir =~/(.*)\/$/){ #format the directory name
        $output_dir=$1;
}
@chromosomes=(); 
if($all_in_one==0){
    open(IN,"$chr_list_file")or die "Can't open file:$!\n";
    @chromosomes=<IN>;
    close(IN);
 }
 else{$chromosomes[0]=$snp_file_name;}
open(LOG,">$output_dir/$source-geno-insert_log.log");
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Starting load process :  $mon/$mday/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
$qh_getSourceid->execute($source);
if($qh_getSourceid->rows<=0){ #source not found in our database
   print STDERR "The specified source was not found in our database\n"; 
   exit(1);
}
else{($source_id)=$qh_getSourceid->fetchrow_array();}
while(@chromosomes>0){   # process file
  $chr=shift(@chromosomes);chomp($chr);$snp_main2load=""; $geno_file=""; $cpg_file="";
  next if($chr ne "X");
  $geno_file="$output_dir/$chr$geno_sufix";
  if(-f $geno_file){ @snp_main=();$count=0;$main_rows=0;$snp_main_file="$output_dir/snp_geno_temp.txt";
      open(SNP,"$geno_file");$header=<SNP>;
      while(<SNP>){
          chomp($_);push(@snp_main,"$_");++$count;
          if($count%500000==0){
             open(MAIN,">$snp_main_file");$main_rows=0;$is_found=0;$has_allele_conflict=0;
             if(MAIN){ 
                foreach my $line(@snp_main){
                 print MAIN "$line\t$is_found\t$has_allele_conflict\n";$main_rows+=1;}
              }close(MAIN);
             if(-f $snp_main_file){
                if($source_id==21 || $source_id==16){
                   $qh_drop_imputed_temp->execute();$qh_create_imputed_temp->execute();
                   $qh_load_imp->execute($snp_main_file);$qh_imp_rowcount->execute();
                   if($qh_imp_rowcount->rows>0){
                      ($load_rows)=$qh_imp_rowcount->fetchrow_array();
                      if($load_rows==$main_rows){ #load was success
                        $qh_update_imp_found->execute();$qh_update_imp_conf_imp->execute();
                        $qh_update_imp_conf_geno->execute();$qh_insert_imp->execute();
                        $qh_insert_imp_conf->execute();$qh_inser_imp_error->execute();
                        print LOG "File snp_imputed: $load_rows of $main_rows loaded\n";
                      }
                    }
                 }
                 else{
                   $qh_drop_geno_temp->execute();$qh_create_geno_temp->execute();
                   $qh_load_geno->execute($snp_main_file);$qh_geno_rowcount->execute();
                   if($qh_geno_rowcount->rows>0){
                      ($load_rows)=$qh_geno_rowcount->fetchrow_array();
                      if($load_rows==$main_rows){ #load was success
                        $qh_update_geno_found->execute();$qh_update_geno_conf_geno->execute();
                        $qh_update_imp_conf_geno->execute();$qh_insert_geno->execute();
                        $qh_insert_geno_conf->execute();$qh_inser_geno_error->execute();
                        print LOG "File snp_genotype: $load_rows of $main_rows loaded\n";
                      }
                    }
                }
              }
            print"$count processed \n";@snp_main=();$main_rows=0; 
          } #($count%1000000==0)
      } #end of while(<SNP>)
       #load last segment 
      if(@snp_main>0){
         open(MAIN,">$snp_main_file");$main_rows=0;$is_found=0;$has_allele_conflict=0;
         if(MAIN){ 
                foreach my $line(@snp_main){
                 print MAIN "$line\t$is_found\t$has_allele_conflict\n";$main_rows+=1;}
          }close(MAIN);
          if(-f $snp_main_file){
                if($source_id==21 || $source_id==16){
                   $qh_drop_imputed_temp->execute();$qh_create_imputed_temp->execute();
                   $qh_load_imp->execute($snp_main_file);$qh_imp_rowcount->execute();
                   if($qh_imp_rowcount->rows>0){
                      ($load_rows)=$qh_imp_rowcount->fetchrow_array();
                      if($load_rows==$main_rows){ #load was success
                        $qh_update_imp_found->execute();$qh_update_imp_conf_imp->execute();
                        $qh_update_imp_conf_geno->execute();$qh_insert_imp->execute();
                        $qh_insert_imp_conf->execute();$qh_inser_imp_error->execute();
                        print LOG "File snp_imputed: $load_rows of $main_rows loaded\n";
                      }
                    }
                 }
                 else{
                   $qh_drop_geno_temp->execute();$qh_create_geno_temp->execute();
                   $qh_load_geno->execute($snp_main_file);$qh_geno_rowcount->execute();
                   if($qh_geno_rowcount->rows>0){
                      ($load_rows)=$qh_geno_rowcount->fetchrow_array();
                      if($load_rows==$main_rows){ #load was success
                        $qh_update_geno_found->execute();$qh_update_geno_conf_geno->execute();
                        $qh_update_imp_conf_geno->execute();$qh_insert_geno->execute();
                        $qh_insert_geno_conf->execute();$qh_inser_geno_error->execute();
                        print LOG "File snp_genotype: $load_rows of $main_rows loaded\n";
                      }
                    }
                }
           }
           @snp_main=();$main_rows=0; 
      }
   } #end of if(-f $geno_file)
 } #end of chromosome loop
 $qh_update_snpMain->execute();
 #$qh_analyze_genotype->execute();
close(SNP);
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
close(LOG); 
print"program completed\n";
 
exit(0);

