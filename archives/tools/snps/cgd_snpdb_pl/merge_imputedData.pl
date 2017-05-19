#!/usr/bin/perl -w

########################## merge_imputedData.pl ########################################################################
# 
# 
# this script will combine the imputed snp data with it corresponding confidence score data
#
# Usage: ./merge_imputedData.pl  snp_file_dir imputed_genotype_dir confidence_score_dir strains_list.txt 
#
# Notes: The script expect: 
#     1. the genotype file to be commas separated and have the following fields:
#         ...,chrom_strand, list of strains_genotype
#     2. the confidence file to be commas separated and have the following fields:
#        ,chrom_strand, list of strains_confidence_scores
#
# also the two input files should not have header lines, and the snp_file.txt should be sorted by snp_id
# 
use DBI;
use Time::localtime;

use vars qw ($opt_h $opt_c $opt_f);
use Getopt::Std;

getopts('hc:f:');
if($opt_h||(!$opt_f || !$opt_c)){
    print <<HELP;

 This script will combine the imputed snp data with it corresponding confidence score data

Usage:

   perl merge_imputedData.pl -c <list_of_chromosome> -f <SNP config file>

Arguments:
   -h  displays this help message
   -c  A file containing a list of chromosome names (one chromosome per line)  
   -f  SNP config file

Example: perl merge_imputedData.pl  -f /scratch/data/snps/imputed/Szatkiewicz_et_al._2008/snp_source.config
HELP
exit;
}
my $user ='lnh';my $pass  ='lucie98';
$dbname ='cgd_snpdb';
$host ='cgd-dev.jax.org';

$pwd=`pwd`;
chomp($pwd);print "$pwd\n";
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
$get_source_id="select source_id from snp_source where source_name=?";
my $qh_getSourceid=$dbh->prepare($get_source_id);
my $get_strains="select strain_id from snp_strain_synonym where synonym_name=?";
my $qh_getStrainId=$dbh->prepare($get_strains);

my $get_strain="select strain_id from snp_strain where strain_name=?";
my $qh_getStrain=$dbh->prepare($get_strain);
my $get_max_strainid="select max(strain_id) from snp_strain";
my $qh_getMaxStrainid=$dbh->prepare($get_max_strainid);
my $get_strain_source="select* from snp_strain_by_source where strain_id=? and source_id=?";
my $qh_get_strain_source=$dbh->prepare($get_strain_source);
my $insertStrain="insert ignore into snp_strain values(?,?)";
my $qh_insertStrain=$dbh->prepare($insertStrain);
my $insertStrainSyn="insert ignore into snp_strain_synonym values(?,?)";
my $qh_insertStrainSyn=$dbh->prepare($insertStrainSyn);
my $insertStrainBySource="insert ignore into snp_strain_by_source values(?,?)";
my $qh_insertStrainBySource=$dbh->prepare($insertStrainBySource);

$qh_drop_geno=$dbh->prepare("drop temporary table if exists geno_load_temp");
$qh_drop_geno_temp=$dbh->prepare("drop temporary table if exists geno_temp");
$qh_drop_geno_conf=$dbh->prepare("drop temporary table if exists geno_conf_load_temp");

 $query="create temporary table if not exists geno_load_temp";
 $query.=" (row_line text, index(row_line(12)))";
my $qh_create_geno_load_temp=$dbh->prepare($query);
$query="create temporary table if not exists geno_temp";
 $query.=" (id int primary key auto_increment,row_line text, index(row_line(12)))";
my $qh_create_geno_temp=$dbh->prepare($query);

$query="create temporary table if not exists geno_conf_load_temp";
 $query.=" (accession varchar(250),row_line text, index(accession))";
my $qh_create_geno_conf_temp=$dbh->prepare($query);

my $qh_loadGenotype=$dbh->prepare("load data local infile ? into table geno_load_temp ignore 1 lines");
my $qh_loadGenotype_temp=$dbh->prepare("insert into geno_temp(row_line) select row_line from geno_load_temp");
my $qh_loadGenotype_conf=$dbh->prepare("load data local infile ? into table geno_conf_load_temp");
$qh_getidlist=$dbh->prepare("select id from geno_temp");
$qh_getgeno_count=$dbh->prepare("select max(id) from geno_temp");
$qh_getgeno_conf_count=$dbh->prepare("select count(*) from geno_conf_load_temp");
$qh_getgeno_line=$dbh->prepare("select row_line from geno_temp where id=?");
$qh_getgeno_conf_line=$dbh->prepare("select row_line from geno_conf_load_temp where accession=?");

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
my ($snp_dir,$source,$file_type,$output_dir,$source_id,$hmmConfDir,$hmmFiledDir);my $count =0;
$start_strain_index=-1;$ref_allele_index=-1;$snp_allele_index=-1;$end_strain_index=-1;
$snpzipfile_sufix="";$genozipfile_sufix="";$geno_sufix="";
$snpid_index=-1;$bp_pos_index=-1;$chrom_index=-1;$strand_index=-1;$file_prefix="";$file_sufix="";
$rsid_index=-1;$all_in_one=0;$snp_file_name="";$confidence_score=-1;$skip_one=0;$source_abrev="";
foreach my $line(@variable){# set global variables from the config file
  chomp($line); 
  if($line=~/SNP_SOURCE=(.+)$/){$source=$1;}
  if($line=~/SOURCE_ABREV=(.+)$/){$source_abrev=$1;}
  if($line=~/SNP_BASE_DIR=(.+)$/){$snp_dir=$1;}
  if($line=~/PIPELINE_DIR=(.+)$/){$output_dir=$1;}
  if($line=~/HMM_CONF_DIR=(.+)$/){$hmmConfDir=$1;}
  if($line=~/HMM_GENO_DIR=(.+)$/){$hmmFiledDir=$1;}


  if($line=~/SNPID_INDEX=(.+)$/){$snpid_index=$1;}
  if($line=~/SNP_CHROM_INDEX=(.+)$/){$chrom_index=$1;}
  if($line=~/SNP_BP_POSITION_INDEX=(.+)$/){$bp_pos_index=$1;}
  if($line=~/SNP_STRAND_INDEX=(.+)$/){$strand_index=$1;} 
  if($line=~/SNP_REF_ALLELE_INDEX=(.+)$/){$ref_allele_index=$1;}
  if($line=~/SNP_OTHER_ALLELE_INDEX=(.+)$/){$snp_allele_index=$1;}
  if($line=~/SNP_FIRST_STRAIN_INDEX=(.+)$/){$start_strain_index=$1;}
  if($line=~/SNP_LAST_STRAIN_INDEX=(.+)$/){$end_strain_index=$1;}
  if($line=~/SNP_RSID_INDEX=(.+)$/){$rsid_index=$1;}

  if($line=~/CONFIDENCE_SCORE=(.+)$/){$confidence_score=$1;}
  if($line=~/ALL_IN_ONE=(.+)$/){$all_in_one=$1;}
  if($line=~/SNP_FILE_NAME=(.+)$/){$snp_file_name=$1;}
  if($line=~/STRAIN_INDEX_SKIP_ONE=(.+)$/){$skip_one=$1;}


  if($line=~/SNP_FILE_TYPE=(.+)$/){$file_type=$1;}
  if($line=~/SNP_FILE_PREFIX=(.+)$/){$file_prefix=$1;}
  if($line=~/SNP_FILE_SUFIX=(.+)$/){$file_sufix=$1;}
  if($line=~/GENO_FILE_SUFIX=(.+)$/){$geno_sufix=$1;}
  if($line=~/SNP_ZIPFILE_SUFIX=(.+)$/){$snpzipfile_sufix=$1;}
  if($line=~/GENO_ZIPFILE_SUFIX=(.+)$/){$genozipfile_sufix=$1;}
  if($line=~/STRAIN_INDEX_SKIP_ONE=(.+)$/){$skip_one=$1;}  
  
}
if(($start_strain_index==-1)||($bp_pos_index==-1)||($chrom_index==-1)){
  print STDERR "One of the following indice has a negative value in the config file:start_strain_index,bp_pos_index,chrom_index\n"; 
  exit(1);
}
if(!(-d $snp_dir)|| !(-d $output_dir)){
  print STDERR "One of the following is not a directory:$snp_dir or $output_dir\nCheck the config file settings\n"; 
  exit(1);
}
if($snp_dir =~/(.*)\/$/){ 
        $snp_dir=$1; #format the directory name
}
if($output_dir =~/(.*)\/$/){ #format the directory name
        $output_dir=$1;
}
if($hmmConfDir =~/(.*)\/$/){ 
        $hmmConfDir=$1; #format the directory name
}
if($hmmFiledDir=~/(.*)\/$/){ #format the directory name
        $hmmFiledDir=$1;
}
if(!(-d $hmmConfDir)){
  print STDERR "input hmmconfidence file directory required\n";
}
if(!(-d $hmmFiledDir)){
  print STDERR "hmmfiled file directory required\n";
  exit(1);
}
@chromosomes=(); 
if($all_in_one==0){
    open(IN,"$chr_list_file")or die "Can't open file:$!\n";
    @chromosomes=<IN>;
 }
else{$chromosomes[0]=$snp_file_name;}
open(LOG,">$snp_dir/merge_log.log");
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
#print "source id is : $source_id\n";
while(@chromosomes>0){   # process file
   $chr=shift(@chromosomes);chomp($chr);$chr=~s/chr//;
   $chr="M" if(($chr =~/M/i)&&($source_id==20));
   $conf_zip="$file_prefix".$chr."$snpzipfile_sufix"; $geno_zip="$file_prefix".$chr."$genozipfile_sufix";
   $hmmconffile_zip=$hmmConfDir."/$conf_zip";$hmmfiledfile_zip=$hmmFiledDir."/$geno_zip";
   $hmmconffile=$hmmConfDir."/$file_prefix".$chr."$file_sufix";
   $hmmfiledfile=$hmmFiledDir."/$file_prefix".$chr."$geno_sufix";
   print "Processing $hmmconffile and $hmmfiledfile\n";
   print LOG "Processing $hmmconffile and $hmmfiledfile\n";
   if(!(-f $hmmfiledfile)){
      if(!(-f $hmmfiledfile_zip)){
           print STDERR "The genotype zip file($hmmfiledfile) required\n";next;
      }
      chdir($hmmFiledDir);system("unzip $geno_zip");
      if(!(-f $hmmfiledfile)){print STDERR "The genotype file($hmmfiledfile) required\n";next;}
      chdir($pwd);
   }
   if(!(-f $hmmconffile)){
      if(!(-f $hmmconffile_zip)){
           print STDERR "The confidence zip file($hmmconffile) required\n";next;
      }
      chdir($hmmConfDir);system("unzip $conf_zip");
      if(!(-f $hmmconffile)){print STDERR "The genotype file($hmmconffile) required\n";next;}
      chdir($pwd);
   }
   open(INHMCON,"$hmmconffile");open(INHMFIL,"$hmmfiledfile");$head=<INHMFIL>;$line=<INHMCON>;
  chomp($line);
   my(@linecontent,@header);$delim=",";if($file_type!=2){$delim="\t";}
   @linecontent= split("$delim",$line);@header= split("$delim",$head);
   if(@linecontent!=@header){
       print LOG "===\ngenotype header:\n$line\nFields count : ".@header."\nConfidence Score header\n$head\n";
       print LOG "Fields count:".@linecontent."\n";
      next;}
   open(TP,">$hmmConfDir/temp_load.txt");
   while(<INHMCON>){chomp($_);
     @fields=split("$delim",$_);$fields[0]=~s/\s*//g;print TP "$fields[0]\t$_\n";
    }
   close(TP); close(INHMCON); close(INHMFIL);
   if($end_strain_index==-1){$end_strain_index= @linecontent;--$end_strain_index;}
   if($end_strain_index<$start_strain_index){
        print "Bad strain index $end_strain_index<$start_strain_index\n";next;
    }my %strains=();$i=$start_strain_index;
   open(OUT,">$hmmfiledfile-merged.txt");
   print OUT "$header[0]$delim$header[1]$delim$header[2]$delim";
   print OUT "$header[3]$delim$header[4]$delim$header[5]";
   $i=$start_strain_index; $geno_alleles="";$score="confidence";
   while($i<=$end_strain_index){
      $strain_name=$header[$i];$strain_name=~s/^\s+//;$strain_name=~s/\s+$//;$strain_name=~s/,$//;
      $geno_alleles.="$delim$strain_name$delim$score";++$i; 
   }
   print OUT "$geno_alleles\n";
  #now check if these strains already exists into our database
  ########################################
  while($i<=$end_strain_index){
      $strain_name=$linecontent[$i];$strain_id=0;$strain_name=~s/^\s+//;$strain_name=~s/\s+$//;
      $strain_name=~s/,$//;$qh_getStrain->execute($strain_name);#check first in snp_strain table
      if($qh_getStrain->rows>0){($strain_id)=$qh_getStrain->fetchrow_array();
         $qh_getStrainId->execute($strain_name); #check in synonym table
         if($qh_getStrainId->rows<=0){$qh_insertStrainSyn->execute($strain_id,$strain_name);}
      }
      else{ #new strain or synonymous
         $qh_getStrainId->execute($strain_name); #check first in snp_synonym
         if($qh_getStrainId->rows>0){($strain_id)=$qh_getStrainId->fetchrow_array();} #strain exists
         else{ #insert new strain
            $qh_getMaxStrainid->execute();
            if($qh_getMaxStrainid->rows>0){
               ($strain_id)=$qh_getMaxStrainid->fetchrow_array();++$strain_id;
             }print "$strain_name\t$strain_id needs to be inserted\n";
            $qh_insertStrain->execute($strain_id,$strain_name);$qh_getStrain->execute($strain_name);
            if($qh_getStrain->rows>0){
               ($strain_id)=$qh_getStrain->fetchrow_array();
               $qh_insertStrainSyn->execute($strain_id,$strain_name);
               print "$strain_name\t$strain_id inserted\n";
            }
         }
      }
      $qh_get_strain_source->execute($source_id,$strain_id);
      if($qh_get_strain_source->rows<=0){ ##insert this strain for this source if does not already exists 
         $qh_insertStrainBySource->execute($source_id,$strain_id);
       }
      $strains{$i}=$strain_id; #a dictionary mapping strain index to strain id   
      if($skip_one==0){++$i; }else{$i+=2;}
  }
  $geno_rows=0;$conf_rows=0;$geno_line_count=0; $conf_line_count=0;
  #
  #get the number of lines in geno and confidence score files
  #
  $linecont = `wc -l $hmmfiledfile`;chomp($linecont);if($linecont=~/^(\d+)\s$hmmfiledfile/){$geno_line_count=$1;}
  $linecont = `wc -l $hmmconffile`;chomp($linecont);if($linecont=~/^(\d+)\s$hmmconffile/){$conf_line_count=$1;}
  $qh_drop_geno->execute();$qh_drop_geno_conf->execute();$qh_drop_geno_temp->execute();
  --$geno_line_count;--$conf_line_count;
  $qh_create_geno_temp->execute() or die "geno temp not created:".mysql_error();
  $qh_create_geno_conf_temp->execute()or die "conf temp not created:".mysql_error();
  $qh_create_geno_load_temp->execute();
  $qh_loadGenotype->execute("$hmmfiledfile")or die "geno temp not loaded:".mysql_error();
  $qh_loadGenotype_temp->execute(); #generate ids
  $qh_loadGenotype_conf->execute("$hmmConfDir/temp_load.txt")or die "conf temp not loaded:".mysql_error();
  #get the number of rows loaded
  $qh_getgeno_count->execute();$qh_getgeno_conf_count->execute();
  if($qh_getgeno_count->rows>0){($geno_rows)=$qh_getgeno_count->fetchrow_array();}
  if($qh_getgeno_conf_count->rows>0){($conf_rows)=$qh_getgeno_conf_count->fetchrow_array();}
  if($conf_line_count!=$conf_rows){
      print LOG "\nBad confidence score load: only $conf_rows of $conf_line_count loaded\n";next;
   }
  if($geno_line_count!=$geno_rows){
      print LOG "\nBad genotype load: only $geno_rows of $geno_line_count loaded\n";next;
   }
  $qh_getidlist->execute() or die mysql_error();
  if($qh_getidlist->rows>0){
     while(($id)=$qh_getidlist->fetchrow_array()){ 
         $qh_getgeno_line->execute($id)or die mysql_error();
         if($qh_getgeno_line->rows>0){
           ($geno_line)=$qh_getgeno_line->fetchrow_array();
            chomp($geno_line);@genocontent= split("$delim",$geno_line);
            my $geno_id=$genocontent[0];$geno_id=~s/\s*//g;
             $qh_getgeno_conf_line->execute("$geno_id") or die "Id not found:".mysql_error();
            while(($line)=$qh_getgeno_conf_line->fetchrow_array()){
               @confcontent=split("$delim",$line);$i=$start_strain_index; $geno_alleles="";
               #if($line=~/^\s*$genocontent[0]\s*$delim/){
                if(@genocontent!=@confcontent){
                     print "===\ngenotype line:\n$geno_line\nconfidence score line\n$line\n";
                     print LOG "===\ngenotype line:\n$geno_line\nconfidence score line\n$line\n";next;
                 }
                  print OUT "$genocontent[0]$delim$genocontent[1]$delim$genocontent[2]$delim";
                  print OUT "$genocontent[3]$delim$genocontent[4]$delim$genocontent[5]";
                  while($i<=$end_strain_index){
                        $allele=$genocontent[$i];$allele=~s/^\s+//;$allele=~s/\s+$//;
                        if($confcontent[$i]>=0.9){$score=2;}elsif($confcontent[$i]<=0.6){$score=0;}else{$score=1;}
                        $geno_alleles.="$delim$allele$delim$score";++$i; 
                   }print OUT "$geno_alleles\n";
              #}
            }
         }
     }
   }
 ########################################
  print "Chromosome $chr done\n";close(OUT);
 }
print "program complete\n";
close(IN);
exit(0);
