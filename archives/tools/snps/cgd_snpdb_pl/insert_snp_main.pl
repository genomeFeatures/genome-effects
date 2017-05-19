#!/usr/bin/perl

###########################################################################################################################
## insert_snp_main.pl
#  This script load snp_main,snp_main_allele_conflict,
#  snp_error_log,snp_CpGSites
# Input: a configuration file

#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : June 2012
#   

#   Usage: perl insert_snp_main.pl -c <list_of_chromosome> -f <SNP config file>
#   where:
#    -h  displays this help message
#    -c  A file containing a list of chromosome names (one chromosome per line)  
#    -f  SNP config file
#
###########################################################################################################################
use DBI;
use Time::localtime;

use vars qw ($opt_h $opt_c $opt_f);
use Getopt::Std;

getopts('hc:f:');
if($opt_h||!$opt_f) {
    print <<HELP;

 This script load snp_main,snp_main_allele_conflict,snp_error_log,snp_CpGSites.
 Runs after snpValidator.

Usage:

   perl insert_snp_main.pl -c <list_of_chromosome> -f <SNP config file>

Arguments:
   -h  displays this help message
   -c  A file containing a list of chromosome names (one chromosome per line)  
   -f  SNP config file

Example: perl insert_snp_main.pl -c chromosome_list_file.txt -f snp_base_dir/snp_source.config
HELP
exit;
}

my $user ='lnh';
my $pwd  ='lucie98';

$dbname ='cgdsnpdb';
$host ='cgd.jax.org';

my $db = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pwd);
if(!$db){print  "Could not connect to database :$!\n"; exit;} 
################################### snp_main ##########################################
$qh_getsnpid=$db->prepare("select snpid from snp_main where snpid=?");
$qh_getSourceid=$db->prepare("select source_id from snp_source where source_name=?");
my $query="drop temporary table if exists snp_main_temp";
my $qh_drop_snp_temp = $db->prepare($query);
$query="create temporary table snp_main_temp(snpid int unsigned default 0,
          provided_ref_allele char(1),local_ref_allele char(1),
          snp_allele char(1),is_conflict tinyint default 0,miscount tinyint default 0,
          is_intergenic tinyint default 0, mutation_type tinyint default 0,
          is_CpG_site tinyint default 0,ref_strain_id smallint default 0,source_id tinyint default 0,
          has_allele_contflict tinyint default 0,is_found tinyint default 0,
          index(is_found),index(has_allele_contflict),
          index(snpid),index(snpid,source_id))";
$qh_create_snpMain_temp=$db->prepare($query);
$query="load data local infile ? into table snp_main_temp";
my $qh_load_this = $db->prepare($query);
$query="select count(*) as rowcount from snp_main_temp ";
my $qh_snp_rowcount = $db->prepare($query);

###### TAG new SNPs ######################################################################
$query="update snp_main_temp t, snp_main m set t.is_found=1 where t.snpid=m.snpid";
my $qh_update_snpMain_temp=$db->prepare($query);
$query="insert into snp_main select distinct snpid,local_ref_allele,snp_allele,is_conflict,";
$query.="is_intergenic,mutation_type,is_CpG_site 
          from snp_main_temp where is_found=0 and is_CpG_site<127 ";
my $qh_insert_newSNP=$db->prepare($query);

##### update error log ##################################################################
$query="insert ignore into snp_error_log select distinct snpid,source_id,is_conflict ";
$query.=" from snp_main_temp where is_conflict>0 and is_CpG_site<127";
$qh_insert_snptemp_error=$db->prepare($query);
$query="insert ignore into snp_error_log select distinct snpid,source_id,miscount ";
$query.=" from snp_main_temp where miscount>0 and is_CpG_site<127 order by snpid";
$qh_insertMiss=$db->prepare($query);

####### update snp_main_allele_conflict ##############################################
$query="update snp_main_temp t, snp_main m set t.has_allele_contflict=1 ";
$query.=" where t.snpid=m.snpid and t.snp_allele!= m.snp_allele";
$qh_update_snpMainConflict=$db->prepare($query);
$query="insert ignore into snp_main_allele_conflict select distinct snpid,";
$query.=" source_id,snp_allele from snp_main_temp 
           where has_allele_contflict=1 and is_CpG_site<127 order by snpid";
$qh_insertMainconf=$db->prepare($query);
$query="insert ignore into snp_error_log select distinct snpid,source_id,6 ";
$query.=" from snp_main_temp where has_allele_contflict=1 and is_CpG_site<127";
$qh_insertAllConfErr=$db->prepare($query);

$query="update snp_main_temp set is_conflict=1 where is_conflict>0";
$qh_update_isConflict=$db->prepare($query);
$query="update snp_main_temp set is_conflict=1 where miscount>0";
$qh_update_MisConflict=$db->prepare($query);
$query="update snp_main_temp set is_conflict=1 where has_allele_contflict>0";
$qh_update_AlleleConflict=$db->prepare($query);

$query="update snp_main m, snp_main_temp t set m.is_conflict=1 
        where  m.snpid=t.snpid and t.is_conflict>0";
$qh_update_mainisConflict=$db->prepare($query);

$qh_analyze=$db->prepare("analyze table snp_main,snp_main_allele_conflict,snp_error_log,snp_CpGSites");

#load cgpsite
$qh_drop_cpg=$db->prepare("drop temporary table if exists snp_CpGSites_temp");
$query="create temporary table if not exists snp_CpGSites_temp(";
$query.="snpid int unsigned default 0,first_allele varchar(8) not null,second_allele varchar(8) not null,";
$query.=" is_found tinyint default 0,index(snpid))";
$qh_create_cpg=$db->prepare($query);
$qh_load_cpg=$db->prepare("load data local infile ? into table snp_CpGSites_temp");
$qh_getCount_cpg=$db->prepare("select count(*) as rowcount from snp_CpGSites_temp");
$qh_update_cpg=$db->prepare("update snp_CpGSites_temp t, snp_CpGSites c set t.is_found=1 where t.snpid=c.snpid");
$query="insert into snp_CpGSites select distinct snpid,first_allele,second_allele";
$query.=" from snp_CpGSites_temp where is_found=0";
$qh_insertNewcpg=$db->prepare($query);

###############################
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
  if($line=~/SNP_BASE_DIR=(.+)$/){$snp_dir=$1;}
  if($line=~/PIPELINE_DIR=(.+)$/){$output_dir=$1;}
  if($line=~/ALL_IN_ONE=(.+)$/){$all_in_one=$1;}
  if($line=~/SNP_FILE_NAME=(.+)$/){$snp_file_name=$1;}
  if($line=~/GENO=(.+)$/){$geno_sufix=$1;} 
  if($line=~/TLOAD=(.+)$/){$snp_main_sufix=$1;} 
  if($line=~/CPG=(.+)$/){$cpg_sufix=$1;} 
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
open(LOG,">$output_dir/$source-insert_log.log");
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
  print "Processing $chr;\n";
  $snp_main2load="$output_dir/$chr$snp_main_sufix";
  $cpg_file="$output_dir/$chr$cpg_sufix";
  #print "$chr: snp_main :$snp_main2load "if(-f $snp_main2load);
  #print "--$chr: snp_cpg :$cpg_file " if(-f $cpg_file);
  #print "\n"; next;
  ##load snp_main
  %uniq=();
  if(-f $snp_main2load){ @snp_main=();$count=0;$snp_main_file="$output_dir/snp_main_temp.txt";
     open(SNP,"$snp_main2load");#remove this file header line
     $header=<SNP>; 
    while(<SNP>){chomp($_);
      ($ID,$ref,$local,$snp,$is_conflict,$miscount,$is_int,$mut,$CpG,$strain)=split(/\t/,$_);
      $qh_getsnpid->execute($ID);
      next if($qh_getsnpid->rows>0);
      push(@snp_main,"$_");++$count;
      if($count%500000==0){
         open(MAIN,">$snp_main_file");$main_rows=0;$has_allele_conflict=0;$is_found=0;
         if(MAIN){  
            foreach my $line(@snp_main){
              ($ID,$ref,$local,$snp,$is_conflict,$miscount,$is_int,$mut,$CpG,$strain)=split(/\t/,$line);
              $local=~s/\s+//g;$snp=~s/\s+//g;
              next if($uniq{$ID}{$local}{$snp}==1);
              print MAIN "$line\t$source_id\t$has_allele_conflict\t$is_found\n";$main_rows+=1;
              $uniq{$ID}{$local}{$snp}=1;
            }
          }close(MAIN);
         if(-f $snp_main_file){
              $qh_drop_snp_temp->execute();$qh_create_snpMain_temp->execute();
              $qh_load_this->execute("$snp_main_file");
              $qh_snp_rowcount->execute();
              if($qh_snp_rowcount->rows>0){
                 ($load_rows)=$qh_snp_rowcount->fetchrow_array();
                 if($load_rows==$main_rows){             #load was success
                    $qh_update_snpMain_temp->execute();  #tag new SNPs
                    $qh_insert_newSNP->execute() ;
                    $qh_insert_snptemp_error->execute(); #update error_log if error found (snpValidator check)
                    $qh_insertMiss->execute(); #update error_log if missing fields(snpValidator check)
                    $qh_update_snpMainConflict->execute(); #set the has_allele_conflict flag 
                    $qh_insertMainconf->execute(); #update snp_main_allele_conflict
                    $qh_insertAllConfErr->execute(); #update error_log
                    $qh_update_isConflict->execute();$qh_update_MisConflict->execute();
                    $qh_update_AlleleConflict->execute();
                    $qh_update_mainisConflict->execute();
                  }
               }
           }
           print LOG "snp_main:$count processed \n";$main_rows=0;@snp_main=();
       }#end of $count%1000000==0
     } #while(<SNP>)
     close(SNP);
   } #end of (-f $snp_main2load)
    #insert the last
   if(@snp_main>0){
     open(MAIN,">$snp_main_file");$main_rows=0;$has_allele_conflict=0;$is_found=0;
     if(MAIN){
            foreach my $line(@snp_main){
              ($ID,$ref,$local,$snp,$is_conflict,$miscount,$is_int,$mut,$CpG,$strain)=split(/\t/,$line);
              $local=~s/\s+//g;$snp=~s/\s+//g;
              next if($uniq{$ID}{$local}{$snp}==1);
              print MAIN "$line\t$source_id\t$has_allele_conflict\t$is_found\n";$main_rows+=1;
              $uniq{$ID}{$local}{$snp}=1;
            }
      }close(MAIN);
      if(-f $snp_main_file){
              $qh_drop_snp_temp->execute();$qh_create_snpMain_temp->execute();
              $qh_load_this->execute("$snp_main_file");
              $qh_snp_rowcount->execute();
              if($qh_snp_rowcount->rows>0){
                 ($load_rows)=$qh_snp_rowcount->fetchrow_array();
                 if($load_rows==$main_rows){             #load was success
                    $qh_update_snpMain_temp->execute();  #tag new SNPs
                    $qh_insert_newSNP->execute() ;
                    $qh_insert_snptemp_error->execute(); #update error_log if error found (snpValidator check)
                    $qh_insertMiss->execute(); #update error_log if missing fields(snpValidator check)
                    $qh_update_snpMainConflict->execute(); #set the has_allele_conflict flag 
                    $qh_insertMainconf->execute(); #update snp_main_allele_conflict
                    $qh_insertAllConfErr->execute(); #update error_log
                    $qh_update_isConflict->execute();$qh_update_MisConflict->execute();
                    $qh_update_AlleleConflict->execute();
                    $qh_update_mainisConflict->execute();
                  }
               }
      }
       print LOG "snp_main: last $main_rows inserted \n";
   }
   #now insert the CpG site data
  if(-f $cpg_file){
     @snp_cpg=();$count=0;$snp_cpg_file="$output_dir/snp_cpg_temp.txt";
     open(SNP,"$cpg_file");#remove this file header line
     $header=<SNP>;
     while(<SNP>){
      chomp($_);push(@snp_cpg,"$_");++$count;
      if($count%500000==0){ open(MAIN,">$snp_cpg_file");$main_rows=0;
         if(MAIN){foreach my $line(@snp_cpg){print MAIN "$line\n";$main_rows+=1;}}close(MAIN);
         if(-f $snp_cpg_file){
              $qh_drop_cpg->execute();$qh_create_cpg->execute();$qh_load_cpg->execute($snp_cpg_file);
              $qh_getCount_cpg->execute();
              if($qh_getCount_cpg->rows>0){
                 ($load_rows)=$qh_getCount_cpg->fetchrow_array();
                 if($load_rows==$main_rows){             #load was success
                    $qh_update_cpg->execute();  #tag new SNPs
                    $qh_insertNewcpg->execute();
                  }
               }
           }
           print LOG "snp_CpGSites: $count processed \n";$main_rows=0;@snp_cpg=();
       }#end of $count%500000==0
     } #while(<SNP>)
     #load last segment
    if(@snp_cpg>0){
       open(MAIN,">$snp_cpg_file");$main_rows=0;
       if(MAIN){
            foreach my $line(@snp_cpg){
              print MAIN "$line\n";$main_rows+=1;}
        }close(MAIN);
        if(-f $snp_cpg_file){
              $qh_drop_cpg->execute();$qh_create_cpg->execute();$qh_load_cpg->execute($snp_cpg_file);
              $qh_getCount_cpg->execute();
              if($qh_getCount_cpg->rows>0){
                 ($load_rows)=$qh_getCount_cpg->fetchrow_array();
                 if($load_rows==$main_rows){             #load was success
                    $qh_update_cpg->execute();  #tag new SNPs
                    $qh_insertNewcpg->execute();
                  }
               }
        }
      print LOG "snp_CpGSites: last $main_rows lines processed \n";$main_rows=0;@snp_cpg=();
    }
   }
   #now load genotype
   
 } #end of chromosome loop
$qh_analyze->execute();
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Program Ends:  $mon/$mday/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
close(LOG);
system("rm $output_dir/*_temp_file.txt");
print"program completed\n";
 
exit(0);

