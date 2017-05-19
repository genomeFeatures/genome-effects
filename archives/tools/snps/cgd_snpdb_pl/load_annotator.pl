#!/usr/bin/perl

##########################################################################################################
## load_annotator.pl
#  This script load SNP annotator generated files into the database.
# 
#
#  Tables updated: snp_transcript, snp_aminoacid, snp_main
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : June 2012
#   Usage: perl load_annotator.pl  -f <SNP config file> -c <chr_list_file>
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

 This script loads the annotator generated files data

Usage:

   perl load_annotator.pl -f <SNP config file> -c <chromosome list>

Arguments:
   -h  displays this help message
   -f  SNP config file
   -c  file containing chromosome list

Example: perl load_annotator.pl -c chromosome_list_file.txt -f snp_source.config
HELP
exit;
}
my $user ='lnh';
my $pwd  ='lucie98';

$dbname ='cgd_snpdb';
$host ='cgd-dev.jax.org';

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
$qh_getSourceid=$dbh->prepare("select source_id from snp_source where source_name=?");
############## snp_transcript_temp
my $qh_drop_tx_temp=$dbh->prepare("drop temporary table if exists snp_transcript_temp");
$query="create temporary table snp_transcript_temp(snpid int unsigned default 0,
            _loc_func_key tinyint default 0,strand tinyint default 0,
            transcript_local_id int unsigned default 0,loc_rank smallint default 0,
            gene_id  int unsigned default 0,is_found tinyint default 0,index(is_found),
            index(snpid),index(snpid,_loc_func_key,transcript_local_id),index(transcript_local_id))";
my $qh_create_tx_temp = $dbh->prepare($query);

$query="load data local infile ? into table snp_transcript_temp";
my $qh_load_tx = $dbh->prepare($query);

   $query="update snp_transcript_temp t, snp_transcript p set t.is_found=1 ";
   $query.=" where t.snpid=p.snpid and t._loc_func_key=p._loc_func_key 
             and t.transcript_local_id=p.transcript_local_id";
my $qh_update_tx_found = $dbh->prepare($query);

   $query="update snp_main t, snp_transcript_temp p set t.is_intergenic=0 ";
   $query.=" where t.snpid=p.snpid";
my $qh_update_is_intergenic = $dbh->prepare($query);

   $query="insert ignore into snp_transcript(snpid,_loc_func_key,strand,transcript_local_id,loc_rank,gene_id) ";
   $query.="select distinct snpid,_loc_func_key,strand,transcript_local_id,loc_rank,gene_id
             from snp_transcript_temp where is_found=0 ";
my $qh_insert_tx = $dbh->prepare($query);

my $qh_tx_rowcount = $dbh->prepare("select count(*) as rowcount from snp_transcript_temp");
my $analyze_tx="analyze table snp_transcript_temp,snp_aminoacid";
my $qh_analyze_tx = $dbh->prepare($analyze_tx);

####################### snp_aminoacid
my $qh_drop_amino_temp=$dbh->prepare("drop temporary table if exists snp_aminoacid_temp");
$query="create temporary table snp_aminoacid_temp(snpid int unsigned default 0,
            transcript_local_id int unsigned default 0,_frame_key  tinyint default 0,
            PosInCDS int default 0,PosInProtein int default 0,
            ref_aa char(3),ref_codon char(3), snp_aa char(3),snp_codon char(3),is_found tinyint default 0,
            index(is_found),
            index(snpid),index(snpid,_frame_key,transcript_local_id),index(transcript_local_id))";
my $qh_create_amino_temp = $dbh->prepare($query);

$query="load data local infile ? into table snp_aminoacid_temp";
my $qh_load_amino = $dbh->prepare($query);
my $qh_amino_rowcount = $dbh->prepare("select count(*) as rowcount from snp_aminoacid_temp");
 $query="update snp_aminoacid_temp t, snp_aminoacid p set t.is_found=1 ";
 $query.=" where t.snpid=p.snpid and t._frame_key=p._frame_key 
             and t.transcript_local_id=p.transcript_local_id";
my $qh_update_amino_found = $dbh->prepare($query);
 $query="insert ignore into snp_aminoacid(snpid,transcript_local_id,_frame_key,PosInCDS,PosInProtein,ref_aa,ref_codon,snp_aa,snp_codon) ";
   $query.="select distinct snpid,transcript_local_id,_frame_key,PosInCDS,PosInProtein,ref_aa,ref_codon,snp_aa,snp_codon
             from snp_aminoacid_temp where is_found=0 ";
my $qh_insert_amino = $dbh->prepare($query);

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
$all_in_one=0;$snp_file_name="";$amino_sufix="";$snp_main_sufix="";
$tx_sufix="";
foreach my $line(@variable){# set global variables from the config file
  chomp($line); 
  if($line=~/SNP_SOURCE=(.+)$/){$source=$1;}
  if($line=~/PIPELINE_DIR=(.+)$/){$output_dir=$1;}
  if($line=~/ALL_IN_ONE=(.+)$/){$all_in_one=$1;}
  if($line=~/SNP_FILE_NAME=(.+)$/){$snp_file_name=$1;}
  if($line=~/AMINO=(.+)$/){$amino_sufix=$1;} 
  if($line=~/TRANSCRIPT=(.+)$/){$tx_sufix=$1;} 
}
if(!(-d $output_dir)){
  print STDERR "One of the following is not a directory:$snp_dir or $output_dir.Check the config file settings\n"; 
  exit(1);
}
if($output_dir =~/(.*)\/$/){$output_dir=$1;}
@chromosomes=(); 
if($all_in_one==0){
    open(IN,"$chr_list_file")or die "Can't open file:$!\n";
    @chromosomes=<IN>;close(IN);
 }
 else{$chromosomes[0]=$snp_file_name;}
open(LOG,">$output_dir/$source-annotator-insert_log.log");
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
  $chr=shift(@chromosomes);chomp($chr);$amino_file=""; $tx_file="";
  print "Processing $chr\n";
  $amino_file="$output_dir/$chr$amino_sufix";$tx_file="$output_dir/$chr$tx_sufix";
  @snp_tx=(); @snp_amino=();#load snp_transcript file
  if(-f $tx_file){  @snp_tx=();$count=0;$main_rows=0;$snp_main_file="$output_dir/snp_tx_temp.txt";
      open(SNP,"$tx_file");
      while(<SNP>){
          chomp($_);push(@snp_tx,"$_");++$count;
          if($count%500000==0){open(MAIN,">$snp_main_file");$main_rows=0;$is_found=0;
             if(MAIN){ 
                foreach my $line(@snp_tx){
                 print MAIN "$line\t$is_found\n";$main_rows+=1;}
              }close(MAIN);
             if(-f $snp_main_file){
                   $qh_drop_tx_temp->execute();$qh_create_tx_temp->execute();
                   $qh_load_tx->execute($snp_main_file);$qh_tx_rowcount->execute();
                   if($qh_tx_rowcount->rows>0){
                      ($load_rows)=$qh_tx_rowcount->fetchrow_array();
                      if($load_rows==$main_rows){ #load was success
                        $qh_update_tx_found->execute();
                        $qh_update_is_intergenic->execute();$qh_insert_tx->execute();
                        print LOG "File $tx_file: $load_rows of $main_rows loaded\n";
                      }
                   }
              }
            @snp_tx=();$main_rows=0; 
          } #($count%1000000==0)
      } #end of while(<SNP>)
       #load last segment  of snp_transcript load
      if(@snp_tx>0){
         open(MAIN,">$snp_main_file");$main_rows=0;$is_found=0;
         if(MAIN){ 
                foreach my $line(@snp_tx){print MAIN "$line\t$is_found\n";$main_rows+=1;}
          }close(MAIN);
          if(-f $snp_main_file){
              $qh_drop_tx_temp->execute();$qh_create_tx_temp->execute();
              $qh_load_tx->execute($snp_main_file);$qh_tx_rowcount->execute();
              if($qh_tx_rowcount->rows>0){
                ($load_rows)=$qh_tx_rowcount->fetchrow_array();
                 if($load_rows==$main_rows){ #load was success
                    $qh_update_tx_found->execute();
                    $qh_update_is_intergenic->execute();$qh_insert_tx->execute();
                    print LOG "File $tx_file: $load_rows of $main_rows loaded\n";
                 }
              }
          }@snp_tx=();$main_rows=0; 
      }#end of @snp_tx>0 for the last segment
    close(SNP);
   } #end of if(-f $tx_file)
   if(-f $amino_file){
      @snp_tx=();$count=0;$main_rows=0;$snp_main_file="$output_dir/snp_amino_temp.txt";
      open(SNP,"$amino_file");
      while(<SNP>){
          chomp($_);push(@snp_tx,"$_");++$count;
          if($count%500000==0){open(MAIN,">$snp_main_file");$main_rows=0;$is_found=0;
             if(MAIN){ 
                foreach my $line(@snp_tx){
                 print MAIN "$line\t$is_found\n";$main_rows+=1;}
              }close(MAIN);
             if(-f $snp_main_file){
                   $qh_drop_amino_temp->execute();$qh_create_amino_temp->execute();
                   $qh_load_amino->execute($snp_main_file);$qh_amino_rowcount->execute();
                   if($qh_amino_rowcount->rows>0){
                      ($load_rows)=$qh_amino_rowcount->fetchrow_array();
                      if($load_rows==$main_rows){ #load was success
                        $qh_update_amino_found->execute();
                        $qh_insert_amino->execute();
                        print LOG "File $amino_file: $load_rows of $main_rows loaded\n";
                      }
                   }
              }
            @snp_tx=();$main_rows=0; 
          } #($count%1000000==0)
      } #end of while(<SNP>)
       #load last segment  of snp_transcript load
      if(@snp_tx>0){
         open(MAIN,">$snp_main_file");$main_rows=0;$is_found=0;
         if(MAIN){ 
                foreach my $line(@snp_tx){print MAIN "$line\t$is_found\n";$main_rows+=1;}
          }close(MAIN);
          if(-f $snp_main_file){
                 $qh_drop_amino_temp->execute();$qh_create_amino_temp->execute();
                 $qh_load_amino->execute($snp_main_file);$qh_amino_rowcount->execute();
                 if($qh_amino_rowcount->rows>0){
                   ($load_rows)=$qh_amino_rowcount->fetchrow_array();
                    if($load_rows==$main_rows){ #load was success
                       $qh_update_amino_found->execute();$qh_insert_amino->execute();
                       print LOG "File $amino_file: $load_rows of $main_rows loaded\n";
                     }
                  }
           }
          @snp_tx=();$main_rows=0; 
      }
     close(SNP);
   }#end of if(-f $amino_file)
 } #end of chromosome loop
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
close(LOG); 
$qh_analyze_tx->execute();
print"program completed\n";
 
exit(0);

