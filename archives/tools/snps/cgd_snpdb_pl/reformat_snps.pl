#!/usr/bin/perl

###########################################################################################################################
## reformat_snps.pl
#  This script reformats the original input snps files to 
#  add the snpid local identifier to the SNPs files to be used in the pipeline
#  the program also generates chromosome ids for X (20), Y(21), and Mt(22) ,and UN(23) 
#    a. SNPs Main file contains: snp_localid,b6_allele/snp_allele,error_flag,chromosome,bp_position,strand
#    b. strains.txt file contains:strain_name,strain_id
#    c. source_missing_snp_allel.txt, a file containing SNPs where all the strains have same genotype
#
# Input: a tab-delimitted/commas separated SNPs file with the following fields:
#     1. LOCAL_IDENTIFIER        
#     2. rsNum ?  
#     3. CHROMOSOME      
#     4. position_currentBuild  
#     5. position_previousBuild?    
#     6. Source (only imputed snps, broad,and gnf have this field)?
#     7. strain_list ...
#     8. strain_list  
#     9. strain_list ...
#     n. strain 

#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : April 2012
#   Modification date : April 2012

#   Usage: perl reformat_snps.pl -c <list_of_chromosome> -f <SNP config file>
#
#   Where: 
#   -h  displays this help message
#   -c  A file containing a list of chromosome names (one chromosome per line)  
#   -f  SNP config file
#
###########################################################################################################################
use DBI;
use POSIX;

use vars qw ($opt_h $opt_c $opt_f);
use Getopt::Std;

getopts('hc:f:');
if($opt_h||!$opt_f) { #||!$opt_c
    print <<HELP;

 This script reformats the original SNPs files to add the snp local id
 and tags all snps where the other allele has more than one genotype
Usage:

   perl reformat_snps.pl -c <list_of_chromosome> -f <SNP config file>

Arguments:
   -h  displays this help message
   -c  A file containing a list of chromosome names (one chromosome per line)  
   -f  SNP config file

Example: perl reformat_snps.pl -c chromosome_list_file.txt -f snp_source.config
HELP
exit;
}
my $user ='lnh';
my $pwd  ='lucie98';

$dbname ='cgdsnpdb';
$host ='cgd.jax.org';

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 

$get_source_id="select source_id from snp_source where source_name=?";
$qh_getSourceid=$dbh->prepare($get_source_id);
$get_strains="select strain_id from snp_strain_synonym where synonym_name=?";
$qh_getStrainId=$dbh->prepare($get_strains);
$get_strain="select strain_id from snp_strain where strain_name=?";
$qh_getStrain=$dbh->prepare($get_strain);
$get_max_strainid="select max(strain_id) from snp_strain";
$qh_getMaxStrainid=$dbh->prepare($get_max_strainid);
$get_strain_source="select* from snp_strain_by_source where strain_id=? and source_id=?";
$qh_get_strain_source=$dbh->prepare($get_strain_source);
my $insertStrain="insert ignore into snp_strain values(?,?)";
$qh_insertStrain=$dbh->prepare($insertStrain);
my $insertStrainSyn="insert ignore into snp_strain_synonym values(?,?)";
$qh_insertStrainSyn=$dbh->prepare($insertStrainSyn);
my $insertStrainBySource="insert ignore into snp_strain_by_source values(?,?)";
$qh_insertStrainBySource=$dbh->prepare($insertStrainBySource);
my $get_snpid="select snpid from snp_position where chromosome_id=? and bp_position=?";
my $qh_get_snpid = $dbh->prepare($get_snpid);


########################
my $create_snp_mult="create table if not exists snp_multiAlleles(
           snpid int unsigned default 0,snp_allele varchar(10),source_id smallint default 0,
          foreign key(snpid) references snp_position(snpid),
          foreign key(source_id) references snp_source(source_id),
          index(snpid),index(source_id))";
my $qh_create_snp_mult = $dbh->prepare($create_snp_mult)or die "Couldn't prepare statement: " . $dbh->errstr;


my $load_multiAlleles="load data local infile ? into table snp_multiAlleles ignore 1 lines";
my $qh_load_multiAlleles = $dbh->prepare($load_multiAlleles)or die "Couldn't prepare statement: " . $dbh->errstr;

my $delete_multiAlleles="delete from snp_multiAlleles where source_id=?";
my $qh_delete_multiAlleles = $dbh->prepare($delete_multiAlleles)or die "Couldn't prepare statement: " . $dbh->errstr;

$qh_drop_error_temp=$dbh->prepare("drop temporary table if exists snp_error");
my $create_error="create temporary table if not exists snp_error(";
   $create_error.="snpid int unsigned default 0,source_id smallint default 0,error_flag tinyint default 0,";
   $create_error.="found tinyint default 0,index(snpid),index(source_id),index(error_flag))";
my $qh_create_error_temp=$dbh->prepare($create_error);
my $update="update snp_error s,snp_error_log sl set s.found=1 where s.snpid=sl.snpid and ";
   $update.="s.source_id=sl.source_id and s.error_flag=sl.error_flag";
$qh_update_error_temp=$dbh->prepare($update);

my $insert_error_temp="insert into snp_error(snpid,source_id,error_flag) ";
   $insert_error_temp.="select distinct snpid,source_id,16 from snp_multiAlleles";
my $qh_insert_error_temp=$dbh->prepare($insert_error_temp);


my $insert_error="insert into snp_error_log(snpid,source_id,error_flag) ";
   $insert_error.="select distinct snpid,source_id,16 from snp_error where found=0";
my $qh_insert_error=$dbh->prepare($insert_error);

my $get_multiAlleles_rowcount="select count(*) as rowcount from snp_multiAlleles where source_id=? ";
my $qh_multiAlleles_rowcount = $dbh->prepare($get_multiAlleles_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;
#######################
my $analyze="analyze table snp_multiAlleles,snp_error_log";
my $qh_analyze = $dbh->prepare($analyze)or die "Couldn't prepare statement: " . $dbh->errstr;

#
#used for cases the snp_allele was not specified
# get_snp_allele($start_strain_index,$ref_allele_index,$line,$file_type)
sub get_snp_allele {
 my ($start_strain_index,$b6_index,$line,$file_type)= @_;
 my $snp_allele="N";@linecontent=();
 if($file_type==2){@linecontent= split(",",$line); # commas separated file
  }else{@linecontent= split(/\t/,$line);}
 if(@linecontent>0){
      my $b6_allele=$linecontent[$b6_index];$b6_allele =~ s/\s*//g;
      $i=$start_strain_index;
      while($i<@linecontent){
          if($linecontent[$i]=~/$b6_allele/i){++$i;}
          else{
            if(($linecontent[$i]=~/A|G|T|C/i)){$snp_allele=$linecontent[$i];last;}  
             ++$i;
          }
      }
   }
 return $snp_allele;
}
sub compl_allele{
 my($allele)= @_;chomp($allele);
 if($allele=~/t/i){$allele="A";}
 elsif($allele=~/a/i){$allele="T";}
 elsif($allele=~/g/i){$allele="C";}
 elsif($allele=~/c/i){$allele="G";}
 return $allele;
}
sub get_b6_index {
  my ($line,$file_type)= @_;my $b6_index=0;my $i=0;my  @linecontent;
  if($file_type==2){@linecontent= split(",",$line); # commas separated file
  }
  else{@linecontent= split(/\t/,$line); # default file type
  }
  if(@linecontent>0){
      while($i<@linecontent){
         if($linecontent[$i]=~/C57BL\/6J/i){ $b6_index=$i;last;}
         ++$i;
      }
   }
  return $b6_index;
}
sub get_str_index {
 my ($line,$file_type,$strain)= @_; my $b6_index=0;my $i=0;my  @linecontent;
 if($file_type==2){@linecontent= split(",",$line); # commas separated file
  }
  else{@linecontent= split(/\t/,$line); # default file type
  }
  $strain=~s/^\s+//;$strain=~s/\s+$//;
  if(@linecontent>0){
      while($i<@linecontent){
         if($linecontent[$i]=~/^\s*$strain\s*$/i){$b6_index=$i;last;}
         ++$i;
      }
   }
  return $b6_index;
}# get command line arguments 
$chr_list_file=$opt_c;$config_file=$opt_f;
#now validate the config file and get global variables
if(!(-f $config_file)){
  print STDERR "Bad configuration file name :$config_file\n";exit(1);
}
open(CONF,"$config_file");
if(!CONF){
  print STDERR "Could not open the configuration file $config_file:$!\n";exit(1);
}
@config=<CONF>;@variable=grep(/=/,@config);
my ($snp_dir,$source,$file_type,$output_dir,$source_id);my $count =0;
$start_strain_index=-1;$ref_allele_index=-1;$snp_allele_index=-1;$end_strain_index=-1;
$snpid_index=-1;$bp_pos_index=-1;$chrom_index=-1;$strand_index=-1;$file_prefix="";$file_sufix="";
$rsid_index=-1;$all_in_one=0;$snp_file_name="";$confidence_score=-1;$skip_one=0;$source_abrev="";
$alleles_index=-1;
#CONFIDENCE_SCORE,ALL_IN_ONE,SNP_FILE_NAME,STRAIN_INDEX_SKIP_ONE
#SOURCE_ABREV=SangerUNC

foreach my $line(@variable){# set global variables from the config file
  chomp($line); 
  if($line=~/SNP_SOURCE=(.+)$/){$source=$1;}
  if($line=~/SOURCE_ABREV=(.+)$/){$source_abrev=$1;}
  if($line=~/SNP_BASE_DIR=(.+)$/){$snp_dir=$1;}
  if($line=~/PIPELINE_DIR=(.+)$/){$output_dir=$1;}

  if($line=~/SNPID_INDEX=(.+)$/){$snpid_index=$1;}
  if($line=~/SNP_CHROM_INDEX=(.+)$/){$chrom_index=$1;}
  if($line=~/SNP_BP_POSITION_INDEX=(.+)$/){$bp_pos_index=$1;}
  if($line=~/SNP_STRAND_INDEX=(.+)$/){$strand_index=$1;} 
  if($line=~/SNP_REF_ALLELE_INDEX=(.+)$/){$ref_allele_index=$1;}
  if($line=~/SNP_OTHER_ALLELE_INDEX=(.+)$/){$snp_allele_index=$1;}
  if($line=~/SNP_ALLELES_INDEX=(.+)$/){$alleles_index=$1;}
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
  
  
}
#Check config file settings
if(($start_strain_index==-1)||($bp_pos_index==-1)||($chrom_index==-1)){
  print STDERR "One of the fields indice has a negative value \n"; 
  exit(1);
}
if(!(-d $snp_dir)|| !(-d $output_dir)){
  print STDERR "One of the following is not a directory:$snp_dir or $output_dir.Check the config file settings\n"; 
  exit(1);
}
$qh_getSourceid->execute($source);
if($qh_getSourceid->rows<=0){ #source not found in our database
   print STDERR "The specified source was not found in our database\n"; 
   exit(1);
}
else{($source_id)=$qh_getSourceid->fetchrow_array();}
#print "The source id of $source is $source_id\n";exit(0);
my($line,@linecontent,$chr,$chr_id);
if($snp_dir =~/(.*)\/$/){ 
        $snp_dir=$1; #format the directory name
}
if($output_dir =~/(.*)\/$/){ #format the directory name
        $output_dir=$1;
}
@chromosomes=(); 
if($all_in_one==0){
    open(IN,"$chr_list_file")or die "Can't open file:$!\n";
    @chromosomes=<IN>;
 }
 else{$chromosomes[0]=$snp_file_name;}

$missing_snp_allele_file=$output_dir."/".$source."_missing_snp_allele.txt";

open(CONT,">$missing_snp_allele_file");  # create the LOG output file
$snp_main_logfile=$output_dir."/multi_allele_load.log";
open(LOG,">$snp_main_logfile"); 
while(@chromosomes>0){   # process file
   $chr=shift(@chromosomes);chomp($chr);
  if($all_in_one==0){
     if($chr=~/^0(\d)$/){$chr_id=$1;} 
     else{$chr_id=$chr;}  $chr_temp=$chr;
     $chr_id=20 if($chr =~/X/i);$chr_id=21 if($chr =~/Y/i);
     $chr_id=22 if($chr =~/M/i);$chr_id=23 if($chr =~/U/i);
     $snp_file_name="";$chr_temp="M" if(($chr =~/M/i)&&($source_id==20));
     $snp_file_name=$snp_dir."/$file_prefix".$chr_temp."$file_sufix";
  }
  else{$snp_file_name="$snp_dir/$snp_file_name";}
  print "The file name is $snp_file_name\n";
  #next;
  if(!(-f $snp_file_name)){print "The file $snp_file_name does not exists\n";next;}
  open(SNP,"$snp_file_name") or die "Can't open $snp_file_name:$!";
  $line=<SNP>;chomp($line); #get the file header line
  $header=$line;
  my(@linecontent,$alleles,$ref_allele,$snp_allele,%strains);
  if($file_type==2){@linecontent= split(",",$line);} # commas separated file
  else{@linecontent= split(/\t/,$line);}              # default file type
  if($end_strain_index==-1){$end_strain_index= @linecontent;--$end_strain_index;}
  if($end_strain_index<$start_strain_index){
     print "Bad strain index $end_strain_index<$start_strain_index\n";next;
  }
  $i=$start_strain_index;#now check if these strains already exists into our database
  $strain_file="$output_dir/strains.txt";
  $snp_main_file=$output_dir."/chr".$chr."_snp_main_file.txt";
  #next if(-f $snp_main_file);
  $snp_mult_allele="$output_dir/chr".$chr."_snp_w_multiallele_file.txt";
  open(ST,">$strain_file");
  while($i<=$end_strain_index){ $strain=$linecontent[$i];
      $strain_name=$linecontent[$i];$strain_id=0;$strain_name=~s/^\s+//;$strain_name=~s/\s+$//;
      $strain_name=~s/,$//;
      $qh_getStrain->execute($strain_name);#check first in snp_strain table
      if($qh_getStrain->rows>0){
         ($strain_id)=$qh_getStrain->fetchrow_array();
         $qh_getStrainId->execute($strain_name); #check in synonym table
         if($qh_getStrainId->rows<=0){
           $qh_insertStrainSyn->execute($strain_id,$strain_name);
         }
      }
      else{ #new strain or synonymous
         $qh_getStrainId->execute($strain_name); #check first in snp_synonym
         if($qh_getStrainId->rows>0){($strain_id)=$qh_getStrainId->fetchrow_array();} #strain exists
         else{ #insert new strain
            $qh_getMaxStrainid->execute();
            if($qh_getMaxStrainid->rows>0){
               ($strain_id)=$qh_getMaxStrainid->fetchrow_array();
               ++$strain_id;
             }
           # print "$strain_name\t$strain_id needs to be inserted\n";
            $qh_insertStrain->execute($strain_id,$strain_name);$qh_getStrain->execute($strain_name);
            if($qh_getStrain->rows>0){
               ($strain_id)=$qh_getStrain->fetchrow_array();
                $qh_insertStrainSyn->execute($strain_id,$strain_name);
            }
         }
      }
      $qh_get_strain_source->execute($source_id,$strain_id);
      if($qh_get_strain_source->rows<=0){ ##insert this strain for this source if does not already exists 
         $qh_insertStrainBySource->execute($source_id,$strain_id);
       }
     $strains{$i}=$strain_id; #a dictionary mapping strain index to strain id    
     print ST "$strain\t$strain_id\n";
     if($skip_one==0){++$i; }
     else{$i+=2;}
    
   } 
   close(ST);
   if($ref_allele_index<0){$ref_allele_index= get_b6_index($line,$file_type);}
   open(MULT,">$snp_mult_allele");open(MAIN,">$snp_main_file");     
   # create the SNP main output file
   if($file_type==2){print MAIN "snpID,alleles,error_flag,$header\n";}
   else{print MAIN "snpID\talleles\terror_flag\t$header\n";}
   print MULT "snpID\tsnpAllele\tsource_id\n";
   #load existing SNPs from this chromosome into memory.
   $count=0;
   while(<SNP>){
      chomp($_);$count +=1;$id="";my @snpid;$rsid="";$snp_accession="";$bp_pos;$local_id=0;
      $line=$_;
      if($file_type==2){@linecontent= split(",",$_);} # commas separated file
      else{@linecontent= split(/\t/,$_);} # default file type
      if(@linecontent>0){
         $chr_id=0;$local_b6allele;$local_snpallele;$alleles="";$b6_al="";$snp_al="";
         if($alleles_index>=0){
            $alleles=$linecontent[$alleles_index];$alleles=~s/^\s+//;$alleles=~s/\s+$//;
            if($alleles=~/^(A|C|T|G){1}\W(A|C|T|G){1}$/i){$b6_al=$1;$snp_al=$2;}
         }
         else { $b6_al = $linecontent[$ref_allele_index];
           if($snp_allele_index<0){
              $snp_al=get_snp_allele($start_strain_index,$ref_allele_index,$line,$file_type);}
           else{$snp_al=$linecontent[$snp_allele_index];}
         }
         $bp_pos=$linecontent[$bp_pos_index];$chr=$linecontent[$chrom_index];$chr=~s/^chr//i;
         if($chr=~/^0(\d)$/){$chr_id=$1;} 
         else{$chr_id=$chr;}  $chr_temp=$chr;$local_id=0;
          $chr_id=20 if($chr =~/X/i);$chr_id=21 if($chr =~/Y/i);
          $chr_id=22 if($chr =~/M/i);$chr_id=23 if($chr =~/U/i);
         $chr_temp="M" if(($chr =~/M/i)&&($source_id==20));
         if($bp_pos=~/\d+/){#get the local id of this SNP
            $qh_get_snpid->execute($chr_id,$bp_pos);
            if($qh_get_snpid->rows>0){
               ($local_id)=$qh_get_snpid->fetchrow_array();
            }
            if(!($snp_al=~/A|G|T|C/i)||!($b6_al=~/A|G|T|C/i)){print CONT "$line\n";next;}
            $local_b6allele="-"; $local_snpallele="-";$error_flag=0;
            %allele_map=();#set flag to detect cases of genotypes with multi alleles
            if($snp_al=~/,/){ @all=split(",",$snp_al); $error_flag=16;
                 while(@all>0){
                    $geno_allele=shift(@all);$geno_allele=~s/\s+//g;$allele_map{$geno_allele}=0;
                 }
            }#generate genotype file
            $i=$start_strain_index;$socre=-1;
            while($i<=$start_strain_index){ 
                  $allele=$linecontent[$i];$strain_id=$strains{$i};if($skip_one>0){$score=$linecontent[$i+1];}
                  $allele=~s/\s+//g;$allele_map{$allele}+=1 if(!($allele=~/$b6_al/i)&&($allele=~/A|G|T|C/i));
                  if($skip_one==0){++$i;}else{$i+=2;}
            }
            $alleles="$b6_al/$snp_al";
            if($error_flag==16){ #get the allele found in most strains
                 $pop_allele=""; $max_count=0;
                 for my $allele(sort keys(%allele_map)){
                     if($max_count==0){
                        $pop_allele=$allele; $max_count=$allele_map{$allele};
                     }
                     else{
                        if($allele_map{$allele}>$max_count){$pop_allele=$allele; $max_count=$allele_map{$allele};}
                      }
                 }
                $alleles="$b6_al/$pop_allele";
             }
             if($file_type==2){print MAIN "$local_id,$alleles,$error_flag,$_\n";}
             else{print MAIN "$local_id\t$alleles\t$error_flag\t$_\n";}
             print MULT "$local_id\t$snp_al\t$source_id\n" if($error_flag==16);
             if($count%500000==0){print"$count processed \n" ;}
           } #end of bp_pos check
        } #end of @linecontent>0 check
    } #end of this chromosome SNP file loop
    close(SNP);close(MAIN);close(MULT);
    #now load and update error log if we have case with multiple genotype
    if(-f $snp_mult_allele){
         $linecont = `wc -l $snp_mult_allele`;chomp($linecont);
         if($linecont=~/^(\d+)\s$snp_mult_allele/){$linecont=$1;
            if($linecont>1){
               $qh_create_snp_mult->execute();$qh_delete_multiAlleles->execute($source_id);
               $qh_load_multiAlleles->execute("$snp_mult_allele");
               $qh_multiAlleles_rowcount->execute($source_id);
               if($qh_multiAlleles_rowcount->rows>0){
                  ($row_count)=$qh_multiAlleles_rowcount->fetchrow_array();
                  if(($linecont-$row_count)==1){
                     $qh_drop_error_temp->execute();$qh_create_error_temp->execute();
                     $qh_insert_error_temp->execute();$qh_update_error_temp->execute();
                     $qh_insert_error->execute();
                    }
                  else{print LOG "$snp_mult_allele was not loaded properly.\n";
                       print LOG "File lines:$linecont- table rows:$row_count\n";
                       $qh_delete_multiAlleles->execute($source_id);
                   }
               }
             }
         }
      $qh_analyze->execute();
    }
    print "Done with chromosome $chr\n"; 
 } #end of chromosome loop
close(IN);close(CONT);close(LOG);
print"program completed\n";
 
exit(0);

