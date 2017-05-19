#!/usr/bin/perl

use DBI;
use Time::localtime;

#*****************************************************************
# This script downloads all the annotation tables for each
# organism found in ucsc_geneannotation_date congig file generated 
# from getCurrentOrgansims.pl script
# 
# Input : a config file (ucsc_geneannotations_xx-xx-xxxx.txt) with the following fields
#       1.organism_group,
#       2.organims_name, 
#       3.ucsc_db, 
#       4.annotation tables list separated by ":"
#
# Output : a config file downloaded_annotations.log with the following fields
#       1.organism_group,
#       2.organims_name, 
#       3.annotation table name, 
#       4.downloaded annotation file name
 
# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : June 2010
# Usage: perl download_ucsc_annotations.pl configfile 
#**************************************************************
#


$config=shift||"/scratch/data/ucsc/ucsc_load_process/configuration/prog_generated_config.txt";
chomp($config);
open(CONF,"$config");
if(!CONF){
   print"Usage: perl run_update_pipeline.pl config file\n";
   exit(1);
}
@filecontent=<CONF>; my($configpath,);
my $path="/scratch/data/ucsc"; # path to gene annotations
my $load_base_path="/scratch/data/ucsc/ucsc_load_process";
$host="genome-mysql.cse.ucsc.edu"; $user="genome";
my $configFile;
while(@filecontent>0){
   $line=shift(@filecontent);chomp($line);
   my($variable,$base_path)=split(",",$line);
   if($variable=~/ucsc_load_base/){
      $load_base_path=$base_path;
   }
   elsif($variable=~/ucsc_org_base/){
      $path=$base_path;
   }
   elsif($variable=~/annot_config/){
      $configFile=$base_path;
   }
   elsif($variable=~/ucsc_db_host/){
      $host=$base_path;
   }
   elsif($variable=~/ucsc_db_usr/){
      $user=$base_path;
   }
}
close(CONF);
my $connection="mysql -h$host -u$user -A";

my $configpath="$load_base_path/configuration";
my $downloadedAnnot="$configpath/downloaded_annotations.log"; #stores the list downloaded tables
my $mainConfig ="$configpath/prog_generated_config.txt";
my $processlog="$configpath/data_load_log.log";
open(LOG,">>$processlog");
if(LOG){
  $tm = localtime;
  # my ($mday, $mon, $yday) = ($tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Starting download_ucsc_annotations.pl program : $mday/$mon/$yday @ $hour:$min:$sec\n";
  print LOG "\n*************************************************************\n";
  chomp($configFile);
  if(!$configFile){
      print LOG "download_ucsc_annotations.pl usage was incorect, missing the config file argument\n";
      close(LOG);exit(1);
   }
   open(IN,"$configFile");open(DN,">$downloadedAnnot");
   if(IN){ @organisms=<IN>;
      foreach my $line(@organisms){
         if($line){ chomp($line);
            my($organism_group,$organism,$db,$tables)=split(",",$line);
            chomp($organism);
            $organism=lc($organism);$organismDir=$organism;
            $organismDir=~s/\s+/-/g; 
            if(!(-d "$path/$organism_group")){mkdir("$path/$organism_group",0777)||print LOG "$!\n";}
            if(!(-d "$path/$organism_group/$organismDir")){ 
               mkdir("$path/$organism_group/$organismDir",0777)||print LOG "$!\n";
            }
            if((-d "$path/$organism_group")&& (-d "$path/$organism_group/$organismDir")){
               if($organism){
                  @filelist=`ls $path/$organism_group/$organismDir`;
                  if(@filelist>0){  # under the archive 
                     if(!(-d "$path/$organism_group/$organismDir/archive")){
                           mkdir("$path/$organism_group/$organismDir/archive",0777)||print LOG "$!\n";
                      }@archivelist=`ls $path/$organism_group/$organismDir/archive`;
                      if(@archivelist>0){
                         $commands=("rm $path/$organism_group/$organismDir/archive/*");system($commands);
                      }
                      $commands=("mv $path/$organism_group/$organismDir/$db* $path/$organism_group/$organismDir/archive/");
                      system($commands);
                  }
                  @tablelist=split(":",$tables);
                  while(@tablelist>0){
                       $table=shift(@tablelist);
                       if($table){
                          if(!($table=~/refFlatSquishMulti/)){
                              #### will remove the row below
                              $file="$db"."_".$table."_"."$mon-$yday";
                              $temp_file="$path/$organism_group/$organismDir/archive/$file";
                              $filename="$path/$organism_group/$organismDir/$file"; 
                              if(-f "$path/$organism_group/$organismDir/archive/$file"){
                                 $lines=`wc -l $temp_file`;
                                 $linecount=0;$lines=~s/^\s+//; #check if the file is empty
                                 if($lines=~/^\s*(\d+)\s+$temp_file/){$linecount=$1;}
                                 system "mv $path/$organism_group/$organismDir/archive/$file $filename" if($linecount>0);
                              }
                              if(!(-f $filename)){$fields="*";$where="";
                                 if($table eq "gbCdnaInfo"){
                                    $fields="acc as accession,version,moddate,direction,source,library,mrnaClone,";
                                    $fields.="sex,tissue,development,cell,description,geneName,author,gi,mol";
                                    $where=" where type=\"EST\"";
                                 }
                                if(!($table=~/sex|source|library|mrnaClone|tissue|development|cell|description|geneName|author|gbCdnaInfo/)){
                                    if(!($organism=~/Mouse|Human/)||(($organism=~/Mouse|Human/)&& !($table=~/ensGene/))){
                                       $commands=("$connection --quick -e 'use $db; select $fields from $table  $where' > $filename");
                                       system($commands); print "$table $db downloaded\n";
                                    }
                                }
                              }
                           } #end of  if(!($table=~/refFlatSquishMulti/)){
                           if(!($table=~/sex|source|library|mrnaClone|tissue|development|cell|description|geneName|author|gbCdnaInfo/)){
                               if(!($organism=~/Mouse|Human/)||(($organism=~/Mouse|Human/)&& !($table=~/ensGene/))){
                                 print DN "$organism_group,$organism,$db,$table,$filename\n";
                               }
                           }
                      }#end of if($table){
                      else{
                          $filename="$path/$organism_group/$organismDir/$table";
                          print DN "$organism_group,$organism,$db,$table,$filename\n";
                      }
                 } #end of while(@tablelist>0){
               } #end of if($organism)
            } #end of path validation
            else{print LOG "$path/$organism_group/$organismDir bad path: $!\n";}
         } #end of line validation -- if($line){
      } #end of foreach organism -- foreach my $line(@organisms){
    } #end of config file input validation-- if(IN){
    else{print LOG "Failed to open $configFileconfig file: $!\n";}
 } #end of log file validation

open(OUT,">>$mainConfig");
if(OUT){
  print OUT "ucsc_table_list,$downloadedAnnot\n";
}
close(OUT);
close(DN);
  $tm = localtime;
 my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);

  print LOG "\n*************************************************************\n";
  print LOG "download_ucsc_annotations.pl program complete: $mday/$mon/$yday @ $hour:$min:$sec\n";
  print LOG "\n*************************************************************\n";
close(LOG);



