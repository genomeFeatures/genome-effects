#!/usr/bin/perl

use DBI;
use Net::FTP;
use File::stat;
use Time::localtime;

#****************************************************************************************
# This script downloads ESTs from ftp://hgdownload.cse.ucsc.edu/goldenPath/XXX/bigZips
# This is a compressed file containing ESTs in GenBank. UCSC updates this sequence data once a
# week via automatic GenBank updates.
# This script is ran after load_geneannotations.pl
# 
# Input : a config file (current_organism_version) 
#         with the following fields
#       1.organism_group,
#       2.organims_name, 
#       3.ucsc_db, 
#       4.version
# Download path: ftp://hgdownload.cse.ucsc.edu/goldenPath/
# ftp to hgdownload.cse.ucsc.edu
# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : December 2011
# Usage: perl program configfile 
#**************************************************************
#

$config=shift||"/scratch/data/ucsc/ucsc_load_process/configuration/Configuration";
chomp($config);open(CONF,"$config");
if(!CONF){
   print"Usage: perl getESTsFromUCSC.pl config_file\n";
   exit(1);
}
@filecontent=<CONF>;
close(CONF);
my $est_path="/raid/seq/mirror_dbEST/ucsc"; # path to genome data
my $load_base_path="/scratch/data/ucsc/ucsc_load_process";

my  $genome_path="/raid/seq";
my %config;
while(@filecontent>0){
   $line=shift(@filecontent);
   chomp($line);
   my($global_variable,$data)=split(",",$line);
   if($data){
     $config{"$global_variable"}="$data";
   }
}
# 
# Set up UCSC data path and ftp login info
# hostname, user, pass and remote dir
#
$ftp_host  = $config{"ucsc_ftp_host"};$ftp_path = $config{"ucsc_ftp_file_path"};
$est_file=$config{"ucsc_est"};
$load_base_path=$config{"ucsc_load_base"};
$user   = $config{"ucsc_user"};$pass = $config{"ucsc_pass"};
#
#Set graber db info
#
$graber_db=$config{"graber_db_dev"};
$graber_db_host=$config{"graber_db_hostdev"};
$graber_db_usr=$config{"graber_db_usr"};
$graber_db_pass=$config{"graber_db_pass"};

my $dbh = DBI->connect("DBI:mysql:$graber_db:$graber_db_host",$graber_db_usr,$graber_db_pass);
my $getOrg="select organism_id,organism_name from organism";  #get current organisms list
my $getOgrv="select organism_version_id from organism_version where ucsc_db=?";

my $createESTsTemp="create temporary table if not exists sequence_temp(";
   $createESTsTemp.="accession char(12),flag tinyint,sequence text,index(accession))";
my $creategbCdnaInfo= "CREATE temporary TABLE IF NOT EXISTS gbCdnaInfo_temp(
                          accession char(12),version smallint unsigned default 0,moddate date default '0000-00-00',
                          direction enum('5','3','0'),source int unsigned default 0,
                          library int unsigned default 0,mrnaClone int unsigned default 0,
                          sex int unsigned default 0,tissue int unsigned default 0,development int unsigned default 0,
                          cell int unsigned default 0,description int unsigned default 0,
                          geneName int unsigned default 0,author int unsigned default 0,gi int unsigned default 0,
                          mol  enum('DNA','RNA','ds-RNA','ds-mRNA','ds-rRNA','mRNA','ms-DNA','ms-RNA','rRNA','scRNA','snRNA','snoRNA','ss-DNA','ss-RNA','ss-snoRNA','tRNA','cRNA','ss-cRNA','ds-cRNA','ms-rRNA'),index(accession),index(gi))";
$createTempTable="create temporary table if not exists est_info_temp (";
$createTempTable.="id int  default 0,name text,crc int unsigned default 0,index(id))";


my $dropTempTable=" drop temporary table if exists est_info_temp";
my $dropEstTemp="drop temporary table if exists sequence_temp";
my $dropgbCdnaTemp=" drop temporary table if exists gbCdnaInfo_temp";

my $loadTempTable="load data infile ? into table est_info_temp ignore 1 lines";
my $loadEstTable="load data infile ? into table sequence_temp";
my $loadgbCdnaTable="load data infile ? into table gbCdnaInfo_temp ignore 1 lines";

my $insertSexTemp=" insert into sex select ?,id,name,crc from est_info_temp where id>0";
my $insertlibTemp=" insert into library select ?,id,name,crc from est_info_temp where id>0";
my $insertSourceTemp=" insert into source select ?,id,name,crc from est_info_temp where id>0";
my $inserttissueTemp=" insert into tissue select ?,id,name,crc from est_info_temp where id>0";
my $insertdevTemp=" insert into development select ?,id,name,crc from est_info_temp where id>0";
my $insertcellTemp=" insert into cell select ?,id,name,crc from est_info_temp where id>0";
my $insertdescTemp=" insert into description select ?,id,name,crc from est_info_temp where id>0";
my $insertauhorTemp=" insert into author select ?,id,name,crc from est_info_temp where id>0";
my $insertgeneNameTemp=" insert into geneName select ?,id,name,crc from est_info_temp where id>0";
my $insertmrnaTemp=" insert into mrnaClone select ?,id,name,crc from est_info_temp where id>0";

$qh_insertSexTemp=$dbh->prepare($insertSexTemp);
$qh_insertlibTemp=$dbh->prepare($insertlibTemp);
$qh_insertSourceTemp=$dbh->prepare($insertSourceTemp);
$qh_inserttissueTemp=$dbh->prepare($inserttissueTemp);
$qh_insertdevTemp=$dbh->prepare($insertdevTemp);
$qh_insertcellTemp=$dbh->prepare($insertcellTemp);
$qh_insertdescTemp=$dbh->prepare($insertdescTemp);
$qh_insertauhorTemp=$dbh->prepare($insertauhorTemp);
$qh_insertgeneNameTemp=$dbh->prepare($insertgeneNameTemp);
$qh_insertmrnaTemp=$dbh->prepare($insertmrnaTemp);

my $insertESTsTemp=" insert into sequence select ?,accession,flag,sequence from sequence_temp where flag!=-99";
my $updateEstsTemp=" update sequence_temp t,sequence s set t.flag=-99 where t.accession=s.qName and s.organism_id=?";

my $insertgbCdnaInfo=" insert into gbCdnaInfo select ?,accession,version,moddate, direction,source,";
   $insertgbCdnaInfo.=" library,mrnaClone,sex,tissue,development,cell,description,";
   $insertgbCdnaInfo.="geneName ,author,gi,mol from gbCdnaInfo_temp where version!=-99";
my $updategbCdnaInfo=" update gbCdnaInfo_temp t,gbCdnaInfo s set t.version=-99 where t.gi=s.gi and s.organism_version_id=?";
my $update_datalog= "update data_load_test set db_line_count=? where organism_version_id=? and table_name=?";

my $getEstRowCount="select count(*) as estcount from sequence_temp";
my $getCdnaRowCount="select count(*) as estcount from gbCdnaInfo_temp";
my $getTempRowCount="select count(*) as estcount from est_info_temp";


my $qh_getEstRowCount       = $dbh->prepare($getEstRowCount)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_getTempRowCount       = $dbh->prepare($getTempRowCount)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_getCdnaRowCount       = $dbh->prepare($getCdnaRowCount)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_orglist           = $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_getOgrv           = $dbh->prepare($getOgrv)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_createESTsTemp    = $dbh->prepare($createESTsTemp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_creategbCdnaInfo  = $dbh->prepare($creategbCdnaInfo)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_createTempTable  = $dbh->prepare($createTempTable)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_insertESTsTemp    = $dbh->prepare($insertESTsTemp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insertFromTemp    = $dbh->prepare($insertFromTemp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insertgbCdnaInfo  = $dbh->prepare($insertgbCdnaInfo)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insertTerms        = $dbh->prepare($insertTerms)or die "Couldn't prepare statement: " . $dbh->errstr; 

my $qh_deleteExistingRows= $dbh->prepare($deleteExistingRows)or die "Couldn't prepare statement: " . $dbh->errstr; 

my $qh_dropTempTable    = $dbh->prepare($dropTempTable)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_dropEstTemp    = $dbh->prepare($dropEstTemp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_dropgbCdnaTemp    = $dbh->prepare($dropgbCdnaTemp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_loadTempTable    = $dbh->prepare($loadTempTable)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_loadEstTable    = $dbh->prepare($loadEstTable)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_loadgbCdnaTable    = $dbh->prepare($loadgbCdnaTable)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_updateEstsTemp= $dbh->prepare($updateEstsTemp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_updategbCdnaInfo=$dbh->prepare($updategbCdnaInfo)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_update_datalog =$dbh->prepare($update_datalog)or die "Couldn't prepare statement: " . $dbh->errstr;

my $configpath        ="$load_base_path/configuration";
my $downloadedConfig ="$configpath/downloaded_annotations.log";
my $processlog        ="$configpath/data_load_log.log";

my %orgmap=();
#get the list of all organisms
  $qh_orglist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
  while(($org_id,$org_name)=$qh_orglist->fetchrow_array()){
        $org_name=lc($org_name);$org_name=~s/^\s+//;$org_name=~s/\s+$//;
        $orgmap{"$org_name"}=$org_id;
  }
#get list of term types

open(LOG,">>$processlog");
if(LOG){
  $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Starting getESTsFromUCSC.pl program :  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
  #load the downloaded config file into a map
  #$organismVersion{"$db"}{"gbCdnaInfo"};
  open(IN,"$downloadedConfig") or die ("Bad config file :$!\n");
  @configContent=<IN>;
  while(@configContent>0){
     $line=shift(@configContent);chomp($line);
     ($organism_group,$thisorganism,$db,$table,$filename)=split(",",$line);
     $organismVersion{"$db"}{"$table"}=$filename;
  }
  close(IN);
  #only download ESTs for the latest version of a given organism
  @organismv_list=`grep all_est $downloadedConfig| cut -d ',' -f '1 2 3 4 5' |sort|uniq`;
  @organism_list=`grep all_est $downloadedConfig| cut -d ',' -f 2 |sort|uniq`;
 # print @organismv_list;
  foreach my $organism(@organism_list){
    chomp($organism);  $db="";$organism_group="";@versionsCount=grep(/$organism/,@organismv_list);
    while($line=shift(@versionsCount)){chomp($line);
         if($line=~/,/){($organism_group,$thisorganism,$dbversion,$table,$filename)=split(",",$line);
            if($organism eq $thisorganism){
               if($db eq ""){$db=$dbversion;}
               else{ $db=~/(\d+)\s*$/;$prev=$1;$dbversion=~/(\d+)\s*$/;$curv=$1;$db=$dbversion if($prev<$curv);}
            }
          }
     }
    if($db ne ""){
      $organism=lc($organism);
      $organismDir=lc($organism); $organismDir=~s/^\s+//; $organismDir=~s/\s+$//;$organismDir=~s/\s+/-/;
      $genome_dir="$genome_path/$organism_group/$organismDir/$db";
      $est_path="$genome_dir/est";
      if(!(-d "$est_path")){mkdir("$est_path",0777)|| print LOG "$!\n";}
      $est_temp="$est_path/temp";
      if(-d "$est_temp"){system("rm -rf $est_temp");
        mkdir("$est_temp",0777)|| print LOG "$!\n";
       }
      else{mkdir("$est_temp",0777)|| print LOG "$!\n";}
      print "Processing $organism\n";
      # we will do the following for each organism:
      # download the est zipped main file, unzip it ,generate the files,then load
      # the file into sequence table
      $ftp_dir="$ftp_path/$db/bigZips";
      chdir($est_temp); #we will run ftp transfer from the est/temp directory
      $ftp = Net::FTP->new("$ftp_host", Debug => 0);
      $ftp->login($user,$pass);$ftp->cwd($ftp_dir);$ftp->binary;
      $ftp->get("$est_file"); #download fresh copy
      $ftp->close;
      print "$db $est_file downloaded\n";
      #now check the status return by the download attemp
      $command="gunzip -t -v $est_file 2> status_file.txt";system($command);$status_id=0; $status_name="Failed";
      open(ST,"status_file.txt"); $ucsc_est_file=$config{"ucsc_est_file"};
      if(ST){@status_list=<ST>; 
            while(@status_list>0){$line=shift(@status_list);
                   if($line=~/\s+OK/){$status_name="OK"; 
                      if(-f $ucsc_est_file){system("rm $ucsc_est_file");}
                       system("gunzip $est_file");
                      if(-f $ucsc_est_file){system("mv $ucsc_est_file ../");last;}
                   }
                   else{$command="rm $est_file"; system($command);}
            }
       }
       close(ST); chdir($est_path);#we will run fa2dat from est directory to generate a tab-delimited file
      #if(-f $ucsc_est_file){
      q{    $command="perl $genome_path/fa2dat.pl -i $ucsc_est_file -o est.dat -m -t";
          system($command);$temp_file="est_temp_file.txt";
          #now generate this file into a database format
          open(IN,"est.dat");open(OUT,">$temp_file");$linecount=0;
          if(OUT){
             while($line=<IN>){
                chomp($line);($accession,$seq)=split("\t",$line);
                if($accession=~/>(.+)\s+(\d+)\s*$/){$accession=$1; $flag=$2;
                    ++$linecount;
                   print OUT "$accession\t$flag\t$seq\n";
                }
             }
          }
          else{
             print LOG "Could not open $temp_file output file: ESTs not loaded for $organism build $db\n";
             next;
          }
          close(IN);close(OUT);
          system("rm est.dat");  #now load this into the data
          $path=`pwd`;chomp($path);$temp_file="$path/$temp_file"; 
          };
          #if(-f $temp_file){  #load ests sequence data $qh_dropEstTemp,$qh_insertESTsTemp ,$qh_createESTsTemp
          q{   $qh_dropEstTemp->execute() or die mysql_error();
              $qh_createESTsTemp->execute()  or die mysql_error();    
              $qh_loadEstTable->execute($temp_file)  or die mysql_error();
              $qh_getEstRowCount->execute()  or die mysql_error();
             if(($estcount)=$qh_getEstRowCount->fetchrow_array()){
                 
                 if($linecount==$estcount){
                    $org_name=lc($organism);$org_name=~s/^\s+//;$org_name=~s/\s+$//;
                    $org_id= $orgmap{"$org_name"};
                    $qh_updateEstsTemp->execute($org_id);
                    $qh_insertESTsTemp->execute($org_id);
                    print "ESTs loaded for $organism build $db : file line count=$linecount loaded= $estcount \n";
                 }else{
                    print LOG "Bad table load: ESTs not loaded for $organism build $db\n";
                    print LOG "We have $linecount in file vs $estcount loaded\n";
                 }
             } 
            };
             @tableList=("sex","source","library","mrnaClone","tissue","development","cell","description","geneName","author");
             $qh_getOgrv->execute($db);
             if(($organism_version_id)=$qh_getOgrv->fetchrow_array()){
                $temp_file=$organismVersion{"$db"}{"gbCdnaInfo"};
                $qh_dropgbCdnaTemp->execute() or die mysql_error();$rowcount=0;
                $qh_creategbCdnaInfo->execute() or die mysql_error();
                $qh_loadgbCdnaTable->execute($temp_file)  or die mysql_error();
                $qh_getCdnaRowCount->execute()  or die mysql_error();
                if(($rowcount)=$qh_getCdnaRowCount->fetchrow_array()){
                   $lines=`wc -l $temp_file`;$linecount=0;$lines=~s/^\s+//;
                   if($lines=~/^\s*(\d+)\s+$temp_file/){$linecount=$1;}
                   --$linecount; $name="gi";$type_id=$term_types{"gbCdnaInfo"};
                   if($linecount==$rowcount){
                      $qh_updategbCdnaInfo->execute($organism_version_id);
                      $qh_insertgbCdnaInfo->execute($organism_version_id);
                      $qh_update_datalog->execute($rowcount,$organism_version_id,"gbCdnaInfo");
                      print "gbCdnaInfo loaded for $organism build $db\n";
                    }
                   else{
                      print LOG "Bad table load: $temp_file not loaded for $organism build $db\n";
                      print LOG "We have $linecount in file vs $rowcount loaded\n";
                   }
                 }

               $qh_dropgbCdnaTemp->execute() or die mysql_error();
               q{
               foreach my $table(@tableList){  
                     $temp_file=$organismVersion{"$db"}{"$table"};
                     $qh_dropTempTable ->execute() or die mysql_error();
                     $qh_createTempTable->execute() or die mysql_error();
                     $qh_loadTempTable->execute($temp_file)  or die mysql_error();
                     $qh_getTempRowCount->execute()  or die mysql_error();
                     if(($rowcount)=$qh_getTempRowCount->fetchrow_array()){
                         $lines=`wc -l $temp_file`;$linecount=0;$lines=~s/^\s+//;
                         if($lines=~/^\s*(\d+)\s+$temp_file/){$linecount=$1;}
                         --$linecount; $name="name";
                         if($linecount==$rowcount){
                           if($table eq "sex"){
                            $query="update est_info_temp t, sex s set t.id=0 where s.organism_version_id=$organism_version_id";
                            $query.=" and s.id=t.id"; $dbh->do($query);
                            $qh_insertSexTemp->execute($organism_version_id);
                            }
                           elsif($table eq "library"){
                            $query="update est_info_temp t, library s set t.id=0 where s.organism_version_id=$organism_version_id";
                            $query.=" and s.id=t.id"; $dbh->do($query);
                            $qh_insertlibTemp->execute($organism_version_id);
                           }
                           elsif($table eq "source"){
                              $query="update est_info_temp t, source s set t.id=0 where s.organism_version_id=$organism_version_id";
                              $query.=" and s.id=t.id"; $dbh->do($query);
                              $qh_insertSourceTemp->execute($organism_version_id);
                           }
                           elsif($table eq "tissue"){
                              $query="update est_info_temp t, tissue s set t.id=0 where s.organism_version_id=$organism_version_id";
                              $query.=" and s.id=t.id"; $dbh->do($query);
                              $qh_inserttissueTemp->execute($organism_version_id);
                            }
                           elsif($table eq "development"){
                              $query="update est_info_temp t,development s set t.id=0 where s.organism_version_id=$organism_version_id";
                              $query.=" and s.id=t.id"; $dbh->do($query);
                              $qh_insertdevTemp->execute($organism_version_id);
                           }
                           elsif($table eq "cell"){
                              $query="update est_info_temp t, cell s set t.id=0 where s.organism_version_id=$organism_version_id";
                              $query.=" and s.id=t.id"; $dbh->do($query);
                              $qh_insertcellTemp->execute($organism_version_id);}
                           elsif($table eq "description"){
                              $query="update est_info_temp t,description s set t.id=0 where s.organism_version_id=$organism_version_id";
                              $query.=" and s.id=t.id"; $dbh->do($query);
                              $qh_insertdescTemp->execute($organism_version_id);
                            }
                           elsif($table eq "author"){
                              $query="update est_info_temp t, author s set t.id=0 where s.organism_version_id=$organism_version_id";
                              $query.=" and s.id=t.id"; $dbh->do($query);
                              $qh_insertauhorTemp->execute($organism_version_id);
                           }
                           elsif($table eq "geneName"){
                             $query="update est_info_temp t, library s set t.id=0 where s.organism_version_id=$organism_version_id";
                             $query.=" and s.id=t.id"; $dbh->do($query);
                             $qh_insertgeneNameTemp->execute($organism_version_id);
                           }
                           elsif($table eq "mrnaClone"){
                             $query="update est_info_temp t, mrnaClone s set t.id=0 where s.organism_version_id=$organism_version_id";
                             $query.=" and s.id=t.id"; $dbh->do($query);
                             $qh_insertmrnaTemp->execute($organism_version_id);
                            }
                           $qh_update_datalog->execute($rowcount,$organism_version_id,$table);
                         }
                         else{
                           print LOG "Bad table load: $temp_file not loaded for $organism build $db\n";
                           print LOG "We have $linecount in file vs $rowcount loaded\n";
                        }
                     } #end of load check
                    $qh_dropTempTable->execute() or die mysql_error();
                 } #end of foreach
               };
            } #end of organism version check
         # } #end of temp_file check
         # else{ print LOG "$temp_file does not exist: ESTs not loaded for $organism build $db\n";}
       } # end of if(-f $ucsc_est_file)
       #if(-d "$est_temp"){system("rm -rf $est_temp");}
   } # end of if($db ne "")
  } # end of foreach my $organism
#}
#optimize tables index

$command="mysql -h$graber_db_host -u$graber_db_usr -p$graber_db_pass -e 'use $graber_db;";
$command.=" ANALYZE TABLE cell,development,description,geneName,library,mrnaClone,sex,sequence,author,source,tissue;'";
system($command);


 close(LOG);
  print "\n Program complete\n";
  exit(0);

