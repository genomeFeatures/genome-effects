#!/usr/bin/perl

#*************************************************************************
# This script downloads ESTs from ftp.ncbi.nih.gov/repository/dbEST
# Author :       Lucie Hutchins Bioinformatic Software Engineer
# Date created: 10/2010
#
q{ dbEST.reports.000000.nnn.gz
   This is a compressed file containing the full reports for all publicly
   available entries in dbEST. 
   nnn is a volume number. 
   volume is defined by range of id_est (1-1000000,1000001-2000000,...)
   volume is updated only when at least one record changed or we have a 
   change in dumping software, otherwise volume from the previous dump 
   remains here.
};
#*************************************************************************

use Net::FTP;
use File::stat;
use Time::localtime;

$tm = localtime;
my ($mday, $mon, $yday) = ($tm->mday, ($tm->mon)+1, ($tm->year)+1900);
my($config_path)=shift(ARGV);
if(!$config_path){
   print "You must provide the path to Graber Lab dbEST mirror (/raid/seq/mirror_dbEST/config)\n";
   exit(1);
}

$config_path=~s/\/\s*$//; #remove the trailing back slash if exists
$config_file="$config_path/est_config.txt";
if(!(-f $config_file)){
   print"The config file est_config.txt is missing under $config_path\n";
   exit(1);
}
open(CONF,"$config_file");
if(!CONF){
   print"Could not open the config file $config_file :$!\n";
   exit(1);
}
@config_filecontent=<CONF>;
close(CONF);
my %config;
while(@config_filecontent>0){
   $line=shift(@config_filecontent);
   chomp($line);
   my($global_variable,$data)=split(",",$line);
   if($data){
     $config{"$global_variable"}="$data";
   }
}

# 
# Set up NCBI data path and ftp login info
# hostname, user, pass and remote dir
#
$host   = $config{"ncbi_host"};  $user   = $config{"ncbi_user"};
$rmtdir = $config{"ncbi_rmtdir"};$pass   = $config{"ncbi_pass"};
$report_prefix=$config{"ncbi_report_prefix"};

#
# Set up Graber Lab dbEST mirror data path
#
$working_dir = $config{"gb_working_dir"};    # path to local EST data
$report_dir  = $config{"gb_report_dir"};     # directory storing ESTs report files
$log_dir     = $config{"gb_log_dir"};        # set up a log file for this update
$temp_dir    = $config{"gb_temp_dir"}; 
$config_dir  = $config{"gb_config_dir"}; 
$lastmodified_since=$config{"lastmodified_since"};
#$gb_user=$config{"gb_db_user"};$gb_pass=$config{"gb_db_pass"};
#$gb_db=$config{"gb_db"};$gb_host=$config{"gb_db_host"};

# format the directory paths
$working_dir =~s/\/\s*$//; $report_dir=~s/\/\s*$//;$log_dir=~s/\/\s*$//;
if(!(-d $working_dir)){
   if(!mkdir($working_dir,0777)){
       print "Could not create $working_dir :$!\n";
       exit(1);
   }
}
if(!(-d $log_dir)){
   if(!mkdir($log_dir,0777)){
       print "Could not create $log_dir :$!\n";
       exit(1);
   }
}
my $update_logfile="$log_dir/update_log_"."$mon$mday$yday".".txt";
my $files2loadIndb="$log_dir/files2loadIndb_"."$mon$mday$yday".".txt";
#my $lastmodified_since= "$log_dir/last_modified.txt";
my %lastmodified_hash; #stores the mapping of each file and its last modified time
open(LOG,">$update_logfile");
if(!LOG){ 
   print "Could not opened the log file $update_logfile: $!\n";
   exit(1);
}
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "NCBI download started on $yday-$mon-$mday Time: $hour:$min:$sec\n";
open(DL,">$files2loadIndb");
if(!DL){ 
   print LOG "Could not opened the log file $files2loadIndb: $!\n";
   exit(1);
}
if(!(-d $report_dir)){
   if(!mkdir($report_dir,0777)){
       print LOG "Could not create $report_dir :$!\n";
       exit(1);
   }
}

if(!(-d $temp_dir)){
   if(!mkdir($temp_dir,0777)){
       print LOG "Could not create $temp_dir :$!\n";
       exit(1);
   }
}
chdir($report_dir); #we will run ftp transfer from the reports directory
#get the list of all the dbEST.reports.000000.nnn.gz from NCBI
$ftp = Net::FTP->new("$host", Debug => 0);
$ftp->login($user,$pass);$ftp->cwd($rmtdir);
$ftp->binary;@rls = $ftp->ls;$ftp->close;
if(@rls<=0){
   print LOG "Could not establish a data transfer chanel with NCBI:$!\n";
   exit(1);
}
@report_files= grep(/$report_prefix/,@rls); #get all the zipped report into and array
if(@report_files<=0){
  print LOG "No reports found NCBI:$!\n";
  exit(1);
}
print LOG "List of report files found NCBI:\n";
%current_lastmodified=();
if(-f $lastmodified_since){
    open(AL,"$lastmodified_since");
    if(AL){
          @lastmodif=<AL>;
          while(@lastmodif>0){
                $line=shift(@lastmodif); chomp($line);($report,$time)=split(",",$line);
                $current_lastmodified{"$report"}="$time";
          }
     }
   close(AL);
}
#
#For each report from NCBI,
#  get its last modified time and
#      check if file exits in mirror
#      set download flag to 1 if file exists and has new modifications
#      set downlod flag to 1 if file is new
#      dowload file if flag set
#      update download status and log files
#
#
open(LM,">$lastmodified_since");
while(@report_files>0){
  $reportfile=shift(@report_files);chomp($reportfile);
  $ftp = Net::FTP->new("$host", Debug => 0);
  $ftp->login($user,$pass);$ftp->cwd($rmtdir);$ftp->binary; $file_exists=0;$todownload=0;
  $lastmodified_time=$ftp->mdtm($reportfile); #get the time this file was last modified
  if(exists($current_lastmodified{"$reportfile"})){ #file exists then check if hasen't been modified since last download
     if($current_lastmodified{"$reportfile"} ne "$lastmodified_time"){
        $todownload=1; $command="mv $reportfile $temp_dir/$reportfile";system($command); 
        #save the old copy in case the download goes wrong
     }
  }
  else{$todownload=1;} #new file
  if($todownload){
      $ftp->get("$reportfile"); #download fresh copy
      $ftp->close;print "$reportfile downloaded\n";
      #now check the status return by the download attemp
      $command="gunzip -t -v $reportfile 2> status_file.txt";system($command);$status_id=0; $status_name="Failed";
      $toLoad=0; #flag set to 1 if file is to be loaded in our database
      open(ST,"status_file.txt");
      if(ST){@status_list=<ST>; 
        while(@status_list>0){$line=shift(@status_list);
              if($line=~/\s+OK/){$status_name="OK";
                 print DL "$reportfile\n";
                 print LM "$reportfile,$lastmodified_time\n" if(LM);
                 if(exists($current_lastmodified{"$reportfile"})){
                     $command="rm $temp_dir/$reportfile"; system($command);
                 }
               }
               else{ #download attempt failed, remove the failed file then get back the old copy 
                   $command="rm $reportfile"; system($command);
                   if(exists($current_lastmodified{"$reportfile"})){
                      $command="mv $temp_dir/$reportfile $reportfile";
                      system($command);
                      $lastmodified_time=$current_lastmodified{"$reportfile"};
                      print LM "$reportfile,$lastmodified_time\n" if(LM);
                   }
               }
         }
         close(ST);
       }
       print LOG "File: $reportfile,Download_status: $status_name\n";
   }
   else{
      print LM "$reportfile,$lastmodified_time\n" if(LM);
   }  
}
#now update the config file to include both the log file and the file2download files
open(CONF,">$config_file");
if(CONF){
   while(($global_variable,$variable_value)=each %config){
       print CONF "$global_variable,$variable_value\n";
   }
   print CONF "update_log_file,$update_logfile\n";
   print CONF "files2loadIndb,$files2loadIndb\n";
   print LOG "ConfigFile_status: Updated\n";
}
close(CONF);
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "NCBI download endeded on $yday-$mon-$mday Time: $hour:$min:$sec\n";
close(LOG);

#now remove all temp files and directories
if( -f "$report_dir/status_file.txt"){
   $command="rm $report_dir/status_file.txt"; system($command);
}
if(-d "$temp_dir"){
   $command="rm -rf $temp_dir"; system($command);
}
print "Program complete\n";
exit(0);

