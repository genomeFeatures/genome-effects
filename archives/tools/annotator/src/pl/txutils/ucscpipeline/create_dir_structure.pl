#!/usr/bin/perl

use Net::FTP;
use File::stat;
use Time::localtime;

#*****************************************************************
# This script  create a local directory structure for the organisms
#

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
# Date : January 2011
# Modified: August 2012
# Usage: perl program configfile 
#**************************************************************
#
$config=shift||"/scratch/data/ucsc/ucsc_load_process/configuration/Configuration";
chomp($config);open(CONF,"$config");
if(!CONF){
   print"Usage: perl create_dir_structure.pl config file\n";
   exit(1);
}
@filecontent=<CONF>;
close(CONF);
my $genome_path="/data/seq"; # path to genome data
my $load_base_path="/scratch/data/ucsc/ucsc_load_process";
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
$report_prefix=$config{"ucsc_file_prefix"};$genome_path= $config{"genome_path"};
$load_base_path=$config{"ucsc_load_base"};$ucsc_file_prefix=$config{"ucsc_file_prefix"};
$user   = $config{"ucsc_user"};$pass = $config{"ucsc_pass"};

my $configpath="$load_base_path/configuration";
my $mainConfig ="$configpath/current_organism_version";
my $progConfig="$configpath/prog_generated_config.txt";
my $processlog="$configpath/data_load_log.log";

open(LOG,">>$processlog");

if(LOG){
  $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Starting create_dir_structure.pl program :  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
  chomp($progConfig);
  if(!$progConfig){
      print LOG "create_dir_structure.pl usage was incorect, missing the config file argument\n";
      close(LOG);
      exit(1);
   }
  open(IN,"$mainConfig");my %annotmap;
  if(!IN){
      print  "Bad config file argument\n";
      close(IN);
      exit(1);
   }
 #get 
  #print  "Bad config file argument\n"; exit(0);
  %orgProcessed=(); %organisms_list=();
  @organisms=<IN>;
  foreach my $line(@organisms){
      if($line){chomp($line);
          my($organism_group,$organism,$db,$version)=split(",",$line);chomp($organism);
          $organisms_list{$organism}{$db}="$line";
      }
  }
 for my $organism(%organisms_list){
     $version=0;$currentdb="";
     for my $db(%{$organisms_list{"$organism"}}){
         $line=$organisms_list{$organism}{$db};
         if($line){
            chomp($line); my($organism_group,$organism,$db,$dbversion)=split(",",$line);
             if($db=~/(\d+)\s*$/){
                  if($version==0){$currentdb=$db;$version=$1;}
                  else{
                      if($version<$1){
                         $currentdb=$db;$version=$1;}
                  }
             }else {$currentdb=$db;}
            chomp($organism);$organism=lc($organism);$organismDir=lc($organism);
            $organismDir=~s/\s+/-/g; $genome_dir="$genome_path/$organism_group/$organismDir/$db";
            $organism_base="$genome_path/$organism_group/$organismDir"; $genome_mol_dir="$genome_dir/mol";
            $genome_data_dir="$genome_dir/dat"; $genome_temp="$genome_dir/temp";
            if(!(-d "$genome_path/$organism_group")){mkdir("$genome_path/$organism_group",0777)|| print LOG "$!\n";}
            if(!(-d "$genome_path/$organism_group/$organismDir")){ 
               mkdir("$genome_path/$organism_group/$organismDir",0777)|| print LOG "$!\n";
            }
            if(!(-d "$genome_dir")){ 
               mkdir("$genome_dir",0777)||print LOG "$!\n";
            }
            if(!(-d "$genome_data_dir")){ 
                    mkdir("$genome_data_dir",0777)|| print LOG "$!\n";
                }
            if(!(-d "$genome_mol_dir")){ 
                    mkdir("$genome_mol_dir",0777)|| print LOG "$!\n";
                }
            #NOW get list of files under dat/ directory
           opendir MYDIR,"$genome_data_dir";
           @dirContent = grep /\.dat$/,readdir MYDIR; #Genome already downloaded for this build
           close(MYDIR); 
           opendir MYDIR,"$genome_mol_dir";
           @molContent = grep /\.fa$/,readdir MYDIR; #Genome already downloaded for this build
           close(MYDIR);                             
           next if((@dirContent>0)&&(@molContent>0));
           print  "Processing $organism_group,$organism,$db\n";
           # we will do the following for each orgaism:
           # download the genome, generate the .dat files, generate symbolic link current_genome
          if(!(-d $genenome_temp)){mkdir("$genome_temp",0777);}
          chdir($genome_temp); #we will run ftp transfer from the reports directory
                    #get all t from NCBI
          #$wkdir=`pwd`; chomp($wkdir);
          #print "running ftp from $wkdir\n";
          $ftp_dir="$ftp_path/$db/chromosomes";
          $ftp = Net::FTP->new("$ftp_host", Debug => 0) ||print LOG "could not connect\n";
          $ftp->login($user,$pass) ||print LOG"could not login $user,$pass :".$ftp->message()."\n";
          $ftp->cwd($ftp_dir) ||print LOG "Bad path $ftp_dir\n";
          $ftp->binary;@rls = $ftp->ls;$ftp->close;
          if(@rls<=0){
                    print LOG "File transfer not successful for $ftp_dir:$!\n";next;
           }
           @report_files= grep(/$report_prefix/,@rls); #get all the zipped report into and array
           if(@report_files<=0){
                print LOG "No reports found $ftp_dir:$!\n";next;
            }
            $filecount=0;
            while(@report_files>0){
                   $reportfile=shift(@report_files);chomp($reportfile);
                   $ftp = Net::FTP->new("$ftp_host", Debug => 0);
                   $ftp->login($user,$pass);$ftp->cwd($ftp_dir);$ftp->binary;
                   $ftp->get("$reportfile"); #download fresh copy
                   $ftp->close;
                   ++$filecount; #print "Downloading $reportfile\n";
           }
           $command="gunzip *.gz";
           system($command);
           chdir($genome_dir); #we will run fa2dat from main genome directory
           if($filecount>0){
              $command="perl $genome_path/fa2dat.pl -i temp -o dat";
              system($command);
              system("mv temp/* mol/");
            }
           system("rm -rf temp");          
        }
     }
     chdir($organism_base); ## now create a symbolic link for current version of this organism
     $wkdir=`pwd`;
     #print "The working directory is $wkdir\n";
     if($version>0){ #assumption is that the latest version is the first on the list
        $symbfile="$organism_base/current_genome";
        if(-d "$symbfile"){system("rm $symbfile");#print "$symbfile removed\n";
         }
        system("ln -s $currentdb current_genome");  
          
      } 
     #print "The current version of $organism is $version -- $currentdb\n" if($currentdb=~/(\d+)\s*$/);
   }
 }
close(IN);
 $tm = localtime;
 my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);

  print LOG "\n************************************************************************************\n";
  print LOG "create_dir_structure.pl program complete: $mday/$mon/$yday @ $hour:$min:$sec\n";
  print LOG "\n************************************************************************************\n";
 close(LOG)


