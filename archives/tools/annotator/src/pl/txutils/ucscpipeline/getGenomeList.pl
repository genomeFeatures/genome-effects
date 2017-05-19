#!/usr/bin/perl

use Net::FTP;
use File::stat;
use Time::localtime;

#*****************************************************************
# This script  generates the list of organisms that have 
# genome data
#

# This script is ran after create_dir_structure.pl
# 
# Input : a config file (current_organism_version) 
#         with the following fields
#       1.organism_group,
#       2.organims_name, 
#       3.ucsc_db, 
#       4.version
# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : January 2011
# Usage: perl program configfile 
#**************************************************************
#
$config=shift||"/scratch/data/ucsc/ucsc_load_process/configuration/Configuration";
chomp($config);open(CONF,"$config");
if(!CONF){
   print"Config file $config missing\n";
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
# 
# hostname, user, pass and remote dir
#

$genome_path= $config{"genome_path"};
$load_base_path=$config{"ucsc_load_base"};
my $gorganisms_w_genome_data="$genome_path/organisms_w_chromosome_data.txt";
my $galaxy_file="$genome_path/galaxy_genomeList.txt";
my $configpath="$load_base_path/configuration";
my $mainConfig ="$configpath/current_organism_version";
my $otherConfig ="$configpath/other_organism_version.txt";
my $processlog="$configpath/data_load_log.log";

open(LOG,">>$processlog");
if(LOG){
  $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Starting getGenomeList.pl program :  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
  chomp($mainConfig);
  if(!$mainConfig){
      print LOG "getGenomeList.pl usage was incorect, missing the config file argument\n";close(LOG); exit(1);
   }
  open(IN,"$mainConfig");
  open(OT,"$otherConfig");
  if(!IN){
      print  LOG "Bad config file argument\n"; close(IN);exit(1);
   } 
  @others=();
  if(OT){@others=<OT>; close(OT);}
  open(OUT,">$gorganisms_w_genome_data");
  open(OUT2,">$galaxy_file");
  print OUT2 "#Harvested from getGenomeList.pl\n";
  print OUT "Organism group,organsim name,version,path2genomedata\n";
  #print  "Bad config file argument\n"; exit(0);
   @organisms=<IN>;
   foreach my $line(@organisms){
         if($line){
            chomp($line);
            my($organism_group,$organism,$db,$version)=split(",",$line);
            chomp($organism);
            $organism=lc($organism);$organismDir=lc($organism);
            $genome_dir="$genome_path/$organism_group/$organismDir/$db";
            $organism_base="$genome_path/$organism_group/$organismDir";
            $organismDir=~s/\s+/-/g; 
            $genome_data_dir="$genome_dir/dat"; $genome_temp="$genome_dir/temp";
            #NOW get list of files under dat/ directory
           if( -d $genome_temp){
                 system("rm -rf $genome_temp");
            }
           if(-d $genome_data_dir){
               opendir MYDIR,"$genome_data_dir";
               @dirContent = grep /\.dat$/,readdir MYDIR; # collect all fasta files
               close(MYDIR);                               # try to update the first file
               if(@dirContent>0){
                  print OUT "$organism_group,$organism,$db,$genome_data_dir\n" if(OUT);
                  print OUT2 "$db\t$organism - $db\n";
               }  
            }              
        }
     }
    foreach my $line(@others){
         if($line){
            chomp($line);
            my($organism_group,$organism,$db,$version)=split(",",$line);
            chomp($organism);
            $organism=lc($organism);$organismDir=lc($organism);
            $genome_dir="$genome_path/$organism_group/$organismDir/$db";
            $organism_base="$genome_path/$organism_group/$organismDir";
            $organismDir=~s/\s+/-/g; 
            $genome_data_dir="$genome_dir/dat"; $genome_temp="$genome_dir/temp";
            #NOW get list of files under dat/ directory
           if( -d $genome_temp){
                 system("rm -rf $genome_temp");
            }
           if(-d $genome_data_dir){
               opendir MYDIR,"$genome_data_dir";
               @dirContent = grep /\.dat$/,readdir MYDIR; # collect all fasta files
               close(MYDIR);                               # try to update the first file
               if(@dirContent>0){
                  print OUT "$organism_group,$organism,$db,$genome_data_dir\n" if(OUT);
                  print OUT2 "$db\t$organism - $db\n";
               }  
            }              
        }
     }
  close(OUT);close(OUT2);
 }
close(IN);
print "Program complete\n";
$tm = localtime;
 my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);

  print LOG "\n************************************************************************************\n";
  print LOG "getGenomeList.pl program complete: $mday/$mon/$yday @ $hour:$min:$sec\n";
  print LOG "\n************************************************************************************\n";
 close(LOG)


