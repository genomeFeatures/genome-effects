#!/usr/bin/env perl

use DBI;
use Time::localtime;
#************************************************************
# This script runs the pipeline
#############################################################

$scriptsBase=`pwd`; chomp($scriptsBase);$scriptsBase=~s/^\s*|\s*$|\/$//;

$config=shift || "$scriptsBase/db_svn/docs/Configuration";
chomp($config);
open(CONF,"$config");
if(!CONF){
   print"Usage: perl run_update.pl config file\n";
   exit(1);
}
@filecontent=<CONF>;

print "$scriptsBase ---\n";
#set default programs path
$getCurrentOrganisms="$scriptsBase/getCurrentOrganisms.pl";
$download_ucsc_annotations="$scriptsBase/download_ucsc_annotations.pl";
$load_geneannotations="$scriptsBase/loadAll_annotations.pl";
$load_ests="$scriptsBase/getESTsFromUCSC.pl";
$load_Static_tables="$scriptsBase/load_Static_tables.pl";
$create_dir_structure="$scriptsBase/create_dir_structure.pl";
$getgenomeList="$scriptsBase/getGenomeList.pl";
while(@filecontent>0){
   $line=shift(@filecontent);chomp($line);
   my($variable,$path)=split(",",$line);
   if($path=~/getCurrentOrganisms/){
      $getCurrentOrganisms=$path;
   }
   elsif($path=~/download_ucsc_annotations/){
      $download_ucsc_annotations=$path;
   }
   elsif($path=~/load_geneannotations/){
      $load_geneannotations=$path;
   }
   elsif($path=~/load_Static_tables/){
      $load_Static_tables=$path;
   }
   elsif($path=~/create_dir_structure/){
      $create_dir_structure=$path;
   }
   elsif($path=~/getGenomeList/){
      $getgenomeList=$path;
   }
   elsif($path=~/getESTsFromUCSC/){
      $load_est=$path;
   }
   
}
close(CONF);
print "Collecting data from the ucsc site\n";
#get current organisms info from ucsc browser
#if(-f $getCurrentOrganisms){
# print "Processing $getCurrentOrganisms\n";
# system("perl $getCurrentOrganisms $config"); #done
#}
#else{print "$getCurrentOrganisms does not exist -"; exit(0);}

#print"Downloading data from ucsc database\n";
# now download annotations
#system("perl $download_ucsc_annotations");

#load new ids
print "loading static tables\n";
system("perl $load_Static_tables");

#print "running load_geneannotations\n";
#load the annotation int our db
#system("perl $load_geneannotations");

#print "running load_est\n";
#system("perl $load_ests");
#now create the expected directory structure
# and download the genome data
##print "Creating directory structure\n";
##system("perl $create_dir_structure");
#now generate the master lists
##print "Generating list of organisms with downloaded genome data\n";
##system("perl $getgenomeList");
print "Program complete\n";






