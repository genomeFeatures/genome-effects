#!/usr/bin/perl

use DBI;
use Time::localtime;
#************************************************************
# This script runs the pipeline
#############################################################

$config=shift || "/scratch/data/ucsc/ucsc_load_process/configuration/Configuration";
chomp($config);
open(CONF,"$config");
if(!CONF){
   print"Usage: perl run_update_pipeline.pl config file\n";
   exit(1);
}
@filecontent=<CONF>;
#set default programs path
$getCurrentOrganisms="/scratch/data/ucsc/ucsc_load_process/src/pl/getCurrentOrganisms.pl";
$download_ucsc_annotations="/scratch/data/ucsc/ucsc_load_process/src/pl/download_ucsc_annotations.pl";

$generate_annotationsformat="/scratch/data/ucsc/ucsc_load_process/src/pl/generate_annotationsformat.pl";
$load_geneannotations="/scratch/data/ucsc/ucsc_load_process/src/pl/load_geneannotations.pl";
$load_ests="/scratch/data/ucsc/ucsc_load_process/src/pl/getESTsFromUCSC.pl";
$load_Static_tables="/scratch/data/ucsc/ucsc_load_process/src/pl/load_Static_tables.pl";
$create_dir_structure="/scratch/data/ucsc/ucsc_load_process/src/pl/create_dir_structure.pl";
$getgenomeList="/scratch/data/ucsc/ucsc_load_process/src/pl/getGenomeList.pl";
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
system("perl $getCurrentOrganisms $config");
print"Downloading data from ucsc database\n";
# now download annotations
system("perl $download_ucsc_annotations");
#load new ids
print "loading static tables\n";
system("perl $load_Static_tables ");
print "running load_geneannotations\n";
#load the annotation int our db
system("perl $load_geneannotations");
#print "running load_est\n";
#system("perl $load_ests");
#now create the expected directory structure
# and download the genome data
print "Creating directory structure\n";
system("perl $create_dir_structure");
#now generate the master lists
print "Generating list of organisms with downloaded genome data\n";
system("perl $getgenomeList");
print "Program complete\n";






