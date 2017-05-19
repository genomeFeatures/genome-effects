#!/usr/bin/perl

use Time::localtime;
#************************************************************
# This script calls load_annotation.pl for every annotation
# of a given organism version
#############################################################

$config=shift || "/scratch/data/ucsc/db_svn/docs/downloaded_annotations.log";
chomp($config);
$load_annotations="/scratch/data/ucsc/db_svn/src/pl/pl/load_annotations.pl";
open(CONF,"$config");
@filecontent=<CONF>;
open(LOG,">graber_tx_load-log.txt");
if(LOG){
  $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Starting load process :  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
  $count=1;
  while(@filecontent>0){
     $line=shift(@filecontent);chomp($line);$dir="";
     next if(!($line=~/mouse,|human,/i));
     my($orgg,$org,$oversion,$gene_prediction,$file)=split(",",$line);
     if($file=~/(.+)\/$oversion/){$dir=$1;}
    next if($gene_prediction=~/chromInfo|author|cell|development|ensGene|description/i);
    next if($gene_prediction=~/gbCdnaInfo|geneName|library|mrnaClone|sex|sequence|source|tissue/i);
    next if(($dir eq "")||($file eq "")); $start=0;$end=1; #default 
    $command="perl $load_annotations -d $dir -f $file -a $gene_prediction -v $oversion -s 0 -e 1";
    print LOG "loading $oversion -- $gene_prediction\n";system($command);
    print "$oversion -- $gene_prediction loaded\n";
  }
}
 $tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
close(LOG);
print "Program complete\n";

