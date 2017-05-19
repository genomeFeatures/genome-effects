#!/usr/bin/perl

use vars qw($opt_h $opt_f);
use Getopt::Std;
getopts('hf:');
open(IN,"$opt_f");open(OUT,">tx_uniq2.txt");

%map=();%expectedMap=();
while(<IN>){chomp($_);@fields=split("\t",$_);
 #@exons=split(",",$fields[2]);
 # next if(@exons<=1);
 #$fields[2]=~s/^(\d+-)//;$fields[2]=~s/(-\d+)$//;
 #@ecount=split(",",$fields[2]); next if(@ecount<=2);
 $fields[2]=~s/,$//;
 $line= join "\t",@fields;$line2= join ".",@fields;
 next if($line2=~/0,0/);
 if(!exists($map{"$line2"})){print OUT "$line\n";$map{"$line2"}=1;}}
close(IN);
print "Program complete\n";
open(IN,"mm10-GRCm38-ensGene-partialTranscripts.bed");
while(<IN>){
  chomp($_);@fields=split("\t",$_);$line2= join ".",@fields;
  $expectedMap{"$line2"}=1;
}$fount=0;
print "mm10-GRCm38-ensGene-partialTranscripts.bed has total tx=".keys(%expectedMap)."\n";
while(($key,$val)=each(%map)){
  ++$found if(exists($expectedMap{"$key"}));
}
print "Total found from graph prog:$found\n";
exit;
