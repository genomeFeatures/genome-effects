#!/usr/bin/perl

my($file_name)=@ARGV;
chomp($file_name);$out="$file_name-out.txt";
open(OUT,">$out");
print OUT "CDSLen\tTranscriptCount\n";
open(IN,"$file_name") or die "$!\n";
$header=<IN>; %cds=();
while(<IN>){ chomp($_);
 ($Name,$chrom,$strand,$txStart,$txEnd,$cdsStarts,$cdsEnds,$exonCount,$exonStarts,$exonEnds,$Name2,$sources)=split("\t",$_);
 @cdsstarts=split(",",$cdsStarts);@cdsends=split(",",$cdsEnds);
 for($i=0;$i<@cdsstarts;++$i){ $len=$cdsends[$i]-$cdsstarts[$i];
    if(!exists($cds{"$len"})){$cds{"$len"}=1;}else{$cds{"$len"}+=1;}} 
}
for my $len(sort keys(%cds)){
  print OUT "$len\t".$cds{"$len"}."\n";
}
print "Program complete\n";
