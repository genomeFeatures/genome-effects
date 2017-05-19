#!/usr/bin/perl
#############################################################################
# this script downloads and and stores vcf files from NCBI
#ftp://ftp.ncbi.nih.gov/snp/organisms/
#####################################################
my($file)=@ARGV;
chomp($file); open(IN,"$file");
$out_file="test-mismatch.txt";open(IN,"$file");open(OUT,">$out_file");
if(!IN){ print "File eror: $! --\n"; exit;}
#get header line: chrom\tstrand\ttxStart\ttxEnd\texonCount\texonStarts\texonEnds\tName2
$header=<IN>; chomp($header); @fields=split(/\t/,$header);%headerFields=();
for my $i (0 .. $#fields){$fieldname=lc($fields[$i]);$headerFields{"$fieldname"}=$i;}
$geneIndex=$headerFields{"name2"};
print "The header line has:".keys(%headerFields)."\n";
while(<IN>){chomp($_); split /\t/;
 if(scalar(@_)!=scalar(@fields)){ #missing gene 
    print OUT join("\t",@_)."\tFieldCount:".scalar(@_)."\n";
  }
}
print "Program complete\n";
