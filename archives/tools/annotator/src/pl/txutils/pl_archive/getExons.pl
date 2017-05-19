#!/usr/bin/perl
####
# This script generate exons and exon junctions files from 
# from a genePred file with annotations
# The data is to test the getConnected components graph program
#
###
use vars qw ($opt_h $opt_f);
use Getopt::Std;
getopts('hf:');
if($opt_h||(!$opt_f)) {
print <<HELP; 
    Usage: ./getExons.pl -f annotations_file_name 
HELP
exit;
}
#
# this function returns two lists: exonJunctions list and internal exons list
#
sub getJunctions {
   ($junction_ref,$exon_ref,$v1_end_ref,$v2_start_ref,$exonStarts,$exonEnds)=@_;
   if($exonStarts=~/,/){ #$starts=".,";$ends="";
      @estarts=split(",",$exonStarts);@ends=split(",",$exonEnds);
      if(@estarts>1){ #there is more than one exon
         ${$v1_end_ref}=$ends[0];${$v2_start_ref}=$estarts[@ends-1];
         for my $i(1 .. $#estarts){
             $v1_end=$ends[$i-1]; $v2_start=$estarts[$i];$v2_end=$ends[$i];
             ${$junction_ref}{"$v1_end"}=$v2_start;  #$ends.="$v1_end,";
             ${$exon_ref}{"$v2_start"}=$v2_end if($i<(@ends-1));
         }
      }
   }
}
open(IN,"$opt_f");$header=<IN>;%fieldsIndex=();chomp($header);
#name	chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	name2	source
@fields=split(/\t/,$header);
for my $i(0..$#fields){$field=$fields[$i];$fieldsIndex{"$field"}=$i;}
$chromindex=$fieldsIndex{"chrom"};$strandindex=$fieldsIndex{"strand"};
$estartsIndex=$fieldsIndex{"exonStarts"};$endsIndex=$fieldsIndex{"exonEnds"};
open(JUNC,">mm10-GRCm38-ensGene-exonJunctions.bed");open(EX,">mm10-GRCm38-ensGene-InternalExons.bed");
open(TX,">mm10-GRCm38-ensGene-partialTranscripts.bed");
print JUNC "track	name=mm10 ensembl junctions	visibility=\n";
print EX "track	name='mm10 ensembl Internal Exons'  description=\n";
%uniqJunc=();%uniqEx=();
while(<IN>){ chomp($_);@fields=split("\t",$_);
   $exonStarts=$fields[$estartsIndex];$exonsEnds=$fields[$endsIndex];
   %junctions=();%exons=();$firstv1_end=0;$lsatv2_start=0; $tx="";$chrom=$fields[$chromindex];
   $strand=$fields[$strandindex];
   getJunctions(\%junctions,\%exons,\$firstv1_end,\$lastv2_start,$exonStarts,$exonsEnds);
   next if(keys(%junctions)<=0); 
   while(($v1_end,$v2_start)=each(%junctions)){$v1_end-=50,$v2_start+=50;
    next if($uniqJunc{"$chrom:$v1_end:$v2_start:$strand"});
    print JUNC "$chrom\t$v1_end\t$v2_start\t.\t.\t$strand\n";$uniqJunc{"$chrom:$v1_end:$v2_start:$strand"}=1;
   }
   for my $v_start(sort keys(%exons)){$v_end=$exons{"$v_start"};
        next if($uniqEx{"$chrom:$v_start:$v_end"});
       print EX "$chrom\t$v_start\t$v_end\n";$tx.="$v_start-$v_end,";$uniqEx{"$chrom:$v_start:$v_end"}=1;
   }
   print TX "$chrom\t$strand\t$firstv1_end,$tx$lastv2_start\n";
}
print "Program complete -- The transcript is: $tx\n";
exit;
