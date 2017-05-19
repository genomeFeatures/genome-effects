#!/usr/bin/perl
#
#Project: Elizabeth.Snyder@jax.org Postdoctoral associate
# Goal: filter out lines in gtf file that do not have strands
#
use vars qw ($opt_h $opt_f); use Getopt::Std;
getopts('hf:');
if($opt_h||!($opt_f)){
  print "Script usage: ./filterGtf.pl -f annotation_filename\n";
  print "Where annotation_filename is the name of the input file (absolute path)\n";
  exit(0);
}
if(!(-f $opt_f)){print "$opt_f does not exist\n"; exit(0);}
else{
   open(W,">$opt_f-withStrand.gtf"); open(N,">$opt_f-noStrand.gtf");open(IN,"$opt_f");
   if(!W ||!N){ print "Error - could not create file : $!\n";exit(0); }
   else{ if(IN){ while(<IN>){ chomp($_);
                @fieds=split("\t",$_);$fieds[6]=~s/^\s+//; $fieds[6]=~s/\s+$//;
               if($fieds[6] eq "."){print N "$_\n";}else{print W "$_\n";}
    }}}
}
print "Program complete\n";
exit(0);

