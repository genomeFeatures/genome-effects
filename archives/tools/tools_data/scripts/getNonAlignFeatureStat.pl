#!/usr/bin/perl

#######################################################################
## getNonAlignFeatureStat.pl
#  Generates statistics on non-aligned portions of
#  ESTs that have some kind of alignment on the genome.

#  These include:
#  1. NON-aligned sequence length
#  2. non-aligned sequence,3.non-aligned sequence type (5'end,3'end)
#  4. non-aligned sequence nucleic acid count (T:2,C:2,A:3,G:4)
#  5. runs of Ts length 
#  6. runs of As length

#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : January, 2012
#   Usage: ./getNonAlignFeatureStat.pl  -h db_host_name -d db_name
#
# seq_type=1 -> 5'end, seq_type=2 -> 3'end
# db table desing: organism,accession,qStart,qEnd,tSart,
#   tEnd,sequence,seq_type,T_count,C_count,G_count,A_count,t_strand
#   runOfTs,runOfAs

#######################################################################


use DBI;
use vars qw ($opt_h $opt_d $opt_g $opt_f $opt_o $opt_v);
use Getopt::Std;

getopts('hg:d:f:o:v:');
if($opt_h) {
    print <<HELP;

 This script generates statistics of non-aligned ends of
 every EST in graber transcript database that has an alignment.
 EST  alignments were downloaded from UCSC and are updated
 on a monthly basis.If the user specified the alignment file,
 that file will be used instead

Usage:

   perl getNonAlignFeatureStat.pl [-g <db host name>][-d <dbName>][-f <alignment file>][-o <output_filename>]
Arguments:
   -h  displays this help message
   -g  The name of the database server (harlequin, or demon) -> optional, default is demon
   -d  The the name of the database-> optional,default is graber_transcriptdb
   -f  The name of a tab-delimitted file containing the alignment data -> optional
       default alignments are generated from a table in the database
       The file should have at least the following field/labels
       qName,qStart,and qEnd
   -o  Output file name (optional) default stdout
   -v Organism version (mm9,mm8, hg18,...) optional, default all 

Example1: perl getNonAlignFeatureStat.pl -h demon.jax.org -d graber_transcriptdb
Example2: perl getNonAlignFeatureStat.pl -d graber_transcriptdb
Example3: perl getNonAlignFeatureStat.pl -h demon.jax.org 
Example4: perl getNonAlignFeatureStat.pl -v hg19 
Example5: perl getNonAlignFeatureStat.pl -f myESTsAlignmentFile.txt -o result.txt
.....
HELP
exit;
}
my $dbname ='graber_transcriptdb';
my $host ='demon';my $user ='lnh';my $pwd  ='lucie';

my($inputFile,$outputFile,$organism);
$inputFile="";$outputFile="";$organism="";
if($opt_g){$host=$opt_g;}if($opt_d){$dbname=$opt_d;} 
if($opt_f){$inputFile=$opt_f;}if($opt_o){$outputFile=$opt_o;}
if($opt_v){$organism=$opt_v;}

my $getChr ="select c.chromosome_id,c.chromosome_name from chromosome c,chromosome_data cd where organism_version_id=? ";
    $getChr.=" and c.chromosome_id=cd.chromosome_id";

my $getOrganisms="select o.organism_name,o.organism_id,ov.organism_version_id,ov.ucsc_db ";
   $getOrganisms.=" from organism o, organism_version ov where o.organism_id=ov.organism_id";

my $getAlignment="select qName,qSize,strand,qStart,qEnd,tStart,tEnd";
   $getAlignment.=" from est_align where organism_version_id=? and chromosome_id=? ";
my $getSequence="select qName,sequence from sequence where organism_id=? and qName=?";

my $dbh = DBI->connect("DBI:mysql:$dbname:$host",$user,$pwd)or die "Can't connect to database!\n";

my $getChr_h= $dbh->prepare($getChr);
my $getOrganisms_h= $dbh->prepare($getOrganisms);
my $getAlignment_h= $dbh->prepare($getAlignment);
my $getSequence_h= $dbh->prepare($getSequence);

my %organisms=(); 
# get organism map
if(!$opt_f){ #for each organism, get
   $getOrganisms_h->execute();
   if($getOrganisms_h->rows>0){
       while(($organism_name,$organism_id,$organism_version_id,$ucsc_db)=$getOrganisms_h->fetchrow_array()){
            if($ucsc_db eq "mm9"){
             open(OUT,">$ucsc_db-ESTs.stat");  # create output file
             if(OUT){
                print OUT "organism_version_id\tAccession\tqStart\tqEnd\ttSart\ttEnd\ttStrand\tqSize\tactual-qSize\thas_same_qSize\t";
                print OUT "five_sequence\tfiv_seq_len\tfive_T_count\tfive_A_count\tfive_maxRunOfTs\tfive_maxRunOfAs\t";
                print OUT "three_sequence\tthree_seq_len\tthree_T_count\tthree_A_count\tthree_maxRunOfTs\tthree_maxRunOfAs\n";
                $getChr_h->execute($organism_version_id);
                if($getChr_h->rows>0){ my %estmap=();
                   while(($chromosome_id,$chromosome_name)=$getChr_h->fetchrow_array()){
                      print "Processing $ucsc_db-$chromosome_name\n";
                      #now get the alignment data for this chromosome
                      $getAlignment_h->execute($organism_version_id,$chromosome_id);
                      if($getAlignment_h->rows>0){
                         while(($qname,$qSize,$strand,$qStart,$qEnd,$tStart,$tEnd)=$getAlignment_h->fetchrow_array()){
                              $line_f="$organism_version_id\t$qname\t$qStart\t$qEnd\t$tStart\t$tEnd\t$strand\t";
                              $line_f.="$qSize\t";
                              next if(($qStart<=5)&&(($qSize-$qEnd)>=5));
                              if(!exists($estmap{"$qname"})){
                                  $getSequence_h->execute($organism_id,$qname); #get this EST sequence
                                  if($getSequence_h->rows>0){
                                     ($qname,$sequence)=$getSequence_h->fetchrow_array();
                                     $estmap{"$qname"}=$sequence;
                                   }
                              }
                              #print "-- has sequence -- $line_f\n";
                              if(exists($estmap{"$qname"})){ 
                                $sequence=$estmap{"$qname"};
                                $actual_qSize=length($sequence);
                                $has_sameQsize=1;  $has_sameQsize=0 if($qSize!=$actual_qSize);
                                $fivePrime_end=substr($sequence,0,$qStart);
                                $threePrime_end=substr($sequence,$qEnd, $actual_qSize-$qEnd);
                                $fivePrime_seq_len=length($fivePrime_end);
                                $threePrime_seq_len=length($threePrime_end);
                                $five_T_count=0;$five_A_count=0;$five_maxRunsOfTs=0;$five_maxRunsOfAs=0;
                                $three_T_count=0;$three_A_count=0;$three_maxRunsOfTs=0;$three_maxRunsOfAs=0;
                                while ($fivePrime_end =~/(T)/ig) {$five_T_count+=1;}
                                while ($fivePrime_end =~/(A)/ig) {$five_A_count+=1;}
                                while ($threePrime_end =~/(T)/ig) {$three_T_count+=1;}
                                while ($threePrime_end =~/(A)/ig) {$three_A_count+=1;}
                                # get max runs of A/T
                               @matches=($threePrime_end =~/((A)+)/ig);
                               foreach $match(@matches){
                                      $three_maxRunsOfAs=length($match) if($three_maxRunsOfAs<length($match));
                               }
                               @matches=();@matches=($threePrime_end =~/((T)+)/ig);
                               foreach $match(@matches){
                                  $three_maxRunsOfTs=length($match) if($three_maxRunsOfTs<length($match));
                                }
                               @matches=();@matches=($fivePrime_end =~/((A)+)/ig);
                               foreach $match(@matches){
                                    $five_maxRunsOfAs=length($match) if($five_maxRunsOfAs<length($match));
                                 }
                               @matches=();@matches=($fivePrime_end =~/((T)+)/ig);
                               foreach $match(@matches){
                                    $five_maxRunsOfTs=length($match) if($five_maxRunsOfTs<length($match));
                                }
                                $line_f.="$actual_qSize\t$has_sameQsize\t";
                                $line_f.="$fivePrime_end\t$fivePrime_seq_len\t$five_T_count\t$five_A_count\t";
                                $line_f.="$five_maxRunsOfTs\t$five_maxRunsOfAs\t";
                                $line_f.="$threePrime_end\t$threePrime_seq_len\t$three_T_count\t$three_A_count\t";
                                $line_f.="$three_maxRunsOfTs\t$three_maxRunsOfAs\n"; 
                                print OUT "$line_f";  
                            }
                           # else{ #missing est sequence
                           #    print "$qName on $ucsc_db chromosome $chromosome_name is missing the sequence\n";
                          #  }
                         }# end of while getAlignment
                      } #end of getAlignment rows>0
                   }#end of chromosome list
                 undef  %estmap;
               }# end of getChr rows>0
             }# end of if(OUT)
             last;
            }
          }#END of organism list
      }# end of getOrganisms rows>0
} #end of opt_f check

$dbh->finish;
$dbh->disconnect;
print"program completed\n";
            
