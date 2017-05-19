#!/usr/bin/perl
#***********************************************************************************************
# This script loads the gene annotation data info of an
# organism,a given assembly built, the annotation name into our database
# Tables affected: transcript,exon,transcript_by_annotation,
#                  transcript_exon, transcript_translation,gene_by_annotation
#                  gene_prediction_by_organism
#
#Usage:
#   perl load_Annotations.pl -d output_directory -f filename -a annotation -v organism_version [-s 0][-e 1]

#Arguments:
#   -h  displays this help message
#   -d  Annotations directory (base directory for annotations)
#   -f  gene annotations input file
#   -a  gene prediction source name . Example: ensGene,refGene,vegaGene,...mgiGene
#   -v  Organism version (mm9,mm8, hg19,...) default mm9
#   -s  specifies whether or not the source uses zero-base or 1-base for feature start (exon, transcript,cds)
#       1 -> 1-base; 0-> zero-base (default)
#   -e  specifies whether or not the feature annotation use zero-base or 1-base for feature start (exon, transcript,cds)
#        0-> zero-base; 1 -> 1-base(default)
#  

#Examples:
#cmd: ./load_Annotations.pl -d /scratch/data/downloads/mgi -f MGI.gff3 -a mgiGene -s 1 -e 1

# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : April 2012
#
#Note:
q{I reviewed  UCSC ,MGI,Ensembl, and Ebi (miRNA) cordinates and here is what I found:
UCSC  internal database representations of coordinates always have a zero-based start and a one-based end.They add 1 to the start before displaying coordinates in the Genome Browser. Therefore, they appear as one-based start, one-based end in the graphical display.

MGI GFF3: has the start and end of the feature  in 1-based integer coordinates.
Ensembl: has the start and end of the feature  in 1-based integer coordinates
Ebi : has the start and end of the feature  in 1-based integer coordinates (assumption based on ensembl)
Uniprot: ?
DOmiceRNAseqRUM (RNASEQ) : has  a zero-based start and a one-based end

Knowing this makes it easier to integreate new datasets into our transcript database and adjust the cordinates to UCSC standard.
The reason being UCSC  has more organisms than any other source we currently load.
So this is what I intend to do:
For every new transcript source other than UCSC, set every starts and ends as needed.
This will include: transcript, exon, and CDS
};
#********************************************************************************************
## set default variables
$dbname="graber_transcriptdb";
$host="harlequin"; $user="lnh";$pass="lucie98";
$oversion="mm9";$gene_prediction="mgiGene";
$start_base=0; $end_base=1;
use DBI;use vars qw ($opt_h $opt_f $opt_d $opt_a $opt_v $opt_s $opt_e);
use Getopt::Std;use LWP 5.64;
#use Time::localtime;
use Time::localtime;
my $browser = LWP::UserAgent->new;
   $browser->timeout(10);$browser->env_proxy;
   $browser->agent('Mozilla/5.0');

getopts('hd:f:a:v:s:e:');
if(($opt_h)|| !($opt_f||$opt_d||$opt_a)){
    print <<HELP;

 This script loads a given gene annotation data
 into Graber Transcript database

Usage:
   perl load_Annotations.pl -d base_directory -f filename -a annotation -v organism_version [-s 0][-e 1]

Arguments:
   -h  displays this help message
   -d  Annotations directory (base directory for annotations)
   -f  gene annotations input file
   -a  gene prediction source name . Example: ensGene,refGene,vegaGene,...mgiGene
   -v  Organism version (mm9,mm8, hg19,...) default mm9
   -s  specifies whether or not the source uses zero-base or 1-base for feature start (exon, transcript,cds)
       1 -> 1-base; 0-> zero-base (default)
   -e  specifies whether or not the source uses zero-base or 1-base for feature end(exon, transcript,cds)
        0-> zero-base; 1 -> 1-base(default)  

Examples:
cmd: ./load_Annotations.pl -d /scratch/data/downloads/mgi -f /scratch/data/downloads/mgi/mm9-mgiGene.txt -a mgiGene -v mm9 -s 1 -e 1
cmd: ./load_Annotations.pl -d /scratch/data/downloads/rnaSeq -f /scratch/data/downloads/rnaSeq/mm9_DOmice_RUMtx_20110609.txt -a DOmiceRNAseqRUM -v mm9 -s 0 -e 1

HELP
exit;
}
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

########### addon ends #########
my $chrId="select chromosome_id from chromosome where chromosome_name=?";
my $qh_chrId=$dbh->prepare($chrId);

my $getOrg="select organism_version_id from organism_version where ucsc_db=?";   #get current organisms list 
my $qh_orgvlist= $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;
my $drop_ex_temp="drop temporary table if exists exon_temp";
#print E "$exon_id\t$chromosome\t$strand\t$exstart\t$exend\t$orgv_id\n";
my $qh_drop_ex_temp = $dbh->prepare($drop_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $create_ex_temp="create temporary table exon_temp(ex_id int unsigned default 0,
          chr_id mediumint unsigned default 0,strand char,ex_start int unsigned default 0,
          ex_end int unsigned default 0,organism_version_id smallint default 0,plus_strand char default '',minus_strand char default '',
         index(ex_id),index(organism_version_id),index(chr_id,strand,ex_start,ex_end))";
my $qh_create_ex_temp = $dbh->prepare($create_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $load_ex_temp="load data local infile ? into table exon_temp (ex_id,chr_id,strand,ex_start,ex_end,organism_version_id)";
my $qh_load_ex_temp = $dbh->prepare($load_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $update_ex_temp="update exon_temp t, exon m ";
   $update_ex_temp.="set t.ex_id=m.exon_id ";
   $update_ex_temp.=" where t.organism_version_id=m.organism_version_id and ";
   $update_ex_temp.=" t.chr_id =m.chromosome_id and t.strand=m.strand ";
   $update_ex_temp.=" and t.ex_start=m.exon_start and t.ex_end=m.exon_end "; 
my $qh_update_ex_temp = $dbh->prepare($update_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $update_ex_temp="update exon_temp t, exon m ";
   $update_ex_temp.="set t.plus_strand ='+' ";
   $update_ex_temp.=" where t.organism_version_id=m.organism_version_id and ";
   $update_ex_temp.=" t.chr_id =m.chromosome_id";
   $update_ex_temp.=" and t.ex_start=m.exon_start and t.ex_end=m.exon_end and m.strand='+' "; 
my $qh_update_ex_plustemp = $dbh->prepare($update_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $update_ex_temp="update exon_temp t, exon m ";
   $update_ex_temp.="set t.minus_strand ='-' ";
   $update_ex_temp.=" where t.organism_version_id=m.organism_version_id and ";
   $update_ex_temp.=" t.chr_id =m.chromosome_id";
   $update_ex_temp.=" and t.ex_start=m.exon_start and t.ex_end=m.exon_end and m.strand='-' "; 
my $qh_update_ex_minustemp = $dbh->prepare($update_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $get_ex_rowcount="select count(*) as rowcount from exon_temp ";
my $qh_get_ex_rowcount = $dbh->prepare($get_ex_rowcount);
my $qh_get_novelex_rowcount = $dbh->prepare("select count(*) as rowcount from exon_temp where ex_id=0");
my $qh_get_novelex_plusrowcount = $dbh->prepare("select count(*) as rowcount from exon_temp where plus_strand='+' and minus_strand=''");
my $qh_get_novelex_plusplusrowcount = $dbh->prepare("select count(*) as rowcount from exon_temp where plus_strand='+' and strand='+' and minus_strand=''");

my $qh_get_novelex_minusrowcount = $dbh->prepare("select count(*) as rowcount from exon_temp where minus_strand='-' and plus_strand=''");
my $qh_get_novelex_minusPlusrowcount = $dbh->prepare("select count(*) as rowcount from exon_temp where plus_strand='+' and minus_strand='-'");

my %orgmap=();$opt_o=$opt_d;%sample_map=();
#get the list of all organisms
my @annots=();$rows=""; $gene_prediction=$opt_a;$text_file=$opt_f;
if($opt_v){$oversion=$opt_v;}
open(LOG,">$opt_d/$oversion-$gene_prediction-log.txt");
open(DEB,">$opt_d/$oversion-$gene_prediction.issues.txt");
if(LOG){
  $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  
  $prediction_id=0;$orgv_id=0;$more=1; 
  print LOG "\n*************************************************************\n";
  print LOG "File name $opt_f: Starting load process :  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
  $qh_orgvlist->execute($oversion) or die "Can't execute query: " . $dbh->errstr . "\n";
  $update_table="$oversion-$gene_prediction";
  ($orgv_id)=$qh_orgvlist->fetchrow_array();
  if(-f "$text_file"){ $rowcount=0;$linecont=0;
         my %exonmap=(); my %transcriptmap=();#WE will generate unique local id for tx and exons
         $transcript="$opt_d/tempfile-$gene_prediction-transcript.txt";
         $exon="$opt_d/tempfile-$gene_prediction-exon.txt"; open(E,">$exon") or die "$exon :$!\n";
         open(FH,"$text_file") or die "Bad file : $!\n";
         print "Processing $text_file\n";
         if(FH){
           # $header=<FH>;chomp($header);
          $line="";
            ################# Get index of fields
            $chromindex=0;$strandindex=1;$txStartindex=2;$txEndindex=3;$count=1;
            while($line=<FH>){ chomp($line);
                ################################################################
                @content=split("\t",$line); ++$count;
                print "$count proccessed\n" if($count%10000==0);
                $name="";$strand="";$txStart=0;$txEnd=0;$cdsStart=0;$cdsEnd=0;$exonCount=0;
                $chrom=$content[$chromindex]if($chromindex>=0 && $chromindex<@content);
                $strand=$content[$strandindex] if($strandindex>=0 && $strandindex<@content);
                $exonStart=$content[$txStartindex] if($txStartindex>=0 && $txStartindex<@content);
                $exonEnd=$content[$txEndindex] if($txEndindex>=0 && $txEndindex<@content);
                $chr_id=0; $chr=$chrom;$chr=~s/^\s*//;$chr=~s/\s*$//;$chr=~s/chr//i;$strand=~s/\s*//g;$exonStart=~s/\s*//g;
                $exonEnd=~s/\s*//g;$chr=~s/\s*//g;$chr="M" if($chr=~/^M$|^MT$/i); $chr="U" if($chr=~/^U$|^UN$/i); 
                if(!exists($chrom_map{$chr})){       
                    $qh_chrId->execute($chr);
                    if($qh_chrId->rows>0){
                      ($chr_id)= $qh_chrId->fetchrow_array();$chrom_map{$chr}=$chr_id;
                    }
                 }else{$chr_id=$chrom_map{$chr};} 
                if(($chr_id<=0)|| ($exonEnd<=0)){
                     print DEB "issue bad chrom id at index $chromindex or txEnd at index $txEndindex: $line\n";next;}
                $ex_id=0;
                if(($exonStart=~/^\s*\d+\s*$/)&&($exonEnd=~/^\s*\d+\s*$/)){
                    if($opt_s ==1){$exonStart-=1;}#adjust exonStart to zero-base cordrinate
                    if($opt_e !=1){$exonEnd+=1;}  #adjust exonEnd to one-base cordrinate
                    $exonmap{$chr_id}{$strand}{$exonStart}{$exonEnd}=$ex_id;
                 } #end of @startExons==@endExons
             } #end of while(<FH>){
         } #end of if(FH){
      close(FH); close(DEB);
      #now display exons 
       if(keys(%exonmap)){
          $exon_id=0;
          for my $chromosome(sort keys(%exonmap)){
            for my $strand (sort keys (%{$exonmap{$chromosome}})){
             for my $exstart(sort keys(%{$exonmap{$chromosome}{$strand}})){
                 for my $exend(sort keys(%{$exonmap{$chromosome}{$strand}{$exstart}})){
                     $exon_id=$exonmap{$chromosome}{$strand}{$exstart}{$exend};
                     print E "$exon_id\t$chromosome\t$strand\t$exstart\t$exend\t$orgv_id\n";
                  }
             }}
          }
       }
       #close generated files
      close(E);print "$exon file generated\n";
      $linecont = `wc -l $exon`;chomp($linecont); #load exon 
      if($linecont=~/^(\d+)\s$exon/){$linecont=$1;
         if($linecont>0){
           $qh_drop_ex_temp->execute();$qh_create_ex_temp->execute();
           $qh_load_ex_temp->execute("$exon");$qh_get_ex_rowcount->execute();
           ($rowcount)=$qh_get_ex_rowcount->fetchrow_array();$novelCount=0;
           print LOG "Unique exons(chromosome,strand,exStart, exEnd): $rowcount of $linecont total were loaded -\n";
           $pluscount=0;$minuscount=0;$minusAndPlus=0;
            if($linecont==$rowcount){
               $qh_update_ex_temp->execute();$qh_get_novelex_rowcount->execute();
               ($novelCount)=$qh_get_novelex_rowcount->fetchrow_array();
              
               $qh_update_ex_plustemp->execute();$qh_get_novelex_plusrowcount->execute();
               ($pluscount)=$qh_get_novelex_plusrowcount->fetchrow_array();
               $qh_update_ex_minustemp->execute();$qh_get_novelex_minusrowcount->execute();
               ($minuscount)=$qh_get_novelex_minusrowcount->fetchrow_array();
               $qh_get_novelex_minusPlusrowcount->execute();
               ($minusAndPlus)=$qh_get_novelex_minusPlusrowcount->fetchrow_array();
            }
            print LOG "Of $rowcount, $novelCount exons were novel\n";
            $existingcount=$rowcount-$novelCount;
            print LOG "Of $existingcount existing exons:\n$pluscount exons were found only on the '+' strand\n ";
            print LOG "$minuscount exons were found only on the '-' strand\n $minusAndPlus exons were found on both strands\n ";
         }
       }
      #clean generated files
      system("rm $opt_d/tempfile-*");   
  } #end of if(-f $test_file){
  $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
 close(LOG);
} #END OF if(LOG)
print "Porgram complete\n";
exit(0);

