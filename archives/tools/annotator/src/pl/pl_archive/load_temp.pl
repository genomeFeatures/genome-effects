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
#print T "$transcript_id\t$chromosome\t$strand\t$txstart\t$txend\t$orgv_id\n"
my $drop_tx_temp="drop temporary table if exists transcript_temp";
my $qh_drop_tx_temp = $dbh->prepare($drop_tx_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
# print T "$transcript_id\t$chromosome\t$strand\t$txstart\t$txend\t$orgv_id\n";
my $create_tx_temp="create temporary table transcript_temp(tx_id int unsigned default 0,
          chr_id mediumint unsigned default 0,strand char(1),tx_start int unsigned default 0,
          tx_end int unsigned default 0,organism_version_id smallint default 0,
         index(tx_id),index(organism_version_id),index(chr_id,strand,tx_start,tx_end))";
my $qh_create_tx_temp = $dbh->prepare($create_tx_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $load_tx_temp="load data local infile ? into table transcript_temp";
my $qh_load_tx_temp = $dbh->prepare($load_tx_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $update_tx_temp="update transcript_temp t, transcript m ";
   $update_tx_temp.="set t.tx_id=m.transcript_id ";
   $update_tx_temp.=" where t.organism_version_id=m.organism_version_id and ";
   $update_tx_temp.=" t.chr_id =m.chromosome_id and t.strand=m.strand ";
   $update_tx_temp.=" and t.tx_start=m.tx_start and t.tx_end=m.tx_end "; 
my $qh_update_tx_temp = $dbh->prepare($update_tx_temp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $get_tx_rowcount="select count(*) as rowcount from transcript_temp ";
my $qh_get_tx_rowcount = $dbh->prepare($get_tx_rowcount);
my $qh_get_noveltx_rowcount = $dbh->prepare("select count(*) as rowcount from transcript_temp where tx_id=0");

my $drop_ex_temp="drop temporary table if exists exon_temp";
#print E "$exon_id\t$chromosome\t$strand\t$exstart\t$exend\t$orgv_id\n";
my $qh_drop_ex_temp = $dbh->prepare($drop_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $create_ex_temp="create temporary table exon_temp(ex_id int unsigned default 0,
          chr_id mediumint unsigned default 0,strand char(1),ex_start int unsigned default 0,
          ex_end int unsigned default 0,organism_version_id smallint default 0,
         index(ex_id),index(organism_version_id),index(chr_id,strand,ex_start,ex_end))";
my $qh_create_ex_temp = $dbh->prepare($create_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $load_ex_temp="load data local infile ? into table exon_temp";
my $qh_load_ex_temp = $dbh->prepare($load_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $update_ex_temp="update exon_temp t, exon m ";
   $update_ex_temp.="set t.ex_id=m.exon_id ";
   $update_ex_temp.=" where t.organism_version_id=m.organism_version_id and ";
   $update_ex_temp.=" t.chr_id =m.chromosome_id and t.strand=m.strand ";
   $update_ex_temp.=" and t.ex_start=m.exon_start and t.ex_end=m.exon_end "; 
my $qh_update_ex_temp = $dbh->prepare($update_ex_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $get_ex_rowcount="select count(*) as rowcount from exon_temp ";
my $qh_get_ex_rowcount = $dbh->prepare($get_ex_rowcount);
my $qh_get_novelex_rowcount = $dbh->prepare("select count(*) as rowcount from exon_temp where ex_id=0");

my %orgmap=();$opt_o=$opt_d;%sample_map=();%sample_dates=*);
#get the list of all organisms
my @annots=();$rows=""; $gene_prediction=$opt_a;$text_file=$opt_f;
if($opt_v){$oversion=$opt_v;}
open(LOG,">>$opt_d/$oversion-$gene_prediction-log.txt");
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
         $exon="$opt_d/tempfile-$gene_prediction-exon.txt"; 
         open(T,">$transcript")or die "$transcript :$!\n";open(E,">$exon") or die "$exon :$!\n";
         open(FH,"$text_file") or die "Bad file : $!\n"; $lcount=0;
         print "Processing $text_file\n";
         if(FH){
            $header=<FH>;chomp($header);$line="";@fields=split("\t",$header);
           #Sample_name 	 chrom 	 strand 	 txStart 	 txEnd 	 exonCount 	 exonStarts 	 exonEnds 	 Sample_date
           print "The header line has ".@fields." fields\n";
            if(@fields>0){$index=0;
                while(@fields>0){$field=shift(@fields);
                      chomp($field);$field=~s/\s*$//;$field=~s/^\s*//;$line.=",$field>$index"; ++$index;}
            }
            print "The fields line is $line\n";
            #exit(0);
            ################# Get index of fields
            $nameindex=-1;$chromindex=-1;$strandindex=-1;$txStartindex=-1;$exonSizesindex=-1;
            $txEndindex=-1;$cdsStartindex=-1;$cdsEndindex=-1;$exonCountindex=-1;$exonStartsindex=-1; 
            $exonEndsindex=-1;$name2index=-1;$exonFramesindex=-1;$matchesindex=-1;$misMatchesindex=-1;
            $qNameindex=-1;$qSizeindex=-1;$qStartindex=-1;$qEndindex=-1;$sizePolyAindex=-1;
            $revSizePolyAindex=-1;$signalPosindex=-1;$revSignalPosindex=-1;$chromStartindex=-1;
            $sample_nameindex=-1; $sample_dateindex=-1;
            $chromEndindex=-1;$proteinIDindex=-1;$has_block_size=0;%chrom_map=();
            if($line=~/,chrom>(\d+)/i){$chromindex=$1;}if($line=~/,name>(\d+)/i){$nameindex=$1;}
            if($line=~/,strand>(\d+)/i){$strandindex=$1;}if($line=~/,txStart>(\d+)/i){$txStartindex=$1;}
            if($line=~/,txEnd>(\d+)/i){$txEndindex=$1;}if($line=~/,cdsStart>(\d+)/i){$cdsStartindex=$1;}
            if($line=~/,cdsEnd>(\d+)/i){$cdsEndindex=$1;}if($line=~/,exonCount>(\d+)/i){$exonCountindex=$1;}
            if($line=~/,exonStarts>(\d+)/i){$exonStartsindex=$1;}if($line=~/,exonEnds>(\d+)/i){$exonEndsindex=$1;}
            if($line=~/,name2>(\d+)/i){$name2index=$1;}if($line=~/,exonFrames>(\d+)/i){$exonFramesindex=$1;}
            if($line=~/,chromStart>(\d+)/i){$txStartindex=$1;}if($line=~/,chromEnd>(\d+)/i){$txEndindex=$1;}
            if($line=~/,proteinID>(\d+)/i){$proteinIDindex=$1;}if($line=~/,matches>(\d+)/i){$matchesindex=$1;}
            if($line=~/,misMatches>(\d+)/i){$misMatchesindex=$1;}if($line=~/,tName>(\d+)/i){$chromindex=$1;}
            if($line=~/,qName>(\d+)/i){$nameindex=$1;}if($line=~/,tStart>(\d+)/i){$txStartindex=$1;}
            if($line=~/,tEnd>(\d+)/i){$txEndindex=$1;}if($line=~/,qStart>(\d+)/i){$qStartindex=$1;}
            if($line=~/,qEnd>(\d+)/i){$qEndindex=$1;}if($line=~/,qSize>(\d+)/i){$qSizeindex=$1;}
            if($line=~/,sizePolyA>(\d+)/i){$sizePolyAindex=$1;}if($line=~/,revSizePolyA>(\d+)/i){$revSizePolyAindex=$1;}
            if($line=~/,signalPos>(\d+)/i){$signalPosindex=$1;}if($line=~/,revSignalPos>(\d+)/i){$revSignalPosindex=$1;}
            if($line=~/,blockSizes>(\d+)/i){$exonSizesindex=$1; $has_block_size=1;}
            $count=1;
            while($line=<FH>){ chomp($line);
                ################################################################
                @content=split("\t",$line); ++$count;
                 print "$count proccessed\n" if($count%100000==0);
                $name="";$strand="";$txStart=0;$txEnd=0;$cdsStart=0;$cdsEnd=0;$exonCount=0;
                $exonStarts="";$exonEnds="";$name2="";$exonFrames="";$proteinID="";$score=0;
                $exonSizes="";$matches=-1;$misMatches=-1;$qStart=-1;$qEnd=-1;$qSize=-1;$sizePolyA=-1;
                $revSizePolyA=-1;$signalPos=-1;$revSignalPos=-1;$trans_id=0;$ex_id=0;
                $chrom=$content[$chromindex]if($chromindex>=0 && $chromindex<@content);
                $name=$content[$nameindex] if($nameindex>=0 && $nameindex<@content);
                $strand=$content[$strandindex] if($strandindex>=0 && $strandindex<@content);
                $txStart=$content[$txStartindex] if($txStartindex>=0 && $txStartindex<@content);
                $txEnd=$content[$txEndindex] if($txEndindex>=0 && $txEndindex<@content);
                $cdsStart=$content[$cdsStartindex] if($cdsStartindex>=0 && $cdsStartindex<@content);
                $cdsEnd=$content[$cdsEndindex] if($cdsEndindex>=0 && $cdsEndindex<@content);
                $exonCount=$content[$exonCountindex] if($exonCountindex>=0 && $exonCountindex<@content);
                $exonStarts=$content[$exonStartsindex] if($exonStartsindex>=0 && $exonStartsindex<@content);
                $exonEnds=$content[$exonEndsindex] if($exonEndsindex>=0 && $exonEndsindex<@content);
                $name2=$content[$name2index] if($name2index>=0 && $name2index<@content);
                $exonFrames=$content[$exonFramesindex] if($exonFramesindex>=0 && $exonFramesindex<@content);
                $proteinID=$content[$proteinIDindex] if($proteinIDindex>=0 && $proteinIDindex<@content);
                $matches=$content[$matchesindex] if($matchesindex>=0 && $matchesindex<@content);
                $misMatches=$content[$misMatchesindex] if($misMatchesindex>=0 && $misMatchesindex<@content);
                $qStart=$content[$qStartindex] if($qStartindex>=0 && $qStartindex<@content);
                $qEnd=$content[$qEndindex] if($qEndindex>=0 && $qEndindex<@content);
                $qSize=$content[$qSizeindex] if($qSizeindex>=0 && $qSizeindex<@content);
                $sizePolyA=$content[$sizePolyAindex] if($sizePolyAindex>=0 && $sizePolyAindex<@content);
                $revSizePolyA=$content[$revSizePolyAindex] if($revSizePolyAindex>=0 && $revSizePolyAindex<@content);
                $signalPos=$content[$signalPosindex] if($signalPosindex>=0 && $signalPosindex<@content);
                $revSignalPos=$content[$revSignalPosindex] if($revSignalPosindex>=0 && $revSignalPosindex<@content);
                $exonSizes=$content[$exonSizesindex] if($exonSizesindex>=0 && $exonSizesindex<@content);
                $chr_id=0; $chr=$chrom;$chr=~s/^\s*//;$chr=~s/\s*$//;$chr=~s/chr//i;$strand=~s/\s*//g;$txStart=~s/\s*//g;
                $txEnd=~s/\s*//g;$chr=~s/\s*//g;
                $chr="M" if($chr=~/^M$|^MT$/i); $chr="U" if($chr=~/^U$|^UN$/i); 
                if(!exists($chrom_map{$chr})){       
                    $qh_chrId->execute($chr);
                    if($qh_chrId->rows>0){
                      ($chr_id)= $qh_chrId->fetchrow_array();$chrom_map{$chr}=$chr_id;
                    }
                 }
                else{$chr_id=$chrom_map{$chr};} 
                if(($chr_id<=0)|| ($txEnd<=0)){
                     print DEB "issue bad chrom id at index $chromindex or txEnd at index $txEndindex: $line\n";next;}
                if(($strand ne "+")&& ($strand ne "-")){print DEB "issue bad strand at index $strandindex: $line\n";next;}
                if($opt_s ==1){$txStart-=1;}#adjust txStart to zero-base cordrinate
                if($opt_e !=1){$txEnd+=1;}#adjust txEnd to one-base cordrinate
                #display this transcript
                @startExons=split(",",$exonStarts);$tx_id=$name;
                @endExons=split(",",$exonEnds);$exon_frame=-1;$exist=0;
                @exonFrames=split(",",$exonFrames); $trans_id=0;$ex_id=0;
                if(@startExons==@endExons){
                    while(@startExons>0){
                        $exonStart=shift(@startExons);$exon_frame=-1;
                        if(@endExons>0){$exonEnd=shift(@endExons);$exon_frame=shift(@exonFrames);}
                        if(($exonStart=~/^\s*\d+\s*$/)&&($exonEnd=~/^\s*\d+\s*$/)){
                          if($opt_s ==1){$exonStart-=1;}#adjust exonStart to zero-base cordrinate
                          if($opt_e !=1){$exonEnd+=1;}  #adjust exonEnd to one-base cordrinate
                          $exonmap{$chr_id}{$strand}{$exonStart}{$exonEnd}=$ex_id;
                        }
                     } #end of @startExons>0
                 } #end of @startExons==@endExons
                 ##################################################
                 $transcriptmap{$chr_id}{$strand}{$txStart}{$txEnd}=$trans_id;
             } #end of while(<FH>){
         } #end of if(FH){
      close(FH); close(DEB);
     # print "Generating tx keys\n";  
     #now generate and display transcripts
     if(keys(%transcriptmap)>0){
        print "Generating tx keys\n"; 
        $transcript_id=0;$exists=0; $exon_id=0;
        for my $chromosome(sort keys(%transcriptmap)){
           for my $strand (sort keys (%{$transcriptmap{$chromosome}})){
                for my $txstart(sort keys(%{$transcriptmap{$chromosome}{$strand}})){
                     for my $txend(sort keys(%{$transcriptmap{$chromosome}{$strand}{$txstart}})){
                           print T "$transcript_id\t$chromosome\t$strand\t$txstart\t$txend\t$orgv_id\n"; 
                     }
                }
            }
         }
       }
       #now display exons 
       if(keys(%exonmap)){
          $exon_id=0;
         for my $chromosome(sort keys(%exonmap)){
             for my $strand (sort keys (%{$exonmap{$chromosome}})){
                 for my $exstart(sort keys(%{$exonmap{$chromosome}{$strand}})){
                      for my $exend(sort keys(%{$exonmap{$chromosome}{$strand}{$exstart}})){
                           $exon_id=$exonmap{$chromosome}{$strand}{$exstart}{$exend};
                           if(($exend>0)&&($exon_id<=0)){
                               print E "$exon_id\t$chromosome\t$strand\t$exstart\t$exend\t$orgv_id\n";
                            }
                       }
                  }
              }
         }
       }
       #close generated files
      close(T);close(E);
       $linecont = `wc -l $transcript`;chomp($linecont); #load transcript 
      if($linecont=~/^(\d+)\s$transcript/){$linecont=$1;
         if($linecont>0){
           $qh_drop_tx_temp->execute();$qh_create_tx_temp->execute();
           $qh_load_tx_temp->execute("$transcript");$qh_get_tx_rowcount->execute();
           ($rowcount)=$qh_get_tx_rowcount->fetchrow_array();$novelCount=0;
           if($linecont==$rowcount){
             $qh_update_tx_temp->execute();$qh_get_noveltx_rowcount->execute();
               ($novelCount)=$qh_get_noveltx_rowcount->fetchrow_array();
           } 
         print LOG "Unique transcripts: $rowcount of $linecont total were loaded - $novelCount transcripts were novel\n";
         }
      }
      $linecont = `wc -l $exon`;chomp($linecont); #load exon 
      if($linecont=~/^(\d+)\s$exon/){$linecont=$1;
         if($linecont>0){
           $qh_drop_ex_temp->execute();$qh_create_ex_temp->execute();
           $qh_load_ex_temp->execute("$exon");$qh_get_ex_rowcount->execute();
           ($rowcount)=$qh_get_ex_rowcount->fetchrow_array();$novelCount=0;
            if($linecont==$rowcount){
               $qh_update_ex_temp->execute();$qh_get_novelex_rowcount->execute();
               ($novelCount)=$qh_get_novelex_rowcount->fetchrow_array();
               
            }print LOG "Unique exons: $rowcount of $linecont total were loaded - $novelCount exons were novel\n";
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

