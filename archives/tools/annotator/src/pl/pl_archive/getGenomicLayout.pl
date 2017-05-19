#!/usr/bin/perl

###################################################################################################
## getGenomicLayout.pl
#  This script generates gene annotation that overlap with
#  a given probe
#
# Input: a tab-delimitted/commas separated SNPs file with the following fields:
#     1. CHR      
#     2. START
#     3. END 
#     4. STRAND
#     5. LOCAL_IDENTIFIER  
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : April 2012

#   Usage: perl getProbeAnnotations.pl -f <input_probe_file_name> -o <output file name> -v <organism version> -a <gene prediction name>
#
#   Where: -f is a file containing probe cordinates to the genome 
#          -o  is the output file name (default standard out)
#          -v  is the ucsc organism version (default mm9)
#          -a  gene prediction name (default mgi)
#          -l  display the list of current gene predictions for the specified organism version
#
###################################################################################################
use DBI;
use POSIX;

use vars qw ($opt_h $opt_l $opt_s $opt_f $opt_t $opt_c $opt_o $opt_v $opt_a);
use Getopt::Std;

$mb=1000000;
getopts('hlsCf:t:c:o:v:a:');
if($opt_h||((!$opt_f)&&(!$opt_l)&&(!$opt_s)&&(!$opt_t)&&(!$opt_c))) { #||!$opt_c
    print <<HELP;

 This script generates gene annotation that overlap with a given probe
Usage:
   perl getProbeAnnotations.pl -f <input_probe_file_name> -o <output file name> -v <organism version> -a <gene prediction name>
Arguments:
   -h  displays this help message
   -l  displays the list of current gene predictions for the specified organism version
   -s  displays the list of current organism versions found in Graber Transcript Database
   -f is a file containing genome range cordinates
   -t is a query term. gene symbol, transcript accession id,gene accession id,...
   -c is a genomic region. Format Chrom:Start:End:Strand (example -c 1:6000000:6000300:+) 
   -o  is the output file name (default standard out)(optional)
   -v  is the ucsc organism version (default mm9)(optional)
   -a  gene prediction name (default mgiGene)(optional)
   -C  returns gene result one isoform per row(default returns aggregate transcript and exons)

Example: perl getProbeAnnotations.pl -v mm9 -l
The above will display all the gene predictions for mm9

Example: perl getProbeAnnotations.pl -f probe_filetxt -o probe_file_geneAnnotation.txt -v mm9 -a ensGene
The above will use Ensembl annotations

HELP
exit;
}

my $user ='pup';
my $pwd  ='puppass';

$dbname ='graber_transcriptdb';
$host ='demon.jax.org';
$organism_version="mm9";
$colapse=0;if($opt_C){$colapse=1;}
#my $dbh = DBI->connect("DBI:mysql:$dbname:$host:$port",$user,$pwd);
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
my $qh_setDefault=$dbh->prepare("");
my $getPosTranscripts="select tx_start,tx_end,transcript_id from transcript where organism_version_id=? ";
   $getPosTranscripts.=" and chromosome_id=? and strand=? and ? between tx_start and tx_end order by tx_start";
my $qh_getPosTranscripts = $dbh->prepare($getPosTranscripts)or die "Couldn't prepare statement: " . $dbh->errstr;

my $getOverlap ="select tx_start,tx_end,transcript_id from transcript where organism_version_id=?";
   $getOverlap .=" and chromosome_id=? and strand=? and tx_end >= ? and tx_start<=? order by tx_start";
my $qh_Overlaplist=$dbh->prepare($getOverlap)or die "Couldn't prepare statement: " . $dbh->errstr;  

my $getTxCord="select organism_version_id,chromosome_id,strand,tx_start,tx_end 
         from transcript where transcript_id=? and organism_version_id=?";
my $qh_getTxCord = $dbh->prepare($getTxCord)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getExCord="select organism_version_id,chromosome_id,strand,exon_start,exon_end 
          from exon where exon_id=? and organism_version_id=?";
my $qh_getExCord = $dbh->prepare($getExCord)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getOverlapExon="select exon_start,exon_end,exon_id from exon where exon_id=? and organism_version_id=? ";
   $getOverlapExon.=" and chromosome_id=? and strand=? and exon_end >= ? and exon_start<=?  order by exon_start";
my $qh_getOverlapExon = $dbh->prepare($getOverlapExon)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getOverlapE="select exon_start,exon_end,exon_id from exon where organism_version_id=? ";
   $getOverlapE.=" and chromosome_id=? and strand=? and exon_end >= ? and exon_start<=?  order by exon_start";
my $qh_getOverlapE = $dbh->prepare($getOverlapE)or die "Couldn't prepare statement: ".$dbh->errstr;


#select* from gene_by_annotation where gene_name="Cfl1" and gene_prediction_id=18;
my $getGeneTx="select transcript_id from gene_by_annotation where gene_name=? and gene_prediction_id=?";
my $qh_getGeneTx = $dbh->prepare($getGeneTx)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTxGene="select gene_name from gene_by_annotation where transcript_id=? and gene_prediction_id=?";
my $qh_getTxGene = $dbh->prepare($getTxGene)or die "Couldn't prepare statement: ".$dbh->errstr;


my $getTxExid="select exon_id from transcript_exon where transcript_id=? and gene_prediction_id=?";
my $qh_getTxExid = $dbh->prepare($getTxExid)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getExidTx="select transcript_id from transcript_exon where exon_id=? and gene_prediction_id=?";
my $qh_getExidTx = $dbh->prepare($getExidTx)or die "Couldn't prepare statement: ".$dbh->errstr;


my $getTranscriptTx="select transcript_id from transcript_by_annotation where transcript_name=? and gene_prediction_id=?";
my $qh_getTranscriptTx = $dbh->prepare($getTranscriptTx)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getChr="select distinct t.chromosome_id, c.chromosome_name from transcript t, chromosome c 
             where organism_version_id=? and t.chromosome_id=c.chromosome_id";
my $qh_getChr = $dbh->prepare($getChr)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getSingleChr="select chromosome_id from chromosome where chromosome_name=?";
my $qh_getSingleChr = $dbh->prepare($getSingleChr)or die "Couldn't prepare statement: ".$dbh->errstr;


my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $qh_getOrgV = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getOrgUcsc="select ucsc_db  from organism_version where organism_version_id=?";
my $qh_getOrgUcsc = $dbh->prepare($getOrgUcsc)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getAllOrgV="select distinct o.organism_name,ov.organism_version_id from organism o,organism_version ov
                 where o.organism_id=ov.organism_id order by organism_name";
my $qh_getAllOrgV = $dbh->prepare($getAllOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getGenePred="select gene_prediction_name from gene_prediction_by_organism op, gene_prediction p
                   where op.organism_version_id =? and op.gene_prediction_id=p.gene_prediction_id 
                   and gene_prediction_name not in ('all_est','all_mrna','estOrientInfo')";
my $qh_getGenePred = $dbh->prepare($getGenePred)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getPredId="select gene_prediction_id from gene_prediction where gene_prediction_name =?";
my $qh_getPredId = $dbh->prepare($getPredId)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTermId="select term,term_type_id from query_terms where term=?";
my $qh_getTermId = $dbh->prepare($getTermId)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getMGIgene="select mgi_symbol from mgi_genes where mgi_id=?";
my $qh_getMGIgene = $dbh->prepare($getMGIgene)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTxAccession="select transcript_name from transcript_by_annotation where transcript_id=? and gene_prediction_id=?";
my $qh_getTxAccession = $dbh->prepare($getTxAccession)or die "Couldn't prepare statement: ".$dbh->errstr;
#getTxExons($tx_id,$org_vid,$pred_id,$chrom_start,$chrom_end){

  $query=" select distinct t.exon_id,t.exon_start,t.exon_end ";
  $query.=" from exon t,transcript_exon a ";
  $query.=" where  a.transcript_id=? and a.organism_version_id=?  ";
  $query.=" and a.gene_prediction_id=? and a.transcript_id=t.transcript_id ";
  $query.=" and a.organism_version_id=a.organism_version_id and t.exon_end >= ? and t.exon_start<=?  order by exon_start";
my $qh_getTxExons=$dbh->prepare($query);
#set defaults
$organism_version=$opt_v if($opt_v);$org_vid=0;$organism_version=~s/\s+//g;
$qh_getOrgV->execute("$organism_version");
if($qh_getOrgV->rows>0){($org_vid)=$qh_getOrgV->fetchrow_array();}
$gene_prediction="";$pred_id=0;%predictions=();
if($opt_a){$gene_prediction=$opt_a;} 
else{ #get the list of annotation sources for this organism version
       #and set the default annotation for this organism
    $qh_getGenePred->execute($org_vid);
    if($qh_getGenePred->rows>0){
       while(($gene_prediction)=$qh_getGenePred->fetchrow_array()){
            $gene_prediction=~s/\s*//g;$predictions{$gene_prediction}=1;
       }
    }
}
if($predictions{"mgiGene"}==1){$gene_prediction="mgiGene";}
elsif($predictions{"knownGene"}==1){$gene_prediction="knownGene";}
elsif($predictions{"ensGene"}==1){$gene_prediction="ensGene";}
elsif($predictions{"refGene"}==1){$gene_prediction="refGene";}
elsif($predictions{"xenoRefGene"}==1){$gene_prediction="xenoRefGene";}
$qh_getPredId->execute("$gene_prediction");if($qh_getPredId->rows>0){($pred_id)=$qh_getPredId->fetchrow_array();}
if($pred_id<=0){
  print "************************\nNote: We do not have $organism_version data from the specified annotation ($gene_prediction)\n";
  print "Please run the script with the -l option to display the current annotations list.\n";
  exit(0);
}
#################################################################################
sub getGenomeAnnotation{
    $annotation="";
   ($chrom_id,$strand,$chrom_start,$chrom_end,$pred_id,$org_vid,$line,$current_txid)= @_;
   $qh_Overlaplist->execute($org_vid,$chrom_id,$strand,$chrom_start,$chrom_end);
   $region_start=0;$region_end=0; %txmap=();%genes=();$gene_names=""; $tx="";$ex_starts="--";
   $ex_ends="--";
   if($qh_Overlaplist->rows>0){
      while(($tx_start,$tx_end,$transcript_id)=$qh_Overlaplist->fetchrow_array()){
            next if(exists($txmap{$transcript_id}));
            next if(($current_txid>0) &&($current_txid!=$transcript_id));
            $qh_getTxAccession->execute($transcript_id,$pred_id);
            next if($qh_getTxAccession->rows<=0);
            if(($region_start>$tx_start)||($region_start==0)){$region_start=$tx_start;}
            if($region_end<$tx_end){$region_end=$tx_end;}
            
            ($name)=$qh_getTxAccession->fetchrow_array();
            if($tx eq ""){$tx.="$name";}else{$tx.=",$name";}
                  #get genes list
            $qh_getTxGene->execute($transcript_id,$pred_id);
            if($qh_getTxGene->rows>0){
               while(($gene_symbol)=$qh_getTxGene->fetchrow_array()){
                   next if(exists($genes{$gene_symbol}));
                   if($gene_names eq ""){$gene_names.="$gene_symbol";}else{$gene_names.=",$gene_symbol";}
                   $genes{$gene_symbol}=1;
               }
             }
             $txmap{$transcript_id}=1;
       }
     }
     $qh_getOverlapE->execute($org_vid,$chrom_id,$strand,$chrom_start,$chrom_end);
     if($qh_getOverlapE->rows>0){ 
        if($qh_getOverlapE->rows>0){
           $current_start=0;$current_end=0;$ex_start=0;$ex_end=0;
           while(($ex_start,$ex_end,$exon_id)=$qh_getOverlapE->fetchrow_array()){
                  $qh_getExidTx->execute($exon_id,$pred_id);%etxmap=();
                  next if($qh_getExidTx->rows<=0);#this exon is not part of this gene annotation
                  while(($this_txid)=$qh_getExidTx->fetchrow_array()){$etxmap{$this_txid}=1;}
                  if($current_txid>0){
                     next if(!exists($etxmap{$current_txid}));}
                  if($current_start==0){$current_start=$ex_start;$current_end=$ex_end;}
                  else{
                     if(($ex_start>=$current_start)&&($ex_start<=$current_end)){
                        if($current_end<$ex_end){$current_end=$ex_end;}
                     }
                     else{ #this is a new genomic region, there is no overlap with previous exon
                        if($ex_starts eq "--"){$ex_starts="$current_start";}else{$ex_starts.=",$current_start";}
                        if($ex_ends eq "--"){$ex_ends="$current_end";}else{$ex_ends.=",$current_end";}
                        $current_start=$ex_start;$current_end=$ex_end;
                     }
                  }
            }
             #now add the last exon
            if($ex_starts eq "--"){
                if($current_start>0){
                  $ex_starts="$current_start";$ex_ends="$current_end";
                }
            }
            else{
                if($current_start>0){
                  $ex_starts.=",$current_start";$ex_ends.=",$current_end" ;
                }
            }
          }
     }
   if(keys(%txmap)>0){
      $annotation= "$line\t$gene_names\t";
      $annotation.= "$region_start\t$region_end\t$tx\t$ex_starts\t$ex_ends\n";
   }
 return $annotation;
}
############################################################################
sub getIsoformAnnotation{
    $annotation="";
   ($chrom_id,$strand,$chrom_start,$chrom_end,$pred_id,$org_vid,$line,$current_txid)= @_;
   $qh_Overlaplist->execute($org_vid,$chrom_id,$strand,$chrom_start,$chrom_end);
   $region_start=0;$region_end=0; %txmap=();$gene_names=""; $tx="";
   if($qh_Overlaplist->rows>0){
      while(($tx_start,$tx_end,$transcript_id)=$qh_Overlaplist->fetchrow_array()){
            next if(exists($txmap{$transcript_id}));
            next if(($current_txid>0) &&($current_txid!=$transcript_id));
            $qh_getTxAccession->execute($transcript_id,$pred_id);
            next if($qh_getTxAccession->rows<=0);
            ($tx)=$qh_getTxAccession->fetchrow_array();
            #get genes list
            $qh_getTxGene->execute($transcript_id,$pred_id);$gene_names="";%genes=();
            if($qh_getTxGene->rows>0){
               while(($gene_symbol)=$qh_getTxGene->fetchrow_array()){
                   next if(exists($genes{$gene_symbol}));
                   if($gene_names eq ""){$gene_names.="$gene_symbol";}else{$gene_names.=",$gene_symbol";}
                   $genes{$gene_symbol}=1;
               }
             }
            ##get the exons of this transcript
            ###==========================
             $qh_getOverlapE->execute($org_vid,$chrom_id,$strand,$chrom_start,$chrom_end);
             $current_start=0;$current_end=0;$ex_start=0;$ex_end=0;
             $ex_starts="--";$ex_ends="--";
             while(($ex_start,$ex_end,$exon_id)=$qh_getOverlapE->fetchrow_array()){
                  $qh_getExidTx->execute($exon_id,$pred_id);%etxmap=();
                  next if($qh_getExidTx->rows<=0);#this exon is not part of this gene annotation
                  while(($this_txid)=$qh_getExidTx->fetchrow_array()){$etxmap{$this_txid}=1;}
                  if($current_txid>0){
                     next if(!exists($etxmap{$current_txid}));}
                  next if(!exists($etxmap{$transcript_id}));
                  if($current_start==0){$current_start=$ex_start;$current_end=$ex_end;}
                  else{
                     if(($ex_start>=$current_start)&&($ex_start<=$current_end)){
                        if($current_end<$ex_end){$current_end=$ex_end;}
                     }
                     else{ #this is a new genomic region, there is no overlap with previous exon
                        if($ex_starts eq "--"){$ex_starts="$current_start";}else{$ex_starts.=",$current_start";}
                        if($ex_ends eq "--"){$ex_ends="$current_end";}else{$ex_ends.=",$current_end";}
                        $current_start=$ex_start;$current_end=$ex_end;
                     }
                  }
            }
            if($ex_starts eq "--"){if($current_start>0){$ex_starts="$current_start";$ex_ends="$current_end";}}
            else{if($current_start>0){$ex_starts.=",$current_start";$ex_ends.=",$current_end";}}
            $annotation.= "$line\t$gene_names\t";
            $annotation.= "$tx_start\t$tx_end\t$tx\t$ex_starts\t$ex_ends\n";
            ###==========================
       }
     }
 return $annotation;
}
############################################################################
if($opt_f){
   if(!(-f $opt_f)){print "File name $opt_f does not exist.\n";exit(1);}
   @header=`head -1 $opt_f`;$has_header=0;$count=0;%fieldmap=();@field_names=();
   $token="chrom\tchromStart\tchromEnd\tstrand";
   if(@header>0){
      if($header[0]=~/strand/i){$has_header=1;chomp($header[0]);$i=0;
         @field_names=split("\t",$header[0]);$token=$header[0];
         while($i<@field_names){$fieldmap{lc($field_names[$i])}=$i;++$i;}
      }
   }
   %chrom_map=();$qh_getChr->execute($org_vid);
   if($qh_getChr->rows>0){
       while(($chr_id,$chr_name)=$qh_getChr->fetchrow_array()){$chrom_map{$chr_name}=$chr_id;}
    }
   open(IN,"$opt_f")or die "Could not open $opt_f: $!\n";
   $head= "***********----------------- **************************\n";
   $head.= "Generating annotations using Graber Transcrip Database.\n";
   $head.= "Your data file $opt_f will be annotated using :\n";
   $head.= "Organism version =$organism_version \nGene Prediction = $gene_prediction\n";
   $head.= "***********----------------- **************************\n";
   $head.= "$token\tgene(s)\t";
   $head.= "Annotation_start\tAnnotation_end\tTranscript(s)\tExon_starts\tExon_ends\n";
   if(IN){
      if($opt_o){open(OUT,">$opt_o")or die "Could not open $opt_o: $!\n";}
      if(-f $opt_o){print OUT "$head";}
      else{print"$head";}
      $chrom;$chrom_start;$chrom_end;$strand;$count=0;
      while($line=<IN>){
        chomp($line);@fields=split("\t",$line);++$count; 
        next if(($count==1)&&($has_header==1));
        if($has_header==0){ #if the file does not have a header, the default is
           if(@fields>=4){   #column1= chrom, column2=chromStart,column3=chromEnd,column4=strand
              $chrom=$fields[0];$chrom_start=$fields[1];$chrom_end=$fields[2];$strand=$fields[3];
           }
         }else{
             next if(@field_names!=@fields);
             $chrom=$fields[$fieldmap{"chrom"}];$chrom_start=$fields[$fieldmap{"chromStart"}];
             $chrom_end=$fields[$fieldmap{"chromEnd"}];$strand=$fields[$fieldmap{"strand"}];
         }
         next if(($chrom_start<=0)||($chrom_end<=0)||!($strand=~/\+|-/));
         $chrom=~s/chr//i;$chrom=~s/ch//i;$chrom_id=$chrom_map{$chrom};$tx_id=0;
         #convert cordinates into BP if entered as MBp
         $chrom_end=$chrom_end*$mb if($chrom_end=~/\d+\.\d+/);
         $chrom_start=$chrom_start*$mb if($chrom_start=~/\d+\.\d+/);
         if(!$colapse){
           $annotations=getGenomeAnnotation($chrom_id,$strand,$chrom_start,$chrom_end,$pred_id,$org_vid,$line,$tx_id);
         }
         else{
           $annotations=getIsoformAnnotation($chrom_id,$strand,$chrom_start,$chrom_end,$pred_id,$org_vid,$line,$tx_id);
         }
         if(-f $opt_o){print OUT "$annotations";}else{print "$annotations";}
      }
    close(IN);
   }
  
}
elsif($opt_c){
  ($chrom,$chrom_start,$chrom_end,$strand)=split(":",$opt_c);
  if(($chrom_start<=0)||($chrom_end<=0)||!($strand=~/\+|-/)){
     print "Bad for the genomic region format\n It should be: -c Chrom:ChromStart:ChromEnd:Strand\n";
     print "For example: -c X:50000000:50000458:+ \n";
  }
  else{
      #get all the transcripts that overlap the region
     $chrom=~s/chr//i;$chrom=~s/ch//i;
     $qh_getSingleChr->execute($chrom);
     if($qh_getSingleChr->rows>0){
        ($chrom_id)=$qh_getSingleChr->fetchrow_array();
        $head= "***********----------------- **************************\n";
        $head.= "Generating annotations using Graber Transcrip Database.\n";
        $head.= "Your data $opt_c will be annotated using :\n";
        $chrom_end=$chrom_end*$mb if($chrom_end=~/\d+\.\d+/);#convert cordinates into BP if entered as MBp
        $chrom_start=$chrom_start*$mb if($chrom_start=~/\d+\.\d+/);
        $head.= "Chromosome=$chrom; ChromStart=$chrom_start bp; ChromEnd=$chrom_end bp; Strand=$strand\n";
        $head.= "Organism version =$organism_version \nGene Prediction = $gene_prediction\n";
        $head.= "***********----------------- **************************\n";
       $head.= "chrom\tchromStart\tchromEnd\tstrand\tgeneNames\t";
       $head.= "Annotation_start\tAnnotation_end\ttranscriptNames\texonStarts\texonEnds\n";
        
        print "$head";$line="$chrom\t$chrom_start\t$chrom_end\t$strand";$tx_id=0;
        if($colapse==0){
           $annotations=getGenomeAnnotation($chrom_id,$strand,$chrom_start,$chrom_end,$pred_id,$org_vid,$line,$tx_id);
         }
         else{
           $annotations=getIsoformAnnotation($chrom_id,$strand,$chrom_start,$chrom_end,$pred_id,$org_vid,$line,$tx_id);
         }
        print "$annotations";
     }
  }
}
elsif($opt_t){ #query by term, get this term type
  $qh_getTermId->execute("$opt_t");
  if($qh_getTermId->rows>0){
    print "***********----------------- **************************\n";
    print "Generating annotations using Graber Transcrip Database.\n";
    print "Your query:$opt_t,will be annotated using :\n";
    print "Organism version =$organism_version \nGene Prediction = $gene_prediction\n";
    print "***********----------------- **************************\n";
    print "chrom\tchromStart\tchromEnd\tstrand\tgeneNames\t";
    print "Annotation_start\tAnnotation_end\ttranscriptNames\texonStarts\texonEnds\n";
     while(($term,$term_type_id)=$qh_getTermId->fetchrow_array()){
        $mgi_symbol=$term;
        if($term_type_id==4){
           $qh_getMGIgene->execute("$term");
           if($qh_getMGIgene->rows>0){
              ($mgi_symbol)=$qh_getMGIgene->fetchrow_array();
            }
        }
         if(($term_type_id==2)||($term_type_id==4)){ ##the term is a gene symbol
            $qh_getGeneTx->execute($mgi_symbol,$pred_id) ;
            if($qh_getGeneTx->rows>0){
               $chr_id=0;$gene_start=0;$gene_end=0;$gstrand="";
               while(($tx_id)=$qh_getGeneTx->fetchrow_array()){#get the cordidates of this transcript 
                   $qh_getTxCord->execute($tx_id,$org_vid);
                   if($qh_getTxCord->rows>0){
                        ($organism_version_id,$chrom_id,$strand,$chrom_start,$chrom_end)=$qh_getTxCord->fetchrow_array();
                        if(($gene_start==0)||($gene_start>$chrom_start)){$gene_start=$chrom_start;}
                        if(($gene_end==0)||($gene_end<$chrom_end)){$gene_end=$chrom_end;}
                         $chr_id=$chrom_id;$gstrand=$strand;
                   }
               }
              $line="$chr_id\t$gene_start\t$gene_end\t$gstrand";$tx_id=0;
              if($colapse==0){
                $annotations=getGenomeAnnotation($chr_id,$strand,$gene_start,$gene_end,$pred_id,$org_vid,$line,$tx_id);
              }
              else{
                $annotations=getIsoformAnnotation($chr_id,$strand,$gene_start,$gene_end,$pred_id,$org_vid,$line,$tx_id);
               }
              print "$annotations";
            }
         }
         elsif($term_type_id==3){##the term is a transcript
             $qh_getTranscriptTx->execute($term,$pred_id);
             if($qh_getTranscriptTx->rows>0){
                ($tx_id)=$qh_getTranscriptTx->fetchrow_array();
                 $qh_getTxCord->execute($tx_id,$org_vid);
                if($qh_getTxCord->rows>0){
                   ($organism_version_id,$chrom_id,$strand,$chrom_start,$chrom_end)=$qh_getTxCord->fetchrow_array();
                   $line="$chrom_id\t$chrom_start\t$chrom_end\t$strand";
                  if($colapse==0){
                      $annotations=getGenomeAnnotation($chrom_id,$strand,$chrom_start,
                                                      $chrom_end,$pred_id,$org_vid,$line,$tx_id);
                   }
                  else{
                      $annotations=getIsoformAnnotation($chrom_id,$strand,$chrom_start,
                                                     $chrom_end,$pred_id,$org_vid,$line,$tx_id);
                   }
                   print "$annotations";
                }
             }   
         }
     }
  }
}
else{
  if($opt_l&& !$opt_s){ #display gene prediction listing for this organism version
    $qh_getGenePred->execute($org_vid);
    if($qh_getGenePred->rows>0){
      print "*****************************************\n";
      print "Gene Predictions for $organism_version\n";
      print "*****************************************\n"; $i=1;
      while(($gene_prediction)=$qh_getGenePred->fetchrow_array()){
            if(($gene_prediction=~/gene/i)||($gene_prediction=~/refFlatSquishMulti/i)){
                print "$i.) $gene_prediction\n";++$i;
            }
      }
   }
 }
 if($opt_s){#display all the organism version found in Graber Transcript database
   $qh_getAllOrgV->execute();
   if($qh_getAllOrgV->rows>0){
      print "*****************************************\n";
      print "Organisms found in Graber Transcript Database\n";
      print "*****************************************\n"; $i=1;
      while(($organism,$id)=$qh_getAllOrgV->fetchrow_array()){
            $qh_getOrgUcsc->execute($id);
            if($qh_getOrgUcsc->rows>0){
               ($uscs_db)=$qh_getOrgUcsc->fetchrow_array();
               if($opt_l){
                  $qh_getGenePred->execute($id);
                  if($qh_getGenePred->rows>0){
                     while(($gene_prediction)=$qh_getGenePred->fetchrow_array()){
                           if(($gene_prediction=~/gene/i)||($gene_prediction=~/refFlatSquishMulti/i)){
                               print "$i.) $organism ($uscs_db)\t $gene_prediction\n";
                          }
                      }
                  }
                  ++$i;
               }
               else{print "$i.) $organism ($uscs_db)\n";++$i;}
            }
       }
    }
  }
}
print"program completed\n";
 
exit(0);

