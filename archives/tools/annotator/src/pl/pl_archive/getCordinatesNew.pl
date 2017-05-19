#!/usr/bin/perl

###################################################################################################
## getCordinates.pl
#  This script generates feature cordinates of the specified term (gene, transcript, genomic region)
#  for a given organism version and gene prediction source. The default organism version is mm9
#  and the default gene prediction is mm9-ensGene.
#  The features are stored into the database using zero-based cordinates standard for feature starts
#  and 1-based cordinate standard for feature ends
#    
# Note: The output of this program can be used by other programs, for example snpAnnotator, getTxFastaSeq 
# Output: either a genePred format (transcript query) 
#         or a gtf format (feature querie, ex, int,...)
#
# 1. a gtf file has the following fields
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
# <seqname> is the chromosome
# <source> is the gene prediction source name
# <feature> is feature type (exon, CDS)
# [attributes] a list of attributes (gene_id;transcript_id;gene_name;feature_id; feature_rank;...)
# 2. a genePred file has the following fields
# name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tid\tname2\texonFrames\n

# Features coordinates in the result are generated in 1-based cordinates standard for <start>, <end>
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : March 2013

#   Usage: perl getCordinates.pl -t <feature type> [-T term]
#           -v <organism version> [-a <gene prediction list>] [-o <outputDir>] [-d db_host]
# Note: 
# 1. If the option -T is not specified, the program will generate all the feature coordinates
#   for the specified organism version and feature type

# 2. if the -a option is not specified, features from all the annotation sources are generated
#  the features are displayed by annotation source

###################################################################################################
use DBI;
use POSIX;
use XML::Simple;

# create object
$xml = XML::Simple->new (ForceArray => 1);

use vars qw ($opt_h $opt_s $opt_o $opt_t $opt_v $opt_a $opt_l $opt_f $opt_d $opt_F $opt_T $opt_g);
use Getopt::Std;

getopts('hlsgo:t:v:a:f:d:FT:');
if($opt_h||(!$opt_v&&!$opt_s&&!$opt_l&&!$opt_g)) {
    print <<HELP;

This script generates the coordinates  of transcripts/features for the specified organism version and annotation source.
The features are stored into the database using the zero-based cordinates standard for feature starts
  and 1-based cordinate standard for feature ends
    
Note: The output of this program can be used by other programs, for example snpAnnotator, getTxFastaSeq 
      result coordinates are in 1-base. Both starts and ends

Output: either a genePred format (transcript query) or a gtf format (feature queries)

Usage: perl getCordinates.pl -t <feature type> -v <organism version>[-T pax6] [-a <gene prediction list>] 
           [-o <outputDir>][-f outfilename]

Example: perl getCordinates.pl -v mm9 -l
The above will display all the gene prediction sources for mm9 currently in graber_transcript database

Example: perl getCordinates.pl -s
The above will display all the organism versions currently in graber_transcript database

Example: perl getCordinates.pl -t tx -o features -v mm10 -a GRCm38-ensGene
  The above will use Ensembl annotations for mm10 to generate all the transcripts
  in gene prediction format and store the file (mm10-GRCm38-ensGene-transcripts.txt) under features/

Example: perl getCordinates.pl -v mm10 -T pax6 -f pax6-transcripts.txt
  The above will use all pax6 annotations for mm10 and generate transcripts data
  in gene prediction format and store the file (pax6-transcripts.txt) under the current working directory

Example: perl getCordinates.pl -t tx -o features -v mm10 -a GRCm38-ensGene -F 
The above will use Ensembl annotations for mm10 to generate all the transcripts
in gtf format and store the file (mm10-GRCm38-ensGene-transcripts.gtf) under features/


Example: perl getCordinates.pl -t tx -v mm10 -a GRCm38-ensGene -o features -f mm10-ensemblTransctipts.txt
The above will use Ensembl annotations for mm10 to generate all the transcripts
in gene prediction format and store the file (mm10-ensemblTransctipts.txt) under features/. 

Example: perl getCordinates.pl -t tx -v mm10 -a GRCm38-ensGene -o features -f mm10-ensemblTransctipts.gff -F gtf
The above will use Ensembl annotations for mm10 to generate all the transcripts
in gtf format and store the file (mm10-ensemblTransctipts.gtf) under features/. 

Arguments:
   -h  displays this help message
  
   -t  feature type: 
       tx     -> Transcripts - default (gene prediction format or gtf format)
       ex     -> Exons      (gtf format)
       fex    -> First Exons (gtf format)
       lex    -> last exons  (gtf format)
       tex    -> terminal exons (5' and 3' exons) (gtf format)
  -T  query term (gene name, transcript, genomic region) 
  -v  is the ucsc organism version (default mm9)
  -a  optional, a commas separated list of gene prediction names to use
  -o  optional, is the output directory, default current working directory
  -f  optional, output file name (relative to the specified output directory)
  -F  optional, output format, default genePred format (only used for the -t tx )
  -l  display the list of current gene predictions for the specified organism version
  -s  display all the organism versions currently in graber_transcript database
  -g  display all the organism versions that have genome sequence (chromosome)
 
HELP
exit;
}

#set defaults
my $user ='pup';my $pwd  ='puppass';$dbname ='graber_transcriptdb';$host ='harlequin.jax.org';
if($opt_d){$host=$opt_d;}
$webservice="http://demon.jax.org/transcriptdb/web-services/index.php";
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
##################################################
# get hits for the search term
##################################################
$query="select distinct g.transcript_id,transcript_name,t.gene_prediction_id 
        from gene_by_annotation g, transcript_by_annotation t 
        where gene_name =? and g.organism_version_id=?  and g.transcript_id=t.transcript_id
        and  g.organism_version_id=t.organism_version_id";
$qh_getGenTxHits=$dbh->prepare($query);
$query="select distinct transcript_id,transcript_name,gene_prediction_id 
        from transcript_by_annotation where transcript_name =? and organism_version_id=? ";
$qh_getTxHits=$dbh->prepare($query);
$query="select distinct t.transcript_id,transcript_name,t.gene_prediction_id  
          from gene_by_annotation g, mgi_genes m, transcript_by_annotation t ";
$query .=" where m.mgi_id =? and m.mgi_symbol=g.gene_name and g.organism_version_id=?
           and g.transcript_id=t.transcript_id and  g.organism_version_id=t.organism_version_id ";
$qh_getMgiHits=$dbh->prepare($query);

$query="select distinct t.transcript_id,t.transcript_name,t.gene_prediction_id 
        from gene_by_annotation g,cc_founders_genes c,transcript_by_annotation t 
        where  c.gene_id =? and c.gene_name=g.gene_name and g.gene_prediction_id=c.gene_prediction_id 
        and g.organism_version_id=? and g.transcript_id=t.transcript_id and  g.organism_version_id=t.organism_version_id  ";
$qh_getCcGenHits=$dbh->prepare($query);
$query="select distinct t.transcript_id,t.transcript_name,t.gene_prediction_id 
        from transcript_by_annotation t,cc_founders_transcripts c
        where  c.transcript_name=? and t.organism_version_id=? and c.gene_prediction_id = t.gene_prediction_id 
        and    t.transcript_name=c.transcript_id ";
$qh_getCcTxHits=$dbh->prepare($query);

$query=" select transcript_id from transcript where organism_version_id=? ";
$query.=" and chromosome_id=?  and strand=(?) and tx_end >=? and tx_start<=? ";
$qh_getCordHits= $dbh->prepare($query);
$query=" select distinct t.transcript_id,transcript_name, gene_prediction_id 
         from transcript t, transcript_exon e where t.organism_version_id=? ";
$query.="and chromosome_id=? and tx_end >=? and tx_start<=?  and t.organism_version_id=e.organism_version_id 
         and t.transcript_id=e.transcript_id";
$qh_getCordAllHits= $dbh->prepare($query);
$query="select distinct transcript_id,transcript_name,gene_prediction_id 
        from transcript_exon where organism_version_id=? and gene_prediction_id in (?)";
$qh_getTxByAnnot=$dbh->prepare($query);
$query="select distinct transcript_id,transcript_name,gene_prediction_id 
        from transcript_exon where organism_version_id=? and gene_prediction_id!=3";
$qh_getTxByAllAnnot=$dbh->prepare($query);

$qh_getChr=$dbh->prepare("select chromosome_id from chromosome where chromosome_name=?");
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$qh_getPredName=$dbh->prepare("select gene_prediction_name from gene_prediction where gene_prediction_id=?");

$query="select distinct gene_prediction_name,g.gene_prediction_id from gene_prediction g,transcript_exon t";
$query.=" where  t.organism_version_id=? and t.transcript_id=? and transcript_name=? ";
$query.=" and t.gene_prediction_id=g.gene_prediction_id";
$qh_getTxPredId=$dbh->prepare($query);

$query="select distinct transcript_id,transcript_name, gene_prediction_id from transcript_exon
       where organism_version_id=? and gene_prediction_id in (?)";
$qh_getTxByAnnot=$dbh->prepare($query);
$query="select distinct transcript_id,transcript_name,gene_prediction_id from transcript_exon 
       where organism_version_id=? and gene_prediction_id!=3";
$qh_getTxByAllAnnot=$dbh->prepare($query);

my $getTxcord="select chromosome_name,strand,tx_start,tx_end from transcript e, chromosome c";
   $getTxcord.=" where transcript_id=? and organism_version_id=? and e.chromosome_id=c.chromosome_id";
   $qh_getTxcord = $dbh->prepare($getTxcord)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTxExons="select distinct chromosome_name,exon_start,exon_end,strand,exon_frame from transcript_exon t,exon e,chromosome c ";
   $getTxExons.=" where  transcript_id=? and transcript_name=? and t.organism_version_id=? and gene_prediction_id in (?)";
   $getTxExons.=" and t.organism_version_id=e.organism_version_id and t.exon_id=e.exon_id and e.chromosome_id=c.chromosome_id";
   $getTxExons.=" order by exon_start asc";
my $qh_getTxExons = $dbh->prepare($getTxExons)or die "Couldn't prepare statement: ".$dbh->errstr;
my $qh_getExonCord = $dbh->prepare("select exon_start,exon_end from exon where exon_id in (?) order by exon_start");

my $qh_getCCgene=$dbh->prepare("select distinct gene_prediction_id, gene_id, gene_name from cc_founders_genes where gene_prediction_id in (?)");
my $qh_getCCtx=$dbh->prepare("select distinct transcript_id,transcript_name,protein_id from cc_founders_transcripts where gene_prediction_id=?");

my $getTxExons="select distinct e.exon_id,exon_start,exon_end,strand from transcript_exon t,exon e ";
   $getTxExons.=" where  transcript_id=? and transcript_name=? and t.organism_version_id=? and gene_prediction_id!=3";
   $getTxExons.=" and t.organism_version_id=e.organism_version_id and t.exon_id=e.exon_id ";
   $getTxExons.=" order by exon_start asc";
my $qh_getTxExonsAll = $dbh->prepare($getTxExons)or die "Couldn't prepare statement: ".$dbh->errstr;

$query="select cdsStart,cdsEnd from transcript_translation ";
$query.=" where transcript_id in (?) and transcript_name in(?) and organism_version_id=? and gene_prediction_id in (?) order by cdsStart asc";
my $qh_getCDS=$dbh->prepare($query);

my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $qh_getOrgV = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getOrgUcsc="select ucsc_db  from organism_version where organism_version_id=?";
my $qh_getOrgUcsc = $dbh->prepare($getOrgUcsc)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTxGene="select distinct gene_name from gene_by_annotation 
                where organism_version_id=? and transcript_id=? and gene_prediction_id in (?)";
my $qh_getTxGene = $dbh->prepare($getTxGene)or die "Couldn't prepare statement: ".$dbh->errstr;
my $qh_getTxGeneAll = $dbh->prepare("select distinct gene_name from gene_by_annotation where organism_version_id=? and transcript_id=?");

%ccgenemap=();

sub displayAnnotSources{
 ($org_version)=@_;$url="$webservice?v=$org_version&l=1";
  $xml_file="annot_list.xml";if(-e "annot_list.xml"){system("rm $xml_file");}
  system("wget -q -O $xml_file '$url'");
   $xml = XML::Simple->new (ForceArray => ['prediction','name']);
  if(-e $xml_file){$data =$xml->XMLin("$xml_file");
    foreach my $source(@{$data->{prediction}}){print $source->{name}[0] . "\n";}
  }if(-e "annot_list.xml"){system("rm $xml_file");}
}
sub displayOrganismVersions{
  ($type)=@_;
  $url="$webservice?$type";$xml_file="annot_list.xml";if(-e "annot_list.xml"){system("rm $xml_file");}
  system("wget -q -O $xml_file '$url'");
  $xml = XML::Simple->new (ForceArray =>1);
  if(-e $xml_file){$data =$xml->XMLin("$xml_file");
    foreach my $organism(@{$data->{organism}}){print $organism->{name}[0] ."\t".$organism->{version}[0]. "\n";}
  }if(-e "annot_list.xml"){system("rm $xml_file");}
}
############################
sub getGenomicRegion{
   ($term)=@_;($chrom,$posData)=split(":",$term);$chromosome_id=0;$chrom_start=0;$chrom_end=0;
  if($posData){$posData=~s/\s*//g;
    if($posData=~/^(\d+.*)-(\d+.*)$/){
       $chrom_start=$1;$chrom_end=$2;$chrom_start=~s/\s+//g;$chrom_end=~s/\s+//g;
       $chrom=~s/CHR//i;$qh_getChr->execute($chrom);($chromosome_id)=$qh_getChr->fetchrow_array() if($qh_getChr->rows>0);
      if($chromosome_id>0){ $mb=1000000;
         if(($chrom_start>0)&&($chrom_end>0)){ #convert cordinates into BP if needed
             $chrom_end *=$mb if($chrom_end=~/\d+\.\d+/);$chrom_start*=$mb if($chrom_start=~/\d+\.\d+/);
         }
      }}
  }
  return "$chromosome_id:$chrom_start:$chrom_end";
}
##returns either the first exon or the last exon
sub getLastTerminalExon{
   ($tx_id,$tx,$org_vid,$prediction_id,$opt)=@_;
   $qh_getTxExons->execute($tx_id,$tx,$org_vid,$prediction_id);$termExon="";
   $minExonStart=0;$maxExonStart=0;$minExon="";$maxExon="";$tstrand ="";
   $firstExon="";$lastExon="";
   while(($chrom,$exstart,$exend,$strand,$frame)=$qh_getTxExons->fetchrow_array()){
        $tstrand =$strand if($tstrand eq "");$exstart+=1;
        if(($minExonStart==0)||($minExonStart>$exstart)){
            $minExonStart=$exstart;$minExon="$chrom,$exstart,$exend,$strand";}
        if(($maxExonStart==0)||($maxExonStart<$exstart)){
            $maxExonStart=$exstart;$maxExon="$chrom,$exstart,$exend,$strand";}
   } 
  $lastExon=($tstrand eq "+")?$maxExon:$minExon;
  $firstExon=($tstrand eq "+")?$minExon:$maxExon;
  $termExon="$firstExon:$lastExon";
 return $termExon;
}
## returns feature list 
## <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
sub getFeatures{
   ($tx_id,$tx,$org_vid,$prediction_id,$opt)=@_;
   $qh_getPredName->execute($prediction_id);($annot)=$qh_getPredName->fetchrow_array();
   $qh_getTxExons->execute($tx_id,$tx,$org_vid,$prediction_id);$exoncount=$qh_getTxExons->rows; if($exoncount<=0){return $rows;}
   #$minStart=0;$maxStart=0;$minCord="";$maxCord="";
   $prank=1; $nrank=$exoncount;$rank=0; $rows=""; $featureType="";
   $qh_getTxGene->execute($org_vid,$tx_id,$prediction_id);($gene_id)=$qh_getTxGene->fetchrow_array();
   $gene_name="";$strand="";$qh_getTxGene->execute($org_vid,$tx_id,$prediction_id); $tx_start+=1;  #adjust tx_start to 1-base
   if($prediction_id>24&&$prediction_id<=40){$gene_name=",$gene_id";$gene_id=$ccgenemap{"$prediction_id"}{"$gene_name"};}
   $attributes="gene_id=$gene_id$gene_name; transcript_id=$tx; exon_count=$exoncount;";$tstrand="";
   $featureType=($opt eq "lex")?"3'Exon":$featureType;$featureType=($opt eq "fex")?"5'Exon":$featureType;
   $featureType=($opt eq "ex")?"exon":$featureType;
   if(($opt eq "lex")||($opt eq "fex")||($opt eq "tex")){
     ($firstExon,$lastExon)=split(":",getLastTerminalExon($tx_id,$tx,$org_vid,$prediction_id,$opt));
     ($chrom,$fexstart,$fexend,$strand)=split(",",$firstExon);($chrom,$lexstart,$lexend,$strand)=split(",",$lastExon);
     if($opt eq "tex"){
       $rows="$chrom\t$annot\t3'Exon\t$lexstart\t$lexend\t.\t$strand\t.\t$attributes\n";
       $rows.="$chrom\t$annot\t5'Exon\t$fexstart\t$fexend\t.\t$strand\t.\t$attributes\n";
      }
      else{
       $rows="$chrom\t$annot\t$featureType\t$lexstart\t$lexend\t.\t$strand\t.\t$attributes\n" if($opt eq "lex");
       $rows="$chrom\t$annot\t$featureType\t$fexstart\t$fexend\t.\t$strand\t.\t$attributes\n" if($opt eq "fex");
      }
   }      
   else{
     #get the cds data for this transcript
     $cdsStart=0;$cdsEnd=0;
     $qh_getCDS->execute($tx_id,$tx,$org_vid,$prediction_id);($cdsStart,$cdsEnd)=$qh_getCDS->fetchrow_array();
     $cdsStart=($cdsStart>0)?$cdsStart+1:$cdsStart;
      while(($chrom,$exstart,$exend,$strand,$frame)=$qh_getTxExons->fetchrow_array()){
             $exstart+=1;$frame=($frame==-99)?".":$frame;$frame=($frame eq "")?".":$frame;
             $frame=($prediction_ids>24&&$prediction_ids<=40)?".":$frame;
             $tstrand =$strand if($tstrand eq "");
             if($exend>0){
                if($strand eq "+"){$rank=$prank;++$prank;}else{$rank=$nrank;--$nrank;}
                $rows.="$chrom\t$annot\texon\t$exstart\t$exend\t.\t$strand\t$frame\t$attributes exon_rank=$rank;\n";
                #CHECK CDS
                if(($cdsStart>0)&&($cdsEnd>0)){
                    if(($cdsStart>=$exstart)&&($exend>$cdsStart)){$exstart=$cdsStart;}
                    elsif(($exstart<$cdsEnd)&&($exend>=$cdsEnd)){$exend=$cdsEnd;}
                    $rows.="$chrom\t$annot\tCDS\t$exstart\t$exend\t.\t$strand\t$frame\t$attributes exon_rank=$rank;\n";
                 }
             }
        }
   }# end of  $qh_gettxexon->rows>0 
 return $rows;  
 }
sub getTranscriptNew{
 ($tx_id,$tx,$org_vid,$prediction_id)=@_;
 $qh_getTxcord->execute($tx_id,$org_vid);#get the coordinates of this transcript
 ($chrom_name,$txstrand,$tx_start,$tx_end)=$qh_getTxcord->fetchrow_array();
 $qh_getTxAll->execute($tx_id,$org_vid);%tx_map=();$prev_tx="";$exids="";
 while(($tx,$exid)=$qh_getTxAll->fetchrow_array()){
     if($prev_tx eq ""){$exids="$exid";$prev_tx=$tx;}
     elsif($prev_tx ne "$tx"){ #this is a new accession id
        if(exists($tx_map{"$exids"})){$tx_map{"$exids"}.=",$prev_tx";}
        else{$tx_map{"$exids"}="$prev_tx";}$exids="$exid";$prev_tx=$tx;
      }else{$exids.=",$exid";}
 }$row="";
 if(exists($tx_map{"$exids"})){$tx_map{"$exids"}.=",$prev_tx";}
 else{$tx_map{"$exids"}="$prev_tx";}
 while(($tx_exons,$txs)=each(%tx_map)){
    next if(!($tx_exons=~/\d+/));
    $cdsStarts="";$cdsEnds=""; $exonStarts="";$exonEnds="";
    @tx_count=split(",",$txs);$tx_string="";
    foreach $tx(@tx_count){if($tx_string eq ""){$tx_string="'$tx'";}else{$tx_string.=",'$tx'";}}
    $query="Select distinct cdsStart,cdsEnd from transcript_translation ";
    $query.=" where transcript_id=$tx_id and organism_version_id=$org_vid and transcript_name in ($tx_string) order by cdsStart asc";

    @exon_count=split(",",$tx_exons);@tx_count=split(",",$txs);
    $cdsStarts="";$cdsEnds="";$exonStarts="";$exonEnds="";$gene="";
    $qh_getCDS->execute($tx_id,$txs,$org_vid)or die $qh_getCDS->errstr;
    while(($cdsStart,$cdsEnd)=$qh_getCDS->fetchrow_array()){
      $cdsStart+1;$cdsStarts=($cdsStarts eq "")? "$cdsStart":"$cdsStarts,$cdsStart";
      $cdsEnds=($cdsEnds eq "")?"$cdsEnd":"$cdsEnds,$cdsEnd";
    }#get exon cordinates 
    $qh_getExonCord->execute($tx_exons);
    while(($exon_start,$exon_end)=$qh_getExonCord->fetchrow_array()){$exon_start+=1;
       if($exonStarts eq ""){$exonStarts="$exon_start";$exonEnds="$exon_end";}else{$exonStarts.=",$exon_start";$exonEnds.=",$exon_end";}
    }
    $qh_getTxGeneAll->execute($org_vid,$tx_id);
    while(($gen)=$qh_getTxGeneAll->fetchrow_array()){
       if($prediction_id>24&&$prediction_id<=40){
          if($ccgenemap{"$prediction_id"}{"$gen"}){$gen=$ccgenemap{"$prediction_id"}{"$gen"}.",$gen";}
       }
       if($gene eq ""){$gene=$gen;}else{$gene.=";$gen";}
    } #
    if((@tx_count>1)&&($qh_getCDS->rows>1)){
       $qh_getPredList->execute($tx_id,$txs); $annot="";
       while(($pred_name)=$qh_getPredList->fetchrow_array()){
           if($annot eq ""){$annot="$pred_name";}else{$annot.=",$pred_name";}}
     # print TS "name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tname2\tsource\n";
       $row="$txs\t$chrom_name\t$txstrand\t$tx_start\t$tx_end\t$cdsStarts\t$cdsEnds\t".@exon_count."\t$exonStarts\t$exonEnds";
       $row.="\t$gene\t$annot\n";
    }
  }
  return $row;
 }
sub getTranscript{
 ($tx_id,$tx,$org_vid,$prediction_id)=@_;
 $qh_getTxcord->execute($tx_id,$org_vid);#get the coordinates of this transcript
 ($chrom_name,$txstrand,$tx_start,$tx_end)=$qh_getTxcord->fetchrow_array();
 $gene="";$exonStarts="";$exonEnds="";$cdsStarts="0";$cdsEnds="0";$annot="";
 $qh_getTxExons->execute($tx_id,$tx,$org_vid,$prediction_id);$qh_getPredName->execute($prediction_id);
 ($annot)=$qh_getPredName->fetchrow_array();$qh_getCDS->execute($tx_id,$tx,$org_vid,$prediction_id);
 $exoncount=$qh_getTxExons->rows; 
 while(($chromosome_name,$exstart,$exend,$strand)=$qh_getTxExons->fetchrow_array()){
        $exstart+=1; #adjust exstart to 1-base
        if($exonStarts eq ""){$exonStarts=$exstart;$exonEnds=$exend;
        }else{$exonStarts.=",$exstart";$exonEnds.=",$exend";}
 }
 $cdsStarts=0;$cdsEnds=0;
 $qh_getCDS->execute($tx_id,$tx,$org_vid,$prediction_id);($cdsStart,$cdsEnd)=$qh_getCDS->fetchrow_array();
 $cdsStart=($cdsStart>0)?$cdsStart+1:$cdsStarts;
 $cdsEnd=($cdsEnd>0)?$cdsEnd:$cdsEnds;
 $qh_getTxGene->execute($org_vid,$tx_id,$prediction_id);
 $tx_start+=1;  #adjust tx_start to 1-base
 while(($gen)=$qh_getTxGene->fetchrow_array()){
       if($prediction_id>24&&$prediction_id<=40){
          if($ccgenemap{"$prediction_id"}{"$gen"}){$gen=$ccgenemap{"$prediction_id"}{"$gen"}.",$gen";}
       }
       if($gene eq ""){$gene=$gen;}else{$gene.=";$gen";}
  } #
  $row="$tx\t$chrom_name\t$txstrand\t$tx_start\t$tx_end\t$cdsStart\t$cdsEnd\t$exoncount\t$exonStarts\t$exonEnds\t$gene\t$annot\n";
  return $row;
 }
############################
$organism_version="";$organism_version=$opt_v if($opt_v);$org_vid=0;
if($organism_version eq ""&& !($opt_s)&&!($opt_g)){
  print "Must specify the organism version. For example -v mm9 .Please run the script with the -h option for usage\n";
  exit(0);
}$organism_version=~s/\s+//g;
if($opt_s){# distinct o.organism_name,ov.organism_version_id ,ucsc_db as version
    $type="s=1";
    print "Organism\tVersion\n";displayOrganismVersions($type);
}
elsif($opt_l){ #display list of annotations for this organism version
   print "Annotation sources for $organism_version\n"; displayAnnotSources($organism_version);
}
elsif($opt_g){ #display list of organism versions with genome data
    $type="g=1";
    print "Organism\tVersion\n";displayOrganismVersions($type);
}
else{
 $qh_getOrgV->execute("$organism_version");if($qh_getOrgV->rows>0){($org_vid)=$qh_getOrgV->fetchrow_array();}
 $type=(!$opt_t)?"tx":$opt_t;
 $gene_prediction="";$pred_id=0;%predictions=();$prediction_ids="";
 if($opt_a){$gene_prediction=$opt_a;@predictions=split(",",$gene_prediction);
   foreach my $pred(@predictions){ #generate gene prediction ids for selected gene predisctions
       $qh_getPredId->execute($pred);if($qh_getPredId->rows>0){
       ($pred_id)=$qh_getPredId->fetchrow_array();if($prediction_ids eq ""){$prediction_ids="$pred_id";}
        else{$prediction_ids.=",$pred_id";}
       }
   }
 } # default to all annotations
 #set the result path
 $file_name="";$type=~s/\/\s*$//g;$feature_type="Transcript";
 if($opt_f){$file_name=$opt_f;}
 elsif($type=~/tx/){
    $file_name="$organism_version-$gene_prediction-transcripts.txt";
    $file_name="$organism_version-$gene_prediction-transcripts.gtf" if($opt_F);
  }else{$file_name="$organism_version-$gene_prediction-features.gtf";}

 if($opt_o){$opt_o=~s/\/\s*$//;$file_name="$opt_o/$file_name";}
 open (TS,">$file_name") or die "$!\n"; $count=0;my $tx_qh;
 if($opt_T){ ## a hit on gene_by annotation table #check if genomic region or mgi id
    $term=$opt_T;
    if($term=~/mgi:/){$qh_getMgiHits->execute($opt_T,$org_vid);$tx_qh=$qh_getMgiHits;} #mgi id?
    else{
       ($chromosome_id,$chrom_start,$chrom_end)=split(":",getGenomicRegion($term));
       if(($chromosome_id>0)&&($chrom_start>0)&&($chrom_end>0)){ #genome region query?
           $tx_qh=$qh_getCordAllHits;$tx_qh->execute($org_vid,$chromosome_id,$chrom_start,$chrom_end);
       }else{## a hit on gene_by annotation table
         $qh_getGenTxHits->execute($opt_T,$org_vid);$tx_qh=$qh_getGenTxHits;
         if($tx_qh->rows<=0){ ## a hit on transcript_by annotation table
            $qh_getTxHits->execute($opt_T,$org_vid);$tx_qh=$qh_getTxHits;
            if($tx_qh->rows<=0){ ## a hit on cc_founders_transcript table (transcript_name)
               $qh_getCcTxHits->execute($opt_T,$org_vid);$tx_qh=$qh_getCcTxHits;
               if($tx_qh->rows<=0){ ## a hit on cc_founders_genes table (gene_id)
                  $qh_getCcGenHits->execute($opt_T,$org_vid);$tx_qh=$qh_getCcGenHits;
               }
            }
          }
       }
    }
 }else{
   if(($prediction_ids ne "")||($prediction_ids>0)){$qh_getTxByAnnot->execute($org_vid,$prediction_ids);$tx_qh=$qh_getTxByAnnot;}
   else{$qh_getTxByAllAnnot->execute($org_vid);$tx_qh=$qh_getTxByAllAnnot;}
 }## load extra gene info if CC founders 
 if($type=~/tx/){
   if(!$opt_F){ print TS "name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tname2\tsource\n";
   }else{print TS "chrom\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n";}
 }else{print TS "chrom\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n";}
 if($tx_qh->rows>0){$count=0; %txmap=();
    while(($tx_id,$tx,$prediction_id)=$tx_qh->fetchrow_array()){
        if($prediction_ids ne ""){if($prediction_ids=~/^\s*(\d+)\s*$/){next if($1!=$prediction_id);}
           next if(!($prediction_ids=~/$prediction_id/));
         }
        if($prediction_id>24&&$prediction_id<=40){
          if(!exists($ccgenemap{"$prediction_id"})){$qh_getCCgene->execute($prediction_id);
             while(($gene_prediction_id, $gene_id, $gene_name )=$qh_getCCgene->fetchrow_array()){ 
               $ccgenemap{"$gene_prediction_id"}{"$gene_name"}=$gene_id;}
          }
        }
        if((!$opt_F)&&($type=~/tx/)){
            if(!exists($txmap{"$tx_id"})){
               $row=getTranscript($tx_id,$tx,$org_vid,$prediction_id);print TS "$row"if($row ne "");
               $txmap{"$tx_id"}=1;
            }
        }
        else{$row=getFeatures($tx_id,$tx,$org_vid,$prediction_id,$type);print TS "$row"if($row ne "");}
    }
  }
 print "$file_name generated\n";
}#end of query by -t option
print"program completed\n";
exit(0);

