#!/usr/bin/perl

###################################################################################################
## getCordinates.pl -- current version
#This script generates the coordinates  of transcripts/features for the specified organism version and annotation source.
#  If an anotation term (gene name, transcript name , genomic region) is spefified 
#  for a given organism version, the program returns all features/transcripts associated with
#  the specified term .

#  The features are stored into the database using zero-based cordinates standard for feature starts
#  and 1-based cordinate standard for feature ends
#    
# Note: The output of this program can be used by other programs, for example snpAnnotator, getTxFastaSeq 
# Outputs: either a genePred format (transcript query) 
#         or a gtf format (feature queries, ex, int,...)

## Features coordinates in gtf file are generated in 1-based start standard and 1-base <end>
## Features coordinates in genePrediction file are generated in 0-based start standard and 1-base <end>

# 1.a gtf file has the following fields
# GTF (Gene Transfer Format) is a refinement to GFF that tightens the specification. The first eight GTF fields are the same as GFF. The group field has been expanded into a list of attributes. Each attribute consists of a type/value pair. Attributes must end in a semi-colon, and be separated from any following attribute by exactly one space.
# The attribute list must begin with the two mandatory attributes:
# gene_id value -       A globally unique identifier for the genomic source of the sequence.
# transcript_id value - A globally unique identifier for the predicted transcript. 
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
# <seqname> is the chromosome
# <source> is the gene prediction source name
# <feature> is feature type (exon, CDS)
# [attributes] a list of attributes (gene_id ""; transcript_id ""; gene_name ""; feature_id ""; feature_rank "";...)
# 
# 2. a genePred file has the following fields
# name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tid\tname2\texonFrames\n
#
#
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : August 2013
###################################################################################################
use DBI;
use POSIX;
use XML::Simple;

# create object
$xml = XML::Simple->new (ForceArray => 1);

use vars qw ($opt_h $opt_s $opt_o $opt_t $opt_v $opt_a $opt_l $opt_f $opt_d $opt_F $opt_T $opt_g,$opt_p);
use Getopt::Std;

getopts('hlsgo:t:v:a:f:d:FT:p:');
if($opt_h||(!$opt_v&&!$opt_s&&!$opt_l&&!$opt_g)) {
    print <<HELP;

This script generates the coordinates  of transcripts/features for the specified organism version and annotation source.
The features are stored into the database using the zero-based cordinates standard for feature starts
and 1-based cordinate standard for feature ends

## Features coordinates in gtf file are generated in 1-based <start> standard and 1-base <end>
## Features coordinates in genePrediction file are generated in 0-based <start> standard and 1-base <end>
    
Note: The output of this program can be used by other programs, for example snpAnnotator, getTxFastaSeq 
      result coordinates are  0-based cordinates standard for feature starts 
      and 1-based cordinate standard for feature ends

Output: either a genePred format (transcript query) or a gtf format (feature queries)

Usage: perl getCordinates.pl -t <feature type> -v <organism version> [-T pax6] [-a <gene prediction list>] 
           [-o <outputDir>][-f outfilename] [-F <result format>]
Examples:

Example: perl getCordinates.pl -v mm9 -l
The above will display all the gene prediction sources for mm9 currently in graber_transcript database

Example: perl getCordinates.pl -s
The above will display all the organism versions currently in graber_transcript database

Example: perl getCordinates.pl -o features -v mm10 -a GRCm38-ensGene
  The above will use the current Ensembl annotations for mm10 to generate all the transcripts
  in gene prediction format and store the file (mm10-GRCm38-ensGene-transcripts.txt) under features/

Example: perl getCordinates.pl -o features -v mm10 -a GRCm38-ensGene -F 
The above will use current Ensembl annotations for mm10 to generate all the transcripts
in gtf format and store the file (mm10-GRCm38-ensGene-transcripts.gtf) under features/

Example: perl getCordinates.pl -v mm10 -a GRCm38-ensGene  -f mm10-ensemblTransctipts.gff -F -p chr
The above will use Ensembl annotations for mm10 to generate all the transcripts
in gtf format and store the file (mm10-ensemblTransctipts.gtf). The chromosome filed will have the "chr" prefix

Example: perl getCordinates.pl -v mm10 -T pax6 -f pax6-transcripts.txt
  The above will use pax6 annotations from all the annotation sources for mm10 and generate transcripts data
  in gene prediction format and store the file (pax6-transcripts.txt) under the current working directory

Example: perl getCordinates.pl -v mm10 -T pax6 -f pax6-transcripts.gtf -F
  The above will use pax6 annotations from all the annotation sources for mm10 and generate transcripts data
  in gtf format and store the file (pax6-transcripts.gtf) under the current working directory


Arguments:
   -h  displays this help message
  
   -t  feature type: 
       tx     -> Transcripts - default (gene prediction format or gtf format)
       ex     -> Exons      (gtf format)
       fex    -> First Exons (gtf format)
       lex    -> last exons  (gtf format)
       tex    -> terminal exons (5' and 3' exons) (gtf format)
  -T  query term (gene name, transcript, genomic region[format chrom:start-end]) 
  -v  is the ucsc organism version (default mm9)
  -a  optional, a commas separated list of gene prediction names to use
  -o  optional, is the output directory, default current working directory
  -f  optional, output file name (relative to the specified output directory)
  -F  optional, output format, default genePred format (only used for the -t tx )
  -l  display the list of current gene predictions for the specified organism version
  -s  display all the organism versions currently in graber_transcript database
  -g  display all the organism versions that have genome sequence (chromosome)
  -p  append prefix to the chromosome field
 
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
$qh_getChr=$dbh->prepare("select chromosome_id from chromosome where chromosome_name=?");
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$qh_getPredName=$dbh->prepare("select gene_prediction_name from gene_prediction where gene_prediction_id=?");

$query="select distinct gene_prediction_name,g.gene_prediction_id from gene_prediction g,transcript_exon t";
$query.=" where  t.organism_version_id=? and t.transcript_id=? and transcript_name=? ";
$query.=" and t.gene_prediction_id=g.gene_prediction_id";
$qh_getTxPredId=$dbh->prepare($query);

$query="select distinct transcript_id,transcript_name,gene_prediction_id from transcript_by_annotation
       where organism_version_id=? and gene_prediction_id=?";
$qh_getTxByAnnot=$dbh->prepare($query);
$query="select distinct transcript_id,transcript_name,gene_prediction_id from transcript_by_annotation 
       where organism_version_id=? and gene_prediction_id!=3";
$qh_getTxByAllAnnot=$dbh->prepare($query);

my $getTxcord="select chromosome_name,strand,tx_start,tx_end from transcript e, chromosome c";
   $getTxcord.=" where transcript_id=? and organism_version_id=? and e.chromosome_id=c.chromosome_id";
   $qh_getTxcord = $dbh->prepare($getTxcord)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTxExons="select distinct chromosome_name,exon_start,exon_end,strand,exon_frame 
                from transcript_exon t,exon e,chromosome c ";
   $getTxExons.=" where transcript_id=? and transcript_name=? and t.organism_version_id=? and gene_prediction_id in (?)";
   $getTxExons.=" and t.organism_version_id=e.organism_version_id and
                  t.exon_id=e.exon_id and e.chromosome_id=c.chromosome_id";
   $getTxExons.=" order by exon_start asc";
my $qh_getTxExons = $dbh->prepare($getTxExons)or die "Couldn't prepare statement: ".$dbh->errstr;

my $qh_getExCount=$dbh->prepare("select count(distinct exon_id) from transcript_exon where transcript_id=? and transcript_name=? and organism_version_id=? and gene_prediction_id in (?)");

my $qh_getCCgene=$dbh->prepare("select distinct gene_prediction_id, gene_id, gene_name from cc_founders_genes where gene_prediction_id in (?)");
my $qh_getCCtx=$dbh->prepare("select distinct transcript_id,transcript_name,protein_id from cc_founders_transcripts where gene_prediction_id=?");

my $getTxExons="select distinct e.exon_id,exon_start,exon_end,strand from transcript_exon t,exon e ";
   $getTxExons.=" where  transcript_id=? and transcript_name=? and t.organism_version_id=? and gene_prediction_id!=3";
   $getTxExons.=" and t.organism_version_id=e.organism_version_id and t.exon_id=e.exon_id ";
   $getTxExons.=" order by exon_start asc";
my $qh_getTxExonsAll = $dbh->prepare($getTxExons)or die "Couldn't prepare statement: ".$dbh->errstr;
my $qh_getExCountAll=$dbh->prepare("select count(distinct exon_id) from transcript_exon where transcript_id=? and transcript_name=? and organism_version_id=? and gene_prediction_id !=3");

$query="select cdsStart,cdsEnd from transcript_translation ";
$query.=" where transcript_id in (?) and transcript_name in(?) and organism_version_id=? and gene_prediction_id in (?) order by cdsStart asc";
my $qh_getCDS=$dbh->prepare($query);

my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $qh_getOrgV = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getOrgUcsc="select ucsc_db  from organism_version where organism_version_id=?";
my $qh_getOrgUcsc = $dbh->prepare($getOrgUcsc)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getTxGene="select distinct gene_name from gene_by_annotation g, transcript_by_annotation t
                where g.organism_version_id=? and t.transcript_name=? and t.gene_prediction_id in (?)
                and t.transcript_id=g.transcript_id and t.gene_prediction_id=t.gene_prediction_id ";
my $qh_getTxGene = $dbh->prepare($getTxGene)or die "Couldn't prepare statement: ".$dbh->errstr;
my $qh_getTxGeneAll = $dbh->prepare("select distinct gene_name from gene_by_annotation where organism_version_id=? and transcript_id=?");
###### this section should be a function call of the webservice ##########################################
my $getAllOrgV="select distinct o.organism_name,ucsc_db as version from organism o,organism_version ov
                 where o.organism_id=ov.organism_id order by organism_name";
my $qh_getAllOrgV = $dbh->prepare($getAllOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;
my $getGenePred="select gene_prediction_name from gene_prediction_by_organism op, gene_prediction p
                   where op.organism_version_id =? and op.gene_prediction_id=p.gene_prediction_id 
                   and gene_prediction_name not in ('all_est','all_mrna','estOrientInfo')";
my $qh_getGenePred = $dbh->prepare($getGenePred)or die "Couldn't prepare statement: ".$dbh->errstr;
$query="select distinct sample_name from rnaseq_sample s, rnaseq_transcript t ";
$query.=" where t.transcript_name=? and gene_prediction_id=? and organism_version_id=? ";
$query.=" and t.sample_id=s.sample_id";
$qh_getSampleName=$dbh->prepare($query);
$query="select distinct g2.gene_name from gene_by_annotation g, gene_by_annotation g2
        where g.gene_name=? and g.organism_version_id=? and g.gene_prediction_id=? and
              g.organism_version_id=g2.organism_version_id and g.gene_prediction_id!=g2.gene_prediction_id 
              and g.transcript_id=g2.transcript_id  and g.gene_prediction_id=41
              and g2.gene_prediction_id=22 ";

$qh_getmgi=$dbh->prepare($query); %ccgenemap=(); %txDup=();
sub displayAnnotSources{
 ($org_version)=@_;$url="$webservice?v=$org_version&l=1";
  $xml_file="annot_list.xml";if(-e "annot_list.xml"){system("rm $xml_file");}
  system("wget -q -O $xml_file '$url'");
  if(-e $xml_file){$data =$xml->XMLin("$xml_file");
    foreach my $source(@{$data->{prediction}}){print $source->{name}[0] . "\n";}
  }if(-e "annot_list.xml"){system("rm $xml_file");}
}
sub displayOrganismVersions{
  ($type)=@_;
  $url="$webservice?$type";$xml_file="annot_list.xml";if(-e "annot_list.xml"){system("rm $xml_file");}
  system("wget -q -O $xml_file '$url'");
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
        $tstrand =$strand if($tstrand eq "");#$exstart+=1;
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
   ($prediction_ids,$tx_qh,$opt,$org_vid,$outfileh,$chrom_prefix)=@_;
  open(my $outfileh,">$file_name") or die "$!\n";
  %dups=();
  while(($tx_id,$tx,$prediction_id)=$tx_qh->fetchrow_array()){
       if($prediction_ids ne ""){
            if($prediction_ids=~/^\s*(\d+)\s*$/){next if($1!=$prediction_id);}next if(!($prediction_ids=~/$prediction_id/));
        }if($prediction_id>24&&$prediction_id<=40){
            if(!exists($ccgenemap{"$prediction_id"})){$qh_getCCgene->execute($prediction_id);
               while(($gene_prediction_id, $gene_id, $gene_name )=$qh_getCCgene->fetchrow_array()){ 
                     $ccgenemap{"$gene_prediction_id"}{"$gene_id"}=$gene_name;
        }}}
  
   $prank=1; $nrank=$exoncount;$rank=0;$rows=""; $featureType="";$gene_id="";$gene_name="";$strand="";
   $annot=""; $qh_getTxExons->execute($tx_id,$tx,$org_vid,$prediction_id);$exoncount=$qh_getTxExons->rows;
   next if($qh_getTxExons->rows<=0);
   $qh_getPredName->execute($prediction_id);($annot)=$qh_getPredName->fetchrow_array();
   $qh_getTxGene->execute($org_vid,$tx,$prediction_id);
   while(($gene_name)=$qh_getTxGene->fetchrow_array()){
          $gene_id=$gene_name if(($gene_name=~/^ENS/)||($gene_name=~/^rsd_g/));
    }
   $gene_id=($gene_id eq "")?$gene_name:$gene_id;
   if($gene_id=~/^ENS|rsd_g/){
     $qh_getmgi->execute($gene_id,$org_vid,$prediction_id);
     if($qh_getmgi->rows>0){($gene_name)=$qh_getmgi->fetchrow_array();}
   }
   if($prediction_id>24&&$prediction_id<=40){
     $gene_name=$ccgenemap{"$prediction_id"}{"$gene_id"} if($gene_name eq "");}
  
   $attributes="gene_id \"$gene_id\";";
   $attributes.=" gene_name \"$gene_name\";" if(($gene_name ne "$gene_id")&&($gene_name ne ""));
   $attributes.=" transcript_id \"$tx\"; exon_count \"$exoncount\";";
   next if(defined $dups{"$attributes"});
   $dups{"$attributes"}=1;
 #  print $outfileh "$attributes\n"; next;
   $tstrand="";$featureType=($opt eq "lex")?"3'Exon":$featureType;
   $featureType=($opt eq "fex")?"5'Exon":$featureType;$featureType=($opt eq "ex")?"exon":$featureType;$sample_id="";
   #get sample ids to fill the source field
   if($tx=~/^rsd_/){ $qh_getSampleName->execute($tx,$prediction_id,$org_vid);
       while(($sample)=$qh_getSampleName->fetchrow_array()){
          if($sample_id eq ""){$sample_id=$sample;}else{$sample_id.=",$sample";}}
   }$annot=($sample_id eq "")?$annot:$sample_id;$rows="";
   if(($opt eq "lex")||($opt eq "fex")||($opt eq "tex")){
     ($firstExon,$lastExon)=split(":",getLastTerminalExon($tx_id,$tx,$org_vid,$prediction_id,$opt));
     ($chrom,$fexstart,$fexend,$strand)=split(",",$firstExon);($chrom,$lexstart,$lexend,$strand)=split(",",$lastExon);
      $lexstart+=1;$fexstart+=1;
     $chrom=($chrom_prefix ne "")?"$chrom_prefix$chrom":$chrom;
     if($opt eq "tex"){
       $rows="$chrom\t$annot\t3'Exon\t$lexstart\t$lexend\t.\t$strand\t.\t$attributes\n";
       $rows.="$chrom\t$annot\t5'Exon\t$fexstart\t$fexend\t.\t$strand\t.\t$attributes\n";
      }else{
       $rows="$chrom\t$annot\t$featureType\t$lexstart\t$lexend\t.\t$strand\t.\t$attributes\n" if($opt eq "lex");
       $rows="$chrom\t$annot\t$featureType\t$fexstart\t$fexend\t.\t$strand\t.\t$attributes\n" if($opt eq "fex");
      }
     print $outfileh "$rows" if($rows ne "");
   }else{ #get the cds data for this transcript
      $cdsStart=0;$cdsEnd=0;$qh_getCDS->execute($tx_id,$tx,$org_vid,$prediction_id);
     ($cdsStart,$cdsEnd)=$qh_getCDS->fetchrow_array(); $cdsStart=($cdsStart>0)?$cdsStart:$cdsStart;
      $prank=1;$nrank=$exoncount;#get the coordinates of this transcript
      $qh_getTxcord->execute($tx_id,$org_vid);($chrom_name,$txstrand,$tx_start,$tx_end)=$qh_getTxcord->fetchrow_array();
      $tx_start+=1;$rows="";$chrom_name=($chrom_prefix ne "")?"$chrom_prefix$chrom_name":$chrom_name;
      $rows.="$chrom_name\t$annot\ttranscript\t$tx_start\t$tx_end\t.\t$txstrand\t.\t$attributes;\n";
      while(($chrom,$exstart,$exend,$strand,$frame)=$qh_getTxExons->fetchrow_array()){
         $frame=($frame==-99)?".":$frame;$frame=($frame eq "")?".":$frame;
         $frame=($prediction_ids>24&&$prediction_ids<=40)?".":$frame;$tstrand =$strand if($tstrand eq "");
         $chrom=($chrom_prefix ne "")?"$chrom_prefix$chrom":$chrom;
         if($exend>0){
            if($strand eq "+"){$rank=$prank;++$prank;}else{$rank=$nrank;--$nrank;}
            $exstart+=1;
            $rows.="$chrom\t$annot\texon\t$exstart\t$exend\t.\t$strand\t$frame\t$attributes exon_rank \"$rank\";\n";
            if(($cdsStart>0)&&($cdsEnd>0)){
                if(($cdsStart>=$exstart)&&($exend>$cdsStart)){$exstart=$cdsStart;}
                elsif(($exstart<$cdsEnd)&&($exend>=$cdsEnd)){$exend=$cdsEnd;}
                $cdsStart+=1;
                $rows.="$chrom\t$annot\tCDS\t$exstart\t$exend\t.\t$strand\t$frame\t$attributes exon_rank \"$rank\";\n";
          }}
         #print $outfileh "$rows" if($rows ne "");
       }
      print $outfileh "$rows" if($rows ne "");
   }
  } 
}

sub getTranscript{
 ($prediction_ids,$tx_qh,$org_vid,$chrom_prefix)=@_;
 $count=0;
 while(($tx_id,$tx,$prediction_id)=$tx_qh->fetchrow_array()){
       if($prediction_ids ne ""){
            if($prediction_ids=~/^\s*(\d+)\s*$/){next if($1!=$prediction_id);}
            next if(!($prediction_ids=~/$prediction_id/));
        }
        if($prediction_id>24&&$prediction_id<=40){
           if(!exists($ccgenemap{"$prediction_id"})){$qh_getCCgene->execute($prediction_id);
               while(($gene_prediction_id, $gene_id, $gene_name )=$qh_getCCgene->fetchrow_array()){ 
                      $ccgenemap{"$gene_prediction_id"}{"$gene_id"}=$gene_name;}
            }
        }
   $qh_getTxcord->execute($tx_id,$org_vid);#get the coordinates of this transcript
   ($chrom_name,$txstrand,$tx_start,$tx_end)=$qh_getTxcord->fetchrow_array();
   $gene="";$exonStarts="";$exonEnds="";$cdsStarts="0";$cdsEnds="0";$annot="";
   $qh_getTxExons->execute($tx_id,$tx,$org_vid,$prediction_id);$qh_getPredName->execute($prediction_id);
   ($annot)=$qh_getPredName->fetchrow_array();$sample_id="";
   if($tx=~/^rsd_/){ $qh_getSampleName->execute($tx,$prediction_id,$org_vid);
       while(($sample)=$qh_getSampleName->fetchrow_array()){
          if($sample_id eq ""){$sample_id=$sample;}else{$sample_id.=",$sample";}}
   }$annot=($sample_id eq "")?$annot:$sample_id;$rows="";
   $exoncount=$qh_getTxExons->rows; 
   while(($chromosome_name,$exstart,$exend,$strand)=$qh_getTxExons->fetchrow_array()){
        #$exstart+=1; #adjust exstart to 1-base
        if($exonStarts eq ""){$exonStarts=$exstart;$exonEnds=$exend;
        }else{$exonStarts.=",$exstart";$exonEnds.=",$exend";}
   }$cdsStarts=0;$cdsEnds=0;
   $qh_getCDS->execute($tx_id,$tx,$org_vid,$prediction_id);($cdsStart,$cdsEnd)=$qh_getCDS->fetchrow_array();
   #$cdsStarts+=1; #adjust cdsStarts to 1-base
   $cdsStart=($cdsStart>0)?$cdsStart:$cdsStarts;$cdsEnd=($cdsEnd>0)?$cdsEnd:$cdsEnds;
   $qh_getTxGene->execute($org_vid,$tx,$prediction_id);$gene_name="";$gene_id="";
   $row="";
   while(($gene_name)=$qh_getTxGene->fetchrow_array()){
          $gene_id=$gene_name if(($gene_name=~/^ENS/)||($gene_name=~/^rsd_g/));
    }$gene_id=($gene_id eq "")?$gene_name:$gene_id;
    if($prediction_id>24&&$prediction_id<=40){$gene_name=$ccgenemap{"$prediction_id"}{"$gene_id"} if($gene_name eq "");}
   $gene_name=$gene_id if($gene_name eq "");
   ++$count; print "$count processed\n" if($count%10000==0);
   $chrom_name=($chrom_prefix ne "")?"$chrom_prefix$chrom_name":$chrom_name;
   $txDup{"$chrom_name\t$txstrand\t$tx_start\t$tx_end"}{"$exoncount\t$exonStarts\t$exonEnds"}.="$tx=$gene_name=$annot=$cdsStart=$cdsEnd;";
 }

}
############################
$organism_version="";$organism_version=$opt_v if($opt_v);$org_vid=0;
if($organism_version eq ""&& !($opt_s)&&!($opt_g)){
  print "Must specify the organism version. For example -v mm9 .Please run the script with the -h option for usage\n";
  exit(0);
}$organism_version=~s/\s+//g;
if($opt_s){# distinct o.organism_name,ov.organism_version_id ,ucsc_db as version
    $type="s=1";print "Organism\tVersion\n";displayOrganismVersions($type);
}elsif($opt_l){ #display list of annotations for this organism version
   print "Annotation sources for $organism_version\n"; displayAnnotSources($organism_version);
}elsif($opt_g){ #display list of organism versions with genome data
    $type="g=1";
    print "Organism\tVersion\n";displayOrganismVersions($type);
}else{
 $qh_getOrgV->execute("$organism_version");if($qh_getOrgV->rows>0){($org_vid)=$qh_getOrgV->fetchrow_array();}
 $type=(!$opt_t)?"tx":$opt_t;
 $gene_prediction="";$pred_id=0;%predictions=();$prediction_ids="";
 if($opt_a){$gene_prediction=$opt_a;@predictions=split(",",$gene_prediction);
   foreach my $pred(@predictions){ #generate gene prediction ids for selected gene predisctions
       $qh_getPredId->execute($pred);if($qh_getPredId->rows>0){
       ($pred_id)=$qh_getPredId->fetchrow_array();
        if($prediction_ids eq ""){$prediction_ids="$pred_id";}
        else{$prediction_ids.=",$pred_id";}
       }
   }
  } # default to all annotations
 #set the result path
 $file_name="";$type=~s/\/\s*$//g;$feature_type="Transcript";
 $chrom_prefix=($opt_p)?$opt_p:"";
 if($opt_f){$file_name=$opt_f;}
 elsif($type=~/tx/){
    $file_name="$organism_version-$gene_prediction-transcripts.txt";
    $file_name="$organism_version-$gene_prediction-transcripts.gtf" if($opt_F);
  }else{$file_name="$organism_version-$gene_prediction-features.gtf";}
 if($opt_o){$opt_o=~s/\/\s*$//;$file_name="$opt_o/$file_name";}
 $count=0;my $tx_qh;
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
 open(TS,">$file_name") or die "$!\n"; 
 if($type=~/tx/){ #gtf file does not have a header
   if(!$opt_F){ 
        print TS "name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tname2\tsource\n";
   } #else{print TS "chrom\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n";}
 }# else{print TS "chrom\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n";}
 if($tx_qh->rows>0){$count=0;
   if((!$opt_F)&&($type=~/tx/)){ 
       print "we have a total of ".$tx_qh->rows." transcripts\n";
       getTranscript($prediction_ids,$tx_qh,$org_vid,$chrom_prefix);
       while(($key,$value)=each(%txDup)){ #this gets rid of duplicate rows
          while(($cord,$feature)=each(%{$value})){ 
               @features=split(";",$feature);%tranc=();%genes=();%annots=();%cdsStarts=();%cdsEnds=();
               foreach $item (@features){
                  ($tx,$gene,$annot,$cdsStart,$cdsEnd)=split("=",$item);
                   $tranc{"$tx"}=1;$genes{"$gene"}=1 if($gene ne "");$annots{"$annot"}=1;
                  if($cdsStart>0){$cdsStarts{"$cdsStart"}=1;$cdsEnds{"$cdsEnd"}=1;}
               }
              $tx=join(",",keys %tranc);$gene=join(",",keys %genes);$annot=join(",",keys %annots);
             $cdsStart=join(",",keys %cdsStarts);$cdsEnd=join(",",keys %cdsEnds);
             if(keys(%cdsStarts)==0){$cdsStart=0;$cdsEnd=0;}
             print TS "$tx\t$key\t$cdsStart\t$cdsEnd\t$cord\t$gene\t$annot\n";
            #print htis transcript
         }
      }
    }
   else{
     print "we have a total of ".$tx_qh->rows." transcripts\n";
     close(TS);
     getFeatures($prediction_ids,$tx_qh,$type,$org_vid,$file_name,$chrom_prefix);
    }
 }
 #now generate file
 print "$file_name generated\n";
}#end of query by -t option
print"program completed\n";
exit(0);

