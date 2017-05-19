#!/usr/bin/perl

###################################################################################################
## getGtfCordinates.pl
#  This script generates feature cordinates of the specified term (gene, transcript, genomic region)
#  for a given organism version and gene prediction source. The default gene prediction is mm9-ensGene.
#  The features are stored into the database using zero-based cordinates standard for feature starts
#  and 1-based cordinate standard for feature ends
#    
#  
# Output: a gtf file with the following fields
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
# <seqname> is the chromosome
# <source> is the gene prediction source name
# <feature> is feature type (exon, CDS)
# [attributes] a list of attributes (gene_id;transcript_id;gene_name;feature_id; feature_rank;...)

#<start>, <end>
#Integers. <start> must be less than or equal to <end>. Sequence numbering starts at 1, 
#         so these numbers should be between 1 and the length of the relevant sequence, inclusive.
## Features coordinates are generated in 1-based cordinates standard for <start>, <end>

#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : March 2013

#   Usage: perl getGtfCordinates.pl -o <outputFile> -v <organismVersion> -a <gene prediction list> -t <query>
#
#   Where: 
#          -o  is the output file name (default standard out)            
#          -v  is the ucsc organism version (default mm9)
#          -a  a commas separated list of gene prediction names (default knowGene,ensGene)
#          -t  query term a commas separated list of query terms      
#          -s  displays the list of current organism versions 
#          -l  display the list of current gene predictions for the specified organism version

###################################################################################################
use DBI;
use POSIX;

use vars qw ($opt_h $opt_s $opt_o $opt_t $opt_v $opt_a $opt_l);
use Getopt::Std;

getopts('hlso:t:v:a:');
if($opt_h||(!$opt_o&&!$opt_s&&!$opt_l)) {
    print <<HELP;

#  This script generates feature cordinates of the specified term (gene, transcript, genomic region)
#  for a given organism version and gene prediction source. The default organism version is mm9
#  and the default gene prediction is mm9-ensGene.
#  The features are stored into the database using zero-based cordinates standard for feature starts
#  and 1-based cordinate standard for feature ends
#    
#  
# Output: a gtf file with the following fields
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes]
# <seqname> is the chromosome
# <source> is the gene prediction source name
# <feature> is feature type (exon, CDS)
# [attributes] a list of attributes (gene_id;transcript_id;gene_name;feature_id; feature_rank;...)
## Features coordinates are generated in 1-based cordinates standard for <start>, <end>

#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : March 2013

#   Usage: perl getGtfCordinates.pl -o <outputFile> -v <organismVersion> -a <gene prediction list> -t <query>
#
#   Where: 
#          -o  is the output file name (default standard out)            
#          -v  is the ucsc organism version (default mm9)
#          -a  a commas separated list of gene prediction names (default knowGene,ensGene)
#          -t  query term a commas separated list of query terms       
#          -s  display the list of current organism versions
#          -l  display the list of current gene predictions for the specified organism version

Example: perl getGtfCordinates.pl -v mm9 -l
The above will display all the gene predictions for mm9

Example: perl getGtfCordinates.pl -t Pax6 -o features -v mm9 -a mm9-ensGene,knownGene 
The above will use Ensembl  and knownGene annotations to generate all the transcripts
for mm9 and store the result file under features/. 

HELP
exit;
}

#set defaults
my $user ='pup';my $pwd  ='puppass';$dbname ='graber_transcriptdb';$host ='harlequin.jax.org';
$flanking=200; $organism_version="mm9";$mb=1000000;
if($opt_d){$host=$opt_d;}
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
##################################################
# get hits for the search term
##################################################
$query="select distinct transcript_id from gene_by_annotation where organism_version_id=?  and gene_name =?";
$qh_getGenTxHits=$dbh->prepare($query);
$query="select distinct transcript_id from transcript_by_annotation where organism_version_id=?  and transcript_name =?";
$qh_getTxHits=$dbh->prepare($query);
$query="select distinct transcript_id from gene_by_annotation g, mgi_genes m ";
$query .=" where g.organism_version_id=? and m.mgi_geneid =? and m.mgi_symbol=g.gene_name ";
$qh_getMgiHits=$dbh->prepare($query);
$query="select distinct transcript_id from cc_founders_transcripts where gene_prediction_id in (?) and transcript_name =?";
$qh_getCcTxHits=$dbh->prepare($query);
$query="select distinct transcript_id from gene_by_annotation g,cc_founders_genes c 
        where g.gene_prediction_id in (?) and g.gene_prediction_id=c.gene_prediction_id and c.gene_id =? and c.gene_name=g.gene_name";
$qh_getCcGenHits=$dbh->prepare($query);
$query=" select transcript_id from transcript where organism_version_id=? ";
$query.=" and chromosome_id=?  and strand=(?) and tx_end >=? and tx_start<=? ";
$qh_getCordHits= $dbh->prepare($query);
$query=" select transcript_id from transcript where organism_version_id=? ";
$query.=" and chromosome_id=? and tx_end >=? and tx_start<=? ";
$qh_getCordAllHits= $dbh->prepare($query);
$query="select distinct transcript_id from transcript_exon where organism_version_id=? and gene_prediction_id in (?)";
$qh_getTxByAnnot=$dbh->prepare($query);
$query="select distinct transcript_id from transcript_exon where organism_version_id=? and gene_prediction_id!=3";
$qh_getTxByAllAnnot=$dbh->prepare($query);


$qh_getChr=$dbh->prepare("select chromosome_id from chromosome where chromosome_name=?");
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$query="select distinct gene_prediction_name,g.gene_prediction_id from gene_prediction g,transcript_exon t";
$query.=" where  t.organism_version_id=? and t.transcript_id=? and transcript_name=? ";
$query.=" and t.gene_prediction_id=g.gene_prediction_id";
$qh_getTxPredId=$dbh->prepare($query);

my $getTxExons="select distinct e.exon_id,exon_start,exon_end  from transcript_exon t,exon e ";
   $getTxExons.=" where  transcript_id=? and transcript_name=? and t.organism_version_id=? and gene_prediction_id in (?)";
   $getTxExons.=" and t.organism_version_id=e.organism_version_id and t.exon_id=e.exon_id ";
   $getTxExons.=" order by exon_start asc";
my $qh_getTxExons = $dbh->prepare($getTxExons)or die "Couldn't prepare statement: ".$dbh->errstr;
$query="select count(distinct exon_id) from transcript_exon where transcript_id=? 
        and transcript_name=? and organism_version_id=? and gene_prediction_id in (?)"
my $qh_getExCount=$dbh->prepare($query);

$query="select cdsStart,cdsEnd from transcript_translation ";
$query.=" where transcript_id=? and transcript_name =? and organism_version_id=? and gene_prediction_id in (?) order by cdsStart asc";
my $qh_getCDS=$dbh->prepare($query);

my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $qh_getOrgV = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getOrgUcsc="select ucsc_db  from organism_version where organism_version_id=?";
my $qh_getOrgUcsc = $dbh->prepare($getOrgUcsc)or die "Couldn't prepare statement: ".$dbh->errstr;

my $getAllOrgV="select distinct o.organism_name,ucsc_db as version from organism o,organism_version ov
                  where o.organism_id=ov.organism_id order by organism_name";
my $qh_getAllOrgV = $dbh->prepare($getAllOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;
my $getGenePred="select gene_prediction_name from gene_prediction_by_organism op, gene_prediction p
                   where op.organism_version_id =? and op.gene_prediction_id=p.gene_prediction_id 
                   and gene_prediction_name not in ('all_est','all_mrna','estOrientInfo')";
my $qh_getGenePred = $dbh->prepare($getGenePred)or die "Couldn't prepare statement: ".$dbh->errstr;
my $getTxGene="select distinct gene_name from gene_by_annotation 
                where organism_version_id=? and transcript_id=? and gene_prediction_id in (?)";
my $qh_getTxGene = $dbh->prepare($getTxGene)or die "Couldn't prepare statement: ".$dbh->errstr;
$organism_version=$opt_v if($opt_v);$org_vid=0;$organism_version=~s/\s+//g;
$qh_getOrgV->execute("$organism_version");if($qh_getOrgV->rows>0){($org_vid)=$qh_getOrgV->fetchrow_array();}

sub getGenomicRegion($term){
  ($chrom,$posData)=split(":",$term);$chromosome_id=0;$chrom_start=0;$chrom_end=0;
  ($chrom_start,$chrom_end)=split("-",$posData);$chrom_start=~s/\s+//g;$chrom_end=~s/\s+//g;
   $chrom=~s/CHR//i;$qh_getChr->execute($chrom);
   ($chromosome_id)=$qh_getChr->fetchrow_array() if($qh_getChr->rows>0);
   if($chromosome_id>0){ $mb=1000000;
      if(($chrom_start>0)&&($chrom_end>0)){ #convert cordinates into BP if needed
        $chrom_end *=$mb if($chrom_end=~/\d+\.\d+/);$chrom_start*=$mb if($chrom_start=~/\d+\.\d+/);
      }
   }
  return "$chromosome_id:$chrom_start:$chrom_end";
}
if($opt_s){ # distinct o.organism_name,ov.organism_version_id ,ucsc_db as version
  if($opt_l){
    print "$organism_version gene prediction sources:\n";
    $qh_getGenePred->execute($org_vid);
    if($qh_getGenePred->rows>0){
       while(($gene_prediction_name)=$qh_getGenePred->fetchrow_array()){print "$gene_prediction_name\n";}}
  }
  else{ 
    $qh_getAllOrgV->execute();
    if($qh_getAllOrgV->rows>0){ print "Organism\tVersion\n";
      while(($organism,$version)=$qh_getAllOrgV->fetchrow_array()){print "$organism\t$version\n";}}
  }
}
elsif($opt_l){
   print "$organism_version gene prediction sources:\n";
   $qh_getGenePred->execute($org_vid);
   if($qh_getGenePred->rows>0){
      while(($gene_prediction_name)=$qh_getGenePred->fetchrow_array()){print "$gene_prediction_name\n";}}
}
else{
 $gene_prediction="";$pred_id=0;%predictions=();$prediction_ids="";
 if($opt_a){$gene_prediction=$opt_a;
    @predictions=split(",",$gene_prediction);
   foreach my $pred(@predictions){ #generate gene prediction ids for selected gene predisctions
     $qh_getPredId->execute($pred);
     if($qh_getPredId->rows>0){
       ($pred_id)=$qh_getPredId->fetchrow_array();
       if($prediction_ids eq ""){$prediction_ids="$pred_id";}
       else{$prediction_ids.=",$pred_id";}
     }
   }
 } # default all annotations
 $file_name="";$opt_o=~s/\/\s*$//;$opt_t=~s/\/\s*$//g;$file_name="$opt_o/$organism_version-$gene_prediction.gtf";
 print "$file_name will be generated\n"; exit(0);
 open (TS,">$file_name") or die "$!\n"; $count=0;
 my $tx_qh;
   if($opt_t){ ## a hit on gene_by annotation table
      $tx_qh=$qh_getGenTxHits;$tx_qh->execute($org_vid,$opt_t);
      if($tx_qh->rows<=0){ ## a hit on transcript_by annotation table
         $tx_qh=$qh_getTxHits;$tx_qh->execute($org_vid,$opt_t);
         if($tx_qh->rows<=0){ ## a hit on mgi_gene table (mgi geneid)
            $tx_qh=$qh_getMgiHits;$tx_qh->execute($org_vid,$opt_t);
            if($tx_qh->rows<=0){ ## a hit on cc_founders_transcript table (transcript_name)
               $tx_qh=$qh_getCcTxHits;$tx_qh->execute($pred_id,$opt_t);
                if($tx_qh->rows<=0){ ## a hit on cc_founders_genes table (gene_id)
                   $tx_qh=$qh_getCcGenHits;$tx_qh->execute($pred_id,$opt_t);
                   if($tx_qh->rows<=0){ ## a hit on genomic region
                      ($chromosome_id,$chrom_start,$chrom_end)=split(":",getGenomicRegion($opt_t));
                      if(($chromosome_id>0)&&($chrom_start>0)&&($chrom_end>0)){
                          $tx_qh=$qh_getCordAllHits;$tx_qh->execute($chromosome_id,$chrom_start,$chrom_end);
                      }
                   }
                }
            }
         }
      }
   }
   else{
     if($prediction_ids ne ""){$tx_qh=$qh_getTxByAnnot; $tx_qh->execute($org_vid,$prediction_ids);}
     else{$tx_qh=$qh_getTxByAllAnnot;$tx_qh->execute($org_vid);}
   }
   if($tx_qh->rows>0){
       print TS "seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute\n";
       while(($tx_id)=$tx_qh->fetchrow_array()){
            $qh_getTxcord->execute($tx_id,$org_vid); #get the coordinates of this transcript
            ($chromosome_name,$strand,$tx_start,$tx_end)=$qh_getTxcord->fetchrow_array();
            $qh_getTxAccession->execute($org_vid,$tx_id,$prediction_ids);
            $tx="";$gene="";$exonStarts="";$exonEnds="";$cdsStarts="";$cdsEnds="";
            $qh_getTxExons->execute($tx_id,$org_vid,$prediction_ids);
            $qh_getCDS->execute($tx_id,$org_vid,$prediction_ids);
            while(($exid,$exstart,$exend)=$qh_getTxExons->fetchrow_array()){
                 $exstart+=1; #adjust tx_start to 1-base
                 if($exonStarts eq ""){$exonStarts=$exstart;$exonEnds=$exend;
                 }else{$exonStarts.=",$exstart";$exonEnds.=",$exend";}
             }
           while(($cdsStart,$cdsEnd)=$qh_getCDS->fetchrow_array()){
                 $cdsStart+=1; #adjust tx_start to 1-base
                 if($cdsStarts eq ""){$cdsStarts=$cdsStart;$cdsEnds=$cdsEnd;
                 }else{$cdsStarts.=",$cdsStart";$cdsEnds.=",$cdsEnd";}
             }while(($transcript)=$qh_getTxAccession->fetchrow_array()){
                 if($tx eq ""){$tx=$transcript;}else{$tx.=",$transcript";}}
            $qh_getTxGene->execute($org_vid,$tx_id,$prediction_ids);#adjust tx_start to 1-base
            $tx_start+=1;
            while(($gen)=$qh_getTxGene->fetchrow_array()){
                if($gene eq ""){$gene=$gen;}else{$gene.=",$gen";}} #
         print  TS  "$tx\t$chromosome_name\t$tx_start\t$tx_end\t$strand\t$cdsStarts\t$cdsEnds\t$exonStarts\t$exonEnds\ttranscript\t$gene\n";
         ++$count; 
      }
   }
 print "$file_name generated\n";
}#end of query by -t option
print"program completed\n";
exit(0);

