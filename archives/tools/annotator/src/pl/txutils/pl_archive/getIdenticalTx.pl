#!/usr/bin/perl

###################################################################################################
## getIdenticalTx.pl
#  For a given organism version, this scripts returns a tally of
#  identical transcripts that differ by coding region.
# Note:
#    We say that two or more transcripts,from the same source or different sources,
#        with different accession ids are identicals iff
#    1. they share the same genomic location (chr,strand, txStart,txEnd)
#    2.  and they have the same number of exons 
#    3.  and all their exons are the same (same exon starts , and same exon ends
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : April 2013

#   Usage: perl getCordinates.pl -v <organism version> [-a <gene prediction list>] [-o <outputDir>]
#
#   Where: 
#          -o  is the output file name (default standard out)
#          -v  is the ucsc organism version (default mm9)
#          -a  a commas separated list of gene prediction names (default knowGene,ensGene)
#          -l  display the list of current gene predictions for the specified organism version
#          -s  display the list of current organism versions
# Note: 
###################################################################################################
use DBI;
use POSIX;
use XML::Simple;
use Time::localtime;
# create object

use vars qw ($opt_h $opt_l $opt_s $opt_o $opt_v $opt_a);
use Getopt::Std;

getopts('hlso:v:a:');
if($opt_h||(!$opt_v&&!$opt_s&&!$opt_l)) {
    print <<HELP;
######################################################################
#  For a given organism version, this scripts returns the tally of
#  identical transcripts that differ by coding region.
# Note:
#    We say that two or more transcripts with different accession ids are identicals iff
#    1. they share the same genomic location (chr,strand, txStart,txEnd)
#    2. and they have the same number of exons 
#    3. and all their exons are the same (same exon starts , and same exon ends  
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : April 2013

#   Usage: perl getCordinates.pl -v <organism version> [-a <gene prediction list>] [-o <outputDir>]
#
#   Where: 
#          -o  is the output file name (default standard out)
#          -v  is the ucsc organism version (default mm9)
#          -a  a commas separated list of gene prediction names (default knowGene,ensGene)
#          -l  display the list of current gene predictions for the specified organism version
           -h  displays this help message
           -s  displays organism versions
 
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
$query="select distinct transcript_name,exon_id from transcript_exon where transcript_id=? 
        and organism_version_id=? order by transcript_name,exon_id ";
$qh_getTxAll=$dbh->prepare($query);
$query="select distinct transcript_name,exon_id from transcript_exon where transcript_id=? 
        and organism_version_id=? and gene_prediction_id=? order by transcript_name,exon_id ";
$qh_getTxByAnnot=$dbh->prepare($query);

$query="Select distinct gene_prediction_id from transcript_by_annotation where transcript_id=?";
$qh_txPred=$dbh->prepare($query);

$query="select transcript_id ,tx_start,tx_end,chromosome_name,strand 
        from transcript t, chromosome c where organism_version_id=? and t.chromosome_id=c.chromosome_id";
$qh_getAllCordHits= $dbh->prepare($query);
$query="select t.transcript_id ,tx_start,tx_end,chromosome_name,strand 
        from transcript t, transcript_by_annotation ta,chromosome c where ta.organism_version_id=? 
        and ta.gene_prediction_id=? and ta.transcript_id=t.transcript_id 
        and ta.organism_version_id=t.organism_version_id and t.chromosome_id=c.chromosome_id";
$qh_getPredidCordHits= $dbh->prepare($query);

$qh_getChr=$dbh->prepare("select chromosome_name from chromosome where chromosome_id=?");
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$qh_getPredName=$dbh->prepare("select gene_prediction_name from gene_prediction where gene_prediction_id=?");

$query="select distinct gene_prediction_name from gene_prediction g, transcript_by_annotation ta
        where ta.transcript_id=? and ta.transcript_name in (?) and ta.gene_prediction_id=g.gene_prediction_id";
$qh_getPredList=$dbh->prepare($query);

$query="Select distinct cdsStart,cdsEnd from transcript_translation ";
$query.=" where transcript_id=? and transcript_name=? and organism_version_id=? and gene_prediction_id=? order by cdsStart asc";
my $qh_getPredCDS=$dbh->prepare($query);

my $qh_getOrgV = $dbh->prepare("select organism_version_id from organism_version where ucsc_db=?");
my $qh_getOrgUcsc = $dbh->prepare("select ucsc_db  from organism_version where organism_version_id=?");
my $getTxGene="select distinct gene_name from gene_by_annotation 
                where organism_version_id=? and transcript_id=? and gene_prediction_id in (?)";
my $qh_getTxGene = $dbh->prepare($getTxGene)or die "Couldn't prepare statement: ".$dbh->errstr;
my $qh_getTxGeneAll = $dbh->prepare("select distinct gene_name from gene_by_annotation where organism_version_id=? and transcript_id=?");


sub displayAnnotSources{
 ($org_version)=@_;$url="$webservice?v=$org_version&l=1";
  $xml_file="annot_list.xml";if(-e "annot_list.xml"){system("rm $xml_file");}
  system("wget -q -O $xml_file '$url'");
  $xml = XML::Simple->new (ForceArray => ['prediction','name']);
  if(-e $xml_file){$data =$xml->XMLin("$xml_file");
    foreach my $source(@{$data->{prediction}}){print $source->{name}->[0]. "\n";}
  }if(-e "annot_list.xml"){system("rm $xml_file");}
}
sub displayOrganismVersions{
  ($type)=@_;
  $url="$webservice?$type";$xml_file="annot_list.xml";if(-e "annot_list.xml"){system("rm $xml_file");}
  system("wget -q -O $xml_file '$url'");
  $xml = XML::Simple->new (ForceArray =>1);
  if(-e $xml_file){ $data =$xml->XMLin("$xml_file");
    foreach my $organism(@{$data->{organism}}){print $organism->{name}[0] ."\t".$organism->{version}[0]. "\n";}
  }if(-e "annot_list.xml"){system("rm $xml_file");}
}
sub getTranscript{
 ($tx_id,$org_vid,$tx_start,$tx_end,$chrom_name,$txstrand)=@_;
 %tx_map=();$prev_tx="";$exids="";
 #I could change this to first get all the accessions associated with txid
 #then for each accession get exon
 $qh_txPred->execute($tx_id);
 while(($pred_id)=$qh_txPred->fetchrow_array()){
   $qh_getTxByAnnot->execute($tx_id,$org_vid,$pred_id);
   while(($tx,$exid)=$qh_getTxByAnnot->fetchrow_array()){#index transcripts that share same set of exons
     if($prev_tx eq ""){$exids="$exid";$prev_tx=$tx;}elsif($prev_tx ne "$tx"){ #this is a new accession id
        $tx_map{"$exids"}{"$prev_tx"}=1;
        $exids="$exid";$prev_tx=$tx;
      }else{$exids.=",$exid";}
   }
 }$row=""; $tx_map{"$exids"}{"$prev_tx"}=1;
 $qh_getTxGeneAll->execute($org_vid,$tx_id);$gene="";
 while(($gen)=$qh_getTxGeneAll->fetchrow_array()){ #get gene list of this transcript region
       if($prediction_id>24&&$prediction_id<=40){
          if($ccgenemap{"$prediction_id"}{"$gen"}){$gen=$ccgenemap{"$prediction_id"}{"$gen"}.",$gen";}
       }if($gene eq ""){$gene=$gen;}else{$gene.=",$gen";}
  }$row="";
  foreach $tx_exons(keys(%tx_map)){ next if(!($tx_exons=~/\d+/));
    @tx_count=sort(keys(%{$tx_map{"$tx_exons"}})); $cdsStarts="";$cdsEnds=""; $exonStarts="";$exonEnds="";
    $tx_string="";foreach $tx(@tx_count){if($tx_string eq ""){$tx_string="'$tx'";}else{$tx_string.=",'$tx'";}}
    $query="Select distinct cdsStart,cdsEnd from transcript_translation ";
    $query.=" where transcript_id=$tx_id and organism_version_id=$org_vid and transcript_name in ($tx_string) order by cdsStart asc";
    $cds_ary_ref=$dbh->selectall_arrayref($query);$cds_count=0;
    foreach $cds(@{$cds_ary_ref}){($cdsStart,$cdsEnd)=@{$cds}; $cdsStart+1; next if(($cdsEnd-$cdsStart)<=1);
         $cdsStarts=($cdsStarts eq "")? "$cdsStart":"$cdsStarts,$cdsStart";$cdsEnds=($cdsEnds eq "")?"$cdsEnd":"$cdsEnds,$cdsEnd";
         ++$cds_count;
    }@exon_count=split(",",$tx_exons);next if((@tx_count<=1)||((@tx_count>1)&&($cds_count<=1)));
    $query="select exon_start,exon_end from exon where exon_id in \($tx_exons\) order by exon_start";
    $exon_ary_ref=$dbh->selectall_arrayref($query);
    foreach $ex(@{$exon_ary_ref}){
      ($exon_start,$exon_end)=@{$ex};$exon_start+=1;
      if($exonStarts eq ""){ $exonStarts="$exon_start";$exonEnds=$exon_end;}
      else{$exonStarts.=",$exon_start";$exonEnds.=",$exon_end";}
    }$annot="";
   $query="select distinct gene_prediction_name from gene_prediction g, transcript_by_annotation ta ";
   $query.="where ta.transcript_id=$tx_id and ta.transcript_name in ($tx_string) and ta.gene_prediction_id=g.gene_prediction_id";
   $pred_arr_ref=$dbh->selectall_arrayref($query);
   foreach $pred(@{$pred_arr_ref}){($pred_name)=@{$pred};
      if($annot eq ""){$annot="$pred_name";}else{$annot.=",$pred_name";}}
   $tx_string=~s/'//g;
    $row.="$tx_string\t$chrom_name\t$txstrand\t$tx_start\t$tx_end\t$cdsStarts\t$cdsEnds\t".@exon_count."\t$exonStarts\t$exonEnds";
    $row.="\t$gene\t$annot\n";
 }
  return $row;
 }
############################
$organism_version=($opt_v)?$opt_v:"";$org_vid=0;
if($organism_version eq ""&& !($opt_s)&&!($opt_g)){
  print "Must specify the organism version. For example -v mm9 .Please run the script with the -h option for usage\n";
  exit(0);
}$organism_version=~s/\s+//g;
if($opt_s){# distinct o.organism_name,ov.organism_version_id ,ucsc_db as version
    $type="s=1";
    print "Organism\tVersion\n";displayOrganismVersions($type);
}elsif($opt_l){ #display list of annotations for this organism version
   print "Annotation sources for $organism_version\n"; displayAnnotSources($organism_version);
}
elsif($opt_g){ #display list of organism versions with genome data
    $type="g=1";
    print "Organism\tVersion\n";displayOrganismVersions($type);
}else{
 $qh_getOrgV->execute("$organism_version");
 if($qh_getOrgV->rows>0){($org_vid)=$qh_getOrgV->fetchrow_array();}
 $gene_prediction="";$pred_id=0;%predictions=();$prediction_ids="";
 if($opt_a){$gene_prediction=$opt_a;@predictions=split(",",$gene_prediction);
    foreach my $pred(@predictions){ #generate gene prediction ids for selected gene predisctions
       $qh_getPredId->execute($pred);
       if($qh_getPredId->rows>0){($pred_id)=$qh_getPredId->fetchrow_array();
          if($prediction_ids eq ""){$prediction_ids="$pred_id";}
          else{$prediction_ids.=",$pred_id";}
       }
    }
  } # default to all annotations
 #set the result path
 $file_name="$organism_version-$gene_prediction-IdenticalTxWithDiffCDS.txt";$type=~s/\/\s*$//g;$feature_type="Transcript";
 if($opt_f){$file_name=$opt_f;}
 if($opt_o){$opt_o=~s/\/\s*$//;$file_name="$opt_o/$file_name";}
 open (TS,">$file_name") or die "$!\n"; $count=0;my $tx_qh;
 print TS "Name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStarts\tcdsEnds\texonCount\texonStarts\texonEnds\tName2\tsources\n";
 if(($prediction_ids ne "")||($prediction_ids>0)){
   $qh_getPredidCordHits->execute($org_vid,$prediction_ids);$tx_qh=$qh_getPredidCordHits;}
 else{$qh_getAllCordHits->execute($org_vid);$tx_qh=$qh_getAllCordHits;}
 if($tx_qh->rows>0){ $rows="";
  $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
       ($tm->mon)+1, ($tm->year)+1900);
  print "\n*************************************************************\n";
  print "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print "\n*************************************************************\n";
  $count=0;
   while(($tx_id ,$tx_start,$tx_end,$chromosome_name,$strand)=$tx_qh->fetchrow_array()){
       $rows=getTranscript($tx_id,$org_vid,$tx_start,$tx_end,$chromosome_name,$strand);
      print TS $rows;++$count; print "$count processed \n" if($count%100000==0);
   }
  $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
       ($tm->mon)+1, ($tm->year)+1900);
  print "\n*************************************************************\n";
  print "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print "\n*************************************************************\n";
 }
 print "$file_name generated\n";
}#end of query by -t option
print"program completed\n";
exit(0);

