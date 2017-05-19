#!/usr/bin/perl

#############################################################################
# This script generates aggregate transcript units for 
#  a given organism version [and strain] 

#
# Note: 
# Aggregate transcript regions are organism version-strain specific
# Aggregate transcript regions are generated after all annotations data of 
# a given organism version/strain have been loaded into the database 
# 
# The process:
# 1. sort transcript by chr, strand,start,end
# 2. generate agregate transcript with associated transcript ids list
# 3. generate a db record
#

#  
##############################################################################
use DBI;use vars qw ($opt_h $opt_d); use Getopt::Std;
use Time::localtime; 

use Cwd;
$dbname="graber_transcriptdb"; $host="harlequin"; $user="lnh";$pass="lucie98";
$dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
getopts('hv:a:');
if(($opt_h)|| !($opt_v)){
    print <<HELP;
    **************************************

    Script usage: ./generateAggregate.pl -v organismVersion(s) [-a genePrediction]

HELP
exit;
}
my $dir = getcwd; chomp($dir);#Current working directory
#return strain specific annotations 
$query="select distinct transcript_name from transcript_by_annotation where transcript_id=? ";
$query.=" and gene_prediction_id=?";
$qh_getCCTxid=$dbh->prepare($query);
$query="select distinct gene_name from gene_by_annotation where transcript_id=? ";
$query.=" and gene_prediction_id =?";
$qh_getCCGene=$dbh->prepare($query);

#Exclude strain specific annotations, also exclude, all_est,and xenorefGene
$query="select distinct transcript_name from transcript_by_annotation where transcript_id=? ";
$query.=" and (gene_prediction_id <24 or gene_prediction_id >40)";
$query.=" and gene_prediction_id not in (1,5)";
$qh_getDefaultTxid=$dbh->prepare($query);
$query="select distinct gene_name from gene_by_annotation where transcript_id=? ";
$query.=" and (gene_prediction_id <24 or gene_prediction_id >40)";
$query.=" and gene_prediction_id not in (1,5)";
$qh_getDefaultGene=$dbh->prepare($query);

$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
$qh_getOrgV = $dbh->prepare("select organism_version_id from organism_version where ucsc_db=?");

#get chromosome list
$query="select distinct chromosome_id,strand from transcript where organism_version_id=?";
$qh_getChr=$dbh->prepare($query);

#strain specific
$query="select distinct tx_start,tx_end,ga.transcript_id ";
$query.=" from transcript t, transcript_by_annotation ga";
$query.=" where gene_prediction_id=? and ga.organism_version_id=? ";
$query.=" and ga.organism_version_id=t.organism_version_id ";
$query.=" and  t.chromosome_id=? and t.strand=? ";
$query.=" and ga.transcript_id=t.transcript_id  order by tx_start ";
$qh_getCCRegioncoord=$dbh->prepare($query);

#organism version specific, exclude strain specific annotations
$query="select distinct tx_start,tx_end,ga.transcript_id ";
$query.=" from transcript t, transcript_by_annotation ga";
$query.=" where ga.organism_version_id=? ";
$query.=" and ga.organism_version_id=t.organism_version_id ";
$query.=" and t.chromosome_id=? and t.strand=? ";
$query.=" and  ga.transcript_id=t.transcript_id ";
$query.=" and (gene_prediction_id <24 or gene_prediction_id >40) ";
$query.=" and gene_prediction_id not in (1,5)";
$query.=" order by tx_start ";
$qh_getRegioncoord=$dbh->prepare($query) or die "error -";

###########################################
## transcript units load
$query="create table if not exists transcript_units(";
$query.=" organism_version_id smallint default 0,gene_prediction_id tinyint unsigned default 0,";
$query.=" chromosome_id mediumint unsigned default 0, strand char,";
$query.=" region_start int unsigned default 0, region_end int default 0,";
$query.=" transcript_list text,gene_list text,";
$query.=" index(organism_version_id),index(chromosome_id),index(gene_prediction_id),";
$query.=" index(organism_version_id,strand,region_start))";
$qh_createTranscriptUnit=$dbh->prepare($query);

#
# Default 
#
$query="delete from transcript_units where organism_version_id=? and gene_prediction_id=?";
$qh_deleteCurrentTxUnits=$dbh->prepare($query);
#

$query="load data local infile ? into table transcript_units ignore 1 lines";
$qh_loadTranscriptUnits=$dbh->prepare($query);

$query="select count(*) from transcript_units where organism_version_id=? and gene_prediction_id=?";
$qh_getTxUnitsRowCount=$dbh->prepare($query);

##################################################################################################
# generates aggregate transcript regions for the specified organism version [and strain]
# Each region has coordinates and list of transcript ids that overlap the region
# aggregate regions are indexed by $chromosome_id:$strand:$current_start
# coordinates are sorted by chrom-strand-start 
#
# Process regions by chrom-strand
#
sub genAggRegions{
 ($region_map_ref,$org_vid,$prediction_id,$chromosome_id,$strand)=@_;
 $qh_fd=($prediction_id>0)?$qh_getCCRegioncoord:$qh_getRegioncoord;
 if($prediction_id>0){$qh_fd->execute($prediction_id,$org_vid,$chromosome_id,$strand);}
 else{$qh_fd->execute($org_vid,$chromosome_id,$strand);}
 %chromMap=();$count=0;
 while(($tx_start,$tx_end,$transcript_id)=$qh_fd->fetchrow_array()){
       $chromMap{$tx_start}{$tx_end}=$transcript_id;++$count;}#now process unique regions
 $current_start=0;$current_end=0;$tx_list="";$more=0;$count=0;
 for my $tx_start(sort { $a <=> $b} keys(%chromMap)){
        $tx_id="";@ends=sort {$a <=> $b} keys(%{$chromMap{$tx_start}});$end_index=@ends;
        --$end_index;$tx_end=$ends[$end_index];
       ++$count;
        foreach $tend(@ends){
                $txid=$chromMap{"$tx_start"}{"$tend"};
                $tx_id=($tx_id eq "")?$txid:"$tx_id,$txid";
        }
        if($current_start==0){$current_start=$tx_start;$current_end=$tx_end;$tx_list="$tx_id";}
        else{
             if($tx_start>$current_end){#this is a new genomic region, no overlap with previous region
                ${$region_map_ref}{"$current_start-$current_end"}{"cord"}=$current_end;
                ${$region_map_ref}{"$current_start-$current_end"}{"tx"}=$tx_list;
                $current_start=$tx_start;$current_end=$tx_end;$tx_list="$tx_id";$more=1;
             }
             else{ 
                 $current_end=$tx_end if($current_end<$tx_end); #adjust end coordinate when needed
                 $tx_list.=",$tx_id";$more=1;
       }}
 }##store last
 if($more){
    ${$region_map_ref}{"$current_start-$current_end"}{"cord"}=$current_end;
    ${$region_map_ref}{"$current_start-$current_end"}{"tx"}=$tx_list;
  }
}
#
######### Main program starts here #############################
#
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
$gene_prediction=($opt_a)?$opt_a:"";$orgv_id=0;$prediction_id=0;
if($gene_prediction ne ""){ 
   $qh_getPredId->execute($gene_prediction);($prediction_id)=$qh_getPredId->fetchrow_array();
}
@orgvs=split(",",$opt_v);
foreach $orgv(@orgvs){
 $qh_getOrgV->execute($orgv);($orgv_id)=$qh_getOrgV->fetchrow_array(); 
 if($orgv_id<=0){print "Bad organism version $orgv does not exists in the database \n"; exit;}

 print "\n\n*************************************************************\n";
 print "* Program Started $mon $mday, $yday @ $hour:$min:$sec\n";
 $annotation=($prediction_id>0)?" - $gene_prediction":"";
 print "* Generating unique transcript regions for $orgv $annotation\n";
 print "*************************************************************\n\n";
 $temp_file="tempfile-$orgv-$gene_prediction.txt";
 open(MAIN,">$temp_file");$orgversion="";%novel_gene_map=();
 print MAIN "organism_v\tgene_prediction\tchrom\tstrand\tregionStart\tregionEnd\ttranscriptsList\tgenesList\n";

#
#Generate new aggregate transcript regions by chromosome/strand
#
 $qh_getChr->execute($orgv_id);
 while(($chromosome_id,$strand)=$qh_getChr->fetchrow_array()){
       print "Processing $chromosome_id,$strand\n";
      # next if($chromosome_id!=18);
       %uniq_region_map=();&genAggRegions(\%uniq_region_map,$orgv_id,$prediction_id,$chromosome_id,$strand);
       $qh_tx=($prediction_id>0)?$qh_getCCTxid:$qh_getDefaultTxid;
       $qh_gen=($prediction_id>0)?$qh_getCCGene:$qh_getDefaultGene;
       print "A total of ".keys(%uniq_region_map)." transcrips for $org_vid,$prediction_id\n";
       while(($key,$data)=each(%uniq_region_map)){
              ($start,$end)=split("-",$key);
               @tx_list=split(",",$uniq_region_map{"$key"}{"tx"}); %tx_map=();%gene_map=();
              
              while(@tx_list>0){$tx_id=shift(@tx_list);
                    if($prediction_id>0){
                       $qh_tx->execute($tx_id,$prediction_id);$qh_gen->execute($tx_id,$prediction_id);
                    }else{$qh_tx->execute($tx_id); $qh_gen->execute($tx_id);}
                    while(($tx_name)=$qh_tx->fetchrow_array()){$tx_map{"$tx_name"}=1;}
                    while(($gen_name)=$qh_gen->fetchrow_array()){$gene_map{"$gen_name"}=1;}
               }
               $tx_list=join(",",sort keys(%tx_map)); $gene_list=join(",",sort keys(%gene_map));
                print MAIN "$orgv_id\t$prediction_id\t$chromosome_id\t$strand\t$start\t$end\t$tx_list\t$gene_list\n";
      }
 }
 close(MAIN);
#
#now load db generated file into the database
#
$linecont = `wc -l $temp_file`;chomp($linecont); 
if($linecont=~/^(\d+)\s$temp_file/){$linecont=$1;
  if($linecont>1){--$linecont;
     $qh_createTranscriptUnit->execute();
     $qh_deleteCurrentTxUnits->execute($orgv_id,$prediction_id);
     $qh_loadTranscriptUnits->execute($temp_file);
     $qh_getTxUnitsRowCount->execute($orgv_id,$prediction_id);
     ($rowCount)=$qh_getTxUnitsRowCount->fetchrow_array();
     print  " file loaded -$temp_file\t$linecont\t$rowCount\n";   
   }#else{print "There are no new genes to load\n";}
 }
 #system("rm $temp_file");
}
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "Program ended  -- $mon $mday, $yday @ $hour:$min:$sec\n";



