#!/usr/bin/perl
########################################################################
## This script collects all the repeats for a given organism version
## from Ensembl and load them into graber_transcriptdb.repeats
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : January 2013
#
#  Input : gets a list of all the ensembl schema currently in our database
#  Output: a tab-delimited file with all the genes and their descriptions
#
## set default variables
$dbname="graber_transcriptdb";
$host="harlequin"; $user="lnh";$pass="lucie98";
$basedir="/scratch/data/downloads/ensembl";

use DBI;use vars qw ($opt_h $opt_o);
use Getopt::Std;use LWP 5.64;
use Time::localtime;
$inputfile="";$names="";$outputfile="";@geneList=();
getopts('ho:');
if(($opt_h)||(!$opt_o)){
    print <<HELP;
This script collects all the repeats for a given organism version
from Ensembl and load them into graber_transcriptdb.repeats

Usage:
   perl download_ensRepeats.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  Output directory

Examples:
cmd: ./download_ensRepeats.pl -o $basedir

HELP
exit;
}
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

my $ens_schema="select db_name, schema_version from external_db where db_name like '%ensGene'";
  $qh_schema=$dbh->prepare($ens_schema);
## Get data to load into cgd_transcript table *
$get_repeats="select g.repeat_feature_id ,s.name as chrom,
                     g.seq_region_start as txStart,g.seq_region_end as txEnd,
                     (case when g.seq_region_strand=1 then '+' else '-' end) as strand,
                     repeat_start, repeat_end,g.repeat_consensus_id,repeat_name,repeat_type,
                     if((g.seq_region_strand= 1),g.seq_region_start+g.repeat_start-1,
                         g.seq_region_end-g.repeat_end+1)as repStart,
                    if((g.seq_region_strand=1),g.seq_region_start+g.repeat_end-1,g.seq_region_end-g.repeat_start+1)as repEnd
                    from repeat_feature g,repeat_consensus ge,seq_region s ,coord_system c
                    where ge.repeat_consensus_id=g.repeat_consensus_id and  g.seq_region_id=s.seq_region_id
                    and s.coord_system_id=c.coord_system_id and c.rank=1 
                    and  s.seq_region_id=g.seq_region_id";
$get_repeat_name="select distinct repeat_name from repeat_consensus";
$get_repeat_type="select distinct repeat_type from repeat_consensus";

my $create_rep_temp = "create temporary table ensembl_repeats_temp(
          organism_version_id smallint default 0,repeat_feature_id int default 0,
          chromosome_id int default 0,seq_region_start int unsigned default 0,seq_region_end int unsigned default 0,
          strand char,repeat_start smallint unsigned default 0, repeat_end smallint unsigned default 0,
          repeat_consensus_id int unsigned default 0,repeat_name_id smallint unsigned default 0,repeat_type_id tinyint default 0,
          repChrStart  int unsigned default 0,repChrEnd  int unsigned default 0, found tinyint default 0,
          index(organism_version_id),index(repeat_feature_id),index(chromosome_id,repChrStart))";
my $create_repex_temp="create temporary  table ensembl_repeats_exon_temp(
 organism_version_id smallint default 0,
 repeat_feature_id int unsigned default 0,
 exon_id int unsigned default 0,found tinyint default 0,
 foreign key(repeat_feature_id) references repeats(repeat_feature_id),
 foreign key(exon_id) references exon(exon_id),index(organism_version_id),
 index(repeat_feature_id),index(exon_id))";

my $getChr="select distinct t.chromosome_id, c.chromosome_name from transcript t, chromosome c 
             where organism_version_id=? and t.chromosome_id=c.chromosome_id";
my $qh_getChr = $dbh->prepare($getChr);

my $getOverlapE="select exon_id from exon where organism_version_id=? ";
   $getOverlapE.=" and chromosome_id=? and strand=? and exon_end >= ? and exon_start<=?  order by exon_start";
my $qh_getOverlapE = $dbh->prepare($getOverlapE);

my $qh_drop_rep_temp =$dbh->prepare("drop temporary table if exists ensembl_repeats_temp");
my $qh_create_rep_temp = $dbh->prepare($create_rep_temp);
my $qh_load_rep =$dbh->prepare("load data local infile ? into table ensembl_repeats_temp ignore 1 lines");
my $qh_get_rep_rowcount=$dbh->prepare("select count(*) as rowcount from ensembl_repeats_temp ");

my $qh_drop_repex_temp =$dbh->prepare("drop temporary table if exists ensembl_repeats_exon_temp");
my $qh_create_repex_temp = $dbh->prepare($create_repex_temp);
my $qh_load_repex =$dbh->prepare("load data local infile ? into table ensembl_repeats_exon_temp ignore 1 lines");
my $qh_get_repex_rowcount=$dbh->prepare("select count(*) as rowcount from ensembl_repeats_exon_temp");

my $update_rep= "update ensembl_repeats_temp t, repeats m set found=1 
                 where t.organism_version_id=m.organism_version_id
                 and t.repeat_feature_id=m.repeat_feature_id";
my $insert_rep = "insert into repeats(organism_version_id ,repeat_feature_id ,
                  chromosome_id,seq_region_start ,seq_region_end,strand,repeat_start,repeat_end,
                  repeat_consensus_id,repeat_name_id,repeat_type_id,repeatChrStart,repeatChrEnd) ";
 $insert_rep.="select distinct organism_version_id ,repeat_feature_id ,
             chromosome_id,seq_region_start ,seq_region_end,strand,repeat_start,repeat_end,
             repeat_consensus_id,repeat_name_id,repeat_type_id,repChrStart,repChrEnd from ensembl_repeats_temp where found=0 ";

my $update_repex= "update ensembl_repeats_exon_temp t, repeats_exon m set found=1 
                 where t.organism_version_id=m.organism_version_id
                 and t.repeat_feature_id=m.repeat_feature_id and t.exon_id=m.exon_id";

my $insert_repex = "insert into repeats_exon(organism_version_id ,repeat_feature_id,exon_id)";
   $insert_repex .= "select distinct organism_version_id ,repeat_feature_id ,exon_id 
                     from ensembl_repeats_exon_temp where found=0";
 
my $qh_update_rep = $dbh->prepare($update_rep)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insert_rep = $dbh->prepare($insert_rep)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_update_repex = $dbh->prepare($update_repex)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insert_repex = $dbh->prepare($insert_repex)or die "Couldn't prepare statement: " . $dbh->errstr;

my $getOrg="select organism_version_id from organism_version where ucsc_db=?";   #get current organisms list 
my $qh_orgvlist= $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;

############################## load repeats name and type ##################################
$qh_drop_repname=$dbh->prepare("drop temporary table if exists repeat_name_temp");
$qh_create_repname=$dbh->prepare("create temporary table repeat_name_temp(repeat_name varchar(255),found tinyint default 0)");
$qh_load_repname=$dbh->prepare("load data local infile ? into table repeat_name_temp ignore 1 lines");
$qh_get_repname_rowcount=$dbh->prepare("select count(*) as rowcount from repeat_name_temp");
$qh_update_repname=$dbh->prepare("update repeat_name_temp t, repeat_name m set found=1 where t.repeat_name=m.repeat_name");
$qh_insert_repname = $dbh->prepare("insert into repeat_name(repeat_name) select distinct repeat_name from repeat_name_temp where found=0");
$qh_get_repname=$dbh->prepare("select repeat_name_id,repeat_name from repeat_name");

$qh_drop_reptype=$dbh->prepare("drop temporary table if exists repeat_type_temp");
$qh_create_reptype=$dbh->prepare("create temporary table repeat_type_temp(repeat_type varchar(255),found tinyint default 0)");
$qh_load_reptype=$dbh->prepare("load data local infile ? into table repeat_type_temp ignore 1 lines");
$qh_get_reptype_rowcount=$dbh->prepare("select count(*) as rowcount from repeat_type_temp");
$qh_update_reptype=$dbh->prepare("update repeat_type_temp t, repeat_type m set found=1 where t.repeat_type=m.repeat_type");
$qh_insert_reptype = $dbh->prepare("insert into repeat_type(repeat_type) select distinct repeat_type from repeat_type_temp where found=0");
$qh_get_reptype=$dbh->prepare("select repeat_type_id,repeat_type from repeat_type");
####################################### Queries ##########################################

##***  *** e2.seq_region_end-tl.seq_start+1
$host="ensembldb.ensembl.org";$port="5306";$user="anonymous";
$basedir="/scratch/data/downloads/ensembl";
%versionmap=("mm9"=>"mm9","mm10"=>"GRCm38","GRCm38"=>"mm10","hg18"=>"hg18","hg19"=>"hg19");
$qh_schema->execute();if($opt_o){$basedir=$opt_o;}
while(($gene_prediction, $eschema)=$qh_schema->fetchrow_array()){
   $dbh2 = DBI->connect("DBI:mysql:$eschema:$host:$port",$user);
   ($oversion,$token)=split("-",$gene_prediction);
   $oversion=$versionmap{"$oversion"};$qh_orgvlist->execute($oversion);
   ($orgvid)=$qh_orgvlist->fetchrow_array();%chrom_map=();$qh_getChr->execute($orgvid);
   if($qh_getChr->rows>0){
       while(($chr_id,$chr_name)=$qh_getChr->fetchrow_array()){$chrom_map{$chr_name}=$chr_id;}
    }
  
    
   print "Processing $orgvid,$oversion,$gene_prediction, $eschema\n";
   if(!(-d $basedir)){mkdir $basedir, 0777;}if(!(-d "$basedir/$oversion")){mkdir "$basedir/$oversion", 0777;}
   if(!(-d "$basedir/$oversion/$eschema")){mkdir "$basedir/$oversion/$eschema", 0777;}
   open(LOG,">$basedir/$oversion/$eschema/repeats-log.txt");$tm = localtime;
   $filename="$basedir/$oversion/$eschema/$gene_prediction-repeats.txt";
   $rep_exon_file="$basedir/$oversion/$eschema/$gene_prediction-repeats_exon.txt";
    ### get repeat name and type 
   $qh_get_repeat_name= $dbh2->prepare($get_repeat_name);$qh_get_repeat_name->execute();
   $qh_get_repeat_type= $dbh2->prepare($get_repeat_type);$qh_get_repeat_type->execute();
   $repeat_name_file="$basedir/$oversion/$eschema/repeat_name.txt";
   $repeat_type_file="$basedir/$oversion/$eschema/repeat_type.txt";
   %repeatname=();%repeattype=();
   if($qh_get_repeat_name->rows>0){
      open(NAME,">$repeat_name_file");print NAME "repeat_name\texists\n";
      while(($repeat_name)= $qh_get_repeat_name->fetchrow_array()){print NAME "$repeat_name\t0\n"; }
      close(NAME);
      $qh_drop_repname->execute();$qh_create_repname->execute();
      $qh_load_repname->execute($repeat_name_file);
      $qh_get_repname_rowcount->execute(); $qh_update_repname->execute();$qh_insert_repname ->execute(); $qh_get_repname->execute();
      if($qh_get_repname->rows>0){
       while(($type_id,$type_name)=$qh_get_repname->fetchrow_array()){$repeatname{"$type_name"}=$type_id;}
     }
   }
   if($qh_get_repeat_type->rows>0){
     open(NAME,">$repeat_type_file");print NAME "repeat_type\texists\n";
     while(($repeat_type)= $qh_get_repeat_type->fetchrow_array()){print NAME "$repeat_type\t0\n"; 
     }close(NAME);
     $qh_drop_reptype->execute();$qh_create_reptype->execute();$qh_load_reptype->execute($repeat_type_file);
     $qh_get_reptype_rowcount->execute();$qh_update_reptype->execute();$qh_insert_reptype->execute();$qh_get_reptype->execute();
     if($qh_get_reptype->rows>0){
       while(($type_id,$type_name)=$qh_get_reptype->fetchrow_array()){$repeattype{"$type_name"}=$type_id;}
     }
   }
   $qh_get_repeats = $dbh2->prepare($get_repeats);
   #open(REP,">$rep_exon_file");
   #if(!(-f $filename)){
   open(OUT,">$filename"); #data will be loaded into ensembl_genes table
   if(OUT){
      $header="organism_version_id\trepeat_feature_id\tchromosome_id\t";
      $header.="seq_region_start\tseq_region_end\tstrand\trepeat_start\trepeat_end\t";
      $header.="repeat_consensus_id\trepeat_name\trepeat_type\trepStart\trepEnd\tfound\n";
      print OUT $header;
      $qh_get_repeats->execute();
      my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
      print LOG "\n*************************************************************\n";
      print LOG "Starting load process :  $mday/$mon/$yday @ $hour:$min:$sec \n";
      print LOG "\n*************************************************************\n";
      if($qh_get_repeats->rows>0){
       while(($repeat_id,$chr,$txStart,$txEnd,$strand,$rstart,$rend,$c_id,$rname,$rtype,$repStart,$repEnd)= $qh_get_repeats->fetchrow_array()){
            $chr_id=$chrom_map{$chr};$rnameid=$repeatname{"$rname"};$rtypeid=$repeattype{"$rtype"};
            print OUT "$orgvid\t$repeat_id\t$chr_id\t$txStart\t$txEnd\t$strand\t$rstart\t$rend\t$c_id\t$rnameid\t$rtypeid\t$repStart\t$repEnd\t0\n";
         }
      }
   }close(OUT);
 # }
   #now get overlap exons
   #a repeat row looks like: 
   #organism_version_id	repeat_feature_id	chromosome_id	seq_region_start	seq_region_end	strand	   repeat_start	repeat_end	  repeat_consensus_id	repeat_class	repeat_type	found
  #61	5833	20	85819919	85820358	-	5871	6279	31	LINE/L1	Type I Transposons/LINE	0
  #repeat start and end are relative to the consensus but not the genome. So I need to compute 
 q{ print "getting overlap exons\n";
   if(-f $filename){
     open(IN,"$filename");$header=<IN>;
     if(REP){
      $header="organism_version_id\trepeat_feature_id\texon_id\tfound\n";
      print REP $header;
      while(<IN>){ chomp($_);
       ($orgvid,$repeat_id,$chr_id,$txStart,$txEnd,$strand,$rstart,$rend,$c_id,$rname,$rtype,$repStart,$repEnd,$found)=split("\t",$_);
        #org_vid,$chrom_id,$strand,$chrom_start,$chrom_end
        $qh_getOverlapE->execute($orgvid,$chr_id,$strand,$repStart,$repEnd);
        if($qh_getOverlapE->rows>0){
           while(($exon_id)=$qh_getOverlapE->fetchrow_array()){print REP "$orgvid\t$repeat_id\t$exon_id\t0\n";}
        }
      }
     } close(REP);
   }
   };
   ###################################################################
   #LOAD repeats
   print "Loading $filename -- check data first\n";
   $linecont = `wc -l $filename`;chomp($linecont);
   if($linecont=~/^(\d+)\s$filename/){$linecont=$1;
      if($linecont>0){
         $qh_drop_rep_temp->execute()or die "bad query ".mysql_error()."\n";
         $qh_create_rep_temp->execute()or die "bad query ".mysql_error()."\n";
         $qh_load_rep->execute("$filename") or die "bad query ".mysql_error()."\n";
         $qh_get_rep_rowcount->execute();($rowcount)=$qh_get_rep_rowcount->fetchrow_array();
         ####################3
         --$linecont;
         if($linecont==$rowcount){$qh_update_rep->execute();$qh_insert_rep->execute();
           print "Processed $linecont lines  and loaded $rowcount lines\n";}
         else{ 
             print LOG "$filename Load failed:only $rowcount of $linecont lines were loaded\n";
               #the bulk load did not go well
         }
     }
   } 
   #LOAD repeats-exon
  q{ print "Loading $rep_exon_file\n";
   $linecont = `wc -l $rep_exon_file`;chomp($linecont);
   if($linecont=~/^(\d+)\s$rep_exon_file/){$linecont=$1;
      if($linecont>0){
         $qh_drop_repex_temp->execute()or die "bad query ".mysql_error()."\n";
         $qh_create_repex_temp->execute()or die "bad query ".mysql_error()."\n";
         $qh_load_repex->execute("$rep_exon_file") or die "bad query ".mysql_error()."\n";
         $qh_get_repex_rowcount->execute();($rowcount)=$qh_get_repex_rowcount->fetchrow_array();
         ####################3
         --$linecont;
         if($linecont==$rowcount){$qh_update_repex->execute();$qh_insert_repex->execute();
           print "Processed $linecont lines  and loaded $rowcount lines\n";}
         else{ 
             print LOG "$filename Load failed:only $rowcount of $linecont lines were loaded\n";
               #the bulk load did not go well
         }
     }
   } 
 };
   ###################################################################
   $tm = localtime;
   my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
  ($tm->mon)+1, ($tm->year)+1900);
   print LOG "\n*************************************************************\n";
   print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
   print LOG "\n*************************************************************\n";
   close(LOG);
 
  $dbh2->disconnect;
 ###########################################################
}
 
print "Porgram complete\n";
exit(0);

