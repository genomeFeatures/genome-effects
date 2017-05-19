#!/usr/bin/perl
#########################################################################################
## This script collects all the gene annotations for a given organism version
## from Ensembl and load them into graber_transcriptdb
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : November 2012
#
#  Input : a directory path where to process files
#  Output: a tab-delimited file with all the annotations(gene, transcripts,exon,cds)
#  and their cordinates
#
## set default variables
## Ensembl release 71 is scheduled for release in April 2013.
## We expect this release to include, among other things:
## check esnsembl blog http://www.ensembl.info/
##############################################################################################
$dbname="graber_transcriptdb";$host="harlequin"; $user="lnh";$pass="lucie98";
$basedir="/scratch/data/downloads/ensembl";
#$oversion="mm9"; $gene_prediction="mm9-ensGene";

use DBI;use vars qw ($opt_h $opt_a $opt_o $opt_v $opt_e);
use Getopt::Std;use LWP 5.64; use Time::localtime;
use POSIX;

$inputfile="";$names="";$outputfile="";@geneList=();
getopts('ho:a:v:e:');
if(($opt_h)||!$opt_a||!$opt_v){
    print <<HELP;
 This script collects all the gene annotations for a given organism version
 from Ensembl and load them into graber_transcriptdb

Usage:
   perl download_ensAnnotations.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  required,Output directory
   -v  required,organism version
   -a  required,gene prediction name (example mm9-ensGene)
   -e  ensembl schema version

Examples:
cmd: ./download_ensAnnotations.pl -o $basedir -a $gene_prediction -v $oversion [-e 68 ]

HELP
exit;
}
sub load_scriptConfig{
  ($config_file,$configMapRef)=@_;@scripts=`cat $config_file`;
   while(@scripts){
         $line=shift @scripts; chomp($line);($name,$script)=split("=",$line);${$configMapRef}{"$name"}=$script;
   }
}
my $dir = getcwd; chomp($dir);
#
# Load scripts path into data structure, assumption is the config file with 
# perl scripts is under the working directory
#
%configMap=();load_scriptConfig("$dir/script_config.txt",\%configMap);
$script_base=$configMap{"sript_base"};$pl_base=$configMap{"gtx_pl"};
$parse_ensAnnot_pl=$configMap{"parseAndLoad_ensAnnotations"};
$load_annotations="$script_base/$pl_base/$parse_ensAnnot_pl";
#**** ensGenes ***************************
my $drop_ens_temp = "drop temporary table if exists ensembl_genes_temp";
my $create_ens_temp = "create temporary table ensembl_genes_temp(gene_name varchar(255),
          biotype varchar(100),organism_version_id smallint default 0,description text,index(organism_version_id))";
my $load_this = "load data local infile ? into table ensembl_genes_temp ignore 1 lines";
my $get_ens_genes_rowcount="select count(*) as rowcount from ensembl_genes_temp ";
my $insert_ensgene = "insert ignore into ensembl_genes(gene_name,biotype,organism_version_id,description) ";
   $insert_ensgene.="select distinct gene_name,biotype,organism_version_id,description from ensembl_genes_temp";

##***  *** e2.seq_region_end-tl.seq_start+1

$get_cds="select distinct t.stable_id as name,
       if((t.seq_region_strand= 1),e.seq_region_start+tl.seq_start-1,
          e2.seq_region_end-tl.seq_end+1)as cdsStart,
      if((t.seq_region_strand=1),e2.seq_region_start+tl.seq_end-1,e.seq_region_end-tl.seq_start+1)as cdsEnd
    from seq_region s,transcript t,translation tl,exon e,exon e2
    where s.coord_system_id=1 and   s.seq_region_id = t.seq_region_id
    and   t.transcript_id=tl.transcript_id 
    and   tl.start_exon_id = e.exon_id and tl.end_exon_id = e2.exon_id";

## Get data to load into cgd_transcript table *
$get_transcript="select distinct g.stable_id as name,s.name as chrom,
                 g.seq_region_start as txStart,g.seq_region_end as txEnd,
                 (case when g.seq_region_strand=1 then '+' else '-' end) as strand,
                  ge.stable_id as name2,ge.biotype,ge.description 
                 from transcript g, gene ge,seq_region s ,coord_system c
                 where g.gene_id=ge.gene_id and s.coord_system_id=c.coord_system_id and c.rank=1 
                 and  g.seq_region_id=s.seq_region_id";


## get data to load into exon_desc *
#### load ensembl transcript exons table ****
$get_exons="select distinct t.stable_id as name,e.seq_region_start as exonStart,";
$get_exons.="e.seq_region_end as exonEnd,phase, e.stable_id as exonid";
$get_exons.=" from exon_transcript et,transcript t,exon e,seq_region s,coord_system c ";
$get_exons.=" where et.transcript_id=t.transcript_id and et.exon_id= e.exon_id ";
$get_exons.=" and e.seq_region_id= s.seq_region_id";
$get_exons.=" and s.coord_system_id=c.coord_system_id and c.rank=1 ";
$get_exons.=" order by et.transcript_id, e.seq_region_start";


$mysqlcmd="mysql -hensembldb.ensembl.org -P 5306 -uanonymous -e 'show databases'";
@databases=`$mysqlcmd`; 
############ TO DO ###############################################################
# I need to create a map between UCSC organism versions                          #
# and ensembl organism schema versions                                           #
# That way, given UCSC version, I can easily get the corresponding               #
# latest ensembl schema version for that organism.                               #
# for eaxmple : mm9 -> *_core_xx_37.. mm10 -> *_core_xx_38**                     #
##################################################################################
%versionmap=("mm9"=>"37","mm10"=>"38","GRCm38"=>"38","hg18"=>"36","hg19"=>"37");
%schema_core=("mm9"=>"mus_musculus_core_","mm10"=>"mus_musculus_core_","GRCm38"=>"mus_musculus_core_",
                "hg18"=>"homo_sapiens_core_","hg19"=>"homo_sapiens_core_");

$oversion=$opt_v;$gene_prediction=$opt_a;@mus=();
$versionid=$versionmap{"$oversion"};$core=$schema_core{"$oversion"};
#set the right ensembl schema to connect to
if($opt_e){$token="$core$opt_e"."_$versionid";@mus=grep(/$token/,@databases);}
else{ @mus=grep(/$core\d+_$versionid/,@databases);}$dbid=0;$eschema="";
foreach $db(@mus){chomp($db);$db=~/$core(\d+)_$versionid/;
  if($opt_e){if($opt_e==$1){$eschema=$db; last;}}
  if($dbid<$1){$eschema=$db;}
}
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);
$host="ensembldb.ensembl.org";$port="5306";$user="anonymous";
my $dbh2 = DBI->connect("DBI:mysql:$eschema:$host:$port",$user);

my $qh_drop_ens_temp = $dbh->prepare($drop_ens_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_create_ens_temp = $dbh->prepare($create_ens_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_load_this = $dbh->prepare($load_this)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_get_ens_genes_rowcount = $dbh->prepare($get_ens_genes_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insert_ensgene = $dbh->prepare($insert_ensgene)or die "Couldn't prepare statement: " . $dbh->errstr;

## get the date of the last update
my $get_lastUpdate="select schema_version from external_db where db_name=?";
my $qh_lastUpdate= $dbh->prepare($get_lastUpdate);
## update the external_db
my $Update_schema="update external_db set schema_version= ? where db_name=?";
my $qh_Update_schema= $dbh->prepare($Update_schema);

my $qh_insertNewdb=$dbh->prepare("Insert into gene_prediction(gene_prediction_name) values(?)");
my $qh_getNewdbId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name=?");
my $qh_insertExternaldb=$dbh->prepare("insert into external_db values(?,?,?)");
my $qh_deleteSchema=$dbh->prepare("delete from external_db where db_name=?");
my $qh_get_exons= $dbh2->prepare($get_exons);
my $qh_get_cds= $dbh2->prepare($get_cds);my $qh_get_transcript = $dbh2->prepare($get_transcript);

my $getOrg="select organism_version_id from organism_version where ucsc_db=?";   #get current organisms list 
my $qh_orgvlist= $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;
#*************************************
$basedir=$opt_o if($opt_o);$basedir=~s/\/\s*$//;
if(!(-d $basedir)){mkdir $basedir, 0777;}
#### this is where I need to set the loop to collect entries from the database
if(!(-d "$basedir/$oversion")){mkdir "$basedir/$oversion", 0777;}
$qh_lastUpdate->execute($gene_prediction);
if($qh_lastUpdate->rows>0){($current_version)=$qh_lastUpdate->fetchrow_array();}
$qh_orgvlist->execute($oversion);($orgvid)=$qh_orgvlist->fetchrow_array();
if($eschema eq $current_version){ print "No new updates ($eschema eq $current_version) \n";}
else
 { 
  if(!(-d "$basedir/$oversion/$eschema")){mkdir "$basedir/$oversion/$eschema", 0777;}
  $qh_getNewdbId->execute($gene_prediction);
  if($qh_getNewdbId->rows<=0){
     $qh_insertNewdb->execute($gene_prediction);
     $qh_getNewdbId->execute($gene_prediction);
     ($gene_prediction_id)=$qh_getNewdbId->fetchrow_array();
     $qh_deleteSchema->execute($gene_prediction); 
     $qh_insertExternaldb->execute($gene_prediction_id,$gene_prediction,$eschema);
  }
  else{
    ($gene_prediction_id)=$qh_getNewdbId->fetchrow_array();
    if($current_version eq ""){$qh_insertExternaldb->execute($gene_prediction_id,$gene_prediction,$eschema);}
    elsif($current_version ne $eschema){$qh_Update_schema->execute($eschema,$gene_prediction);}
  }$more=1; my @cds=();$rows="";@txs=();@exons=();
 #load cds and exons into memory
 $qh_get_cds->execute();
 if($qh_get_cds->rows>0){while(@row=$qh_get_cds->fetchrow_array()){push(@cds,join(",",@row));}}
 $qh_get_exons->execute();
 if($qh_get_exons->rows>0){while(@row=$qh_get_exons->fetchrow_array()){push(@exons,join(",",@row));}}
 $qh_get_transcript->execute();
 open(LOG,">$eschema.txt");$tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Starting load process :  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
 #check when this file was last updated and compare with what we have in the database
 $main_file="$basedir/$oversion/$eschema/$gene_prediction.txt"; $i=0;  %gene_list=();
 $gene_filename="$basedir/$oversion/$eschema/$gene_prediction-genesDesc.txt";
 #get the list of transcript, gene, exons ids of this organism version block by block
 open(OUT,">$main_file");open(GEN,">$gene_filename");
 if(OUT){
    print OUT "name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\t";
    print OUT "exonFrames\texonIds\tname2\n"; 
    print GEN "gene\tbiotype\torganism_version_id\tdescription\n";
    %gene_list=();
    if($qh_get_transcript->rows>0){$i=0;
       while(($name,$chrom,$txStart,$txEnd,$strand,$name2,$biotype,$desc)= $qh_get_transcript->fetchrow_array()){
          @tokens=grep(/$name/,@exons);$cdsStart=0;$cdsEnd=0;$exonCount=0;$exonStarts="";$exonEnds="";
          $exonFrames="";$exonsIds="";
          if(@tokens>0){
            while(@tokens>0){$line=shift(@tokens);($tx_name,$exstart,$exend,$frame,$exonid)=split(",",$line);
               if($name eq $tx_name){$exonStarts.="$exstart,";$exonEnds.="$exend,";$exonFrames.="$frame,";$exonsIds.="$exonid,";}
             }
          }# now get the CDS cordinates of this transcript
         @cds_cord=grep(/$name/,@cds);
         if(@cds_cord>0){while(@cds_cord>0){$line=shift(@cds_cord);($tx_name,$cdsStart,$cdsEnd)=split(",",$line);}} 
         $gene_list{"$name2"}=$desc;
         print OUT "$name\t$chrom\t$strand\t$txStart\t$txEnd\t$cdsStart\t$cdsEnd\t$exonCount\t";
         print OUT "$exonStarts\t$exonEnds\t$exonFrames\t$exonsIds\t$name2\n";++$i;
         print GEN "$name2\t$biotype\t$orgvid\t$desc\n";
         print LOG " $i transcripts proccessed\n" if($i%10000==0);
       }
    }close(OUT); close(GEN);
    #LOAD gene table
    ###################################################################
    #get the number of lines in the file
   $linecont = `wc -l $gene_filename`;chomp($linecont);
   if($linecont=~/^(\d+)\s$gene_filename/){$linecont=$1;
      if($linecont>0){
         $qh_drop_ens_temp->execute()or die "bad query ".mysql_error()."\n";
         $qh_create_ens_temp->execute()or die "bad query ".mysql_error()."\n";
         $qh_load_this->execute("$gene_filename") or die "bad query ".mysql_error()."\n";
         $qh_get_ens_genes_rowcount->execute();($rowcount)=$qh_get_ens_genes_rowcount->fetchrow_array();
         --$linecont;
         if($linecont==$rowcount){$qh_insert_ensgene->execute();
           print LOG "$gene_filename loaded: $rowcount of $linecont lines were loaded\n";}
         else{print LOG "$gene_filename Load failed:only $rowcount of $linecont lines were loaded\n";}
     }
   } 
  ###################################################################
  #load annotations into the database
   $dir="$basedir/$oversion/$eschema";
   $command="perl $load_annotations -d $dir -f $main_file -a $gene_prediction -v $oversion -s 1 -e 1";
   system($command); 
 }
 $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
       ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
 close(LOG);
}
print "Porgram complete\n";
exit(0);

