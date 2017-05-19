#!/usr/bin/perl
########################################################################
## This script collects all the gene descriptions for a given organism version
## from Ensembl and load them into graber_transcriptdb.ensemblGenes
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
This script collects all the gene descriptions for a given organism version
from Ensembl and load them into graber_transcriptdb.ensemblGenes

Usage:
   perl download_ensGeneDesc.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  Output directory

Examples:
cmd: ./download_ensGeneDesc.pl -o $basedir

HELP
exit;
}
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

my $ens_schema="select db_name, schema_version from external_db where db_name like '%ensGene'";
  $qh_schema=$dbh->prepare($ens_schema);

my $drop_ens_temp = "drop temporary table if exists ensembl_genes_temp";
my $create_ens_temp = "create temporary table ensembl_genes_temp(gene_name varchar(255),
          biotype varchar(100),organism_version_id smallint default 0,description text,index(organism_version_id))";
my $load_this = "load data local infile ? into table ensembl_genes_temp ignore 1 lines";
my $delete_this= "delete from ensembl_genes where organism_version_id=?";
my $get_ens_genes_rowcount="select count(*) as rowcount from ensembl_genes_temp ";
my $insert_ensgene = "insert into ensembl_genes(gene_name,biotype,organism_version_id,description) ";
   $insert_ensgene.="select distinct gene_name,biotype,organism_version_id,description from ensembl_genes_temp";

my $qh_drop_ens_temp = $dbh->prepare($drop_ens_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_create_ens_temp = $dbh->prepare($create_ens_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_load_this = $dbh->prepare($load_this)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_get_ens_genes_rowcount = $dbh->prepare($get_ens_genes_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_delete_this = $dbh->prepare($delete_this)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insert_ensgene = $dbh->prepare($insert_ensgene)or die "Couldn't prepare statement: " . $dbh->errstr;
my $getOrg="select organism_version_id from organism_version where ucsc_db=?";   #get current organisms list 
my $qh_orgvlist= $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;
## Get data to load into cgd_transcript table *
$get_genes="select distinct ge.stable_id ,ge.biotype,ge.description 
                 from gene ge,seq_region s ,coord_system c
                 where s.coord_system_id=c.coord_system_id and c.rank=1 
                 and  s.seq_region_id=ge.seq_region_id";

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
   ($orgvid)=$qh_orgvlist->fetchrow_array();
   print "Processing $orgvid,$oversion,$gene_prediction, $eschema\n";
   if(!(-d $basedir)){mkdir $basedir, 0777;}if(!(-d "$basedir/$oversion")){mkdir "$basedir/$oversion", 0777;}
   if(!(-d "$basedir/$oversion/$eschema")){mkdir "$basedir/$oversion/$eschema", 0777;}
   open(LOG,">$basedir/$oversion/$eschema/genes-log.txt");$tm = localtime;
   $filename="$basedir/$oversion/$eschema/$gene_prediction-genesDesc.txt";
   $qh_get_genes = $dbh2->prepare($get_genes);
   open(OUT,">$filename"); #data will be loaded into ensembl_genes table
   if(OUT){
      print OUT "gene\tbiotype\torganism_version_id\tdescription\n";
      $qh_get_genes->execute();
      my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
      print LOG "\n*************************************************************\n";
      print LOG "Starting load process :  $mday/$mon/$yday @ $hour:$min:$sec \n";
      print LOG "\n*************************************************************\n";
      print "generating genes file $filename\n";open(OUT,">$filename");
     if($qh_get_genes->rows>0){
        while(($name,$biotype,$desc)= $qh_get_genes->fetchrow_array()){
            print OUT "$name\t$biotype\t$orgvid\t$desc\n";
        }
     }
   }
   close(OUT);
   ###################################################################
    #get the number of lines in the file
   $linecont = `wc -l $filename`;chomp($linecont);
   if($linecont=~/^(\d+)\s$filename/){$linecont=$1;
      if($linecont>0){
         $qh_drop_ens_temp->execute()or die "bad query ".mysql_error()."\n";
         $qh_create_ens_temp->execute()or die "bad query ".mysql_error()."\n";
         $qh_load_this->execute("$filename") or die "bad query ".mysql_error()."\n";
         $qh_get_ens_genes_rowcount->execute();($rowcount)=$qh_get_ens_genes_rowcount->fetchrow_array();
         --$linecont;
         if($linecont==$rowcount){$qh_delete_this->execute($orgvid);$qh_insert_ensgene->execute();
           print "Processed $linecont lines  and loaded $rowcount lines\n";}
         else{ 
             print LOG "$filename Load failed:only $rowcount of $linecont lines were loaded\n";
               #the bulk load did not go well
         }
     }
   } 
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

