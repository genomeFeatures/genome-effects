#!/usr/bin/perl
########################################################################
## This script collects all the protein domain annotations
## from http://pfam.sanger.ac.uk/protein/P26367?output=xml
#  and load them into graber_transcriptdb.pfamDomains table
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : February 2012
#
#  Input : a directory path where to process files
#  Output: download, unzip and parse each genome into
#    a tab-delimited file with all the protein domains and their aa cordinates 
#    on the genome
#
# Data sources: ftp://ftp.uniprot.org
#**********************************************************************
## set default variables
$dbname="graber_transcriptdb";
$host="demon"; $user="lnh";$pass="lucie98";
use DBI;use vars qw ($opt_h $opt_f $opt_o);
use Getopt::Std;use LWP 5.64;
use Time::localtime;

use XML::Simple;
use Data::Dumper;

# create object
$xml = new XML::Simple (KeyAttr=>[]);


$organism=""; # default organism
$inputfile="";$names="";$outputfile="";@geneList=();
my $browser = LWP::UserAgent->new;
   $browser->timeout(10);$browser->env_proxy;
   $browser->agent('Mozilla/5.0');

getopts('ho:f:');
if(($opt_h)||(!$opt_f && !$opt_o)){
    print <<HELP;

 This script collects all the protein domain annotations
 from http://pfam.sanger.ac.uk/protein/xxx?output=xml
  and load them into graber_transcriptdb.pfamDomains table

Usage:
   perl download_pfamDomains.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  Output directory
  

Examples:
cmd: ./download_pfamDomains.pl -o /scratch/data/downloads/proteinDomains

HELP
exit;
}
chdir($opt_o);
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);
my $getOrg="select distinct o.organism_id,organism_sc_name,organism_tax_id from organism o,";   #get current organisms list 
   $getOrg.=" proteinDomains p where o.organism_id=p.organism_id";
my $qh_orglist    = $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;
my $getOrg="select organism_id,organism_sc_name,organism_tax_id from organism";   #get protein list 

my $create_table="create table if not exists pfamDomains(
  organism_id SMALLINT default 0,
  uniprot_protein_id varchar(25) not null,
  pfams_id varchar(25) not null,
  pfam_name varchar(100),
  pfam_type varchar(25) not null,
  feature_aa_start smallint default 0,
  feature_aa_end smallint default 0,
  index(organism_id),index(uniprot_protein_id),
  index(pfam_id),index(pfam_type),FOREIGN KEY(organism_id) REFERENCES organism(organism_id)
)ENGINE=MyISAM";
my $qh_create_table= $dbh->prepare($create_table);

my $load_this="load data local infile ? into table pfamDomains";
my $qh_load_this = $dbh->prepare($load_this)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_get_pd=$dbh->prepare("select distinct uniprot_protein_id from proteinDomains where organism_id=? ");
   #and pfams like '%PF%'");
my $qh_get_pf=$dbh->prepare("select* from pfamDomains where organism_id=? and uniprot_protein_id=?");

my $get_rowcount="select count(*) as rowcount from pfamDomains where organism_id=?";
my $qh_get_rowcount = $dbh->prepare($get_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;

my $delete_pf="delete from pfamDomains where organism_id=?";
my $qh_delete_pf = $dbh->prepare($delete_pf)or die "Couldn't prepare statement: " . $dbh->errstr;

my $analyze="analyze table pfamDomains ";
my $qh_analyze = $dbh->prepare($analyze)or die "Couldn't prepare statement: " . $dbh->errstr;

my %files_hd=();
open(LOG,">pfam_log.log");
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Starting load process :  $mday/$mon/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
#collect a
#get the list of all organisms

$qh_orglist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
$url="http://pfam.sanger.ac.uk/protein/";$protein="";$sufix="?output=xml"; $count=0;
while(($org_id,$organism_sc_name,$organism_tax_id)=$qh_orglist->fetchrow_array()){
   $organism_sc_name=~s/^\s+//;$organism_sc_namee=~s/\s+$//;
   $organism_tax_id=~s/^\s+//;$organism_tax_id=~s/\s+$//; $filename="pfam_load.txt";
   print "Processing $org_id,$organism_sc_name,$organism_tax_id\n";
   next if($org_id !=38); #only load mouse
   open(OUT,">$filename");
   $qh_get_pd->execute($org_id);$line_count=0;
   if($qh_get_pd->rows>0){
      while(($pids)=$qh_get_pd->fetchrow_array()){
         @ids=split(";",$pids);
         while(@ids>0){$protein="";
             $protein=shift(@ids);$count+=1;$protein=~s/^\s+//;$protein=~s/\s+$//;
             if($protein){
                #check if this protein is already loaded
                $qh_get_pf->execute($org_id,$protein);
                next if($qh_get_pf->rows>0);
                $temp_url="$url$protein$sufix"; $xml_file="$protein$sufix";
                eval{
                   system("wget $temp_url -q");
                };print LOG  $@ if $@;
                if(-f $xml_file){
                   #print LOG "pROCESSING $xml_file\n";
                   @html_test=`grep "<html>" $xml_file`;
                   next if(@html_test>0);
                   $data = $xml->XMLin("$xml_file");
                   if(ref($data->{entry}->{matches}->{match}) ne 'ARRAY'){
                     eval{
                       $pfam=$data->{entry}->{matches}->{match};
                       $pfam_accession=$pfam->{accession};$pfam_name=$pfam->{id};
                       $pfam_type=$pfam->{type};$pfam_start=$pfam->{location}->{start};
                       $pfam_end=$pfam->{location}->{end};
                       print OUT "$org_id\t$protein\t$pfam_accession\t$pfam_name\t";
                       print OUT "$pfam_type\t$pfam_start\t$pfam_end\n";$line_count+=1;
                     };print LOG  $@ if $@;
                   }
                   else{
                   foreach $pfam(@{$data->{entry}->{matches}->{match}}){
                     eval{
                       $pfam_accession=$pfam->{accession};$pfam_name=$pfam->{id};
                       $pfam_type=$pfam->{type};$pfam_start=$pfam->{location}->{start};
                       $pfam_end=$pfam->{location}->{end};
                       print OUT "$org_id\t$protein\t$pfam_accession\t$pfam_name\t";
                       print OUT "$pfam_type\t$pfam_start\t$pfam_end\n";$line_count+=1;
                     };print LOG  $@ if $@;
                   }}
                  system("rm $opt_o/$xml_file");
                }
               if($count %200==0){print "$count processed\n"; sleep(3);}
             } #end of if($protein)
         } #end of while(@ids>0)
      } # end of while(($pids)
   } #end of $qh_get_pd->rows>0
  if(-f "$filename"){
     $qh_create_table->execute();$qh_delete_pf->execute($org_id);
     $qh_load_this->execute("$filename");$row_count=0;
     $qh_get_rowcount->execute($org_id);
     if($qh_get_rowcount->rows>0){($row_count)=$qh_get_rowcount->fetchrow_array();}
     print LOG "$organism_sc_name:$line_count lines total $row_count loaded\n";
  }
}
$qh_analyze->execute();$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
close(LOG); 
print "Porgram complete\n";
exit(0);

