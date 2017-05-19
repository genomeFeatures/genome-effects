#!/usr/bin/perl
########################################################################
## This script collects ensembl proteins for a given organism
## from ftp://ftp.ensembl.org/pub/release-66/fasta/mus_musculus/pep/Mus_musculus.NCBIM37.66.pep.all.fa.gz
#  and load them into graber_transcriptdb.proteinEnsembl table
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
#    a fasta file with all the protein 
#
# Data sources: ftp://ftp.ensembl.org
#**********************************************************************
## set default variables
$dbname="graber_transcriptdb";
$host="demon"; $user="lnh";$pass="lucie98";
use DBI;use vars qw ($opt_h $opt_f $opt_o);
use Getopt::Std;use LWP 5.64;
#use Time::localtime;
use POSIX;

$organism=""; # default organism
$inputfile="";$names="";$outputfile="";@geneList=();
my $browser = LWP::UserAgent->new;
   $browser->timeout(10);$browser->env_proxy;
   $browser->agent('Mozilla/5.0');

getopts('ho:f:');
if(($opt_h)||(!$opt_f && !$opt_o)){
    print <<HELP;

 This script collects ensembl proteins for a given organism
 from ftp://ftp.ensembl.org/pub/release-66/fasta/mus_musculus/pep/Mus_musculus.NCBIM37.66.pep.all.fa.gz
 and load them into graber_transcriptdb.proteinEnsembl table

Usage:
   perl download_ensemblProtein.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  Output directory
  

Examples:
cmd: ./download_ensemblProtein.pl -o /scratch/data/downloads/ensemblProteins

HELP
exit;
}
chdir($opt_o);
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);
my $getOrg="select organism_id,organism_sc_name,organism_tax_id from organism";   #get current organisms list 
my $qh_orglist    = $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;

my $getOrgv="select organism_version_id,coordinate_build from organism_version";   #get current organisms list 
my $qh_orgvlist    = $dbh->prepare($getOrgv)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_dropSeq = $dbh->prepare("drop table if exists proteinEnsembl");
# print OUT "$orgv_id\t$chr\t$tx_start\t$tx_end\t$protein\t$transcript\t$gene\t$seq";
my $create_Seq="create table if not exists proteinEnsembl(
    organism_version_id SMALLINT default 0,
    chromosome varchar(15),transcript_start int unsigned default 0,
    transcript_end int unsigned default 0,
    protein_id varchar(25) not null,
    transcript varchar(25) not null,
    gene varchar(25) not null,
    sequence text,
    FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
    index(organism_version_id),index(protein_id),index(transcript),index(gene)
  )ENGINE=MyISAM";
$qh_create_Seq= $dbh->prepare($create_Seq);
my $load_seq="load data local infile ? into table proteinEnsembl";
my $qh_load_seq = $dbh->prepare($load_seq)or die "Couldn't prepare statement: " . $dbh->errstr;
my $get_seq_rowcount="select count(*) as rowcount from proteinEnsembl ";
my $qh_get_seq_rowcount = $dbh->prepare($get_seq_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;

my $delete_seq_log="delete from feature_load_log where table_name='proteinEnsembl'";
my $qh_delete_seq_log = $dbh->prepare($delete_seq_log)or die "Couldn't prepare statement: " . $dbh->errstr;

my $insert_seq_log="insert into feature_load_log(table_name,segment_name,file_line_count,db_line_count,moddate)
                  values('proteinEnsembl',?,?,?,concat(CURDATE(),':',CURTIME()))";
my $qh_insert_seq_log = $dbh->prepare($insert_seq_log)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_analyze_seq = $dbh->prepare("analyze table proteinEnsembl")or die "Couldn't prepare statement: " . $dbh->errstr;

#collect a
$zip_link="ftp://ftp.ensembl.org/pub/release-66/fasta/mus_musculus/pep/Mus_musculus.NCBIM37.66.pep.all.fa.gz";
$zip_name="Mus_musculus.NCBIM37.66.pep.all.fa.gz";$text_file="Mus_musculus.NCBIM37.66.pep.all.fa"; 
#if(-f $zip_name){system("mv $zip_name archive");}
#if(-f $text_file){system("mv $text_file archive");}
if( !(-f $text_file)&&!(-f $zip_name)){system("wget $zip_link");}
if((-f $zip_name)&& !(-f $text_file)){system("gunzip $zip_name");}
$more=1; my @protein=();$rows="";$sequence_file="$opt_o/sequence_load.txt";
if(-f $text_file){ #get the specified fields
   $abuild="NCBI37";$chr="";$tx_start=0;$tx_end=0;$gene="";$transcript="";$protein_id="";
   open(FH,"$text_file"); $lcount=0;open(OUT,">sequence_load.txt");
   $first=<FH>;chomp($first);if($first=~/chromosome:(.+):(.+):(\d+):(\d+):/){$abuild=$1;$chr=$2;$tx_start=$3;$tx_end=$4;}
   if($abuild=~/NCBI\w(\d+)/){$abuild="NCBI$1";}
   $qh_orgvlist->execute();
   if($qh_orgvlist->rows>0){
      while(($orgv_id,$build)=$qh_orgvlist->fetchrow_array()){last if($build=~/$abuild/);}
    }
   close(FH);
   open(FH,"$text_file");
   @protein=();$prev_line="";
   while(<FH>){ chomp($_); 
     if(!($_=~/^>/)){push(@protein,$_);}
     else{
        if($prev_line eq ""){$prev_line=$_;}
        else{ # we have a sequence
             $seq= join "",@protein;$seq=~s/\s*//g;
             if($prev_line=~/chromosome:(.+):(.+):(\d+):(\d+):/){$abuild=$1;$chr=$2;$tx_start=$3;$tx_end=$4;}
             if($prev_line=~/\sgene:(.+)\s+transcript:/){$gene=$1;} 
             if($prev_line=~/\stranscript:(.+)\s+gene_biotype:/){$transcript=$1;}
             if($prev_line=~/^>(.+)\s+pep:/){$protein_id=$1;}
             print OUT "$orgv_id\t$chr\t$tx_start\t$tx_end\t$protein_id\t$transcript\t$gene\t$seq\n";
             $prev_line=$_;@protein=();
        }
     }
  }
}
close(OUT);
if(-f "sequence_load.txt"){
  $qh_dropSeq->execute();$qh_create_Seq->execute();
  open(IN,"sequence_load.txt");$rows="";$index=0;
  while(<IN>){ chomp($_);$rows.="$_\n"; ++$index;
     if($index%50000==0){
         $filename="sql_seq_load_segment.txt";open(OUT,">$filename");
         if(OUT){print OUT "$rows";$rows="";}close(OUT);
         if(-f "$filename"){
              $filename="$opt_o/$filename"; #get the number of lines in the file
              $linecont = `wc -l $filename`;chomp($linecont);
              if($linecont=~/^(\d+)\s$filename/){$linecont=$1;
                 if($linecont>0){ $qh_load_seq->execute($filename)or die "Couldn't prepare statement: " . $dbh->errstr;
                 }
               } 
          } 
      }   
   }
   if($rows ne ""){
       $filename="sql_seq_load_segment.txt";open(OUT,">$filename");
       if(OUT){print OUT "$rows";$rows="";}close(OUT);
         if(-f "$filename"){$filename="$opt_o/$filename";
              $linecont = `wc -l $filename`;chomp($linecont);
              if($linecont=~/^(\d+)\s$filename/){$linecont=$1;
                 if($linecont>0){ 
                    $qh_load_seq->execute("$filename")or die "Couldn't prepare statement: ".$dbh->errstr;
                 }
             } 
       } 
    }
   #get the total number of lines of this segment that were inserted
   $qh_get_seq_rowcount->execute();
  ($rowcount)= $qh_get_seq_rowcount->fetchrow_array();
   #now insert this info into the log
   $filename="$opt_o/sequence_load.txt";
   $linecont = `wc -l $filename`;chomp($linecont);
   if($linecont=~/^(\d+)\s$filename/){$linecont=$1;
      $qh_delete_seq_log->execute();
      $qh_insert_seq_log->execute("$filename",$linecont,$rowcount);
    }
   close(IN);
 }
$qh_analyze_seq->execute();
print "Porgram complete\n";
exit(0);

