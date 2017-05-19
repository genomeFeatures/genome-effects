#!/usr/bin/perl

########################################################################
## This script collects all the miRNA annotations
# http://www.ebi.ac.uk/enright-srv/microcosm/htdocs/targets/v5/
#http://www.ebi.ac.uk/enright-srv/microcosm/htdocs/targets/v5/info.html
## from http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/download.pl
#  and load them into graber_transcriptdb.miRNAs table
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
#    a tab-delimited file with all the miRNA and their cordinates 
#    on the genome
q{ Documentation from EBI
About MicroCosm Targets

MicroCosm is a web resource developed by the Enright Lab at the EMBL-EBI containing computationally predicted targets for microRNAs across many species. The miRNA sequences are obtained from the miRNA Registry and most genomic sequence from EnsEMBL. We aim to provide the most up-to-date and accurate predictions of miRNA targets and hence this resource will be updated regularly. 
Finding miRNA Target Sites

We currently use the miRanda algorithm to identify potential binding sites for a given miRNA in genomic sequences. The current version uses dynamic programing alignment to identify highly complementary sites which are scored between 0 and 100, where 0 represents no complementarity and 100 complete complementary. The algorithm uses a weighted scoring system and rewards complementarity at the 5' end of the microRNA. Currently we demand strict complementarity at this so-called seed region in accordance with recent publications, this in practice means that we throw away alignments where more than one base in this region is not complementary to a target site. Target sites selected in this fashion are passed through the Vienna RNA folding routines in order to estimate their thermodynamic stability.

Finally every potential target site in a 3'UTR detected is checked to see whether the site is conserved in orthologous transcripts from other species. Currently we demand that for a site to be conserved it must be detected at the same position in a cross-species orthologous UTR alignment by an miRNA of the same family. Each target must be conserved in at least two species for inclusion in the database (With the exception of Human and Chimp whose sequences are too similar).

The entire process of assembling miRNAs, genomic sequence, cross species UTR alignments and miRanda analysis is performed in parallel on a high-performance compute cluster according to the protocol described below.

------------------------------------
http://epigenie.com/pursuit-of-the-moving-mirna-target/
1. miRanda
2. On TargetScan(S)
3. PicTar
Structure Before Sequence
Unlike the algorithms discussed above, thermodynamics based algorithms like DIANA-microT and RNAHybrid place more emphasis on target structure than seed complementarity.
1. DIANA
2. RNAHybrid
3. Arrival of the Departures
4. STarMir
5. RNA22
};
# Note: The cordinates 
#
# Data sources: http://www.ebi.ac.uk
# Date: March 27,2013 
# Note: Changing source to http://www.mirbase.org/ftp.shtml instead of from http://www.ebi.ac.uk/
#       the reason being mirbase has updates and organism versions
# link ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz
#**********************************************************************
## set default variables
$dbname="graber_transcriptdb";
$host="harlequin"; $user="lnh";$pass="lucie98";
use DBI;use vars qw ($opt_h  $opt_o);
use Getopt::Std;use LWP 5.64;
#use Time::localtime;
use POSIX;

$organism=""; # default organism
$inputfile="";$names="";$outputfile="";@geneList=();
my $browser = LWP::UserAgent->new;
   $browser->timeout(10);$browser->env_proxy;
   $browser->agent('Mozilla/5.0');

getopts('ho:');

if(($opt_h)||!($opt_o)){
    print <<HELP;

 This script collects all the miRNA annotations
 from http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/download.pl

Usage:

   perl download_miRNA.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  Output directory
  

Examples:
cmd: ./download_miRNA.pl -o /scratch/data/downloads/miRNAs

HELP
exit;
}
chdir($opt_o);$pwd=`pwd`;
#print "running prog from $pwd\n";
#exit(0);
$domain_temp="mirna.txt";
$url="http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/download.pl";
$response = $browser->get($url, ":content_file" => $domain_temp);
open(IN,"$domain_temp");if(IN){@filecontent=<IN> ;}close(IN);
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user,$pass);
my $getChr ="select chromosome_id,chromosome_name from chromosome";              # get current chromosomes list
my $getOrg="select organism_id,organism_name,organism_sc_name from organism";   #get current organisms list 
my $delete_this="delete from miRNAs where organism_id=?";
my $loadfile= "load data local infile ? into table miRNAs ignore 1 lines";
my $getOrgV= "select organism_version_id,ucsc_db from organism_version where organism_id=?";
my $getTx_id =" select t.transcript_id,ta.transcript_name from transcript_by_annotation ta, transcript t ";
   $getTx_id .=" where t.organism_version_id=? and ta.transcript_id=t.transcript_id ";

my $qh_chrlist    = $dbh->prepare($getChr)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_orglist    = $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_delete    = $dbh->prepare($delete_this)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_load    = $dbh->prepare($loadfile)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_orgv = $dbh->prepare($getOrgV)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_tx_id = $dbh->prepare($getTx_id)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_rowCount=$dbh->prepare("select count(*) from miRNAs where organism_id=?");

my $delete_log="delete from feature_load_log where table_name=?";
my $qh_delete_log = $dbh->prepare($delete_log)or die "Couldn't prepare statement: " . $dbh->errstr;


my $insert_log="insert into feature_load_log(table_name,segment_name,file_line_count,db_line_count,moddate)
                  values(?,?,?,?,concat(CURDATE(),':',CURTIME()))";
my $qh_insert_log = $dbh->prepare($insert_log)or die "Couldn't prepare statement: " . $dbh->errstr;

my(%chrmap,%orgmap);
  # get the list of all the chromosomes
$qh_chrlist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
while(($chr_id, $chr_name)=$qh_chrlist->fetchrow_array()){
   $chr_name=lc($chr_name);$chr_name=~s/^\s+//;$chr_name=~s/\s+$//;
   $chrmap{"$chr_name"}=$chr_id;
}
#get the list of all organisms
$qh_orglist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
while(($org_id,$common_name,$scientific_name)=$qh_orglist->fetchrow_array()){
   $common_name=~s/^\s+//;$common_name=~s/\s+$//;$scientific_name=~s/^\s+//;
   next if($org_id==58);$scientific_name=~s/\s+$//;$orgmap{"$scientific_name"}="$org_id";
}
#collect all the text zip files
$dom="";
while(@filecontent>0){
   $line=shift(@filecontent);chomp($line);
   if($line=~/<h1>miRBase Targets Download<\/h1>/){
      $line=~s/^\s+//; $line=~s/\s+$//;$dom=$line;$more=1;
      while(($more)&&(@filecontent>0)){
          $line=shift(@filecontent);chomp($line); $line=~s/^\s+//; $line=~s/\s+$//;
          if($line=~/<\/table>/){$dom.=$`;last;}
          else{$dom.=$line;}
      }
     last;
   }
}
@files=split("<TR",$dom);
if(!(-d temp)){mkdir("temp",0777)||print LOG "$!\n";}
#system("mv *.v5.* archive");
foreach my $file(@files){
    chomp($file); $organism_name="";$zip_link="";$zip_name="";$text_file="";
    @tds=split("<\/TD>",$file);
    while(@tds>0){
        $td=shift(@tds);
        if($td=~/<I>(.+)<\/I>/){$organism_name=$1;}
        if($td=~/HREF="(.+\/v5\/(.+))">TXT<\/A/i){$zip_link=$1;$zip_name=$2;}
    }next if(!exists($orgmap{"$organism_name"}));
    if(!(-f $zip_name)&&($zip_link ne "")){
       system("wget -q $zip_link");if(-f $zip_name){system("unzip -qu $zip_name");}  
    }if($zip_name=~/arch\.(.+)\.zip/){$text_file=$1;}
    #print "Processing $organism_name\n";
    if(-f $text_file){ #get the specified fields
        open(IN,"$text_file");
        if(IN){
            $more=1;$header="";$version="";
            ##source-version miRanda 3.0
            ##created on:2007-10-31
            while($more){ $header=<IN>;chomp($header);if($header=~/^##GROUP\s+/){$more=0;}}
            open(OUT,">$text_file-trim.txt");chomp($pwd);
            print OUT "organism_id\tfeature_name\tCHR\tchr_id\tSTART\tEND\tSTRAND\t";
            print OUT "transcript_local_id\tTRANSCRIPT_AC\tEXTERNAL_NAME\n";
            $filename="$pwd/$text_file-trim.txt";@fields=split("\t",$header);
            $org_id=$orgmap{"$organism_name"};%field_index=();$org_version_id=0;$current_v=0;
            #get current version of this organism
            $qh_orgv->execute($org_id) or die "Can't execute query: " . $dbh->errstr . "\n";
            while(($organism_version_id,$ucsc_db)=$qh_orgv->fetchrow_array()){
              $ucsc_db=~s/\s*$//;
              if($ucsc_db=~/(\d+)$/){
                 if($current_v==0){$current_v=$1;$org_version_id=$organism_version_id;}
                 else{if($current_v<$1){$current_v=$1;$org_version_id=$organism_version_id;}}
               }
            }$qh_tx_id->execute($org_version_id);%transcripMap=();
            while(($tx_id,$tx_name)=$qh_tx_id->fetchrow_array()){
                 $tx_name=lc($tx_name);$transcripMap{"$tx_name"}=$tx_id;}
            #now load the transcript map of this organism version
            $i=0;while($i<@fields){$field_index{$fields[$i]}=$i;++$i;}
            #print "$filename -- $organism_name has ".@fields." fields\n";
            while($line=<IN>){
                   chomp($line);
                   @field_content=split("\t",$line);
                   next if(@field_content!=@fields);
                   $name=$field_content[$field_index{"SEQ"}];$CHR=$field_content[$field_index{"CHR"}];
                   $START=$field_content[$field_index{"START"}];$END=$field_content[$field_index{"END"}];
                   #Adjust start and end cordinates : asumption is that EBI and Ensembl use the same
                   #cordinate system 1-base 
                   $START-=1; #adjust only Start to zero-base cordrinate
                   $STRAND=$field_content[$field_index{"STRAND"}];
                   $TRANSCRIPT_ID=$field_content[$field_index{"TRANSCRIPT_ID"}];
                   $EXTERNAL_NAME=$field_content[$field_index{"EXTERNAL_NAME"}];
                   $CHR=lc($CHR); $CHR="m" if($CHR eq "mt");$CHR="u" if($CHR eq "unkn");
                   $chr_id=0; $tx_id=0;$chr_id=$chrmap{"$CHR"} if(exists($chrmap{"$CHR"}));
                   $STRAND="-" if($STRAND ne "+");$tx_name=lc($TRANSCRIPT_ID);
                   if(exists($transcripMap{"$tx_name"})){$tx_id=$transcripMap{"$tx_name"};}
                   print OUT "$org_id\t$name\t$CHR\t$chr_id\t$START\t$END\t$STRAND\t";
                   print OUT "$tx_id\t$TRANSCRIPT_ID\t$EXTERNAL_NAME\n";
            }close(OUT);
            #NOW LOAD this file into the database
            $linecont = `wc -l $filename`;chomp($linecont); #load transcript 
            if($linecont=~/^(\d+)\s$filename/){$linecont=$1;--$linecont;
               $qh_delete->execute($org_id);
               $qh_load->execute("$filename"); ##########
               $qh_rowCount->execute($org_id);
               $qh_delete_log->execute($filename);
               ($rowcount)=$qh_rowCount->fetchrow_array();
               $update_table="$organism_name-miRNA";
               $qh_insert_log->execute($update_table,$filename,$linecont,$rowcount);
               #clean generated files
               system("rm $text_file");system("rm *.zip");
            }  
       }
    }  
 }
print "Program complete\n";
exit(0);

