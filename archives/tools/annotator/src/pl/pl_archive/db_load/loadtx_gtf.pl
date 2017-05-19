#!/usr/bin/perl
########################################################################
## This script parses Al's gtf files to generated genes
##  (gene_name, gene_id, biotype),inserts
## new chromosomes into the database, and generates
#  the the genePrediction file with the associated gene
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : February 2013
#
#  Input : gets a list of all the ensembl schema currently in our database
#  Output: a tab-delimited file with all the genes and their descriptions
#
## set default variables
use DBI;
$dbname="graber_transcriptdb";
$host="harlequin"; $user="lnh";$pass="lucie98";
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

my $getThisChr="select chromosome_id from chromosome where chromosome_name=?";
my $qh_gtchrid     = $dbh->prepare($getThisChr)or die "Couldn't prepare statement: " . $dbh->errstr;
my $insertChr    ="insert ignore into chromosome(chromosome_name) values(?)";
my $qh_insertchr   = $dbh->prepare($insertChr)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_getPredid= $dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name=?");

my $drop_ens_temp = "drop temporary table if exists ensembl_genes_temp";
my $create_ens_temp = "create temporary table ensembl_genes_temp(gene_prediction_id tinyint unsigned,
                      transcript_id varchar(100),transcript_name varchar(255),protein_id varchar(100),
                      index(gene_prediction_id),index(transcript_id),index(protein_id))";
my $load_this = "load data local infile ? into table ensembl_genes_temp ignore 1 lines";
my $delete_this= "delete from cc_founders_transcripts where gene_prediction_id=?";
my $get_ens_genes_rowcount="select count(*) as rowcount from ensembl_genes_temp ";
my $insert_ensgene = "insert into cc_founders_transcripts(gene_prediction_id,transcript_id,transcript_name,protein_id) ";
   $insert_ensgene.="select distinct gene_prediction_id,transcript_id,transcript_name,protein_id from ensembl_genes_temp";

my $qh_drop_ens_temp = $dbh->prepare($drop_ens_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_create_ens_temp = $dbh->prepare($create_ens_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_load_this = $dbh->prepare($load_this)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_get_ens_genes_rowcount = $dbh->prepare($get_ens_genes_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_delete_this = $dbh->prepare($delete_this)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insert_ensgene = $dbh->prepare($insert_ensgene)or die "Couldn't prepare statement: " . $dbh->errstr;


## get the gtf file
$oversion="mm9";%chrmap=();
@strains=("129S1","A_J","C57BL","CAST","NOD","NZO","PWK","WSB","129P2","129S5","AKR","BALB","C3H","CBA","DBA","LP_J");
foreach $strain(@strains){
  $gtffile="$strain/$strain.gtf"; $dir=$strain;$strain=($strain eq "A_J")? "AJ":$strain; 
  $strain=($strain eq "LP_J")? "LPJ":$strain; 
  open(LOG,">>$strain/$oversion-$oversion-$strain"."Transcripts-log.txt");
  $founder_genes="$dir/$oversion-$strain"."TranscriptDesc.txt";%txmap=();
  $annotation="$oversion-$strain"."Gene";$predid=0;
  $qh_getPredid->execute($annotation);
  if($qh_getPredid->rows>0){($predid)=$qh_getPredid->fetchrow_array();}
  print "Processing $annotation \t$predid \n";
  if(-f $gtffile){
    @genes=`cut -f 9 $gtffile|sort|uniq`;%proteinmap=();%tx_map=();
    while(@genes>0){
       $gene=shift(@genes);chomp($gene);@fields=split(";",$gene); @txs=grep(/transcript_id/,@fields); 
       @proteins=grep(/protein_id/,@fields);
       if($txs[0]=~/transcript_id\s+"(.+)"/){
            $tx_id=$1;$tx_id=~s/\s+//g;$tx_map{"$tx_id"}=$gene;
            if($proteins[0]=~/protein_id\s+"(.+)"/){$proteinmap{"$tx_id"}=$1;}
        }
    }
    open(GE,">$founder_genes"); print GE "gene_prediction_id\ttranscript_id\ttranscript_name\tprotein_id\n";
    while(($tx_id,$tx_line)=each(%tx_map)){
         @row=split(";",$tx_line); $protein_id="";$transcript_name="";
         foreach my $line(@row){
            $transcript_name=$1 if($line=~/transcript_name\s+"(.+)"/);
            #$protein_id=$1 if($line=~/protein_id\s+"(.+)"/);
         }
        $protein_id = $proteinmap{"$tx_id"};
        print GE "$predid\t$tx_id\t$transcript_name\t$protein_id\n" if($transcript_name ne "");
    }
  }
 print "Loading  $founder_genes\n";#now load genes into 
 if(-f $founder_genes){
    $linecont = `wc -l $founder_genes`;chomp($linecont);
    if($linecont=~/^(\d+)\s$founder_genes/){$linecont=$1;
       if($linecont>0){ print "Loading $founder_genes \n"; --$linecont;
          $qh_drop_ens_temp->execute();$qh_create_ens_temp->execute();$qh_load_this->execute($founder_genes);
          $qh_get_ens_genes_rowcount->execute();($rowcount)=$qh_get_ens_genes_rowcount->fetchrow_array();
          if($rowcount==$linecont){
             $qh_delete_this->execute($predid); $qh_insert_ensgene->execute();
             print LOG "$founder_genes loaded. $rowcount of $linecont lines were loaded\n";
           }
          else{ print LOG "$founder_genes Load failed:only $rowcount of $linecont lines were loaded\n";}
        }
    }
  }
}
print "complete\n";


