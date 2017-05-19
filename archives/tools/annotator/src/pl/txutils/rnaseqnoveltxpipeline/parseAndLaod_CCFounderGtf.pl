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
                      gene_id varchar(100),gene_name varchar(255),biotype varchar(100),index(gene_prediction_id))";
my $load_this = "load data local infile ? into table ensembl_genes_temp ignore 1 lines";
my $delete_this= "delete from cc_founders_genes where gene_prediction_id=?";
my $get_ens_genes_rowcount="select count(*) as rowcount from ensembl_genes_temp ";
my $insert_ensgene = "insert ignore into cc_founders_genes(gene_prediction_id,gene_id,gene_name,biotype) ";
   $insert_ensgene.="select distinct gene_prediction_id,gene_id,gene_name,biotype from ensembl_genes_temp";

my $qh_drop_ens_temp = $dbh->prepare($drop_ens_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_create_ens_temp = $dbh->prepare($create_ens_temp)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_load_this = $dbh->prepare($load_this)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_get_ens_genes_rowcount = $dbh->prepare($get_ens_genes_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_delete_this = $dbh->prepare($delete_this)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_insert_ensgene = $dbh->prepare($insert_ensgene)or die "Couldn't prepare statement: " . $dbh->errstr;


## get the gtf file
$oversion="mm9";%chrmap=();
#@strains=("129S1","A_J","C57BL","CAST","NOD","NZO","PWK","WSB","129P2","129S5","AKR","BALB","C3H","CBA","DBA","LP_J");
@strains=("129S1","A_J","CAST","NOD","NZO","PWK","WSB","C57BL");
foreach $strain(@strains){
  $gtffile="$strain/$strain.gtf"; $txt_file="$strain/$strain.txt";$dir=$strain;
  $strain=($strain eq "A_J")? "AJ":$strain; $strain=($strain eq "LP_J")? "LPJ":$strain; 
  $db_file="$dir/$oversion-$strain"."Gene.txt";
  open(LOG,">>$strain/$oversion-$oversion-$strain"."Gene-log.txt");
  $founder_genes="$dir/$oversion-$strain"."GeneDesc.txt";%txmap=();
  $annotation="$oversion-$strain"."Gene";$predid=0;
  $qh_getPredid->execute($annotation);
  if($qh_getPredid->rows>0){($predid)=$qh_getPredid->fetchrow_array();}
  print "Processing $annotation \t$predid \n";
  if(-f $txt_file){
     @chromosomes=`cut -f 2 $txt_file|sort|uniq`;
     foreach my $chr(@chromosomes){
       $chr="M" if($chr=~/MT/);$chr=~s/\s+//g;$chr=~s/chr//i;$chr=~s/ch//i;$qh_gtchrid->execute($chr);
       if($qh_gtchrid->rows<=0){$qh_insertchr->execute($chr);$qh_gtchrid->execute($chr);}
       ($chr_id)=$qh_gtchrid->fetchrow_array();$chrmap{$chr}=$chr_id;
     } #now get the gene colums 
    @genes=`cut -f 9 $gtffile|sort|uniq`;%genemap=();%tx_map=();
    while(@genes>0){
       $gene=shift(@genes);@fields=split(";",$gene); @txs=grep(/transcript_id/,@fields); 
       if($txs[0]=~/transcript_id\s+"(.+)"/){$tx_id=$1;$tx_id=~s/\s+//g;$tx_map{"$tx_id"}=$gene;}
    }
    open(IN,"$txt_file"); open(OUT,">$db_file");open(GE,">$founder_genes");
    print GE "gene_prediction_id\tgene_id\tgene_name\tbiotype\n";
    print OUT "name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tname2\n";
    while(<IN>){chomp($_); $tx_line=$_;
         @row=split(/\t/,$tx_line); $tx_id=$row[0]; $tx_id=~s/\s+//g; 
         $gene_id="";$gene_name="";$biotype=""; $gene=$tx_map{"$tx_id"};$tx="";
         if($gene=~/gene_id\s+"(.+)";\s+transcript_id\s+"(.+)";\s+gene_name\s+"(.+)";\s+gene_biotype\s+"(.+)"/){
            $gene_id=$1;$tx=$2; $gene_name=$3;$biotype=$4;
            $gene_id=$` if($gene_id=~/"/);$gene_name=$` if($gene_name=~/"/);$biotype=$` if($biotype=~/"/);
            $tx=$` if($tx=~/"/);
            $genemap{"$gene_id"}="$gene_name=$biotype"; print OUT "$tx_line\t$gene_id,$gene_name\n";
         }
    }print "$db_file generated\n"; #display gene file
    while(($gene_id,$data)=each(%genemap)){
      ($gene_name,$biotype)=split(/=/,$data);print GE "$predid\t$gene_id\t$gene_name\t$biotype\n";}
 }
 #now load genes into 
 if(-f $founder_genes){
    $linecont = `wc -l $founder_genes`;chomp($linecont);
    if($linecont=~/^(\d+)\s$founder_genes/){$linecont=$1;
       if($linecont>0){ print "Loading $founder_genes \n"; --$linecont;
          $qh_drop_ens_temp->execute();$qh_create_ens_temp->execute();$qh_load_this->execute($founder_genes);
          $qh_get_ens_genes_rowcount->execute();($rowcount)=$qh_get_ens_genes_rowcount->fetchrow_array();
          if($rowcount==$linecont){
             #$qh_delete_this->execute($predid); 
             $qh_insert_ensgene->execute();
             print LOG "$founder_genes loaded. $rowcount of $linecont lines were loaded\n";
           }
         else{ print LOG "$founder_genes Load failed:only $rowcount of $linecont lines were loaded\n";}
        }
    }
 }
 #now load the transcript
 $cmd="perl pl/graber_txdb/load_annotations.pl -d $dir -f $db_file -a $annotation -v $oversion -s 0 -e 1";
 system($cmd);
 
}
print "complete\n";


