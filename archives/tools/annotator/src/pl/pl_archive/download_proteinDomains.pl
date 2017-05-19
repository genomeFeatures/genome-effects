#!/usr/bin/perl
########################################################################
## This script collects all the protein  annotations
## from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
## The UniProt Knowledgebase has been created from Swiss-Prot, TrEMBL and PIR-PSD.
## It consists of two parts, one containing fully manually annotated records
## and another one with computationally analysed records awaiting full manual
## annotation. The two sections will be referred to as the Swiss-Prot 
## Knowledgebase and TrEMBL Protein Database, respectively. PIR-PSD release 80.0 
## of 31-Dec-2004 has been fully integrated into these sections. This was the 
## last release of PIR-PSD. 

##Note: The file contains the four-weekly updates of the UniProt Knowledgebase,
#       consisting of UniProtKB/Swiss-Prot (fully annotated curated entries) 

#  and load them into graber_transcriptdb.proteinDomains table
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
$host="harlequin"; $user="lnh";$pass="lucie98";
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

 This script collects all the protein domain annotations
 from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
 and load them into graber_transcriptdb.proteinDomains table

Usage:
   perl download_proteinDomains.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  Output directory
  

Examples:
cmd: ./download_proteinDomains.pl -o /scratch/data/downloads/proteinDomains

HELP
exit;
}
chdir($opt_o);
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);
my $getOrg="select organism_id,organism_sc_name,organism_tax_id from organism";   #get current organisms list 
my $qh_orglist    = $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;

my $qh_dropSeq = $dbh->prepare("drop table if exists proteinSequence");
my $create_Seq="create table if not exists proteinSequence(
    organism_id SMALLINT default 0,
    uniprot_protein_id varchar(255) not null,
    isoform varchar(50),
    sequence_len smallint default 0,
    sequence text,
    FOREIGN KEY(organism_id) REFERENCES organism(organism_id),
    index(organism_id),index(uniprot_protein_id)
  )ENGINE=MyISAM";
$qh_create_Seq= $dbh->prepare($create_Seq);
my $load_seq="load data local infile ? into table proteinSequence";
my $qh_load_seq = $dbh->prepare($load_seq)or die "Couldn't prepare statement: " . $dbh->errstr;
my $get_seq_rowcount="select count(*) as rowcount from proteinSequence ";
my $qh_get_seq_rowcount = $dbh->prepare($get_seq_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;

my $delete_seq_log="delete from feature_load_log where table_name='proteinSequence'";
my $qh_delete_seq_log = $dbh->prepare($delete_seq_log)or die "Couldn't prepare statement: " . $dbh->errstr;


my $insert_seq_log="insert into feature_load_log(table_name,segment_name,file_line_count,db_line_count,moddate)
                  values('proteinSequence',?,?,?,concat(CURDATE(),':',CURTIME()))";
my $qh_insert_seq_log = $dbh->prepare($insert_seq_log)or die "Couldn't prepare statement: " . $dbh->errstr;



my $qh_droptable = $dbh->prepare("drop table if exists  proteinDomains");
my $create_table="create table if not exists  proteinDomains(
  organism_id SMALLINT default 0,
  uniprot_protein_id varchar(255) not null,
  uniprot_protein_name varchar(255) not null,
  protein_evidence varchar(255) not null,
  pfams varchar(255) not null,
  common_gene_name varchar(255) not null,
  other_transcript_name varchar(255) not null,
  other_protein_id varchar(255) not null,
  other_gene_name  varchar(255) not null,
  feature_type varchar(50),
  feature_aa_start smallint default 0,
  feature_aa_end smallint default 0,
  feature_name varchar(100),
  index(organism_id),index(feature_type),index(uniprot_protein_id),index(other_transcript_name),
  index(other_protein_id),index(other_gene_name),index(pfams),index(common_gene_name),
  FOREIGN KEY(organism_id) REFERENCES organism(organism_id)
)ENGINE=MyISAM";
my $qh_create_table= $dbh->prepare($create_table);
my $load_this="load data local infile ? into table  proteinDomains";
my $qh_load_this = $dbh->prepare($load_this)or die "Couldn't prepare statement: " . $dbh->errstr;

my $get_rowcount="select count(*) as rowcount from  proteinDomains ";
my $qh_get_rowcount = $dbh->prepare($get_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;

my $delete_log="delete from feature_load_log where table_name=' proteinDomains'";
my $qh_delete_log = $dbh->prepare($delete_log)or die "Couldn't prepare statement: " . $dbh->errstr;


my $insert_log="insert into feature_load_log(table_name,segment_name,file_line_count,db_line_count,moddate)
                  values(' proteinDomains',?,?,?,concat(CURDATE(),':',CURTIME()))";
my $qh_insert_log = $dbh->prepare($insert_log)or die "Couldn't prepare statement: " . $dbh->errstr;

my $analyze="analyze table  proteinDomains ";
my $qh_analyze = $dbh->prepare($analyze)or die "Couldn't prepare statement: " . $dbh->errstr;
my $qh_analyze_seq = $dbh->prepare("analyze table proteinSequence");


my %orgmap=();
#get the list of all organisms
$qh_orglist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
while(($org_id,$organism_sc_name,$organism_tax_id)=$qh_orglist->fetchrow_array()){
   $organism_sc_name=~s/^\s+//;$organism_sc_namee=~s/\s+$//;
   $organism_tax_id=~s/^\s+//;$organism_tax_id=~s/\s+$//;
   next if($org_id==58);
   $orgmap{$organism_tax_id}="$org_id,$organism_sc_name";
}
my %files_hd=();
#collect a
$zip_link="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz";
$zip_name="uniprot_sprot.dat.gz";$text_file="uniprot_sprot.dat"; #$test_file="test_uniprot_sprot.dat";
if(-f $zip_name){system("mv $zip_name archive");}
if(-f $text_file){system("mv $text_file archive");}
if( !(-f $text_file)&&!(-f $zip_name)){system("wget -q $zip_link");}
if((-f $zip_name)&& !(-f $text_file)){system("gunzip -q $zip_name");}
print "$text_file generated \n";
$more=1; my @domain=();$rows=""; 
$sql_file="$opt_o/uniprot_sprot_load.txt"; $index=0; 
$sequence_file="$opt_o/uniprot_sprot_sequence_load.txt";
if((-f $text_file)){ #get the specified fields
 #  system("head -41000 $text_file > $test_file");
   open(FH,"$text_file"); $lcount=0; open(OUT,">$sql_file");
   open(SEQ,">$sequence_file");
   while(<FH>){ chomp($_); 
     if(!($_=~/^\/\//)){push(@domain,"$_");}
     else{     
          ##get supporting fields for this protein  
          @taxa=grep(/OX\s/,@domain);if(@taxa<=0){@domain=();next;}
          $refseq=""; $ensembl="";$NCBI_TaxID=0;$GeneName="";$isoform="";
          $seq_len=0;$seq="";@interprot=();@ensembl=();@refSeq=();
          $uniprotProteinid="";$ProteinName="";$ProteinEvidence="===";@feature=();
          $pfams="";$genename="";#CC         IsoId=P15379-14; Sequence=Displayed;
          @isof=grep(/Sequence=Displayed/,@domain);
          if(@isof>0){
            if($isof[0]=~/IsoId=(.+);\s*Sequence=Displayed/){$isoform=$1;}
          }
          while(@taxa){$line=shift(@taxa);chomp($line);
             if($line=~/NCBI_TaxID=(\d+);/){$NCBI_TaxID=$1;last;}
          }
          #check if we have this organism in our database
          if(!exists($orgmap{$NCBI_TaxID})){@domain=();next;}
          #DR   Ensembl; ENSMUST00000050487; ENSMUSP00000057836; ENSMUSG00000024610.
          #DR   Ensembl; ENSMUST00000097563; ENSMUSP00000095171; ENSMUSG00000024610.
          #DR   Ensembl; ENSMUST00000167610; ENSMUSP00000126688; ENSMUSG00000024610.
          @ensembl= grep(/Ensembl/,@domain);$is_gene=0;
          @interprot=grep(/InterPro;/,@domain);
          @refSeq=grep(/RefSeq/,@domain);
          #if(@ensembl<=0){ #DR   GermOnline; ENSMUSG00000005087; Mus musculus.
          #   if(
          #  @ensembl= grep(/ENSMUSG/,@domain);
          # }
          ($org_id,$organism_sc_name)=split(",",$orgmap{$NCBI_TaxID});
          while(@domain>0){
             $line=shift(@domain);chomp($line);$line=~s/^\s*//;$line=~s/\s*$//;
             @tokens=split(/\s/,$line);
             if(@tokens>0){ $line1=$tokens[0];$line1=~s/^\s*//;$line=~s/\s*$//;
                if($line1=~/^AC$/){$i=1; 
                   while($i<@tokens){$uniprotProteinid.=$tokens[$i];++$i;}
                 }
                if($line1=~/^DE$/){
                     if($line=~/^DE\s+RecName: Full=(.+);/){$ProteinName=$1;}
                 }
                if($line1=~/^PE$/){
                    if($line=~/^PE\s+.+:\s*(.+);/){$ProteinEvidence=$1;}
                }
                if($line1=~/^FT$/){push(@feature,$line);} 
                if($line1=~/^DR$/){
                   if($line=~/Pfam;\s+(.+)\s+\d+\./){$pfams.=":$1";}
                 }
                 if($line1=~/^GN$/){
                     if($line=~/Name(s)?=(.+);/){$genename.="$2;";}
                 }
              }
             #get the sequence data SQ   SEQUENCE   778 AA;
             if($line=~/SQ\s+SEQUENCE\s+(\d+)\s*AA;/){
                $seq_len=$1;$seq="";
                while(@domain>0){
                   $line=shift(@domain);chomp($line);$line=~s/\s*//g;
                   $seq.=$line;
                }
             }
          }
          #YP_001816784.1; NC_010563
          $data= "$org_id\t$uniprotProteinid\t$ProteinName\t$ProteinEvidence\t$pfams\t$genename";
          $computed_len=length($seq);$i=0;$rows="";
          print SEQ "$org_id\t$uniprotProteinid\t$isoform\t$seq_len\t$seq\n";
          if(@ensembl>0){
             foreach $row(@ensembl){$row=~s/^\s*//;
                next if(!($row=~/^DR\s+/));
                ($token,$tx_id,$prot_id,$gene)=split(";",$row);$gene=~s/\.$//;
                if(@feature>0){
                   $i=0;
                   while($i<@feature){
                         $line=$feature[$i]; chomp($line);$line=~s/^\s*//;$line=~s/\s*$//;++$i;
                         if($line=~/^FT\s+(.+)\s+(\d+)\s+(\d+)\s+(.+)$/){
                            @fields=split(/\s+/,$line);
                            if(@fields>=4){
                               $type=$fields[1]; $start= $fields[2];$end=$fields[3];$name=$4;
                               $rows.="$data\t$tx_id\t$prot_id\t$gene\t$type\t$start\t$end\t$name\n";
                            }
                         }
                     }
                 }
              }    
          }
          if(@refSeq>0){$i=0;
            foreach $row(@refSeq){$row=~s/^\s*//;
                next if(!($row=~/^DR\s+/));
                ($token,$prot_id,$tx_id)=split(";",$row);$gene=$tx_id;
                if(@feature>0){$i=0;
                   while($i<@feature){
                         $line=$feature[$i]; chomp($line);$line=~s/^\s*//;$line=~s/\s*$//;++$i;
                         if($line=~/^FT\s+(.+)\s+(\d+)\s+(\d+)\s+(.+)$/){
                            @fields=split(/\s+/,$line);
                            if(@fields>=4){
                               $type=$fields[1]; $start= $fields[2];$end=$fields[3];$name=$4;
                               $rows.="$data\t$tx_id\t$prot_id\t$gene\t$type\t$start\t$end\t$name\n";
                            }
                         }
                     }
                 }
              }    
          }
          if(@interprot>0){$i=0;
            foreach $row(@interprot){next if(!($row=~/^DR\s+/));
                ($token,$prot_id,$tx_id)=split(";",$row);$gene=$tx_id;
                if(@feature>0){$i=0;
                   while($i<@feature){
                         $line=$feature[$i]; chomp($line);$line=~s/^\s*//;$line=~s/\s*$//;++$i;
                         if($line=~/^FT\s+(.+)\s+(\d+)\s+(\d+)\s+(.+)$/){
                            @fields=split(/\s+/,$line);
                            if(@fields>=4){
                               $type=$fields[1]; $start= $fields[2];$end=$fields[3];$name=$4;
                               $rows.="$data\t$tx_id\t$prot_id\t$gene\t$type\t$start\t$end\t$name\n";
                            }
                         }
                     }
                 }
              }    
          }
          print OUT "$rows";
        }
  }
}
#close(OUT);close(SEQ);
print "domain file generated --\n";
if(-f "$sql_file"){
  $qh_droptable->execute();$qh_create_table->execute();
  open(IN,"$sql_file");$rows="";$index=0;
  while(<IN>){ chomp($_);$rows.="$_\n"; ++$index;
     if($index%50000==0){
         $filename="sql_load_segment.txt";open(OUT,">$filename");
         if(OUT){
            print OUT "$rows";$rows="";
         }
         close(OUT);
         if(-f "$filename"){
              $filename="$opt_o/$filename"; #get the number of lines in the file
              $linecont = `wc -l $filename`;chomp($linecont);
              if($linecont=~/^(\d+)\s$filename/){$linecont=$1;
                if($linecont>0){ $qh_load_this->execute($filename)or die "Couldn't prepare statement: " . $dbh->errstr;
                   print "Loaded new segment\n";
                 }
               } 
          } 
      }   
   }
   if($rows ne ""){
       $filename="sql_load_segment.txt";open(OUT,">$filename");
       if(OUT){print OUT "$rows";$rows="";}close(OUT);
         if(-f "$filename"){$filename="$opt_o/$filename";
              $linecont = `wc -l $filename`;chomp($linecont);
              if($linecont=~/^(\d+)\s$filename/){$linecont=$1;
                 if($linecont>0){ 
                    $qh_load_this->execute("$filename")or die "Couldn't prepare statement: ".$dbh->errstr;
                    print "Loaded last segment\n";
                 }
             } 
       } 
    }
   #get the total number of lines of this segment that were inserted
   $qh_get_rowcount->execute();
  ($rowcount)= $qh_get_rowcount->fetchrow_array();
   #now insert this info into the log
   #$filename="$sql_file";
   $linecont = `wc -l $sql_file`;chomp($linecont);
   if($linecont=~/^(\d+)\s$sql_file/){$linecont=$1;
      $qh_delete_log->execute();
      $qh_insert_log->execute("$sql_file",$linecont,$rowcount);
      print "Log table updated\n";
    }
   close(IN);
 }
$qh_analyze->execute();
print "Loding sequence table\n";
if(-f $sequence_file){
  $qh_dropSeq->execute();$qh_create_Seq->execute();
  open(IN,"$sequence_file");$rows="";$index=0;
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
                    print "Loding last segment of sequence table\n";
                 }
             } 
       } 
    }
   #get the total number of lines of this segment that were inserted
   $qh_get_seq_rowcount->execute();
  ($rowcount)= $qh_get_seq_rowcount->fetchrow_array();
   #now insert this info into the log
   $filename="$opt_o/sequence_load.txt";
   $linecont = `wc -l $sequence_file`;chomp($linecont);
   if($linecont=~/^(\d+)\s$sequence_file/){$linecont=$1;
      $qh_delete_seq_log->execute();
      $qh_insert_seq_log->execute("$sequence_file",$linecont,$rowcount);
    }
   close(IN);
 }
$qh_analyze_seq->execute();

print "Porgram complete\n";


