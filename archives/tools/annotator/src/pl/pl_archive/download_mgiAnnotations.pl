#!/usr/bin/perl
########################################################################
## This script collects all the gene annotations
## from ftp://ftp.informatics.jax.org/pub/mgigff/MGI.gff3.gz
#  and load them into graber_transcriptdb.mgi table
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : February 2012
#   Modified: March 2013
#  Input : a directory path where to process files
#  Output: download, unzip and parse the file into
#    a tab-delimited file with all the annotations (gene, transcripts,exon,cds)
#    and their cordinates and relationship. mgiGene coordinates are in 1-base 
#
# Problem: The CDS segment is associated to the transcript parent but not the exon.
#    there is no direct association between the exon and the CDS feature
# making it harder to associate a given CDS segment to the corresponding exon.
# This is not an issue for single exon-transcript but if a transcript has multiple
# exons, to establish the association exon-cds_segment, you have to
# first collect all the exons associated to a given transcript then
# for every cds, check the exon that coordinates overlap this cds coordinates
#
# also you have to compute the cds start and stop 
# Note MGI updates this file every 4th of the month
# So I will set a cron job that runs every 7th of the month
#
# Data sources: ftp://ftp.informatics.jax.org/pub/mgigff/MGI.gff3.gz
#**********************************************************************
## set default variables
$dbname="graber_transcriptdb";$host="harlequin"; $user="lnh";$pass="lucie98";
$oversion="mm9"; $gene_prediction="mgiGene";
$load_Annotations="~/work/projects/graber_transcriptdb/src/pl/load_annotations.pl";
use DBI;use vars qw ($opt_h $opt_f $opt_o);
use Getopt::Std;use LWP 5.64;
#use Time::localtime;
use Time::localtime;
$organism=""; # default organism
$inputfile="";$names="";$outputfile="";@geneList=();
my $browser = LWP::UserAgent->new;
   $browser->timeout(10);$browser->env_proxy;
   $browser->agent('Mozilla/5.0');

getopts('ho:f:');
if(($opt_h)||(!$opt_o)){
    print <<HELP;

 This script collects all the gene annotations
 from ftp://ftp.informatics.jax.org/pub/mgigff/MGI.gff3.gz
  and load them into graber_transcriptdb.mgi table

Usage:
   perl download_proteinDomains.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  Output directory
  

Examples:
cmd: ./download_mgiAnnotations.pl -o /scratch/data/downloads/mgi

HELP
exit;
}
chdir($opt_o);
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

## get the date of the last update
my $get_lastUpdate = "select schema_version from external_db where db_name=?";
my $qh_lastUpdate  = $dbh->prepare($get_lastUpdate)or die "Couldn't prepare statement: " . $dbh->errstr;

my $drop_mgi_temp = "drop temporary table if exists mgi_genes_temp";
my $qh_drop_mgi_temp = $dbh->prepare($drop_mgi_temp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $create_mgi_temp = "create temporary table mgi_genes_temp(mgi_id varchar(25),
          mgi_symbol varchar(100),biotype varchar(100),found tinyint default 0,
         index(mgi_id),index(mgi_symbol),index(biotype))";
my $qh_create_mgi_temp = $dbh->prepare($create_mgi_temp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $load_this = "load data local infile ? into table mgi_genes_temp";
my $qh_load_this = $dbh->prepare($load_this)or die "Couldn't prepare statement: " . $dbh->errstr;

my $update_this= "update mgi_genes_temp t, mgi_genes m set found=1 where t.mgi_id=m.mgi_id";
my $qh_update_this = $dbh->prepare($update_this)or die "Couldn't prepare statement: " . $dbh->errstr;

my $insert_mgigene = "insert into mgi_genes(mgi_id,mgi_symbol,biotype) ";
   $insert_mgigene.="select distinct mgi_id,mgi_symbol,biotype from mgi_genes_temp where found=0 ";
my $qh_insert_mgigene = $dbh->prepare($insert_mgigene)or die "Couldn't prepare statement: " . $dbh->errstr;

my $get_mgi_genes_rowcount="select count(*) as rowcount from mgi_genes_temp ";
my $qh_get_mgi_genes_rowcount = $dbh->prepare($get_mgi_genes_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;

## update the external_db
my $Update_schema="update external_db set schema_version= ? where db_name=?";
my $qh_Update_schema= $dbh->prepare($Update_schema);

my $qh_insertNewdb=$dbh->prepare("Insert into gene_prediction(gene_prediction_name) values(?)");
my $qh_getNewdbId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name=?");
my $qh_insertExternaldb=$dbh->prepare("insert into external_db values(?,?,?)");
my $qh_deleteSchema=$dbh->prepare("delete from external_db where db_name=?");


$zip_link="ftp://ftp.informatics.jax.org/pub/mgigff/MGI.gff3.gz";$zip_name="MGI.gff3.gz";$text_file="MGI.gff3";
if(-f $zip_name){system("mv $zip_name archive");}
if(-f $text_file){system("mv $text_file archive");}
if(!(-f $text_file)&&!(-f $zip_name)){system("wget -q $zip_link");}
if((-f $zip_name)&& !(-f $text_file)){system("gunzip -q $zip_name");}
$more=1; my @annots=();$rows="";
open(LOG,">mgi_log.txt");$tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Starting load process :  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
 #check when this file was last updated and compare with what we have in the database
 @heads=`head -30 $text_file`;$newupdate="";$current_version="";
 $head=join("", @heads);
  if($head=~/MGI database:\s*#\s+Last updated:\s*(.+,\s?\d+)\s*#/){$newupdate=$1;}
  if($head=~/Genome build:\s+(.+)-C57BL/){$version=$1;
       if($version=~/GRCm38/){$oversion="mm10"; $gene_prediction="GRCm38-mgiGene";}
  }
######################################################
$qh_lastUpdate->execute($gene_prediction);
if($qh_lastUpdate->rows>0){($current_version)=$qh_lastUpdate->fetchrow_array();}
if($newupdate eq $current_version){ print "No new updates ($newupdate eq $current_version) \n";}
else
 { 
   $qh_getNewdbId->execute($gene_prediction); #get the id of the annotation source if exists
   if($qh_getNewdbId->rows<=0){               #new version, insert
      $qh_insertNewdb->execute($gene_prediction);$qh_getNewdbId->execute($gene_prediction);
      ($gene_prediction_id)=$qh_getNewdbId->fetchrow_array();
      $qh_deleteSchema->execute($gene_prediction); 
      $qh_insertExternaldb->execute($gene_prediction_id,$gene_prediction,$newupdate);
   }
  else{($gene_prediction_id)=$qh_getNewdbId->fetchrow_array();
     if($current_version eq ""){$qh_insertExternaldb->execute($gene_prediction_id,$gene_prediction,$newupdate);}
     elsif($current_version ne $newupdate){$qh_Update_schema->execute($newupdate,$gene_prediction);}
  }
  $main_file="$opt_o/$oversion-$gene_prediction.txt"; $index=0;  %gene_list=();
  if(-f $text_file){ 
         #get the list of transcript, gene, exons ids of this organism version block by block
         open(FH,"$text_file"); $lcount=0; open(OUT,">$main_file");
         open(GE,">mgi_genes.txt"); #data will be loaded into mgi_genes table
         if(OUT){
            print OUT "name\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\t";
            print OUT "exonFrames\tname2\tmgi_id\ttranscriptType\tbiotype\n";  
            while(<FH>){ chomp($_);$_=~s/^\s+//; next if($_=~/^#\s|^##gff-/); 
              if(!($_=~/^#{3}/)){ push(@annots,"$_"); }
              else{ ##get data for this block
                %transcripts=();@genes=();@genes=grep(/gene/,@annots);
                while(@genes>0){ $line=shift(@genes);
                       ($chr,$db,$type,$start,$end,$score,$strand,$phase,$attr)=split("\t",$line);
                       if($attr=~/ID=MGI:MGI:(\d+);Name=(.+);Dbxref.+;bioType=(.+)/){
                         $gene_list{"MGI:$1"}{"name"}=$2; $gene_list{"MGI:$1"}{"biotype"}=$3;}
                }
                next if(keys(%gene_list)<=0); #skip cases we do not have gene annotation
                @transcript=grep(/transcript|mRNA/,@annots);
                $name="";$chrom="";$strand="";$txStart=0;$txEnd=0;$cdsStart=0;$cdsEnd=0;
                $exonCount=0;$exonStarts="";$exonEnds="";$score=0;$name2="";$exonFrames=""; 
                next if(@transcript<=0);  #skip cases we do not have transcript annotation
                while(@transcript>0){
                   $line=shift(@transcript);
                   ($chr,$db,$ttype,$tstart,$tend,$score,$tstrand,$phase,$attr)=split("\t",$line);
                     if($attr=~/ID=$db:(.+);Parent=.+;Dbxref=MGI:MGI:(\d+)/){
                         $cds_start=0;$cds_end=0;$exon_count=0;
                         $tx_id=$1;$mgi_id="MGI:$2";$exonStarts="";$exonEnds="";$score;$exonFrames="";
                         $name2=$gene_list{"$mgi_id"}{"name"}; $biotype=$gene_list{"$mgi_id"}{"biotype"};
                         $transcripts{$tx_id}{"cds"}=();@tx_data=grep(/$tx_id/,@annots);
                         @exons=grep(/exon/,@tx_data); @cds=grep(/CDS/,@tx_data);
                         while($line=shift(@cds)){ #load all the cds of this gene
                            ($chr,$db,$type,$start,$end,$score,$strand,$phase,$attr)=split("\t",$line);
                            if($attr=~/Parent=$db:(.+);/){$transcripts{$1}{"cds"}{$start}="$end,$phase";}
                         }
                         while($line=shift(@exons)){
                               ($chr,$db,$type,$estart,$eend,$score,$strand,$phase,$attr)=split("\t",$line);
                             if($attr=~/Parent=$db:(.+);/){
                                if($1 eq "$tx_id"){ #only get the exons associate with this transcript
                                     $exonStarts.="$estart,"; $exonEnds.="$eend,"; $frame="-1";++$exon_count;
                                     if(keys(%{$transcripts{$tx_id}{"cds"}})>0){
                                        for my $cstart(sort keys(%{$transcripts{$tx_id}{"cds"}})){
                                            ($cend,$cphase)=split(",",$transcripts{$tx_id}{"cds"}{$cstart});
                                            if($cds_start==0){$cds_start=$cstart;}
                                            else{if($cds_start>$cstart){$cds_start=$cstart;}}
                                            if($cds_end==0){$cds_end=$cend;}
                                            else{if($cds_end<$cend){$cds_end=$cend;}}
                                            if(($cstart>=$estart)&&($cend<=$eend)){$frame=$cphase;}
                                        } #end of cds loop
                                     } #end of cds check
                                     $exonFrames.="$frame,";
                                  }#end of  if($1 eq "$tx_id"){ 
                              } #end of if($attr=~/Parent=$db:(.+);/){
                        } #end of while($line=shift(@exons)){
                        #display this transcript
                      $chr=~s/^\s*//;$chr=~s/\s*$//;
                      print OUT "$tx_id\t$chr\t$strand\t$tstart\t$tend\t$cds_start\t$cds_end\t";
                      print OUT "$exon_count\t$exonStarts\t$exonEnds\t$exonFrames\t$name2\t$mgi_id\t$ttype\t$biotype\n";
                    } #end of  if($attr=~/ID=$db
                } #end of while(shift(@transcript)){
              @annots=();
             } #end of this gene
           } #end of file loop  while(<FH>){
        } #end of if if(OUT){
       close(OUT); 
      #now generate gene list
      for my $mgi_id(sort keys(%gene_list)){
         $gene=$gene_list{"$mgi_id"}{"name"}; $biotype=$gene_list{"$mgi_id"}{"biotype"};
         print GE "$mgi_id\t$gene\t$biotype\t0\n" if(GE);
      }
      close(GE) if(GE);
      $filename="$opt_o/mgi_genes.txt";
      #get the number of lines in the file
      $linecont = `wc -l $filename`;chomp($linecont);
      if($linecont=~/^(\d+)\s$filename/){$linecont=$1;
         if($linecont>0){
            $qh_drop_mgi_temp->execute();$qh_create_mgi_temp->execute();
            $qh_load_this->execute("$filename") or die "bad query ".mysql_error()."\n";
            $qh_get_mgi_genes_rowcount->execute();($rowcount)=$qh_get_mgi_genes_rowcount->fetchrow_array();
            if($linecont==$rowcount){$qh_update_this->execute();$qh_insert_mgigene->execute();
               print LOG "$filename Loaded  $rowcount of $linecont lines were loaded\n";}
            else{ 
               print LOG "$filename Load failed:only $rowcount of $linecont lines were loaded\n";
               #the bulk load did not go well
             }
          }
       } 
       #now call the database update script
       $cmd=" perl $load_Annotations -d /scratch/data/downloads/mgi -v $oversion -f $main_file -a $gene_prediction -s 1 -e 1";
       system($cmd);  
     } #end of if(-f $test_file){
 } #end of if($newupdate eq $dbupdate){
 $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
  print LOG "\n*************************************************************\n";
  print LOG "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print LOG "\n*************************************************************\n";
 close(LOG);
print "Porgram complete\n";
exit(0);

