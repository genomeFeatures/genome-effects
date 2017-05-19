#!/usr/bin/perl

######################################################################
## This script gets SNP annotations from the current version
## of ensembl database . The organism is mus_musculus
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : September 2011
#   Usage: ./getEnsemblvariations.pl 
#
#  Output:
#
# Assumption: the asumption is that ensembl database namimg convention is
#              mus_musculus_variation_&_&c
#  Where & is one or more digits and c is a character
# for example : mus_musculus_variation_62_37o
#
# Tables to load:
q{
  cgd_ensembl_variations
};
#######################################################################
## set default variables
$data_dir=`pwd`;                        #default data directory
$local_db="cgd_snpdb37";                #default lookingglass snp database name
$localhost="lookingglass.jax.org";      #deafult db server
$luser="lnh";$lpwd="sql4lucie";        
$sql_scripts_path="/home/lnh/ensembl_load/cgd_snpdb_sql";


use DBI;
use vars qw ( $opt_h $opt_o);
use Getopt::Std;
use Time::localtime;
getopts('hi:o:');
if($opt_h) {
    print <<HELP;

 This script gets SNP annotations from the current version
 of ensembl database . The organism is mus_musculus

Usage:

   perl getEnsemblvariations.pl -d <path to data directory>

Arguments:
   -h  displays this help message
   -o  path to data directory (optional)

HELP
exit;
}
####################################### Queries ##########################################
 
#get ensembl variations
$get_variations="select (case when s.name='X' then 20 when s.name='Y' then 21 when s.name='MT' then 22 else s.name end)as 
           chromosome_id,seq_region_start as bp_position,(case when seq_region_strand=1 then '+' else '-' end)as strand,
           variation_name as accession_id,tv. allele_string as alleles,tv.consequence_types as function_flass,tv.cds_start posInCds, 
          tv.translation_start as posInProtein,tv.codon_allele_string as codons,tv.pep_allele_string as aminoacids,
          feature_stable_id as tx_accession_id from transcript_variation tv, variation_feature v, seq_region s  
         where  v.seq_region_id=? and v.seq_region_id=s.seq_region_id and v.variation_feature_id=tv.variation_feature_id
          and tv.allele_string=v.allele_string";
#cds_start >0 and feature_stable_id="ENSMUST00000169862" and


$seq_regionQuery="select s.name,s.seq_region_id from seq_region s, coord_system c where s.coord_system_id=c.coord_system_id 
                    and c.rank=1";


########################################
if($opt_o){$data_dir=$opt_o;} 

# set up some flags
$bad_localdb_con=0;$bad_remotedb_con=0;
#local db connection
my $dbh2 = DBI->connect("DBI:mysql:$local_db:$localhost",$luser,$lpwd);
if(!$dbh2){
    print  "Could not connect to the local database $local_db:$!\n"; exit;
}

$host="ensembldb.ensembl.org";
$port="5306";
$user="anonymous";

$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print "Starting SNP download :  $mday/$mon/$yday @ $hour:$min:$sec \n";

$qh_getEnsemblSchema=$dbh2->prepare("select schema_version from cgd_externaldb where db_name='ENSEMBL'");
$qh_getEnsemblSchema->execute();$db_version="";
if(($db_version)=$qh_getEnsemblSchema->fetchrow_array()){
   # check if data directory is provided by the user
  # print "The current version is $db_version\n";
   $variation_db=$db_version; $variation_db=~s/core/variation/i;
   #print "The SNP database is $variation_db\n"; 
   #ensembl db connection
   my $dbh = DBI->connect("DBI:mysql:$variation_db:$host:$port",$user);
   if(!$dbh){
       print  "Could not connect to Ensembl database :$!\n"; exit;
   } 
   $qh_getvariations=$dbh->prepare($get_variations);
   $data_dir=~s/\s+$//;$data_dir=~s/\/$//;
   $organism_data_dir="$data_dir/$db_version"; 
   $qh_getVariationCount= $dbh2->prepare("select count(*) from cgd_ensembl_variations where chromosome_id=?");
   $qh_updVariation= $dbh2->prepare("UPDATE cgd_annot_load_test set db_table_line_count=? where db_table_name=?");

   if(!(-d "$organism_data_dir")){ mkdir("$organism_data_dir",0777);}
   if(-d "$organism_data_dir"){
         #generate seq_region_file
         $seq_file="$organism_data_dir/ensemblSeqregionList.txt";
         $db_comand="mysql --user=$user --port=$port -h $host -e";
         $command="$db_comand 'use $db_version;$seq_regionQuery'> $seq_file";
         if(!(-f $seq_file)){
             print "we will generate the chromosome file\n";
             system($command);
         }
         #print "okido \n";
       
        #now open the process the chromosome file and download SNPs for each chromosome
        open(CHR,"$seq_file");
       if(CHR){
           $header=<CHR>; $count=0;
           open(OUT, ">$sql_scripts_path/ensemblVariations_load.sql");
           if(OUT){
               print OUT "USE $local_db;\n\n";
                    ##now insert the schema
               open(SCHEMA,"$sql_scripts_path/ensemblVariations_schema.sql");
               while($line=<SCHEMA>){
                      chomp($line);print OUT "$line\n";
               }
               close(SCHEMA);
           }
           close(OUT);
           $comand=" mysql -h$localhost -u$luser -p$lpwd < $sql_scripts_path/ensemblVariations_load.sql";
           system($comand);
           while($line=<CHR>){
               chomp($line);
               my($chr,$seq_region_id)=split("\t",$line);
               #get the snps for this chromosome
               #last if($count>0); ++$count ;
               if($seq_region_id=~/\d+/){
                 $chr_id=$chr; $chr_id=20 if($chr=~/X/i);$chr_id=21 if($chr=~/y/i);
                 $chr_id=22 if($chr=~/m/i);$chr_file="$organism_data_dir/chr$chr"."_SNPs.txt";
                 open(OUT, ">$sql_scripts_path/ensemblVariations_load.sql");
                 if(OUT){
                    print OUT "USE $local_db;\n\n";
                 }
                 if(!(-f $chr_file)){
                      $qh_getvariations->execute($seq_region_id);$item_count=$qh_getvariations->rows;
                      if($item_count>0){
                          open(CF,">$chr_file");$item_count=0;            
                          while(($chromosome_id,$bp_position,$strand,$accession_id,$alleles,
                                  $function_flass,$posInCds,$posInProtein,$codons,$aminoacids,$tx_accession_id)=
                                  $qh_getvariations->fetchrow_array()){
                                  if(CF){
                                     $posInCds=0 if(!$posInCds);  $posInProtein=0 if(!$posInProtein);
                                     $codons="---" if($codons eq "");$aminoacids="---" if($aminoacids eq "");
                                     $frame=0;
                                     if($codons=~/^\s*(\w{3})\/(\w{3})\s*$/){
                                        $ref_codon=$1; $c_codon=$2;
                                        if(substr($ref_codon,0,1) ne substr($c_codon,0,1)){
                                           $frame=1;
                                        }
                                        elsif(substr($ref_codon,1,1) ne substr($c_codon,1,1)){
                                           $frame=2;
                                        }
                                        else{$frame=3;}
                                     }
                                     $data="$chromosome_id\t$bp_position\t$strand\t$accession_id\t";
                                     $data.="$alleles\t$function_flass\t$posInCds\t$posInProtein\t";
                                     $data.="$codons\t$aminoacids\t$tx_accession_id\t$frame\n";
                                     if($alleles=~/^\s*\w{1}\/\w{1}\s*$/){
                                        print CF "$data"; ++$item_count;
                                     }
                                   }
                            } #end of while SNPs
                            close(CF);
                         } #end of if $item_count>0
                 }
                 else{
                     #get line count of this file
                      $command="wc -l $chr_file";
                      $item_count=`$command`;
                      if($item_count=~/(\d+)/){
                        $item_count=$1;
                      }
                  }
                  $table_name="cgd_ensembl_variations $chr";
                  if(OUT){
                      if((-f "$chr_file")){
                          print OUT "LOAD DATA LOCAL INFILE '$chr_file' into table cgd_ensembl_variations;\n";
                          print OUT "DELETE FROM cgd_annot_load_test WHERE db_table_name='$table_name';\n\n";
                          print OUT "INSERT INTO cgd_annot_load_test VALUES('$chr_file',$item_count,'$table_name',0);\n\n";       
                      }
                     close(OUT);
                   }
                  #now load data into the database
                  $comand=" mysql -h$localhost -u$luser -p$lpwd < $sql_scripts_path/ensemblVariations_load.sql";
                 q{ system($comand);
                 
                  $qh_getVariationCount->execute($chr_id);
                  if(($rowCount)=$qh_getVariationCount->fetchrow_array()){
                       $qh_updVariation->execute($rowCount,"$table_name");
                  }
                  };
                } #end of seq_region id validity
                 print "Chromosome $chr SNPs loaded\n"; #now validate the updates (table loads)
             } #end of while chromosome list
        } #end of chromosome file valid
  print "SNPs annotations downloaded\n";
  #################################### SNP annotator tables ###############################
 }
}
$command="mysql -h$localhost -u$luser -p$lpwd -e 'use $local_db;";
$command.=" ANALYZE TABLE cgd_ensembl_variations;'";
system($command);
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print "Program complete at :  $mday/$mon/$yday @ $hour:$min:$sec \n";
exit(0);
           
 

