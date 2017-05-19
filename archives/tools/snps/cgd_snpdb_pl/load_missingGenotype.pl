#!/usr/bin/perl
use vars qw($opt_h $opt_f $opt_c);
use Getopt::Std;use DBI;
getopts('hf:c:');

if($opt_h||!$opt_f||!$opt_c) { #||!$opt_c
    print <<HELP;

 This script loads missing genotype data

Usage:

   perl load_missingGenotype.pl -f <config_file> -c <chrom_list>

Arguments:
   -h  displays this help message
   -f  SNP config file
   -c  commas separated list of chromosomes

Example: perl load_missingGenotype.pl -f  docs/Yang2011_snp_source.config -c X,Y
HELP
exit;
}
my $user ='lnh';
my $pwd  ='sql4lucie';
$dbname ='cgdsnpdb';
$host ='cgd.jax.org';
$dbname="cgdsnpdb";
$dbh=DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user,$pwd);
if(!$dbh){print mysql_error(); exit;}
#I need to load the strain list and index strain ids in the header
$qh_getStrainID=$dbh->prepare("select strain_id from snp_strain_synonym where synonym_name=?");
open(IN,"$opt_f");
@config=<IN>;@variable=grep(/=/,@config);
my ($snp_dir,$source,$output_dir,$source_id);my $source_id=0;
$all_in_one=0;$local_idIndex=-1;$alleles_index=-1;
$snp_file_name="";$geno_sufix="";$snp_main_sufix="";
$snp_allele_index=-1;$ref_allele_index=-1;$ref_strainID=0;
$first_strainIndex=-1;$chrom_index=-1; $pos_index=-1;
foreach my $line(@variable){# set global variables from the config file
  chomp($line);
  if($line=~/SNP_SOURCE=(.+)$/){$source=$1;}
  if($line=~/SOURCE_ID=(.+)$/){$source_id=$1;}
  if($line=~/PIPELINE_DIR=(.+)$/){$output_dir=$1;}
  if($line=~/SNP_MAIN_FIRST_STRAIN_INDEX=(.+)$/){$first_strainIndex=$1;}
  if($line=~/SNPID_LOCAL_INDEX=(.+)$/){$local_idIndex=$1;}
  if($line=~/SNP_ALLELES_INDEX=(.+)$/){$alleles_index=$1;}
  if($line=~/SNP_NAIN_CHROM_INDEX=(.+)$/){$chrom_index=$1;}
  if($line=~/REF_STRAIN_ID=(.+)$/){$ref_strainID=$1;}
  if($line=~/SNP_NAIN_BP_POSITION_INDEX=(.+)$/){$pos_index=$1;}
  if($line=~/SNP_REF_ALLELE=(.+)$/){$ref_allele_index=$1;}
  if($line=~/OTHER_ALLELE=(.+)$/){$snp_allele_index=$1;}
  if($line=~/ALL_IN_ONE=(.+)$/){$all_in_one=$1;}
  if($line=~/SNP_MAIN_SUFIX=(.+)$/){$snp_file_name=$1;}
  if($line=~/GENO=(.+)$/){$geno_sufix=$1;}
}
close(IN);
%complement=("A"=>"T","T"=>"A","C"=>"G","G"=>"C");@chroms=split(",",$opt_c);
%strainsIndex=();%fieldsIndex=();
while(@chroms>0){$chr=shift(@chroms);
  if(-f "$output_dir/chr$chr$snp_file_name"){
     open(IN,"$output_dir/chr$chr$snp_file_name");$line=<IN>; chomp($line);@header= split(",",$line);
     $geno_file="$chr$geno_sufix";
     open(OUT,">$output_dir/$geno_file")or die "Bad fiel: $!\n";
     print OUT "snpID\tsource_id\tstraind_id\tgenotype_allele\tref_strain_id\n";
     for my $i(0 .. $#header){ $strain_id=0;$field=$header[$i];
       if($i>=8){
          $strain_name=$header[$i];$strain_name=~s/\s*$//;$qh_getStrainID->execute($strain_name);
          if($qh_getStrainID->rows>0){($strain_id)=$qh_getStrainID->fetchrow_array();}
          $strainsIndex{"$i"}=$strain_id;
       }
       $fieldsIndex{"$field"}=$i;#print "$i\t$strain_id\t--".$header[$i]."--\n";
     }
    while(<IN>){chomp($_); @fields=split ",",$_;
      $local_id=$fields[$local_idIndex];$alleles=$fields[$alleles_index];
      $local_ref_allele="";$local_snp_allele=""; $snp_allele=uc($fields[$snp_allele_index]);
      if($alleles=~/^(\w{1})\/(\w{1})/){
         $local_ref_allele=uc($1);$local_snp_allele=uc($2);}
      $ref_allele=uc($fields[$ref_allele_index]);next if($local_ref_allele eq "");
      $ref_allele=~s/^\s*//;$ref_allele=~s/\s*$//;$snp_allele=~s/^\s*//;$snp_allele=~s/\s*$//;
      $reverse=0;
      next if($local_id<=0);
      if($ref_allele ne $local_ref_allele){if($complement{"$ref_allele"} eq $local_ref_allele){$reverse=1;}
      }
      for my $i(8 .. $#fields){#snpID	source_id	straind_id	genotype_allele	ref_strain_id
           $genotype_allele=uc($fields[$i]);$genotype_allele=$1 if($genotype_allele=~/^(\w{1})/);
           $genotype_allele=($reverse==1)?$complement{"$genotype_allele"}:$genotype_allele;
           print OUT "$local_id\t$source_id\t".$strainsIndex{"$i"}."\t$genotype_allele\t$ref_strainID\n";
       }
    }
   close(OUT);
  }
}
exit(0);
