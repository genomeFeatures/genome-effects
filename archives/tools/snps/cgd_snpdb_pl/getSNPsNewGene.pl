#!/usr/bin/perl
################################################################
q{
This script generates SNP annotations for selected strains

Author: Lucie Hutchins, Scientific Software Engineer
Date : October 2011

};
################################################################

use DBI;
use vars qw ($opt_h $opt_s $opt_c);
use Getopt::Std;
use POSIX;
getopts('hs:c:');
if($opt_h||!$opt_s ||!$opt_c) {
    print <<HELP;

 This script generates SNPs annotation for provided strains

Usage:

   perl getSNPsNews.pl -s <strain_file> -c <snp_file>]

Arguments:
   -h  displays this help message
   -s  The file containing name(s) of selected strains (multiple names are separated with a commas
   -g  A flag to specified to only include SNPs on transcript
   -d  A flag to filter rows based on genotype similarity :
        0-> display all cases ; 
        1-> only display cases where the genotype is different between selected strains
        2-> only display cases where the genotype is the same between selected strains
   -p  name of the snp provider (source)
   -c  file containing snp positions format: chrom:position

Example: perl getSNPsNew.pl -s strain_filename -c snp_file_list

HELP
exit;
}


my $dbname ='cgd_snpdb';
my $host ='cgd-dev';
my $user ='pup';
my $pwd  ='puppass';

my $dbh = DBI->connect("DBI:mysql:$dbname:$host",$user,$pwd)
                            or die "Can't connect to database!\n";
$source_id=21;$imputed_genotype_table="snp_imputed"; $genotype_table="snp_genotype";
# get chromosomes list
my $qh_getChr =$dbh->prepare("select chromosome_id from snp_chromosome where chromosome_name=?");
my $qh_getChr_name =$dbh->prepare("select chromosome_name from snp_chromosome where chromosome_id=?");
#now get the strains listing for this source
$qh_getStrainId=$dbh->prepare("select strain_id from snp_strain_synonym where synonym_name=?");
$get_strain="select strain_id from snp_strain where strain_name=?";
$qh_getStrain=$dbh->prepare($get_strain);
$get_max_strainid="select max(strain_id) from snp_strain";
$qh_getMaxStrainid=$dbh->prepare($get_max_strainid);
$get_strain_source="select* from snp_strain_by_source where strain_id=? and source_id=?";
$qh_get_strain_source=$dbh->prepare($get_strain_source);
my $insertStrain="insert ignore into snp_strain values(?,?)";
$qh_insertStrain=$dbh->prepare($insertStrain);
my $insertStrainSyn="insert ignore into snp_strain_synonym values(?,?)";
$qh_insertStrainSyn=$dbh->prepare($insertStrainSyn);

# get the SNPs list if all SNPs
my $qh_snpId  =$dbh->prepare("select snpid from snp_position where chromosome_id=? and bp_position=?");
my $qh_snpIdRaonge =$dbh->prepare("select snpid from snp_position where chromosome_id=? and bp_position between ? and ?");
my $qh_getGeneSnp=$dbh->prepare("select distinct snpid from snp_transcript where gene_id=?");
my $qh_snpPos=$dbh->prepare("select chromosome_id,bp_position from snp_position where snpid=?");
   $query="select concat(ref_allele,'/',snp_allele) as allele,is_intergenic,mutation_type,is_CpG_site";
   $query.= " from snp_main where snpid=?";
my $qh_snpAllele=$dbh->prepare($query);
# get the SNPs list for SNPs located on transcript
$query="select t._loc_func_key,description ,strand,transcript_local_id,";
$query.=" loc_rank from snp_transcript t, snp_loc_func s ";
$query.=" where snpid=? and gene_id=?  and t._loc_func_key=s._loc_func_key";
my $qh_snpTx=$dbh->prepare($query);
my $qh_getTx=$dbh->prepare("select accession_id from cgd_transcripts_desc where transcript_local_id=?");
my $qh_getGene=$dbh->prepare("select accession_id from cgd_genes_desc where gene_id=?");
#now get the list of accessions for this SNP  
my $qh_getAccession =$dbh->prepare("select accession_id from snp_accession where snpid=? and snpid_id_type in(1,3) order by snpid_id_type limit 1");
my $qh_getSource=$dbh->prepare("select source_id from snp_source where source_name=?"); 
my $query="select _frame_key,PosInCDS,PosInProtein,ref_aa,ref_codon,snp_aa,snp_codon";
 $query.=" from snp_aminoacid where snpid=? and transcript_local_id=?";
 $qh_getaa=$dbh->prepare($query);
#Now get the genotype of this snp
my $getGenotype= " select distinct strain_id,genotype_allele from $genotype_table ";
   $getGenotype.=" where snpid=? ";
my $imgetGenotype= " select distinct strain_id,genotype_allele from $imputed_genotype_table ";
   $imgetGenotype.=" where snpid=? ";
$qh_getGeno=$dbh->prepare($getGenotype);
$qh_getImpGeno=$dbh->prepare($imgetGenotype);
my $qh_getMGI =$dbh->prepare("select marker_symbol from cgd_genes_ensembl_mgi where gene_id=?");
my $qh_getMGI_geneid =$dbh->prepare("select gene_id from cgd_genes_ensembl_mgi where marker_symbol=?");

$strainList_file="";$chrom_file="";
if($opt_c){$chrom_file=$opt_c;}if($opt_s){$strainList_file=$opt_s;} 
$strainList="";#$strainList=join(",",@strains);
if(-f $strainList_file){open(IN,"$strainList_file");
  if(IN){while($line=<IN>){chomp($line);push(@strains,$line);}}
}%strainMap=();@strains_index=();$i=0;
while(@strains>0){
 $strain=shift(@strains);$qh_getStrainId->execute($strain);$strain_id=0;
 if($qh_getStrainId->rows>0){($strain_id)=$qh_getStrainId->fetchrow_array();}
 $strainMap{$strain}=$strain_id;push(@strains_index,$strain);
}
if($opt_p){$getSource->execute($opt_p);$source_id=$getSource->fetchrow_array();}
$filename="$chrom_file-annotations.csv";open(IN,"$chrom_file");open(OUT ,">$filename");
$notfound="$chrom_file-notFount.csv";open(NOT,">$notfound");
if(!IN ||!OUT){print "bad file: $!\n";exit(1);}
print OUT "SNPID,Chromosome,NCBI37_Bp_position,Alleles,Mutation_type,Is_CpG_site,";
print OUT "Gene,Transcript,Function_class,aminoacidInfo,Feature_rank,Minor_allele%,Missing_genotype%";
foreach $strain (@strains_index){print OUT ",$strain";}
print OUT "\n";
%chrom_map=();$chrom_id=0;
while($gene=<IN>){
  chomp($gene);$qh_getMGI_geneid->execute($gene);
  if($qh_getMGI_geneid->rows<=0){print NOT "$gene\n";next;}
  ($gene_id)=$qh_getMGI_geneid->fetchrow_array();
  $qh_getGeneSnp->execute($gene_id);next if($qh_getGeneSnp->rows<=0);
  while(($snpid)=$qh_getGeneSnp->fetchrow_array()){
        %geno_map=();%allele_map=();$geno_data="";
        $qh_snpAllele->execute($snpid);
        $qh_snpPos->execute($snpid);
        ($chromosome_id,$position)=$qh_snpPos->fetchrow_array();
        $qh_getChr_name->execute($chromosome_id);
        ($chrom)=$qh_getChr_name->fetchrow_array();
        $qh_getAccession->execute($snpid);($snp_accession)= $qh_getAccession->fetchrow_array();
        #get the genotype
        $qh_getImpGeno->execute($snpid);%strain_count=();%has_geno=();%miss_geno=();
        if($qh_getImpGeno->rows>0){
          while(($strain_id,$geno_allele)=$qh_getImpGeno->fetchrow_array()){
             if($geno_allele=~/A|T|C|G/i){ $has_geno{$strain_id}=1;
                if(!exists($allele_map{uc($geno_allele)})){$allele_map{uc($geno_allele)}{$strain_id}=1;}
                else{$allele_map{uc($geno_allele)}{$strain_id}+=1;}
             }
             else{$miss_geno{$strain_id}=1;}
             $strain_count{$strain_id}=1;
             if(!exists($geno_map{$strain_id})){
                $geno_map{$strain_id}=$geno_allele if($geno_allele=~/A|T|C|G/i);}
             else{$geno_map{$strain_id}.=":$geno_allele" if($geno_allele=~/A|T|C|G/i);}
          }
        }
        $qh_getGeno->execute($snpid);#check if some strains w missing genotype in imputed table
                                     #have a genotype in the genotype table 
        %genoTable_map=();
        if($qh_getGeno->rows>0){
          while(($strain_id,$geno_allele)=$qh_getGeno->fetchrow_array()){
             if($geno_allele=~/A|T|C|G/i){ $has_geno{$strain_id}=1;
                if(!exists($allele_map{uc($geno_allele)})){$allele_map{uc($geno_allele)}{$strain_id}=1;}
                else{$allele_map{uc($geno_allele)}{$strain_id}+=1;}
             }
             else{#missing genotype strains
               $miss_geno{$strain_id}=1;
             }
             $strain_count{$strain_id}=1;
             if(!exists($genoTable_map{$strain_id})){
                 $genoTable_map{$strain_id}=$geno_allele if($geno_allele=~/A|T|C|G/i);}
             else{$genoTable_map{$strain_id}.=":$geno_allele" if($geno_allele=~/A|T|C|G/i);}
          }
        }
        $str_count=keys(%strain_count);
         $min=0;$max=0;$miss=0;
        while(($strain_id,$exists)=each(%miss_geno)){
           if(!exists($has_geno{$strain_id})){$miss+=1;}
        }
        while(($geno,$strains_count)=each(%allele_map)){
          if($min==0||$min>keys(%{$allele_map{$geno}})){$min=keys(%{$allele_map{$geno}});}
          if($max==0||$max<keys(%{$$allele_map{$geno}})){$max=keys(%{$allele_map{$geno}});}
        }
        #$miss=$str_count-($max+$min);
        $min="$min/$str_count";$miss="$miss/$str_count";
        $allele="";
        $mut_type="Transition";$CpG_site="No";
        ($allele,$is_intergenic,$mutation_type,$is_CpG_site)=$qh_snpAllele->fetchrow_array();
        ($ref_allele,$other_allele)=split(/\//,$allele);
        foreach $strain (@strains_index){
                 $strain_id=$strainMap{$strain};
                if($geno_map{$strain_id}){$geno_data.=",".$geno_map{$strain_id};}
                elsif($strain_id==7){$geno_data.=",$ref_allele";}
                elsif($genoTable_map{$strain_id}){$geno_data.=",".$genoTable_map{$strain_id};}
                else{$geno_data.=",N";}
        }
        if($allele ne ""){
           $mut_type="Transversion" if($mutation_type==2);
           $CpG_site="Yes" if($is_CpG_site==1);$description="Intergenic";$loc_rank="-1";
           #get snp location
            if($is_intergenic==0){ $qh_snpTx->execute($snpid,$gene_id);
               while(($loc_id,$description ,$strand,$transcript_local_id,$loc_rank)=$qh_snpTx->fetchrow_array()){
                 $qh_getTx->execute($transcript_local_id);($tx)=$qh_getTx->fetchrow_array();
                 $codon="--";
                 if(($loc_id>0)&&($loc_id!=3)&&($loc_id!=4)){
                    $qh_getaa->execute($snpid,$transcript_local_id)or die mysql_error();
                    ($frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon)=$qh_getaa->fetchrow_array();
                    $codon="PosInCondon:$frame-PosInCDS:$PosInCDS-codonPosInProtein:$PosInProtein";
                    $codon.="-refAminoAcid:$ref_aa/$ref_codon-otherAminoAcid:$snp_aa/$snp_codon";
                  }
                  print OUT "$snp_accession,$chrom,$position,$allele,$mut_type,$CpG_site,$gene,";
                  print OUT "$tx,$description,$codon,$loc_rank,$min,$miss$geno_data\n";
                }
            }
           else{$gene="--";$tx="--";
                print OUT "$snp_accession,$chrom,$position,$allele,$mut_type,$CpG_site,$gene,";
                print OUT "$tx,$description,--,$loc_rank,$min,$miss$geno_data\n";
            }
       
         }
     }
}

 print "Program complete\n";
 exit(0);
      
