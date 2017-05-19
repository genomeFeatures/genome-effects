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
$source_id=16;$genotype_table="snp_imputed"; 
# get chromosomes list
my $qh_getChr =$dbh->prepare("select chromosome_id from snp_chromosome where chromosome_name=?");
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
my $qh_snpIdRange =$dbh->prepare("select snpid from snp_position where chromosome_id=? and bp_position between ? and ?");
my $qh_snpPos=$dbh->prepare("select chromosome_id,bp_position from snp_position where snpid=?");
   $query="select concat(ref_allele,'/',snp_allele) as allele,is_intergenic,mutation_type,is_CpG_site";
   $query.= " from snp_main where snpid=?";
my $qh_snpAllele=$dbh->prepare($query);
# get the SNPs list for SNPs located on transcript
$query="select t._loc_func_key,description ,strand,transcript_local_id,";
$query.=" loc_rank,gene_id from snp_transcript t, snp_loc_func s ";
$query.=" where snpid=? and t._loc_func_key=s._loc_func_key";
my $qh_snpTx=$dbh->prepare($query);
my $qh_getTx=$dbh->prepare("select accession_id from cgd_transcripts_desc where transcript_local_id=?");
my $qh_getGene=$dbh->prepare("select accession_id from cgd_genes_desc where gene_id=?");
#now get the list of accessions for this SNP  
my $qh_getAccession =$dbh->prepare("select accession_id from snp_accession where snpid=? and snpid_id_type in(1,3) order by snpid_id_type limit 1");
my $qh_getSource=$dbh->prepare("select source_id from snp_source where source_name=?"); 
my $query="select _frame_key,PosInCDS,PosInProtein,ref_aa,ref_codon,snp_aa,snp_codon";
 $query=" from snp_aminoacid where snpid=? and transcript_local_id=?";
 $qh_getaa=$dbh->prepare($query);
#Now get the genotype of this snp
my $getGenotype= " select strain_id,genotype_allele from $genotype_table ";
   $getGenotype.=" where snpid=? and source_id=?";
$qh_getGeno=$dbh->prepare($getGenotype);
my $qh_getMGI =$dbh->prepare("select marker_symbol from cgd_genes_ensembl_mgi where gene_id=?");

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
 print "Processing $strain -- $strain_id\n";
}
exit(0);
if($opt_p){$getSource->execute($opt_p);$source_id=$getSource->fetchrow_array();}
$filename="$chrom_file-annotations.csv";open(IN,"$chrom_file");open(OUT ,">$filename");
if(!IN ||!OUT){print "bad file: $!\n";exit(1);}
print OUT "SNPID\tChromosome\tNCBI37_Bp_position\tAlleles\tMutation_type\tIs_CpG_site\t";
print OUT "Gene\tTranscript\tFunction_class\taminoacidInfo\tFeature_rank\tMinor_allele%\tMissing_genotype%";
foreach $strain (@strains_index){print OUT "\t$strain";}
print OUT "\n";
%chrom_map=();$chrom_id=0;
while($snpdata=<IN>){
  chomp($snpdata);($chrom,$position)=split(":",$snpdata);
  $chrom=~s/chr//i; $chrom=~s/\s//g;$chrom_id=0;
  if(!exists($chrom_map{$chrom})){
     $qh_getChr->execute($chrom);
     ($chrom_id)=$qh_getChr->fetchrow_array();
     $chrom_map{$chrom}=$chrom_id;
   }else{$chrom_id=$chrom_map{$chrom};}
  if($chrom_id>0){
     $qh_snpId->execute($chrom_id,$position);
     if($qh_snpId->rows<=0){ print "$chrom_id,$position not found\n";}
     else{%geno_map=();%allele_map=();$geno_data="";
       ($snpid)=$qh_snpId->fetchrow_array();$qh_snpAllele->execute($snpid);
        $qh_getAccession->execute($snpid);($snp_accession)= $qh_getAccession->fetchrow_array();
        #get the genotype
        $qh_getGeno->execute($snpid,$source_id);%strain_count=();
        if($qh_getGeno->rows>0){
          while(($strain_id,$geno_allele)=$qh_getGeno->fetchrow_array()){
             if(!exists($allele_map{uc($geno_allele)})){$allele_map{uc($geno_allele)}=1;}
             else{$allele_map{uc($geno_allele)}+=1;}$strain_count{$strain_id}=1;
             if(!exists($geno_map{$strain_id})){$geno_map{$strain_id}=$geno_allele;}
             else{$geno_map{$strain_id}=":$geno_allele";}
          }
        }
        $str_count=keys(%strain_count);
        $min=0;$max=0;$miss=0;
        while(($geno,$strain_count)=each(%allele_map)){
          if($min==0||$min>$strain_count){$min=$strain_count;}
          if($max==0||$max<$strain_count){$max=$strain_count;}
        }
        $miss=$str_count-($max+$min);$min="$min/$str_count";$miss="$miss/$str_count";
        foreach $strain (@strains_index){
                 $strain_id=$strainMap{$strain};
                if($geno_map{$strain_id}){$geno_data.=",".$geno_map{$strain_id};}
                else{$geno_data.=",-";}
        }
        if($qh_snpAllele->rows>0){ $mut_type="Transition";$CpG_site="No";
           ($allele,$is_intergenic,$mutation_type,$is_CpG_site)=$qh_snpAllele->fetchrow_array();
           $mut_type="Transversion" if($mutation_type==2);
           $CpG_site="Yes" if($is_CpG_site==1);$description="Intergenic";$loc_rank="-1";
           #get snp location
            if($is_intergenic==0){ $qh_snpTx->execute($snpid);
               while(($loc_id,$description ,$strand,$transcript_local_id,$loc_rank,$gene_id)=$qh_snpTx->fetchrow_array()){
                 $gene="";
                 $qh_getTx->execute($transcript_local_id);($tx)=$qh_getTx->fetchrow_array();
                 $qh_getMGI->execute($gene_id);($gene)=$qh_getMGI->fetchrow_array() if($qh_getMGI->rows>0);
                 if($gene eq ""){$qh_getGene->execute($gene_id);($gene)=$qh_getGene->fetchrow_array();}
                 $codon="--";
                 if(($loc_id>0)&&($loc_id!=3)&&($loc_id!=4)){
                    $qh_getaa->execute($snpid,$transcript_local_id);
                    ($frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon)=$qh_getaa->fetchrow_array();
                    $codon="PosInCondon:$frame-PosInCDS:$PosInCDS-codonPosInProtein:$PosInProtein";
                    $condon.="-refAminoAcid:$ref_aa/$ref_codon-otherAminoAcid:$snp_aa/$snp_codon";
                  }
                  #if($chrom_id==11&&$position==36841359){
                    print OUT "$snp_accession,$chrom,$position,$allele,$mut_type,$CpG_site,$gene,";
                    print OUT "$tx,$description,$codon,$loc_rank,$min,$miss$geno_data\n";
                 #}
                }
            }
           else{$gene="--";$tx="--";
                print OUT "$snp_accession,$chrom,$position,$allele,$mut_type,$CpG_site,$gene,";
                print OUT "$tx,$description,--,$loc_rank,$min,$miss$geno_data\n";
            }
       
         }
     }
  }
}

 print "Program complete\n";
 exit(0);
      
