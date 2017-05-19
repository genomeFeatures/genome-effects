#!/usr/bin/perl
################################################################
q{
Author: Lucie Hutchins, Scientific Software Engineer
Date : August 2012

Project Description:
I would like to take every unique exon in the mouse transcript database and intersect it with the SNPs for the 8 founder strains (7 when you consider that B6 is the reference).
I would like the output formatted in this specific way:
Exon ID<tab>All transcripts that include this exon in a semi-colon separated list<tab>chromosome<tab>start<tab>stop<tab>strand<tab>SNP list<end of line>
The SNP list I would like formatted like this:  
position:strain=base,strain=base;position:strain=base,strain=base
So for example, if an exon had a SNP at position 35678987 that was different in A/J and C3H/HeJ, and another SNP at position 35679045 that was different in PWK and CAST, the SNP list would look like this:

35678987:A/J=A,C3H/HeJ=A;35679045:PWK=C,CAST=C
I'd like this as soon as possible, as I want it for some new work with Elissa on DO data.  Let me know if it's harder than I think or if it isn't clear.

};
################################################################

use DBI;
use vars qw ($opt_h $opt_s $opt_v $opt_g $opt_t);
use Getopt::Std;
use POSIX;
getopts('hs:v:g:t:');
if($opt_h||(!$opt_s&&!$opt_g)) {
    print <<HELP;

 This script generates SNPs annotation that overlap a given genomic region (exon,..)
 for provided strains list

Usage:

   perl getSNPsOnFeatures.pl -s <strain_file> -v <organism version> -t <feature type>

Arguments:
   -h  displays this help message
   -s  A list containing name(s) of selected strains (multiple names are separated with a commas)
   -g  A flag specifying the strain group (only use this option if the -s option is not used) :
         1-> Classical strains ; 
         2-> Wild Derived
         3-> Wild
         4-> CC Founders
   -t  A flag to specify the feature type:
         ex -> exons
         tx -> transcripts

Example: perl getSNPsOnFeatures.pl -s strainlist -v mm9 -g 4 -t ex

HELP
exit;
}
my $dbname ='cgd_snpdb';my $host ='cgd-dev.jax.org';
my $user ='pup';my $pwd  ='puppass';
my $g_dbname ='graber_transcriptdb';my $g_host ='demon.jax.org';
$g_dbname="graber_transcriptdb";$g_host="harlequin.jax.org"; #$user="lnh";$pass="lucie98";

$organism_version="mm9";  # default organism version
$strain_group=4;         # CC Founders is the default group
$feature_type="ex";     # exons are the default features
$ref_strain_id=7;       # B6 is the reference strain
my $dbh = DBI->connect("DBI:mysql:$dbname:$host",$user,$pwd)or die "Can't connect to database!\n";
my $gdbh = DBI->connect("DBI:mysql:$g_dbname:$g_host",$user,$pwd)or die "Can't connect to database!\n";

######################### TRANSCRIPT DB QUERIES ##################################################
# get chromosomes list
$query="select distinct chromosome_name, c.chromosome_id from chromosome c, chromosome_data d ";
$query.="where d.organism_version_id=? and d.chromosome_id=c.chromosome_id";
my $gqh_getChr =$gdbh->prepare($query);
#get exon list
$query="select exon_id,strand,exon_start,exon_end from exon where organism_version_id=? and chromosome_id=?";
my $gqh_getExons=$gdbh->prepare($query);
#get transcript list
$query="select transcript_id,strand,tx_start,tx_end from transcript where organism_version_id=? and chromosome_id=?";
my $gqh_getTranscripts=$gdbh->prepare($query);

#get gene list for this transcript
$query="select distinct gene_name from gene_by_annotation where transcript_id =? and organism_version_id=? ";
my $gqh_getTxGenes=$gdbh->prepare($query);

$query="select distinct transcript_id from transcript_exon where exon_id=? and organism_version_id=? and gene_prediction_id in(?)";
my $gqh_getExonsTxids=$gdbh->prepare($query);
$query="select distinct transcript_name from transcript_by_annotation ";
$query.=" where transcript_id in(?) and organism_version_id=? and gene_prediction_id in(?) order by transcript_name";
my $gqh_getTxAccessions=$gdbh->prepare($query);
#get organism version id
my $getOrgV="select organism_version_id from organism_version where ucsc_db=?";
my $gqh_getOrgV = $gdbh->prepare($getOrgV)or die "Couldn't prepare statement: ".$dbh->errstr;

######################### CGD QUERIES ############################################################
$imputed_genotype_table="snp_imputed"; $genotype_table="snp_genotype";
# get chromosomes list
my $qh_getChr =$dbh->prepare("select chromosome_id from snp_chromosome where chromosome_name=?");
my $qh_getChr_name =$dbh->prepare("select chromosome_name from snp_chromosome where chromosome_id=?");
#get the strain group
$query="select distinct strain_name,s.strain_id from snp_strain s, snp_strain_by_group sp ";
$query.=" where sp.group_id=? and sp.strain_id=s.strain_id ";
my $qh_getGroupStrains=$dbh->prepare($query);
$qh_snpAllele=$dbh->prepare("select ref_allele from snp_main where snpid=?");
#get the strain id given the strain name
$qh_getStrainId=$dbh->prepare("select strain_id from snp_strain_synonym where synonym_name=?");
$qh_getStrain=$dbh->prepare("select strain_id from snp_strain where strain_name=?");
# get the SNPs list if all SNPs
my $qh_snpIdRaonge =$dbh->prepare("select snpid, bp_position from snp_position where chromosome_id=? and bp_position between ? and ?");
#get the genotype of this snp
my $getGenotype= " select distinct strain_id,genotype_allele from $genotype_table ";
   $getGenotype.=" where snpid=? and genotype_allele in ('A','T','C','G')";
my $imgetGenotype= "select distinct strain_id,genotype_allele from $imputed_genotype_table ";
   $imgetGenotype.=" where snpid=? and genotype_allele in ('A','T','C','G')";
$qh_getGeno=$dbh->prepare($getGenotype);$qh_getImpGeno=$dbh->prepare($imgetGenotype);
#############################################

$organism_version=$opt_v if($opt_v);$org_vid=0;$organism_version=~s/\s+//g;
if($organism_version ne "mm9"){print "We only have the SNPs data for mm9\n";exit(1);}
$gqh_getOrgV->execute("$organism_version");
if($gqh_getOrgV->rows>0){($org_vid)=$gqh_getOrgV->fetchrow_array();}
@strains=();if($opt_g){$strain_group=$opt_g;}
$strainList_file="";$chrom_file="";$i=0;%strainIdMap=();%strainGenoMap=();
if($opt_s){@strains=split(",",$opt_s);
  while($strain=shift(@strains)){
        $qh_getStrainId->execute($strain);$strain_id=0;
        if($qh_getStrainId->rows>0){($strain_id)=$qh_getStrainId->fetchrow_array();}
        $strainIdMap{$strain_id}=$strain;
   }   
}
if(!$opt_s){
  $qh_getGroupStrains->execute($strain_group);
  while(($strain,$strain_id)=$qh_getGroupStrains->fetchrow_array()){$strainIdMap{$strain_id}=$strain;}
}if($opt_t){$feature_type=$opt_t;}
$filename="$organism_version-$feature_type-SNPs.txt";
$genefile="$organism_version-$feature_type-SNPs-transcript-genes-table.txt";
#get the list of chromosome of this organism version
$gqh_getChr->execute($org_vid);
if($gqh_getChr->rows>0){
   open(OUT ,">$filename");$feature="exonID"; if(!OUT){print "bad file: $!\n";exit(1);}
   if($feature_type=~/tx/i){$feature="transcriptID";
      print OUT "transcriptList\tgeneLists\tchromosome\tchromStart\tchromEnd\tstrand\n";
    }else{ print OUT "$feature\ttranscriptList\tchromosome\tchromStart\tchromEnd\tstrand\tSNP_list\n";}
   open(OUT2 ,">$genefile");print OUT2 "transcript\tgenes\n";$count=1;
  while(($chrom,$chr_id)=$gqh_getChr->fetchrow_array()){
     next if($chr_id>22);#the SNP database does not have chr_id>22
     if($chrom=~/^\s*M\s*$/i){$chrom="MT";}if($chrom=~/^\s*U\s*$/i){$chrom="UN";}
     #get the features of this chromosome
     if($feature_type=~/ex/i){ %tx_map=();
        $gqh_getExons->execute($org_vid,$chr_id);
        while(($exon_id,$strand,$feature_start,$feature_end)=$gqh_getExons->fetchrow_array()){
             #get all the transcripts that contain this exon id
             $gqh_getExonsTxids->execute($exon_id,$org_vid,"19");$tex_ids="";$transcriptAcc="";$maintx="";
             if($gqh_getExonsTxids->rows>0){
                while(($tx_id)=$gqh_getExonsTxids->fetchrow_array()){$gene_list="";$transcriptAcc="";$this_tx="";
                     $gqh_getTxAccessions->execute($tx_id,$org_vid,"19");
                     while(($thistx_id)=$gqh_getTxAccessions->fetchrow_array()){
                         if($transcriptAcc eq ""){$transcriptAcc="$thistx_id";}else{$transcriptAcc.=",$thistx_id";}}
                     if($maintx eq ""){$maintx=$transcriptAcc;}else{$maintx.=",$transcriptAcc";}
                     if(!exists($tx_map{$tx_id})){
                        $gqh_getTxGenes->execute($tx_id,$org_vid);
                        while(($gene)=$gqh_getTxGenes->fetchrow_array()){
                           if($gene_list eq ""){$gene_list="$gene";}else{$gene_list.=",$gene";}}
                        print OUT2 "$transcriptAcc\t$gene_list\n";$tx_map{$tx_id}=1;
                     }
                 }
                $qh_snpIdRaonge->execute($chr_id,$feature_start,$feature_end);$snplist="";
                while(($snpid,$bp_position)=$qh_snpIdRaonge->fetchrow_array()){$qh_snpAllele->execute($snpid);
                     if($qh_snpAllele->rows>0){($ref_allele)=$qh_snpAllele->fetchrow_array();
                       #get the genotype of selected strain #
                       $qh_getImpGeno->execute($snpid);%has_geno=();%allele_map=();
                       if($qh_getImpGeno->rows>0){
                         while(($strain_id,$geno_allele)=$qh_getImpGeno->fetchrow_array()){
                            if(exists($strainIdMap{$strain_id})){
                              $has_geno{$strain_id}=1;$allele_map{uc($geno_allele)}{$strain_id}=1;}
                         }
                       }#now check if all the selected strains have a genotype 
                       if(keys(%strainIdMap)!=keys(%has_geno)){$qh_getGeno->execute($snpid);
                         if($qh_getGeno->rows>0){
                            while(($strain_id,$geno_allele)=$qh_getGeno->fetchrow_array()){ 
                              if(exists($strainIdMap{$strain_id})){$allele_map{uc($geno_allele)}{$strain_id}=1;}
                            }
                          }
                       }#get all trains with different genotype than reference strain
                       $ref_allele=uc($ref_allele);%strainGenoMap=();
                       $ref_strain=$strainIdMap{$ref_strain_id};
                       if(keys(%allele_map)>0){
                          while(($geno_allele)=each(%allele_map)){
                                while(($strain_id)=each(%{$allele_map{$geno_allele}})){
                                    if(!($geno_allele=~/$ref_allele/i)){$strainGenoMap{$strain_id}="$geno_allele";}}
                           }$snpdata="";
                           while(($strain_id)=each(%strainGenoMap)){
                               $this_strain_geno=$strainGenoMap{$strain_id};$this_strain_name=$strainIdMap{$strain_id};
                               if($snpdata eq ""){
                                  $snpdata="$bp_position:$ref_strain=$ref_allele,$this_strain_name=$this_strain_geno";}
                               else{$snpdata.=",$this_strain_name=$this_strain_geno";}
                           }
                           if($snpdata ne ""){if($snplist eq ""){$snplist="$snpdata";}else{$snplist.=";$snpdata";}}
                       }
                    }# end of  if($qh_snpAllele->rows>0)
                 } #END of while(($snpid,$bp_position)=$qh_snpIdRaonge->
                 print OUT "$exon_id\t$maintx\t$chrom\t$feature_start\t$feature_end\t$strand\t$snplist\n";
                # ++$count; last if($count>2000);
               }# end of if($gqh_getExonsTxids->rows>0)
        } #end of  while(($exon_id,$strand,$feature_start,$feature_end)=$gqh_getExons->
     } #end of  if($feature_type=~/ex/i){
  }
}
print "Program complete\n";
exit(0);
