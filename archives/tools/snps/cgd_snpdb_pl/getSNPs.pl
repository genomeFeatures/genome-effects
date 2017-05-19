#!/usr/bin/perl
################################################################
q{
This script generates SNP annotations for selected strains

Author: Lucie Hutchins, Scientific Software Engineer
Date : October 2011

};
################################################################

use DBI;
use vars qw ($opt_h $opt_s $opt_r $opt_g $opt_d $opt_p $opt_c);
use Getopt::Std;

getopts('hs:r:gd:p:c:');
if($opt_h||!$opt_s ) {
    print <<HELP;

 This script generates SNPs annotation between provided strains

Usage:

   perl getSNPs.pl -s <list_of_strains_separated_w_commas> [-r <ref_strain_name>] [-g ] [-d <0,1,2>] [-p <strain_provider>][-c <chromosome>]

Arguments:
   -h  displays this help message
   -s  The name(s) of selected strains (multiple names are separated with a commas
   -r  The reference strain name (Black 6 (C57BL/J6) is the default strain)
   -g  A flag to specified to only include SNPs on transcript
   -d  A flag to filter rows based on genotype similarity :
        0-> display all cases ; 
        1-> only display cases where the genotype is different between selected strains
        2-> only display cases where the genotype is the same between selected strains
   -p  name of the snp provider (source)
   -c  chromosome name

Example: perl getSNPs.pl -s DBA/1J,DBA/2J
HELP
exit;
}


my $dbname ='cgdsnpdb';
my $host ='cgd';
my $user ='pup';
my $pwd  ='puppass';

my $dbh = DBI->connect("DBI:mysql:$dbname:$host",$user,$pwd)
                            or die "Can't connect to database!\n";

$genotype_table="snp_genotype"; $imputed_genotype="snp_imputed";
# get chromosomes list
my $getChr ="select distinct chromosome_id,chromosome_name from snp_chromosome";

#now get the strains listing for this source
$get_strains="select strain_id from snp_strain_synonym where synonym_name=?";
$qh_getStrainId=$dbh->prepare($get_strains);
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
my $snplistQuery  ="select snpid,bp_position,black6_allele,snp_allele "; 
   $snplistQuery .=" from snp_main where chromosome_id=? and is_public=1 ";

# get the SNPs list for SNPs located on transcript
my $txsnplistQuery  ="select snpid,bp_position,black6_allele,snp_allele "; 
   $txsnplistQuery .=" from snp_main where chromosome_id=? and is_public=1 and is_intergenic=0";

#now get the list of accessions for this SNP  
my $getAccession =" select accession_id from snp_accession where snpid=? order by snpid_id_type";
  
#Now get the genotype of this snp
my $getGenotype= " select genotype_allele from snp_genotype ";
   $getGenotype.=" where snpid=? and strain_id=? and genotype_allele in('A','T','C','G')";

#Now get the genotype of this snp from this source
my $getSourceGenotype= " select genotype_allele from snp_genotype ";
   $getSourceGenotype.=" where snpid=? and strain_id=? and source_id=? and genotype_allele in('A','T','C','G')";

#Now get the genotype of this snp
my $getImpGenotype= " select genotype_allele from snp_imputed";
   $getImpGenotype.=" where snpid=? and strain_id=? ";
#now get the sources of this SNP
my $getSNPSource= " select distinct source_name from snp_source s, snp_by_source sc ";
   $getSNPSource.=" where sc.snpid=? and sc.source_id=s.source_id order by source_name "; 

#get the exon of this transcript containing the SNP
my $getExondata="select accession_id,exon_start,exon_end,exon_rank from cgd_transcript_exons te, ";
   $getExondata.="cgd_exons_desc ed where transcript_local_id=? and exon_start<=? ";
   $getExondata.=" and exon_end>=? and te.exon_id=ed.exon_id";

#now get tx accession 
my $getTxAccession =" select accession_id from cgd_transcripts_desc where transcript_local_id=? ";
#now get mgi symbol 
my $getMgiAccession =" select marker_symbol from cgd_genes_ensembl_mgi where gene_id=? ";
my $getSource="select source_id from snp_source where source_name=?";
my $isThisSNPsource="select source_id from snp_by_source where snpid=? and source_id=?";

my $getSNPLocation="select t._loc_func_key,description ,strand,transcript_local_id,gene_id ";
    $getSNPLocation.=" from snp_transcript t, snp_loc_func f where snpid=?  ";  
    $getSNPLocation.=" and t._loc_func_key=f._loc_func_key order by gene_id"; 

#get the transcript(s) data containing this SNP
my $getTxdata="select t.transcript_start,t.transcript_end,gd.accession_id, gd.gene_type  ";
   $getTxdata.=" from cgd_transcripts t,cgd_genes_desc gd where t.transcript_local_id=? and t.gene_id=? ";
   $getTxdata.="  and t.gene_id=gd.gene_id";

my $getCoding=" select _frame_key,PosInCDS,PosInProtein,black6_aa,b6_codon,snp_aa,snp_codon ";
    $getCoding.=" from snp_aminoacid where snpid=? and transcript_local_id=?";

my $genotypeBySource= $dbh->prepare($getSourceGenotype);
my $isSNPInSource  =$dbh->prepare($isThisSNPsource)or die "Couldn't prepare statement: " . $dbh->errstr;
my $getSource  =$dbh->prepare($getSource)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpLoc     = $dbh->prepare($getSNPLocation)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpCoding = $dbh->prepare($getCoding)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpTx      = $dbh->prepare($getTxdata)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpMgi     = $dbh->prepare($getMgiAccession)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpTxAc     = $dbh->prepare($getTxAccession)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpExon     = $dbh->prepare($getExondata)or die "Couldn't prepare statement: " . $dbh->errstr;

my $snpchr   = $dbh->prepare($getChr)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpstr   = $dbh->prepare($get_strains)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snplist  = $dbh->prepare($snplistQuery)or die "Couldn't prepare statement: " . $dbh->errstr;
my $txsnplist  = $dbh->prepare($txsnplistQuery)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpac    = $dbh->prepare($getAccession)or die "Couldn't prepare statement: " . $dbh->errstr;

my $snpgeno  = $dbh->prepare($getGenotype)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpImpgeno  = $dbh->prepare($getImpGenotype)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpSource  = $dbh->prepare($getSNPSource)or die "Couldn't prepare statement: " . $dbh->errstr;

my($strainList,$ref_strain);
$ref_strain="C57BL/6J"; #default strain
$chrom="";$geno_filter=0;
if($opt_c){$chrom=$opt_c;}
if($opt_s){$strainList=$opt_s;} 
if($opt_r){$ref_strain=$opt_r;}
if($opt_d){$geno_filter=$opt_d;}

#$strain1="DBA/1J"; $strain2="DBA/2J";

%strainMap=();
@strains=split(",",$strainList);
@strains_index=();$i=0;
while(@strains>0){
 $strain=shift(@strains);
 $snpstr->execute("$strain");
 $strainMap{$strain}=$snpstr->fetchrow_array();
 push(@strains_index,$strain);
}
$source_id=0;
if($opt_p){
   $getSource->execute($opt_p);
   $source_id=$getSource->fetchrow_array();
}
#exit(0);
# get the list of all the chromosomes
$snpchr->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
$scount=0;
if($snpchr->rows > 0){ #process each chromosome 
   while(($chromosome_id,$chrom_name)= $snpchr->fetchrow_array()){
       #my($chromosome_id,$chrom_name)=@row;
       $chrom=~s/\s*//;$chrom_name=~s/\s*//;
       if($chrom ne ""){next if(lc($chrom) ne lc($chrom_name));}
       $filename="chr$chrom_name"."_SNPs.txt";
       print "Processing  $chrom_name\n";
       #get the snp list
       $theList= $snplist;
       if($opt_g){$theList=$txsnplist;}
       $theList->execute($chromosome_id) or die "Can't execute query: " . $dbh->errstr . "\n";
       if($theList->rows>0){
           open(OUT ,">$filename") or die "can't open $filename\n";
           print OUT "SNPID\tChromosome\tNCBI37_Bp_position\tAlleles\t";
           print OUT "Exonid\tExon_rank\tExon_start\tExon_end\tTranscript_id\tTranscript_start\t";
           print OUT "Transcript_end\tgene_strand\tEnsembl_geneid\tMgi_symbol\tGene_type\t";
           print OUT "snpFunction_class\treference_codon/Consensus_codon\treference_aminoacid/Consensus_aminoacid\t";
           print OUT "snp_frame\tsnpPositionInCds\tcodonPositionInProtein\tSource(s)";
          foreach $strain (@strains_index){print OUT "\t$strain";}
          if(!$opt_d){ print OUT "\tsameGenotype\n";}
          else{print OUT "\n";}
          $count=0;
           while(($snpid,$bp_position,$black6_allele,$snp_allele)=$theList->fetchrow_array()){ 
             #get the snp data    
            if($source_id>0){ next if(!($isSNPInSource->execute($snpid,$source_id)));}
             next if(!($snpac->execute($snpid)));
             $accession_id="";
             while($id=$snpac->fetchrow_array()){
               if(($id=~/^\s*rs/)||($id=~/JAX/i)){
                  if($accession_id eq ""){$accession_id="$id";}
                  else{$accession_id .=",$id";}
                }
             }
             #get the source(s0 of this SNP
             $source=""; 
             if($opt_p){
                 $source=$opt_p;
             }
             else{ $sc=0;$snpSource->execute($snpid);
                   while(($source_name)=$snpSource->fetchrow_array()){
                        if($sc==0){$source=$source_name;++$sc;}
                        else{$source.=",$source_name";}
                   }
             }
             #get SNP location
             #get the genotype of every strain on the list
             $genodata="";$i=0;$samegeno="Yes";$allele="";$prevallele="";++$count;$mis=0;
             ###########################################################################
             foreach $strain (@strains_index){
                 $allele="";
                 if($source_id>0){
                    $genotypeBySource->execute($snpid,$strainMap{$strain},$source_id);
                    if($genotypeBySource->rows>0){$allele=($genotypeBySource->fetchrow_array())[0];}
                 }
                 else{
                    $snpImpgeno->execute($snpid,$strainMap{$strain})or die "Can't execute query: " . $dbh->errstr . "\n";
                    if($snpImpgeno->rows>0){ #this strain has genotype data in the imputed genotype table for this SNP
                       $allele=($snpImpgeno->fetchrow_array())[0];
                    }
                    else{ #a miss in the imputed genotye table, check the main genotype table
                       $snpgeno->execute($snpid,$strainMap{$strain})or die "Can't execute query: " . $dbh->errstr . "\n";
                       if($snpgeno->rows>0){$allele=($snpgeno->fetchrow_array())[0];}
                    }
                }
                if(($allele ne "")&&($prevallele ne "")&&(lc($prevallele) ne lc($allele))){$samegeno="No";}
                $prevallele=$allele;$genodata.="\t$allele";if($allele eq ""){$mis=1;}
             }
             if($opt_d==1){ #only display cases with different genotype
               next if($samegeno eq "Yes");
             }
             elsif($opt_d==2){ #only display cases with same genotype
               next if($samegeno eq "No");
             }
             next if($mis);
             next if(!($snpLoc->execute($snpid)));
             if($snpLoc->rows>0){ # snp located on transcript
                while(($_loc_func_key,$functionClass,$gene_strand,$transcript_local_id,$gene_id)=$snpLoc->fetchrow_array()){
                   #get tx data;
                   if($gene_strand==-1){$gene_strand="-";}
                   else{$gene_strand="+";}
                   $snpTx->execute($transcript_local_id,$gene_id);
                   ($transcript_start,$transcript_end,$ensembl_geneid,$gene_type)=$snpTx->fetchrow_array();
                   #get ensembl tx id, get mgi gene symbol
                   $snpMgi->execute($gene_id);$mgi_gene="na";$tx_id="na";
                   if($snpMgi->rows>0){($mgi_gene)=$snpMgi->fetchrow_array();}
                   $snpTxAc->execute($transcript_local_id);
                   if($snpTxAc->rows>0){($tx_id)=$snpTxAc->fetchrow_array();}
                   if(($_loc_func_key<=0)&&($_loc_func_key!=-3)){
                      print OUT "$accession_id\t$chrom_name\t$bp_position\t$black6_allele/$snp_allele\t";                
                      print OUT "na\tna\tna\tna\t$tx_id\t$transcript_start\t";
                      print OUT "$transcript_end\t$gene_strand\t$ensembl_geneid\t$mgi_gene\t$gene_type\t";
                      print OUT "$functionClass\tna\tna\tna\tna\tna\t";
                      print OUT " $source\t$genodata";
                      if(!$opt_d){print OUT "\t$samegeno\n";}
                      else{print OUT "\n";}
                   }
                   else{ #snp locate on exon 
                      #get exon data
                      $snpExon->execute($transcript_local_id,$bp_position,$bp_position);
                      ($Exonid,$Exon_start,$Exon_end,$Exon_rank)=$snpExon->fetchrow_array();             
                      print OUT "$accession_id\t$chrom_name\t$bp_position\t$black6_allele/$snp_allele\t";                
                      print OUT "$Exonid\t$Exon_rank\t$Exon_start\t$Exon_end\t$tx_id\t$transcript_start\t";
                      print OUT "$transcript_end\t$gene_strand\t$ensembl_geneid\t$mgi_gene\t$gene_type\t";
                      if(($_loc_func_key>0)&&($_loc_func_key!=3)&&($_loc_func_key!=4)){
                          #get coding info if snp coding
                          $snpCoding->execute($snpid,$transcript_local_id);
                          ($snp_frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,
                              $snp_aa,$snp_codon)=$snpCoding->fetchrow_array();
                           print OUT "$functionClass\t$ref_codon/$snp_codon\t$ref_aa/$snp_aa\t";
                           print OUT "$snp_frame\t$PosInCDS\t$PosInProtein\t";
                           print OUT " $source\t$genodata";
                           if(!$opt_d){print OUT "\t$samegeno\n";}
                           else{print OUT "\n";}
                       }
                      else{ #SNP ON non coding exon
                          print OUT "$functionClass\tna\tna\tna\tna\tna\t";
                          print OUT " $source\t$black6_allele\t$genodata";
                          if(!$opt_d){print OUT "\t$samegeno\n";}
                          else{print OUT "\n";}
                       }
                   }
                
                }
              }
             else{ #intergenic snp
                 print OUT "$accession_id\t$chrom_name\t$bp_position\t$black6_allele/$snp_allele\t";                
                 print OUT "na\tna\tna\tna\tna\tna\tna\tna\tna\tna\tna\t";
                 print OUT "Intergenic\tna\tna\tna\tna\tna\t$source\t$genodata";
                 if(!$opt_d){print OUT "\t$samegeno\n";}
                 else{print OUT "\n";}
             }
             ###########################################################################  
            # print OUT "$snpdata\t$genodata\t$samegeno\n";
             
           }
           close(OUT);
           #now compress this file;
           $zipfile="$filename"."."."tar"."."."gz";
           $command="tar -cvzf $zipfile $filename";
           system($command); if(-d "temp"){system("mv $filename temp");}
          ++$scount; 
      }
   }
 }

 print "Program complete\n";
 exit(0);
      
