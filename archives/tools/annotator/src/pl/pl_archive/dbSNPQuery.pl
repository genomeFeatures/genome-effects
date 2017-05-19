#!/usr/bin/perl

######################################################################
## This script generates SNPs file by chromosome
## from the snp database
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : April 2011
#   last update : May 2013
#   
#
#####################################################################
use DBI;
use vars qw ($opt_h $opt_d $opt_s);
use Getopt::Std;
use POSIX;
getopts('hd:s:');
if($opt_h) {
    print <<HELP;

 This script generates SNPs file by chromosome for the specified SNP source.
 The chromosome file names have the format: chrxx"."_SourceDate_snps.txt"
  for example (chr1_Imputed2012-snps.txt)

Usage:

   perl dbSNPQuery.pl [-d resultPath] [-s strainsList]

Arguments:
   -h  displays this help message
   -d  A directory where to store data. Default working directory
   -s  user specified strain(s) or a file containing strains (one per line)

Example: perl dbSNPQuery.pl [-d resultPath] [-s strainsFileName]

HELP
exit;
}
my $dbname ='cgdsnpdb';my $host ='cgd';my $user ='pup';my $pwd  ='puppass';
my $dbh = DBI->connect("DBI:mysql:$dbname:$host",$user,$pwd)
                            or die "Can't connect to database! $!\n";
my %strainmap;my $source_id= 21; #Imputed - Jeremy R. Wang et al. 2012 , default imputed source
# get chromosomes list
my $getChr ="select chromosome_id, chromosome_name from snp_chromosome where chromosome_id<22 ";
#now get the strains listing for this source
my $getStrains="select synonym_name,s.strain_id from snp_strain_synonym s, snp_strain_by_source sc ";
   $getStrains .=" where sc.source_id=? and sc.strain_id=s.strain_id order by synonym_name ";
# get the SNPs list for the imputed source
my $snplistQuery  ="select m.snpid,chromosome_id,bp_position ";
   $snplistQuery .=" from snp_position m, snp_by_source s";
   $snplistQuery .=" where chromosome_id=? and s.source_id=? and m.snpid=s.snpid ";

my $snpdataQuery  ="select distinct ref_allele,snp_allele,is_CpG_site,mutation_type";
   $snpdataQuery .=" from snp_main where snpid=? ";
#now get the list of accessions for this SNP  
my $getAccession =" select accession_id from snp_accession where snpid=? order by snpid_id_type limit 1";

#Now get the Imputed genotype of this snp
my $getGenotype= " select i.strain_id,genotype_allele from snp_imputed i, snp_strain s";
   $getGenotype.=" where i.snpid=? and source_id=? and i.strain_id=s.strain_id order by strain_name";

my $snpchr   = $dbh->prepare($getChr)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpstr   = $dbh->prepare($getStrains)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snplist  = $dbh->prepare($snplistQuery)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpdata  = $dbh->prepare($snpdataQuery)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpac    = $dbh->prepare($getAccession)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpgenoImp  = $dbh->prepare($getGenotype)or die "Couldn't prepare statement: " . $dbh->errstr;
#Now get the Imputed genotype of this snp
my $getGenotype= " select i.strain_id,genotype_allele from snp_genotype i, snp_strain s";
   $getGenotype.=" where i.snpid=? and source_id=? and i.strain_id=s.strain_id order by strain_name";
my $snpgeno  = $dbh->prepare($getGenotype)or die "Couldn't prepare statement: " . $dbh->errstr;
$year=2012; %source_name=(20=>"DiversityArrayYang2011",21=>"ImputedJeremyR.Wang2012");
$source_id=20; $source=$source_name{$source_id};
#
#get the list of specified strains into a map
#
%userStrains=();
if($opt_s){
  if(-f $opt_s){open(IN,"$opt_s");while(<IN>){chomp($_);$userStrains{"$_"}=0;}}
  else{@strains=split(",",$opt_s);
       foreach $str(@strains){chomp($str);$userStrains{"$str"}=0;}
   }
}
#print "The source name is $source\n"; exit;
my @row=(); my @strainIndex=();%strmap=();$straindata="";%indexMap=();
#index strain positions
#three dictionaries:
# strain_name/strain_id, strain_id/index, strain_id/strain_name
$i=0; $snpstr->execute($source_id) or die "Can't execute query: ".$dbh->errstr . "\n";
while(($strain_name,$strain_id)=$snpstr->fetchrow_array()){
     if(exists($userStrains{"$strain_name"})){ 
        $indexMap{"$strain_id"}=$i;$straindata.=",$strain_name";++$i;
        push(@strainIndex,"$strain_id:$strain_name");
        $userStrains{"$strain_name"}=$strain_id;$strmap{"$strain_id"}=$strain_name;
      }
}
print "Total user strains from this $source : ".keys(%indexMap)."\n";
# get the list of all the chromosomes
$snpchr->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
if($snpchr->rows > 0){ #process each chromosome 
   ## here I need to filter user's strains
   while(($chromosome_id,$chromosome_name) = $snpchr->fetchrow_array()){
       #next if($chromosome_id!=20);
       $filename=($opt_d)?"$opt_d/chr$chromosome_name"."_$source-snps.txt":"chr$chromosome_name"."_$source-snps.txt";
       print "Processing  $chromosome_name\n"; 
       # now for each item category, open the data get kenroy listing
       open(OUT ,">$filename") or die "can't open $filename\n";
       print OUT "SNPID,MutationType,Alleles,Build37BpPos,IsCpGSite,MinorAllele,MinorAlleleFreq,StrainWithGenoCount";
       print OUT "$straindata\n";  #get the snp list
       $snplist->execute($chromosome_id,$source_id) or die "Can't execute query: " . $dbh->errstr . "\n";
       while(($snpid,$chromosome_id,$bp_position)=$snplist->fetchrow_array){ 
           #get the snp data
           $snpdata->execute($snpid);next if($snpdata->rows<=0);
           ($black6_allele,$snp_allele,$is_CpG_site,$mutation_type)=$snpdata->fetchrow_array();
           $snpac->execute($snpid);($accession_id)=$snpac->fetchrow_array();
           $black6_allele=uc($black6_allele);$snp_allele=uc($snp_allele);
           $mutation=($mutation_type==2)?"Transversion":"Transition";$IsCpGSite=($is_CpG_site==1)?"Yes":"No";
           $snpLine="$accession_id,$mutation,$black6_allele/$snp_allele,$bp_position,$IsCpGSite";$snpgenoh="";
           #get the genotype of this snp
           if($source_id==21||$source_id==16){ $snpgenoImp->execute($snpid,$source_id);$snpgenoh=$snpgenoImp;}
           else{$snpgeno->execute($snpid,$source_id);$snpgenoh=$snpgeno;}
           next if($snpgenoh->rows<=0);
           $genodata=""; my %allelemap=(); $allelemap{$black6_allele}=0;$allelemap{$snp_allele}=0;%genolist=();
           while(($strain_id,$geno_allele)=$snpgenoh->fetchrow_array()){
                 $geno_allele=uc($geno_allele);
                if(exists($indexMap{"$strain_id"})){
                 $index=$indexMap{"$strain_id"};$genolist{"$index"}{"$geno_allele"}=1;
                 $allelemap{$geno_allele}+=1;}
           } #compute minor allele
           @genos=();
           foreach $index(sort keys(%genolist)){
               $genoalleles="";
               while(($thisallele,$alcount)=each(%{$genolist{"$index"}})){
                      $genoalleles=($genoalleles eq "")?$thisallele:":$thisallele";}
                $genos[$index]=$genoalleles;
           }
           $genodata=join(",",@genos);$strcount=0;
           $b6_count=$allelemap{$black6_allele};$snp_count=$allelemap{$snp_allele}; 
           while(($geno_al,$count)=each(%allelemap)){$strcount+=$count if($geno_al=~/A|C|T|G/i);}
           next if($strcount<=0);
           $maf=($snp_count>$b6_count)? sprintf("%.1f",($b6_count/$strcount)*100):sprintf("%.1f",($snp_count/$strcount)*100);
           $minorallele=($snp_count>$b6_count)?$black6_allele:$snp_allele; $snpLine .=",$minorallele,$maf,$strcount,$genodata";        
           print OUT "$snpLine\n";   
       }
       close(OUT); #now compress this file;
       $zipfile="$filename"."."."tar"."."."gz";
       $command="tar -cvzf $zipfile $filename"; system($command); 
       if(-f $filename){system("rm $filename");}
   }
 }
 print "Program complete\n";
 exit(0);
           
 

