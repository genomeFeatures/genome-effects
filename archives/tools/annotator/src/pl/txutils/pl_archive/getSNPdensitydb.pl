#!/usr/bin/perl
use vars qw($opt_h $opt_f);
use DBI;
use Getopt::Std;
use Time::localtime;
getopts('hf:');

$mbp="1000000";
$window="250000";$kb="250000";
$database="cgdsnpdb";$host="cgd";$pwd="puppass";$usr="pup";
$db=DBI->connect("DBI:mysql:database=$database;host=$host",$usr,$pwd);
if(!$opt_f){
print <<HELP;
   Usage: ./program -f filename -c cordinates_file
   Where filename is a file containing gene intervals and cordinates_file contains
   interval cordinates
  Note: the ensembl gene sometimes is not on the same chromosome as the interval markers 
        Narayanan said it was expected.
HELP
exit;
}
################# database connection and queries setings 
$qh_truncateSNPtemp=$db->prepare("truncate temporary table if exists snp_allele_temp");
$query="create temporary table if not exists snp_allele_temp(snpid int unsigned primary key,";
$query.="ref_allele char(1),snp_allele char(1))";
$qh_createSNPtemp=$db->prepare($query);
$query="insert into snp_allele_temp select p.snpid, ref_allele,snp_allele from snp_position p, snp_main m";
$query.=" where chromosome_id=? and bp_position between ? and ? and p.snpid=m.snpid";
$insert_tempSNP=$db->prepare($query);
$query="select p.snpid, ref_allele,snp_allele from snp_position p, snp_main m";
$query.=" where chromosome_id=? and bp_position between ? and ? and p.snpid=m.snpid";
$get_snpList=$db->prepare($query);

$query="select t.strain_id,genotype_allele from snp_imputed t,";
$query.=" (select strain_id from snp_strain_by_group where group_id=4) as s ";
$query.=" where snpid=? and t.strain_id=s.strain_id and source_id=21 ";
$get_geno=$db->prepare($query);
$get_chrom=$db->prepare("select chromosome_id,chromosome_name from snp_chromosome");
$getb6_allele=$db->prepare("select ref_allele,snp_allele from snp_main where snpid=?");
$query="select strain_name, s.strain_id from snp_strain_by_group p, snp_strain s";
$query.=" where p.group_id=4 and p.strain_id=s.strain_id order by strain_name";
$get_ccFounders=$db->prepare($query); $get_ccFounders->execute();

##################### reset strain genotype counts 
sub resetCount {
 ($indexMap_ref,$rowMap_ref)=@_;
 while(($index,$pattern)=each(%{$indexMap_ref})){${$rowMap_ref}{"$pattern"}{"count"}=0;}
  ${$rowMap_ref}{"all"}{"count"}=0;
} 
sub getChunks{
 ($chr_id,$start,$thisend,$indexMap_ref,$strainsIdPattern_ref,$rowMap_ref)=@_;
 $get_snpList->execute($chr_id,$start,$thisend);
 if($get_snpList->rows>0){$i=0;
    while(($snpid,$ref_allele,$snp_allele)=$get_snpList->fetchrow_array()){++$i;
         @genos=(0,0,0,0,0,0,0,0);$get_geno->execute($snpid);
         $ref_allele=uc($ref_allele);$snp_allele=uc($snp_allele);
         %allele_freq=();#ref_allele stores the ref_allele and snp_allele frequencies for a given snp
         $allele_freq{"$ref_allele"}{"count"}=0; $allele_freq{"$ref_allele"}{"strains"}="";
         $allele_freq{"$snp_allele"}{"count"}=0;$allele_freq{"$snp_allele"}{"strains"}="";
         if($get_geno->rows>0){ #process this snp
           $strain_id=7;$pattern=${$strainsIdPattern_ref}{"$strain_id"};$index=${$rowMap_ref}{"$pattern"}{"index"};
           $allele_freq{"$ref_allele"}{"strains"}="$strain_id,"; $allele_freq{"$ref_allele"}{"count"}+=1;
           while(($strain_id,$genotype_allele)=$get_geno->fetchrow_array()){
             $genotype_allele=uc($genotype_allele);
             $allele_freq{"$genotype_allele"}{"strains"}.="$strain_id,"; 
             $allele_freq{"$genotype_allele"}{"count"}+=1;
           }#now set the pattern of this genotype and update counts only if all 8 cc have known geno
           if($allele_freq{"$ref_allele"}{"count"}>0 && $allele_freq{"$snp_allele"}{"count"}>0){
              ${$rowMap_ref}{"all"}{"count"}+=1;  #at least one strain is different, update total count
              if($allele_freq{"$ref_allele"}{"count"}==1){ #b6 is the only one different
                 $strain_id=$allele_freq{"$ref_allele"}{"strains"};$strain_id=~s/,//;
                 $pattern=${$strainsIdPattern_ref}{"$strain_id"}; 
                 ${$rowMap_ref}{"$pattern"}{"count"}+=1;#update count where this strain is the only one with varation
              }elsif($allele_freq{"$snp_allele"}{"count"}==1){#all other strains but this strain have same geno as b6 
                 $strain_id=$allele_freq{"$snp_allele"}{"strains"};$strain_id=~s/,//;
                 $pattern=${$strainsIdPattern_ref}{"$strain_id"}; 
                 ${$rowMap_ref}{"$pattern"}{"count"}+=1;
              }
           } #else all the strains carry the reference or the other allele
      } #end of if($get_geno->rows>0){
   }
  }#end of if($get_snpList->rows>0)
}
################################################################################## 
# parse region data
# for each SNP in the region, get the strain genotype alleles
# if all strain have same genotype, don't don't inlcude this SNP in the count
# if at least one of the cc founders is different, include this in the total count
# if only this cc is different, increase this cc count
#
#################################################################################
sub getRegionData {
 ($chr_id,$start,$end,$indexMap_ref,$strainsIdPattern_ref,$rowMap_ref)=@_;
 resetCount($indexMap_ref,$rowMap_ref);$range=$end-$start; ++$count;$more=1;
 while($more){
      $newEnd=(($start+$window)>$end)?$end:($start+$window);
      getChunks($chr_id,$start,$newEnd,$indexMap_ref,$strainsIdPattern_ref,$rowMap_ref);
      $start=$newEnd+1; $more=0 if($newEnd==$end);
  }
}
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
       ($tm->mon)+1, ($tm->year)+1900);
  print "\n*************************************************************\n";
  print "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print "\n*************************************************************\n";
@strainsIndex=();$i=0;%strainsIdPattern=(); #maps strain id to pattern
#create a structure that detects genotype patterns from the 8 CC founders genotype row 
# indexMap maps strain index to a pattern
%indexMap=("0"=>"10000000","1"=>"01000000","2"=>"00100000","3"=>"00010000",
           "4"=>"00001000","5"=>"00000100","6"=>"00000010","7"=>"00000001");
%rowMap=();%chrommap=();$get_chrom->execute();
if($get_chrom->rows>0){
  while(($chr_id,$chr_name)=$get_chrom->fetchrow_array()){
    $chr_name=~s/^\s*//;$chr_name=~s/\s*$//;$chrommap{"$chr_name"}=$chr_id;}
}
if($get_ccFounders->rows>0){
   while(($strain_name,$strain_id)=$get_ccFounders->fetchrow_array()){
        $strainsIndex[$i]=$strain_name;
        $pattern=$indexMap{"$i"};$rowMap{"$pattern"}{"strain_name"}=$strain_name;
        $rowMap{"$pattern"}{"strain_id"}=$strain_id;$rowMap{"$pattern"}{"index"}=$i;
        $strainsIdPattern{"$strain_id"}=$pattern;
        ++$i;
  }
}
open(OUT,">$opt_f-SNPs.txt");open(OUT2,">log-log.txt");
#get the header line
#ensID	marker	marker.chr	marker.loc	lod	lodDropL	lodDropR	lower_Name	lower_Chromosome	lower_Position	lower_Avg.cM	upper_Name	upper_Chromosome	upper_Position	upper_Avg.cM
open(IN,"$opt_f"); $header=<IN>; chomp($header);
print OUT "$header\tregionStart\tregionEnd\tTotal\t".join("\t",@strainsIndex)."\n";$count=0;
while(<IN>){chomp($_);
  @fields=split("\t",$_);
  if(@fields==15){ $line= join("\t",@fields);#++$count;
     $lodDropL=$fields[5]; $lodDropR=$fields[6];$lodDropL=~s/^\s*//;$lodDropL=~s/\s*$//;
     $chr=$fields[2];$chr=~s/^\s*//;$chr=~s/\s*$//;$lodDropR=~s/^\s*//;$lodDropR=~s/\s*$//;
     $chrl=$fields[8];$region_start=$fields[9]*$mbp;$chrr=$fields[12];$region_end=$fields[13]*$mbp;
     $chrl=~s/\s*$//;$chrl=~s/^\s*//;$chrr=~s/\s*$//;$chrr=~s/^\s*//;
     if(($chr ne "$chrl")||($chr ne "$chrr")){
        print OUT2 join("\t",@fields)."\t-- $chrl:$region_start\t$chrr:$region_end\n";next;}
     $chr_id=$chrommap{"$chr"};
     getRegionData($chr_id,$region_start,$region_end,\%indexMap,\%strainsIdPattern,\%rowMap);
    if($rowMap{"all"}{"count"}>0){
       print OUT "$line\t$region_start\t$region_end\t"; @genos=(0,0,0,0,0,0,0,0);
       while(($index,$pattern)=each(%indexMap)){$genos[$index]= $rowMap{"$pattern"}{"count"};}
       print OUT $rowMap{"all"}{"count"}."\t".join("\t", @genos)."\n"; 
       print "$count processed\n" if($count%100000==0);
    }
  }
 
}
 $tm = localtime;
  my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
       ($tm->mon)+1, ($tm->year)+1900);
  print "\n*************************************************************\n";
  print "Program Ends:  $mday/$mon/$yday @ $hour:$min:$sec \n";
  print "\n*************************************************************\n";
print "Program complete\n";


