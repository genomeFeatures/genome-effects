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
#steps:
q{
1. load interval into db intervalid, chr,start,end,line
2. get list of chromosomes from interval table
3. For each chrom, get list of SNPs then proceed
};
#####################################################################
use DBI;
use vars qw ($opt_h $opt_d $opt_f);
use Getopt::Std;
use POSIX;
use Time::localtime;
$mbp=1000000;
getopts('hd:f:');
if($opt_h) {
    print <<HELP;

 This script generates SNPs file by chromosome for the specified SNP source.
 The chromosome file names have the format: chrxx"."_SourceDate_snps.txt"
  for example (chr1_Imputed2012-snps.txt)

Usage:

   perl dbSNPQuery.pl [-d resultPath]

Arguments:
   -h  displays this help message
   -d  A directory where to store data. Default working directory

Example: perl dbSNPQuery.pl [-d resultPath]

HELP
exit;
}
my $dbname ='cgdsnpdb';my $host ='cgd';my $user ='lnh';my $pwd  ='sql4lucie';
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user,$pwd)
                            or die "Can't connect to database!\n";

$qh_trunc=$dbh->prepare("truncate table snps_temp");
$query="create temporary table if not exists snps_temp(snpid int unsigned primary key,";
$query.="chrom tinyint default 0,Pos int unsigned default 0,alleles varchar(4),";
$query.="ref_Count int default 0,other_Count int default 0,uniqStrid smallint default 0,strains varchar(25),";
$query.="index(chrom,Pos),index(Pos))";
$qh_create=$dbh->prepare($query);
$qh_load=$dbh->prepare("load data local infile ? into table snps_temp ignore 1 lines");
$get_region=$dbh->prepare("select* from snps_temp where chrom=? and Pos between ? and ?");


$query="create temporary table if not exists intervals(interval_id int unsigned primary key auto_increment,";
$query.="chrom_id tinyint default 0, start int unsigned default 0, end int unsigned default 0,line text,index(chrom_id),";
$query.="index(start),index(end))";
$qh_createIntervals=$dbh->prepare($query);
$load_intervals=$dbh->prepare("load data local infile ? into table intervals (chrom_id,start,end,line)");
$get_rowsCount=$dbh->prepare("select count(*) from intervals");

$qh_truncIntervSnp=$dbh->prepare("truncate table snp_interval");
$query="create temporary table if not exists snp_interval(interval_id int unsigned,snpid int unsigned, ";
$query.="ref_allele char(1),snp_allele char(1),index(interval_id),index(snpid))";
$create_intsnpmap=$dbh->prepare($query)or die "Couldn't prepare statement: " . $dbh->errstr;
$query="insert into snp_interval select interval_id,p.snpid,ref_allele,snp_allele from intervals i,snp_position p,snp_main m ";
$query.="where i.chrom_id=? and i.chrom_id=p.chromosome_id and p.bp_position between start and end and p.snpid=m.snpid";
$load_intervalSnps=$dbh->prepare($query)or die "Couldn't prepare statement: " . $dbh->errstr;
$get_snpcount=$dbh->prepare("select count(*) from snp_interval");


$get_intervL=$dbh->prepare("select line from intervals where interval_id=?");

my %strainmap;my $source_id= 21; #Imputed - Jeremy R. Wang et al. 2012 , default imputed source
# get chromosomes list
#now get cc founder strains listing 
my $getStrains="select strain_name,s.strain_id from snp_strain s, snp_strain_by_group sc ";
   $getStrains .=" where group_id=4 and sc.strain_id=s.strain_id order by strain_name ";
# get the SNPs list for the imputed source
my $snplistQuery  ="select m.snpid,chromosome_id,bp_position,interval_id ";
   $snplistQuery .=" from snp_position m,intervals i, snp_by_source s";
   $snplistQuery .=" where chromosome_id=? and chromosome_id=chrom_id ";
   $snplistQuery .=" and bp_position between start and end"; 
   $snplistQuery .=" and s.source_id=? and m.snpid=s.snpid ";


#Now get the Imputed genotype of this snp
my $getGenotype= " select i.strain_id,genotype_allele from snp_imputed i, snp_strain s";
   $getGenotype.=" where snpid=? and source_id=? and i.strain_id in(1,2,7,8,16,129,213,214) 
                   and i.strain_id=s.strain_id order by strain_name";
$query="select distinct chrom_id,chromosome_name from intervals i, snp_chromosome c ";
$query.=" where i.chrom_id=c.chromosome_id";
my $snpchr   = $dbh->prepare($query);
my $snpstr   = $dbh->prepare($getStrains)or die "Couldn't prepare statement: " . $dbh->errstr;

my $snplist  = $dbh->prepare($snplistQuery)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpdata  = $dbh->prepare($snpdataQuery)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpac    = $dbh->prepare($getAccession)or die "Couldn't prepare statement: " . $dbh->errstr;
my $snpgenoImp  = $dbh->prepare($getGenotype)or die "Couldn't prepare statement: " . $dbh->errstr;
$year=2012; %source_name=(20=>"DiversityArrayYang2011",21=>"ImputedJeremyR.Wang2012");
$source_id=21; $source=$source_name{$source_id};
#print "The source name is $source\n"; exit;
# get the list of all the chromosomes
#load intervals
open(IN,"$opt_f");if(!IN){print "$!\n";exit;}open(LOG,">$opt_f-log.txt");
open(OUT,">$opt_f-intervals.txt");$count=0;
$region_headerline=<IN>;chomp($region_headerline);
  while(<IN>){chomp($_); @fields=split("\t",$_);
     if(@fields==15){ $line= join(",",@fields);++$count;
        $lodDropL=$fields[5]; $lodDropR=$fields[6];$lodDropL=~s/^\s*//;$lodDropL=~s/\s*$//;
        $chr=$fields[2];$chr=~s/^\s*//;$chr=~s/\s*$//;$lodDropR=~s/^\s*//;$lodDropR=~s/\s*$//;
        $chrl=$fields[8];$region_start=$fields[9]*$mbp;$chrr=$fields[12];$region_end=$fields[13]*$mbp;
        $chrl=~s/\s*$//;$chrl=~s/^\s*//;$chrr=~s/\s*$//;$chrr=~s/^\s*//;
        if(($chr ne "$chrl")||($chr ne "$chrr")){
            print LOG join("\t",@fields)."\t-- $chrl:$region_start\t$chrr:$region_end\n";next;}
          $chr_id =($chr=~/x/i)?20:$chr;$chr_id =($chr=~/y/i)?21:$chr;
         print OUT "$chr_id\t$region_start\t$region_end\t$line\n";
     } 
  }
close(OUT); print "Program complete\n";
$qh_createIntervals->execute();$load_intervals->execute("$opt_f-intervals.txt");
$get_rowsCount->execute();($rowcount)=$get_rowsCount->fetchrow_array();
print "Total from file:$count\n";print "Total from db:$rowcount\n";
if($count!=$rowcount){
  print "load failed\n";exit;
}
#index strain positions
$snpstr->execute() or die "Can't execute query: ".$dbh->errstr . "\n";
while(($strain_name,$strain_id)=$snpstr->fetchrow_array()){ push(@strainIndex,"$strain_id:$strain_name");}
$straindata="";%indexMap=();@strains=();
for my $i (0 .. $#strainIndex) {
      ($strain_id,$strain_name)=split(":",$strainIndex[$i]);$straindata.=",$strain_name";
       $indexMap{"$strain_id"}=$i;$strains[$i]=$strain_name;
     # print "$strain_id,$strain_name,$i\n";
 } 
#exit;
$straindata=~s/^,//; %chrmap=();%regionMap=();
my @strainIndex=();$snpchr->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
if($snpchr->rows > 0){ #process each chromosome 
  $count=0;
  while(($chromosome_id,$chromosome_name) =$snpchr->fetchrow_array()){ #
     # $create_intsnpmap->execute();$qh_truncIntervSnp->execute();
     # $load_intervalSnps->execute($chromosome_id);
     #$get_snpcount->execute();($rowCount)=$get_snpcount->fetchrow_array();
      $tm = localtime;
      my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
       ($tm->mon)+1, ($tm->year)+1900);
     print "\n*************************************************************\n";
     print "chrom $chromosome_name started:  $mday/$mon/$yday @ $hour:$min:$sec \n";
     print "\n*************************************************************\n";
     $snplist->execute($chromosome_id,$source_id) or die "Can't execute query: " . $dbh->errstr . "\n";
     $Totalsnps=$snplist->rows; 
     print "$chromosome_name SNPs count: $Totalsnps\n";
       $tm = localtime;
      my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
       ($tm->mon)+1, ($tm->year)+1900);
     print "\n*************************************************************\n";
     print "chrom $chromosome_name ended:  $mday/$mon/$yday @ $hour:$min:$sec \n";
     print "\n*************************************************************\n";
     last;
     $tm = localtime;
     open(SNP,">$opt_f-chr$chromosome_name-SNPdb.txt");
     print SNP "$region_headerline\tTotal\t".join("\t",@strains)."\n";
    
     $filename="chr$chromosome_id"."_$source-snpsdb.txt";
     print "Processing  $chromosome_name\n"; $count=0;%snpmap=();
     #get intervals for this chr
       open(OUT ,">$filename") or die "can't open $filename\n";
       print OUT "SNPID\tchrom\tPos\talleles\tref_alleleCount\tother_alleleCount\tuniqStrid\t";
       print OUT "\t$straindata\n";$count=0;  #get the snp list
       
       while(($interval_id,$snpid,$black6_allele,$snp_allele)=$snplist->fetchrow_array){ 
           #get the snp data
           $black6_allele=uc($black6_allele);$snp_allele=uc($snp_allele);
           $snpLine="$snpid\t$chromosome_id\t$bp_position\t$black6_allele/$snp_allele";$snpgenoh="";
           #get the genotype of this snp
           $snpgenoImp->execute($snpid,$source_id);next if($snpgenoImp->rows<=0);
           $strain_id=7;$index=$indexMap{"$strain_id"};
           %allele_freq=();#ref_allele stores the ref_allele and snp_allele frequencies for a given snp
           $allele_freq{"$black6_allele"}{"strains"}="$strain_id,"; $allele_freq{"$black6_allele"}{"count"}=1;
           $genodata=""; my %allelemap=(); $allelemap{$black6_allele}=0;$allelemap{$snp_allele}=0;@genolist=(0,0,0,0,0,0,0,0);
           $genolist[$index]="$black6_allele";
           while(($strain_id,$geno_allele)=$snpgenoImp->fetchrow_array()){
                  $index=$indexMap{"$strain_id"}; $geno_allele=uc($geno_allele);
                  $genolist[$index]=$geno_allele;
                  $allele_freq{"$geno_allele"}{"strains"}.="$strain_id,";$allele_freq{"$geno_allele"}{"count"}+=1;
           } #compute minor allele
           next if(keys(%allele_freq)<=1); #all cc founders have the same genotype
           next if(($allele_freq{"$black6_allele"}{"count"}+$allele_freq{"$snp_allele"}{"count"})!=8);
           ++$count;print "$count processed \n" if($count%500000==0);
           $genodata=join(",",@genolist);$uniqId="-";
           if($allele_freq{"$black6_allele"}{"count"}==1){$uniqId=$allele_freq{"$black6_allele"}{"strains"};}
           elsif($allele_freq{"$snp_allele"}{"count"}==1){$uniqId=$allele_freq{"$snp_allele"}{"strains"};}$uniqId=~s/,//;
           $snpLine .="\t".$allele_freq{"$black6_allele"}{"count"}."\t".$allele_freq{"$snp_allele"}{"count"};
           $snpLine.="\t$uniqId\t$genodata"; $ref_count=$allele_freq{"$black6_allele"}{"count"};
           $snp_count=$allele_freq{"$snp_allele"}{"count"};
           $snpmap{"$bp_position"}="$ref_count,$snp_count,$uniqId";      
           print OUT "$snpLine\n"; #last if($count%10000==0);
       }
       close(OUT);
       $tm = localtime;
      print "$chromosome_name- Total SNPs:$Totalsnps; TotalDifferent: ".keys(%snpmap)."\n";
      print "Total regions:".keys(%{$regionMap{"$chromosome_id"}})."\n";
      my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
       ($tm->mon)+1, ($tm->year)+1900);
      print "\n*************************************************************\n";
      print "parsing $chromosome_name from db ended:  $mday/$mon/$yday @ $hour:$min:$sec \n";
      print "\n*************************************************************\n";
      #now load the database
      $qh_create->execute();$qh_trunc->execute();$qh_load->execute($filename);
       #now for every interval on this chromosome, get the region stats
      $newCount=0;
      while(($linenumber,$line)=each(%{$regionMap{"$chromosome_id"}})){
          @fields=split("\t",$line);$region_start=$fields[9]*$mbp;$region_end=$fields[13]*$mbp;
          @genolist=(0,0,0,0,0,0,0,0);$total=0; $newCount+=1;
          $get_region->execute($chromosome_id,$region_start,$region_end);
          while(($snpid,$chrom,$Pos,$alleles,$ref_count,$snp_count,$uniqStrid,$strains)=
                 $get_region->fetchrow_array()){$total+=1; 
              if(($ref_count==1)||($snp_count==1)){$index=$indexMap{"$uniqStrid"};$genolist[$index]+=1;}
          }
          print "$newCount regions processed \n" if($newCount%100==0);
          if($total>0){print SNP "$line\t$total\t".join ("\t",@genolist)."\n";}
          #last if($newCount>100);
       }
      $tm = localtime;
      my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, 
     ($tm->mon)+1, ($tm->year)+1900);
     print "\n*************************************************************\n";
     print "parsing regions on $chromosome_name ended:  $mday/$mon/$yday @ $hour:$min:$sec \n";
     print "\n*************************************************************\n";
   }
 }
 print "Program complete\n";
 exit(0);
           
 

