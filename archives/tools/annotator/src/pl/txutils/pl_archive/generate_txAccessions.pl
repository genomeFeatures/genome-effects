#!/usr/bin/perl

#############################################################################
# This script porcesses novel transcripts data to
# 1. assign/generate rsd numbres 
# 2. generate a config file to be use later in the pipeline 
#    format: group,organism,version,annotation,genePrediction file

# 3. generate genePrediction file to load in db
# 4. generate rnaseq_transcript file to load in db
# 5. Load generated files into db
#
# The process:
#  get the current maximum of rsd_number from the database
#  Every novel transcript is given a rsd accession id (rsd_number)
#  the accession id of novel transcripts have the following format:rsd_[number]
#  where number is a positive integer value generated when a novel transcript is detected.
#  The following steps are taken when processing novel transcripts:
#  1. Get the database current max rsd_[number] (the number part of the accession id) into a global variable(maxNumber)
#  2. Create a data structure (dbrsdmap)to store database current novel transcripts 
#     and associated rsd numbers --> this is a dictionary
#     the assumption is that memory size is not an object. If this becomes a concern, we can
#     have this check done on the fly. In order word we can adjust the granularity of the map key.
#     --> this is empty if there is no novel transcripts in db
#   
#  3. Foreach strain:
#    i. Index the new rnaSeq file into data structures to facilitate detection of identical transcripts   
#    ii. foreach novel transcript (txNew) from generated data structure, do:
#       a. assign the rsd number
#          a.1 if(Exists(dbrsdmap{txNew})) -> thisAccession=dbrsdmap{txNew}
#          a.2 else{++maxNumber;  thisAccession=_rsd_[maxNumber];dbrsdmap{txNew}=thisAccession
#       b. generate organismversion-strain-sampleID-thisAccession record 
#       c. generate genePrediction record 
#  
##############################################################################
use DBI;use vars qw ($opt_h $opt_d); use Getopt::Std;use Time::localtime; use POSIX; ## set default variables
use Cwd;
$dbname="graber_transcriptdb"; $host="harlequin"; $user="lnh";$pass="lucie98";
$dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
getopts('hd:');
if(($opt_h)|| !($opt_d)){
    print <<HELP;
   print  "Run the script with the -h option to display usage\n"; 
HELP
exit;
}
my $dir = getcwd; chomp($dir);
print "Current working directory is $dir\n";
#get the current max of rsd_number
$query="select MAX(cast(substring(transcript_name from locate('_rsd_',transcript_name)+5) as unsigned)) ";
$query.=" from transcript_by_annotation WHERE transcript_name like '_rsd_%'";
$qh_getCurrentMaxRsdNumber= $dbh->prepare($query);

$query="select transcript_id from transcript t, chromosome c where organism_version_id=? and c.chromosome_name=? 
        and c.chromosome_id=t.chromosome_id and t.strand=? and t.tx_start=? and t.tx_end=? ";
$qh_getCordHits= $dbh->prepare($query);

my $getTxExons="select distinct gene_prediction_id,transcript_name, exon_start,exon_end from transcript_exon t,exon e  ";
   $getTxExons.=" where  transcript_id=? and t.organism_version_id=? and t.exon_id=e.exon_id ";
   $getTxExons.=" order by gene_prediction_id,transcript_name, exon_start asc";
my $qh_getTxExons = $dbh->prepare($getTxExons);
$qh_getPredId=$dbh->prepare("select gene_prediction_id from gene_prediction where gene_prediction_name in (?)");
my $qh_getOrgV = $dbh->prepare("select organism_version_id from organism_version where ucsc_db=?");
my $qh_getSampleid=$dbh->prepare("select sample_id from rnaseq_sample where sample_name=?");
my $qh_insertSampleid=$dbh->prepare("insert into rnaseq_sample(sample_name,sample_desc) values(?,?)");
#GET NOVEL transcripts accession ids
$query="select distinct transcript_id,exon_start,exon_end ";
$query.=" from transcript_exon t, exon e where transcript_name =? and t.organism_version_id=? ";
$query.=" and t.organism_version_id=e.organism_version_id ";
$query.=" and t.exon_id=e.exon_id order by exon_start ";
$qh_getexons=$dbh->prepare($query);
$query="select distinct transcript_name from transcript_by_annotation
        where organism_version_id=? and transcript_name like 'rsd_%' ";
$qh_getTxRsds=$dbh->prepare($query);
$query="select tx_start,tx_end,chromosome_name, strand from transcript t, chromosome c where transcript_id=? ";
$query.=" and organism_version_id=? and c.chromosome_id=t.chromosome_id ";
$qh_getTxCord=$dbh->prepare($query);

#load rnaseq_transcript table
$qh_deleteNovel=$dbh->prepare("delete from rnaseq_transcript where organism_version_id=? and gene_prediction_id=?");
$qh_loadRnaSeqTx=$dbh->prepare("load data local infile ? into table rnaseq_transcript ignore 1 lines");
$query="select count(*) from rnaseq_transcript where organism_version_id=? and gene_prediction_id=?";
$qh_getRowCount=$dbh->prepare($query);
## Binary search implementation. on a sorted sequencial array, getting the min and max are O(1) constant 
## I can use this for sorted array, sorted hash keys (cpp map)
sub binarySearch{ ($array_ref,$token,$index)=@_;
 $min=0;$max=@{$array_ref}-1; $value=-1; ##implementing binary search
 while($min<=$max){  $mid=$min+ceil(($max-$min)/2);
  if($keys[$min]==$token){$value=@{$array_ref}[$min];$$index=$min; last;}
  elsif($keys[$max]==$token){$value=@{$array_ref}[$max];$$index=$max;last;}
  elsif($keys[$mid]==$token){$value=@{$array_ref}[$mid];$$index=$mid;last;}
  else{ if($keys[$mid]<$token){$min=$mid+1;}else{$max=$mid-1;}}#compute the new range to search
 } return $value;
}

#######################################################################
# loadCurrentRsd: index current novel transcripts into a data structure
# in memory to facilitate tx-rsd mapping
sub loadCurrentRsd{
  my($rsd_hash_ref,$organism_version_id)=@_; $qh_getTxRsds->execute($organism_version_id);
  if($qh_getTxRsds->rows>0){
    while(($tx_name)=$qh_getTxRsds->fetchrow_array()){
       $qh_getexons->execute($tx_name,$organism_version_id);$exonStarts=""; $exonEnds="";$tx_id=0;
       while(($transcript_id,$exon_start,$exon_end)=$qh_getexons->fetchrow_array()){
            if($exonStarts eq ""){$exonStarts="$exon_start"; $exonEnds="$exon_end";$tx_id=$transcript_id;}
            else{$exonStarts.=",$exon_start"; $exonEnds.=",$exon_end";}
       }
       if($tx_id>0){
         $qh_getTxCord->execute($tx_id,$organism_version_id);
         ($txstart,$txend,$chrom,$strand)=$qh_getTxCord->fetchrow_array();
          $key="$chrom-$strand-$txstart-$txend";
          ${$rsd_hash_ref}{"$organism_version_id"}{"$key"}{"$exonStarts:$exonEnds"}=$tx_name;
       }
    }
  }
}$qh_getCurrentMaxRsdNumber->execute();my $currentMaxRsdNumber=0;
if($qh_getCurrentMaxRsdNumber->rows){$currentMaxRsdNumber=$qh_getCurrentMaxRsdNumber->fetchrow_array();}
$currentMaxRsdNumber=$currentMaxRsdNumber<=0?0:$currentMaxRsdNumber;$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
$token=50000;$index=-1;
#$value=binarySearch(\@keys,$token,\$index); print "Found $token - $value at position $index\n ";
@files=`ls $opt_d`; open(LOG,">novel_transcript.log");
print LOG "Max number=$currentMaxRsdNumber Program Started $mon $mday, $yday @ $hour:$min:$sec\n";
#filename format: organismVersion.annotationSource.smapleId.(.+)?.txt  example: mm9.mm9-WSBGene.WSB_11304.txt
%novelMap=();$local_counter=$currentMaxRsdNumber; $tx_accession_prefix="rsd_"; %strainFiles=();
#index files first by strains then by samples
foreach $file(@files){chomp($file); 
   if($file=~/(.+)\.(.+)\.(.+)(\..+)?\.txt/){$strainFiles{"$2"}{"$3"}=$file;$orgversion=$1;}
 }
 while(($prediction,$samplefiles)=each(%strainFiles)){
    open(MAIN,">tempfile-$prediction.txt");$orgversion="";
    print MAIN "name\tchrom\tstrand\ttxStart\ttxEnd\texonCount\texonStarts\texonEnds\tname2\n";
    open(db,">tempfile-$prediction-rnaseq_transcript.txt");
    print db "organism_version_id\tgene_prediction_id\ttranscript_name\tsample_id\n";
   print "Processing $prediction \n";  
    while(($sample,$file)=each(%{$samplefiles})){
     if($file=~/(.+)\.(.+)\.(.+)(\..+)?\.txt/){
      #get/generate database local ids for organism version,annotation source, and sample name 
      $sample=$3;$qh_getSampleid->execute($sample);$orgversion=$1;$desc=$4;
      if($qh_getSampleid->rows<=0){$qh_insertSampleid->execute($sample,$desc);$qh_getSampleid->execute($sample);}
      ($sample_id)=$qh_getSampleid->fetchrow_array();$prediction_id=0;$orgv_id=0; 
      $qh_getOrgV->execute($1);($orgv_id)=$qh_getOrgV->fetchrow_array(); $prediction=$2; $prediction=~s/B6Gene/ensGene/;
      #check if you need to load the rsd for this organism
      if($currentMaxRsdNumber>0){if(!exists($novelMap{"$orgv_id"})){
       print "Loading existing rsd ids\n";loadCurrentRsd(\%novelMap,$orgv_id);}}
      $qh_getPredId->execute($prediction);($prediction_id)=$qh_getPredId->fetchrow_array();
      open(IN,"$opt_d/$file");
      if(IN){     
         $header=<IN>;chomp($header);$line="";@fields=split("\t",$header);
         if(@fields){$index=0;
            while(@fields>0){$field=shift(@fields);chomp($field);
               if($line eq ""){$line="$field>$index";}else{$line.=",$field>$index";} ++$index;}
         }# get fields index from header line
         $chromindex=-1;$strandindex=-1;$txStartindex=-1;$txEndindex=-1;$geneIndex=-1;
         $exonCountindex=-1;$exonStartsindex=-1;$exonEndsindex=-1;
         if($line=~/,?chrom>(\d+)/i){$chromindex=$1;}if($line=~/,strand>(\d+)/i){$strandindex=$1;}
         if($line=~/,txStart>(\d+)/i){$txStartindex=$1;}if($line=~/,txEnd>(\d+)/i){$txEndindex=$1;}
         if($line=~/,exonCount>(\d+)/i){$exonCountindex=$1;}
         if($line=~/,exonStarts>(\d+)/i){$exonStartsindex=$1;}if($line=~/,exonEnds>(\d+)/i){$exonEndsindex=$1;}
         if($line=~/,name2>(\d+)/i){$geneIndex=$1;}@rowIndexmap=(); %txStartmap=();
         #create indexes- generate unique transcript regions;
         #Indexes are created to detect identical transcripts accross samples/strains --
         while(<IN>){ chomp($_);
           split /\t/;$txstart=$_[$txStartindex];$gene=$_[$geneIndex];
           $rowIndexmap[$.]=join ":",($_[$exonStartsindex],$_[$exonEndsindex]);
           if(!exists($txStartmap{"$txstart"})){
              $txStartmap{"$txstart"}=join ":",($_[$txEndindex],$_[$exonCountindex],$_[$chromindex],$_[$strandindex],$.,$gene);}
           else{
            $txStartmap{"$txstart"}.=";".join ":",($_[$txEndindex],$_[$exonCountindex],$_[$chromindex],$_[$strandindex],$.,$gene);}
         }#end of index creation
         $txCoordHit=0;$total=0;$overlap_total=0;
        while(($txstart,$line)=each(%txStartmap)){ split /;/,$line;
            while(@_){ $newline=shift(@_); ++$total; ++$overlap_total if($line=~/;/);
               ($txend,$excount,$chrom,$strand,$linenumber,$gene)=split(":", $newline);
                $chrom=~s/chr//i;$chrom=~s/ch//i;$chrom=uc($chrom);$chrom="M" if($chrom eq "MT");
                $qh_getCordHits->execute($orgv_id,$chrom,$strand,$txstart,$txend);
                if($qh_getCordHits->rows>0){$txCoordHit+=1;}
                if(!$currentMaxRsdNumber){ #first time to create _rsd_ ids no need to check exons
                    $exons=$rowIndexmap[$linenumber];$key="$chrom-$strand-$txstart-$txend";$accession="";
                    #identical transcripts across strain/smaple share the same access
                    if(!exists($novelMap{"$orgv_id"}{"$key"}{"$exons"})){#create a new accession id
                      ++$local_counter;$accession="$tx_accession_prefix$local_counter";
                      $novelMap{"$orgv_id"}{"$key"}{"$exons"}="$accession";
                    }else{$accession=$novelMap{"$orgv_id"}{"$key"}{"$exons"};}
                   print db "$orgv_id\t$prediction_id\t$accession\t$sample_id\n";($exstarts,$exends)=split(":",$exons);
                   print MAIN "$accession\t$chrom\t$strand\t$txstart\t$txend\t$excount\t$exstarts\t$exends\t$gene\n";
                }
            }
     
        }#end of while(($txstart,$line)=each(%txStartmap)){  
     } # end of  if(IN){  
    }# end of if ($file=~exp)
   }#  while(($sample,$file)=each(%{$samplefiles}))
   close(db); close(MAIN);
   #load rnaseq_transcript and the annotation 
   $rnaseqTxfile="$dir/tempfile-$prediction-rnaseq_transcript.txt";$annotFile="$dir/tempfile-$prediction.txt";
   $linecont = `wc -l $rnaseqTxfile`;chomp($linecont); #load ensembl exon ids
   if($linecont=~/^(\d+)\s$rnaseqTxfile/){$linecont=$1;
      if($linecont>0){--$linecont;
         $qh_deleteNovel->execute($orgv_id,$prediction_id);
         $qh_loadRnaSeqTx->execute($rnaseqTxfile);
         $qh_getRowCount->execute($orgv_id,$prediction_id);($rowCount)=$qh_getRowCount->fetchrow_array();
         print LOG " $rnaseqTxfile\t$linecont\t$rowCount\n";   
       } 
   }
  #load annotations
  $cmd="perl pl/graber_txdb/load_annotations.pl -d $dir -f $annotFile -a $prediction -v $orgversion -s 0 -e 1";
  system($cmd);system("rm tempfile-*");         
 }
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "Program ended Total= $total -- $mon $mday, $yday @ $hour:$min:$sec\n";



