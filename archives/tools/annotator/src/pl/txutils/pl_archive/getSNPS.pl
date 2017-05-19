#!/usr/bin/perl

q{
This script is to generate a test SNPs set for the annotator
We will include SNPs from all function classes : Intergenic, Intronic,UTR, and coding
 The sample will contain 100 cases of each
};
use DBI;
use Getopt::Std;
use vars qw($opt_h,$opt_f,$opt_v);
getopts("hf:v:");

if(!$opt_f){
print <<HELP;
   Usage: ./program -f filename 
   Where filename is the results file name
HELP
exit;
}

$user="pup"; $upass="puppass";$database="cgdsnpdb";$host="cgd";
$db=DBI->connect("DBI:mysql:database=$database;host=$host",$user,$upass);

$query="select p.snpid,chromosome_id,bp_position,ref_allele,snp_allele,mutation_type,is_CpG_site from snp_position p, snp_main m";
$query.=" where is_intergenic=1  and p.snpid=m.snpid limit 100";
$get_snpInterGenicList=$db->prepare($query);

$query="select p.snpid,chromosome_id,bp_position,ref_allele,snp_allele,mutation_type,is_CpG_site from snp_position p, snp_main m";
$query.=" where p.snpid=? and p.snpid=m.snpid ";
$get_snpData=$db->prepare($query);

$get_intronic=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=0 limit 100");
$get_3utrplus=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=3 and  strand=1 limit 50");
$get_3utrminus=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=3 and  strand!=1 limit 50");
$get_5utrplus=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=4 and  strand=1 limit 50");
$get_5utrminus=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=4 and  strand!=1 limit 50");

$get_codingsynplus=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=1 and  strand=1 limit 50");
$get_codingsynminus=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=1 and  strand!=1 limit 50");
$get_codingNonsynplus=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=2 and  strand=1 limit 50");
$get_codingNonsynminus=$db->prepare("select distinct snpid from snp_transcript where _loc_func_key=2 and  strand!=1 limit 50");
$query="select distinct accession_id,_frame_key,PosInCDS,PosInProtein,ref_aa,ref_codon,snp_aa,snp_codon ";
$query.=" from snp_aminoacid a,cgd_transcripts_desc t where a.snpid=? and a.transcript_local_id=t.transcript_local_id ";
$get_codon=$db->prepare($query);

$query="select distinct accession_id from snp_transcript a,cgd_transcripts_desc t";
$query.="  where a.snpid=? and a._loc_func_key=? and a.strand=? a.transcript_local_id=t.transcript_local_id ";
$get_tx=$db->prepare($query);

#get the intergenic, Intronic, and UTR SNPs
$get_snpInterGenicList->execute();open(OUT,">$opt_f");
print OUT "snpid\tchromosome\tposition\tref\talt\tmutation_type\tis_CpG_site\tclass\t";
print OUT "name,_frame_key,PosInCDS,PosInProtein,ref_aa,ref_codon,snp_aa,snp_codon\n";
if($get_snpInterGenicList->rows>0){
  while(($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpInterGenicList->fetchrow_array()){
      print OUT "$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\tIntergenic\t.\n";
  }
}
$get_intronic->execute();
if($get_intronic->rows>0){
  while(($snpid)=$get_intronic->fetchrow_array()){
      $get_snpData->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         print OUT "$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\tIntronic\t.\n";
      }
  }
}
$get_3utrplus->execute();
if($get_3utrplus->rows>0){
  while(($snpid)=$get_3utrplus->fetchrow_array()){
      $get_snpData->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         
         print OUT "$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\t3'UTR:+\t.\n";
      }
  }
}
$get_3utrminus->execute();
if($get_3utrminus->rows>0){
  while(($snpid)=$get_3utrminus->fetchrow_array()){
      $get_snpData->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         print OUT "$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\t3'UTR:-\t.\n";
      }
  }
}
$get_5utrplus->execute();
if($get_5utrplus->rows>0){
  while(($snpid)=$get_5utrplus->fetchrow_array()){
      $get_snpData->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         print OUT "$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\t5'UTR:+\t.\n";
      }
  }
}
$get_5utrminus->execute();
if($get_5utrminus->rows>0){
  while(($snpid)=$get_5utrminus->fetchrow_array()){
      $get_snpData->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         print OUT "$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\t5'UTR:-\t.\n";
      }  
  }
}
$get_codingsynplus->execute();
if($get_codingsynplus->rows>0){
  while(($snpid)=$get_codingsynplus->fetchrow_array()){
      $get_snpData->execute($snpid);$data=""; $get_codon->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         $data="$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\tCodingSyn:+";
      }
      while(($tx,$frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon)=$get_codon->fetchrow_array()){
        print OUT "$data\t$tx,$frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon\n";
      }
  }
}
$get_codingsynminus->execute();
if($get_codingsynminus->rows>0){
  while(($snpid)=$get_codingsynminus->fetchrow_array()){
      $get_snpData->execute($snpid);$data=""; $get_codon->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         $data="$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\tCodingSyn:-";
      }
     while(($tx,$frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon)=$get_codon->fetchrow_array()){
        print OUT "$data\t$tx,$frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon\n";
      }
  }
}

$get_codingNonsynplus->execute();
if($get_codingNonsynplus->rows>0){
  while(($snpid)=$get_codingNonsynplus->fetchrow_array()){
      $get_snpData->execute($snpid);$data=""; $get_codon->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         $data="$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\tCodingNonSyn:+";
      }
      while(($tx,$frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon)=$get_codon->fetchrow_array()){
        print OUT "$data\t$tx,$frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon\n";
      }
  }
}
$get_codingNonsynminus->execute();
if($get_codingNonsynminus->rows>0){
  while(($snpid)=$get_codingNonsynminus->fetchrow_array()){
      $get_snpData->execute($snpid);$data=""; $get_codon->execute($snpid);
      if($get_snpData->rows>0){
         ($snpid,$chromosome_id,$bp_position,$ref_allele,$snp_allele,$mutation_type,$is_CpG_site)=$get_snpData->fetchrow_array();
         $data="$snpid\t$chromosome_id\t$bp_position\t$ref_allele\t$snp_allele\t$mutation_type\t$is_CpG_site\tCodingNonSyn:-";
      }
     while(($tx,$frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon)=$get_codon->fetchrow_array()){
        print OUT "$data\t$tx,$frame,$PosInCDS,$PosInProtein,$ref_aa,$ref_codon,$snp_aa,$snp_codon\n";
      }
  }
}





