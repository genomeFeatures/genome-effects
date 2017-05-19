#!/usr/bin/perl

######################################################################
## This script gets gene annotations from the current version
## of ensembl database . The organism is mus_musculus
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : April 2011
#   Usage: ./getEnsemblAnnotations.pl 
#
#  Output:
#
# Assumption: the asumption is that ensembl database namimg convention is
#              mus_musculus_core_&_&c
#  Where & is one or more digits and c is a character
# for example : mus_musculus_core_62_37o
#
# Tables to load:
q{
| cgd_genes                | 
| cgd_genes_desc           | 
| cgd_genes_ensembl_mgi    | 
| cgd_genes_ensembl_entrez | 
| cgd_ensembl_protein      | 

| cgd_externaldb           | 
| cgd_transcript_exons     | 
| cgd_transcripts          | 
| cgd_transcripts_desc     | 
  cgd_exons_desc
};
#######################################################################
## set default variables
$organism_prefix="mus_musculus_core_"; # default organism
$data_dir=`pwd`;                        #default data directory
$local_db="cgd_snpdb";                #default lookingglass snp database name
$localhost="cgd-dev.jax.org";      #deafult db server
$luser="lnh";$lpwd="lucie98";        
$sql_scripts_path="/home/lnh/ensembl_load/cgd_snpdb_sql";

$organism_file="/data/annotations/ensembl/organism_list.txt"; #list of organisms
$organism_path="/data/annotations/ensembl/";

$host="ensembldb.ensembl.org";
$port="5306";
$user="anonymous";

use DBI;
use vars qw ( $opt_h $opt_o $opt_i $opt_d $opt_l);
use Getopt::Std;

sub getLatestVersion(){
 $db_comand="mysql --user=$user --port=$port -h $host -e";
 #get the most rescent versionm of mus_musculus_core
 $comand= "$db_comand 'show databases;'"; @dblist=`$comand`;
 #now get of the latest version of mus_musculus_core_
 @mus_musculus= grep(/$organism_prefix/,@dblist);
 $eschema=0; $erelease=0;$db_version="";
 %releasemap;
 foreach my $version(@mus_musculus){ #sort the releases into a map
    chomp($version);
    if($version=~/$organism_prefix(\d+)_(\d+)/){
       $releasemap{$1}{$2}="$version";
    }
 }
 # now get the most current release and the corresponding version
 foreach $ensembl_schema(keys %releasemap){
          $eschema=$ensembl_schema if($ensembl_schema>$eschema);
 }
 ## now get the latest release of this version
 foreach $ensembl_release(keys %{$releasemap{$eschema}}){
          $erelease=$ensembl_release if($ensembl_release>$erelease);
 }
 $eschema=$releasemap{$eschema}{$erelease};

 return $eschema;
}
getopts('hi:o:d:l:');
if($opt_h) {
    print <<HELP;

 This script gets gene annotations from the current version
 of ensembl database . The organism is mus_musculus

Usage:

   perl getEnsemblAnnotations.pl [-i <ensembl_db_prefix (mus_musculus_core_)> -d <path to data directory>

Arguments:
   -h  displays this help message
   -i  ensembl_db_prefix for example mus_musculus_core_ (optional)
   -o  path to data directory (optional)
   -d  cgd SNP database name (cgd_snpdb, or cgdsnpdb) (optional)
   -l  cgd snp database server (cgd-dev.jax.org default) (optional)

HELP
exit;
}
####################################### Queries ##########################################
 ##Get data to load cgd_genes table (all ensembl genes)****
   $get_genes="select distinct g.gene_id, (case when s.name='X' then 20 
                                 when s.name='Y' then 21
                                 when s.name='MT' then 22 else s.name end)as chromosome_id,
     g.seq_region_strand,g.seq_region_start, g.seq_region_end
     from gene g,seq_region s ,coord_system c where g.seq_region_id=s.seq_region_id
     and s.coord_system_id=c.coord_system_id and c.rank=1";
## Get data to load cgd_genes_desc table***
  $get_gene_desc="select distinct g.gene_id,g.stable_id as ensembl_accession_id,g.biotype as gene_type,10
     from gene g,seq_region s,coord_system c
     where g.seq_region_id=s.seq_region_id
     and s.coord_system_id=c.coord_system_id and c.rank=1";

$seq_regionQuery="select s.name,s.seq_region_id from seq_region s, coord_system c where s.coord_system_id=c.coord_system_id and c.rank=1";

## Get data to load into cgd_ensembl_mgi *****
$get_mgi_gene=" select distinct x.dbprimary_acc as mgi_id,x.display_label as marker_symbol,
                  ge.description as marker_name, ge.gene_id
                 from  xref x,gene ge, external_db e,seq_region s,coord_system c 
                 where x.external_db_id=e.external_db_id
                 and e.db_name='MGI' AND x.xref_id=ge.display_xref_id 
                 and ge.seq_region_id= s.seq_region_id
                 and s.coord_system_id=c.coord_system_id and c.rank=1";

## Get data to load into cgd_ensembl_entrez *****
$get_entrez_gene='select distinct x.dbprimary_acc as entrez_id,x.display_label as marker_symbol, 
                 x.description as marker_name,ge.gene_id 
                 from  external_db e,xref x,object_xref o,gene ge,seq_region s,coord_system c 
                 where e.db_name="EntrezGene" and e.external_db_id=x.external_db_id
                 and  x.xref_id=o.xref_id and o.ensembl_object_type="Gene" and o.ensembl_id=ge.gene_id 
                 and ge.seq_region_id= s.seq_region_id and s.coord_system_id=c.coord_system_id and c.rank=1';

#get ensembl variations
$get_variations="select (case when s.name='X' then 20 when s.name='Y' then 21 when s.name='MT' then 22 else s.name end)as 
           chromosome_id,seq_region_start as bp_position,(case when seq_region_strand=1 then '+' else '-' end)as strand,
           variation_name as accession_id,tv.allele_string as alleles,tv.consequence_types as function_flass,tv.cds_start posInCds, 
          tv.translation_start as posInProtein,tv.codon_allele_string as codons,tv.pep_allele_string as aminoacids
          from transcript_variation tv, variation_feature v, seq_region s  where  tv.variation_feature_id=v.variation_feature_id
          and tv.allele_string=v.allele_string and tv.consequence_types=v.consequence_type and v.seq_region_id=s.seq_region_id ";
#cds_start >0 and feature_stable_id="ENSMUST00000169862" and
##*** load ensembl protein table ***

$get_protein="select distinct t.transcript_id,tl.stable_id as protein_id,
       seq_start,start_exon_id,seq_end,end_exon_id,
       if((t.seq_region_strand= 1),e.seq_region_start+seq_start-1,
          e.seq_region_end-seq_start+1)as cds_start,
      if((t.seq_region_strand=1),e2.seq_region_start+seq_end-1,e2.seq_region_end-seq_end-1)as cds_stop
    from seq_region s,transcript t,translation tl,exon e,exon e2
    where s.coord_system_id=1 and   s.seq_region_id = t.seq_region_id
    and   t.transcript_id=tl.transcript_id 
    and   tl.start_exon_id = e.exon_id and tl.end_exon_id = e2.exon_id";

## Get data to load into cgd_transcript table *
$get_transcript=" select transcript_id as transcript_local_id,gene_id,(case when s.name='X' then 20 
                                 when s.name='Y' then 21 when s.name='MT' then 22 else s.name end)as chromosome_id,
                    seq_region_start as transcript_start,seq_region_end as transcript_end,
                    (case when seq_region_strand=1 then '+' else '-' end) as strand
                    from transcript g,seq_region s ,coord_system c
                    where s.coord_system_id=c.coord_system_id and c.rank=1 and  g.seq_region_id=s.seq_region_id";


## get data to load into transcript_desc *
$get_transcript_desc="select transcript_id,stable_id from transcript g,seq_region s ,coord_system c
                    where g.seq_region_id=s.seq_region_id
                    and s.coord_system_id=c.coord_system_id and c.rank=1";

## load ensembl transcript exons table ****
$get_transcript_exon="select distinct t.transcript_id,e.exon_id,(case when s.name='X' then 20 
                                 when s.name='Y' then 21
                                 when s.name='MT' then 22 else s.name end)as chromosome_id,
              e.seq_region_start as exon_tart,e.seq_region_end as exon_end,
             t.rank as exon_rank, e.phase as start_phase,e.end_phase
            from exon_transcript t,exon e,seq_region s,coord_system c
            where s.coord_system_id=c.coord_system_id and c.rank=1 and s.seq_region_id = e.seq_region_id
            and   e.exon_id= t.exon_id";

## get data to load into exon_desc *
$get_exon_desc="select exon_id,stable_id from exon g,seq_region s,coord_system c
                    where g.seq_region_id=s.seq_region_id
                     and s.coord_system_id=c.coord_system_id and c.rank=1";

## generate tables for snp annotator -- exon.table.txt
$get_exon_table="select et.transcript_id,et.exon_id,stable_id,e.seq_region_id,e.seq_region_start,";
$get_exon_table.="e.seq_region_end,phase,end_phase,et.rank,tc.seq_region_strand";
$get_exon_table.=" from exon_transcript et,exon e,seq_region s,coord_system c ,transcript tc";
$get_exon_table.=" where tc.transcript_id=et.transcript_id and et.exon_id= e.exon_id ";
$get_exon_table.=" and e.seq_region_id= s.seq_region_id";
$get_exon_table.=" and s.coord_system_id=c.coord_system_id and c.rank=1";

## load exon_transcript.table.txt file
$get_exon_tx_table="select et.transcript_id,et.rank,et.exon_id from exon_transcript et,transcript t ,seq_region s,coord_system c";
$get_exon_tx_table.=" where et.transcript_id= t.transcript_id ";
$get_exon_tx_table.=" and t.seq_region_id= s.seq_region_id";
$get_exon_tx_table.=" and s.coord_system_id=c.coord_system_id and c.rank=1";
$get_exon_tx_table.=" order by et.transcript_id,et.rank  ";

## load gene table
$get_gene_table="select g.gene_id,g.seq_region_id,g.seq_region_start,";
$get_gene_table.= " g.seq_region_end,g.stable_id";
$get_gene_table.=" from gene g,seq_region s,coord_system c";
$get_gene_table.=" where g.seq_region_id=s.seq_region_id";
$get_gene_table.=" and s.coord_system_id=c.coord_system_id and c.rank=1";


## load translation table
$get_translation_table="select tl.transcript_id,tl.seq_start,tl.start_exon_id,tl.seq_end,tl.end_exon_id ";
$get_translation_table.=" from translation tl,transcript t,seq_region s ,coord_system c";
$get_translation_table.=" where tl.transcript_id= t.transcript_id";
$get_translation_table.=" and  t.seq_region_id=s.seq_region_id ";
$get_translation_table.=" and s.coord_system_id=c.coord_system_id and c.rank=1";
## load transcript table
$get_tx_table="select t.transcript_id,gene_id,stable_id, t.seq_region_id,seq_region_start,";
$get_tx_table.=" t.seq_region_end,t.seq_region_strand,t.biotype,s.name";
$get_tx_table.=" from transcript t,seq_region s,coord_system c";
$get_tx_table.=" where t.seq_region_id= s.seq_region_id";
$get_tx_table.=" and s.coord_system_id=c.coord_system_id and c.rank=1";

$db_version=getLatestVersion();
#local db connection
my $dbh2 = DBI->connect("DBI:mysql:database=$local_db;host=$localhost;mysql_local_infile=1",$luser,$lpwd);
if(!$dbh2){
    print  "Could not connect to the local database $local_db:$!\n"; exit;
}
#ensembl db connection
my $dbh = DBI->connect("DBI:mysql:$db_version:$host:$port",$user);
if(!$dbh){
    print  "Could not connect to Ensembl database :$!\n"; exit;
} 

 my $qh_get_exon_table= $dbh->prepare($get_exon_table);
   my $qh_get_exon_tx_table= $dbh->prepare($get_exon_tx_table);
   my $qh_get_gene_table= $dbh->prepare($get_gene_table);
   my $qh_get_translation_table= $dbh->prepare($get_translation_table);
   my $qh_get_tx_table= $dbh->prepare($get_tx_table);
  
   my $qh_get_genes= $dbh->prepare($get_genes);
   my $qh_get_gene_desc= $dbh->prepare($get_gene_desc);
   my $qh_get_mgi_gene= $dbh->prepare($get_mgi_gene);
   my $qh_get_entrez_gene   = $dbh->prepare($get_entrez_gene);
   my $qh_get_transcript = $dbh->prepare($get_transcript);
   my $qh_get_transcript_desc= $dbh->prepare($get_transcript_desc);
   my $qh_get_exon_desc= $dbh->prepare($get_exon_desc);
   my $qh_get_transcript_exon= $dbh->prepare($get_transcript_exon);
   my $qh_get_protein= $dbh->prepare($get_protein);
   
   $delete_genes=$dbh2->prepare("delete from cgd_genes");
   $delete_genes_desc=$dbh2->prepare("delete from cgd_genes_desc");
   $delete_genes_ens_mgi=$dbh2->prepare("delete from cgd_genes_ensembl_mgi");
   $delete_genes_ens_entr=$dbh2->prepare("delete from cgd_genes_ensembl_entrez");
   $delete_prot=$dbh2->prepare("delete from cgd_ensembl_protein");
   $delete_tx=$dbh2->prepare("delete from cgd_transcripts");
   $delete_tx_ex=$dbh2->prepare("delete from cgd_transcript_exons");
   $delete_tx_desc=$dbh2->prepare("delete from cgd_transcripts_desc");
   $delete_ex_desc=$dbh2->prepare("delete from cgd_exons_desc");
   
   $load_genes=$dbh2->prepare("load data local infile ? into table cgd_genes ignore 1 lines");
   $load_genes_desc=$dbh2->prepare("load data local infile ? into table cgd_genes_desc  ignore 1 lines");
   $load_genes_ens_mgi=$dbh2->prepare("load data local infile ? into table cgd_genes_ensembl_mgi  ignore 1 lines");
   $load_genes_ens_entr=$dbh2->prepare("load data local infile ? into table cgd_genes_ensembl_entrez  ignore 1 lines");
   $load_prot=$dbh2->prepare("load data local infile ? into table cgd_ensembl_protein ignore 1 lines");
   $load_tx=$dbh2->prepare("load data local infile ? into table cgd_transcripts ignore 1 lines");
   $load_tx_ex=$dbh2->prepare("load data local infile ? into table cgd_transcript_exons ignore 1 lines");
   $load_tx_desc=$dbh2->prepare("load data local infile ? into table cgd_transcripts_desc ignore 1 lines");
   $load_ex_desc=$dbh2->prepare("load data local infile ? into table cgd_exons_desc ignore 1 lines");

   $qh_getGeneCount= $dbh2->prepare(" select count(*) from cgd_genes");
   $qh_getGeneDescCount=$dbh2->prepare(" select count(*) from cgd_genes_desc");
   $qh_getGeneMgiCount=$dbh2->prepare(" select count(*) from cgd_genes_ensembl_mgi");
   $qh_getGeneEntrezCount=$dbh2->prepare(" select count(*) from cgd_genes_ensembl_entrez");
   $qh_getGeneProteinCount=$dbh2->prepare(" select count(*) from cgd_ensembl_protein");
   $qh_getTxCount=$dbh2->prepare(" select count(*) from cgd_transcripts");
   $qh_getTxDescCount=$dbh2->prepare(" select count(*) from cgd_transcripts_desc");
   $qh_getTxExonCount=$dbh2->prepare(" select count(*) from cgd_transcript_exons");
   $qh_getExonDescCount=$dbh2->prepare(" select count(*) from cgd_exons_desc");
  $qry="DELETE from cgd_annot_load_test where db_table_name in('cgd_genes','cgd_genes_desc','cgd_genes_ensembl_mgi',";
  $qry.="'cgd_genes_ensembl_entrez','cgd_ensembl_protein','cgd_transcripts','cgd_transcripts_desc',";
  $qry.="'cgd_exons_desc','cgd_transcript_exons')";
  $qh_delete_log=$dbh2->prepare($qry);
  $qh_load_log=$dbh2->prepare("LOAD DATA LOCAL INFILE ? INTO TABLE cgd_annot_load_test");
  $qh_update_schema=$dbh2->prepare("UPDATE cgd_externaldb set schema_version=? where db_name='ENSEMBL'");
  $qh_update_date=$dbh2->prepare("UPDATE cgd_externaldb set update_time=concat(CURDATE(),':',CURTIME()) where schema_version=?");
########################################
if($opt_i){$organism_prefix = $opt_i;} 
if($opt_o){$data_dir=$opt_o;} 
if($opt_d){$local_db=$opt_d;}
if($opt_l){$localhost=$opt_l;}
# set up some flags
$bad_localdb_con=0;$bad_remotedb_con=0;
$qh_getEnsemblSchema=$dbh2->prepare("select schema_version from cgd_externaldb where db_name='ENSEMBL'");
$qh_getEnsemblSchema->execute();$old_version="";
if(($old_version)=$qh_getEnsemblSchema->fetchrow_array()){
   if($old_version eq $db_version){
      print "Update not needed, database annotations in synch with Ensembl\n";exit;
    }
    else{print "Update needed,local db uses $old_version vs $db_version for Ensembl\n";}
   # check if data directory is provided by the user
   $data_dir=~s/\s+$//;$data_dir=~s/\/$//;
   $organism_data_dir="$data_dir/$db_version";
   print "Processing $organism_data_dir,\n";
   $checksum="111111111"; #9 bits corresponding to the nine files of gene annotations
   $validate_download=""; # a string of 9 bits where a bit is set when the corresponding file is loaded
   if(!(-d "$organism_data_dir")){ mkdir("$organism_data_dir",0777);}
   $seq_file="$organism_data_dir/ensemblSeqregionList.txt";
   $command="$db_comand 'use $db_version;$seq_regionQuery'> $seq_file";
   if(!(-f $seq_file)){system($command);}
   open(LOG,">$organism_data_dir/ensemblAnnot_load.log");$data="";
   # get gene annotations
   if(!(-f "$organism_data_dir/genes_table.txt")){
     open(ANT,">$organism_data_dir/genes_table.txt");
     if(ANT){ print ANT "gene_id\tchromosome_id\tgene_strand\tgene_start\tgene_end\n";
         $qh_get_genes->execute();$item_count=0;
         while(($gene_id,$chromosome_id,$gene_strand,$gene_start,$gene_end)=$qh_get_genes->fetchrow_array()){
              print ANT "$gene_id\t$chromosome_id\t$gene_strand\t$gene_start\t$gene_end\n";
             ++$item_count;
          }
          $validate_download.="1" if($item_count>0);
     }else{$validate_download.="0";}close(ANT);
   }else{ #get line count of this file
        $command="wc -l $organism_data_dir/genes_table.txt";$item_count=`$command`;
        if($item_count=~/(\d+)/){$item_count=$1;}--$item_count;
        if($item_count>0){$validate_download.="1";}
        else{$validate_download.="0";}
   } 
  if(OUT){
      if((-f "$organism_data_dir/genes_table.txt")){
           $delete_genes->execute();$load_genes->execute("$organism_data_dir/genes_table.txt");
           $qh_getGeneCount->execute();
           if($qh_getGeneCount->rows>0){
              ($row_count)=$qh_getGeneCount->fetchrow_array();
               print LOG "$organism_data_dir/genes_table.txt\t$item_count\tcgd_genes\t$row_count\n";
           }
       }
   }print "Genes table downloaded\n";
   # now get gene accessions
  if(!(-f "$organism_data_dir/")){open(ANT,">$organism_data_dir/genes_desc_table.txt");
     if(ANT){print ANT "gene_id\taccession_id\tgene_type\tsource_id\n";
       $qh_get_gene_desc->execute();$item_count=0;
        while(($gene_id,$accession_id,$gene_type,$source_id)=$qh_get_gene_desc->fetchrow_array()){
             print ANT "$gene_id\t$accession_id\t$gene_type\t$source_id\n";
             ++$item_count;
        }
        $validate_download.="1" if($item_count>0);
     }else{$validate_download.="0";}close(ANT);
  }else{
    #get line count of this file
    $command="wc -l $organism_data_dir/genes_desc_table.txt";
    $item_count=`$command`;
    if($item_count=~/(\d+)/){$item_count=$1;}--$item_count;
    if($item_count>0){$validate_download.="1";}
    else{$validate_download.="0";}
  }
  if(OUT){
        if((-f "$organism_data_dir/genes_desc_table.txt")){
           $delete_genes_desc->execute();$load_genes_desc->execute("$organism_data_dir/genes_desc_table.txt");
            $qh_getGeneDescCount->execute();
           if($qh_getGeneDescCount->rows>0){
              ($row_count)=$qh_getGeneDescCount->fetchrow_array();
               print LOG "$organism_data_dir/genes_desc_table.txt\t$item_count\tcgd_genes_desc\t$row_count\n";
           }
        }
  }print "Genes desc table downloaded\n";
   # now get cgd_genes_ensembl_mgi
  if(!(-f "$organism_data_dir/genes_mgi_table.txt")){
      open(ANT,">$organism_data_dir/genes_mgi_table.txt");
     if(ANT){print ANT "mgi_geneid\tmarker_symbol\tmarker_name\tgene_id\n";
        $qh_get_mgi_gene->execute(); $item_count=0;
        while(($mgi_geneid,$marker_symbol,$marker_name,$gene_id)=$qh_get_mgi_gene->fetchrow_array()){
             print ANT "$mgi_geneid\t$marker_symbol\t$marker_name\t$gene_id\n";
              ++$item_count;
        }
       $validate_download.="1" if($item_count>0);
     }else{$validate_download.="0";}close(ANT);
  }else{$command="wc -l $organism_data_dir/genes_mgi_table.txt";
       $item_count=`$command`;
       if($item_count=~/(\d+)/){$item_count=$1;}--$item_count;
       if($item_count>0){$validate_download.="1";}
       else{$validate_download.="0";}
  }
  if(OUT){
    if((-f "$organism_data_dir/genes_mgi_table.txt")){
         $delete_genes_ens_mgi->execute();$load_genes_ens_mgi->execute("$organism_data_dir/genes_mgi_table.txt");
         $qh_getGeneMgiCount->execute();
         if($qh_getGeneMgiCount->rows>0){
            ($row_count)= $qh_getGeneMgiCount->fetchrow_array();
             print LOG "$organism_data_dir/genes_mgi_table.txt\t$item_count\tcgd_genes_ensembl_mgi\t$row_count\n";
         }
     }
   }print "MGI Genes table downloaded\n";
   if(!(-f "$organism_data_dir/genes_entrez_table.txt")){
      open(ANT,">$organism_data_dir/genes_entrez_table.txt");
     if(ANT){print ANT "entrez_geneid\tmarker_symbol\tmarker_name\tgene_id\n";
        $qh_get_entrez_gene->execute();$item_count=0;
        while(($mgi_geneid,$marker_symbol,$marker_name,$gene_id)=$qh_get_entrez_gene->fetchrow_array()){
             print ANT "$mgi_geneid\t$marker_symbol\t$marker_name\t$gene_id\n";
             ++$item_count;
        }
        $validate_download.="1" if($item_count>0);
     }else{$validate_download.="0";}close(ANT);
  }else{$command="wc -l $organism_data_dir/genes_entrez_table.txt";
       $item_count=`$command`;
       if($item_count=~/(\d+)/){$item_count=$1;}
       --$item_count;
       if($item_count>0){$validate_download.="1";}
       else{$validate_download.="0";}
  }
  if(OUT){
        if((-f "$organism_data_dir/genes_entrez_table.txt")){
           $delete_genes_ens_entr->execute();$load_genes_ens_entr->execute("$organism_data_dir/genes_entrez_table.txt");
           $qh_getGeneEntrezCount->execute();
           if($qh_getGeneEntrezCount->rows>0){
             ($row_count)=$qh_getGeneEntrezCount->fetchrow_array();
              print LOG "$organism_data_dir/genes_entrez_table.txt\t$item_count\tcgd_genes_ensembl_entrez\t$row_count\n";
           }
        }
   }print "ENTREZ Genes table downloaded\n";
  if(!(-f "$organism_data_dir/genes_protein_table.txt")){
     open(ANT,">$organism_data_dir/genes_protein_table.txt");
     if(ANT){print ANT "tx_id\tprotein\tseq_start\tstart_ex_id\tseq_end\tend_ex_id\tcdstart\tcdse\n";
        $qh_get_protein->execute();$item_count=0;
        while(($tx_id,$protein,$seq_start,$start_ex_id,$seq_end,$end_ex_id,$cdstart,$cdse)=$qh_get_protein->fetchrow_array()){
             print ANT "$tx_id\t$protein\t$seq_start\t$start_ex_id\t$seq_end\t$end_ex_id\t$cdstart\t$cdse\n";
             ++$item_count;
        }$validate_download.="1" if($item_count>0);
      }else{$validate_download.="0";}close(ANT);
  }else{$command="wc -l $organism_data_dir/genes_protein_table.txt";
       $item_count=`$command`;
       if($item_count=~/(\d+)/){$item_count=$1;}--$item_count;
       if($item_count>0){$validate_download.="1";}
       else{$validate_download.="0";}
  }
  if(OUT){
        if((-f "$organism_data_dir/genes_protein_table.txt")){
           $delete_prot->execute();$load_prot->execute("$organism_data_dir/genes_protein_table.txt");
           $qh_getGeneProteinCount->execute();
           if($qh_getGeneProteinCount->rows>0){
              ($row_count)=$qh_getGeneProteinCount->fetchrow_array();
               print LOG "$organism_data_dir/genes_protein_table.txt\t$item_count\tcgd_ensembl_protein\t$row_count\n";
           }
        }
   } print "Protein table downloaded\n";
  if(!(-f "$organism_data_dir/transcript_table.txt")){
     open(ANT,">$organism_data_dir/transcript_table.txt");
     if(ANT){print ANT "tx_id\tgene_id\tchr\ttx_start\ttx_end\tstrand\n";
        $qh_get_transcript->execute();$item_count=0;
        while(($tx_id,$gene_id,$chrom_id,$tx_start,$tx_end,$strand)=$qh_get_transcript->fetchrow_array()){
             print ANT "$tx_id\t$gene_id\t$chrom_id\t$tx_start\t$tx_end\t$strand\n";
             ++$item_count;
        }
        $validate_download.="1" if($item_count>0);
     }else{$validate_download.="0";}close(ANT);
  }else{$command="wc -l $organism_data_dir/transcript_table.txt";
       $item_count=`$command`;
       if($item_count=~/(\d+)/){$item_count=$1;}--$item_count;
       if($item_count>0){$validate_download.="1";}
       else{$validate_download.="0";}
  }
  if(OUT){
        if((-f "$organism_data_dir/transcript_table.txt")){
           $delete_tx->execute();$load_tx->execute("$organism_data_dir/transcript_table.txt");
           $qh_getTxCount->execute();
           if($qh_getTxCount->rows>0){
             ($row_count)=$qh_getTxCount->fetchrow_array();
             print LOG "$organism_data_dir/transcript_table.txt\t$item_count\tcgd_transcripts\t$row_count\n";
           }
        }
   }print "Transcripts table downloaded\n";
  if(!(-f "$organism_data_dir/transcript_desc_table.txt")){
     open(ANT,">$organism_data_dir/transcript_desc_table.txt");
     if(ANT){print ANT "tx_id\ttx_accession\n";
        $qh_get_transcript_desc->execute();$item_count=0;
        while(($tx_id,$tx)=$qh_get_transcript_desc->fetchrow_array()){
             print ANT "$tx_id\t$tx\n";++$item_count;
        }
        $validate_download.="1" if($item_count>0);
     }else{$validate_download.="0";}close(ANT);
  }else{$command="wc -l $organism_data_dir/transcript_desc_table.txt";
       $item_count=`$command`;
       if($item_count=~/(\d+)/){$item_count=$1;}--$item_count;
       if($item_count>0){$validate_download.="1";}
       else{$validate_download.="0";}
  }
 if(OUT){
       if((-f "$organism_data_dir/transcript_desc_table.txt")){
           $delete_tx_desc->execute();$load_tx_desc->execute("$organism_data_dir/transcript_desc_table.txt");
           $qh_getTxDescCount->execute();
           if($qh_getTxDescCount->rows>0){
              ($row_count)=$qh_getTxDescCount->fetchrow_array();
               print LOG "$organism_data_dir/transcript_desc_table.txt\t$item_count\tcgd_transcripts_desc\t$row_count\n";
           }
        }
   }print "Transcript desc table downloaded\n";
  if(!(-f "$organism_data_dir/exon_desc_table.txt")){
     open(ANT,">$organism_data_dir/exon_desc_table.txt");
     if(ANT){print ANT "exon_id\texon_accession\n";
        $qh_get_exon_desc->execute();$item_count=0;
        while(($ex_id,$ex_accession)=$qh_get_exon_desc->fetchrow_array()){
             print ANT "$ex_id\t$ex_accession\n";++$item_count;
        }
        $validate_download.="1" if($item_count>0);
     }else{$validate_download.="0";}close(ANT);
  }else{$command="wc -l $organism_data_dir/exon_desc_table.txt";
       $item_count=`$command`;
       if($item_count=~/(\d+)/){$item_count=$1;}--$item_count;
        if($item_count>0){$validate_download.="1";}
        else{$validate_download.="0";}
  }
 if(OUT){
        if((-f "$organism_data_dir/exon_desc_table.txt")){
          $delete_ex_desc->execute();$load_ex_desc->execute("$organism_data_dir/exon_desc_table.txt");
          $qh_getExonDescCount->execute();
          if($qh_getExonDescCount->rows>0){
             ($row_count)=$qh_getExonDescCount->fetchrow_array();
              print LOG "$organism_data_dir/exon_desc_table.txt\t$item_count\tcgd_exons_desc\t$row_count\n";
           }
        }
   }print "exon_desc_table table downloaded\n";
  if(!(-f "$organism_data_dir/transcript_exon_table.txt")){
     open(ANT,">$organism_data_dir/transcript_exon_table.txt");
     if(ANT){print ANT "tx_id\texon_id\tchr\texon_start\texon_end\texon_rank\tstart_phase\tend_phase\n";
        $qh_get_transcript_exon->execute();$item_count=0;
        while(($tx_id,$ex_id,$chr,$ex_start,$ex_end,$ex_rank,$start_phase,$end_phase)=$qh_get_transcript_exon->fetchrow_array()){
             print ANT "$tx_id\t$ex_id\t$chr\t$ex_start\t$ex_end\t$ex_rank\t$start_phase\t$end_phase\n";
             ++$item_count;
        }
        $validate_download.="1" if($item_count>0);
     }else{$validate_download.="0";}close(ANT);
  }else{
       #get line count of this file
       $command="wc -l $organism_data_dir/transcript_exon_table.txt";
       $item_count=`$command`;
       if($item_count=~/(\d+)/){$item_count=$1;}--$item_count;
       if($item_count>0){$validate_download.="1";}
       else{$validate_download.="0";}
  }
  if(OUT){
        if((-f "$organism_data_dir/transcript_exon_table.txt")){
            $delete_tx_ex->execute();$load_tx_ex->execute("$organism_data_dir/transcript_exon_table.txt");
            $qh_getTxExonCount->execute();
            if($qh_getTxExonCount->rows>0){
               ($row_count)= $qh_getTxExonCount->fetchrow_array();
                print LOG "$organism_data_dir/transcript_exon_table.txt\t$item_count\tcgd_transcript_exons\t$row_count\n";
             } 
        }
   }
  print "Transcript exons table downloaded\n";
  #################################### SNP annotator tables ###############################
  $text="$db_version\n";
  open(SH,">$organism_data_dir/ensembleSchema.txt") or die ("Bad file :$!\n");
  if(SH){print SH "$text"; close(SH);}
  close(LOG);
  # now update the database
 if($checksum eq $validate_download){ #the annotations were downloaded without issues
    print "Updating the database\n";
    $qh_delete_log->execute();
    $qh_load_log->execute("$organism_data_dir/ensemblAnnot_load.log");
    $qh_update_schema->execute($db_version);
    $qh_update_date->execute($db_version);
 }
 else{print "Database update fails\n";}
}
print "Program complete\n";
exit(0);
           
 

