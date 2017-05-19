#!/usr/bin/perl

######################################################################
## This script load the search terms into the search engine
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
q{
7 | Ensembl gene id          | cgd_genes_desc           | accession_id  |
8 | Ensembl transc id        | cgd_transcripts_desc     | accession_id  |
9 | MGI gene id              | cgd_genes_ensembl_mgi    | mgi_geneid    |
10 | MGI gene symbol          | cgd_genes_ensembl_mgi    | marker_symbol |
17 | Entrez gene id           | cgd_genes_ensembl_entrez | entrez_geneid |
18 | Entrez gene symbol       | cgd_genes_ensembl_entrez | marker_symbol |

19 | Ensembl Protein          | cgd_ensembl_protein      | protein_id    |
11 | SNP accession id         | snp_accession            | accession_id  |
};
use DBI;
use Time::localtime;
$local_db="cgd_snpdb";                #default lookingglass snp database name
$localhost="cgd-dev.jax.org";      #deafult db server
$luser="lnh";$lpwd="lucie98";        

#local db connection
my $dbh = DBI->connect("DBI:mysql:database=$local_db;host=$localhost;mysql_local_infile=1",$luser,$lpwd);
if(!$dbh){print  "Could not connect to the local database $local_db:$!\n"; exit;}
$query="create temporary table cgd_genes_desc_temp(accession_id varchar(30),found tinyint default 0)";
$qh_create_genedesc_temp=$dbh->prepare($query);
$query="insert into cgd_genes_desc_temp(accession_id) select distinct accession_id from cgd_genes_desc";
$qh_insert_genedesc_temp=$dbh->prepare($query);
$query="update cgd_genes_desc_temp t,cgd_items i set t.found=1 where t.accession_id=i.item and items_type_id =7";
$qh_update_genedesc_temp=$dbh->prepare($query);
$query="insert into cgd_items select accession_id,7 from cgd_genes_desc_temp where found=0";
$qh_load_ens_gene=$dbh->prepare($query);


$query="create temporary table cgd_transcripts_desc_temp(accession_id varchar(30),found tinyint default 0)";
$qh_create_transcripts_desc_temp=$dbh->prepare($query);
$query="insert into cgd_transcripts_desc_temp(accession_id) select distinct accession_id from cgd_transcripts_desc";
$qh_insert_transcripts_desc_temp=$dbh->prepare($query);
$query="update cgd_transcripts_desc_temp t,cgd_items i set t.found=1 where t.accession_id=i.item and items_type_id =8";
$qh_update_transcripts_desc_temp=$dbh->prepare($query);
$query="insert into cgd_items select accession_id,8 from cgd_transcripts_desc_temp where found=0";
$qh_load_ens_tx=$dbh->prepare($query);

$query="create temporary table cgd_genes_ensembl_mgi_temp(";
$query.="marker_symbol varchar(30),mgi_geneid varchar(30),found tinyint default 0)";
$qh_create_genes_ensembl_mgi_temp=$dbh->prepare($query);
$query="insert into cgd_genes_ensembl_mgi_temp(marker_symbol,mgi_geneid) ";
$query.=" select distinct marker_symbol,mgi_geneid from cgd_genes_ensembl_mgi";
$qh_insert_genes_ensembl_mgi_temp=$dbh->prepare($query);
$query="update cgd_genes_ensembl_mgi_temp t,cgd_items i set t.found=1 where t.marker_symbol=i.item and items_type_id =10";
$qh_update_mgi_gene_temp=$dbh->prepare($query);
$query="insert into cgd_items select marker_symbol,10 from cgd_genes_ensembl_mgi_temp where found=0";
$qh_load_mgi_gene=$dbh->prepare($query);

$qh_update_mgifound=$dbh->prepare("update cgd_genes_ensembl_mgi_temp set found=0");
$query="update cgd_genes_ensembl_mgi_temp t,cgd_items i set t.found=1 where t.mgi_geneid=i.item and items_type_id =9";
$qh_update_mgi_id_temp=$dbh->prepare($query);
$query="insert into cgd_items select mgi_geneid,9 from cgd_genes_ensembl_mgi_temp where found=0";
$qh_load_mgi_id=$dbh->prepare($query);

$query="create temporary table cgd_genes_ensembl_entrez_temp(";
$query.="marker_symbol varchar(30),entrez_geneid varchar(30),found tinyint default 0)";
$qh_create_genes_ensembl_entrez_temp=$dbh->prepare($query);
$query="insert into cgd_genes_ensembl_entrez_temp(marker_symbol,entrez_geneid) ";
$query.=" select distinct marker_symbol,entrez_geneid from cgd_genes_ensembl_entrez";
$qh_insert_genes_ensembl_entrez_temp=$dbh->prepare($query);
$query="update cgd_genes_ensembl_entrez_temp t,cgd_items i set t.found=1 where t.marker_symbol=i.item and items_type_id =18";
$qh_update_genes_ensembl_entrez_temp=$dbh->prepare($query);
$query="insert into cgd_items select marker_symbol,18 from cgd_genes_ensembl_entrez_temp where found=0";
$qh_load_entr_gene=$dbh->prepare($query);


$qh_update_entrezfound=$dbh->prepare("update cgd_genes_ensembl_entrez_temp set found=0");
$query="update cgd_genes_ensembl_entrez_temp t,cgd_items i set t.found=1 where t.entrez_geneid=i.item and items_type_id =17";
$qh_update_entrez_id_temp=$dbh->prepare($query);
$query="insert into cgd_items select entrez_geneid,17 from cgd_genes_ensembl_entrez_temp where found=0";
$qh_load_entr_id=$dbh->prepare($query);

$query="create temporary table cgd_ensembl_protein_temp(protein_id varchar(30),found tinyint default 0)";
$qh_create_ensembl_protein_temp=$dbh->prepare($query);
$query="insert into cgd_ensembl_protein_temp(protein_id) select distinct protein_id from cgd_ensembl_protein";
$qh_insert_ensembl_protein_temp=$dbh->prepare($query);
$query="update cgd_ensembl_protein_temp t,cgd_items i set t.found=1 where t.protein_id=i.item and items_type_id =19";
$qh_update_ensembl_protein_temp=$dbh->prepare($query);
$query="insert into cgd_items select protein_id,19 from cgd_ensembl_protein_temp where found=0";
$qh_load_prot_id=$dbh->prepare($query);

$query="create temporary table snp_accession_temp(accession_id varchar(30),found tinyint default 0)";
$qh_create_snp_accession_temp=$dbh->prepare($query);
$query="insert into snp_accession_temp(accession_id) select distinct accession_id from snp_accession";
$qh_insert_snp_accession_temp=$dbh->prepare($query);
$query="update snp_accession_temp t,cgd_items i set t.found=1 where t.accession_id=i.item and items_type_id =11";
$qh_update_snp_accession_temp=$dbh->prepare($query);
$query="insert into cgd_items select accession_id,11 from snp_accession_temp where found=0";
$qh_load_snp_id=$dbh->prepare($query);

#load genes
open(LOG,">engine_load.log");
print "Loading CGD search engine\n";
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Starting load process :  $mon/$mday/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
$qh_create_genedesc_temp->execute();
$qh_insert_genedesc_temp->execute();$qh_update_genedesc_temp->execute();
$qh_load_ens_gene->execute();
print LOG "Ensembl Genes loaded\n";
$qh_create_transcripts_desc_temp->execute();$qh_insert_transcripts_desc_temp->execute();
$qh_update_transcripts_desc_temp->execute();
$qh_load_ens_tx->execute();
print LOG "Ensembl Transcripts loaded\n";
$qh_create_genes_ensembl_mgi_temp->execute();$qh_insert_genes_ensembl_mgi_temp->execute();
$qh_update_mgi_gene_temp->execute();
$qh_load_mgi_gene->execute();
print LOG "MGI Genes loaded\n";
$qh_update_mgifound->execute();$qh_update_mgi_id_temp->execute();
$qh_load_mgi_id->execute();
print LOG "MGI Gene IDs loaded\n";
$qh_create_genes_ensembl_entrez_temp->execute();
$qh_insert_genes_ensembl_entrez_temp->execute();$qh_update_genes_ensembl_entrez_temp->execute();
$qh_load_entr_gene->execute();
print LOG "Entrez Genes loaded\n";
$qh_update_entrezfound->execute();$qh_update_entrez_id_temp->execute();
$qh_load_entr_id->execute();
print LOG "Entrez Gene IDs loaded\n";
$qh_create_ensembl_protein_temp->execute();
$qh_insert_ensembl_protein_temp->execute();$qh_update_ensembl_protein_temp->execute();
$qh_load_prot_id->execute();
print LOG "Ensembl Proteins loaded\n";
$qh_create_snp_accession_temp->execute();$qh_insert_snp_accession_temp->execute();
$qh_update_snp_accession_temp->execute();
$qh_load_snp_id->execute();
print LOG "SNP accession ids loaded\n";
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Program Ends:  $mon/$mday/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
close(LOG); 
print "Program complete\n";



