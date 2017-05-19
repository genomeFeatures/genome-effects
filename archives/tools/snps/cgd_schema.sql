use cgd_snpdb_schema;

/****************************************************************************
 Author: Lucie Hutchins, Scientific Software Engineer
 Date :  May 2012
 Project: CGD SNP integration database
 
 Description: runs QC on a given SNP set and tags all the issues found,
              assign the gene implication, mutation type, CpG sites, and more

******************************************************************************
************************ *************/

drop table if exists snp_error_log;
drop table if exists cgd_annot_load_test;
drop table if exists cgd_mouse_codonFrequency;
drop table if exists cgd_externaldb;
drop table if exists cgd_blosum62NCBI;

drop table if exists cgd_ensembl_protein;
drop table if exists cgd_genes_ensembl_entrez;
drop table if exists cgd_genes_ensembl_mgi;
drop table if exists cgd_transcript_exons;
drop table if exists cgd_transcripts_desc;
drop table if exists cgd_genes_desc;
drop table if exists cgd_ensembl_variations;
drop table if exists snp_strain_synonym;
drop table if exists snp_strain_by_group;

drop table if exists cgd_items;

drop table if exists snp_accession;
drop table if exists snp_CpGSites;
drop table if exists snp_by_source;

drop table if exists snp_transcript;
drop table if exists snp_aminoacid;

drop table if exists snp_imputed;
drop table if exists snp_genotype;
drop table if exists snp_genotype_conflict;

drop table if exists snp_main; 
drop table if exists snp_multiAlleles;
drop table if exists snp_main_allele_conflict;

/**************** Primary key *************************/
drop table if exists cgd_aminoacid;
drop table if exists cgd_mutation_type;
drop table if exists snp_source;
drop table if exists snp_chromosome;
drop table if exists snp_frame;
drop table if exists snp_id_type;
drop table if exists snp_loc_func;
drop table if exists snp_strain_groups;
drop table if exists snp_strain;
drop table if exists snp_error_flag;
drop table if exists cgd_items_type;
drop table if exists cgd_transcripts;
drop table if exists cgd_exons_desc;
drop table if exists cgd_genes;
drop table if exists snp_position;
drop table if exists cgd_organism_codonCount;
/**** Tables creation starts here  ********************************************
primary key tables first
/******************************************************************************
 cgd_aminoacid:  stores aminoacid mapping between different annotations
 and the chemical characteristics of each amino acid
*******************************************************************************/
create table if not exists cgd_aminoacid(
  name  varchar(20) not NULL,
  three_letter_abrev char(3) not NULL,
  one_letter_abrev  char(1) PRIMARY KEY,
  acidity_polarity  varchar(25) ,
  hydropathy_index  float(2,1))ENGINE=InnoDB;
/******************************************************************************
cgd_mutation_type: Static table to store mutation type map (transition, transversion)
*******************************************************************************/
create table if not exists cgd_mutation_type(
 mutation_type tinyint primary key,
 mutation_name varchar(15)
)ENGINE=InnoDB;
/******************************************************************************
snp_source: Static table to store SNP provider info
*******************************************************************************/
create table if not exists snp_source(
 source_id smallint primary key,
 source_name  varchar(50),
 website   varchar(150),
 address  varchar(100),
 telephone  varchar(20),
 fax varchar(20),
 pi varchar(50),
 contact_person varchar(50)
)ENGINE=InnoDB;
/******************************************************************************
snp_chromosome: Static table to generate local id for every chromosome
*******************************************************************************/
create table if not exists snp_chromosome(
 chromosome_id tinyint not null,
 chromosome_name char(2),
 chr_size int,
 primary key(chromosome_id)
)ENGINE=InnoDB;
/******************************************************************************
snp_frame: Static table to generate local id for every possible
 SNP position within a codon [1,2,3]
*******************************************************************************/
create table if not exists snp_frame(
 _frame_key tinyint primary key,
 description varchar(30)
)ENGINE=InnoDB;
/******************************************************************************
snp_id_type: Static table to generate local id for every possible
 types of SNP accession rsid,ssid,submitterid,..)
*******************************************************************************/
create table if not exists snp_id_type(
 snpid_type_id tinyint not null,
 snpid_type varchar(20),
 primary key(snpid_type_id)
)ENGINE=InnoDB;
/******************************************************************************
snp_loc_func: Static table to generate local id for every possible
 types of SNP location (exon,intron, 3'utr,5'utr,...)
*******************************************************************************/
create table if not exists snp_loc_func(
 _loc_func_key tinyint primary key,
 description varchar(30)
)ENGINE=InnoDB;
/******************************************************************************
snp_strain_groups: Static table to generate a unique id for every strain group
(CC Founders, Wild, Wild Derived, Classical,..)
*******************************************************************************/
create table if not exists snp_strain_groups(
 group_id tinyint primary key auto_increment,
 group_name varchar(50)
)ENGINE=InnoDB;
/******************************************************************************
snp_strain: Static table to generate local id for strain_name
*******************************************************************************/
create table if not exists snp_strain(
 strain_id smallint primary key,
 strain_name varchar(40)
)ENGINE=InnoDB;
/******************************************************************************
snp_error_flag: asigns a local id to different types of error found during
 our QC step
*******************************************************************************/
create table if not exists snp_error_flag(
 error_code tinyint primary key,
 description varchar(100)
)ENGINE=InnoDB;
/******************************************************************************
cgd_items_type: Static table to store different type of search terms
     (MGI, Ensembl,....)
*******************************************************************************/
create table if not exists cgd_items_type(
 items_type_id tinyint  PRIMARY KEY,
 type_name varchar(25) NOT NULL,
 item_table_name varchar(50),
 item_key_name varchar(30)
)ENGINE=InnoDB;
/*********************************************************************************
 local ensembl genes

**********************************************************************************/
CREATE TABLE IF NOT EXISTS cgd_genes( 
    gene_id int unsigned NOT NULL,
    chromosome_id tinyint NOT NULL,
    gene_strand tinyint default 0,
    gene_start int unsigned default 0,
    gene_end int unsigned default 0,
    primary key(gene_id),
    index(chromosome_id),index(chromosome_id,gene_strand,gene_start),
    foreign key(chromosome_id) references snp_chromosome(chromosome_id)
    )ENGINE=InnoDB;
/*******************************************************************************
  gene transcript mapping foreign key(chromosome_id) references snp_chromosome(chromosome_id)  ON DELETE CASCADE,
********************************************************************************/
CREATE TABLE IF NOT EXISTS cgd_transcripts(
    transcript_local_id  int unsigned NOT NULL,
    gene_id int unsigned NOT NULL,
    chromosome_id tinyint NOT NULL,
    transcript_start int unsigned default 0,
    transcript_end int unsigned default 0,
    transcript_strand char(1),
    primary key(transcript_local_id),
    index(gene_id),index(chromosome_id),
    foreign key (gene_id) references cgd_genes (gene_id)   ON DELETE CASCADE,
    index(chromosome_id,transcript_strand,transcript_start)
)ENGINE=InnoDB;
/********************************************************************************
 exon accession ids
*********************************************************************************/
CREATE TABLE IF NOT EXISTS cgd_exons_desc(
  exon_id  int unsigned not null,
  accession_id varchar(30),
  primary key(exon_id),
  index(accession_id))ENGINE=InnoDB;
/***************************************************************************************
 snp_position: generates a unique id for a given chromosome-bp_position-ref_strain_id pair
****************************************************************************************/
create table if not exists snp_position(
 snpid int unsigned primary key auto_increment, 
 chromosome_id tinyint default 0,
 bp_position int unsigned default 0,
 ref_strain_id smallint default 0, 
 foreign key (chromosome_id) references snp_chromosome(chromosome_id),
 foreign key(ref_strain_id) references snp_strain(strain_id),
 index(chromosome_id),index(bp_position),index(chromosome_id,bp_position),
 index(ref_strain_id))ENGINE=InnoDB;
/******************************************************************************
cgd_organism_codonCount: Static table to store organism map
*******************************************************************************/
create table if not exists cgd_organism_codonCount(
 organism_id tinyint PRIMARY KEY auto_increment,
 organism_name varchar(40),
 cds_count  int(11) DEFAULT 0,
 codon_count int(10) unsigned DEFAULT 0
)ENGINE=InnoDB;

/************* end of  primary key tables ********************************************/
/***************************************************************************************
 snp_error_log: SNP to error code mapping
****************************************************************************************/
create table if not exists snp_error_log(
 snpid int unsigned default 0,
 source_id smallint default 0,
 error_flag tinyint default 0,
 primary key(snpid,source_id,error_flag),
 foreign key(snpid) references snp_position(snpid),
 foreign key(source_id) references snp_source(source_id),
 foreign key(error_flag) references snp_error_flag(error_code),
 index(snpid),index(source_id),index(error_flag))ENGINE=InnoDB;
/***************************************************************************************
 Ensembl annotations load testing
****************************************************************************************/
CREATE TABLE IF NOT EXISTS cgd_annot_load_test(
   annotation_file varchar(200), file_line_count int default 0,
   db_table_name varchar(50), 
   db_table_line_count int default 0)ENGINE=InnoDB;

/******************************************************************************
cgd_mouse_codonFrequency: Static table to store codon frequency
by organism
*******************************************************************************/
create table if not exists cgd_mouse_codonFrequency(
 organism_id   tinyint DEFAULT 0,
 codon   char(3) ,
 aa_one_letter_abrev char(1),
 codon_fraction_aa float(2,2),
 codon_fraction_1000  decimal(4,2),
 codon_count  int(11),
 foreign key(aa_one_letter_abrev) references cgd_aminoacid(one_letter_abrev),
 foreign key(organism_id) references cgd_organism_codonCount(organism_id),
 INDEX(codon),index(aa_one_letter_abrev),index(organism_id)
)ENGINE=InnoDB;
/******************************************************************************
cgd_externaldb: Static table to keep the load log
*******************************************************************************/
create table if not exists cgd_externaldb(
 source_id smallint not null default 0,
 db_name varchar(20),
 schema_version varchar(50),
 update_time datetime not null default "0000-00-00 00:00:00",
 foreign key(source_id) references snp_source(source_id)
)ENGINE=InnoDB;

/******************************************************************************
cgd_blosum62NCBI: Static table to store the blosum62NCBI matrix info
*******************************************************************************/
create table if not exists cgd_blosum62NCBI(
 amino_1 char(1) not null,
 amino_2 char(1) not null,
 score  tinyint default 0,
 foreign key(amino_1) references cgd_aminoacid(one_letter_abrev),
 index(amino_1),index(amino_2)
)ENGINE=InnoDB;
/************************* gene annotation supporting tables start here **********
This section will be replace with graber_transcriptdb
/********************************************************************************
 gene accession id mapping
*********************************************************************************/
CREATE TABLE IF NOT EXISTS cgd_genes_desc( 
    gene_id int unsigned,
    accession_id varchar(30),
    gene_type varchar(20),
    source_id smallint,
    foreign key(gene_id) references cgd_genes(gene_id),
    foreign key(source_id) references snp_source(source_id),
    index(gene_id),index(source_id),
    index(accession_id))ENGINE=InnoDB;

/*******************************************************************************
  now load ensembl_mgi mapping
********************************************************************************/
 CREATE TABLE IF NOT EXISTS cgd_genes_ensembl_mgi(
    mgi_geneid varchar(20),
    marker_symbol varchar(30),
    marker_name varchar(250),
    gene_id int unsigned NOT NULL,
    foreign key(gene_id) references cgd_genes(gene_id),
    index(gene_id),index(mgi_geneid), index(marker_symbol))ENGINE=InnoDB;

/*********************************************************************************
load ensembl entrez mapping 
**********************************************************************************/
 CREATE TABLE IF NOT EXISTS cgd_genes_ensembl_entrez(
    entrez_geneid varchar(20),
    marker_symbol varchar(30),
    marker_name varchar(250),
    gene_id int unsigned NOT NULL,
    foreign key(gene_id) references cgd_genes(gene_id),
    index(entrez_geneid),index(gene_id),
    index(marker_symbol))ENGINE=InnoDB;

/*****************************************************************************
 transcript exon mapping 
   foreign key(exon_id) references cgd_exon_desc(exon_id)
 
******************************************************************************/
CREATE TABLE IF NOT EXISTS cgd_transcript_exons(
    transcript_local_id  int unsigned not null,
    exon_id  int unsigned default 0,
    chromosome_id tinyint not null,
    exon_start   int unsigned default 0,
    exon_end   int unsigned default 0,
    exon_rank   tinyint default 0,
    start_phase   tinyint,
    end_phase tinyint,
    index(transcript_local_id),index(chromosome_id),
    index(exon_id),
    foreign key(transcript_local_id) references  cgd_transcripts(transcript_local_id),
    foreign key(chromosome_id) references snp_chromosome(chromosome_id),
    foreign key(exon_id) references cgd_exons_desc(exon_id)
    )ENGINE=InnoDB;
/*****************************************************************************
 transcript accession ids
******************************************************************************/
CREATE TABLE IF NOT EXISTS cgd_transcripts_desc(
    transcript_local_id  int unsigned default 0,
    accession_id varchar(30),
    foreign key(transcript_local_id) references  cgd_transcripts(transcript_local_id),
    index(transcript_local_id),
    index( accession_id))ENGINE=InnoDB;
/*********************************************************************************
 protein - transcript mapping
**********************************************************************************/
CREATE TABLE IF NOT EXISTS cgd_ensembl_protein(
    transcript_local_id  int unsigned NOT NULL,
    protein_id  varchar(30),
    seq_start     smallint default 0,
    start_exon_id  int unsigned NOT NULL,
    seq_end    smallint default 0,
    end_exon_id  int unsigned default 0,
    cds_start  int default 0,
    cds_end   int default 0,
    index(transcript_local_id),index(start_exon_id),index(end_exon_id),index(protein_id),
    foreign key(transcript_local_id) references  cgd_transcripts(transcript_local_id),
    foreign key(start_exon_id) references cgd_exons_desc(exon_id)
    )ENGINE=InnoDB;
/************************* gene annotation supporting tables end here **********
/*****************Strain supporting tables start
*****************************************************************************
snp_strain_synonym: Static table to map strain synonyms
*******************************************************************************/
create table if not exists snp_strain_synonym(
 strain_id smallint default 0,
 synonym_name varchar(200),
 primary key(strain_id,synonym_name),
 foreign key(strain_id) references snp_strain(strain_id),
 index(strain_id)
)ENGINE=InnoDB;
/******************************************************************************
snp_strain_by_group: Static table to map every strain to a strain group
*******************************************************************************/
create table if not exists snp_strain_by_group(
 strain_id smallint default 0,
 group_id tinyint default 0,
 primary key(strain_id,group_id),
 foreign key(strain_id) references snp_strain(strain_id),
 foreign key(group_id) references snp_strain_groups(group_id)
)ENGINE=InnoDB;
/*****************Strain supporting tables end
/**************************************************************************************
 local ensembl variations
***************************************************************************************
CREATE TABLE IF NOT EXISTS cgd_ensembl_variations( 
    chromosome_id tinyint default 0,bp_position int unsigned default 0,
    strand char(1),accession_id varchar(255),
    alleles varchar(20),function_flass varchar(200),posInCds int default 0,
    posInProtein int default 0,codons varchar(100),aminoacids varchar(50),
    tx_accession_id varchar(128),snp_frame tinyint default 0,
    foreign key(chromosome_id) references snp_chromosome(chromosome_id),
   index(chromosome_id), index(chromosome_id,bp_position),index(accession_id),
   index(tx_accession_id))ENGINE=InnoDB;

/*************************** Search Engine supporting  related tables start here ***********/
/***************************************************************************************
 cgd_items search engine
****************************************************************************************/
create table if not exists cgd_items(
 item varchar(50) not null,
 items_type_id tinyint default 0,
 primary key(item,items_type_id),
 foreign key(items_type_id) references cgd_items_type(items_type_id),
 index(item),index(items_type_id))ENGINE=InnoDB;
/**************************** SNP load supporting tables starts here *******************
/***************************************************************************************
 snp_accession: stores accession ids info for every SNP 
****************************************************************************************/
create table if not exists snp_accession(
  accession_id char(20) not null,
  snpid_id_type_id tinyint not null,
  snpid int unsigned default 0,
  source_id smallint default 0,
  primary key(accession_id,snpid,source_id),
  index(snpid),index(source_id),index(snpid_id_type_id),index(accession_id),
  foreign key(snpid) references snp_position(snpid),                            
  foreign key(source_id) references snp_source(source_id),
  foreign key(snpid_id_type_id) references snp_id_type(snpid_type_id)
)ENGINE=InnoDB;
/***************************************************************************************
 snp_by_source: stores SNP - source mapping
****************************************************************************************/
create table if not exists snp_by_source(
  snpid int unsigned default 0,
  source_id smallint default 0,
  primary key(snpid,source_id),
  foreign key(snpid) references snp_position(snpid),                            
  foreign key(source_id) references snp_source(source_id),
  index(snpid),index(source_id))ENGINE=InnoDB;
/***************************************************************************************
 snp_CpGSites: SNP to CpC Site mapping
****************************************************************************************/
create table if not exists snp_CpGSites(
 snpid int unsigned not null default 0,
 first_allele varchar(8) not null,second_allele varchar(8) not null,
 foreign key(snpid) references snp_position(snpid),
 index(snpid))ENGINE=InnoDB;
/***************************************************************************************
 snp_transcript: SNP to transcript mapping
****************************************************************************************/
create table if not exists snp_transcript(
 snpid int unsigned default 0,
 _loc_func_key tinyint default 0,
 strand tinyint default 0,
 transcript_local_id int unsigned default 0,
 gene_id int unsigned default 0,
 foreign key(snpid) references snp_position(snpid),
 foreign key(_loc_func_key) references snp_loc_func(_loc_func_key),
 foreign key(transcript_local_id) references cgd_transcripts(transcript_local_id),
 foreign key(gene_id) references cgd_genes(gene_id),
 index(snpid),index(transcript_local_id),index(_loc_func_key),
 index(strand),index(gene_id))ENGINE=InnoDB;
/***************************************************************************************
 snp_aminoacid: SNP to transcript to codon mapping
****************************************************************************************/
create table if not exists snp_aminoacid(
 snpid int unsigned default 0,
 transcript_local_id int unsigned default 0,
 _frame_key tinyint default 0,
 PosInCDS int default 0,
 PosInProtein int default 0,
 ref_aa char(3) not null,
 ref_codon  char(3) not null,
 snp_aa char(3) not null,
 snp_codon char(3) not null,
 foreign key(snpid) references snp_position(snpid),
 foreign key(transcript_local_id) references cgd_transcripts(transcript_local_id),
 foreign key(_frame_key) references snp_frame(_frame_key),
 index(snpid),index(transcript_local_id)
)ENGINE=InnoDB;
/***************************************************************************************
 snp_genotype: SNP genotype (excluding data from imputed sources)
****************************************************************************************/
create table if not exists snp_genotype(
 snpid int unsigned default 0,
 source_id smallint default 0,
 strain_id  smallint default 0,
 genotype_allele char(1) not null,
 primary key(snpid,source_id,strain_id,genotype_allele),
 foreign key(snpid) references snp_position(snpid),
 foreign key(source_id) references snp_source(source_id),
 foreign key(strain_id) references snp_strain(strain_id),
 index(snpid),index(source_id),
 index(strain_id))ENGINE=InnoDB;
/***************************************************************************************
 snp_imputed: SNP genotype with confidence score(only data from imputed sources)
****************************************************************************************/
create table if not exists snp_imputed(
 snpid int unsigned default 0,
 strain_id  smallint default 0,
 genotype_allele char(1) not null,
 confidence tinyint default 0,
 source_id smallint default 0,
 primary key(snpid,strain_id,genotype_allele,source_id),
 foreign key(snpid) references snp_position(snpid),
 foreign key(source_id) references snp_source(source_id),
 foreign key(strain_id) references snp_strain(strain_id),
 index(snpid),index(source_id),
 index(strain_id))ENGINE=InnoDB;
/***************************************************************************************
 snp_genotype_conflict: keeps SNP genotype conflict between sources
****************************************************************************************/
create table if not exists snp_genotype_conflict(
 snpid int unsigned default 0,
 source_id smallint default 0,
 strain_id  smallint default 0,
 genotype_allele char(1) not null,
 primary key(snpid,source_id,strain_id,genotype_allele),
 foreign key(snpid) references snp_position(snpid),
 foreign key(source_id) references snp_source(source_id),
 foreign key(strain_id) references snp_strain(strain_id),
 index(snpid),index(source_id),
 index(strain_id))ENGINE=InnoDB;

/***************************************************************************************
 snp_main: SNP main table
****************************************************************************************/
create table if not exists snp_main(
 snpid int unsigned default 0, 
 ref_allele char(1) not null, snp_allele char(1) not null,
 is_conflict tinyint default 0,
 is_intergenic tinyint default 1, 
 mutation_type tinyint default 0, 
 is_CpG_site tinyint default 0, 
 foreign key(snpid) references snp_position(snpid), 
 foreign key(mutation_type) references cgd_mutation_type(mutation_type), 
 index(snpid),index(is_CpG_site),index(is_intergenic),index(is_conflict),
 index(mutation_type ))ENGINE=InnoDB;  
/***************************************************************************************
 snp_main_allele_conflict: stores snp_allele conflict between sources
****************************************************************************************/
create table if not exists snp_main_allele_conflict(
   snpid int unsigned default 0,
   source_id smallint default 0, snp_allele char(1) not null,
   primary key(snpid,source_id,snp_allele),
   foreign key(snpid) references snp_position(snpid),                            
   foreign key(source_id) references snp_source(source_id),
   index(snpid),index(source_id))ENGINE=InnoDB;




