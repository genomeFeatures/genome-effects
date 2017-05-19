/*genome_integrationdb

/*******************************************************************************************
 Author:  Lucie Hutchins, Scientific Software Engineer
 Date :   October 2013
 Project: transcripts integration database
 
 Description: Integrate trasncript annotations from various gene predictions and organisms.
              generate unique genomic regions for each transcript (exon-set),
              and exon for all the organisms and annotations
Also provide a repository for novel transcripts generated from
various rnaSeq analysis -
Requirement to add a transcriptome in genome_integrationdb 
1. Have a publication associated with your transcriptome
2. Have the transcriptome file in gene prediction format (gtf format in progress)
              
*******************************************************************************************/
use genomeeffectdb;
/***** transcript core ***********/
drop table if exists transcript_exon;
drop table if exists gene_by_annotation ;
drop table if exists transcript_by_annotation;
drop table if exists transcript_translation ;
drop table if exists rnaseq_transcript; 
drop table if exists transcript;
drop table if exists exon;
drop table if exists unique_regions;
drop table if exists unique_regions_synonym;
/*** features ************/
drop table if exists gene_prediction;
drop table if EXISTS organism_version;
drop table if exists organism;
drop table if exists organism_group;
drop table if exists chromosome;
/****** Data log group starts here *****************************************
/****************************************************************************
Log table, keeps track of when a give annotation was loaded into the database
table_name -> organism_version - annotation 
*****************************************************************************/
create table if not exists feature_load_log(
  table_name varchar(50),segment_name varchar(100),
  file_line_count int unsigned default 0,
  db_line_count  int unsigned default 0,
  moddate  date default "0000-00-00", 
  index(table_name)
)ENGINE=MyISAM;
/****************************************************************************
Log table, keeps track of the schema versions of downloads from db external sources 
(MGI, Ensembl,)
*****************************************************************************/
create table if not exists external_db(
 external_db_id tinyint,
 db_name varchar(50),
 schema_version varchar(100))ENGINE=MyISAM;
/****** Data log group ends here *****************************************/

/****** Core group starts here *****************************************/
/*************************************************************************
  list of all the chromosomes/scafolds found  in UCSC databases.
Note:
  Handy for mapping between chromosome_id and chromosome_name
**************************************************************************/
CREATE TABLE IF NOT EXISTS chromosome(
 chromosome_id MEDIUMINT UNSIGNED primary key auto_increment,
 chromosome_name varchar(255),
 index(chromosome_name)
)ENGINE=InnoDB;

/*************************************************************************
 A table storing different classes of organisms
**************************************************************************/
CREATE TABLE IF NOT EXISTS organism_group(
 organism_group_id tinyint unsigned PRIMARY KEY AUTO_INCREMENT,
 organism_group_name varchar(100),
 index(organism_group_name)
)ENGINE=InnoDB;

/*************************************************************************
  list of all the organisms found  in UCSC 
  
Note: 
  handy for query using either organism common name, organism scientific name,
  or organism taxonomy id
 organism_id is the same as the organism taxonomy id from NCBI
***************************************************************************/
CREATE TABLE IF NOT EXISTS organism(
 organism_id SMALLINT unsigned PRIMARY KEY,
 organism_name varchar(250),
 organism_sc_name VARCHAR(250),
 organism_group_id tinyint unsigned default 0,
 INDEX(organism_group_id),
 index(organism_name),
 FOREIGN KEY(organism_group_id) REFERENCES organism_group(organism_group_id)
 )ENGINE=InnoDB;

CREATE TABLE IF NOT EXISTS organism_version(
 organism_version_id SMALLINT PRIMARY KEY AUTO_INCREMENT,
 organism_id SMALLINT unsigned default 0,
 ucsc_db varchar(100),
 coordinate_build varchar(100),
 index(organism_id), index(organism_version_id,organism_id),
 FOREIGN KEY(organism_id) REFERENCES organism(organism_id)
)ENGINE=InnoDB;
/***************************************************************************
  list of all the gene predictions found  in UCSC for each organism
  ALSO some will be from some new sources
Note:
 handy for mapping between organism and gene predictions
****************************************************************************/
CREATE TABLE IF NOT EXISTS gene_prediction(
 gene_prediction_id tinyint unsigned PRIMARY KEY AUTO_INCREMENT,
 gene_prediction_name varchar(230),
 gene_prediction_desc varchar(255),
 source_link varchar(255),
 index(gene_prediction_name)
)ENGINE=InnoDB;
/***********************************************************
classifies features with found in our integrated database
1. gene
2. transcript
3. repeat
4. snp
5. microRNA
6. microsatelite
7. protein domain
8. protein
9. exon
10. gap
11. QTL
12. Intron
************************************************************/
create table if not exists feature_type(
 feature_type_id smallint unsigned primary key auto_increment,
 feature_type_name varchar(255)
)ENGINE=InnoDB;
insert into feature_type(feature_type_name) values("Gene","Transcript","Exon","SNP","Repeat","Protein","Protein Domain","Microsatelite");
insert into feature_type(feature_type_name) values("miRNA","Intron","UTR","QTL","Gap");
create table if not exists feature(
 feature_id int unsigned primary key auto_increment,
 feature_type_id smallint unsigned default 0,feature_name varchar(255),
 foreign key(feature_type_id) references feature_type(feature_type_id),
 index(feature_type_id),index(feature_name)
)ENGINE=InnoDB;
/****** a table dump loaded from ucsc.uniProt.featureClass*********/
create table if not exists uniprot_featureClass(
 id smallint unsigned primary key ,
 val varchar(255),index(val)
)ENGINE=InnoDB;
/****** a table dump loaded from ucsc.uniProt.featureType*********/
create table if not exists uniprot_featureType(
 id smallint unsigned primary key ,
 val varchar(255),index(val)
)ENGINE=InnoDB;
/****** a table dump loaded from ucsc.uniProt.commonName***/
create table if not exists uniprot_commonName(
 taxon_id int unsigned default 0,
 val varchar(255),
)ENGINE=InnoDB;
/**********************************************************
 Repeats associated tables
 From ensembl schema tables repeat_feature ,repeat_consensus
***********************************************************/
create table if not exists repeat_name(
  repeat_name_id smallint unsigned primary key auto_increment,
  repeat_name varchar(255))ENGINE=InnoDB;
create table if not exists repeat_type(
  repeat_type_id smallint unsigned primary key auto_increment,
  repeat_type varchar(255))ENGINE=InnoDB;;

/*****************************************************************************
 list of all the transcripts
 found  in UCSC, MGI, Enesmebl for each organismVersion/gene_prediction
 This table generates unique id (transcript_id) for a given
 organism_version_id/chromosome_id/strand/chrStar/chrEnd

*****************************************************************************/
CREATE TABLE IF NOT EXISTS transcript(
 transcript_id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
 organism_version_id SMALLINT default 0,
 chromosome_id  MEDIUMINT UNSIGNED default 0,
 strand char not null,
 tx_start  INT UNSIGNED default 0,
 tx_end  INT UNSIGNED default 0,
 INDEX(organism_version_id),
 index(organism_version_id,chromosome_id),
 index(organism_version_id,chromosome_id,strand),
 INDEX(organism_version_id,chromosome_id,strand,tx_start,tx_end),
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(chromosome_id) REFERENCES chromosome(chromosome_id)
)ENGINE=MyISAM;

/**********************************************************************
 mapping between transcript and gene prediction source
 Handy for some types of queries 
***********************************************************************/
CREATE TABLE IF NOT EXISTS transcript_by_annotation(
 organism_version_id SMALLINT default 0,
 transcript_name varchar(255) not null,
 transcript_id INT UNSIGNED default 0,
 gene_prediction_id tinyint unsigned default 0,
 feature_type_id smallint unsigned default 0,
 PRIMARY KEY(transcript_name,transcript_id,gene_prediction_id,feature_type_id),
 index(transcript_name),
 INDEX(transcript_id),
 INDEX(gene_prediction_id),index(feature_type_id)),
 index(organism_version_id,transcript_id,gene_prediction_id),
 FOREIGN KEY(gene_prediction_id) REFERENCES gene_prediction(gene_prediction_id),
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(feature_type_id ) REFERENCES feature_type(feature_type_id),
 FOREIGN KEY(transcript_id) REFERENCES transcript(transcript_id)
)ENGINE=MyISAM;

/************************************************************
 Transcript - Translation(CDS) mapping
 
*************************************************************/
CREATE TABLE IF NOT EXISTS transcript_translation(
 organism_version_id SMALLINT default 0,
 transcript_name varchar(255) not null,
 transcript_id INT UNSIGNED default 0,
 gene_prediction_id tinyint unsigned default 0,
 cdsStart  INT UNSIGNED default 0,
 cdsEnd	  INT UNSIGNED default 0,
 PRIMARY KEY(transcript_name,transcript_id,gene_prediction_id,cdsStart,cdsEnd),
 INDEX(organism_version_id),
 index(transcript_name),
 INDEX(transcript_id),
 INDEX(gene_prediction_id),
 index(transcript_id,gene_prediction_id),
 INDEX(transcript_id,gene_prediction_id,cdsStart,cdsEnd),
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(transcript_id) REFERENCES transcript(transcript_id),
 FOREIGN KEY(gene_prediction_id) REFERENCES gene_prediction(gene_prediction_id)
 )ENGINE=MyISAM;

/*****************************************************************************
 list of all the exons
 found  in UCSC UCSC, MGI, Enesmebl for each organismVersion/gene_prediction
 This table generates unique id (exon_id) for a given
 organism_version_id/chromosome_id/strand/chrStar/chrEnd
*****************************************************************************/
CREATE TABLE IF NOT EXISTS exon(
 exon_id INT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
 organism_version_id SMALLINT default 0,
 chromosome_id  MEDIUMINT UNSIGNED default 0,
 strand char not null,
 exon_start  INT UNSIGNED default 0,
 exon_end  INT UNSIGNED default 0,
 INDEX(organism_version_id),
 index(organism_version_id,chromosome_id),
 index(organism_version_id,chromosome_id,strand),
 INDEX(organism_version_id,chromosome_id,strand,exon_start,exon_end),
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(chromosome_id) REFERENCES chromosome(chromosome_id)
)ENGINE=MyISAM;
/**********************************************************************
 mapping between transcript and exon by gene prediction source
***********************************************************************/
CREATE TABLE IF NOT EXISTS transcript_exon(
 organism_version_id SMALLINT default 0,
 transcript_name varchar(255) not null,
 transcript_id INT UNSIGNED default 0,
 exon_id INT UNSIGNED default 0,
 gene_prediction_id tinyint unsigned default 0,
 exon_frame tinyint default -99,
 PRIMARY KEY(organism_version_id,transcript_name,transcript_id,exon_id,gene_prediction_id,exon_frame),
 FOREIGN KEY(transcript_id) REFERENCES transcript(transcript_id),
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(exon_id) REFERENCES exon(exon_id),
 FOREIGN KEY(gene_prediction_id) REFERENCES gene_prediction(gene_prediction_id),
 INDEX(transcript_id),
 INDEX(exon_id),index(gene_prediction_id),
 index(organism_version_id,transcript_id,gene_prediction_id,exon_id)
)ENGINE=MyISAM;
/**********************************************************************
 mapping between transcript  and gene by gene prediction source
***********************************************************************/
CREATE TABLE IF NOT EXISTS gene_by_annotation(
 organism_version_id SMALLINT default 0,
 gene_name VARCHAR(255),
 transcript_name varchar(255),
 transcript_id INT UNSIGNED default 0,
 gene_prediction_id tinyint unsigned default 0,
 PRIMARY KEY(gene_name,transcript_name,transcript_id,gene_prediction_id),
 INDEX(gene_prediction_id),
 INDEX(transcript_id),
 index(gene_name),
 index(organism_version_id,transcript_id,gene_prediction_id),
 FOREIGN KEY(transcript_id) REFERENCES transcript(transcript_id),
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(gene_prediction_id) REFERENCES gene_prediction(gene_prediction_id)
 )ENGINE=MyISAM;
 
/*******************************************************************
 A table to store extra information from Ensembl gene annotations
 updated by organism_version_id
*******************************************************************/
create table if not exists genes(
  gene_name varchar(255),
  biotype varchar(100),
  description text,
  source_id smallint default 0,
  FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
  index(gene_name),index(biotype),index(organism_version_id)
)ENGINE=MyISAM;

/***************************************************
 A table to store extra information for Ensembl exon
 updated by organism_version_id
****************************************************/
create table if not exists exon_desc(
 exon_id int unsigned default 0,
 organism_version_id smallint default 0,
 accession_id varchar(50), source_id smallint default 0,
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(exon_id) REFERENCES exon(exon_id),
 INDEX(exon_id),index(organism_version_id)
)ENGINE=MyISAM;

/******** Additional Genome Features ************************
/********************************** RNA seq data load ********
source name: rnaSeq isoform quantification and qualification,
From Nazira
************************************************************/
create table if not exists rnaseq_sample(
  sample_id smallint unsigned primary key auto_increment,
  sample_name varchar(50),sample_desc varchar(255),
  sample_tissue varchar(50),sample_age varchar(50),genetic_structure varchar(100));

/*create table if not exists rnaseq_sample_date(
   date_id tinyint unsigned primary key auto_increment,
   date_name varchar(50));
/******* Transcript sample mapping ******************/
create table if not exists rnaseq_transcript(
 organism_version_id smallint default 0,
 gene_prediction_id tinyint unsigned default 0,
 transcript_name varchar(100) not null,
 sample_id smallint unsigned default 0,
 -- date_id tinyint unsigned default 0,
 foreign key(organism_version_id) references organism_version(organism_version_id),
 -- foreign key(transcript_id) references transcript(transcript_id),
 foreign key(sample_id) references rnaseq_sample(sample_id),
 -- foreign key(date_id) references rnaseq_sample_date(date_id),
 FOREIGN KEY(gene_prediction_id) REFERENCES gene_prediction(gene_prediction_id),
 index(organism_version_id),index(gene_prediction_id),index(transcript_name),index(sample_id)
 );

create table if not exists repeats(
  organism_version_id smallint default 0,
  repeat_feature_id int unsigned default 0,
  chromosome_id mediumint unsigned default 0,
  seq_region_start int unsigned default 0,
  seq_region_end int unsigned default 0,
  strand char,
  repeat_start smallint unsigned default 0,
  repeat_end smallint unsigned default 0,
  repeat_consensus_id int unsigned default 0,
  repeat_name_id  smallint unsigned default 0,
  repeat_type_id  smallint unsigned default 0,
  repeatChrStart  int unsigned default 0,
  repeatChrEnd  int unsigned default 0, 
  primary key(organism_version_id,repeat_feature_id),
  foreign key(organism_version_id) references organism_version(organism_version_id),
  foreign key(chromosome_id) references chromosome(chromosome_id),
  foreign key(repeat_name_id) references repeat_name(repeat_name_id),
  foreign key(repeat_type_id) references repeat_type(repeat_type_id),
  index(organism_version_id),index(chromosome_id,repeatChrStart),
  index(repeat_feature_id),index(repeat_consensus_id),index(repeat_type_id)
)ENGINE=MyISAM;
/**********************************************************
 Protein Domain table (UNIPROT-SWISSPROT)-- mapping protein/id/gene trough ensembl
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
***********************************************************/
create table if not exists proteinDomains(
  segment_id SMALLINT default 0,
  organism_id SMALLINT default 0,
  uniprot_protein_id varchar(255) not null,
  uniprot_protein_name varchar(255) not null,
  protein_evidence varchar(255) not null,
  pfams varchar(255) not null,
  common_gene_name varchar(255) not null,
  other_transcript_name varchar(255) not null,
  other_protein_id varchar(255) not null,
  other_gene_name  varchar(255) not null,
  feature_type varchar(50),
  feature_aa_start smallint default 0,
  feature_aa_end smallint default 0,
  feature_name varchar(100),index(segment_id),
  index(organism_id),index(feature_type),index(uniprot_protein_id),index(other_transcript_name),
  index(other_protein_id),index(other_gene_name),index(pfams),index(common_gene_name),
  FOREIGN KEY(organism_id) REFERENCES organism(organism_id)
)ENGINE=MyISAM;
/**********************************************************
 miRNA table--  mapping miRNA trough ensembl/ucsc
http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/download.pl
***********************************************************/
create table if not exists miRNAs(
  organism_id SMALLINT default 0,
  name varchar(50),
  chromosome_name varchar(50),
  chromosome_id SMALLINT default 0,
  chrom_start  INT UNSIGNED default 0,
  chrom_end  INT UNSIGNED default 0,
  strand char not null,
  transcript_id INT UNSIGNED default 0,
  transcript_name varchar(255) not null,
  external_name varchar(255) not null,
  index(transcript_name),
  index(external_name),
  index(organism_id,chromosome_id),
  index(chromosome_id,strand,chrom_start,chrom_end),
  index(transcript_id),
  FOREIGN KEY(organism_id) REFERENCES organism(organism_id),
  FOREIGN KEY(transcript_id) REFERENCES transcript(transcript_id),
  FOREIGN KEY(chromosome_id) REFERENCES chromosome(chromosome_id)
)ENGINE=MyISAM;

/************************************************************************
 organism - est alignment mapping
 UCSC
************************************************************************/
CREATE TABLE IF NOT EXISTS est_align(
  organism_version_id SMALLINT default 0,
  matches int unsigned default 0,
  misMatches int unsigned default 0,
  strand char(1) not null,
  qName char(12) not null,
  qSize int unsigned default 0,
  qStart int unsigned default 0,
  qEnd int unsigned default 0,
  chromosome_id  MEDIUMINT UNSIGNED default 0,
  tStart int unsigned default 0,
  tEnd int unsigned default 0,
  index(organism_version_id),index(organism_version_id,chromosome_id),index(qName),
  index(organism_version_id,chromosome_id,strand,tStart,tEnd),
  FOREIGN KEY(chromosome_id) REFERENCES chromosome(chromosome_id),
  FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id)
)ENGINE=MyISAM;

/************************************************************
 organism - est orientation mapping
 UCSC
*************************************************************/
CREATE TABLE IF NOT EXISTS est_orientation(
  organism_version_id SMALLINT default 0,
  qName char(12) not null,
  chromosome_id  MEDIUMINT UNSIGNED default 0, 
  tStart int unsigned default 0,
  tEnd int unsigned default 0,
  sizePolyA smallint default 0, 
  revSizePolyA smallint default 0,
  signalPos smallint default 0, 
  revSignalPos smallint default 0,
  index(organism_version_id),index(organism_version_id,chromosome_id),index(qName),
  index(organism_version_id,chromosome_id,tStart,tEnd),
  FOREIGN KEY(chromosome_id) REFERENCES chromosome(chromosome_id),
  FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id)
)ENGINE=MyISAM;

create table if not exists novel_genes(
 organism_version_id SMALLINT default 0,
 transcript_name varchar(255) not null,
 novel_gene_name varchar(255) not null,
 gene_prediction_id tinyint unsigned default 0,
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(gene_prediction_id) REFERENCES gene_prediction(gene_prediction_id),
 index(organism_version_id),index(transcript_name),index(novel_gene_name),index(gene_prediction_id)
);
create table if not exists novel_genes_synonyms(
 organism_version_id SMALLINT default 0,
 gene_prediction_id tinyint unsigned default 0,
 novel_gene_name varchar(255) not null,
 synonym_gene_name varchar(255) not null,
 tx_list text,
 FOREIGN KEY(organism_version_id) REFERENCES organism_version(organism_version_id),
 FOREIGN KEY(gene_prediction_id) REFERENCES gene_prediction(gene_prediction_id),
 index(novel_gene_name),index(synonym_gene_name)
);

