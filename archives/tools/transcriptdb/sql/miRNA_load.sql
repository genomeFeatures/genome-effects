/*use graber_transcriptdb_temptable;

/****************************************************************
 Author: Lucie Hutchins, Scientific Software Engineer
 Date :  February 2012
 Project: feature (protein domains, miRNA)
           annotation integration with graber_transcriptdb database
 
1. miRNA
 Description: data is downloaded from 
http://www.ebi.ac.uk/enright-srv/microcosm/cgi-bin/targets/v5/download.pl

2.protein domains
              
*****************************************************************/

use graber_transcriptdb;

/*drop table if exists miRNAs;*/
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

drop table if exists proteinDomains;
create table if not exists proteinDomains(
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

DROP TABLE IF EXISTS feature_load_log;
create table if not exists feature_load_log(
 table_name varchar(50),segment_name varchar(100),
 file_line_count int unsigned default 0,
 db_line_count int unsigned default 0
)ENGINE=MyISAM;
 



