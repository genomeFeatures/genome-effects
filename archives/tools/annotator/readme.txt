/****************************************************************
 The following steps are needed to load/update graber_transcriptdb
 Author : Lucie Hutchins
 Date :   July 2010
 
 Note: this process will run monthly or every two/three months
       to synchronize both our database and ucsc, MGI,and Ensembl
Check this documentation : snap-0.15-manual.pdf 
*****************************************************************/

1. download the current organism/db/build listing from ucsc browser
   run [ getCurrentOrganisms.pl ] script to generate the following config files:
   a.ucsc_organism_version_date.txt commas separated file with the following fields
        1.organism_group,
        2.organims_name, 
        3.ucsc_db, 
        4.organism_version
   b.ucsc_geneAnnotationConfig.txt commas separated file with the following fields
        1.organism_group,
        2.organims_name, 
        3.ucsc_db, 
        4.annotation tables list separated by ":"
   c. Configuration file
      contains list of config files to be used

2. download the annotations tables for each organism
   run: [ download_ucsc_annotations.pl ] 
   and generate downloaded_annotations.log with fields:
   a.organism_group,
   b.organims_name, 
   c.annotation table name, 
   d.downloaded annotation file name

4. generate db ids if new chromosome, organism, organism version,gene annotations,..
   by running load_Static_tables.pl
5. Run db_backup db_backup.pl
6. Run [ load_geneannotation.pl ] to generate the ucsc_dataload.sql, a load script,
   and all the data files to be loaded into the database . Also update:
   a.the chomosome
   b.the chromosome_data
   c.organism_group
   d.organism
   e.organism_version table
   f.gene_anotation table
 
 5. run the ucsc_dataload.sql to load/update database tables.
 
 6. generate dataload stattistics tables

Notes:
We have cases in the database where a given gene from one source
has multiple isoforms and all the isoforms have the same txStart and txEnd
The good part is that they are given the same txId, the bad part is that
the translation table stores the translation data using the txID, now we need to add
the transcript accession id to the translation table to avoid cases where
the same txID has multiple cdsStart and cdsEnds 
For example:
ENSG00000143257 has 6 transcript in hg18 and all these tx have
the same txStart and txEnd but different cdsStart and cdsEnds

Notes:
I have created load scripts to load annotation sources other UCSC into our database.
These source include:
1. Ensembl -> download_ensAnnotations.pl (although UCSC provides ensembl annotations,
  we find it difficult to track the Ensembl schema version we are currently using.
  also if we like to have more than one schema version of a given organism version,
  or additional annotations - for example exons accession id)
  only for mouse at the moment
2. MGI -> download_mgiAnnotations.pl -- this script should run once a month (the 4th) because
   MGI generates the load file every 4th of the month
3. Protein Domains (48 organisms) - download_proteinDomains.pl
  Data from uniprot has the following information 
  +-----------------------+--------------+------+-----+---------+-------+
| Field                 | Type         | Null | Key | Default | Extra |
+-----------------------+--------------+------+-----+---------+-------+
| organism_id           | smallint(6)  | YES  | MUL | 0       |       | 
| uniprot_protein_id    | varchar(255) | NO   | MUL | NULL    |       | 
| uniprot_protein_name  | varchar(255) | NO   |     | NULL    |       | 
| protein_evidence      | varchar(255) | NO   |     | NULL    |       | 
| pfams                 | varchar(255) | NO   | MUL | NULL    |       | 
| common_gene_name      | varchar(255) | NO   | MUL | NULL    |       | 
| other_transcript_name | varchar(255) | NO   | MUL | NULL    |       | 
| other_protein_id      | varchar(255) | NO   | MUL | NULL    |       | 
| other_gene_name       | varchar(255) | NO   | MUL | NULL    |       | 
| feature_type          | varchar(50)  | YES  | MUL | NULL    |       | 
| feature_aa_start      | smallint(6)  | YES  |     | 0       |       | 
| feature_aa_end        | smallint(6)  | YES  |     | 0       |       | 
| feature_name          | varchar(100) | YES  |     | NULL    |       | 
+-----------------------+--------------+------+-----+---------+-------+
including the sequence information (protein sequence)
+----------------+-------------+
| organism_name  | organism_id |
+----------------+-------------+
| A.gambiae      |           1 | 
| A.mellifera    |           2 | 
| C.briggsae     |           4 | 
| C.elegans      |           5 | 
| C.intestinalis |           6 | 
| C.remanei      |           8 | 
| Cat            |           9 | 
| Chicken        |          10 | 
| Chimp          |          11 | 
| Cow            |          12 | 
| D.ananassae    |          13 | 
| D.erecta       |          14 | 
| D.grimshawi    |          15 | 
| D.melanogaster |          16 | 
| D.mojavensis   |          17 | 
| D.persimilis   |          18 | 
| D.sechellia    |          20 | 
| D.simulans     |          21 | 
| D.virilis      |          22 | 
| D.yakuba       |          23 | 
| Dog            |          24 | 
| Elephant       |          25 | 
| Fugu           |          26 | 
| Gorilla        |          28 | 
| Guinea pig     |          29 | 
| Horse          |          30 | 
| Human          |          31 | 
| Lamprey        |          32 | 
| Lancelet       |          33 | 
| Lizard         |          34 | 
| Marmoset       |          35 | 
| Medaka         |          36 | 
| Mouse          |          38 | 
| Opossum        |          39 | 
| Orangutan      |          40 | 
| Panda          |          42 | 
| Pig            |          43 | 
| Platypus       |          44 | 
| Rabbit         |          45 | 
| Rat            |          46 | 
| Rhesus         |          47 | 
| S.cerevisiae   |          48 | 
| S.purpuratus   |          49 | 
| Sea hare       |          50 | 
| Sheep          |          51 | 
| Turkey         |          54 | 
| X.tropicalis   |          55 | 
| Zebrafish      |          57 | 



