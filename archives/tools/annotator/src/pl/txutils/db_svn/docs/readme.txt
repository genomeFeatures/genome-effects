/****************************************************************
 The following steps are needed to load/update graber_transcriptdb
 Author : Lucie Hutchins
 Date :   July 2010
 
 Note: this process will run monthly or every two/three months
       to synchronize both our database and ucsc database
4472 avahi     20   0 20800 2740 1220 S    2  0.0 295:14.02 avahi-daemon
 
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

3. generate data format.For each annotation
  run [ generate_ucsc_annotationformat.pl ]
  to generate a configuration file containing the data fields,
  fieldcount, and field index accross organisms.
  the file annotations_format.log contains the following fields:
  a. annotation_name
  b. organism
  c. organism_group
  d. a list of fieldCount,fields and the corresponding index

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


