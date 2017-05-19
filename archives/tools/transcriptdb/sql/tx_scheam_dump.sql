-- MySQL dump 10.11
--
-- Host: harlequin    Database: graber_transcriptdb
-- ------------------------------------------------------
-- Server version	5.0.96

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `cc_founders_genes`
--

DROP TABLE IF EXISTS `cc_founders_genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `cc_founders_genes` (
  `gene_prediction_id` tinyint(3) unsigned default NULL,
  `gene_id` varchar(100) default NULL,
  `gene_name` varchar(200) default NULL,
  `biotype` varchar(50) default NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `cc_founders_transcripts`
--

DROP TABLE IF EXISTS `cc_founders_transcripts`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `cc_founders_transcripts` (
  `gene_prediction_id` tinyint(3) unsigned default NULL,
  `transcript_id` varchar(100) default NULL,
  `transcript_name` varchar(255) default NULL,
  `protein_id` varchar(100) default NULL,
  KEY `gene_prediction_id` (`gene_prediction_id`),
  KEY `transcript_id` (`transcript_id`),
  KEY `protein_id` (`protein_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chromosome`
--

DROP TABLE IF EXISTS `chromosome`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `chromosome` (
  `chromosome_id` mediumint(8) unsigned NOT NULL auto_increment,
  `chromosome_name` varchar(255) default NULL,
  PRIMARY KEY  (`chromosome_id`),
  KEY `chromosome_name` (`chromosome_name`)
) ENGINE=InnoDB AUTO_INCREMENT=1656059 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `chromosome_data`
--

DROP TABLE IF EXISTS `chromosome_data`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `chromosome_data` (
  `organism_version_id` smallint(6) default '0',
  `chromosome_id` mediumint(8) unsigned default '0',
  `fileName` varchar(255) default NULL,
  `chrom_size` int(10) unsigned default '0',
  KEY `organism_version_id` (`organism_version_id`),
  KEY `chromosome_id` (`chromosome_id`),
  KEY `organism_version_id_2` (`organism_version_id`,`chromosome_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ensembl_exons`
--

DROP TABLE IF EXISTS `ensembl_exons`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ensembl_exons` (
  `exon_id` int(10) unsigned default '0',
  `organism_version_id` smallint(6) default NULL,
  `accession_id` varchar(50) default NULL,
  KEY `exon_id` (`exon_id`),
  KEY `accession_id` (`accession_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ensembl_genes`
--

DROP TABLE IF EXISTS `ensembl_genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ensembl_genes` (
  `gene_name` varchar(255) default NULL,
  `biotype` varchar(100) default NULL,
  `organism_version_id` smallint(6) default '0',
  `description` text,
  KEY `gene_name` (`gene_name`),
  KEY `biotype` (`biotype`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `est_align`
--

DROP TABLE IF EXISTS `est_align`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `est_align` (
  `organism_version_id` smallint(6) default '0',
  `matches` int(10) unsigned default '0',
  `misMatches` int(10) unsigned default '0',
  `strand` char(1) NOT NULL,
  `qName` char(12) NOT NULL,
  `qSize` int(10) unsigned default '0',
  `qStart` int(10) unsigned default '0',
  `qEnd` int(10) unsigned default '0',
  `chromosome_id` mediumint(8) unsigned default '0',
  `tStart` int(10) unsigned default '0',
  `tEnd` int(10) unsigned default '0',
  KEY `organism_version_id` (`organism_version_id`),
  KEY `organism_version_id_2` (`organism_version_id`,`chromosome_id`),
  KEY `qName` (`qName`),
  KEY `organism_version_id_3` (`organism_version_id`,`chromosome_id`,`strand`,`tStart`,`tEnd`),
  KEY `chromosome_id` (`chromosome_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `est_orientation`
--

DROP TABLE IF EXISTS `est_orientation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `est_orientation` (
  `organism_version_id` smallint(6) default '0',
  `qName` char(12) NOT NULL,
  `chromosome_id` mediumint(8) unsigned default '0',
  `tStart` int(10) unsigned default '0',
  `tEnd` int(10) unsigned default '0',
  `sizePolyA` smallint(6) default '0',
  `revSizePolyA` smallint(6) default '0',
  `signalPos` smallint(6) default '0',
  `revSignalPos` smallint(6) default '0',
  KEY `organism_version_id` (`organism_version_id`),
  KEY `organism_version_id_2` (`organism_version_id`,`chromosome_id`),
  KEY `qName` (`qName`),
  KEY `organism_version_id_3` (`organism_version_id`,`chromosome_id`,`tStart`,`tEnd`),
  KEY `chromosome_id` (`chromosome_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `exon`
--

DROP TABLE IF EXISTS `exon`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `exon` (
  `exon_id` int(10) unsigned NOT NULL auto_increment,
  `organism_version_id` smallint(6) default '0',
  `chromosome_id` mediumint(8) unsigned default '0',
  `strand` char(1) NOT NULL,
  `exon_start` int(10) unsigned default '0',
  `exon_end` int(10) unsigned default '0',
  PRIMARY KEY  (`exon_id`),
  KEY `organism_version_id` (`organism_version_id`),
  KEY `organism_version_id_2` (`organism_version_id`,`chromosome_id`),
  KEY `organism_version_id_3` (`organism_version_id`,`chromosome_id`,`strand`),
  KEY `organism_version_id_4` (`organism_version_id`,`chromosome_id`,`strand`,`exon_start`,`exon_end`),
  KEY `chromosome_id` (`chromosome_id`)
) ENGINE=MyISAM AUTO_INCREMENT=9484660 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `external_db`
--

DROP TABLE IF EXISTS `external_db`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `external_db` (
  `external_db_id` tinyint(4) default NULL,
  `db_name` varchar(50) default NULL,
  `schema_version` varchar(100) default NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `feature_load_log`
--

DROP TABLE IF EXISTS `feature_load_log`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `feature_load_log` (
  `table_name` varchar(50) default NULL,
  `segment_name` varchar(100) default NULL,
  `file_line_count` int(10) unsigned default '0',
  `db_line_count` int(10) unsigned default '0',
  `moddate` date default '0000-00-00'
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene_by_annotation`
--

DROP TABLE IF EXISTS `gene_by_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene_by_annotation` (
  `organism_version_id` smallint(6) default '0',
  `gene_name` varchar(255) NOT NULL default '',
  `transcript_id` int(10) unsigned NOT NULL default '0',
  `gene_prediction_id` tinyint(3) unsigned NOT NULL default '0',
  PRIMARY KEY  (`gene_name`,`transcript_id`,`gene_prediction_id`),
  KEY `gene_prediction_id` (`gene_prediction_id`),
  KEY `transcript_id` (`transcript_id`),
  KEY `gene_name` (`gene_name`),
  KEY `organism_version_id` (`organism_version_id`,`transcript_id`,`gene_prediction_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene_prediction`
--

DROP TABLE IF EXISTS `gene_prediction`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene_prediction` (
  `gene_prediction_id` tinyint(3) unsigned NOT NULL auto_increment,
  `gene_prediction_name` varchar(230) default NULL,
  `gene_prediction_desc` varchar(255) default NULL,
  PRIMARY KEY  (`gene_prediction_id`),
  KEY `gene_prediction_name` (`gene_prediction_name`)
) ENGINE=InnoDB AUTO_INCREMENT=62 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `gene_prediction_by_organism`
--

DROP TABLE IF EXISTS `gene_prediction_by_organism`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `gene_prediction_by_organism` (
  `organism_version_id` smallint(6) default '0',
  `gene_prediction_id` tinyint(3) unsigned default '0',
  `load_date` varchar(20) default NULL,
  KEY `organism_version_id` (`organism_version_id`),
  KEY `gene_prediction_id` (`gene_prediction_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `mgi_genes`
--

DROP TABLE IF EXISTS `mgi_genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `mgi_genes` (
  `mgi_id` varchar(25) default NULL,
  `mgi_symbol` varchar(100) default NULL,
  `biotype` varchar(100) default NULL,
  KEY `mgi_id` (`mgi_id`),
  KEY `mgi_symbol` (`mgi_symbol`),
  KEY `biotype` (`biotype`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `miRNAs`
--

DROP TABLE IF EXISTS `miRNAs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `miRNAs` (
  `organism_id` smallint(6) default '0',
  `name` varchar(50) default NULL,
  `chromosome_name` varchar(50) default NULL,
  `chromosome_id` smallint(6) default '0',
  `chrom_start` int(10) unsigned default '0',
  `chrom_end` int(10) unsigned default '0',
  `strand` char(1) NOT NULL,
  `transcript_id` int(10) unsigned default '0',
  `transcript_name` varchar(255) NOT NULL,
  `external_name` varchar(255) NOT NULL,
  KEY `transcript_name` (`transcript_name`),
  KEY `external_name` (`external_name`),
  KEY `organism_id` (`organism_id`,`chromosome_id`),
  KEY `chromosome_id` (`chromosome_id`,`strand`,`chrom_start`,`chrom_end`),
  KEY `transcript_id` (`transcript_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `novel2known_genes`
--

DROP TABLE IF EXISTS `novel2known_genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `novel2known_genes` (
  `novel_gene_name` varchar(255) NOT NULL,
  `known_gene_name` varchar(255) NOT NULL,
  KEY `novel_gene_name` (`novel_gene_name`),
  KEY `known_gene_name` (`known_gene_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `novel_genes`
--

DROP TABLE IF EXISTS `novel_genes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `novel_genes` (
  `organism_version_id` smallint(6) default '0',
  `transcript_name` varchar(255) NOT NULL,
  `novel_gene_name` varchar(255) NOT NULL,
  `gene_prediction_id` tinyint(3) unsigned default '0',
  KEY `organism_version_id` (`organism_version_id`),
  KEY `transcript_name` (`transcript_name`),
  KEY `novel_gene_name` (`novel_gene_name`),
  KEY `gene_prediction_id` (`gene_prediction_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `novel_genes_synonyms`
--

DROP TABLE IF EXISTS `novel_genes_synonyms`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `novel_genes_synonyms` (
  `organism_version_id` smallint(6) default '0',
  `gene_prediction_id` tinyint(3) unsigned default '0',
  `novel_gene_name` varchar(255) NOT NULL,
  `synonym_gene_name` varchar(255) NOT NULL,
  `tx_list` text,
  KEY `organism_version_id` (`organism_version_id`),
  KEY `gene_prediction_id` (`gene_prediction_id`),
  KEY `novel_gene_name` (`novel_gene_name`),
  KEY `synonym_gene_name` (`synonym_gene_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `organism`
--

DROP TABLE IF EXISTS `organism`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `organism` (
  `organism_id` smallint(6) NOT NULL auto_increment,
  `organism_name` varchar(250) default NULL,
  `organism_sc_name` varchar(250) default NULL,
  `organism_tax_id` smallint(6) default '0',
  `organism_group_id` tinyint(3) unsigned default '0',
  PRIMARY KEY  (`organism_id`),
  KEY `organism_group_id` (`organism_group_id`),
  KEY `organism_name` (`organism_name`),
  CONSTRAINT `organism_ibfk_1` FOREIGN KEY (`organism_group_id`) REFERENCES `organism_group` (`organism_group_id`)
) ENGINE=InnoDB AUTO_INCREMENT=64 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `organism_group`
--

DROP TABLE IF EXISTS `organism_group`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `organism_group` (
  `organism_group_id` tinyint(3) unsigned NOT NULL auto_increment,
  `organism_group_name` varchar(100) default NULL,
  PRIMARY KEY  (`organism_group_id`),
  KEY `organism_group_name` (`organism_group_name`)
) ENGINE=InnoDB AUTO_INCREMENT=7 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `organism_version`
--

DROP TABLE IF EXISTS `organism_version`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `organism_version` (
  `organism_version_id` smallint(6) NOT NULL auto_increment,
  `organism_id` smallint(6) default '0',
  `ucsc_db` varchar(100) default NULL,
  `coordinate_build` varchar(100) default NULL,
  PRIMARY KEY  (`organism_version_id`),
  KEY `organism_id` (`organism_id`),
  KEY `organism_version_id` (`organism_version_id`,`organism_id`),
  CONSTRAINT `organism_version_ibfk_1` FOREIGN KEY (`organism_id`) REFERENCES `organism` (`organism_id`)
) ENGINE=InnoDB AUTO_INCREMENT=100 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `pfamDomains`
--

DROP TABLE IF EXISTS `pfamDomains`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pfamDomains` (
  `organism_id` smallint(6) default '0',
  `uniprot_protein_id` varchar(25) NOT NULL,
  `pfam_id` varchar(25) NOT NULL,
  `pfam_name` varchar(100) default NULL,
  `pfam_type` varchar(25) NOT NULL,
  `feature_aa_start` smallint(6) default '0',
  `feature_aa_end` smallint(6) default '0',
  KEY `organism_id` (`organism_id`),
  KEY `uniprot_protein_id` (`uniprot_protein_id`),
  KEY `pfam_id` (`pfam_id`),
  KEY `pfam_type` (`pfam_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `proteinDomains`
--

DROP TABLE IF EXISTS `proteinDomains`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `proteinDomains` (
  `organism_id` smallint(6) default '0',
  `uniprot_protein_id` varchar(255) NOT NULL,
  `uniprot_protein_name` varchar(255) NOT NULL,
  `protein_evidence` varchar(255) NOT NULL,
  `pfams` varchar(255) NOT NULL,
  `common_gene_name` varchar(255) NOT NULL,
  `other_transcript_name` varchar(255) NOT NULL,
  `other_protein_id` varchar(255) NOT NULL,
  `other_gene_name` varchar(255) NOT NULL,
  `feature_type` varchar(50) default NULL,
  `feature_aa_start` smallint(6) default '0',
  `feature_aa_end` smallint(6) default '0',
  `feature_name` varchar(100) default NULL,
  KEY `organism_id` (`organism_id`),
  KEY `feature_type` (`feature_type`),
  KEY `uniprot_protein_id` (`uniprot_protein_id`),
  KEY `other_transcript_name` (`other_transcript_name`),
  KEY `other_protein_id` (`other_protein_id`),
  KEY `other_gene_name` (`other_gene_name`),
  KEY `pfams` (`pfams`),
  KEY `common_gene_name` (`common_gene_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `proteinSequence`
--

DROP TABLE IF EXISTS `proteinSequence`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `proteinSequence` (
  `organism_id` smallint(6) default '0',
  `uniprot_protein_id` varchar(255) NOT NULL,
  `isoform` varchar(50) default NULL,
  `sequence_len` smallint(6) default '0',
  `sequence` text,
  KEY `organism_id` (`organism_id`),
  KEY `uniprot_protein_id` (`uniprot_protein_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `query_terms`
--

DROP TABLE IF EXISTS `query_terms`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `query_terms` (
  `term` varchar(255) default NULL,
  `term_type_id` tinyint(4) default NULL,
  KEY `term` (`term`),
  KEY `term_type_id` (`term_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `repeat_name`
--

DROP TABLE IF EXISTS `repeat_name`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `repeat_name` (
  `repeat_name_id` smallint(5) unsigned NOT NULL auto_increment,
  `repeat_name` varchar(255) default NULL,
  PRIMARY KEY  (`repeat_name_id`)
) ENGINE=MyISAM AUTO_INCREMENT=1851 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `repeat_type`
--

DROP TABLE IF EXISTS `repeat_type`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `repeat_type` (
  `repeat_type_id` smallint(5) unsigned NOT NULL auto_increment,
  `repeat_type` varchar(255) default NULL,
  PRIMARY KEY  (`repeat_type_id`)
) ENGINE=MyISAM AUTO_INCREMENT=13 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `repeats`
--

DROP TABLE IF EXISTS `repeats`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `repeats` (
  `organism_version_id` smallint(6) NOT NULL default '0',
  `repeat_feature_id` int(10) unsigned NOT NULL default '0',
  `chromosome_id` mediumint(8) unsigned default '0',
  `seq_region_start` int(10) unsigned default '0',
  `seq_region_end` int(10) unsigned default '0',
  `strand` char(1) default NULL,
  `repeat_start` smallint(5) unsigned default '0',
  `repeat_end` smallint(5) unsigned default '0',
  `repeat_consensus_id` int(10) unsigned default '0',
  `repeat_name_id` smallint(5) unsigned default '0',
  `repeat_type_id` smallint(5) unsigned default '0',
  `repeatChrStart` int(10) unsigned default '0',
  `repeatChrEnd` int(10) unsigned default '0',
  PRIMARY KEY  (`organism_version_id`,`repeat_feature_id`),
  KEY `repeat_name_id` (`repeat_name_id`),
  KEY `organism_version_id` (`organism_version_id`),
  KEY `chromosome_id` (`chromosome_id`,`repeatChrStart`),
  KEY `repeat_feature_id` (`repeat_feature_id`),
  KEY `repeat_consensus_id` (`repeat_consensus_id`),
  KEY `repeat_type_id` (`repeat_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `rnaseq_sample`
--

DROP TABLE IF EXISTS `rnaseq_sample`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `rnaseq_sample` (
  `sample_id` smallint(5) unsigned NOT NULL auto_increment,
  `sample_name` varchar(50) default NULL,
  `sample_desc` varchar(255) default NULL,
  PRIMARY KEY  (`sample_id`)
) ENGINE=MyISAM AUTO_INCREMENT=286 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `rnaseq_transcript`
--

DROP TABLE IF EXISTS `rnaseq_transcript`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `rnaseq_transcript` (
  `organism_version_id` smallint(6) default '0',
  `gene_prediction_id` tinyint(3) unsigned default '0',
  `transcript_name` varchar(100) NOT NULL,
  `sample_id` smallint(5) unsigned default '0',
  KEY `organism_version_id` (`organism_version_id`),
  KEY `gene_prediction_id` (`gene_prediction_id`),
  KEY `transcript_name` (`transcript_name`),
  KEY `sample_id` (`sample_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `term_types`
--

DROP TABLE IF EXISTS `term_types`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `term_types` (
  `term_type_id` tinyint(4) NOT NULL auto_increment,
  `term_key_name` varchar(20) default NULL,
  `table_name` varchar(50) default NULL,
  `description` varchar(100) default NULL,
  PRIMARY KEY  (`term_type_id`)
) ENGINE=InnoDB AUTO_INCREMENT=33 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `transcript`
--

DROP TABLE IF EXISTS `transcript`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `transcript` (
  `transcript_id` int(10) unsigned NOT NULL auto_increment,
  `organism_version_id` smallint(6) default '0',
  `chromosome_id` mediumint(8) unsigned default '0',
  `strand` char(1) NOT NULL,
  `tx_start` int(10) unsigned default '0',
  `tx_end` int(10) unsigned default '0',
  PRIMARY KEY  (`transcript_id`),
  KEY `organism_version_id` (`organism_version_id`),
  KEY `organism_version_id_2` (`organism_version_id`,`chromosome_id`),
  KEY `organism_version_id_3` (`organism_version_id`,`chromosome_id`,`strand`),
  KEY `organism_version_id_4` (`organism_version_id`,`chromosome_id`,`strand`,`tx_start`,`tx_end`),
  KEY `chromosome_id` (`chromosome_id`)
) ENGINE=MyISAM AUTO_INCREMENT=4122924 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `transcript_by_annotation`
--

DROP TABLE IF EXISTS `transcript_by_annotation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `transcript_by_annotation` (
  `organism_version_id` smallint(6) default '0',
  `transcript_name` varchar(255) NOT NULL,
  `transcript_id` int(10) unsigned NOT NULL default '0',
  `gene_prediction_id` tinyint(3) unsigned NOT NULL default '0',
  PRIMARY KEY  (`transcript_name`,`transcript_id`,`gene_prediction_id`),
  KEY `transcript_name` (`transcript_name`),
  KEY `transcript_id` (`transcript_id`),
  KEY `gene_prediction_id` (`gene_prediction_id`),
  KEY `organism_version_id` (`organism_version_id`,`transcript_id`,`gene_prediction_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `transcript_exon`
--

DROP TABLE IF EXISTS `transcript_exon`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `transcript_exon` (
  `organism_version_id` smallint(6) NOT NULL default '0',
  `transcript_name` varchar(255) NOT NULL,
  `transcript_id` int(10) unsigned NOT NULL default '0',
  `exon_id` int(10) unsigned NOT NULL default '0',
  `gene_prediction_id` tinyint(3) unsigned NOT NULL default '0',
  `exon_frame` tinyint(4) NOT NULL default '-99',
  PRIMARY KEY  (`organism_version_id`,`transcript_name`,`transcript_id`,`exon_id`,`gene_prediction_id`,`exon_frame`),
  KEY `transcript_id` (`transcript_id`),
  KEY `exon_id` (`exon_id`),
  KEY `gene_prediction_id` (`gene_prediction_id`),
  KEY `organism_version_id` (`organism_version_id`,`transcript_id`,`gene_prediction_id`,`exon_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `transcript_translation`
--

DROP TABLE IF EXISTS `transcript_translation`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `transcript_translation` (
  `organism_version_id` smallint(6) default '0',
  `transcript_name` varchar(255) NOT NULL,
  `transcript_id` int(10) unsigned NOT NULL default '0',
  `gene_prediction_id` tinyint(3) unsigned NOT NULL default '0',
  `cdsStart` int(10) unsigned NOT NULL default '0',
  `cdsEnd` int(10) unsigned NOT NULL default '0',
  PRIMARY KEY  (`transcript_name`,`transcript_id`,`gene_prediction_id`,`cdsStart`,`cdsEnd`),
  KEY `organism_version_id` (`organism_version_id`),
  KEY `transcript_name` (`transcript_name`),
  KEY `transcript_id` (`transcript_id`),
  KEY `gene_prediction_id` (`gene_prediction_id`),
  KEY `transcript_id_2` (`transcript_id`,`gene_prediction_id`),
  KEY `transcript_id_3` (`transcript_id`,`gene_prediction_id`,`cdsStart`,`cdsEnd`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2014-01-26 15:20:48
