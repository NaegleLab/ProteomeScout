-- MySQL dump 10.13  Distrib 5.5.28, for debian-linux-gnu (i686)
--
-- Host: localhost    Database: ptmscout_dev
-- ------------------------------------------------------
-- Server version	5.5.28-0ubuntu0.12.04.2

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

CREATE DATABASE ptmscout;
USE ptmscout;

--
-- Table structure for table `GO`
--

DROP TABLE IF EXISTS `GO`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `GO` (
  `id` int(10) unsigned zerofill NOT NULL AUTO_INCREMENT,
  `aspect` enum('F','P','C') DEFAULT NULL,
  `GO` varchar(10) NOT NULL DEFAULT '',
  `term` text NOT NULL,
  `date` datetime NOT NULL,
  `version` varchar(10) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `uniqueEntry` (`aspect`,`GO`),
  KEY `GO` (`GO`)
) ENGINE=InnoDB AUTO_INCREMENT=45551 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `GO_hierarchy`
--

DROP TABLE IF EXISTS `GO_hierarchy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `GO_hierarchy` (
  `parent_id` int(10) unsigned NOT NULL,
  `child_id` int(10) unsigned NOT NULL,
  KEY `parent_id` (`parent_id`),
  KEY `child_id` (`child_id`),
  CONSTRAINT `GO_hierarchy_ibfk_1` FOREIGN KEY (`parent_id`) REFERENCES `GO` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `GO_hierarchy_ibfk_2` FOREIGN KEY (`child_id`) REFERENCES `GO` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `MS`
--

DROP TABLE IF EXISTS `MS`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MS` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `peptide` varchar(150) COLLATE latin1_bin NOT NULL DEFAULT '',
  `experiment_id` int(10) unsigned NOT NULL DEFAULT '0',
  `protein_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`) USING BTREE,
  UNIQUE KEY `uniqueMSEntry` (`peptide`,`experiment_id`,`protein_id`),
  KEY `FK_MS_Experiment` (`experiment_id`),
  KEY `FK_MS_protein` (`protein_id`),
  CONSTRAINT `FK_MS_protein` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`),
  CONSTRAINT `MS_ibfk_1` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=431598 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `MS_modifications`
--

DROP TABLE IF EXISTS `MS_modifications`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MS_modifications` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `MS_id` int(10) unsigned NOT NULL DEFAULT '0',
  `peptide_id` int(10) unsigned NOT NULL DEFAULT '0',
  `modification_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `FK_phosphopep_id` (`peptide_id`),
  KEY `FK_MS_id` (`MS_id`),
  KEY `modification_id` (`modification_id`),
  CONSTRAINT `MS_modifications_ibfk_1` FOREIGN KEY (`MS_id`) REFERENCES `MS` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `MS_modifications_ibfk_2` FOREIGN KEY (`modification_id`) REFERENCES `PTM` (`id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=649460 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `PTM`
--

DROP TABLE IF EXISTS `PTM`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `PTM` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(100) NOT NULL,
  `position` enum('anywhere','core','c-terminal','n-terminal') DEFAULT NULL,
  `accession` varchar(10) DEFAULT NULL,
  `target` varchar(1) DEFAULT NULL,
  `mono_mass_diff` float DEFAULT NULL,
  `avg_mass_diff` float DEFAULT NULL,
  `parent_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`,`accession`),
  KEY `parent_id` (`parent_id`),
  CONSTRAINT `PTM_ibfk_3` FOREIGN KEY (`parent_id`) REFERENCES `PTM` (`id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=4694 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `PTM_keywords`
--

DROP TABLE IF EXISTS `PTM_keywords`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `PTM_keywords` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `PTM_id` int(10) unsigned NOT NULL,
  `keyword` varchar(100) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `PTM_id_2` (`PTM_id`,`keyword`),
  KEY `PTM_id` (`PTM_id`),
  CONSTRAINT `PTM_keywords_ibfk_1` FOREIGN KEY (`PTM_id`) REFERENCES `PTM` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=2615 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `PTM_taxonomy`
--

DROP TABLE IF EXISTS `PTM_taxonomy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `PTM_taxonomy` (
  `PTM_id` int(10) unsigned NOT NULL,
  `taxon_id` int(10) unsigned NOT NULL,
  KEY `PTM_id` (`PTM_id`),
  KEY `taxon_id` (`taxon_id`),
  CONSTRAINT `PTM_taxonomy_ibfk_1` FOREIGN KEY (`PTM_id`) REFERENCES `PTM` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `PTM_taxonomy_ibfk_2` FOREIGN KEY (`taxon_id`) REFERENCES `taxonomy` (`node_id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ambiguity`
--

DROP TABLE IF EXISTS `ambiguity`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `ambiguity` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `peptide` varchar(100) NOT NULL DEFAULT '',
  `protein_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  UNIQUE KEY `uniqueConstraint_key` (`peptide`,`protein_id`),
  KEY `FK_ambiguity_1` (`protein_id`),
  CONSTRAINT `ambiguity_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=3223 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `data`
--

DROP TABLE IF EXISTS `data`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `data` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` enum('data','stddev') NOT NULL DEFAULT 'data',
  `units` varchar(20) NOT NULL,
  `run` varchar(20) NOT NULL DEFAULT 'average',
  `label` varchar(45) NOT NULL DEFAULT '',
  `priority` int(10) unsigned NOT NULL DEFAULT '0',
  `value` float DEFAULT NULL,
  `NA` tinyint(1) NOT NULL DEFAULT '0',
  `MS_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  KEY `FK_data_MS` (`MS_id`),
  CONSTRAINT `data_ibfk_1` FOREIGN KEY (`MS_id`) REFERENCES `MS` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=147124 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment`
--

DROP TABLE IF EXISTS `experiment`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` text NOT NULL,
  `author` text,
  `date` datetime NOT NULL,
  `description` text,
  `contact` varchar(45) DEFAULT NULL,
  `PMID` int(10) unsigned DEFAULT NULL,
  `URL` text,
  `published` tinyint(1) NOT NULL DEFAULT '0',
  `ambiguity` tinyint(1) NOT NULL DEFAULT '0',
  `export` tinyint(1) unsigned NOT NULL DEFAULT '0',
  `experiment_id` int(10) unsigned DEFAULT NULL,
  `dataset` text NOT NULL,
  `submitter` varchar(50) NOT NULL DEFAULT 'ptmscout_admin@mit.edu',
  `volume` int(11) DEFAULT NULL,
  `page_start` varchar(10) DEFAULT NULL,
  `page_end` varchar(10) DEFAULT NULL,
  `journal` varchar(45) DEFAULT NULL,
  `publication_month` enum('','january','february','march','april','may','june','july','august','september','october','november','december') DEFAULT NULL,
  `publication_year` int(4) unsigned DEFAULT NULL,
  `public` int(1) unsigned NOT NULL DEFAULT '1',
  `status` enum('preload','loading','loaded','error') NOT NULL DEFAULT 'preload',
  `submitter_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `public` (`public`),
  KEY `submitter_id` (`submitter_id`),
  KEY `ready` (`status`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `experiment_ibfk_2` FOREIGN KEY (`submitter_id`) REFERENCES `users` (`id`) ON DELETE SET NULL ON UPDATE CASCADE,
  CONSTRAINT `experiment_ibfk_3` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE SET NULL ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=1266 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment_condition`
--

DROP TABLE IF EXISTS `experiment_condition`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment_condition` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` enum('cell','tissue','drug','stimulus','environment') NOT NULL,
  `value` varchar(100) NOT NULL,
  `experiment_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `experiment_condition_ibfk_1` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=29527 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment_error`
--

DROP TABLE IF EXISTS `experiment_error`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment_error` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `experiment_id` int(10) unsigned NOT NULL,
  `line` int(10) unsigned NOT NULL,
  `accession` varchar(45) NOT NULL,
  `peptide` varchar(150) NOT NULL,
  `message` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `experiment_error_ibfk_1` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=343 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `expression`
--

DROP TABLE IF EXISTS `expression`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `expression` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probeset_id` varchar(45) NOT NULL DEFAULT '',
  `genechip` enum('gnf1h','gnf1m','HG-U133A') NOT NULL DEFAULT 'gnf1h',
  `species_id` int(10) unsigned NOT NULL,
  `name` text NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `probeset_id` (`probeset_id`),
  KEY `species_id` (`species_id`),
  CONSTRAINT `expression_ibfk_1` FOREIGN KEY (`species_id`) REFERENCES `species` (`id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=69793 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `expression_acc`
--

DROP TABLE IF EXISTS `expression_acc`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `expression_acc` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` enum('gene_symbol','refseq','uniprot','alias') NOT NULL,
  `value` varchar(45) NOT NULL,
  `probeset_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  KEY `value` (`value`),
  KEY `probeset` (`probeset_id`),
  CONSTRAINT `expression_acc_ibfk_1` FOREIGN KEY (`probeset_id`) REFERENCES `expression` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=424781 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `expression_collection`
--

DROP TABLE IF EXISTS `expression_collection`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `expression_collection` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(50) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB AUTO_INCREMENT=4 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `expression_samples`
--

DROP TABLE IF EXISTS `expression_samples`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `expression_samples` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probeset_id` int(10) unsigned NOT NULL,
  `collection_id` int(10) unsigned NOT NULL,
  `tissue_id` int(10) unsigned NOT NULL,
  `value` float NOT NULL,
  PRIMARY KEY (`id`),
  KEY `probeset_id` (`probeset_id`,`collection_id`,`tissue_id`),
  KEY `collection_id` (`collection_id`),
  KEY `tissue_id` (`tissue_id`),
  CONSTRAINT `expression_samples_ibfk_1` FOREIGN KEY (`collection_id`) REFERENCES `expression_collection` (`id`) ON UPDATE CASCADE,
  CONSTRAINT `expression_samples_ibfk_2` FOREIGN KEY (`tissue_id`) REFERENCES `expression_tissue` (`id`) ON UPDATE CASCADE,
  CONSTRAINT `expression_samples_ibfk_3` FOREIGN KEY (`probeset_id`) REFERENCES `expression` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=6715679 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `expression_tissue`
--

DROP TABLE IF EXISTS `expression_tissue`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `expression_tissue` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(50) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB AUTO_INCREMENT=197 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `invitations`
--

DROP TABLE IF EXISTS `invitations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `invitations` (
  `invited_email` varchar(50) NOT NULL,
  `experiment_id` int(10) unsigned NOT NULL,
  `inviting_user_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`invited_email`,`experiment_id`),
  KEY `user_id_index` (`inviting_user_id`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `invitations_ibfk_1` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `invitations_ibfk_2` FOREIGN KEY (`inviting_user_id`) REFERENCES `users` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `peptide`
--

DROP TABLE IF EXISTS `peptide`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `peptide` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `pep_aligned` varchar(15) CHARACTER SET latin1 NOT NULL DEFAULT '',
  `site_pos` int(10) unsigned NOT NULL DEFAULT '0',
  `site_type` char(1) CHARACTER SET latin1 NOT NULL DEFAULT '',
  `protein_domain_id` int(10) unsigned DEFAULT NULL,
  `protein_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`) USING BTREE,
  KEY `FK_peptide_protein` (`protein_id`),
  KEY `protein_domain_id` (`protein_domain_id`),
  CONSTRAINT `peptide_ibfk_4` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `peptide_ibfk_3` FOREIGN KEY (`protein_domain_id`) REFERENCES `protein_domain` (`id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=237084 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `peptide_predictions`
--

DROP TABLE IF EXISTS `peptide_predictions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `peptide_predictions` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `source` varchar(40) CHARACTER SET latin1 COLLATE latin1_bin NOT NULL DEFAULT 'scansite',
  `value` varchar(20) NOT NULL DEFAULT '',
  `score` float DEFAULT '0',
  `peptide_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  UNIQUE KEY `UNIQUE_pepId_source_value` (`peptide_id`,`source`,`value`),
  KEY `FK_peptide_prediction` (`peptide_id`),
  CONSTRAINT `peptide_predictions_ibfk_1` FOREIGN KEY (`peptide_id`) REFERENCES `peptide` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=508904 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `permissions`
--

DROP TABLE IF EXISTS `permissions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `permissions` (
  `user_id` int(10) unsigned NOT NULL,
  `experiment_id` int(10) unsigned NOT NULL,
  `access_level` enum('view','owner') NOT NULL DEFAULT 'view',
  KEY `user_id` (`user_id`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `permissions_ibfk_1` FOREIGN KEY (`user_id`) REFERENCES `users` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `permissions_ibfk_2` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein`
--

DROP TABLE IF EXISTS `protein`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `sequence` text NOT NULL,
  `species_id` int(10) unsigned NOT NULL,
  `acc_gene` varchar(30) DEFAULT NULL,
  `name` varchar(100) NOT NULL DEFAULT '',
  `date` datetime NOT NULL,
  PRIMARY KEY (`id`),
  KEY `species_id` (`species_id`),
  KEY `sequence` (`sequence`(20)),
  CONSTRAINT `protein_ibfk_1` FOREIGN KEY (`species_id`) REFERENCES `species` (`id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=38417 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_GO`
--

DROP TABLE IF EXISTS `protein_GO`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_GO` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `protein_id` int(10) unsigned NOT NULL DEFAULT '0',
  `GO_id` int(10) unsigned NOT NULL DEFAULT '0',
  `date` datetime NOT NULL,
  PRIMARY KEY (`id`),
  KEY `protein_id` (`protein_id`),
  KEY `GO_id` (`GO_id`),
  CONSTRAINT `protein_GO_ibfk_2` FOREIGN KEY (`GO_id`) REFERENCES `GO` (`id`) ON UPDATE CASCADE,
  CONSTRAINT `protein_GO_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=837022 DEFAULT CHARSET=utf8;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_acc`
--

DROP TABLE IF EXISTS `protein_acc`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_acc` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` varchar(30) DEFAULT NULL,
  `value` varchar(45) NOT NULL DEFAULT '',
  `protein_id` int(10) unsigned NOT NULL DEFAULT '0',
  `out_of_date` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  UNIQUE KEY `UNIQUE_Entry_forProtein` (`protein_id`,`value`,`type`) USING BTREE,
  KEY `FK_acc_protein` (`protein_id`),
  CONSTRAINT `protein_acc_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=226803 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_domain`
--

DROP TABLE IF EXISTS `protein_domain`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_domain` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `label` varchar(45) NOT NULL DEFAULT '',
  `start` int(10) unsigned NOT NULL DEFAULT '0',
  `stop` int(10) unsigned NOT NULL DEFAULT '0',
  `p_value` float NOT NULL DEFAULT '0',
  `source` varchar(45) NOT NULL DEFAULT 'pfam',
  `params` text,
  `protein_id` int(10) unsigned NOT NULL DEFAULT '0',
  `version` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  KEY `FK_pfam_protein` (`protein_id`),
  CONSTRAINT `protein_domain_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=109490 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_expression`
--

DROP TABLE IF EXISTS `protein_expression`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_expression` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `protein_id` int(10) unsigned zerofill NOT NULL DEFAULT '0000000000',
  `probeset_id` varchar(45) NOT NULL DEFAULT '',
  PRIMARY KEY (`id`),
  UNIQUE KEY `Unique` (`protein_id`,`probeset_id`),
  KEY `FK_protein_expression_1` (`protein_id`),
  KEY `FK_protein_expression_ann` (`probeset_id`),
  CONSTRAINT `FK_protein_expression_ann` FOREIGN KEY (`probeset_id`) REFERENCES `expression` (`probeset_id`),
  CONSTRAINT `protein_expression_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=32591 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `session`
--

DROP TABLE IF EXISTS `session`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `session` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `user_id` int(10) unsigned NOT NULL,
  `data_file` varchar(100) NOT NULL,
  `load_type` enum('new','reload','append','extension') NOT NULL,
  `parent_experiment` int(10) unsigned DEFAULT NULL,
  `change_description` text NOT NULL,
  `units` varchar(20) NOT NULL,
  `stage` enum('config','metadata','confirm','complete') NOT NULL DEFAULT 'config',
  `experiment_id` int(10) unsigned DEFAULT NULL,
  `date` datetime NOT NULL,
  PRIMARY KEY (`id`),
  KEY `user_id` (`user_id`,`parent_experiment`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `session_ibfk_1` FOREIGN KEY (`user_id`) REFERENCES `users` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `session_ibfk_2` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=424 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `session_columns`
--

DROP TABLE IF EXISTS `session_columns`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `session_columns` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `session_id` int(10) unsigned NOT NULL,
  `type` enum('hidden','data','stddev','accession','peptide','species','modification','run','none') NOT NULL DEFAULT 'none',
  `label` varchar(45) NOT NULL,
  `column_number` int(10) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  KEY `session_id` (`session_id`),
  CONSTRAINT `session_columns_ibfk_1` FOREIGN KEY (`session_id`) REFERENCES `session` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=5779 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `species`
--

DROP TABLE IF EXISTS `species`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `species` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(100) NOT NULL,
  `taxon_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`),
  KEY `taxon_id` (`taxon_id`),
  CONSTRAINT `species_ibfk_2` FOREIGN KEY (`taxon_id`) REFERENCES `taxonomy` (`node_id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=105 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `taxonomy`
--

DROP TABLE IF EXISTS `taxonomy`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `taxonomy` (
  `node_id` int(10) unsigned NOT NULL,
  `kingdom` varchar(1) NOT NULL,
  `name` varchar(100) NOT NULL,
  `strain` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`node_id`),
  KEY `kingdom` (`kingdom`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `users`
--

DROP TABLE IF EXISTS `users`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `users` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `username` varchar(30) NOT NULL,
  `salted_password` varchar(64) NOT NULL,
  `salt` varchar(10) NOT NULL,
  `name` varchar(50) NOT NULL,
  `email` varchar(50) NOT NULL,
  `institution` varchar(100) NOT NULL,
  `date_created` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `active` int(1) NOT NULL DEFAULT '0',
  `activation_token` varchar(50) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `username` (`username`,`email`)
) ENGINE=InnoDB AUTO_INCREMENT=2446 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2012-11-26 16:54:53
