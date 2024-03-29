-- MySQL dump 10.13  Distrib 5.5.29, for debian-linux-gnu (x86_64)
--
-- Host: localhost    Database: ptmscout
-- ------------------------------------------------------
-- Server version	5.5.29-0ubuntu0.12.04.1

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
) ENGINE=InnoDB AUTO_INCREMENT=138304 DEFAULT CHARSET=utf8;
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
  `query_accession` varchar(45) COLLATE latin1_bin NOT NULL,
  `peptide` varchar(150) COLLATE latin1_bin NOT NULL DEFAULT '',
  `experiment_id` int(10) unsigned NOT NULL DEFAULT '0',
  `protein_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`) USING BTREE,
  KEY `FK_MS_Experiment` (`experiment_id`),
  KEY `FK_MS_protein` (`protein_id`),
  CONSTRAINT `FK_MS_protein` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`),
  CONSTRAINT `MS_ibfk_1` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=1249284 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `MS_ambiguity`
--

DROP TABLE IF EXISTS `MS_ambiguity`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MS_ambiguity` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `locus` varchar(30) NOT NULL,
  `alt_accession` varchar(20) NOT NULL DEFAULT '',
  `ms_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `uniqueConstraint_key` (`alt_accession`,`ms_id`),
  KEY `FK_MS` (`ms_id`),
  CONSTRAINT `MS_ambiguity_ibfk_1` FOREIGN KEY (`ms_id`) REFERENCES `MS` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=116327 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `MS_annotations`
--

DROP TABLE IF EXISTS `MS_annotations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MS_annotations` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `MS_id` int(10) unsigned NOT NULL,
  `type_id` int(10) unsigned NOT NULL,
  `value` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `MS_id` (`MS_id`),
  KEY `label` (`value`),
  KEY `set_id` (`type_id`),
  CONSTRAINT `MS_annotations_ibfk_2` FOREIGN KEY (`type_id`) REFERENCES `annotations` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `MS_annotations_ibfk_1` FOREIGN KEY (`MS_id`) REFERENCES `MS` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=22201 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `MS_data`
--

DROP TABLE IF EXISTS `MS_data`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `MS_data` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` enum('data','stddev') NOT NULL DEFAULT 'data',
  `units` varchar(20) NOT NULL,
  `run` varchar(20) NOT NULL DEFAULT 'average',
  `label` varchar(45) NOT NULL DEFAULT '',
  `priority` int(10) unsigned NOT NULL DEFAULT '0',
  `value` float DEFAULT NULL,
  `MS_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  KEY `FK_data_MS` (`MS_id`),
  CONSTRAINT `MS_data_ibfk_1` FOREIGN KEY (`MS_id`) REFERENCES `MS` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=252655 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=1788634 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=4710 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=2645 DEFAULT CHARSET=latin1;
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
-- Table structure for table `annotation_sets`
--

DROP TABLE IF EXISTS `annotation_sets`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotation_sets` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` text NOT NULL,
  `owner_id` int(10) unsigned NOT NULL,
  `experiment_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  KEY `owner_id` (`owner_id`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `annotation_sets_ibfk_2` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `annotation_sets_ibfk_1` FOREIGN KEY (`owner_id`) REFERENCES `users` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=5 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `annotations`
--

DROP TABLE IF EXISTS `annotations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `annotations` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `set_id` int(10) unsigned NOT NULL,
  `type` enum('cluster','numeric','nominative') NOT NULL,
  `name` varchar(100) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `one_name_per_set` (`set_id`,`name`),
  KEY `set_id` (`set_id`),
  CONSTRAINT `annotations_ibfk_1` FOREIGN KEY (`set_id`) REFERENCES `annotation_sets` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=1312 DEFAULT CHARSET=latin1;
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
  `experiment_id` int(10) unsigned DEFAULT NULL,
  `dataset` text NOT NULL,
  `volume` int(11) DEFAULT NULL,
  `page_start` varchar(10) DEFAULT NULL,
  `page_end` varchar(10) DEFAULT NULL,
  `journal` varchar(45) DEFAULT NULL,
  `publication_month` enum('','january','february','march','april','may','june','july','august','september','october','november','december') DEFAULT NULL,
  `publication_year` int(4) unsigned DEFAULT NULL,
  `public` int(1) unsigned NOT NULL DEFAULT '1',
  `job_id` int(10) unsigned DEFAULT NULL,
  `submitter_id` int(10) unsigned DEFAULT NULL,
  `modified_residues` varchar(26) NOT NULL,
  `type` enum('compendia','experiment','dataset') NOT NULL DEFAULT 'experiment',
  PRIMARY KEY (`id`),
  KEY `public` (`public`),
  KEY `submitter_id` (`submitter_id`),
  KEY `experiment_id` (`experiment_id`),
  KEY `modified_residues` (`modified_residues`),
  KEY `type` (`type`),
  CONSTRAINT `experiment_ibfk_2` FOREIGN KEY (`submitter_id`) REFERENCES `users` (`id`) ON DELETE SET NULL ON UPDATE CASCADE,
  CONSTRAINT `experiment_ibfk_3` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE SET NULL ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=1419 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment_PTM`
--

DROP TABLE IF EXISTS `experiment_PTM`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment_PTM` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `experiment_id` int(10) unsigned NOT NULL,
  `PTM_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `experiment_id_2` (`experiment_id`,`PTM_id`),
  KEY `experiment_id` (`experiment_id`),
  KEY `PTM_id` (`PTM_id`),
  CONSTRAINT `experiment_PTM_ibfk_1` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `experiment_PTM_ibfk_2` FOREIGN KEY (`PTM_id`) REFERENCES `PTM` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=580 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=29594 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=284490 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiment_subsets`
--

DROP TABLE IF EXISTS `experiment_subsets`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `experiment_subsets` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `owner_id` int(10) unsigned NOT NULL,
  `experiment_id` int(10) unsigned NOT NULL,
  `annotation_set_id` int(10) unsigned NOT NULL,
  `name` varchar(100) NOT NULL,
  `foreground_query` blob NOT NULL,
  `background_query` blob NOT NULL,
  PRIMARY KEY (`id`),
  KEY `name` (`name`),
  KEY `owner_id` (`owner_id`),
  KEY `experiment_id` (`experiment_id`),
  KEY `annotation_set_id` (`annotation_set_id`),
  CONSTRAINT `experiment_subsets_ibfk_3` FOREIGN KEY (`annotation_set_id`) REFERENCES `annotation_sets` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `experiment_subsets_ibfk_1` FOREIGN KEY (`owner_id`) REFERENCES `users` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `experiment_subsets_ibfk_2` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
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
-- Table structure for table `jobs`
--

DROP TABLE IF EXISTS `jobs`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `jobs` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `status` enum('configuration','in queue','started','finished','error') NOT NULL DEFAULT 'configuration',
  `failure_reason` text NOT NULL,
  `stage` varchar(20) NOT NULL,
  `progress` int(10) unsigned NOT NULL,
  `max_progress` int(10) unsigned NOT NULL,
  `status_url` varchar(250) DEFAULT NULL,
  `resume_url` varchar(250) DEFAULT NULL,
  `result_url` varchar(250) DEFAULT NULL,
  `name` varchar(250) NOT NULL,
  `type` enum('load_experiment','load_annotations','load_dataset') NOT NULL,
  `user_id` int(10) unsigned NOT NULL,
  `created` datetime NOT NULL,
  `restarted` datetime DEFAULT NULL,
  `finished` datetime DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `user_id` (`user_id`),
  CONSTRAINT `jobs_ibfk_1` FOREIGN KEY (`user_id`) REFERENCES `users` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=58 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `peptide`
--

DROP TABLE IF EXISTS `peptide`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `peptide` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `scansite_date` datetime DEFAULT NULL,
  `pep_aligned` varchar(15) CHARACTER SET latin1 NOT NULL DEFAULT '',
  `site_pos` int(10) unsigned NOT NULL DEFAULT '0',
  `site_type` char(1) CHARACTER SET latin1 NOT NULL DEFAULT '',
  `protein_domain_id` int(10) unsigned DEFAULT NULL,
  `protein_id` int(10) unsigned NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`) USING BTREE,
  UNIQUE KEY `site_pos` (`site_pos`,`protein_id`),
  KEY `FK_peptide_protein` (`protein_id`),
  KEY `protein_domain_id` (`protein_domain_id`),
  CONSTRAINT `peptide_ibfk_3` FOREIGN KEY (`protein_domain_id`) REFERENCES `protein_domain` (`id`) ON UPDATE CASCADE,
  CONSTRAINT `peptide_ibfk_4` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=3417481 DEFAULT CHARSET=latin1 COLLATE=latin1_bin;
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
  `locus` varchar(30) NOT NULL DEFAULT '',
  `name` varchar(300) NOT NULL DEFAULT '',
  `date` datetime NOT NULL,
  `remove` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`id`),
  KEY `species_id` (`species_id`),
  KEY `sequence` (`sequence`(20)),
  CONSTRAINT `protein_ibfk_1` FOREIGN KEY (`species_id`) REFERENCES `species` (`id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=151050 DEFAULT CHARSET=latin1;
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
  CONSTRAINT `protein_GO_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `protein_GO_ibfk_2` FOREIGN KEY (`GO_id`) REFERENCES `GO` (`id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=1256727 DEFAULT CHARSET=utf8;
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
) ENGINE=InnoDB AUTO_INCREMENT=839261 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=245470 DEFAULT CHARSET=latin1;
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
) ENGINE=InnoDB AUTO_INCREMENT=118718 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_mutations`
--

DROP TABLE IF EXISTS `protein_mutations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_mutations` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `protein_id` int(10) unsigned NOT NULL,
  `acc_id` varchar(30) NOT NULL,
  `mutationType` varchar(25) NOT NULL,
  `location` int(5) unsigned NOT NULL,
  `original` varchar(3) NOT NULL,
  `mutant` varchar(3) NOT NULL,
  `date` datetime NOT NULL,
  `annotation` varchar(256) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `protein_FK` (`protein_id`),
  CONSTRAINT `protein_mutations_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=98209 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_regions`
--

DROP TABLE IF EXISTS `protein_regions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_regions` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `label` varchar(100) NOT NULL,
  `source` enum('predicted','parsed') NOT NULL,
  `start` int(10) unsigned NOT NULL,
  `stop` int(10) unsigned NOT NULL,
  `protein_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  KEY `protein_id` (`protein_id`),
  CONSTRAINT `protein_regions_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=2018 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `protein_scansite`
--

DROP TABLE IF EXISTS `protein_scansite`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `protein_scansite` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `source` varchar(40) NOT NULL DEFAULT 'scansite',
  `value` varchar(20) NOT NULL,
  `score` float NOT NULL,
  `percentile` float NOT NULL,
  `site_pos` int(10) unsigned NOT NULL,
  `protein_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `source` (`source`,`value`,`site_pos`,`protein_id`),
  KEY `protein_id` (`protein_id`),
  CONSTRAINT `protein_scansite_ibfk_1` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=5747591 DEFAULT CHARSET=latin1;
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
  `resource_type` enum('experiment','annotations','dataset') NOT NULL DEFAULT 'experiment',
  `load_type` enum('new','reload','append','extension') NOT NULL,
  `parent_experiment` int(10) unsigned DEFAULT NULL,
  `change_name` text NOT NULL,
  `change_description` text NOT NULL,
  `units` varchar(20) NOT NULL,
  `stage` enum('config','metadata','condition','confirm','complete') NOT NULL DEFAULT 'config',
  `experiment_id` int(10) unsigned DEFAULT NULL,
  `date` datetime NOT NULL,
  PRIMARY KEY (`id`),
  KEY `user_id` (`user_id`,`parent_experiment`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `session_ibfk_1` FOREIGN KEY (`user_id`) REFERENCES `users` (`id`) ON DELETE CASCADE ON UPDATE CASCADE,
  CONSTRAINT `session_ibfk_2` FOREIGN KEY (`experiment_id`) REFERENCES `experiment` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=705 DEFAULT CHARSET=latin1;
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
  `type` enum('hidden','data','stddev','accession','peptide','sites','species','modification','run','none','cluster','numeric','nominative','MS_id') NOT NULL DEFAULT 'none',
  `label` varchar(45) NOT NULL,
  `column_number` int(10) unsigned NOT NULL,
  PRIMARY KEY (`id`),
  KEY `session_id` (`session_id`),
  CONSTRAINT `session_columns_ibfk_1` FOREIGN KEY (`session_id`) REFERENCES `session` (`id`) ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=15032 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `species`
--

DROP TABLE IF EXISTS `species`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `species` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(250) NOT NULL,
  `taxon_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `name` (`name`),
  KEY `taxon_id` (`taxon_id`),
  CONSTRAINT `species_ibfk_2` FOREIGN KEY (`taxon_id`) REFERENCES `taxonomy` (`node_id`) ON UPDATE CASCADE
) ENGINE=InnoDB AUTO_INCREMENT=1700 DEFAULT CHARSET=latin1;
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
  `parent_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`node_id`),
  KEY `kingdom` (`kingdom`),
  KEY `parent_id` (`parent_id`),
  KEY `name_index` (`name`),
  CONSTRAINT `taxonomy_ibfk_1` FOREIGN KEY (`parent_id`) REFERENCES `taxonomy` (`node_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `uniprot_swissprot`
--

DROP TABLE IF EXISTS `uniprot_swissprot`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `uniprot_swissprot` (
  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `accession` varchar(20) NOT NULL,
  `locus` varchar(30) NOT NULL,
  `name` varchar(150) NOT NULL,
  `species` varchar(100) NOT NULL,
  `sequence` text NOT NULL,
  PRIMARY KEY (`id`),
  KEY `species_search` (`species`),
  FULLTEXT KEY `seq_search` (`sequence`)
) ENGINE=MyISAM AUTO_INCREMENT=571965 DEFAULT CHARSET=latin1;
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
  `access_level` enum('reviewer','researcher') NOT NULL DEFAULT 'researcher',
  `expiration` datetime DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `username` (`username`,`email`)
) ENGINE=InnoDB AUTO_INCREMENT=2485 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2013-06-14 15:16:45
