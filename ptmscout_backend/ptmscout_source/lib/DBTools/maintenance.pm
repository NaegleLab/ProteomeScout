use strict;
use warnings;
use DBI;
use DBTools::insertFunctions;
use DBTools::queryTools;
use DBTools::deleteFunctions;
use DBTools::dbEntry;
use DBTools::deleteFunctions;
use GOTools;
use commonTools;
use ncbiDBTools;
use uniprotTools;
# $proteinIdArr = maintainDomains($dbh)
# Runs pfam for all proteins that do not already have a pfam prediction
# Inputs: $dbh - database handle (domains are committed)
# Outputs: $proteinIdArr - ref. to array of protein ids
# Kristen Naegle
# March 29, 2008
sub maintainDomains($){
    my ($dbh) = @_;
#    my $CPU = returnCPU();
    
    my $domainHash = returnDomainHashFromGlobalFile();
    my $proteinIdArr = returnProteinIdsWithNoDomains($dbh);
    handleDomainsForProteinIds($dbh, $proteinIdArr, $domainHash);
    
    return $proteinIdArr;

}

# updateGO($dbh, $REMOVE_IEA);
# Updates GO, first removes all GO terms, loads new files (MUST MANUALLY LOAD OBO FILE), then maintiansGOTerms for all proteins
# Inputs: $dbh - database handle
#         $REMOVE_IEA - whether to run with IEA removed from GO Annotation files
# Kristen Naegle
# August 14, 2009
sub updateGO($$){
    my ($dbh, $REMOVE_IEA) = @_;
    #delete all GO terms
    deleteAllGO($dbh);
    $dbh->commit();

    retrieveAllGOAnnotationFiles($REMOVE_IEA);

    maintainGOTerms($dbh);
    $dbh->commit();
}
 
# maintainGOTerms($dbh)
# Adds GO terms to all proteins in database
# Inputs: $dbh - database handle
# Outputs: none
# Kristen Naegle
# June 12, 2009
sub maintainGOTerms($){
    my($dbh) = @_;
    
    my $sth = $dbh->prepare('SELECT id from protein');
    $sth->execute();
    my $proteinIds = returnArrayOfResultsOnCol($sth, 0);
   # my $ontologyFile = "/home/knaegle/SVN/knaegle/scripts/GO_DATA/gene_ontology.1_2.obo";
    #my ($rootDir, $rootFaile, $ontologyFile) = returnGlobalGOFileInfo();
    #my $domainHash = returnDomainHashFromGlobalFile();
    

    handleGOTermsForProteinIds($dbh, $proteinIds);
    
    $dbh->commit();
}

#cleanAccessions($)
#Removes duplicate accession/protein ids from database. These should NOT be generated. Only insertAcc in handleAcc if accession does not already exist
# Please do not use
# Kristen Naegle
# March 30, 2008
sub cleanAccessions($){
    my ($dbh) = @_;
    
    # find proteins with duplicate accessions
  #  my $sth = $dbh->prepare('select protein.id from acc, protein where acc.protein_id=protein.id and acc.value in ( select value from acc group by value having count(*) > 1)');
    my $sth = $dbh->prepare('select protein_id from acc where type=\'swissprot\' group by protein_id having count(*) > 1');
    $sth->execute();
    my $duplicates = returnArrayOfResultsOnCol($sth, 0);
    #my @proteinIdArr; 
    print "Proteins with duplicates\n";
    foreach my $duplicateProtein (@$duplicates){
	# get the protein id. delete accession and then protein ids
	print $duplicateProtein."\n";
	# delete domains, delete protein, delete accessions
	deleteAccByProteinId($dbh, $duplicateProtein);
 	deleteDomainsByProteinId($dbh, $duplicateProtein);
 	deleteProteinById($dbh, $duplicateProtein);
    }
    print 'There were'.scalar(@$duplicates)."\n";
}

# (\@proteinsRemoved) = removeOrphanedProteins($dbh)
# This removes proteins that aren't linked to phosphopep or ambiguity DOES NOT Commit
# Inputs: $dbh - database handle
# Outputs: $proteinsRemoved - ref to array of protein ids that were removed
# Kristen Naegle
# January 9, 2009
sub removeOrphanedProteins($){
    my ($dbh) = @_;

    # get the proteins that aren't linked to ambiguity, MS or phosphopep
    my $sth = $dbh->prepare('SELECT protein.id from protein left outer join phosphopep on protein.id=phosphopep.protein_id left outer join ambiguity on protein.id=ambiguity.protein_id left outer join MS on MS.protein_id=protein.id where phosphopep.protein_id IS NULL and ambiguity.protein_id is NULL and MS.protein_id is NULL');
    $sth->execute();
    my $count = 0;
    my @proteinsRemoved;
    while(my @row = $sth->fetchrow_array){
	$count += 1;
	my $proteinId = $row[0];

	deleteAccByProteinId($dbh, $proteinId);
	deleteDomainsByProteinId($dbh, $proteinId);
	deleteProteinExpressionByProteinId($dbh, $proteinId);
	deleteProteinGOByProteinId($dbh, $proteinId);
	my $error = deleteProteinById($dbh, $proteinId);
	if(!$error){
	    push @proteinsRemoved, $proteinId;
	}

	

    }
   
    return (\@proteinsRemoved);
}

# (\@peptidesRemoved) = removeOrphanedPhosphopeps($dbh)
# This removes phosphopeps that aren't linked to an experiment (note may have to call removeOrphanedProteins after this as well)
# Inputs: $dbh - database handle
# Outputs: $peptidesRemoved - ref to array of phosphopep ids that were removed
# Kristen Naegle
# January 16, 2009
sub removeOrphanedPhosphopeps($){
    my ($dbh) = @_;
    my $sth = $dbh->prepare('SELECT phosphopep.id from phosphopep left outer join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_phosphopep.phosphopep_id IS NULL');
    $sth->execute();
    my @removed;
    
    
    while(my @row = $sth->fetchrow_array){
	#dependencies are phosphopep_predictions
	my $phosphopep = $row[0];
	deletePhosphopepPredByPepId($dbh, $phosphopep);
	deletePhosphopepByPepId($dbh, $phosphopep);
	push @removed, $phosphopep;

    }
    return \@removed;

}

# updatePfamSitesForPhosphopeps($dbh, $phosphopepIds)
# Given a list of phosphopep Ids, update the pfam_site field. To be run when domain predictions are updated. 
# Inputs: $dbh - database handle
#         $phosphopepIds - ids of phosphopep table
# Kristen Naegle
# June 29, 2009
sub updatePfamSitesForPhosphopeps($$){
    my ($dbh, $phosphopepIds) = @_;

    foreach my $phosphopepId (@$phosphopepIds){
	# get the protein id and the position from phosphopep table
	my $sth = $dbh->prepare('SELECT * from phosphopep where id=?');
	$sth->execute($phosphopepId);
	
	my $results = returnMultipleResultsForExSTH($sth);
	my $result = $results->[0]; #expect single result since id is unique
	my $proteinId = $result->{'protein_id'};
	my $pos = $result->{'site_pos'};
	my $pfamSite = returnPfamSite($dbh, $proteinId, $pos);
	my $updateSTH = $dbh->prepare('UPDATE phosphopep set pfam_site=? where phosphopep.id=?');
	$updateSTH->execute($pfamSite, $phosphopepId);
	
    }


}

#  maintainAmbiguousExperiments($dbh, $rdbh, $LOCAL, $REFSEQ)
# Get all the unique peptide/protein ids from the database (dbh) that belong to experiments with ambiguous loaded.  Then based on whether just LOCAL or REFSEQ or both are ones than force reloading from those databases
# Inputs: $dbh - database handle 
#         $rdbh - a refseq database handle 
#         $LOCAL - 1 if you want to load based on overlap with members in $dbh
#         $REFSEQ - 1 if you want to laod based on overlap with refseq 
# Kristen Naegle
# July 2, 2009
sub maintainAmbiguousExperiments($$$$){
    my($dbh, $rdbh, $LOCAL, $REFSEQ) = @_;
    
    my $FORCE = 1;
    
    # get all experiments that have ambiguity =1 and then run with force
    my $sth = $dbh->prepare('SELECT id from experiment where ambiguity=1');
    $sth->execute();
    my $ids=returnArrayOfResultsOnCol($sth,0);
    foreach my $id (@$ids){
	print "---------INSERTING AMBIGUITY FOR EXPERIMENT ID: $id-------------\n";
	print "DEBUG: In maintain: FORCE: $FORCE\tLOCAL: $LOCAL\tREFSEQ: $REFSEQ\n";
	insertAmbiguousPeptidesForExpId($dbh, $rdbh, $id, $LOCAL, $REFSEQ, $FORCE);
	$dbh->commit();
	
    }
   

}

# Given dbh handle of an update database (which starts out as copy of current database), completely reload experiments and update all third party information).
# 1. Creates a copy of all experiment loading statements
# 2. Clears all tables (that aren't expression tables).  Optimizes tables
# 3. Updates GO files
# 4. Update Pfam files
# 5. Updates Refseq database (MAKE SURE TO DO THIS off-line)
# 6. Re-loads cleaned experiments
# 7. Also off-line - load updated experiments for phospho.ELM, phosphoSite and uniprotKB searches
# Inputs: $dbh - database handle - for update database (don't do production)
# Kristen Naegle
# Nov. 23, 2009
sub createNewVersion($){
    my ($dbh) = @_;

    #1. First get all experiment loading statements. Place them in a globalVar directory for later reloading.
    print "--------Making backup of load statements--------\n";
    my $f = makeBackupOfLoadStatements($dbh); 
    
    #2. Prompt user to continue..make sure they say yes to correct database handle, to updated RefSeq
    print "Have you changed replaced phosphoELM and phosphoSite symbolic links if applicable?";
    my $alteredDB = <STDIN>;
    chomp $alteredDB;
    if($alteredDB ne 'yes'){
	print "Do it then!\n";
	exit;
    }

    print "--------Updating Uniprot Files----\n";
    retrieveAndParseUniprotFileSearches();
    
    print "--------UPDATING REFSEQ--------\n";
    my $CHECK = 0;
    my $files = updateRefSeq($CHECK);
    print "ADDED file numbers @$files\n";
    
    #3. Wipe all tables
    print "------Wiping all tables from database, except expression related-----\n";
   clearAllNonExpressionTables($dbh);
 $dbh->commit();

    #4. Update GO files
    print "------Updating Gene Onotlogy------\n";
    my $REMOVE_IEA = 1;
    updateGO($dbh, $REMOVE_IEA);
     $dbh->commit(); #however, once cleared, this should really only be updating files.

    #5. Update pfam files. -- SAVING this for after Bioperl handles parsing of new HMMER3 format


    #6. Walk through the file of datasets to load
    print "-------RELOADING experiments------\n";
    loadMultipleExperimentFile($f);
}

# loadMultipleExperimentFile($dbh, $file);
# Given a file of command line statementst to call_loadDataFile, load each statement
# Inputs:$dbh - database handle 
#        $file - line by line experiment loading statements. Like that produced in makeBackupOfLoadStatements
# Kristen Naegle
# Nov. 24, 2009
sub loadMultipleExperimentFile($$){
    my ($dbh, $file) = @_;

    open(FI, $file) || die "Can't open $file for reading\n";
    while(defined(my $line = <FI>)){
	#print $line;
	my @line; 

	while ( $line =~ m/\'(.+?)\'|\s*([^\']+?)\s/ ) {
	    
	    my $match = defined($1) ? $1 : $2;
    
	    $line = substr($line, length($`) + length($&));
    
	    #print "\t$match\n";
	    push @line, $match;
    
	} 

	my($command, $script, $dataFile, $name, $author, $description, $contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_link, $submissionEmail, $primaryMods, $journal, $pub_date, $volume, $pages) = @line;# = split(' ', $line);
			
	print "DEBUG: $dataFile\nname:$name\nauthor:$author\ndescription:$description\ncontact:$contact\nPMID:$PMID\nURL:$URL\npublished:$published\nambiguity:$AMBIGUITY\nEXPORT:$EXPORT\n$expId_link\nsubmissionEmail:$submissionEmail\nprimareyMods:$primaryMods\n";


	my ($expId, $accepted, $rejected) = call_loadDataFile($dbh, $dataFile,$name, $author, $description, $contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_link, $submissionEmail, $primaryMods, $journal, $pub_date, $volume, $pages);

	print "Experiment ID: $expId\n";
	print "Accepted: $accepted Rejected:$rejected\n";	

    }

    # only in end maintain entire database
    maintainDBAfterFileLoad($dbh);

    close(FI);

}

# makeBackupOfLoadStatements();
# Makes a back up of all load statements -time stamped and updates the symbolic link in the backpu path loadStatements to reflect latest backup
# Kristen Naegle
# Nov. 23, 2009
sub makeBackupOfLoadStatements($){
    my ($dbh) = @_;
    my $backupDir = $globalVars::BACKUP_PATH;
    my $date = returnTimeStr();
    my @date = split(' ', $date); #get rid of hours minutes, etc.
    $date = $date[0];
    my $f = $backupDir.'loadStatements';
    my $file = $f.".".$date;
    writeLoadsForAllExperiments($dbh, $file);
    `ln -sf $file $f`;
    return $f;
}

# writeLoadsForAllExperiments($dbh, $file);
# For all experiments in database, prints load statements into desired file
# Inputs: $dbh - database handle
#         $file - target file (overwritten) for load statements
# Kristen Naegle
# November 23, 2009
sub writeLoadsForAllExperiments($$){
    my ($dbh, $file) = @_;

    # get all experiment ids
    my $sth = $dbh->prepare('select experiment.id, count(*) as NUM from experiment join MS on MS.experiment_id=experiment.id group by experiment.id order by NUM desc');
    $sth->execute();
    my $expIds = returnArrayOfResultsOnCol($sth,0);
    open(F, ">$file") || die "Can't open $file for writing\n";
    foreach my $expId (@$expIds){
	my $loadStr = outputDataLoadStatementForExpId($dbh, $expId);
	print F $loadStr."\n";

    }
    close(F);
}


# $files = updateRefSeq($CHECK);
# This automatically deletes refseq table, parses refseq directory, and loads all current protein files for mammalian vertebrates.
# Inputs: $CHECK - boolean value to check for gi before load
# Outputs: $files- ref to array of file numbers loaded
# Kristen Naegle
# Nov 24, 2009
sub updateRefSeq($){
    my ($CHECK) = @_;
    my $rdbh = returnRefSeqDBHNOCommit();
    

    # first have to clear refseq 
    clearRefSeqDB($rdbh);
    $rdbh->commit();

    #Get files to load
    my $files = addAllFilesToDB($rdbh, $CHECK);

    return $files;
}

# Global log file tells us what processes were started, what time and when they finished.  If the last process doesn't have a finished time, then we'll run maitenance, like domains, GO Terms etc.  Run this at the start of each load data file.
# Inputs: $dbh -database handle
# Outputs: $MAINTAINED - boolean value that says whetehr or not maintenance was run
# Kristen Naegle
# Nov 24, 2009
sub checkAndMaintainKilledProcess($){
    my($dbh) = @_;
    
    my $KILLED = checkForKilledProcess();

    if($KILLED){
	writeProcess('checkAndMaintainKilledProcess', \@_);

# RUN MAINTENANCE STUFF
	print "-----Maintaining Database after a killed process-----\n";
# # 	print "Maintaining local ambiguity links\n";
# # 	my $LOCAL = 1;
# # 	my $REFSEQ = 1;
# # 	my $rdbh = returnRefSeqDBHNOCommit();
# # 	maintainAmbiguousExperiments($dbh, $rdbh, $LOCAL, $REFSEQ);
# # 	$dbh->commit();
# # 	$rdbh->rollback();
	
# # 	print "Maintaining expression linkages\n";
# # 	my ($numHandled, $numFailed) = maintainExpressionLinkages($dbh);
# # 	print "Handled: $numHandled\tFAILED: $numFailed\n";
# # 	$dbh->commit();
	
	print "Maintaining domains in case of failure\n";
	maintainDomains($dbh);
	
	# Expensive..but probably necessary to make sure all goes ok
	print "Maintaining GO Terms\n";
	maintainGOTerms($dbh);

	writeProcessFinish();
    }
}

# $predIds = maintainPELMKinaseAnnotations($dbh);
# Calls handlePELMKinaseAnnotations  wiht the global symbolically linked pelm file
# Inputs: $dbh - database handle
# Outputs: $predIs - ref. to array of prediction ids inserted
# Kristen Naegle
# Nov. 24, 2009
sub maintainPELMKinaseAnnotations($){
    my ($dbh) = @_;
    my $pelmRoot = $globalVars::DATA_PATH;
    my $pelmFile = $pelmRoot.'phosphoELM';
    my $predIds = handlePELMKinaseAnnotations($dbh, $pelmFile);
    return $predIds;


}

1;

