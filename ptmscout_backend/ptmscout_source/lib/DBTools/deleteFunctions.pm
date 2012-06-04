use strict;
use warnings;
use DBTools::insertFunctions;


# deleteExperiment($dbh, $expId)
# Deletes all MS specific portions of an experiment from database
# Deletes: All MSIds - including their data and MSPhosphopep entries
# Inputs: $dbh - database handle
#         $expId - id of experiment to delete 
# Kristen Naegle
# March 23, 2008
sub deleteExperiment($$){
    my($dbh, $expId) = @_;

    # Leave proteinsalone

    # GET all MS_ids in an experiment
    my $MS_ids = returnMSIdsForExpId($dbh, $expId);
    if(scalar(@$MS_ids) > 1 && $MS_ids->[0] != -1){
	my $dataDel = deleteDataByMSId($dbh, $MS_ids);
	my $MS_Pepdel = deleteMSPhosphopepByMSId($dbh, $MS_ids);
	my $MS_del = deleteMSById($dbh, $MS_ids);
    } 
    deleteExpById($dbh, $expId);
    
}

# $count = deleteDataByMSId($dbh, \@MSIds)
# Deletes from data table MS_id's
# Inputs: $dbh - database handle
#         \@MSIds - array of MS Ids to delete
# Outputs: $count - the number deleted
# Kristen Naegle
# March 23, 2008
sub deleteDataByMSId($$){
    my($dbh, $MSIds) = @_;
    my $sth = $dbh->prepare('DELETE FROM data where MS_id=?');
    my $count = 0;
    foreach my $MS_id (@$MSIds){
	$sth->execute($MS_id);
	$count +=1;
    }
    return $count;
}

# $count = deleteMSPhosphopepByMSId($dbh, \@MSIds)
# Deletes from MS_phosphopep table MS_id's
# Inputs: $dbh - database handle
#         \@MSIds - array of MS Ids to delete
# Outputs: $count - the number deleted
# Kristen Naegle
# March 23, 2008
sub deleteMSPhosphopepByMSId($$){
    my($dbh, $MSIds) = @_;
    my $sth = $dbh->prepare('DELETE FROM MS_phosphopep where MS_id=?');
    my $count = 0;
    foreach my $MS_id (@$MSIds){
	$sth->execute($MS_id);
	$count +=1;
    }
    return $count;
}

# $count = deleteMSById($dbh, \@MSIds)
# Deletes from MS by id 
# Inputs: $dbh - database handle
#         \@MSIds - array of MS Ids to delete
# Outputs: $count - the number deleted
# Kristen Naegle
# March 23, 2008
sub deleteMSById($$){
    my($dbh, $MSIds) = @_;
    my $sth = $dbh->prepare('DELETE FROM MS where id=?');
    my $count = 0;
    foreach my $MS_id (@$MSIds){
	$sth->execute($MS_id);
	$count +=1;
    }
    return $count;
}

# deleteExpById($dbh, $ExpId)
# Deletes from experiment by id
# Inputs: $dbh - database handle
#         $expId - single experiment id to delete
# Kristen Naegle
# March 23, 2008
sub deleteExpById($$){
    my($dbh, $expId) = @_;
    my $sth = $dbh->prepare('DELETE FROM experiment where id=?');
    $sth->execute($expId);
}



# deleteProteinById($dbh, $proteinId)
# Deletes protein from protein table
# Inputs: $dbh - database handle
#         $proteinId - id to protein table
# Output: $error - returns error if domains or accessions weren't deleted first
# Kristen Naegle
# March 30, 2008
sub deleteProteinById($$){
    my ($dbh, $proteinId) = @_;
    my $error = 0;
    
    # check for accessions
    my $accArr = returnAccValuesByProteinId($dbh, $proteinId);
    if(scalar(@$accArr) == 1 && $accArr->[0] eq '-1'){
	my $domainIds = returnDomainIdsForProteinId($dbh, $proteinId);
	if(scalar(@$domainIds) == 1 && $domainIds->[0] eq '-1'){
	    my $sth = $dbh->prepare('DELETE FROM protein where id=?');
	    $sth->execute($proteinId);
	 }
	else{
	    $error = 1;
	    handleError('deleteProteinById', 'Attempt to delete protein with existing domains failed', \@_);
	}
    }
    else{
	$error = 1;
	handleError('deleteProteinById', 'Attempt to delete protein with existing accessions failed', \@_);
    }
    
    return $error;
}


# $error = deleteAccByProteinId($dbh, $proteinId);
# Deletes a single accession based on protein_id NEED TO CHANGE THIS
# Inputs: $dbh -database handle
#         $proteinId - (FK) protein Id
# Outputs: nothing
# Kristen Naegle
# March, 28, 2008
sub deleteAccByProteinId($$){
    my ($dbh, $proteinId) = @_;
   
    # check that other accessions exist
    my $sth = $dbh->prepare('Select id from acc where protein_id=?');
    $sth->execute($proteinId);
    my $accIds = returnArrayOfResultsOnCol($sth,0);
    $sth->finish();
    $sth = $dbh->prepare('delete FROM acc where id=?');
    foreach my $accId (@$accIds){
	$sth->execute($accId);
    }
    $sth->finish();
}

# deleteAccByAccId($dbh, $accId);
# Deletes a single accession based on acc.id
# Inputs: $dbh -database handle
#         $acc Id - id of accession table
# Outputs: nothing
# Kristen Naegle
# June 12, 2009
sub deleteAccByAccId($$){
    my ($dbh, $accId) = @_;
   
    
    my $sth = $dbh->prepare('delete FROM acc where id=?');
    $sth->execute($accId);
    
    $sth->finish();
}

# $error = deleteDomainsByProteinId($$)
# Delete all domains that point to protein.id
# Inputs: $dbh - database handle
#         $proteinId - (FK) id to protein table
# Kristen Naegle
# March 28, 2008
sub deleteDomainsByProteinId($$){
    my ($dbh, $proteinId) = @_;
   
    # check that other accessions exist
    my $sth = $dbh->prepare('Select id from domain where protein_id=?');
    $sth->execute($proteinId);
    my $domainIds = returnArrayOfResultsOnCol($sth,0);
    $sth->finish();
    $sth = $dbh->prepare('delete FROM domain where id=?');
    foreach my $id (@$domainIds){
	
	$sth->execute($id);
    }
    $sth->finish();
}

# $error = deleteDomain($dbh, $domainId)
# Delete all domains that point to protein.id
# Inputs: $dbh - database handle
#         $domainId - id of domain table 
# Kristen Naegle
# June 6, 2009
sub deleteDomain($$){
    my ($dbh, $domainId) = @_;
    my $error = 0;
    # check that other accessions exist
    my $sth = $dbh->prepare('delete FROM domain where id=?');
    $sth->execute($domainId);
    $sth->finish();
    return $error;
}

# deleteProteinExpressionByProteinId($dbh, $proteinId)
# Removes entries in protein_expression tables linking to protein_id
# Inputs: $dbh - database handle
#         $proteinId - id to protein table, and entry in protein_expression
# Kristen Naegle
# January 9, 2009
sub deleteProteinExpressionByProteinId($$){
    my ($dbh, $proteinId) = @_;
    my $sth = $dbh->prepare('delete FROM protein_expression where protein_id=?');
    $sth->execute($proteinId);
    $sth->finish();

}

# deleteProteinGOByProteinId($dbh, $proteinId);
# Removes protein_GO entry based on protein id
# Inputs: $dbh - database handle
#         $proteinId - protein_id (FK to protein table);
# January 16, 2009
# Kristen Naegle
sub deleteProteinGOByProteinId($$){
    my ($dbh, $proteinId) = @_;
    my $sth=$dbh->prepare('delete from protein_GO where protein_id=?');
    $sth->execute($proteinId);
    $sth->finish();
}

# deleteAmbiguityByProteinId($dbh, $proteinId);
# Removes an ambiguity entry when deleting by protein Id
# Inputs: $dbh - database handle
#         $proteinId - id to protein table, FK in ambiguity
# Kristen Naegle
# January 9, 2009
sub deleteAmbiguityByProteinId($$){
    my ($dbh, $proteinId) = @_;
    my $sth = $dbh->prepare('delete FROM ambiguity where protein_id=?');
    $sth->execute($proteinId);
    $sth->finish();

}

# $error = deleteProteinRecordsByProteinId($dbh, $proteinId);
# Remove all records linked to protein and the protein.  DOES not allow if protein is linked to a phosphopep
# Records to delete: acc, protein_GO (not yet), domains, protein_expression and protein
# Inputs: $dbh - database handle
#         $proteinId - key to protein table to delete
# Outputs: $error- 0 if no error based on acc/domain existence and 1 otherwise.
# Kristen Naegle
# January 9, 2009
sub deleteProteinRecordsByProteinId($$){
    my ($dbh, $proteinId) = @_;
    deleteAccByProteinId($dbh, $proteinId);
    deleteDomainsByProteinId($dbh, $proteinId);
    deleteProteinExpressionByProteinId($dbh, $proteinId);
    deleteAmbiguityByProteinId($dbh, $proteinId);
    my $error = deleteProteinById($dbh, $proteinId);    
    return $error;
}

# deletePhosphopepPredByPepId($dbh, $phosphopepId)
# Deletes the phosphopep prediction based on the phosphopep_id (FK)
# Inputs: $dbh - database handle
#         $phosphopepId - id in phosphopep and FK to phosphopep_prediction
# Kristen Naegle
# January 16, 2009
sub deletePhosphopepPredByPepId($$){
    my ($dbh, $phosphopepId) = @_;

    my $sth = $dbh->prepare('DELETE from phosphopep_prediction where phosphopep_id=?');
    $sth->execute($phosphopepId);
    $sth->finish();

}

# deletePhosphopepByPepId($dbh, $phosphopepId);
# Deletes the phosphopep based on the id
# Inputs: $dbh - database handle
#         $phosphopepId - id in phosphopep 
# Kristen Naegle
# January 16, 2009
sub deletePhosphopepByPepId($$){
    my($dbh, $phosphopepId) = @_;
    my $sth = $dbh->prepare('DELETE from phosphopep where id=?');
    $sth->execute($phosphopepId);
    $sth->finish();
    

}

# deleteALLGO($dbh)
# Removes all protein_GO and GO entries
# Inputs: $dbh -database handle
# Kristen Naegle
# August 14, 2009
sub deleteAllGO($){
    my($dbh) = @_;
    my $sth=$dbh->prepare('DELETE FROM protein_GO');
    $sth->execute();
    $sth=$dbh->prepare('DELETE FROM GO');
    $sth->execute();
}

# clearAllNonExpressionTables($dbh);
# This clears all tables except expression tables (it does delete protein_expression linkages).  Use this when you want to clear the update database to prepare fresh loading. 
# Inputs: $dbh - database handle !! This should be your update database, not production!
# Outputs: none
# Kristen Naegle
# November 23, 2009
sub clearAllNonExpressionTables($){
    my($dbh) = @_;

    my @tables = ('acc', 'ambiguity', 'data', 'domain', 'protein_GO', 'phosphopep_prediction', 'MS_phosphopep', 'phosphopep', 'MS', 'protein_expression', 'GO', 'experiment', 'protein');

    my $sth;
    my $optimize;

    for(my $i=0; $i<scalar(@tables); $i++){
	my $table = $tables[$i];
	print "Deleting items from $table\n";
	$sth = $dbh->prepare("delete quick from $table");
	$sth->execute();
	$optimize = $dbh->prepare("optimize table $table");
	$optimize->execute();
    }
   
}

1;
