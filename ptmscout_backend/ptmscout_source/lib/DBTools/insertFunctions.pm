use strict;
use warnings;
use DBI;
use DBTools::resultParsing;
use DBTools::queryTools;
use errorHandling;

# $expId = insertExperiment($dbh, $name, $author, $date, $description, $contact, $PMID, $URL, $published, $Ambiguity, $export, $expId_link, $dataFileLocation);
# Inserts experiment into experiment table
# Inputs: $dbh - database handler
#         $fileName - name of experiment
#         $author - author of experiment
#         $date - date of experiment 
#         $description - additional info about experiment
#         $contact - contact info for dataset
#         $PMID - pubmed id if published
#         $URL - URL if exists to dataset/paper
#         $published - boolean value - if data published this is set true.
#         $ambiguity - boolean value if ambiguity should/was loaded for fragments
#         $export - boolean value to say whether a dataset can be exported
#         $expId_Link - expId if this is a modification of an existing experiment
#         $dataFileLocation - location of data file 
#         $submissionEmail - email of person submitting dataset
#         $primaryMods - string of primary modifications (e.g. 'S,T' or 'Y' or 'S,T,Y', 'K')
# Outputs: $id - the id to that experiment
# Kristen Naegle
# March 3, 2008
# Mod. June 24, 2009 - to include ambiguity bit
sub insertExperiment($$$$$$$$$$$$$$$$$$$){
    my ($dbh, $fileName, $author, $date, $description, $contact, $PMID, $URL, $published, $ambiguity, $export, $expId_link, $dataFileLocation, $submissionEmail, $primaryMods, $journal, $pub_date, $volume, $pages) = @_;
    #check to see if that experiment exists
    my $expId = returnExperimentId($dbh, $fileName);
    if($expId < 0){ #it doesn't already exist so insert and return latest id. 
	my $sth = $dbh->prepare('INSERT INTO experiment (name,author, date,description, contact, PMID, URL, published, ambiguity, export, experiment_id, dataset, submitter, primaryModification, volume, pages, journal, pub_date)  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?, ?,?,?,?)'); 

	$sth->execute($fileName, $author, $date, $description, $contact, $PMID, $URL, $published, $ambiguity, $export, $expId_link, $dataFileLocation, $submissionEmail, $primaryMods, $volume, $pages, $journal, $pub_date) or die "Couldn't execute statement: " . $sth->errstr;
	my $lastID_sth = $dbh->prepare('SELECT LAST_INSERT_ID()');
	$lastID_sth->execute or die "Couldn't execute statement: " . $sth->errstr;
	my $lastIDResults = returnArrayOfResultsOnCol($lastID_sth, 0);
	$expId = $lastIDResults->[0];
	$sth->finish;
	$lastID_sth->finish;
    }
    else{ #update the fields that might have changed. 
	my $sth = $dbh->prepare('UPDATE experiment set PMID=?, URL=?, published=?, ambiguity=?, export=?, dataset=?, submitter=? where id=?');
	$sth->execute($PMID, $URL, $published, $ambiguity, $export, $dataFileLocation, $submissionEmail, $expId);
    }
    return $expId;
}


# $id = insertProtein($dbh, $sequence, $species, $acc_gene, $name, $date)
# Inserts protein into protein table
# Inputs: $dbh - database handler
#         $sequence - protein sequence
#         $species - common name of species
#         $acc_gene - entrez gene name
#         $name - long name
#         $date - date of entry Month-Year (02-2008)
# Outputs: $id - the id to that protein
# Kristen Naegle
# March 3, 2008, 
# Modified March 26, 2008 to include date
sub insertProtein($$$$$$){
    my ($dbh, $sequence, $species, $acc_gene, $name, $date) = @_;
    my $proteinId;
    my $sth = $dbh->prepare('INSERT INTO protein (sequence ,species, acc_gene,name, date)  VALUES (?,?,?,?,?)'); 

    $sth->execute($sequence, $species, $acc_gene, $name, $date) or die "Couldn't execute statement: " . $sth->errstr;
    $proteinId = getLastInsertId($dbh);
    $sth->finish;
    return $proteinId;
}

# $id = insertAcc($dbh, $type, $value, $proteinId);
# Inserts accession into accession table
# Inputs: $dbh - database handler
#         $type - type of accession ('gi', 'swissprot', 'uniprot', 'david', 'entrez_protein', 'gene_synonym', 'undefined')
#         $value - accession 
#         $proteinId - (FK) key to protein table id 
# Outputs: $id - the id to that accession
# Kristen Naegle
# March 3, 2008
# Updated Feb 20, 2008 to make sure record doesn't exist before inserting.
# Updated Dec. 4, 2009 to add OUT OF DATE field to indicate an old accession
sub insertAcc($$$$$){
    my ($dbh, $type, $value, $proteinId, $OUT_OF_DATE) = @_;
    my $accId;
    
    #ensure this isn't a duplicate entry
    my %argsAcc;
    my $strType = '"'.$type.'"';
    my $charAllow = 45;
    if(length($value) > $charAllow){
	$value = substr($value, 0, $charAllow);
    }
    my $strValue = '"'.$value.'"';
    $argsAcc{'type'} = $strType;
    $argsAcc{'value'} = $strValue; #cut this off to 45 characters
    $argsAcc{'protein_id'} = $proteinId;
    my $accIdExist = checkForExistence($dbh, '*', 'acc', \%argsAcc);
    
    if(!$accIdExist){
	my $sth = $dbh->prepare('INSERT INTO acc (type, value, protein_id, out_of_date)  VALUES (?,?,?,?)'); 
	# check for FK constraing
	my %argsPep;
	$argsPep{'id'} = $proteinId;
	my $proteinIdExist = checkForExistence($dbh, '*', 'protein', \%argsPep);
	if($proteinIdExist){

	    $sth->execute($type, $value, $proteinId, $OUT_OF_DATE) or die "Couldn't execute statement: " . $sth->errstr;
	    $accId = getLastInsertId($dbh);
	    
	}
	else{
	    $accId = -1;
	    handleError('insertAcc', 'Attempted acc insert where protein_id FK fails\n', \@_);
	}
	$sth->finish;
    }
    else{
	my $sth = $dbh->prepare('SELECT id from acc where type=? and value=? and protein_id=?');
	$sth->execute($type, $value, $proteinId);
	$accId = returnSingleResultOnCol($sth, 0);
    }


    return $accId;
}

# $id = getLastInsertId($dbh)
# returns the id from last insert
# Inputs: $dbh - database handler
# Outputs: $id - id from last insert
# Kristen Naegle
# March 3, 2008
sub getLastInsertId($){
    my $dbh = shift;
    my $lastId;
    my $lastID_sth = $dbh->prepare('SELECT LAST_INSERT_ID()');
    $lastID_sth->execute or die "Couldn't execute statement: " . $lastID_sth->errstr;
    my @row = $lastID_sth->fetchrow_array;
    $lastId = $row[0];
    $lastID_sth->finish;
    return $lastId;
}

# $id = insertDomain($dbh, $label, $start, $stop, $p_value, $source, $params, $protein_id, $version);
# Inserts domain into domain table
# Inputs: $dbh - database handler
#         $label - name of domain
#         $start - starting position of domain
#         $stop - end position of domain
#         $p_value - significance of domain
#         $source - source of prediction ('pfam')
#         $params - params used in that prediction (e.g. pfam_fs, pValueCutoff=x)
#         $proteinId - (FK) key to protein table id 
#         $version - version release number of the domain tools used
# Outputs: $id - the id to that domain
# Kristen Naegle
# March 3, 2008
sub insertDomain($$$$$$$$$){
    my ($dbh, $label, $start, $stop, $p_value, $source, $params, $protein_id, $version) = @_;
    my $domainId;
    my $sth = $dbh->prepare('INSERT INTO domain (label, start, stop, p_value, source, params, protein_id, version) values (?, ?, ?, ?, ?, ?, ?, ?)');

    # check for FK constraing
    my %argsPep;
    $argsPep{'id'} = $protein_id;
    my $proteinIdExist = checkForExistence($dbh, '*', 'protein', \%argsPep);
    if($proteinIdExist){
	$sth->execute($label, $start, $stop, $p_value, $source, $params, $protein_id, $version) or die "Couldn't execute statement: " . $sth->errstr;
	$domainId = getLastInsertId($dbh);
    }
    else{
	handleError('insertDomain', 'Attempt to insert domain where FK protein_id fails\n', \@_);
	$domainId = -1;
    }
    $sth->finish;
    return $domainId;
}


# \@domainIds = insertDomainArr($dbh, \@label, \@start, \@stop, \@p_value, $source, $params, $protein_id, $version);
# # Inserts multiple domains into domain table for a single protein (and assumes the source and params are the same
# Inputs: $dbh - database handler
#         \@label - name of domain
#         \@start - starting position of domain
#         \@stop - end position of domain
#         \@p_value - significance of domain
#         $source - source of prediction ('pfam')
#         $params - params used in that prediction (e.g. pfam_fs, pValueCutoff=x)
#         $proteinId - (FK) key to protein table id 
#         $version - version release number of domain tools used
# Outputs: \@ids - the ids to those domain
# Kristen Naegle
# March 3, 2008
sub insertDomainArr($$$$$$$$$){
    my ($dbh, $label, $start, $stop, $p_value, $source, $params, $protein_id, $version) = @_;
    my @domainId;
    my $sth = $dbh->prepare('INSERT INTO domain (label, start, stop, p_value, source, params, protein_id, version) values (?, ?, ?, ?, ?, ?, ?, ?)');
    for (my $i=0; $i < scalar(@$label); $i++){
	$sth->execute($label->[$i], $start->[$i], $stop->[$i], $p_value->[$i], $source, $params, $protein_id, $version) or die "Couldn't execute statement: " . $sth->errstr;
	push @domainId, getLastInsertId($dbh);
    }
    $sth->finish;
    return \@domainId;
}

# ($MS_phosphopepId, $phosphopepId) = insertPhosphopep($dbh, $pep_tryps, $pep_aligned, $site_pos, $site_type, $pfam_site, $protein_id, $MS_id);
#  Inserts phosphopep into phosphopep table AND MS_phosphopep table
#  Inputs: $dbh - database handler
#          $pep_tryps - idealized tryps phosphopep
#          $pep_aligned - aligned peptide +/-7aa from protein sequence
#          $site_pos - position of site in sequence
#          $site_type - type of phosphorylation (S/T/Y);
#          $pfam_site - label of domain if this site falls in domain (~~~ else)
#          $protein_id - (FK) id into protein table from where alignment came from
#          $MS_id - id of MS line in experiment (to link to in MSPhosphopep table)
#  Outputs: $MS_phosphopepId - id to the MS_phosphopep table entry
#           $phosphopepId - id to the phosphopep table entry
# Kristen Naegle
# March 10, 2008
sub insertPhosphopep($$$$$$$$){
    my ($dbh, $pep_tryps, $pep_aligned, $site_pos, $site_type, $pfam_site, $protein_id, $MS_id) = @_;
    my $phosphopepId;
    my $MS_phosphopepId;
    my $sth = $dbh->prepare('INSERT INTO phosphopep (pep_tryps, pep_aligned, site_pos, site_type, pfam_site, protein_id) values (?, ?, ?, ?, ?, ?)');
    $sth->execute($pep_tryps, $pep_aligned, $site_pos, $site_type, $pfam_site, $protein_id) or die "Couldn't execute statement: " . $sth->errstr;
    $phosphopepId = getLastInsertId($dbh);
   # print "In insertPhosphopep: $phosphopepId\n";
    $MS_phosphopepId = insertMSPhosphopep($dbh, $MS_id, $phosphopepId);
    $sth->finish;
    return ($MS_phosphopepId, $phosphopepId);

}

# $id = insertMSPhosphopep($dbh, $MS_id, $phosphopepId)
# Inserts MS Phosphopep into MS_Phosphopep table
# Inputs: $dbh - database handler
#         $MS_id - (FK) id to MS table 
#         $phosphopepId - (FK) id to phosphopep table
# Outputs: $id - id to MS_phosphopep table
# Kristen Naegle
# March 10, 2008
sub insertMSPhosphopep($$$){
    my ($dbh, $MS_id, $phosphopepId) = @_;
    my $MS_phosphopepId;
    
    # Check for phosphopep and MS id first since both are FK? 
    my %argsPep;
    $argsPep{'id'} = $phosphopepId;
    my $pepIdExist = checkForExistence($dbh, '*', 'phosphopep', \%argsPep);
    my %argsMS;
    $argsMS{'id'} = $MS_id;
    my $MSIdExist = checkForExistence($dbh, '*', 'MS', \%argsMS);
    # Return MS_phosphopepId if record already exists
    my $sth = $dbh->prepare('SELECT id from MS_phosphopep where MS_id=? and phosphopep_id=?');
    $sth->execute($MS_id, $phosphopepId);
    $MS_phosphopepId = returnSingleResultOnCol($sth, 0);
    if($MS_phosphopepId < 0){
					      
	if($MSIdExist && $pepIdExist){
	    my $sth_ref = $dbh->prepare('INSERT INTO MS_phosphopep (MS_id, phosphopep_id) values (?, ?)');
	    $sth_ref->execute($MS_id, $phosphopepId);
	    $MS_phosphopepId = getLastInsertId($dbh);
	    $sth_ref->finish;
	}
	else{
	    $MS_phosphopepId = -1;
	    my $estr = "FATAL: Foreign Key existence problem, MS_idExist=$MSIdExist and phosphopepIDExist=$pepIdExist";
	    handleError('insertMSPhosphopep', $estr, \@_);
	}
    }
    return $MS_phosphopepId;

}

# $id = insertMSLine($dbh, $phosphopep, $experiment_id, $protein_id);
# Inserts MS Phosphopep into MS_Phosphopep table
# Inputs: $dbh - database handler
#         $phosphopep - the phosphorylated peptide detected in experiment
#         $experiment_id - (FK) id to experiment table
#         $protein_id - (FK) id to protein table
# Outputs: $id - id to MS table
# Kristen Naegle
# March 10, 2008
sub insertMSLine($$$$){
    my ($dbh, $phosphopep, $experiment_id, $protein_id) = @_;
    my $MSId;

    my $sth = $dbh->prepare('INSERT INTO MS (phosphopep, experiment_id, protein_id)  VALUES (?,?,?)'); 

    my %argsExp;
    $argsExp{'id'} = $experiment_id;
    my $expIdExist = checkForExistence($dbh, '*', 'experiment', \%argsExp);
    my %argsProtein;
    $argsProtein{'id'} = $protein_id;
    my $proteinIdExist = checkForExistence($dbh, '*', 'protein', \%argsProtein);
    if($proteinIdExist && $expIdExist){
	$sth->execute($phosphopep, $experiment_id, $protein_id) or die "Couldn't execute statement: " . $sth->errstr;
	$MSId = getLastInsertId($dbh);
	$sth->finish;
    }
    else{
	$MSId = -1;
	my $estr = "FATAL: Foreign Key existence problem, protein_idExist=$proteinIdExist and experimentIDExist=$expIdExist";
	handleError('insertMSLine', $estr, \@_);
    }
    return $MSId;
}

# $id = insertPhosphopepPrediction($dbh, $source, $value, $score, $phosphopep_id)
# Inserts Phosphopep prediction into phosphopep_prediction table
# Inputs: $dbh - database handler
#         $source - source (e.g. scansite or pelm)
#         $value - name of predicted interaction (kinase name, protein_SH2 etc.)
#         $score - score if known, can be null
#         $phosphopep_id - (FK) id to phosphopep table
# Outputs: $id - id to phosphopep_prediction table
# Kristen Naegle
# April 6, 2008
sub insertPhosphopepPrediction($$$$$){
    my ($dbh, $source, $value, $score, $phosphopep_id) = @_;
    my $predictionId;

# Check for phosphopep_prediction, based on phosphopep_id, source and value
    my $sthP = $dbh->prepare('select id from phosphopep_prediction where phosphopep_id=? and source=? and value=?');
    $sthP->execute($phosphopep_id, $source, $value);
    $predictionId = returnSingleResultOnCol($sthP, 0);
    if($predictionId == -1){
    

	my $sth = $dbh->prepare('INSERT INTO phosphopep_prediction (source, value, score, phosphopep_id) VALUES (?, ?, ?, ?)');
	
	my %argsExp;
	$argsExp{'id'} = $phosphopep_id;
	my $phosphopepIdExist = checkForExistence($dbh, '*', 'phosphopep', \%argsExp);
	if($phosphopepIdExist){
	    $sth->execute($source, $value, $score, $phosphopep_id) or die "Couldn't execute statement: " . $sth->errstr;
	    $predictionId = getLastInsertId($dbh);
	    $sth->finish;
	}
	else{
	    $predictionId = -1;
	    my $estr = "FATAL: Foreign Key existence problem, phosphopep_idExist=$phosphopepIdExist";
	    handleError('insertPhosphopepPrediction', $estr, \@_);
	}
    }
    else{
	print "Prediction for $source, $value, and $phosphopep_id already existed\n";
    }
    return $predictionId;
}

# $dataID = insertData($dbh, $type, $run, $label, $priority, $value, $MS_id);
# Inserts data into data table
# Inputs: $dbh - database handle
#         $type - enumerated type: 
#         $run - the run name, average is default
#         $label - the label of that measurement (e.g. 0min)
#         $priority - priority of that measurement
#         $value - the value of the measurement (float)
#         $MS_id - (FK) the id for a measurement 
# Outputs: $id - id of the data inserted
# Kristen Naegle
# March 17, 2008
sub insertData($$$$$$$){
    my($dbh, $type, $run, $label, $priority, $value, $MS_id) = @_;
    my $id;
    my $NA = 0;
   
    #handle enumerated type here?
    my %args;
    $args{'id'} = $MS_id;
    my $MSIdExist = checkForExistence($dbh, '*', 'MS', \%args);
    if($MSIdExist){
	my $sth = $dbh->prepare('INSERT INTO data (type, run, label, priority, value, NA, MS_id) VALUES (?, ?, ?, ?, ?, ?, ?)');
	if($value =~ /[a-z]/i){
	    $NA = 1;
	    my $newValue;
	    $value = $newValue; #null
	}
    
	$sth->execute($type, $run, $label, $priority, $value, $NA, $MS_id) or die "Couldn't execute statement: ".$sth->errstr;
	$id = getLastInsertId($dbh);
	$sth->finish;
    }
    else{
	$id = -1;
	my $estr = "FATAL: Foreign Key existence problem, MSID $MS_id does not exist";
	handleError('insertMSLine', $estr, \@_);
    }
    return $id;

}

# $id = insertGO($dbh, $GO, $aspect, $ontologyVersion)
# Inserts GO term and aspect into GO table
# Inputs: $dbh - database handle
#         $GO - gene ontology term
#         $aspect - aspect value, 'C', 'F', or 'P'
#         $ontologyVersion - version of ontology file this is from
# Outputs: $id - id value to table entry inserted
# Kristen Naegle
# April 15, 2008
# June 23, 2009 - added version 
sub insertGO($$$$$){
    my ($dbh, $GO, $term, $aspect, $version) = @_;
    my $id;
    my $sth = $dbh->prepare('INSERT into GO (GO, term, aspect, version) VALUES (?, ?, ?, ?)');
    $sth->execute($GO, $term, $aspect, $version);
    $id = getLastInsertId($dbh);
    return $id;

}

# $id = insertProteinGO($dbh, $protein_id, $GO_id, $annotationVersion)
# Inserts table entry into protein_GO
# Inputs: $dbh - database handle
#         $protein_id (FK) - protein.id
#         $GO_id (FK) - GO.id
#         $annotationVersion - version of annotaiton this record came from
# Outputs: $id - returns -1 if not successful (FK problems), else returns inserted id
# Kristen Naegle
# April 15, 2008
# June 23, 2009 - added verion
sub insertProteinGO($$$$){
    my ($dbh, $protein_id, $GO_id, $version) = @_;
    my $id;
    my %protArgs;
    $protArgs{'id'} = $protein_id;
    my $protIdExist = checkForExistence($dbh, '*', 'protein', \%protArgs);
    my %goArgs;
    $goArgs{'id'} = $GO_id;
    my $goIdExist = checkForExistence($dbh, '*', 'GO', \%goArgs);
    if($protIdExist && $goIdExist){
	my $sth = $dbh->prepare('INSERT into protein_GO (protein_id, GO_id, version) VALUES (?,?,?)');
	$sth->execute($protein_id, $GO_id, $version);
	$id = getLastInsertId($dbh);
	$sth->finish;
	
    }
    else{
	$id = -1;
	my $estr = "FATAL: Foregin Key existence probelm, protein_id $protein_id exist is $protIdExist and GO_id $GO_id exist is $goIdExist";
	handleError('insertProteinGO', $estr, \@_);
    }
    return $id;

}

# $id = insertAmbiguity($dbh, $peptide, $proteinId);
# Insert ambiguity fragment
# Inputs: $dbh - database handle
#         $peptide - upper case peptide fragment
#         $proteinId - id to protein table
# Outputs: $id - ambiguity id, -1 if there was an error (checks for FK existence)
# Kristen Naegle 
# August 6, 2008
sub insertAmbiguity($$$){
    my ($dbh, $peptide, $proteinId) = @_;
    my $id = -1;
    my %protArgs;
    $protArgs{'id'} = $proteinId;
    my $protIdExist = checkForExistence($dbh, '*', 'protein', \%protArgs);
    if($protIdExist){
	my $sth=$dbh->prepare('INSERT into ambiguity (peptide, protein_id) VALUES (?,?)');
	$sth->execute($peptide, $proteinId);
	$id = getLastInsertId($dbh);
    }
    else{
	handleError('insertAmbiguity', 'FATAL: Foreign key existence problem, protein_id $proteinId does not exist', \@_);
    }
    return $id;
}

# $id = insertProteinExpression($dbh, $proteinId, $probesetId);
# Inserts into linker table protein_id and probeset_id
# Inputs: $dbh - database handle
#         $proteinId - FK id to protein table
#         $probeset_id - id to both expression and expression_ann
# Outputs: $id - -1 if failure to insert
# Kristen Naegle
# August 12, 2008
sub insertProteinExpression($$$){
    my ($dbh, $proteinId, $probeId) = @_;
    my $id = -1;
    my %protArgs; 
    $protArgs{'id'} = $proteinId;
    my $protIdExist = checkForExistence($dbh, '*', 'protein', \%protArgs);
    if($protIdExist){
	my $sth = $dbh->prepare('INSERT into protein_expression (protein_id, probeset_id) VALUES (?,?)');
	$sth->execute($proteinId, $probeId);
	$id = getLastInsertId($dbh);
    }
    else{
	handleError('insertProteinExpression', 'FATAL: Foreign key existence problem, protein_id $proteinId does not exist', \@_);
     
    }
    return $id;
}

# updateAccGeneFieldForProteinId($dbh, $proteinId, $acc_gene)
# Update the acc_gene field of protein table
# Inputs: $dbh - database handle
#         $proteinId - id to protein table
#         $acc_gene - acc_gene field of protein table to be set
# Kristen Naegle
# Feb 9, 2009
sub updateAccGeneFieldForProteinId($$$){
    my ($dbh, $proteinId, $acc_gene) = @_;
    
    my $sth = $dbh->prepare('UPDATE protein set acc_gene=? where id=?');
    $sth->execute($acc_gene, $proteinId);
    $sth->finish;
    
}

### -----------------Simple Selects--------------------------- ##

# $str = createSelectString($value, $table, \%args);
# Creates a select statement for using in a prepare statement. Creates select (value) from a single table ($table) where $key=$args{$key}
# Inputs: $value - value to select from table e.g. *
#         $table - the table to select from
#         \%args - hash ref, where keys define the column fields and the value of that $arg{$key} is the value it must be equal to.  More than one key means it will create a multiple and statement;
# Outptus: $str - SELECT $value FROM $table WHERE $key1=$args{$key1} AND $key2=$args{$key2} AND...
# Kristen Naegle
# March 11, 2008
sub createSelectString($$$){
    my ($value, $table, $args) = @_;

    my $str = "SELECT $value FROM $table WHERE ";
    my @keys = keys %$args;
    $str .= $keys[0]."=".$args->{$keys[0]};
    for (my $i=1; $i< scalar(@keys); $i++){
	$str .= " and ".$keys[$i]."=".$args->{$keys[$i]};
    }
    #print "Select Statement: \n $str\n";
    return $str;
}



# \@domainIds = returnDomainIdsForProteinId($dbh, $proteinId)
# Returns all the domains that exist for a particular protein
# Inputs: $dbh - database handler
#         $proteinId (FK in domain table) - id in protein table
# Outputs: \@domains - array of domains for that protein
# Kristen Naegle
# March 9, 2008
sub returnDomainIdsForProteinId($$){
    my ($dbh, $proteinId) = @_;
    my $sth = $dbh->prepare('SELECT * FROM domain where protein_id = ?');
    $sth->execute($proteinId) || die "Couldn't execute statement: ". $sth->errstr."\n";
    my $domains = returnArrayOfResultsOnCol($sth, 0);
    $sth->finish;
    return $domains;
    
}




#$pepId = returnPhosphopepIdOnPepAlignedAndProteinId($dbh, $aligned ,$proteinId)
# returns id to phosphopep table for something with the exact alignment from the same protein
# Inputs: $dbh - database handler
#         $aligned - aligned sequence
#         $proteinID - (FK) id of protein from where alignment was made
# Outputs: $id - id to phosphopep table
# Kristen Naegle
# March 10, 2008
sub returnPhosphopepIdOnPepAlignedAndProteinId($$$){
    my($dbh, $aligned, $proteinId) = @_;
    
    my $sth = $dbh->prepare('SELECT id FROM phosphopep where pep_aligned = ? and protein_id = ?');
    $sth->execute($aligned, $proteinId) or die "Couldn't execute statement: ".$sth->errstr;
    #expect a single value for this
    my($errorCode, $id) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleErrorCode('returnPhosphopepIdOnPepAlignedAndProteinId', 'More than one result for select statement', \@_);}
    $sth->finish;
    return ($errorCode, $id);

}

# $expId = returnExperimentId($dbh, $name);
# returns experiment id based on name
# Inptus: $dbh - database handler
#         $name - name of experiment (must be identical)
# Outputs: $expId - id to experiment table, -1 if it wasn't found, if there was more than one for the same name, returns the first
# Kristen Naegle
# March 2008
sub returnExperimentId($$){
    my ($dbh, $name) = @_;
    my $expID;
    my $sth = $dbh->prepare('SELECT id FROM experiment where name = ?');
    $sth->execute($name);
    my $errorCode;
    ($errorCode, $expID) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleErrorCode('returnExperimentId', 'More than one result for select statement', \@_);}
    $sth->finish;
    return $expID;
}

# $proteinId = returnProteinIdBySeqGene($dbh, $sequence, $gene, $species)
# Returns id to protein table based on the same sequence and geneName
# Inputs: $dbh - datbase handler
#         $sequence - the sequence (since if we don't find accession we need to look for the sequence). 
#         $gene - gene name in entrez
# Outputs: $protein_id - protein id of that accession number (-1 if does not exist in table) and if more than one exists..writes error, and returns first proteinId
# Kristen Naegle
# March 2, 2008
# Modified April 4th, require Species as well
sub returnProteinIdBySeqGene($$$$){
    my($dbh, $sequence, $gene, $species) = @_;
    my $proteinId;
    my $sth;
    if($gene eq 'NULL'){
	$gene =~ s/\'//g;
	$sth = $dbh->prepare('SELECT * FROM protein where sequence = ? and species = ? and acc_gene is NULL');
	$sth->execute($sequence, $species);
    }
    else{
	$sth = $dbh->prepare('SELECT * FROM protein where sequence = ? and species=? and acc_gene= ?');
	$sth->execute($sequence, $species, $gene);
    }

   
    my $errorCode;
    ($errorCode, $proteinId) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleErrorCode('returnProteinIdBySeqGene', 'More than one result for select statement', \@_);}
    $sth->finish;
    return $proteinId;
}

# $proteinId = returnProteinIdBySeq($dbh, $sequence, $gene, $species)
# Returns id to protein table based on the same sequence and species
# Inputs: $dbh - datbase handler
#         $sequence - the sequence (since if we don't find accession we need to look for the sequence). 
#         $species - two word formal name of species
# Outputs: $protein_id - protein id of that accession number (-1 if does not exist in table) and if more than one exists..writes error, and returns first proteinId
# Kristen Naegle
# August 6, 2008
sub returnProteinIdBySeq($$$){
    my($dbh, $sequence, $species) = @_;
    my $proteinId;
    my $sth;
    $sth = $dbh->prepare('SELECT * FROM protein where sequence = ? and species = ?');
    $sth->execute($sequence, $species);
    
    my $errorCode;
    ($errorCode, $proteinId) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleErrorCode('returnProteinIdBySeqGene', 'More than one result for select statement', \@_);}
    $sth->finish;
    return $proteinId;
}

# $accId = returnAccIdByValue($dbh, $accValue);
# return accession Id for an accession value. Returns -1 if it doesn't exist
# Inputs: $dbh - database handle
#         $acc - accession value
# Outputs: $accId - id to that accession, -1 if it doesn't exist in table
# Kristen Naegle
# March 2, 2008
sub returnAccIdByValue($$){
    my($dbh, $acc) = @_;
    my $accId;
    my $sth = $dbh->prepare('SELECT * FROM acc where value = ?');
    $sth->execute($acc);
    my $errorCode;
    ($errorCode, $accId) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleErrorCode('returnAccIdByValue', 'More than one result for select statement', \@_);}
    $sth->finish;
    return $accId;
}

# $accId = returnAccIdByValueProteinId($dbh, $acc, $proteinId);
# return accession Id for an accession value and proteinid. Returns -1 if it doesn't exist
# Inputs: $dbh - database handle
#         $acc - accession value
#         $proteinId - (FK) protein id
# Outputs: $accId - id to that accession, -1 if it doesn't exist in table
# Kristen Naegle
# March 2, 2008
sub returnAccIdByValueProteinId($$$){
    my($dbh, $acc, $proteinId) = @_;
    my $accId;
    my $sth = $dbh->prepare('SELECT * FROM acc where value = ? and protein_id= ?');
    $sth->execute($acc, $proteinId);
    my $errorCode;
    ($errorCode, $accId) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleErrorCode('returnAccIdByValue', 'More than one result for select statement', \@_);}
    $sth->finish;
    return $accId;
}

# $accArr = returnAccArrOfTypeForProtein($dbh, $proteinId, $type);
# Returns an array of accession of the type indicated ('swissprot', 'gene_synonym', etc)
# Inputs: $dbh - database handle
#         $proteinId - (FK) protein id
#         $type - type value in acc. table
# Outputs: $accArr - ref to array of accessions (the values), if no accessions of that type exist for a protein, then first entry is a -1
# Kristen Naegle
# June 18, 2009
sub returnAccArrOfTypeForProtein($$$){
    my($dbh, $proteinId, $type) = @_;
    my $sth = $dbh->prepare("SELECT value FROM acc where protein_id=? and type=?");
    $sth->execute($proteinId, $type);
    my $errorCode;
    my $accArr = returnArrayOfResultsOnCol($sth, 0);
    $sth->finish;
    return $accArr;
}

# $proteinId = returnProteinIdByAcc($dbh, $accId);
# return (FK) protein ID based on accession Id 
# Inputs: $dbh - database handle
#         $acc - value in accession table
# Outputs: $proteinId (FK) - protein_id id to protein table
# Kristen Naegle
# March 2, 2008
sub returnProteinIdByAcc($$){
    my($dbh, $accId) = @_;
    my $proteinId;
    my $sth = $dbh->prepare('SELECT protein_id FROM acc where value = ?');
    $sth->execute($accId);
    my $errorCode;
    ($errorCode, $proteinId) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleError('returnProteinIdByAcc', 'More than one result for select statement', \@_);}
    $sth->finish;
    return $proteinId;
}

# $accArr = returnAccValuesByProteinId($dbh, $proteinId);
# returns an array of accession values (or all accession values) for a protein (based onid to protein table)
# Inputs: $dbh - database handle
#         $proteinId - (FK) id to protein table
# Outputs: $accArr - array of accession for that protein (that currently exist)
# Kristen Naegle
# March 3, 2008
sub returnAccValuesByProteinId($$){
    my($dbh, $proteinId) = @_;
    my $sth = $dbh->prepare('SELECT value FROM acc where protein_id = ?');
    $sth->execute($proteinId) || die  "Couldn't execute statement: " . $sth->errstr;
    my $accArr = returnArrayOfResultsOnCol($sth, 0);
    $sth->finish;
    return $accArr;
}

# $sequence = returnSeqByProteinId($dbh, $proteinId)
# Returns the protein sequence for a protein id
# Inputs: $dbh - database handle
#         $proteinId - (FK) id to protein table
# Outputs: $sequence - protein sequence
# Kristen Naegle
# March 10, 2008
sub returnSeqByProteinId($$){
    my($dbh, $proteinId) = @_;
    my ($errorCode, $sequence);
    my $sth = $dbh->prepare('SELECT sequence FROM protein where id = ?');
    $sth->execute($proteinId) || die "Couldn't execute statement: " . $sth->errstr;
    ($errorCode, $sequence) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleErrorCode('returnSeqByProteinId', 'More than one result for select statement', \@_);}
    $sth->finish;
    return $sequence;

}

# $name = returnNameByProteinId($dbh, $proteinId)
# Returns the protein sequence for a protein id
# Inputs: $dbh - database handle
#         $proteinId - (FK) id to protein table
# Outputs: $name - long protein name
# Kristen Naegle
# March 10, 2008
sub returnNameByProteinId($$){
    my($dbh, $proteinId) = @_;
    my ($errorCode, $name);
    my $sth = $dbh->prepare('SELECT name FROM protein where id = ?');
    $sth->execute($proteinId) || die "Couldn't execute statement: " . $sth->errstr;
    ($errorCode, $name) = returnSingleResultOnCol($sth, 0);
    if($errorCode){handleErrorCode('returnSeqByProteinId', 'More than one result for select statement', \@_);}
    $sth->finish;
    return $name;

}


# $MS_id = returnMSIdByPeptideExperiment($dbh, $phosphopeptide, $$proteinId, expId)
# Returns the MS_id for an MS table with a phosphopeptide and experiment id (-1 otherwise)
# Inputs: $dbh - database handle
#         $phosphopeptide - phosphopeptide 
#         $proteinId - protein id in MS table
#         $expId - (FK) experiment id
# Outputs: $MS_id - id of a MS table entry
# Kristen Naegle
# March 17, 2008
sub returnMSIdByPeptideExperiment($$$$){
    my ($dbh, $pep, $proteinId, $expId) = @_;

    my $MS_id;

    my $sth=$dbh->prepare('SELECT id FROM MS WHERE phosphopep=? AND experiment_id=? and protein_id=?');
    $sth->execute($pep, $expId, $proteinId);
    my($errorCode, $id) = returnSingleResultOnCol($sth, 0);
    if($errorCode){
	handleError('returnMSIdByPeptideExperiment', 'More than one MS_id for peptide and experiment', \@_);
    }
    if(!$errorCode){
	$MS_id = $id;
    }
    else{
	$MS_id = -1;
    }
    return $MS_id;

}


# $MS_id = returnMSIdByPhosphopepExperimentId($dbh, $phosphopeptide,$expId);
# Returns the MS_id for an MS table with a phosphopeptide and experiment id (-1 otherwise)
# Inputs: $dbh - database handle
#         $phosphopeptide - phosphopeptide 
#         $expId - (FK) experiment id
# Outputs: $MS_id - id of a MS table entry
# Kristen Naegle
# July 23, 3009
sub returnMSIdByPhosphopepExperimentId($$$){
    my ($dbh, $pep, $expId) = @_;

    my $MS_id;

    my $sth=$dbh->prepare('SELECT id FROM MS WHERE phosphopep=? AND experiment_id=?');
    $sth->execute($pep, $expId);
    my($errorCode, $id) = returnSingleResultOnCol($sth, 0);
    if($errorCode){
	handleError('returnMSIdByPeptideExperiment', 'More than one MS_id for peptide and experiment', \@_);
    }
    if(!$errorCode){
	$MS_id = $id;
    }
    else{
	$MS_id = -1;
    }
    return $MS_id;

}
# (\@MS_phosphopepIds, \@phosphopepIds) = returnMSPhosphopepTableForMSId($dbh, $MS_id);
# Returns the array of MS_phosphopep ids and phosphopep Ids in MSPhosphopep with MS_id
# Inputs: $dbh - database handle
#         $MS_id - (FK) id to MS table
# Outputs: $MS_phosphopepIds - (FK) ref to array of MS_phosphopepIds
#          $phosphopepIds - (FK) ref to array of phosphopep Ids
# Kristen Naegle
# March 17, 2008
sub returnMSPhosphopepTableForMSId($$){
    my($dbh, $MS_id) = @_;
    
    my $sth = $dbh->prepare('SELECT * FROM MS_phosphopep WHERE MS_id = ?');
    
    $sth->execute($MS_id);
    
    my $MS_phosphopepIds = returnArrayOfResultsOnCol($sth, 0);
    $sth->execute($MS_id);
    my $phosphopepIds = returnArrayOfResultsOnCol($sth, 2);
    return($MS_phosphopepIds, $phosphopepIds);
}
# (\@MS_phosphopepIds, \@phosphopepIds) = returnPhosphopepTableForMSId($dbh, $MS_id)
# Returns the array of MS_phosphopep ids and phosphopep Ids in MSPhosphopep with MS_id
# Inputs: $dbh - database handle
#         $MS_id - (FK) id to MS table
# Outputs: $MS_phosphopepIds - (FK) ref to array of MS_phosphopepIds
#          $phosphopepIds - (FK) ref to array of phosphopep Ids
# Kristen Naegle
# March 17, 2008
sub returnPhosphopepTableForMSId($$){
    my($dbh, $MS_id) = @_;
    
    my $sth = $dbh->prepare('SELECT phosphopep.id, phosphopep.pep_aligned FROM MS_phosphopep join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id WHERE MS_id = ?');
    
    $sth->execute($MS_id);
    
    my $phosphopepIds = returnArrayOfResultsOnCol($sth, 0);
    $sth->execute($MS_id);
    my $phosphopeps = returnArrayOfResultsOnCol($sth, 1);
    return($phosphopepIds, $phosphopeps);
}

# (\@MS_phosphopepIds, \@phosphopepIds) = returnPhosphopepTableSiteSpeciesForMSId($dbh, $MS_id, $siteType, $species)
# Returns the array of MS_phosphopep ids and phosphopep Ids in MSPhosphopep with MS_id
# Inputs: $dbh - database handle
#         $MS_id - (FK) id to MS table
#         $siteType - type of site 'Y' 'S' or 'T'
#         $species - species wanted, for mixed datasets like pELM and phosphosite
# Outputs: $MS_phosphopepIds - (FK) ref to array of MS_phosphopepIds
#          $phosphopepIds - (FK) ref to array of phosphopep Ids
# Kristen Naegle
# March 17, 2008
sub returnPhosphopepTableSiteSpeciesForMSId($$$$){
    my($dbh, $MS_id, $siteType, $species) = @_;
    
    my $sth = $dbh->prepare('SELECT phosphopep.id, phosphopep.pep_aligned FROM MS_phosphopep join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id join protein on protein.id=phosphopep.protein_id WHERE MS_id = ? and site_type=? and species=?');
    
    $sth->execute($MS_id, $siteType, $species);
    
    my $phosphopepIds = returnArrayOfResultsOnCol($sth, 0);
    $sth->execute($MS_id, $siteType, $species);
    my $phosphopeps = returnArrayOfResultsOnCol($sth, 1);
    return($phosphopepIds, $phosphopeps);
}

# $dataId = returnDataIdByMost($dbh, $type, $run, $label, $MS_id);
# Returns an existing id for a data entry, or a -1 if it doesn't match
# Inputs: $dbh - database handle
#         $type - type 
#         $run - run number/description
#         $label - label for data
#         $MS_id - (FK) MS table id
# Outputs: $dataId - id entry to data table that matches input criteria
# Kristen Naegle
# March 18, 2008 
sub returnDataIdByMost($$$$$){
    my ($dbh, $type, $run, $label, $MS_id) = @_;

    my $sth = $dbh->prepare('SELECT id from data where type= ? and run= ? and label =? and MS_id = ?');
    $sth->execute($type, $run, $label, $MS_id);
    my $id = returnSingleResultOnCol($sth, 0);
    return $id;

}

# \@MS_ids = returnMSIdsForExpId($dbh, $expId);
# Returns ref to array of MS ids for a particular experiment id
# Inputs: $dbh - database handle
#         $expId - (FK) experiment id
# Outputs: $MS_ids - ref to array of MS ids, first value is -1 if there were no results
# Kristen Naegle
# March 21, 2008
sub returnMSIdsForExpId($$){
    my ($dbh, $expId) = @_;
    my $sth=$dbh->prepare('SELECT id from MS where experiment_id=?');
    $sth->execute($expId);
    my $MS_ids = returnArrayOfResultsOnCol($sth, 0);
    return $MS_ids;
}

# \@proteinIds = returnProteinIdsForExpId($dbh, $expId)
# Returns unique array of protein ids in an experiment
# Inputs: $dbh - database handle
#         $expId - (FK) experiment Id
# Outputs: \@proteinIds - array of protein ids that belong to that experiment
# Kristen Naegle
# April 15, 2008
sub returnProteinIdsForExpId($$){
    my($dbh, $expId) = @_;
    my $sth = $dbh->prepare('SELECT protein_id from MS where experiment_id = ? group by protein_id');
    $sth->execute($expId);
    my $proteinIds = returnArrayOfResultsOnCol($sth,0);
    return $proteinIds;

}

# $GO_id = returnGOIdForGOTerm($dbh, $GO);
# Returns the GO ID based on the GO term value and aspect (F, C, or P)
# Inputs: $dbh - database handle
#         $GO - GO value
#         $aspect - aspect of GO term to look at 
# Outputs: $GO_id - key to GO table (-1 if doesn't exist)
# Kristen Naegle
# April 15, 2008
sub returnGOIdForGOTerm($$$){
    my ($dbh, $GO, $aspect) = @_;

    my $sth = $dbh->prepare('SELECT id from GO where GO = ? and aspect=?');
    $sth->execute($GO, $aspect);
    my $GO_id = returnSingleResultOnCol($sth, 0);
    return $GO_id;
}

# $protein_GO_id = returnProteinGOId($dbh, $protein_id, $GO_id);
# Returns linkage table id, protein_GO, for protein_id and a GO_id
# Inputs: $dbh -database handle
#         $protein_id - (FK) protein_id or protein.id
#         $GO_id (FK) - GO.id 
# Outputs: $protein_GO_id - id to protein_GO table
# Kristen Naegle
# April 15, 2008
sub returnProteinGOId($$$){
    my ($dbh, $protein_id, $GO_id) = @_;
    
    my $sth = $dbh->prepare('SELECT id from protein_GO where protein_id=? and GO_id=?');
    $sth->execute($protein_id, $GO_id);
    my $protein_GO_id = returnSingleResultOnCol($sth, 0);
    return $protein_GO_id;

}

# $phosphopep_id = returnPhosphopepIdOnAccPepTryps($dbh, $acc, $pepTryps);
# Return the phosphopep_id given an accession value and a tryps peptide
# Inputs: $dbh -database handle
#         $acc - value of accession
#         $pepTryps - trypsinized peptide
# Outputs: $phosphopep_id - id to corresponding phosphopep (-1 if dne)
# Kristen Naegle
# October 13, 2008
sub returnPhosphopepIdOnAccPepTryps($$$){
    my($dbh, $acc, $pepTryps) = @_;

    # find protein id based on sp
    my $proteinId = returnProteinIdByAcc($dbh, $acc);
    if($proteinId == -1){
	handleError('returnPhosphopepIdOnAccPepTryps', "Could not find protein id with acc $acc", \@_);
	return -1;
    }
    
    # find phosphopep_id on protein_id and pep_tryps
    my $sth = $dbh->prepare('select id from phosphopep where protein_id=? and pep_tryps =?');
    $sth->execute($proteinId, $pepTryps);
    my $phosphopep_id = returnSingleResultOnCol($sth, 0);
    return $phosphopep_id;
    

}
# $phosphopep_id = returnPhosphopepIdOnAccPepAligned($dbh, $acc, $pep);
# Return the phosphopep_id given an accession value and a tryps peptide
# Inputs: $dbh -database handle
#         $acc - value of accession
#         $pepTryps - trypsinized peptide
# Outputs: $phosphopep_id - id to corresponding phosphopep (-1 if dne)
# Kristen Naegle
# Nov 5, 2008 - needed to switch from trypsinized peptide to aligned for pelm_kinase annotations
sub returnPhosphopepIdOnAccPepAligned($$$){
    my($dbh, $acc, $pep) = @_;

    # find protein id based on sp
    my $proteinId = returnProteinIdByAcc($dbh, $acc);
    if($proteinId == -1){
	handleError('returnPhosphopepIdOnAccPepAligned', "Could not find protein id with acc $acc", \@_);
	return -1;
    }
    
    # find phosphopep_id on protein_id and pep_tryps
    my $sth = $dbh->prepare('select id from phosphopep where protein_id=? and pep_aligned =?');
    $sth->execute($proteinId, $pep);
    my $phosphopep_id = returnSingleResultOnCol($sth, 0);
    return $phosphopep_id;
    

}

# $phosphopepPrediction_id = returnPhosphoPepPredId($dbh, $source, $value, $phosphopep_id);
# Returns the id to the record for a prediction in phosphopep_prediction table. -1 if record DNE
# Inputs: $dbh - database handle
#         $source - source entry in table
#         $value - value entry in table
#         $phosphopep_id - FK to phosphopep
# Outputs: $id - -1 if record DNE
# Kristen Naegle
# October 14, 2008
sub returnPhosphopPepPredId($$$$){
    my ($dbh, $source, $value, $phosphopep_id) = @_;
    my $sth = $dbh->prepare('SELECT id from phosphopep_prediction where source=? and value=? and phosphopep_id=?');
    $sth->execute($source, $value, $phosphopep_id);
    my $id = returnSingleResultOnCol($sth, 0);
    return $id;
}

# $id = returnAmbIdForFragProtein($dbh, $fragment, $proteinId)
# return id of ambiguity table for a specific fragment and protein
# Inputs: $dbh - database handle
#         $fragment - peptide fragment, capitalized before checked
#         $proteinId - protein id
# Outputs: $id - id of ambiguity table, -1 if DNE
# Kristen Naegle
# August 6, 2008
sub returnAmbIdForFragProtein($$$){
    my ($dbh, $fragment, $proteinId) = @_;
    $fragment = uc($fragment);
    my $sth = $dbh->prepare('SELECT id from ambiguity where peptide=? and protein_id=?');
    $sth->execute($fragment, $proteinId);
    my ($errorCode, $result) = returnSingleResultOnCol($sth, 0);
    return $result;

}

# $hashReslt = returnPhosphopepRecordOnSiteProteinId($dbh, $site_pos, $proteinId);
# Returns the record, in hash form, for phosphopeps that have site and protein as specified. Built intentionally for finding peps with specific phosphorylation sties from alignment work
# Inputs: $dbh - database handle
#         $site_pos - position of site phosphorylation
#         $proteinId - id of protein
# Outputs: $hashResult - ref. to hash with keys 'id', 'pep_aligned', 'site_type', 'pfam_site' .. complains if more than one result, if no results all are -1
# Kristen Naegle
# Nov 5th, 2008
sub returnPhosphopepRecordOnSiteProteinId($$$){
    my ($dbh, $site_pos, $proteinId) = @_;
    
    my $sth = $dbh->prepare('SELECT id, pep_aligned, site_type, pfam_site from phosphopep where site_pos=? and protein_id=?');
    $sth->execute($site_pos, $proteinId);
    #my $hashResult = 
    my @row = $sth->fetchrow_array;
    my %hash;
    if(not defined $row[0]){ 
#	print "setting row\n";
	@row = (-1, -1, -1, -1); 
    }
    $hash{'id'} = $row[0];
    $hash{'pep_aligned'} = $row[1];
    $hash{'site_type'} = $row[2];
    $hash{'pfam_site'} = $row[3];
    if($row[0] != -1){
	@row = $sth->fetchrow_array; 
	if(@row){ 
	    print "ERROR: Found too many phosphopep records on site_pos $site_pos and protein_id $proteinId\n"; 
	}
    }
    return \%hash;
    
}
# ---------------------REFSEQ Datbase -------------------- #

# $sth = returnInsertStatementRefSeq($dbh)
# Prepare a statement handle for an insert gi, xp, name, species, sequence and date are the order of statement
# Inputs: $dbh - refseq database handle
# Outputs: $sth - statement handler
# Kristen Naegle
# April 30, 2008
sub returnInsertStatementRefSeq($){
    my ($dbh) = @_;
    my $sth=$dbh->prepare('INSERT INTO refseq (gi, xp, name, species, sequence, date) VALUES (?, ?, ?, ?, ?, ?)');
    return $sth;
}

# insertRefSeqEntry($sth, $gi, $xp, $name, $species, $sequence, $date)
# Executes statement handle insert (see returnInsertStatementRefSeq)
# Inputs: $sth - statement handle
#         $gi - gi accession
#         $xp - ref seq accession
#         $name - protein name
#         $species - two word name of protein
#         $sequence - protein sequence
#         $date - Month Year 
# Kristen Naegle
# April 30, 2008
sub insertRefSeqEntry($$$$$$$){
    my ($sth, $gi, $xp, $name, $species, $sequence, $date) = @_;
    $sth->execute($gi, $xp, $name, $species, $sequence, $date);

}

# \%hash = returnGIOnFragmentMatch($dbh, $fragment, $species)
# Given a species and trypsinized fragment (or any peptide) return a hash of the values that match in refseq database
# Inputs: $dbh - database handle
#         $fragment - peptide to match to subset of sequence
#         $species - two word name of species
# Outputs: $hash - reference to hash with key being gi and value is pointer to hash with keys: 'sequence', 'name', 'species'
# Kristen Naegle
# April 30, 2008
sub returnGIOnFragmentMatch($$$){
    my ($dbh, $fragment, $species) = @_;
    my %hash;
    my $sth = $dbh->prepare('SELECT gi, sequence, name, species from refseq where sequence REGEXP ? and species REGEXP ?');
    $sth->execute($fragment, $species);
    while (my @row = $sth->fetchrow_array){
	my $gi = $row[0];
	my %newHash;
	$newHash{'sequence'} = $row[1];
	$newHash{'name'} = $row[2];
	$newHash{'species'} = $row[3];
	$hash{$gi} = \%newHash;
    }
    return \%hash;

}

# $alignedHash = returnAlignedPeptidesForMSId($$)
# Returns a hash of hashes of aligned peptides for an MS. Top key is a number and second hash has keys 'phosphopep.aligned', 'phosphopep.type', 'phosphopep.site'
# 
# Kristen Naegle
# June 18, 2008
sub returnAlignedPeptidesForMSId($$){
    my ($dbh, $MSId) = @_;
    my %alignedHash;
    my $sth = $dbh->prepare('Select phosphopep.pep_aligned, phosphopep.site_pos, phosphopep.site_type from MS join MS_phosphopep on MS_phosphopep.MS_id=MS.id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS.id=?');
    $sth->execute($MSId);
    my $count = 0;
    while(my @row=$sth->fetchrow_array){
	$count += 1;
	my %hash;
	($hash{'phosphopep.aligned'}, $hash{'phosphopep.site'}, $hash{'phosphopep.type'}) = @row;
	 $alignedHash{$count} = \%hash;
    }
     
    return \%alignedHash;
}

# \@proteinIds = returnProteinIdsByGeneSpecies($dbh, $acc_gene, $species)
# Return all the protein ids that match a particular gene in a species
# Inputs: $dbh - database handle
#         $acc_gene - protein table entry for gene accession
#         $species - species
# Outputs: $proteinIds - reference to array of protein ids with that gene name and species
# Kristen Naegle
# May 12th, 2009
sub returnProteinIdsByGeneSpecies($$$){
    my ($dbh, $acc_gene, $species) = @_;
    my @proteinIds;

    my $sth = $dbh->prepare('SELECT id from protein where acc_gene = ? and species=?');
    $sth->execute($acc_gene, $species);
    while(my @row = $sth->fetchrow_array){
	push @proteinIds, $row[0];
    }
    return \@proteinIds;
}

# \@proteinIds = returnProteinIdsByGeneSynSpecies($dbh, $geneSynonyms, $species)
# Return all the protein ids that match a particular gene (by synonym) in a species
# Inputs: $dbh - database handle
#         $gene_synonyms - ref to an array of gene synonyms
#         $species - species
# Outputs: $proteinIds - reference to array of protein ids with that gene name and species
# Kristen Naegle
# May 12th, 2009
sub returnProteinIdsByGeneSynSpecies($$$){
    my ($dbh, $geneSynonyms, $species) = @_;
    my @proteinIds;

    my $geneSynStr = makeSelectINString($geneSynonyms, 1); #quoted string
    print "DEBUG: Gene Syns: $geneSynStr and species: $species\n";
    my $sth = $dbh->prepare("SELECT protein.id from acc join protein on acc.protein_id=protein.id where acc.value IN $geneSynStr and species=?");
    $sth->execute($species);
    
    while(my @row = $sth->fetchrow_array){
	push @proteinIds, $row[0];
    }
    return \@proteinIds;
}

# $proteinId = returnProteinIdBySeqSpecies($dbh, $sequence, $species);
# Return a protein id that matches based on sequence and species
# Inputs: $dbh - database handle
#         $sequence - amino acid sequence 
#         $species - species
# Outputs: $proteinId - protein id of protein table with matching species and sequence
# Kristen Naegle
# May 12th, 2009
sub returnProteinIdBySeqSpecies($$$){
    my ($dbh, $sequence, $species) = @_;
    my $proteinId;

    
    my $sth = $dbh->prepare('SELECT id from protein where sequence=? and species=?');
    $sth->execute($sequence, $species);
    
    $proteinId = returnSingleResultOnCol($sth, 0);

    return $proteinId;
}

# ($phosphopep, $proteinId) = returnMSPhosphopepForMSId($dbh, $MSId)
# Return MS.phosphopep based on MS.id 
# Written for translation between databases
# Inputs: $dbh - database handle
#         $MSId - id to MS table
# Outputs: $phosphopep - phosphopep in MS table
#          $proteinId - protein_id in MS table
# Kristen Naegle
# July 23, 3009
sub returnMSPhosphopepProteinIdForMSId($$){
    my ($dbh, $MSId) = @_;

    my $sth = $dbh->prepare('SELECT * from MS where id=?');
    $sth->execute($MSId);
#    my $phosphopep = returnSingleResultOnCol($sth, 0);
    my $results = returnMultipleResultsForExSTH($sth);
    my $hash = $results->[0];
    my $phosphopep = $hash->{'phosphopep'};
    my $proteinId = $hash->{'protein_id'};
    return ($phosphopep, $proteinId) ;

}

1;
