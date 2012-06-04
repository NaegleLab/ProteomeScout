use strict;
use warnings;
use DBI;
use DBTools::resultParsing;
use DBTools::insertFunctions;
use DBTools::enrichment;
use entrezTools;
use errorHandling;


# ($errorCode, $pfamSite) = returnPfamSite($dbh, $proteinId, $pos);
# Returns the pfam site in which a site lays (~~~ , if it doesn't reside in known pfam site)
# Inputs: $dbh - database handler
#         $proteinId - id in protein table
#         $pos - position of site
# Outputs: $errorCode - if more than one domain exists for a site, or no domains exist for that protein
#          $pfamSite - domain where site falls in, should always have value unless error
# Kristen Naegle
# March 8, 2008
sub returnPfamSite($$$){
    my ($dbh, $proteinId, $pos) = @_;
    my ($errorCode, $pfamSite);
    $errorCode = 0;
    my $sth = $dbh->prepare('SELECT label FROM domain where ? >= start and ? <= stop and protein_id = ? and source REGEXP \'PFAM\'');
    $sth->execute($pos, $pos, $proteinId) || die "Couldn't execute statement: " . $sth->errstr;
    my ($errCodeFetch, $label) = returnSingleResultOnCol($sth, 0);

    if($label eq '-1'){
	my $sth_pid = $dbh->prepare('SELECT label FROM domain where protein_id = ? and source REGEXP \'PFAM\'');
	$sth_pid->execute($proteinId) || die "Couldn't execute statement: ".$sth->errstr;
	my @row = $sth_pid->fetchrow_array;
	if(@row){
	    $pfamSite = "~~~";
	}
	else{
	    $errorCode = 1; #Empty..no pfam predictions
	    handleError('returnPfamSite', "No pfam predictions exist for protein Id $proteinId", \@_);
	}
    }
    elsif($errCodeFetch){
	$errorCode = 2; #more than one site..this shouldn't have happened!
	handleError('returnPfamSite', "BAD! More than one pfam prediction for the same region exist for protein Id $proteinId", \@_);
    }
    else{
	$pfamSite = $label;
    }
    return($errorCode, $pfamSite);
}

# $exist = checkForExistence($dbh, $value, $table, $args);
# Checks for the existence of a value in a table according to args
# Inputs: $dbh - database handler
#         $value - value to select from table e.g. *
#         $table - the table to select from
#         \%args - hash ref, where keys define the column fields and the value of that $arg{$key} is the value it must be equal to.  More than one key means it will create a multiple and statement;
# Outputs: $exist - 1 if it existed in table 0 else
# Kristen Naegle
# March 11, 2008
sub checkForExistence($$$$){
    my ($dbh, $value, $table, $args) = @_;
    my $flag;
    my $statement = createSelectString($value, $table, $args);
    my $sth = $dbh->prepare($statement);
    $sth->execute;
    my $array = returnArrayOfResultsOnCol($sth, 0);
    if($array->[0] != -1){
	$flag = 1;
    }
    else{
	$flag = 0;
    }
    return $flag;
}

# \@missingAcc = returnMissingProteins($dbh, \@acc)
# Return the accessions in input acc array that are not in the database
# Inputs: $dbh - database handle
#         $accArr - reference to array of accessions to check
# Outputs: $missingAcc - ref to array of accessions not found in database.
#          \%accHash - ref to hash has keys equal to accessions in database and value equal to protein_id for accession
# Kristen Naegle
# March 24, 2008
sub returnMissingProteins($$){
    my ($dbh, $accArr) = @_;
    my @missAcc;
    my %accHash;
    my $sth= $dbh->prepare('SELECT protein_id FROM acc WHERE value=?');
    foreach my $acc (@$accArr){
	$sth->execute($acc);
	my ($errorCode, $val) = returnSingleResultOnCol($sth, 0);
	if($val == -1){
	    push @missAcc, $acc;
	}
	else{
	    $accHash{$acc} = $val;
	}

    }
    return (\@missAcc, \%accHash);

}

# $proteinIdArr = returnProteinIdsWithNoDomains($dbh);
# Returns array of protein ids that do not have any domain predictions yet
# Inputs: $dbh - database handle
# Outputs: $proteinIdArr - ref. to array of protein ids 
# Kristen Naegle
# March 29, 2008
sub returnProteinIdsWithNoDomains($){
    my ($dbh) = @_;
    my $sth= $dbh->prepare('SELECT protein.id FROM protein LEFT OUTER JOIN domain on domain.protein_id=protein.id where (domain.protein_id is NULL)');
    $sth->execute();
    my $proteinIdArr = returnArrayOfResultsOnCol($sth, 0);
    return $proteinIdArr;
}

# \%hash = returnPhosphopepIdHashWithNoPredictions($dbh, $source)
# Return those phosphopeps that do not have predictions with a particular source
# Inputs: $dbh - database Handle
#         $source - source (e.g. scansite, or pelm)
# Outputs: \%hash - hash with key phosphopep.id and value phosphopep.pep_aligned
# Kristen Naegle
# April 5, 2008 
# Fixed Oct. 29, 2008 - had introduced error in SQL select due to appearance of pelm_kinases in prediction table
sub returnPhosphopepIdHashWithNoPredictions($$){
    my ($dbh, $source) = @_;
    my %hash;
   # my $sth=$dbh->prepare('SELECT phosphopep.id, phosphopep.pep_aligned from phosphopep LEFT OUTER JOIN phosphopep_prediction on phosphopep_prediction.phosphopep_id=phosphopep.id where (phosphopep_prediction.phosphopep_id is NULL or phosphopep_prediction.source NOT REGEXP ?)');
    my $sth = $dbh->prepare('SELECT phosphopep.id, phosphopep.pep_aligned from phosphopep where phosphopep.id NOT IN (select phosphopep_id from phosphopep_prediction where source REGEXP ?)');
    $sth->execute($source);
    my @row;
    while(@row = $sth->fetchrow_array){
	$hash{$row[0]} = $row[1];

    }
    return \%hash;
    
}

# $sth = returnSingleTableFieldSearchSTH($dbh, $table, $field, $expression);
# Returns statement handle for dbh based on SELECT * from $TABLE where $FIELD $XPRESSION ?"
# Inputs: $dbh - database handle
#         $table - the table you are choosing from
#         $field - table column name you are comparing
#         $expression - expression value (e.g. IS IN, REGEXP, REGEXP BINARY, =, etc.)
# Outputs: $sth - statement handle
# Kristen Naegle
# April 2, 2008
sub returnSingleTableFieldSearchSTH($$$$){
    my ($dbh, $table, $field, $expression) = @_;    
    
    my $select = "SELECT * from $table where $field $expression ?";
    print "SELECT Statement: $select\n";
    my $sth = $dbh->prepare($select);
    return $sth;

} 

# (\%hash, $fieldHashRef) = returnResultsSingleTableSearch($dbh, $table, $field, $expression, $valueArr);
# Conglomerate a general table query into a hash result
# Returns statement handle for dbh based on SELECT * from $TABLE where $FIELD $XPRESSION ?"
# Inputs: $dbh - database handle
#         $table - the table you are choosing from
#         $field - table column name you are comparing
#         $expression - expression value (e.g. IS IN, REGEXP, REGEXP BINARY, =, etc.)
#         $valueArr - ref. to an array of values for which the expression is equal too.
# Outputs: \%hash - hash with keys equal to the values of input and value an array of results. $fieldHashRef gives translation from field to array col. ONLY HANDLES ONE RETURN.
# Kristen Naegle
# April/November 2008
sub returnResultsSingleTableSearch($$$$$){
    my ($dbh, $table, $field, $expression, $valueArr) = @_;
    my $sth = returnSingleTableFieldSearchSTH($dbh, $table, $field, $expression);
    my $count = 0;
    my %hash;
#    my %fieldHash;
    my $fieldHashRef;
    if(!(uc($expression) eq 'IN')){
	foreach my $val (@$valueArr){
	    $sth->execute($val);
	    #print "executing with $val\n";
	    my @row = $sth->fetchrow_array;
	    #print "Result @row\n";
	    if(@row){
		push @{$hash{$val}}, @row;
		#if(!$count){ #first time through, set up hash
		if(not defined $fieldHashRef){
		    $fieldHashRef = returnFieldHashFromSTH($sth);
		}
	
	    }
	    else{
		push @{$hash{$val}}, -1;
	    }
	    $count +=1;
	}
	 
	
    } # end if not an IS IN 

    else{
	my $str = '('.$valueArr->[0];
	for (my $i=1; $i < scalar(@$valueArr); $i++){
	    $str .= ','.$valueArr->[$i];
	}
	$str .= ')';
	print "String: $str\n";
	$sth->execute($str);
	$fieldHashRef = returnFieldHashFromSTH($sth);
	my @row = $sth->fetchrow_array;
	push @{$hash{1}}, @row;
	
    }

    return (\%hash, $fieldHashRef);

}

# (\%hash) = returnMultipleResultsSingleTableSearch($dbh, $table, $field, $expression, $valueArr);
# Conglomerate a general table query into a hash result
# Returns statement handle for dbh based on SELECT * from $TABLE where $FIELD $XPRESSION ?"
# Inputs: $dbh - database handle
#         $table - the table you are choosing from
#         $field - table column name you are comparing
#         $expression - expression value (e.g. IS IN, REGEXP, REGEXP BINARY, =, etc.)
#         $valueArr - ref. to an array of values for which the expression is equal too.
# Outputs: \%hash - hash with keys equal to the values of input and value an array of hashes for multiple row effects. The second level of hashes have fields equal to field of interest and value equal to the value of that table column.  
# Kristen Naegle
# April/November 2008 - Does not currently handle IN statements
sub returnMultipleResultsSingleTableSearch($$$$$){
    my ($dbh, $table, $field, $expression, $valueArr) = @_;
    my $sth = returnSingleTableFieldSearchSTH($dbh, $table, $field, $expression);
    my $count = 0;
    my %hash;
    my $fieldHashRef;
    foreach my $val (@$valueArr){
	$sth->execute($val);
	if(not defined $fieldHashRef){
	    $fieldHashRef = returnFieldHashFromSTH($sth);
	}
	#print "executing with $val\n";
	while (my @row = $sth->fetchrow_array){
	    #print "Result: @row\n";
	    my %hR;
	    foreach my $field (keys %$fieldHashRef){
		my $col = $fieldHashRef->{$field};
		#print "$field maps to $col\n";
		$hR{$field} = $row[$col];
	    }
	    push @{$hash{$val}}, \%hR;
	    
	    
	}
	if(not defined $hash{$val}){
	    my %hR;
	    foreach $field (%$fieldHashRef){
		$hR{$field} = -1;
	    }
	    push @{$hash{$val}}, \%hR;
	}
	
    }
    
    
    
    return (\%hash);
    
}

# (\@results) = returnMultipleResultsForExSTH($sth);
# Conglomerate a general table query into a hash result - assumes sth has already been executed
# Inputs: $sth - the executed statement handle
# Outputs: \@results - an array of hash results. Hash has keys equal to fields in search.
# Kristen Naegle
# November 8, 2008
sub returnMultipleResultsForExSTH($){
    my ($sth, ) = @_;
    my @results;
    my $fieldHashRef;
    
    
    $fieldHashRef = returnFieldHashFromSTH($sth);
    
	
    while (my @row = $sth->fetchrow_array){
	
	my %hR;
	foreach my $field (keys %$fieldHashRef){
	    my $col = $fieldHashRef->{$field};
	    #print "$field maps to $col\n";
	    $hR{$field} = $row[$col];
	}
	push @results, \%hR;
	
	    
    }
    if(not @results){
	my %hR;
	foreach my $field (%$fieldHashRef){
	    $hR{$field} = -1;
	}
	push @results, \%hR;
    }
    
    
    return (\@results);
    
}

# \%fieldHash = returnFieldHashFromSTH($sth)
# Use this to translate between a column name and a column value for a $row[$column[ reference. Given an executed statement handle, return a hash with fields equal to the column names and values equal to the column value of row. 
# Inputs: $sth - executed statement handle
# Outputs: $fieldHash -ref. to hash with field names as keys and translation to column values as value
# Kristen Naegle
# Date ?
sub returnFieldHashFromSTH($){
    my ($sth) = @_;
    my %fieldHash;
    my $colNames = $sth->{NAME};
    for (my $i=0; $i < scalar(@$colNames); $i++){
	$fieldHash{$colNames->[$i]} = $i;
    }
    return  \%fieldHash;
}


# \%accHash = returnProtAccHash($dbh, $proteinIdArr);
# Return the accHash for all protein ids passed in, i.e. for every proteinId return the accession for GO fetch
# Inputs: $dbh - databae handle
#         $proteinIdArr - ref to array of protein Ids
# Outputs: \%accHash - keys are proteinId and values are the accession, chooses swissprot first and acc_gene second
# Kristen Naegle
# April 15, 2008
sub returnProtAccHash($$){
    my ($dbh, $proteinIdArr) = @_;
    my %hash;
    my $value = 'value';
    my $table = 'acc';
    my %argsHash;
    $argsHash{'protein_id'} = '?';
    $argsHash{'type'} = "'swissprot'";
    my $str = createSelectString($value, $table, \%argsHash);
    print "Command: $str\n";
    my $sth = $dbh->prepare($str);
    print "Protein ID:      \t Accession\n";
    my %geneArgs;
    $geneArgs{'id'} = '?';
    my $strGene = createSelectString('acc_gene', 'protein', \%geneArgs);
    my $sthGene = $dbh->prepare($strGene);
    foreach my $proteinId (@$proteinIdArr){
	if(not defined $hash{$proteinId}){
	    $sth->execute($proteinId);
	    my $acc = returnSingleResultOnCol($sth, 0);

	    if($acc eq '-1'){
		$sthGene->execute($proteinId);
		my $gene = returnSingleResultOnCol($sthGene, 0);
		#print "   $proteinId\t $gene\n";
		$hash{$proteinId} = $gene;
	    }
	    else{
		$hash{$proteinId} = $acc;
	    }
	}
    }
    return \%hash;
}

# \%accHash = returnGeneListForExperiment($dbh, $expId);
# Return the accHash for an experiment, i.e. for every proteinId return the accession for GO fetch
# Inputs: $dbh - databae handle
#         $expId - experiment Id
# Outputs: \%accHash - keys are proteinId and values are the accession, chooses swissprot first and acc_gene second
# Kristen Naegle
# April 15, 2008
# CLEARLY NOT WORKING
sub returnGeneListForExperiment($$){
    my ($dbh, $experiment_id) = @_;

    # foreach protein in experiment, get GO accessible accession
    my $proteinIds = returnProteinIdsForExpId($dbh, $experiment_id);
    my $accHash = returnProtAccHash($dbh, $proteinIds);
# #    foreach my $protein (@$proteinIds){
# #	print $pro
# #    }
}

#(\%hash) = returnMSInfoByExpId($dbh, $expId);
# return the trypsinized fragments and gi numbers for an experiment id 
# Inputs: $dbh - database handle
#         $expId - (FK) experiment Id in MS table and experiment table
# Outputs: Returns a hash of hashes. Top hash has MS_id as key and second hash has MS.phosphopep, acc.value, protein.name, protein.sequence, protein.species, phosphopep.aligned as keys
# Kristen Naegle
# April 28, 2008        
sub returnMSInfoByExpId($$){
    my ($dbh, $expId) = @_;
    my %hash;
   
    my $sth=$dbh->prepare('SELECT MS.id, MS.phosphopep, acc.value, protein.name, protein.sequence, protein.species, phosphopep.pep_aligned, protein.acc_gene, protein.id FROM MS join protein on protein.id=MS.protein_id join acc on acc.protein_id=protein.id join MS_phosphopep on MS_phosphopep.MS_id=MS.id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where experiment_id=? and acc.type=\'gi\'');
    
    $sth->execute($expId);
    my $i = 0;
    while(my @row= $sth->fetchrow_array){
	my %hashKey;
	$hashKey{'MS.id'} = $row[0];
	$hashKey{'MS.phosphopep'} = $row[1];
	$hashKey{'acc.value'} = $row[2];
	$hashKey{'protein.name'} = $row[3];
	$hashKey{'protein.sequence'} = $row[4];
	$hashKey{'protein.species'} = $row[5];
	$hashKey{'phosphopep.pep_aligned'} = $row[6];
	$hashKey{'protein.acc_gene'} = $row[7];
	$hashKey{'protein.id'} = $row[8];
	my $id = $i;
	$i+=1;
	$hash{$id} = \%hashKey;
    }
    
    return (\%hash);

}

# \%hash = returnPepsForMSIds($dbh, \@MS_ids)
# Return a hash (with keys MS_ids) and values are array of the singly phosphorylated, idealized trypsinizations for that MS.id fragment
# Inputs: $dbh - database handle
#         \@MS_ids - ref. to array of MS_ids 
# Outputs: \%hash - ref to hash, keys are MS_ids and values are array of trypsinized, singly phospho peptides
# Kristen Naegle
# May 2, 2008
sub returnPepsForMSIds($$){
    my ($dbh, $MS_ids) = @_;
    my %hash;
    my $sth = $dbh->prepare('SELECT phosphopep.pep_tryps FROM MS join MS_phosphopep on MS_phosphopep.MS_id=MS.id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS.id=?');
    foreach my $MS_id (@$MS_ids){
	$sth->execute($MS_id);
	my @tryps;
	while (my @row = $sth->fetchrow_array){
	    push @tryps, $row[0];
	}
	push @{$hash{$MS_id}}, @tryps;
    }
    return \%hash;

}

# \%hash = returnAllPepTrypsHits($dbh, $pep, $species)
# Returns a hash of all phosphopep_aligned, protein_id, gene, and names for those proteins which have a trypsinized fragment hit
# Assumes dbh has been loaded with other protein options for that tryps already
# Inputs: $dbh - database handle
#         $pep - the trypsinized peptide
#         $species - the species to limit proteins to
# Outputs: \%hash - keys are an index into a hit number and values are a hash with keys: phophsopep.pep_aligned, protein.id, protein.acc_gene, and protein.name
# Kristen Naegle
# May 1, 2008
sub returnAllPepTrypsHits($$$){
    my ($dbh, $pep, $species) = @_;
    $species = lc($species);
    my %hash;
    my $sth = $dbh->prepare('SELECT pep_aligned, protein.id, acc_gene, name from phosphopep join protein on protein.id=phosphopep.protein_id where pep_tryps = ? and protein.species=?');
    $sth->execute($pep, $species);
    my $i = 0;
    while (my @row = $sth->fetchrow_array){
	my %tHash;
	($tHash{'phosphopep.pep_aligned'}, $tHash{'protein.id'}, $tHash{'protein.acc_gene'}, $tHash{'protein.name'}) = @row;
	$hash{$i} = \%tHash;
	$i += 1;
    }
    return \%hash;
}

# $sth = returnDataForMSIdSTH($dbh);
# create a statement handle for data selection based on variable MS_id
# Inputs: $dbh - database handle
# Outputs: $sth - statement handle - selects data based on condition, time or cell (in type)
# Kristen Naegle
# May 15, 2008
sub returnDataForMSIdSTH($){
    my ($dbh) = @_;
#    my $sth = $dbh->prepare('SELECT type, run, label, value, NA FROM data where MS_id= ? and (type REGEXP \'time\' OR type REGEXP \'condition\' or type REGEXP \'cell\') order by priority');
    my $sth = $dbh->prepare('SELECT type, run, label, value, NA FROM data where MS_id=? order by run, priority');
    my %hash;
    $hash{'type'} = 0;
    $hash{'run'} = 1;
    $hash{'label'} = 2;
    $hash{'value'} = 3;
    $hash{'NA'} = 4;
    return ($sth, \%hash);


}

# (\@data, \@labels) = returnDataForMSId($sth, $MS_id, $cHash);
# Given a statement handle from returnDataForMSIdSTH, find the data vector for an MS_id. Return as hash based on different runs.
# Inputs: $sth - statement handle to find data for an MS_id given an order
#         $MS_id - MS id to execute query on
#         $cHash - another output of returnDataForMSIdSTH - the key value mapping for field namesin the executed statment
# Outputs: \@data - reference to array of data
#          \@labels - reference to array of labels that describe data (in order)
# Kristen Naegle
# May 15, 2008
sub returnDataForMSId($$$){
    my ($sth, $MS_id, $cHash) = @_;
    $sth->execute($MS_id);
    my @data;
    my @labels;
    my %data;
    my %labels;
    while(my @row = $sth->fetchrow_array){
	my $value = $row[$cHash->{'value'}];
	my $run = $row[$cHash->{'run'}];
	if(not defined $data{$run}){
	    $data{$run} = [];
	    $labels{$run} = [];
	}
	#if($value == 0){
	if($row[$cHash->{'NA'}]){
	    $value = 'NaN';
	}
	#}
	push @{$data{$run}}, $value;
	#replace underscores in label with dashes for matlab readability:
	my $labelValue = $row[$cHash->{'label'}];
	$labelValue =~ s/_/-/g;
	chomp $labelValue;
	$labelValue =~ s/\cM//; #this is for grandfathered weirdness prior to dos2unix requirement on data loading.
	my $type = $row[$cHash->{'type'}];
	$type =~ s/_/-/g;
			
#	push @{$labels{$run}}, $row[$cHash->{'type'}]."_".$row[$cHash->{'label'}];
	push @{$labels{$run}}, $type."_".$labelValue;
	#print "DEBUG: new label $type _ $labelValue\n";
			 
    }
    my @runs = keys %labels;
    @labels = @{$labels{$runs[0]}}; 
    return (\%data, \@labels);
}

# (\@data, \@labels) = returnDataForMSIdNoStddev($sth, $MS_id, $cHash);
# Given a statement handle from returnDataForMSIdSTH, find the data vector for an MS_id. Return as hash based on different runs.
# Inputs: $sth - statement handle to find data for an MS_id given an order
#         $MS_id - MS id to execute query on
#         $cHash - another output of returnDataForMSIdSTH - the key value mapping for field namesin the executed statment
# Outputs: \@data - reference to array of data
#          \@labels - reference to array of labels that describe data (in order)
# Kristen Naegle
# January 19, 2010 - mod from returnDataForMSId to do the same, but exclude stddev from consideration. For use in feature enrichment.
sub returnDataForMSIdNoStddev($$$){
    my ($sth, $MS_id, $cHash) = @_;
    $sth->execute($MS_id);
    my @data;
    my @labels;
    my %data;
    my %labels;
    while(my @row = $sth->fetchrow_array){
	my $value = $row[$cHash->{'value'}];
	my $run = $row[$cHash->{'run'}];
	if(not defined $data{$run}){
	    $data{$run} = [];
	    $labels{$run} = [];
	}
	#if($value == 0){
	if($row[$cHash->{'NA'}]){
	    $value = 'NaN';
	}
	#}

	#replace underscores in label with dashes for matlab readability:
	my $labelValue = $row[$cHash->{'label'}];
	$labelValue =~ s/_/-/g;
	chomp $labelValue;
	$labelValue =~ s/\cM//; #this is for grandfathered weirdness prior to dos2unix requirement on data loading.
	my $type = $row[$cHash->{'type'}];
	$type =~ s/_/-/g;
	if($type !~ m/stddev/){
#	push @{$labels{$run}}, $row[$cHash->{'type'}]."_".$row[$cHash->{'label'}];
	    push @{$labels{$run}}, $type."_".$labelValue;
	    push @{$data{$run}}, $value;
	    #print "DEBUG: new label $type _ $labelValue\n";
	}
    }
    my @runs = keys %labels;
    @labels = @{$labels{$runs[0]}}; 
    return (\%data, \@labels);
}



# \@pep_aligned = returnPepAlignedForMSids($dbh, $MSids)
# For an array of MSids, return an array of aligned peptides
# Inputs: $dbh - database handle
#         $MSids - ref to array of MS ids
# Outputs: \@pep_aligned - ref to array of aligned peptides
# June, 2008
# Kristen Naegle
sub returnPepAlignedForMSids($$){
    my ($dbh, $MSids) = @_;
    my @pep_aligned;
    my $sth = $dbh->prepare('SELECT phosphopep.pep_aligned FROM MS join MS_phosphopep on MS_phosphopep.MS_id=MS.id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS.id=?');
    foreach my $id (@$MSids){
	$sth->execute($id);
	my $peps = returnArrayOfResultsOnCol($sth, 0);
	push @pep_aligned, @$peps;
	 
    }
    return \@pep_aligned;
}

# \@fgnd_msIds = returnForegroundBasicMath($dbh, $MSidArr, $dataLabel_num, $dataLabel_denom, $value)
# Return a foreground from an array of MS ids based on whether the numerator data divided by the denominator data is greater than or equal to some value
# Inputs: $dbh - database handle
#         $MSidArr - background from which to pull MIS ids
#         $dataLabel_num - the term for the data label you wish to be numerator
#         $dataLabel_denom - the term for the data label you wish to be the denominator
#         $value - the cutoff value that it should be above
# Outputs: $fgnd_msIds - ref. to array of ms ids for which numerator/denominator >= value
# Kristen Naegle
# July 1, 2008
sub returnForegroundBasicMath($$$$$){
    my($dbh, $MSidArr, $dataLabel_num, $dataLabel_denom, $value) = @_;
    my @fgnd;
    my $sth = $dbh->prepare('SELECT value from data where MS_id=? and label=? order by priority');
    foreach my $msid (@$MSidArr){
	$sth->execute($msid, $dataLabel_num);
	my $num = returnSingleResultOnCol($sth, 0);
	if(not defined($num)){
	    next;
	}
	if($num == -1){
	    handleError('returnForegroundBasicMath', 'Value for dataLabel numerator did not exist, check your label value', \@_);
	    return \@fgnd;
	}
	$sth->execute($msid, $dataLabel_denom);
	my $denom = returnSingleResultOnCol($sth, 0);
	if(not defined($denom)){
	    next;
	}
	if($denom == -1){
	    handleError('returnForegroundBasicMath', "Value for dataLabel denominator did not exist for $msid, check your label value", \@_);
	    return \@fgnd;
	}
	if($denom == 0){
	    next;
	}
	my $ratio = $num/$denom;
	if($ratio >= $value){
	    push @fgnd, $msid;
	}
    }
    return \@fgnd;
}

# $proteinIdArr = returnProteinIdsWithDomainPhosphorylation($dbh, $domainName, \@species)
# This returns the protein ids that have known phosphorylations occuring in the specified domain and when @species is not empty, only for a particular species as specified
# Inputs: $dbh - database handle
#         $domainName - name of domain. This is an exact match, not a regexp
#         $speciesArr - array of species name
# Kristen Naegle
# July 10, 2008
# Mod Nov. 5th, 2008 to include only certain species
sub returnProteinIdsWithDomainPhosphorylation($$$){
    my ($dbh, $domainName, $speciesArr) = @_;
    my $sth;
    if(scalar(@$speciesArr)){
	my $speciesStr =  createSTHArrayVal($speciesArr);
	$sth =  $dbh->prepare("SELECT protein.id from protein join phosphopep on phosphopep.protein_id=protein.id where pfam_site=? and species IN $speciesStr group by protein.id");
    }
    else{
	$sth = $dbh->prepare('SELECT protein.id from protein join phosphopep on phosphopep.protein_id=protein.id where pfam_site=? group by protein.id');

    }
    $sth->execute($domainName);
    my $results = returnArrayOfResultsOnCol($sth, 0);
    return $results;

}



# $proteinIdArr = returnProteinIdsWithDomain($dbh, $domainName)
# This returns the protein ids that have known domain and contain phosphorylation sites on it (not necessarily in domain, but proteins exist in database not linked to known phosphorylation sites
# This will include all species
# Inputs: $dbh - database handle
#         $domainName - name of domain. This is an exact match, not a regexp
# Kristen Naegle
# July 16, 2008
sub returnProteinIdsWithDomain($$$){
    my ($dbh, $domainName, $pval) = @_;
    my $sth = $dbh->prepare('SELECT domain.protein_id from domain join phosphopep on phosphopep.protein_id=domain.protein_id where domain.label=? and p_value <= ? group by domain.protein_id');
    $sth->execute($domainName, $pval);
    my $results = returnArrayOfResultsOnCol($sth, 0);
    
    return $results;

}

# ($startArr, $stopArr) = returnDomainDetailsForProtein($dbh, $proteinId, $domainName)
# returns start and stop position of domain(s) in a protein. Returns -1 for both if domain wasn't found
# Inputs: $dbh - database handle
#         $proteinId - id of protein
#         $domainName - case insensitive, but text exact, name of domain
# Outputs: $startArr - ref. to array of start positions of each domain hit
#          $stopArr - the corresponding stop positions
# Kristen Naegle
# July 10, 2008
sub returnDomainDetailsForProtein($$$$){
    my ($dbh, $proteinId, $domainName, $pval) = @_;
    #print "PVALUE in Domiain: $pval\n";
    my $sth = $dbh->prepare('SELECT start, stop from domain where protein_id=? and p_value<=? and label=? order by start');
    $sth->execute($proteinId, $pval, $domainName);
    my $starts = returnArrayOfResultsOnCol($sth, 0);
    $sth->execute($proteinId, $pval, $domainName);
    my $stops = returnArrayOfResultsOnCol($sth, 1);
    return ($starts, $stops);
}


# ($acc_gene,$species, $site_typeArr, $site_posArr) = returnPhosphoSitesInProteinDomain($dbh, $proteinId, $domainName)
# Return all the phosphorylation sites in a given protein domain
# Inputs: $dbh - database handle
#         $proteinId - protein.id
#         $domainName - name of domain (exact name of domain)
# Outputs: $acc_gene - gene name of protein
#          $species - species 
#          $site_typeArr - array of site types
#          $site_posArr - corresponding array of site positions (rel. to start of protein)
# Kristen Naegle
# July 13, 2008
sub returnPhosphoSitesInProteinDomain($$$){
    my ($dbh, $proteinId, $domainName)= @_;
    
    my $sth = $dbh->prepare('SELECT site_pos, site_type FROM phosphopep join protein on protein.id=phosphopep.protein_id where protein_id=? and pfam_site=? order by site_pos');
    $sth->execute($proteinId, $domainName);
    my $site_posArr = returnArrayOfResultsOnCol($sth, 0);
    $sth->execute($proteinId, $domainName);
    my $site_typeArr = returnArrayOfResultsOnCol($sth, 1);
    my $sth_p = $dbh->prepare('SELECT acc_gene, species from protein where protein.id=?');
    $sth_p->execute($proteinId);
    my $acc_gene = returnSingleResultOnCol($sth_p, 0);
    $sth_p->execute($proteinId);
    my $species = returnSingleResultOnCol($sth_p, 1);
    return ($acc_gene,$species, $site_typeArr, $site_posArr);


}


# $gohash = returnGOHashForProteinId($dbh, $proteinId
# Returns a hash with keys 'F', 'C' and 'P' of GO terms and values (in a hash) for a protein Id
# Inputs: $dbh - database handle
#         $proteinId - protein table id
# Outputs: $hash - ref to hash with top keys 'F' 'C' and 'P' and values are arrays of hashes with keys 'GO' and 'term'
# Kristen Naegle
# July 18, 2008
sub returnGOHashForProteinId($$){
    my ($dbh, $proteinId) = @_;
    my $sth = $dbh->prepare('SELECT GO.GO, GO.term, GO.aspect from protein_GO join GO on protein_GO.GO_id=GO.id where protein_GO.protein_id = ?');
    $sth->execute($proteinId);
    my %hash;
    while(my @row = $sth->fetchrow_array()){
	my %GO;
	my $aspect = $row[2];
	
	$GO{'GO'} = $row[0];
	$GO{'term'} = $row[1];
	if(not defined $hash{$aspect}){
	    @{$hash{$aspect}} = \%GO;
	}
	else{
	    push @{$hash{$aspect}}, \%GO;
	}
    }
    return \%hash;

}

# $hashRef = returnAccHashByTypeFromArr($accArr)
# Converts an accession array into a hash with keys indicating accession type
# Inputs: $accArr - array of accessions, see returnAccValuesByProteinId
# Outputs: $hashRef - reference to hash with keys indicating type and values an array of accession of that type from array
# Kristen Naegle
# July 20, 2008
sub returnAccHashByTypeFromArr($){
    my ($accArr) = @_;
    my %hash;
    foreach my $acc (@$accArr){
	my $type = returnAccType($acc);
	if(not defined $hash{$type}){
	    @{$hash{$type}} = $acc;
	}
	else{
	    push @{$hash{$type}}, $acc;
	}
	
    }
    return \%hash;
}

# $proteinHash = returnProteinDescForProteinId($dbh, $proteinId);
# returns a hash, where keys are protein table entries and values 
# Inputs: $dbh - database handle
#         $proteinId - id to protein table
# Outputs: $proteinHash - hash with keys id, sequence, species, acc_gene, name, and date. If it didn't exist id value is -1
# Kristen Naegle
# July 20, 2008
sub returnProteinDescForProteinId($$){
    my ($dbh, $proteinId) = @_;
    
    my $sth=$dbh->prepare('SELECT id, sequence, species, acc_gene, name, date from protein where id=?');
    $sth->execute($proteinId);
    my @row = $sth->fetchrow_array;
    my %hash;
    if($row[0]){
	($hash{'id'}, $hash{'sequence'}, $hash{'species'}, $hash{'acc_gene'}, $hash{'name'}, $hash{'date'}) = @row;
    }
    else{
	$hash{'id'} = -1;
    }
    return \%hash;
}

# $proteinHash = returnPhosphopepPredForPhosphopepId($dbh, $phosphopep_id)
# returns a hash, where keys are protein table entries and values 
# Inputs: $dbh - database handle
#         $phosphopep_id - FK to pep_prediction, phosphopep id
# Outputs: $predHash - array of hashes with keys id, 
# Kristen Naegle
# November 6, 2008
sub returnPhosphopepPredForPhosphopepId($$){
    my ($dbh, $pepId) = @_;
    my @arr;
    my $sth=$dbh->prepare('SELECT id, source, value, score from phosphopep_prediction where phosphopep_id=?');
    $sth->execute($pepId);
    while (my @row = $sth->fetchrow_array){
	my %hash;
	if($row[0]){
	}
	else{
	    @row = (-1, -1, -1, -1);
	}
	($hash{'id'}, $hash{'source'}, $hash{'value'}, $hash{'score'}) = @row;
	push @arr, \%hash;
    }
	return \@arr;
}


# \@proteinIds = returnProteinsWithoutExpression($dbh)
# Returns all the protein ids that do not appear in protein_expression -- ONLY for mouse and human proteins
# Inputs: $dbh - database handle
# Outputs: $proteinIds - ref. to array of protein Ids
# Kristen Naegle
# August 12, 2008
# Jan 16, 2009 - updated to only look at human and mouse
sub returnProteinsWithoutExpression($){
    my ($dbh) = @_;
    my $sth = $dbh->prepare("SELECT protein.id from protein left outer join protein_expression on protein_expression.protein_id=protein.id WHERE protein_expression.protein_id IS NULL and (species='homo sapiens' or species='mus musculus')");
    $sth->execute();
    my $proteinIds = returnArrayOfResultsOnCol($sth, 0);
    return $proteinIds;
}

# $convHash = convertMSIdToDescription($dbh, $MSIdArr)
# Create a hash that points from the MS id to a gene_site1_site2..etc
# Inputs: $dbh - database handle
#         $MSIdArr - ref. to an array of MS Ids 
# Outputs: $convHash - ref. to hash with keys equal to MS id and value equal to name with gene and sites
# Kristen Naegle
# December 12, 2008
sub convertMSIdToDescription($$){
    my ($dbh, $MSIdArr) = @_;
    my %conv;
    my $sth = $dbh->prepare('SELECT acc_gene, site_type, site_pos from MS join protein on MS.protein_id=protein.id join MS_phosphopep on MS_phosphopep.MS_id=MS.id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS.id=?');
    foreach my $MSId (@$MSIdArr){
	$sth->execute($MSId);
	
	my $result = returnMultipleResultsForExSTH($sth);
	my $name = ($result->[0])->{'acc_gene'};
	foreach my $r (@$result){
	    my $site_type = $r->{'site_type'};
	    my $site_pos = $r->{'site_pos'};
	    $name .= "_".$site_type.$site_pos;
	}
	$conv{$MSId} = $name;
    }
    
    return \%conv;
}

# $hash = returnMSIdInfoHash($dbh, $MSId)
# Given an MS.id, return a hash struct with each of the terms GO.F, GO.C, GO.P scansite_kinase, scansite_bind, domains, and site_domain
# Returns unique domains per a protein. -- think that's correct.  Same with pfam_site, if multiple sites in a peptide, returns only one domain if redundant, otherwise both
# Inputs: $dbh - database handle
#         $MSId - id to MS Table
# Outputs: $hash -ref to hash with fields GO.P, GO.C, GO.F, pfam_site, scansite_bind, scansite_kinase and that points a hash with keys of corresponding labels
# Kristen Naegle
# January 28, 2009
sub returnMSIdInfoHash($$){
    my ($dbh, $MSId) = @_;
    my %hash;
    my $sth;
    my $results;
    # Get GO Terms 
    $sth=$dbh->prepare('SELECT * from MS join protein_GO on MS.protein_id=protein_GO.protein_id join GO on protein_GO.GO_id=GO.id where MS.id=? and aspect=?');
    my @aspects = ('F', 'C', 'P');
    foreach my $aspect (@aspects){
	$sth->execute($MSId, $aspect);
	$results =  returnMultipleResultsForExSTH($sth);
	my $field = $aspect;
	$hash{$field} = reshapeHashFromResults($results, 'term', 'GO');
    }
    # Get Domain Terms
    $sth=$dbh->prepare('SELECT * from MS join domain on domain.protein_id=MS.protein_id where MS.id=?');
    $sth->execute($MSId);
    $results = returnMultipleResultsForExSTH($sth);
    $hash{'domains'} = reshapeHashFromResults($results, 'label', 'p_value');
    
# Get Scansite Terms
    $sth=$dbh->prepare('SELECT * from MS join MS_phosphopep on MS.id=MS_phosphopep.MS_id join phosphopep_prediction on MS_phosphopep.phosphopep_id=phosphopep_prediction.phosphopep_id where MS.id=? and source=?');
    my @sources = ('scansite_kinase', 'scansite_bind', 'pelm_kinase');
    foreach my $source (@sources){
	my $s = $source; #"scansite_".$source;
	$sth->execute($MSId, $s);
	$results = returnMultipleResultsForExSTH($sth);
	$hash{$s} = reshapeHashFromResults($results, 'value', 'score');

    }
 
   # Get Pfam_site terms
    $sth=$dbh->prepare('SELECT * from MS join MS_phosphopep on MS_phosphopep.MS_id=MS.id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS.id=?');
    $sth->execute($MSId);
    $results = returnMultipleResultsForExSTH($sth);
    $hash{'pfam_site'} = reshapeHashFromResults($results, 'pfam_site', 'site_pos');
    


    return \%hash;
}

# $hashReshaped = reshapeHashFromResults($results, $keyField, $valField);
# given an array of hashes from an executed sth, push this into a new hash, reshaping by the key col and the val col - ignore -1 
# Inputs: $results - array of results of hashes from an executed sth for example
#        $keyField - the field of the hashes for whcih you want to be the new key
#        $valField - the corresponding value field that will be set in the new hash
# Kristen Naegle
# January 28, 2009
sub reshapeHashFromResults($$$){
    my ($results, $keyField, $valField) = @_;
    my %hash;
    
    foreach my $result (@$results){
	my $key = $result->{$keyField};
	if(not defined ($hash{$key}) and $key ne '-1'){
	    $hash{$key} = $result->{$valField};
	}
	
    } 
    return \%hash;

}

# $MSIdHash =  returnMSIdsWithMiscleavedForms($dbh, $expId)
# For an experiment id, finds all MS ids that have a cleaved form and a mathing miscleaved form.  See MatlabIO/plotMiscleavages for examples
# Inputs: $dbh - database handle
#         $expId - id to experiment table
# Outputs: $MSIdHash - hash with hit MS.id (the cleaved form) as key and an array of miscleaved MS Ids
# Kristen Naegle
# Feb. 26, 2009
sub returnMSIdsWithMiscleavedForms($$){
    my ($dbh, $expId) = @_;
    my $MSids = returnMSIdsForExpId($dbh, $expId);
#get the phosphopeps from MS 
    my %MSIdHash; 
    foreach my $MSId (@$MSids){ # for every MS id find if there is more than one match
	my $sth = $dbh->prepare('select phosphopep from MS where MS.id=?');
	$sth->execute($MSId);
	my $phosphopep = returnSingleResultOnCol($sth, 0);
	my $sthHit = $dbh->prepare('select MS.* from MS where experiment_id=? and MS.phosphopep REGEXP binary ?');
	$sthHit->execute($expId, $phosphopep);
	my $results = returnMultipleResultsForExSTH($sthHit);
	if(scalar(@$results) > 1){
	    
#	print "$MSId\t$phosphopep\n";
	    my @rs;
	    foreach my $result (@$results){
		#print "$result->{'id'}\t$result->{'phosphopep'}\n";
		if($result->{'id'} != $MSId){
		    push @rs, $result->{'id'};
		}
	    }
     
	    push @{$MSIdHash{$MSId}}, @rs;
	}
	
    }
    return \%MSIdHash;
    
}

# ($proteinId, $accExists) = findProteinByRichSeq($dbh, $acc, $richSeq)
# Determine whether a protein exists in the database, looks by accession types first, then a sequence and species to speed up search
# Inputs: $dbh - database handle
#         $acc - the accession used to retrieve the richSeq record
#         $richSeq - rich seq object of retrieved record
# Outputs: $proteinId - id of protein table that matches , -1 if record DNE
#          $accExists - 1 if accession already exists and 0 otherwise
# Kristen Naegle
# May 13, 2009
sub findProteinByRichSeq($$$){
    my ($dbh, $acc, $richSeq) = @_;
    my ($proteinId, $accExists);
    $accExists = 0;
    
    $proteinId = returnProteinIdByAcc($dbh, $acc); # first look by accession. 
    if($proteinId != -1){
	$accExists = 1;
	print "DEBUG: Found protein based on accession\n";
	return($proteinId, $accExists);
    }
    
    # doesn't exist by the direct accession, need to look by other fields
    my ($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc) = getProteinFromRichSeq($richSeq);

#     # try again by the $primaryAcc that was returned - need to make sure it's not an empty field - see if it matches $acc
#     # Turns out this is bad...when a primaryAcc may not contain isoform specific information... removed Dec. 8, 2009
#   #   if($primaryAcc){
# # 	if($primaryAcc ne $acc){
# # 	    $proteinId = returnProteinIdByAcc($dbh, $primaryAcc);
# # 	    if($proteinId != -1){
# # 		print "DEBUG: Found protein based on rich seq primary accession, different than accession passed in acc:$acc primary:$primaryAcc\n";
# # 		$accExists = 0; #protein exists, but by a different accession
# # 		return($proteinId, $accExists);
# # 	    }
# # 	}
# #     }

    #next try to look up based on gene name and geneSynonyms - limit query on species - added species as an index in protein table to speed it up
    $proteinId = findProteinByGenes($dbh, $species, $gene, $geneSynonyms, $sequence);
    if($proteinId != -1){
	$accExists = 0;
	print "DEBUG: Found protein based on gene name\n";
	return($proteinId, $accExists);
    }
    
    # if that didn't work, finally do the expensive operation and see if a protein of the exact sequence matches for that species.   -  this is the final search
    $proteinId = returnProteinIdBySeqSpecies($dbh, $sequence, $species);
    $accExists = 0;
    if($proteinId > 0){ print "DEBUG: Found protein based on sequence\n";}
    return($proteinId, $accExists);

}


#$proteinId = findProteinByGenes($dbh, $species, $gene, $geneSynonyms, $sequence);
# Find a protein based on acc_gene and geneSynonyms - limit by species
# Speeds up search by returning all the records by genes, and cross-referencing the exact sequence for isoform matches.  Adds $gene to the gene_synonym array in case 
# Inputs: $dbh - database handle 
#         $species - species of protein table
#         $ gene - acc_gene of gene table - but will be considered as synonym as well
#         $geneSynonyms - ref to an array of gene synonyms
#         $sequence - amino acid sequence
# Outputs: $proteinId - id of protein if matches, -1 if no match
# May 12th, 2009
sub findProteinByGenes($$$$$){
    my($dbh, $species, $gene, $geneSynonyms, $sequence) = @_;
    my $proteinId=-1; 
    
    my @newGeneSyn;
    push @newGeneSyn, @$geneSynonyms;
    push @newGeneSyn, $gene;  # in case gene is a geneSynonym..create new gene Syn
    # Look up first by acc_gene - if allows for empty acc_gene just in case..
    if($gene){
	my $proteinIdsGene = returnProteinIdsByGeneSpecies($dbh, $gene, $species);
    #if this isn't empty, do any of them match specific isofrom by sequence?

	$proteinId = matchProteinIdsBySequence($dbh, $proteinIdsGene, $sequence);
    }
    # now repeat for proteinIdsByGeneSynonym- allows for empty geneSynonmys
    if($proteinId == -1 and ($geneSynonyms or $gene)){
	print "DEBUG: gene: $gene and syns: @$geneSynonyms\n";
	my $proteinIdsByGeneSyn = returnProteinIdsByGeneSynSpecies($dbh, \@newGeneSyn, $species);
	
	$proteinId = matchProteinIdsBySequence($dbh, $proteinIdsByGeneSyn, $sequence);
    }

    return $proteinId;
}

# $proteinId = matchProteinIdsBySequence($dbh, $proteinIds, $sequence);
# for an array of protein ids find the one that has an exact sequence (i.e. isoform) match
# assumes proteinIds are only for records of the correct gene and species.
# Inputs: $dbh - database handle
#         $proteinIds - an array of protein ids to check sequences for (again, assumes they are already matched for accession and species)
#         $sequence - amino acid sequence of protein to match
# Outputs: $proteinId - id of matching protein, -1 if no match
# Kristen Naegle
# March 12th, 2009
sub matchProteinIdsBySequence($$$){
    my($dbh, $proteinIds, $sequence) = @_;
    my $proteinId = -1;
    if(!scalar(@$proteinIds)){
	return $proteinId;
    }

    my $idArrStr = makeSelectINString($proteinIds, 0);

#    print "DEBUG: ID ARR STR: $idArrStr\n";
    my $sth = $dbh->prepare("SELECT id from protein where id IN $idArrStr and sequence = ?");
    $sth->execute($sequence);
    my $ids = returnArrayOfResultsOnCol($sth, 0);
    if(scalar(@$ids) > 1){
	print "ERROR IN DATABASE: You have multiple protein entries for the same sequence, found with protein ids: @$proteinIds\n";
    }
    $proteinId = $ids->[0];
    return $proteinId;
}

# $speciesArr = returnArrayOfSpeciesFromProteinIds($dbh, $proteinIdArr);
# REturn the array of unique species in an array of protein ids
# For use in GO Tools - in order to load correct parser hash
# Inputs: $dbh - database handle
#         $proteinIdArr - ref. to an array of protein ids
# Outputs: $speciesArr - ref to array of resulting species that represents protein id arr
# Kristen Naegle
# May 21, 2009
sub returnArrayOfSpeciesFromProteinIds($$){
    my ($dbh, $proteinIds) = @_;
    my $proteinStr = makeSelectINString($proteinIds, 0);
    my $sth = $dbh->prepare("SELECT species from protein where id IN $proteinStr group by species");
    $sth->execute();
    my $speciesArr = returnArrayOfResultsOnCol($sth, 0);
    return $speciesArr;
}

# $id = getLastEntryId($dbh, $table)
# Get the last id entry in a table
# Inputs; $dbh - database handle
#         $table - name of table
# OUtputs: $id - id of last table entry
# Kristen Naegle
# May 26, 2009
sub getLastEntryId($$){
    my ($dbh, $table) = @_;
    
    my $sth = $dbh->prepare("SELECT id from $table order by id desc limit 1");
    $sth->execute();
    my $id = returnSingleResultOnCol($sth,0);
    return $id;
}

# $entries = getEntriesBetweenIds($dbh, $table, $start, $stop);
# Return actual ids that exist in table between $start and $stop ids
# Inputs: $dbh - database handle
#         $table - name of table
#         $start - first id to include
#         $stop - last id to include
# Outputs: $entries - ref. to array of ids that exist in table
# Kristen Naegle
# May 26, 2009
sub getEntriesBetweenIds($$$$){
    my ($dbh, $table, $start, $stop) = @_;

    my @arr = ($start..$stop);
    #print "@arr\n";
    my $selectStr = makeSelectINString(\@arr, 0);
    my $sth = $dbh->prepare("SELECT id from $table where id IN $selectStr");
    $sth->execute();
    my $entries = returnArrayOfResultsOnCol($sth, 0);
    return $entries;

}

# $domainHash = returnDomainDescHashForProteinId($dbh, $proteinId)
# return a hash of domains (with keys equal to the id of the domain) and values equal to a ref to a hash with keys 'id' 'start' 'stop' 'label' 'p_value'
# Inputs: $dbh - database handle
#         $proteinId - foreign key in domain to protein table
#         $key - name to use as reference key (for example 'id' or 'start')
# Outputs: $domainHash - top key is id and bottom kesy are table labels
# Kristen Naegle
# June 6, 2009
sub returnDomainDescHashForProteinId($$$){
    my ($dbh, $proteinId, $key) = @_;
    
    my $sth = $dbh->prepare('SELECT * from domain where protein_id=?');
    $sth->execute($proteinId);
    my $results = returnMultipleResultsForExSTH($sth);
    my %hash;
    foreach my $result (@$results){
	my $keyId = $result->{$key};
	if(defined $hash{$keyId}){
	    print "ERROR: key $key has multiple identical entries for value $keyId\n";
	}
	$hash{$keyId} = $result;
    }
    return \%hash;
}

# ($overlapHash, $domainIdsToRemove) = checkOverlappingDomains($dbh, $proteinId)
# Given a domainDescHash for a single protein, determine if there are overlaps (this is because I need to handle the case when one domain is completely covered by another
#  Inputs: $dbh -database handle
#          $proteinId - id of protein table FK in domain
# Outputs: $overlapHash - hash, with keys equal to a number of sets of domains and points to an array of domain description hashes that overlap eachother
#          $domainIdsToRemove - these are the domains to remove based on the fact that they overlap with another domain with a lower p-value
# Kristen Naegle
# June 6, 2009
sub checkOverlappingDomains($$){
    my ($dbh, $proteinId) = @_;

    my $domainHash = returnDomainDescHashForProteinId($dbh, $proteinId, 'start');
    my %hashOverlap; # keep an hash of things that overlap - each is an array 
    
    my @starts = keys %$domainHash;
    my @sorted = sort {$a <=> $b} @starts;
    my $count = 0;
    for (my $i=0; $i<$#sorted-1; $i++){
	my $hash1 = $domainHash->{$sorted[$i]};
	my $hash2 = $domainHash->{$sorted[$i+1]};
	my $stop = $hash1->{'stop'};
	my $start = $hash2->{'start'};
	if (($start - $stop) < 0){
	    my @arr; 
	    push @arr, $hash1;
	    push @arr, $hash2;
	    $hashOverlap{$count} = \@arr;
	    $count++;
	#    print "ProteinId: $proteinId\n";
	#    print "ERROR: Domains overlap: $hash1->{'label'} and $hash2->{'label'} with ids $hash1->{'id'} and $hash2->{'id'}\n\n";
	}
 
    }

    #return the domain id that should be kept for each overlap 
    my @domainsToRemove;
    foreach my $key (keys %hashOverlap){
	my $arr = $hashOverlap{$key};
	my $pvalIndex = 0;

	my $currHash = $arr->[$pvalIndex];
	my $pvalLow = $currHash->{'p_value'};
	#   print "DEBUG: Low pval to start: $pvalLow\n";
	for(my $i=1; $i <= scalar(@$arr)-1; $i++){
	    $currHash = $arr->[$i];
	    my $pval = $currHash->{'p_value'};
	    if($pval < $pvalLow){
		$pvalIndex = $i;
	    }
	    
	}
	#keep $pvalIndex for $key - push al
	for(my $i=0; $i<=scalar(@$arr); $i++){
	    if($i != $pvalIndex){
		my $h = $arr->[$i];
		push @domainsToRemove,$h->{'id'}; 
	    }
	}

	
    }
    return (\%hashOverlap, \@domainsToRemove);
}

# \@proteinIds = findFragmentsFromProteins($dbh, $peptide, $species);
# Given a peptide fragment and a species, find other proteins in the database that hit that fragment (for use in ambiguity)
# Inputs: $dbh - database handle
#         $peptide - fragment peptide
#         $species - binomial name of species
# Outputs: $proteinIds - ref. to array of protein ids (of the given species) that have the protein fragment in them
# Kristen Naegle
# June 24, 2009
sub findFragmentsFromProteins($$$){
    my ($dbh, $peptide, $species) = @_;
    
    my $sth = $dbh->prepare('SELECT id from protein where sequence REGEXP ? and species=?');
    $sth->execute($peptide, $species);
    my $proteinIds = returnArrayOfResultsOnCol($sth,0);
    return $proteinIds;

    
}

# \@genes = returnGeneListForProteinId($dbh, $proteinId);
# Given a protein id return all genes (acc_gene and all gene_synonyms). If acc_gene is empty than this is not put on that list
# Inputs: $dbh - database handle
#         $proteinId - id of protein table
# Outputs: $genes - ref to array of genes
# Kristen Naegle
# June 26, 2009
sub returnGeneListForProteinId($$){
    my ($dbh, $proteinId) = @_;

    my @genes;

    # get acc_gene 
    my $sth = $dbh->prepare('SELECT acc_gene from protein where id=?');
    $sth->execute($proteinId);
    my $acc_gene = returnSingleResultOnCol($sth, 0);
    if($acc_gene and $acc_gene ne '-1'){
	push @genes, $acc_gene;
    }
    #get gene synonyms
    $sth=$dbh->prepare('SELECT value from acc where protein_id=? and type=\'gene_synonym\'');
    $sth->execute($proteinId);
    my $results = returnArrayOfResultsOnCol($sth, 0);
    if($results->[0] ne '-1'){
	push @genes,@$results;
    }
    return \@genes;
}

# %transHash = translateDBHMSIds($dbh_old, $dbh_new, $expId_old, $expId_new)
# Given two databases and two corresponding experiment ids, return a translation hash from new database MS ids to old database MS ids
# Inputs: $dbh_old - the old database (translate from these)
#         $dbh_new - the new database handle (translate to these)
#         $expId_old - experiment id corresponding to old dbh
#         $expId_new - experiment id corresponding to new dbh
# Outputs: $transHash - ref. to hash with keys equal to old dbh msIds and values equal to new dbh MSIds.
# Kristen Naegle
# July 23, 2009
sub translateDBHMSIds($$$$){
    my ($dbh_old, $dbh_new, $expId_old, $expId_new) = @_;
    #translate an experiment, compare phosphopeps restricting to MS ids.
    my $msIds_new =  returnMSIdsForExpId($dbh_new, $expId_new);
    my %trans;
    foreach my $msId_new (@$msIds_new){
	#get the phosphopeptide
	my ($phosphopep, $proteinId) = returnMSPhosphopepProteinIdForMSId($dbh_new, $msId_new);
    my $map_msId = returnMSIdByPhosphopepExperimentId($dbh_old, $phosphopep,$expId_old);
    $trans{$map_msId} = $msId_new;
    }
    return \%trans;

}

# $logFileArr = returnLogFiles($dbh)
# Returns all current log files used by database
# Inputs: $dbh - database handle
# Outputs: $logFileArr - ref to array of log files. Does not include LOG_PATH which you can find in globalVars.pm
# Kristen Naegle
# November 15, 2009
sub returnLogFiles($){
    my($dbh) = @_;

    my $sth = $dbh->prepare('SELECT errorLog from experiment');
    $sth->execute();
    my $logFiles = returnArrayOfResultsOnCol($sth, 0);
    return $logFiles;

}

# $dataFileArr = returnDataFiles($dbh)
# Returns all current data files used by database
# Inputs: $dbh - database handle
# Outputs: $dataFileArr - ref to array of data files. Does not include DATA_PATH which you can find in globalVars.pm
# Kristen Naegle
# November 15, 2009
sub returnDataFiles($){
    my($dbh) = @_;

    my $sth = $dbh->prepare('SELECT dataset from experiment');
    $sth->execute();
    my $logFiles = returnArrayOfResultsOnCol($sth, 0);
    return $logFiles;

}


# $runArray = returnRunArray($dbh, $expId);
# Return the unique array of runs in an experiment from the data column. When runs don't appear equally, this will return it according to descending order of number of datapoints in a run
# Inputs: $dbh - database handle
#         $expId - id of experiment table under consideration
# outputs: $runArray - ref. to array of run names.
# Kristen Naegle
# Jan 12, 2009
sub returnRunArray($$){
    my($dbh, $expId) = @_;
    
    my @runs;
    my $sth = $dbh->prepare('SELECT run, count(*) as NUM from MS join data on MS.id=data.MS_id where experiment_id=? group by run order by NUM desc');
    $sth->execute($expId);
    my $results = returnMultipleResultsForExSTH($sth);
    for(my $i=0; $i<scalar(@$results); $i++){
	my $h = $results->[$i];
	push @runs, $h->{'run'};
    }
    return \@runs;
}

# \@expIds = returnUnpublishedExpIds($dbh);
# Find all unpubished experiments
# Inputs: $dbh - database handle
# Outputs: $expIds -ref to array of unpublished experiment ids
# Kristen Naegle
# January 28, 2010
sub returnUnpublishedExpIds($){
    my($dbh) = @_;
    
    my $sth = $dbh->prepare('SELECT id from experiment where published=0');
    $sth->execute();
    my $ids = returnArrayOfResultsOnCol($sth,0);
    return $ids;
}

# $uniprot = returnAccByTypeForMSId($dbh, $type, $MSid)
# Get the uniprot number for an MSId 
# Inputs: $dbh - database handle
#         $type - type of accession. e.g. 'gi', 'refseq', 'uniprot'
#         $MSid - ms.id entry
# Outputs: $uniprot - uniprot number
# Kristen Naegle
# January 28, 2010
sub returnAccByTypeForMSId($$$){
    my($dbh, $type, $msId) = @_;

    my $sth = $dbh->prepare('Select value from acc join MS on MS.protein_id=acc.protein_id where type=? and MS.id=?');
    $sth->execute($type, $msId);
    my $accs = returnArrayOfResultsOnCol($sth, 0);
    return $accs->[0];
}


1;
