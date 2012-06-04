use strict;
use warnings;
use DBI;
use DBTools::insertFunctions;
use DBTools::queryTools;
use DBTools::deleteFunctions;
use DBTools::maintenance;
use entrezTools;
use phosphopepTools;
use HMMTools;
use commonTools;
use GOTools;
use scansiteParser;
use expressionTools;


# ($proteinId, $iProteinFlag, $iAccFlag) = handleProteinsByRichSeq($dbh, $acc, $richSeq, $date);
# Handles insertion of protein by:
# 1. Check to see if exact accession exists (if so, return proteinId for it)
# 2. If not, see if protein exists based on sequence (get from entrez) and gene (if it does return proteinId)
# 3. If not, add protein AND Accession entry AND Pfam domain predictions
# Inputs: $dbh - database handler
#         $acc - an accession value
# Outputs: $proteinId - the id for the protein (that either existed or was inserted)
#          $iProteinFlag - flag, 1 if protein was inserted, 0 otherwise
#          $iAccFlag - flag, 1 if accession was inserted, 0 otherwise
# Kristen Naegle
# March 8, 2008
sub handleProteinByRichSeq($$$$){
    my($dbh, $acc, $richSeq, $date) = @_;
    my($proteinId, $iProteinFlag, $iAccFlag, $accId, $accExists);
    $iProteinFlag = 0;
    $iAccFlag = 0;
    my ($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newGI);
    
    ($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newGI) = getProteinFromRichSeq_requireAll($richSeq);
    my @giOld;
    if($errorCode){
	#   print "Error: Accession retrieval had error code $errorCode\n";
	my $t = translateAccGetErrorCode($errorCode);
	print "Error Code ($errorCode) means: $t\n";
	handleError('handleProteinByAcc', $t, \@_);
	if($newGI){
	    push @giOld, $acc;
	    handleError('handleProteinByRichSeq', "REPLACING incoming accession $acc with new accession $newGI\n", \@_);
	  #  print "DEBUG: Found new gi replacement record $newGI\n";
	    #get new richseq and check that there are no errorCodes with it..only attempt one new record for now
	    $richSeq = returnRichSeqFromAcc($newGI);
	    # if($richSeq == -1) {
# 		#RETRY 
# 	    }
	    
	    my $newNewGI; 
	    ($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newNewGI) = getProteinFromRichSeq_requireAll($richSeq);
	    push @giOld, $newGI;
	    while($newNewGI){
		push @giOld, $newNewGI;
		$richSeq = returnRichSeqFromAcc($newNewGI);
		($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc, $newNewGI) = getProteinFromRichSeq_requireAll($richSeq);
		
	    }
	    
	    print "DEBUG: Accessions we followed for a record: \n@giOld\n";
	    printLineToLog( "Accessions we followed for a record starting with $acc: @giOld\n");
	    if($errorCode){
		print "DEBUG: ERORR CODE FOR new gi $newGI is $errorCode\n";
		my $t = translateAccGetErrorCode($errorCode);
		print "Error Code ($errorCode) means: $t\n";
		handleError('handleProteinByAcc', $t, \@_);
		return(-1,0,0);
	    } #else we break and continue on as normal.
	    
	    #if we make to here than protein 
	}
	else{
	    return(-1, 0, 0); #means there was a WARNING that was not a new record
	}
    }
    # if error code is a warning code handle here - if it's a replacement record, get new richSeq
    
    if(@giOld){
	$acc = $giOld[$#giOld]; #it's the newest accession we care about looking for
    }
    ($proteinId, $accExists) = findProteinByRichSeq($dbh, $acc, $richSeq);  #tested - this works
    # if protein id is -1 then insert
    if($proteinId == -1){
	$proteinId = insertProtein($dbh, $sequence, $species, $gene, $name, $date);
	$iProteinFlag = 1;
    }
    if($accExists and $proteinId == -1){
	print "ERROR: Found a problem where I say acc exists but protein doesn't!!! Acc: $acc\n";
    }
    if(!$accExists and $proteinId != -1){
	($accId, $iAccFlag) = handleAcc($dbh, $acc, $proteinId);
	
    }



    # regardless of what was inserted before, handle gene, geneSyn, and primaryAcc;


    if($proteinId != -1){
	# new as of Dec. 4, 2009 - insert all other accessions. 
	my $globalAccTypes = $globalVars::ACC_TYPES;
	my $alternateAccList = returnValidDBAccList($richSeq, $globalAccTypes);
	foreach my $accToInsert (@$alternateAccList){
	    my($accAddtl, $iaccAddlt) = handleAcc($dbh, $accToInsert, $proteinId);
	}

	if($primaryAcc){
	    if(!($primaryAcc eq $acc)){
		my ($accIdP, $iAccFlagP)= handleAcc($dbh, $primaryAcc, $proteinId);
	    }
	    if(@giOld){ #have to insert old gi's so that files can find the proteins once they've been loaded
		foreach my $giOld (@giOld){
		    my ($accIdgi, $iAccFlaggi)= handleAcc($dbh, $giOld, $proteinId, 1);
		}
	    }
	}

	my @genes;
	push @genes, @$geneSynonyms;
	push @genes, $gene;
	print "Inserting Gene Synonyms: @$geneSynonyms\n";
	my $geneAccIds = handleGeneSynonymsForProteinId($dbh, \@genes, $proteinId);  #modify this to send in $gene + $geneSynonyms and make sure it never inserts a geneSyn that's the same as acc_gene
    }



    
    if($iProteinFlag and !$iAccFlag){
	print "ERROR! Protein Id:$proteinId inserted without Accession: $acc\n";
	handleError('handleProteinByRichSeq', "Protein: $proteinId inserted without accession: $acc", \@_);
#	exit;

    }

    return($proteinId, $iProteinFlag, $iAccFlag);
}


# ($accId, $iAccFlag) = handleAcc($dbh, $acc, $proteinId)
# Handles an accession for a protein_id. If accession already exists, returns accId, if it does not then it adds it to database. If undefined type, does not insert
# Inputs: $dbh -database handle
#         $acc - accession number, type derived from this
#         $proteinId (FK) - protein_id in accession table
# Outputs: $accId - id in acc table. Returns -1 if there was an error
#          $iAccFlag - 1 if inserted, 0 else
# Kristen Naegle
# March 27, 2008
# Modified Dec. 4, 2009 - to handle extra optional input variable indicating an accession is out of date 
sub handleAcc(){

    my ($dbh, $acc, $proteinId, $OLD);
    if(scalar(@_) < 3){
	print "FATAL ERROR: Call to handleACC without enough args\n";
    }
    elsif(scalar(@_) == 3){
	($dbh, $acc, $proteinId) = @_;
	$OLD = 0;
    }
    elsif(scalar(@_) == 4){
	($dbh, $acc, $proteinId, $OLD) = @_;
    }
    else{
	print "FATAL ERROR: Call to handleACC without enough args\n";
    }
    my $accId=-1;
    my $iAccFlag = 0;
    my $type = returnAccType($acc);
    if($proteinId != -1){
	if($type eq 'undefined'){
	    handleError('handleAcc', "Accession $acc has undefined type", \@_);
	    $accId = -1;
	}
	else{
	    $accId = returnAccIdByValueProteinId($dbh, $acc, $proteinId);
	    
	    if($accId != -1) {
		return($accId, $iAccFlag);
	    }
	    $accId = insertAcc($dbh, $type, $acc, $proteinId, $OLD);
	    $iAccFlag = 1;
	}
    }
    
    return ($accId, $iAccFlag);
}

# \@accIds =  handleGeneSynonymsForProteinId($dbh, $geneSynonyms, $proteinId);
# Insert an array of gene synonyms for a protein Id, makes sure that gene_synonym is not acc_gene. (insertAcc checks to make sure it doesn't exist)
# Inputs: $dbh - database handle
#         $geneSynonyms - array of gene synonyms to insert
#         $proteinId - id to protein table (FK in acc)
# Outputs: \@accIds - ref. to array of accession Ids. 
# Kristen Naegle
# Feb. 9, 2009
sub handleGeneSynonymsForProteinId($$$){
    my ($dbh, $geneSynonyms, $proteinId) = @_;
    my @accIds;
    
    #get acc_gene for proteinId 

    my $proteinHash = returnProteinDescForProteinId($dbh, $proteinId);
    my $acc_gene = $proteinHash->{'acc_gene'};

    foreach my $gene_syn (@$geneSynonyms){
	if($gene_syn ne $acc_gene and $gene_syn ne ''){
	   # print "Inserting $gene_syn\n";
	    my $accId = insertAcc($dbh, 'gene_synonym', $gene_syn, $proteinId, 0);
	    push @accIds, $accId;
	}
    }
    return \@accIds;
}

# ($errorCode, $domainIds) = handleDomainsForProteinId($dbh, $proteinId, $acc, $pvalueCutoff)
# Inserts pfam domain predictions for a proteinID (only if they don't already exist)
# Inputs: $dbh - database handler
#         $proteinId - id key to protein table
#         $pvalueCutoff - minimum pvalue required for accepting prediction
#         $CPU - number of processors to use for computation
#         $pfamDomainHash - see parse pfam file
# Outputs: $errorCode - nonzero if there was a problem with domain insertion
#          $domainIds - array of domain ids that were inserted for protein, unless error, this should never be empty (default for no predictions, the whole protein is set to "~~~"
# Kristen Naegle
# March 10, 2008
# Updated June 6, 2009 - since there is an error in overlapping domains that are completely incorporated
# Updated June 19, 2009 - to include look for parsed domains instead
sub handleDomainsForProteinId($$$$$){
    my ($dbh, $proteinId, $pvalueCutoff, $CPU, $domainHash) = @_;
    #my ($errorFlag, $dHashRef, $source, $params);
    my $domainIds;
    my $COMPUTE = 0;
    my $errorFlag = 0;
    my $domainsExist = returnDomainIdsForProteinId($dbh, $proteinId);
    if($domainsExist->[0] != -1){
	$domainIds = $domainsExist;
# #	$errorCode = 3; # domains are there, but don't want to flag problem. but print error to file, so we can see there was a problem
	handleWarning('handleDomainsForProteinId', "Domains already existed for protein id=$proteinId", \@_);
	print "DEBUG: Domains already existed for protein id $proteinId\n";
    }
    else{
	my ($label, $start, $stop, $p_value, $pfam_id, $description, $errorParsed, $source, $params, $version);


	# get the accession here instead
	my $accArr = returnAccArrOfTypeForProtein($dbh, $proteinId, 'swissprot');
	if($accArr->[0] eq '-1'){
	    $errorParsed = 1;
	}
	else{
	    foreach my $a (@$accArr){
		($errorParsed, $label, $start, $stop, $p_value, $pfam_id, $description, $source, $params, $version) = returnDomainsForProteinParsed($domainHash, $a);
		if (!$errorParsed){ 
		    print "DEBUG: Success! Found $a in hash for protien id $proteinId\n";
		    last;
		}
	      }	
	}
	
	if($errorParsed){
	    $COMPUTE = 1;
	    print "DEBUG: running computation on protein id $proteinId\n";
	    my $accArr = returnAccValuesByProteinId($dbh, $proteinId);
	    my $acc = $accArr->[0];
	    if($acc eq '-1'){
		handleError('handleDomainsForProteinId', "Protein Id $proteinId does not have accession", \@_);
		$errorFlag = 1;
	    }
	    else{

		($errorFlag, $start, $stop, $label, $p_value, $source, $params, $version) = returnDomainsByPrediction($dbh, $proteinId, $acc, $pvalueCutoff, $CPU);
	
	    }
	}

	
	if(!$errorFlag){
	    if(scalar(@$start)){
		print "DEBUG: Inserting domains!!\n";
		print "DEBUG: Source: $source\tParams: $params\tVersion:$version\n";
		$domainIds = insertDomainArr($dbh, $label, $start, $stop, $p_value, $source, $params, $proteinId, $version)
		}
	    else{
		my $sequence =  returnSeqByProteinId($dbh, $proteinId); 
		my $domainEmpty = insertDomain($dbh, '~~~', 0, length($sequence), 0, $source, $params, $proteinId, $version);
		push @$domainIds, $domainEmpty;

	    }
	    
	

	}
    

    #gross handling to take care of overlapping domains (when one domain completely falls in another one)

# add if COMPUTATIONAL then do this
	if($COMPUTE){
	    my ($hashOverlap, $domainsToRemove) = checkOverlappingDomains($dbh, $proteinId);
	if(scalar(keys %$hashOverlap)){ #need to remove domains
	    #delete domains
	    foreach my $domainId (@$domainsToRemove){
		deleteDomain($dbh, $domainId);
	    }
	    
	}
	    print "Domains found by COMPUTATION\n";
	}
	else{
	    print "domains found by PARSING\n";
	}
	
    }
    return ($errorFlag, $domainIds);
    
}



# $error = writeFastaFileForProteinId($proteinId, $acc, $fastaFile);
# Appends fasta file format to fasta file
# Inputs: $dbh - database handle
#         $proteinId - protein id (id for protein table)
#         $acc - accession 
#         $fastaFile - file to write fasta output to (appends)
# Kristen Naegle
# March 28, 2008
sub writeFastaFileForProteinId($$$$){
    my ($dbh, $proteinId, $acc, $fastaFile) = @_;
    my $error = 0;
    #find accession, name and sequence for proteinId
    my $sequence = returnSeqByProteinId($dbh, $proteinId);
    my $name = returnNameByProteinId($dbh, $proteinId);
    my $accArr = returnAccValuesByProteinId($dbh, $proteinId);
    if($sequence eq '-1' or $name eq '-1' or scalar(@$accArr) ==0){
	handleError('writeFastaFileForProteinId', 'Attempted to write fasta file for nonexisting protein or accession data', \@_);
	$error = 1;
    }
    else{
	my $header = createFastaHeader($accArr, $name);
    
	if(!-e $fastaFile){system("touch $fastaFile");}
	open(FH_98x, ">>$fastaFile") || die "Can't open Fasta File $fastaFile for appending\n";
	print FH_98x "$header\n";
	print "DEBUG: Fasta header $header\n";
	print FH_98x "$sequence\n\n";
	close(FH_98x);
    }
    return $error;
}

# $fastaHeader = createFastaHeader($accArr, $name);
# Returns a fasta header
# Inputs: $accArr - array of accessions for a protein id
#         $name - name for a protein id
# Outputs: $header - the fasta header to write at the top of a fasta file
# Kristen Naegle
# March 29, 2008
sub createFastaHeader($$){
    my ($accArr, $name) = @_;
    my $header = '>';
   foreach my $acc (@$accArr){
       $acc =~ s/^ //; # this was causing a problem with parsing into hmm if started with a space
       $header .= "$acc|";
    }
    $header .= "$name";
    return $header;
}

# ($errorCode, $MS_id, \@MS_phosphopepIds, \@phosphopepIds) = handlePhosphoSites($dbh, $pep, $expId, $proteinId)
# handles and inserts phosphorylation sites into phosphopep for a protein and a peptide from an MS line. If MS_id exists based on the peptide and experiment id, then returns MS_id and doesn't insert anything. 
# Inputs: $dbh - database handler
#         $pep - phosphopeptide from MS line (replaces pX values with x)
#         $expId - experiment table id this is associated with
#         $proteinId - protein table id this peptide is associated with
# Outputs: $errorCode - nonzero if problems with insertions 
#          \@MS_phosphopepIds - array of MS_phosphopep table ids that were generated on insertion
#          \@phosphopepIds - array of phosphopep table ids that were generated on insertion
# Kristen Naegle
# March 11, 2008
sub handlePhosphoSites($$$$){
    my ($dbh, $pep, $expId, $proteinId) = @_;
    my $errorCode = 0;
    my($MS_phosphopepId, $phosphopepId);
    my(@MS_phosphopepIds, @phosphopepIds);
    my $numAA = 7;
    $pep =~ s/pS/s/g;
    $pep =~ s/pT/t/g;
    $pep =~ s/pY/y/g;
    chomp $pep;
    
    my $MS_id = returnMSIdByPeptideExperiment($dbh, $pep, $proteinId, $expId);
#    print "MSID $MS_id\n";

    my ($errorCodeSingle, $phosphoPeps_arr) = returnSinglyPhosphoArray($pep);

    if($errorCodeSingle){
	$errorCode = 1; # NO phosphorylation site found;
	handleError('handlePhosphoSite', 'No phosphorylation site found', \@_);
    }
    elsif($MS_id == -1){
	my $sequence = returnSeqByProteinId($dbh, $proteinId);
	# get pfam array for protein
	my $count = 0;
	foreach my $singlePep (@$phosphoPeps_arr){
	    $count += 1;
	    my $tryps = returnTrypsPep($singlePep);
	    my ($errorCode_align, $pos, $type, $aligned) = returnAlignedandCodeForSinglyPhospho($sequence, $singlePep, $numAA);
	   
	    if(!$errorCode_align){
		if($count <= 1){ # make sure not to insert MS for multiply phospho forms
		    $MS_id = insertMSLine($dbh, $pep, $expId, $proteinId); #don't insert until we know everything will work out
		    if($MS_id == -1){
			$errorCode = 1;
			handleError('handlePhosphoSites', 'insertMSLine did not work', \@_);
			return($errorCode, $MS_id, \@MS_phosphopepIds, \@phosphopepIds);
		    }
		}
		my ($pepError, $pepId) = returnPhosphopepIdOnPepAlignedAndProteinId($dbh, $aligned ,$proteinId);
		if($pepId==-1){
		    my ($errorPfam, $pfam) = returnPfamSite($dbh, $proteinId, $pos);
		    if(!$errorPfam){
			
			($MS_phosphopepId, $phosphopepId) = insertPhosphopep($dbh, $tryps, $aligned, $pos, $type, $pfam, $proteinId, $MS_id);
			push @MS_phosphopepIds, $MS_phosphopepId;
			push @phosphopepIds, $phosphopepId;
			if($MS_phosphopepId == -1){ $errorCode = 1};
		    }
		    else{
			$errorCode = 1; # Error with pfam site
			handleError('handlePhosphoSite', 'Problem with getting pfam prediction where site falls', \@_);
		    }
		    
		}
		else{ #only need to insert a link between existing peptide and this MS id. 
		    $MS_phosphopepId = insertMSPhosphopep($dbh, $MS_id, $pepId);
		    if($MS_phosphopepId == -1){ 
			$errorCode = 1;
		    }
		    push @MS_phosphopepIds, $MS_phosphopepId;
		    push @phosphopepIds, $pepId;
		}
		
	    }
	    else{ #Can't find peptide in sequence
		$errorCode = 1; 
		handleError('handlePhosphoSite', 'Problem during alignment', \@_);
	    }

	}
	
    }

    elsif($MS_id != -1){

	my ($MS_phosphopepIds, $phosphopepIds) = returnMSPhosphopepTableForMSId($dbh, $MS_id);
	@MS_phosphopepIds = @$MS_phosphopepIds;
	@phosphopepIds = @$phosphopepIds;
    }
    
    return($errorCode, $MS_id, \@MS_phosphopepIds, \@phosphopepIds);
}


# \@data_ids = insertData($dbh, $dataHash, $line, $MS_id);
# Assumes you've already checked for existence of MS_id elsewhere, so will insert the data that exists for MS_id according to the current line and $dataHash loaded earlier from calling function/script. 
# Inputs: $dbh - database handle
#         $dataHash - hash of data headers
#         $line - the line being handled
#         $MS_id - (FK) id to MS entry
# Kristen Naegle
# March 17, 2008
sub handleData($$$$$){
    my($dbh, $dataHash, $line, $run, $MS_id) = @_;
    my @ids;
    my $error = 0;
    my @line = split("\t", $line);
    
    foreach my $col (keys %$dataHash){
	my ($type, $label, $priority) = returnDataHashValues($dataHash->{$col});
	#check to make sure we aren't putting in something with the same run number and the same priority..if so then set error.
	my $value = $line[$col]; #check 
	$value =~ s/ //g;
	if($value eq ''){
	    print "Found empty value\n";
	    $value = "NF";
	}
	my ($dataIdE) = returnDataIdByMost($dbh, $type, $run, $label, $MS_id);
	my $id;
	if($dataIdE == -1){
	    $id = insertData($dbh, $type, $run, $label, $priority, $value, $MS_id);
	}
	else{ 
	    $id = $dataIdE;
	}
	if($id != -1){
	    push @ids, $id;
	    
	}
	
	else{
	    $error = 1;
	}
	
	
    }
    
    return ($error, \@ids);
    
}

# \%hash = handleAllProteins($dbh, $accArr)
# Foreach accession in the array, adds the necessary accession and protein if it doesn't already exist
# Implements using batch query to entrez tools for speed
# Inputs: $dbh - database handle
#         $accArr - reference to array of accessions, please uniquify this before function call
# Outputs: $hashRef - reference to hash, with key $acc and array of values [$proteinId, $iProteinFlag, $iAccFlag]
# Kristen Naegle
# March 24, 2008
sub handleAllProteins($$){
    my ($dbh, $accArr) = @_;
    my (@proteinIds, @iProteinFlag, @iAccFlag);
    my %hash;
    my $date = returnMonthYear();
    #my @accArr = keys %$accHash;
    
    #remove isoform versions from accessions and remove accessions of an unknown type?
    my @newAccArr;
    foreach my $acc2 (@$accArr) {
	push @newAccArr, $acc2;
	
    }
    print "DEBUG: uniquifying array\n";
    $accArr = returnUniqueArray(\@newAccArr);

    print "DEBUG: Looking for missing accessions\n";
    my ($missAcc, $accIdHash) = returnMissingProteins($dbh, $accArr);
    my @genPeptAcc;
    my $groupSize = 100;
    my $start = 0;
    my $stop;
    foreach my $acc (@$missAcc){
	my ($errorSeq, $seq, $GENPEPT) = returnSeqFromAcc($acc);
#	if(!$errorSeq){
	if($GENPEPT){
	    push @genPeptAcc, $acc;
	    
	    print "DEBUG: $acc is saved for GenPept\n";
	}
	else{
	    
	    my ($proteinId, $iProteinFlag, $iAccFlag) = handleProteinByRichSeq($dbh, $acc, $seq, $date);
	    print "DEBUG: $acc is being inserted: proteinId: $proteinId\n";
	    $dbh->commit();
	    $hash{$acc} = [$proteinId, $iProteinFlag, $iAccFlag];
	}
    }
    
    #}
    
    print "Retreiving Genpept Stream".scalar(@genPeptAcc)." in groups of $groupSize\n";    
    while(scalar(@genPeptAcc)){
	if(scalar(@genPeptAcc) < $groupSize){
	    $stop = scalar(@genPeptAcc) - 1;
	}
	else{
	    $stop = $groupSize - 1;
	   
	}
	my @miss = @genPeptAcc[$start..$stop];
	splice(@genPeptAcc, $start, $stop-$start+1);
	my $richSeqHash = returnGenPeptQueryByBatch(\@miss);
	foreach my $acc (keys %$richSeqHash){
	    my ($proteinId, $iProteinFlag, $iAccFlag) = handleProteinByRichSeq($dbh, $acc, $richSeqHash->{$acc}, $date);
	    $hash{$acc} = [$proteinId, $iProteinFlag, $iAccFlag];
	}
	$dbh->commit();
    } # END WHILE MISSING GENPEPT ACCESSIONS
    
    foreach my $accInDB (keys %$accIdHash){
	$hash{$accInDB} = [$accIdHash->{$accInDB}, 0, 0];
	
    }
    
    return \%hash;
}
    

# ($error, $MS_id_Inserted) = loadLine($dbh, $line, $dataHash, $DATA, $expId, $accCol, $pepCol, $runCol);
# Loads a line from a data file, requires knowledge of header line through inputs. Loads (only if these don't already exist: phosphopep(s), protein, domain(s), accession(s), MS, data. No LONGER handles protein and domain. This is handled en masse in loadDataFile
# Inputs: $dbh - database handle
#         $line - the data line
#         $dataHash - hash of data found from header
#         $DATA - boolean flag about whether data exists in this file
#         $expId - (Key) id to experiment table
#         $accCol - column number of accession entry
#         $pepCol - column number of trypsinized phosphopeptide
#         $runCol - column number of run (if -1 will set to 'average')
# Outputs: $error - see error log for specific error, but this is true if any error occured
#          $MS_id_Inserted - 0 if MS_id was not inserted and 1 if a new one was inserted.
# Kristen Naegle
# March 18, 2008
# Modified October 16, 2008 - To pass back MS_id_Inserted flag - so we can catch redundanceis in large datasets.
sub loadLine($$$$$$$$){
    my ($dbh, $line, $dataHash, $DATA, $expId, $accCol, $pepCol, $runCol) = @_;
    my $error = 0;
    my $HMM_pvalueCutoff= 0.00001;

    my $pep = returnField($line, $pepCol);
    chomp $pep;
    $pep =~ s/pY/y/g;
    $pep =~ s/pS/s/g;
    $pep =~ s/pT/t/g;
    my $MS_id;
    my $acc = returnField($line, $accCol);
    my $accOld = $acc;
    $acc = removeAccIsoforms($accOld); 
    my $proteinId = returnProteinIdByAcc($dbh, $acc); # should have loaded proteins before this as bulk1
    my $MS_id_Inserted = 0;
    if($proteinId==-1){
	$error = 1;
	handleError('loadLine', "Protein does not exist in database, see earlier errors for acc=$acc with original form $accOld", \@_);
	$MS_id = -1;
    }
    else{     # not a bad protein
	
	my $MS_id_check = returnMSIdByPeptideExperiment($dbh, $pep, $proteinId, $expId);
	my ($errorPPep, $MS_phosphopepIds, $phosphopepIds);
	if($MS_id_check == -1){
	    ($errorPPep, $MS_id, $MS_phosphopepIds, $phosphopepIds) = handlePhosphoSites($dbh, $pep, $expId, $proteinId);
	    $error = $error | $errorPPep;
	    #   }
	    if($MS_id != -1){
		$MS_id_Inserted = 1;
		if($DATA){
		    my $run = returnRun($line, $runCol);
		    
		    my ($errorData, $dataIds) = handleData($dbh, $dataHash, $line, $run, $MS_id);
		    $error = $error|$errorData;
		    
		}
	    }
	    
	}
	else{ #If peptide MS exists, then just add the data
	    if($DATA){
		my $MS_id = $MS_id_check;
		my $run = returnRun($line, $runCol);
		my ($errorData, $dataIds) = handleData($dbh, $dataHash, $line, $run, $MS_id);
		$error = $error|$errorData;
		$MS_id_Inserted = 0;
	    }
	}
    }# end bad protein
    
	return ($error, $MS_id_Inserted);
}

# loadProteinAndDomainsInFile($dbh, $dataFile)
# For a data file.. load all the proteins and domains into dbh
# Inputs: $dbh - database handle
#         $dataFile - tab separated data file
# Outputs: $proteinHash - hash of proteins (keys $acc, values = [$proteinId, $iProteinFlag, $iAccFlag]
# Kristen Naegle
# March 31, 2008
sub loadProteinAndDomainsInFile($$){
    my ($dbh, $dataFile) = @_;

    #my ($dirSource, $fileName, $ext) = returnFileNames($dataFile);
;
    #my $CPU = 6;
    my $CPU = returnCPU();
    my $HMM_pvalueCutoff= 0.00001;

    # Load all proteins
    my $accCol = returnAccCol($dataFile);
    if($accCol == -1){
	handleError('loadDataFile', 'No Valid Accession Column', \@_);
	exit;
    }
    my $accHash = returnUniqColumnHash($dataFile, $accCol);
    my @acc = keys %$accHash;
    #my $accHash = returnAllHash($dataFile);
    print "INSERTING PROTEINS .....\n";
    printLineToLog("------------PROTEIN INSERTION ALL---------------------\n");
    my $proteinHash = handleAllProteins($dbh, \@acc);
    

    # For all proteins with no domain predictions, make domain predictions
    print "INSERTING DOMAINS.... Predicting ...\n";
    printLineToLog('--------------DOMAINS FOR PROTEIN ALL--------\n');
    maintainDomains($dbh);
    return $proteinHash;
}

# loadProteinAndDomainsInFile($dbh, $dataFile)
# For a data file.. load all the proteins and domains into dbh
# Inputs: $dbh - database handle
#         $dataFile - tab separated data file
# Outputs: $proteinHash - hash of proteins (keys $acc, values = [$proteinId, $iProteinFlag, $iAccFlag]
# Kristen Naegle
# March 31, 2008
sub loadProteinsInFile($$){
    my ($dbh, $dataFile) = @_;

    #my ($dirSource, $fileName, $ext) = returnFileNames($dataFile);
;
    #my $CPU = 6;
    my $CPU = returnCPU();
    my $HMM_pvalueCutoff= 0.00001;

    # Load all proteins
    my $accCol = returnAccCol($dataFile);
    if($accCol == -1){
	handleError('loadProteinsInFile', 'No Valid Accession Column', \@_);
	exit;
    }
    my $accHash = returnUniqColumnHash($dataFile, $accCol);
    my @acc = keys %$accHash;
    #my $accHash = returnAccHash($dataFile);
    print "INSERTING PROTEINS .....\n";
    printLineToLog("------------PROTEIN INSERTION ALL---------------------\n");
    my $proteinHash = handleAllProteins($dbh, \@acc);
    

  
    return $proteinHash;
}

# moving the script into this function, which calls loadDataFile, but also time stamps the log and places it in, as well as reporting errors.
# ($expId, $accepted, $rejected) = call_loadDataFile($dbh,$dataFile,$name, $author, $description, $contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_Link, $submissionEmail, $primaryMods);
# Load a data file into the database
# Inputs: $dbh - database handle (assumes opened with no commit)
#         $dataFile - tab separated data file
#         $name - name of experiment in database
#         $author - author of expt
#         $description - description of experiment
#         $contact - contact info for dataset
#         $PMID - pubmed id if published
#         $URL - URL if exists to dataset/paper
#         $published - boolean value - if data published this is set true.
#         $AMBIGUITY - boolean value - 1 if you want to load ambiguous relationships (ONLY FOR tryspinized fragemnts..do not use for datasets with aligned peptides), 0 otherwise.
#         $EXPORT - boolean value, 1 if allowed for export from PTMScout and 0 if not
#         $expId_Link - integer, if this is a different version of a dataset already loaded, than this is the experiment id for that dataset. 0 is default (a non-existant experiment id)
#         $submissionEmail - email of person submitting dataset
#         $primaryMods - list of primary modificaitons in string format
# Outputs: $expId - id to experiment table
#          $accepted - the number accepted
#          $rejected - the number rejected
# Kristen Naegle
# Nov 25, 2009
sub call_loadDataFile($$$$$$$$$$$$$$$$$$){
    my($dbh,$dataFile,$name, $author, $description,$contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_link, $submissionEmail, $primaryMods, $journal, $pub_date, $volume, $pages) = @_;
    
    #maintain domains in case there was a failure in the middle of some previous load leaving domains empty
    checkAndMaintainKilledProcess($dbh); 

    writeProcess('call_loadDataFile', \@_);
    
    flushLog();

    my ($expId, $accepted, $rejected);
    
    ($expId, $accepted, $rejected) = loadDataFile($dbh, $dataFile,$name, $author, $description, $contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_link, $submissionEmail, $primaryMods, $journal, $pub_date, $volume, $pages);
    

# Time stamp log.
    my $newLog = dateStampLogFile();
    print "LOG AT: $newLog\n";   

# enter log into experiment table
#first split to get just the last part
    my @log = split('/', $newLog);
    my $logFile = $log[$#log];
    enterLog($dbh, $expId, $logFile);
    
    $dbh->commit();

    

    writeProcessFinish();
    
    return($expId, $accepted, $rejected);
}

# Maintain ambiguous experiments (local only), expression and pelm Kinase annotations. Can call this after a new dataset load.
# Inputs: $dbh - database handle
# Kristen Naegle
# Dec. 8, 2009
sub maintainDBAfterFileLoad($){
    my($dbh) = @_;
    print "-----Maintaining Database after load-----\n";
    print "Maintaining local ambiguity links\n";
    my $LOCAL = 1;
    my $REFSEQ = 0;
    my $rdbh = returnRefSeqDBHNOCommit();
    maintainAmbiguousExperiments($dbh, $rdbh, $LOCAL, $REFSEQ);
    $dbh->commit();
    $rdbh->rollback();
    
    print "Maintaining expression linkages\n";
    my ($numHandled, $numFailed) = maintainExpressionLinkages($dbh);
    print "Handled: $numHandled\tFAILED: $numFailed\n";
    $dbh->commit();
	
    print "Maintaining phosphoELM Kinase annotations\n";
    my ($predIds) = maintainPELMKinaseAnnotations($dbh);
    
}


# ($expId, $accepted, $rejected) = loadDataFile($dbh,$dataFile,$name, $author, $description, $contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_Link, $submissionEmail, $primaryMods);
# Load a data file into the database
# Inputs: $dbh - database handle (assumes opened with no commit)
#         $dataFile - tab separated data file
#         $name - name of experiment in database
#         $author - author of expt
#         $description - description of experiment
#         $contact - contact info for dataset
#         $PMID - pubmed id if published
#         $URL - URL if exists to dataset/paper
#         $published - boolean value - if data published this is set true.
#         $AMBIGUITY - boolean value - 1 if you want to load ambiguous relationships (ONLY FOR tryspinized fragemnts..do not use for datasets with aligned peptides), 0 otherwise.
#         $EXPORT - boolean value, 1 if allowed for export from PTMScout and 0 if not
#         $expId_Link - integer, if this is a different version of a dataset already loaded, than this is the experiment id for that dataset. 0 is default (a non-existant experiment id)
#         $submissionEmail - email of person submitting dataset
#         $primaryMods - list of primary modificaitons in string format
# Outputs: $expId - id to experiment table
#          $accepted - the number accepted
#          $rejected - the number rejected
# Kristen Naegle
# March 18, 2008 
# Modified March 27, 2008 - to insert all proteins first
# Mod. Jan 9, 2009 to allow more fields
# Mod. July 12, 2009- to include new fields
sub loadDataFile($$$$$$$$$$$$$$$$$$){
    my($dbh,$dataFile,$name, $author, $description,$contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_Link, $submissionEmail, $primaryMods, $journal,$pub_date, $volume, $pages) = @_;
    my $error;
    my $MS_id_Inserted;

    my $date = returnTimeStr;

    my $domainHash = returnDomainHashFromGlobalFile();


    my $firstProteinId = getLastEntryId($dbh, 'protein');
    loadProteinsInFile($dbh, $dataFile);


    my $lastProteinId_1 = getLastEntryId($dbh, 'protein');
    my $insertedProteins_1 = getEntriesBetweenIds($dbh, 'protein', $firstProteinId, $lastProteinId_1); 
#    loadProteinAndDomainsInFile($dbh, $dataFile);
    print "INSERTING DOMAINS.... Predicting ...\n";
    printLineToLog("--------------DOMAINS FOR PROTEINS LOADED--------\n");
    handleDomainsForProteinIds($dbh, $insertedProteins_1, $domainHash);
    
    my @fileArrTemp = split('/', $dataFile);
    my $fileName = $fileArrTemp[$#fileArrTemp];

    my $expId = insertExperiment($dbh, $name, $author, $date, $description, $contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_Link, $fileName, $submissionEmail, $primaryMods, $journal, $pub_date, $volume, $pages);
    $dbh->commit(); # Make sure to commit this before an individual line can fail
    print "Experiment ID: $expId\n";

    
## Column numbers into datafile
    my $accCol = returnAccCol($dataFile);
    #my $pepCol = returnColumnNumber($dataFile, 'pep:tryps');
    my $pepCol = returnPeptideCol($dataFile);
    if($pepCol < 0){
	handleError('loadDataFile', 'Peptide column not found! Make sure to denote by pep in header', \@_);
	exit;
    }

    open(FH_IN, $dataFile) || die "Can't open $dataFile for reading in $0\n";
    my $line = <FH_IN>; # take out header
    my $dataHash = returnDataHash($line);
    my @dataCols = keys %$dataHash;
    my $DATA = 0;
    my $runCol;
    if(scalar(@dataCols) > 0){
	$DATA = 1;
    }
    if($DATA){
	$runCol = returnColumnNumber($dataFile, 'run');
    }
    my $lineCount = 2; # to mimic direct line number in file
    my $rejected = 0;
    my $accepted = 0;

    while(defined($line = <FH_IN>)){
	$error = 0;
	printBreakerLineToLog();
	printLineToLog("Data file: Line # $lineCount\n");
	    ($error, $MS_id_Inserted) = loadLine($dbh, $line, $dataHash, $DATA, $expId, $accCol, $pepCol, $runCol);
	if($error){
	    printLineToLog("REJECTED\n\n");
	    $dbh->rollback();
	    $rejected += 1;
	}
	else{
	    if($MS_id_Inserted){
		printLineToLog("ACCEPTED\n\n");
	    }
	    else{
		printLineToLog("ACCEPTED -- IGNORED\n\n");
	    }
	    $dbh->commit();
	    $accepted +=1;
	}
	$lineCount += 1;
    }
    close(FH_IN);
    print "-----Loading Scansite Predictions -----\n";
    loadScansitePredictions($dbh);
   $dbh->commit();


    if($AMBIGUITY){
	print "-----Loading Ambiguous relationships ----\n";
	my $rdbh = returnRefSeqDBHNOCommit();
	my ($LOCAL, $REFSEQ) = 1;
	my $FORCE = 0;
	insertAmbiguousPeptidesForExpId($dbh, $rdbh, $expId, $LOCAL, $REFSEQ, $FORCE);
	#  print "There were $count ambiguous peptides\n";
	$rdbh->rollback();
	$dbh->commit();
#	print "----- Predicting Domains for Proteins loaded for Ambiguity-----\n";
	#maintainDomains($dbh);
    }

    my $lastProteinId = getLastEntryId($dbh, 'protein');
    if($lastProteinId != $firstProteinId){
	my $insertedProteins = getEntriesBetweenIds($dbh, 'protein', $firstProteinId, $lastProteinId); 
	print "DEBUG: Inserting GO Terms from $firstProteinId to $lastProteinId\n";
    
    
	print "-----Loading GO Terms -----\n";
	printLineToLog("-------GO TERMS ----------\n");
#	my $ontologyFile = "/home/knaegle/SVN/knaegle/scripts/GO_DATA/gene_ontology.1_2.obo";
	handleGOTermsForProteinIds($dbh, $insertedProteins);
	$dbh->commit();
	# For all proteins with no domain predictions, make domain predictions
	print "INSERTING DOMAINS.... Predicting ...\n";
	printLineToLog("--------------DOMAINS FOR PROTEINS LOADED BY AMBIGUITY--------\n");
#    maintainDomains($dbh);
	my $proteinDomainsToPredict = getEntriesBetweenIds($dbh, 'protein', $lastProteinId_1, $lastProteinId);
	handleDomainsForProteinIds($dbh, $proteinDomainsToPredict, $domainHash);
    }
    $dbh->commit();
  
    

    return($expId, $accepted, $rejected);
    
}

# handleDomainsForProteinIds($dbh, $proteinIdArr, $domainHash)
# Given a list of protein ids, predict and insert domains
# Inputs: $dbh - database handle
#         $proteinIdArr - ref. to array of protein ids to make and insert domain predictions for
#         $domainHash - ref. to domain hash (parsed pfam file)
# Kristen Naegle
# May 31st, 2009
# Added parsing pfam file June 22, 2009
sub handleDomainsForProteinIds($$$){
    my ($dbh, $proteinIdArr, $domainHash) = @_;
    my $CPU = returnCPU();
    my $pvalueCutoff = 0.00001;
    foreach my $proteinId (@$proteinIdArr){
	my $accArr = returnAccValuesByProteinId($dbh, $proteinId);
	my $acc = $accArr->[0];
	if(!$acc){
	    print "ERROR: Accession for $proteinId protein is empty!!\n";
	    exit;
	}
	if($acc eq '-1'){
	    handleError('handleDomainsForProteinIds', "Protein Id $proteinId does not have accession", \@_);
	    next;
	}
	else{
	    my ($errorCode, $domainIds) = handleDomainsForProteinId($dbh, $proteinId, $pvalueCutoff, $CPU, $domainHash);
	    $dbh->commit();
	}
    }


}

# enterLog($expId, $logFile)
# After call_loadDataFile, time stamp the errorlog and then use this to enter into experiment table
# Inputs: $dbh - database handle
#         $expId - id of experiment table
#         $logFile - the name of the log file (better be less than 30 characters)
# Kristen Naegle
# May 11th, 2009
sub enterLog($$$){
    my ($dbh, $expId, $logFile) = @_;
    #check to see if there's a truncation
    if(length($logFile) > 30){
	print "ERROR: Your log file $logFile was truncated, limit of 30 characters\n";
    }
    my $sth = $dbh->prepare('update experiment set errorLog=? where id=?');
    $sth->execute($logFile, $expId);
    $sth->finish();
    
}

# loadScansitePredictions($dbh);
# Inserts scansite predictions for any phosphopeps that don't currently have predictions
# Inputs: $dbh- database handle
# Outputs: insertion into phosphopep_prediction
# Kristen Naegle
# April 4, 2008
# Modified Oct. 27, 2008 - To eliminate duplicate records. 
sub loadScansitePredictions($){
    my ($dbh) = @_;
    my $source = 'scansite';
    # get phosphopep id's missing scansite predictions 
    my $missingHash = returnPhosphopepIdHashWithNoPredictions($dbh, $source);
    print "FOUND MISSING SCANSITE PREDICTIONS\n";
    foreach my $id (keys %$missingHash){
	print "\t$id\t$missingHash->{$id}\n";
	my $pep = $missingHash->{$id};
	my $hashRef = returnScansiteHashForDB($pep);
	print "Scansite Predictions\n";
	foreach my $count (keys %$hashRef){
	    my ($source, $value, $score) = @{$hashRef->{$count}};
	    print "\t\t source->$source\tValue->$value\tScore->$score\n";
	    my $prediction_id = insertPhosphopepPrediction($dbh, $source, $value, $score, $id);
	    print "INSERTED with ID: $prediction_id\n";
	}
	# get scansite predictions for pep
    }
    


}



# \%hash = insertGOTerms($dbh, $proteinId, $associations, $aspect, $ontologyVersion, $annotationVersion)
# Inserts Go terms in array of associations for a particular proteinId and aspect type 
# Inputs: $dbh = database handle
#         $proteinId - protein.id 
#         $associations - ref to array of GO terms 
#         $aspect - 'P', 'C' or 'F'
# Outputs: \%hash - ref to hash, key is $association in array and value is array of flags [$iGO, $iProteinGO] -- these are 1 if they needed to be inserted
# Kristen Naegle
# April 15, 2008
sub insertGOTerms($$$$$$$){
    my ($dbh, $proteinId, $associations, $ontology, $aspect, $ontologyVersion, $annotationVersion) = @_;
    my %hash;
    foreach my $GO (@$associations){
	my $GO_id;
	my $protein_go_id;
	# check to see if GO exists in GO table
	my $iGO = 0;
	my $iProteinGO = 0;
	$GO_id = returnGOIdForGOTerm($dbh, $GO, $aspect);
	if($GO_id == - 1){
	    # insert into both tables
	    #my $node = $ontology->nodeFromId($GO);
	    my $term = $ontology->get_term($GO);
	    if($term){
		$GO_id = insertGO($dbh, $GO, $term->name, $aspect, $ontologyVersion);
		print "DEBUG: Inserting $term \t aspect:$aspect \t proteinId: $proteinId\n";
		$protein_go_id = insertProteinGO($dbh, $proteinId, $GO_id, $annotationVersion);
		$iGO = 1;
		$iProteinGO = 1;
	    #print "INSERTED GO and PROTEIN_GO\n";
	    }
	    else{
		handleError('insertGOTerms', 'Tried to insert GO terms for a term that is undefined. Term: $term, GO: $GO', \@_);
	    }
	}
	else{
	    $protein_go_id = returnProteinGOId($dbh, $proteinId, $GO_id);
	    if($protein_go_id == -1){
		
		# insert protein_go linkage 
		$protein_go_id = insertProteinGO($dbh, $proteinId, $GO_id, $annotationVersion);
		if($protein_go_id != -1){
		    $iProteinGO = 1;
		    #print "INSERTED Protein_GO only\n";
		}
	    }
	}
	my @result = ($iProteinGO, $iGO);
	push @{$hash{$GO}}, @result;
    }
 
    return \%hash;
}

# insertAmbiguousPeptidesForExpId($dbh, $refseqDbh, $expId, $LOCAL, $REFSEQ, $FORCE);
# Creates a temporary file for loading into ambigous peptides for alternate proteins given an experiment Id. 
# Need to finish .. # should be fixed seems we won't insert more than the same pep_tryps (also need to load into test database (check accession and peptryps existence)..this assumes production 
# Bases ambiguity on tryps pep hit in refseq that DOES NOT have the exact same sequence
# Inputs: $dbh - database handle
#         $rdbh - refseq database handle
#         $expId - the experiment Id for which to handle
#         $LOCAL - 1 if you want to add ambiguity from local database and 0 if not
#         $REFSEQ - 1 if you want to add ambiguity from refseq database overalp and 0 if not
#         $FORCE - if 0 than don't force renewal of ambiguity check (i.e. if ambiguity entries exist for a peptide fragment, return without updating 
# Kristen Naegle
# May 1, 2008
# Modified August 6, 2008 - to accomadate new table
# Mod June 24, 2009 - to include ambiguity for proteins already in database
sub insertAmbiguousPeptidesForExpId($$$$$$){
    my ($dbh, $rdbh, $expId, $LOCAL, $REFSEQ, $FORCE) = @_;
    my $count = 0;
    print "DEBUG: IN insertAmbiguous FORCE: $FORCE\n";
    my $hashRef = returnMSInfoByExpId($dbh, $expId);
    for my $i (keys %$hashRef){
	my $tempHash = $hashRef->{$i};
	my $id = $tempHash->{'MS.id'};
	my $gi = $tempHash->{'acc.value'};
	my $peptide = uc($tempHash->{'MS.phosphopep'});
	my $proteinId = $tempHash->{'protein.id'};
	print "Handling for $gi\t$peptide\t$proteinId\n";
	handleAmbiguityForPeptide($dbh, $rdbh, $peptide, $proteinId, $LOCAL, $REFSEQ, $FORCE);
    }
    
}


# handleAmbiguityForPeptide($dbh, $rdbh, $peptide, $proteinId);
# Handles ambiguity for a peptide fragment, inserts into protein and ambiguity tables as needed
# Inputs: $dbh - database handle
#         $rdbh- refseq database handle
#         $peptide - the peptide fragment
#         $proteinId - the protein that fragment is currently assigned to
#         $LOCAL - 1 if you want to add ambiguity from local database and 0 if not
#         $REFSEQ - 1 if you want to add ambiguity from refseq database overalp and 0 if not
#         $FORCE - if 0 than don't force renewal of ambiguity check (i.e. if ambiguity entries exist for a peptide fragment, return without updating 
# Kristen Naegle
# August 8, 2008
# Mod June 24, 2009 - to include a search in current database handle
# Mod July 1, 2009 - to include LOCAL, REFSEQ and FORCE so that I can maintain ambiguity and do it in a quick way
sub handleAmbiguityForPeptide($$$$$$$){
    my($dbh, $rdbh, $peptide, $proteinId, $LOCAL, $REFSEQ, $FORCE) = @_;
    
    print "DEBUG: Force: $FORCE\n";

    my $ambId = returnAmbIdForFragProtein($dbh, $peptide, $proteinId);
    if(!$FORCE){
	if($ambId != -1){
	    return;
	}	
    }

    my $proteinHash = returnProteinDescForProteinId($dbh, $proteinId);
    my $species = $proteinHash->{'species'};
    my $date = returnMonthYear(); 
 
    # Look up other proteins that match 
    if($REFSEQ){
	my $matchHash = returnGIOnFragmentMatch($rdbh, $peptide, $species);
	my @accArr = keys %$matchHash;
	my $hashRichSeq = returnGenPeptQueryByBatch(\@accArr);
	foreach my $gi (keys %$matchHash){
	    
	    my $rHash = $matchHash->{$gi};
	    my $richSeq = $hashRichSeq->{$gi};
	    
	    print "Matches: $gi\t$rHash->{'name'}\t$rHash->{'species'}\n"; #multiple gi entries for the same protein so need to handle that 
	    my ($proteinId_GI, $iProteinFlag, $iAccFlag) = handleProteinByRichSeq($dbh, $gi, $richSeq, $date); 
	    print "Insertion:  proteinId: $proteinId_GI proteinFlag:$iProteinFlag acc:$iAccFlag\n";
	    my ($ambId, $iAmb) = handleAmbInsertion($dbh, $peptide, $proteinId_GI);
	    print "Amb Id: $ambId and Insert: $iAmb\n";
	    
	}
	my ($ambIdO, $iAmbO) = handleAmbInsertion($dbh, $peptide, $proteinId);  #adding insertion of the original
	print "Amb Id: $ambIdO and Insert: $iAmbO\n";
	# insert all proteins and fragments (including current) - includes at least one accession for new proteins 
    }

    if($LOCAL){
	#now do it from the current database:
	my $ambInsertedFromDB = handleAmbiguityFromCurrentDatabase($dbh, $peptide, $proteinId, $species);
    }
}

# \@insertedAmbIds = handleAmbiguityFromCurrentDatabase($dbh, $peptide, $proteinId, $species);
# Given a peptide and the database handle to our development databases, insert other ambiguity matches.
# Inputs: $dbh - database handle
#         $peptide - the peptide fragment
#         $proteinId - protein.id
#         $species - binomial name of species of proteinId
# Outputs: \@ambIds - ref to array of ambiguity ids that were inserted
# Kristen Naegle
# June 24, 2009
sub handleAmbiguityFromCurrentDatabase($$$$){
    my ($dbh, $peptide, $proteinId, $species) = @_;
    my @ambIds;
    my $proteinIdsHit = findFragmentsFromProteins($dbh, $peptide, $species);
    foreach my $proteinHit (@$proteinIdsHit){
	my ($ambId, $iAmb) = handleAmbInsertion($dbh, $peptide, $proteinHit);
	push @ambIds, $ambId;
    }
    return \@ambIds;
}

# ($ambId, $iAmb) = handleAmbInsertion($dbh, $peptide, $proteinId);
# inserts ambiguity entry unless it already exists
# Inputs: $dbh - database handle
#         $peptide - the peptide fragment
#         $proteinId - protein.id
# Outputs: $ambId - id to ambiguity table of entry
#          $iAmb - boolean flag, when 1 means ambiguity table was entered
# Kristen Naegle
# August 8, 2008
sub handleAmbInsertion($$$){
    my($dbh, $peptide, $proteinId) = @_;
    my $iAmb = 0;
    $peptide = uc($peptide);
    my $ambId = returnAmbIdForFragProtein($dbh, $peptide, $proteinId);
    if($ambId == -1){
	$ambId = insertAmbiguity($dbh, $peptide, $proteinId);
	$iAmb = 1;
    }
    return($ambId, $iAmb);
}


# Insert PELM kinase predictions into phosphopep_prediction from pELM file
# Inputs: $dbh - database handle
#         $pelm_file - current pelm file of annotations
#         $release_date - year for release date (e.g. 2006, 2008)
# Kristen Naegle
# October 13, 2008
sub handlePELMKinaseAnnotations($$){
    my ($dbh, $pelm_file) = @_;

    my $source = "pelm_kinase";

    my @predIds; 
    my $count = 0;
    #Check for headers and get column numbers that match necessary headers
    my %colHash;
    my @cols = ('acc:sp', 'pmids', 'kin:pelm', 'pep:tryps');
    foreach my $col (@cols){
	$colHash{$col} = returnColumnNumber($pelm_file, $col);
	if($colHash{$col} < 0){
	    print "ERROR: Can't find column $col in $pelm_file\n";
	    exit;
	}
    }
    open(PELM, $pelm_file) || die "Can't open phosphoELM file for reading\n";

    my $line = <PELM>; #cut header
    while(defined($line=<PELM>)){
	chomp $line;
	my @line = split("\t",$line); 
	my $kinases = $line[$colHash{'kin:pelm'}];

	if($kinases){
	    $count += 1;
	    # look up phosphopep.id based on sp acc and pep_tryps
	    my $acc = $line[$colHash{'acc:sp'}];
	    
	    my $pepTryps = $line[$colHash{'pep:tryps'}];
	    my $phosphopep_id = returnPhosphopepIdOnAccPepAligned($dbh, $acc, $pepTryps);
	    if($phosphopep_id != -1){
		my $pmids = $line[$colHash{'pmids'}];
		my @pmids = split(';', $pmids);
		
		my @kinases = split(';', $kinases);
		for(my $i=0; $i<scalar(@kinases); $i++){
		    my $k = $kinases[$i];
		    my $pmid = $pmids[$i];
		    # check for existence of record based on phosphopep_id, source, and value # if doesn't exist, insert into table
		    my $predId = returnPhosphopPepPredId($dbh, $source, $k, $phosphopep_id);
		    if($predId == -1){
			$predId = insertPhosphopepPrediction($dbh, $source, $k, $pmid, $phosphopep_id);
			if($predId == -1){
			    handleError('handlePELMKinaseAnnotations', "Could not insert pep prediction for $source, $k, $pmid, and $phosphopep_id", \@_);
			}
			else{
			    print "FOUND Kinase: @kinases and phosphopep_id: $phosphopep_id\n"; 
			    push @predIds, $predId;
			}
		    }
		    
		}
	    }
	    else{
		handleError('handlePELMKinaseAnnotations', "Could not find a phosphopep_id fitting description for $acc and $pepTryps", \@_);
	    }
	}
    }
    close(PELM);
    print "FOUND $count PELM Kinase annotations\n";
    return \@predIds;

} 

# $string = outputDataLoadStatementForExpId($dbh, $expId);
# Return the string for calling the data load based on info in the experiment table
# Inputs: $dbh -database handle
#         $expId - experiment id
# Outputs: $string - string including perl call_loadDataFile followed by appropriate arguments.
# Kristen Naegle
# Sept. 15, 2009
sub outputDataLoadStatementForExpId($$){
    my ($dbh, $expId) = @_;
    my $sth = $dbh->prepare('SELECT * from experiment where id=?');
    $sth->execute($expId);
    my $results = returnMultipleResultsForExSTH($sth);
    my $hash = $results->[0];
    my $dataDir = $globalVars::DATA_PATH;
    my $string = "perl call_loadDataFile.pl $dataDir$hash->{'dataset'} '$hash->{'name'}' '$hash->{'author'}' '$hash->{'description'}' '$hash->{'contact'}' $hash->{'PMID'} '$hash->{'URL'}' $hash->{'published'} $hash->{'ambiguity'} $hash->{'export'} $hash->{'experiment_id'} '$hash->{'submitter'}' '$hash->{'primaryModification'}' '$hash->{'journal'}' '$hash->{'pub_date'}' '$hash->{'volume'}' '$hash->{'pages'}'";
    return $string;
}

1;
