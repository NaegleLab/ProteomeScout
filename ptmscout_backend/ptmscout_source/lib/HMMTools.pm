use strict;
use warnings;
use Bio::Tools::Hmmpfam;
use SBMLTools;
use entrezTools;
use commonTools;
use DBTools::dbEntry;
use globalVars;


# $HMMHashRef = returnHMMFileAccHash($hmmInput)
# Parse an Hmmpfam for unique accessions
# Inputs: $hmmFile - hmm file output of pfam (example: hmmpfam Pfam_fs fastaInput > pfamOutput
# Outputs: A hash where keys are the fasta names and an array of parsed HMM features
# Kristen Naegle
# September 14, 2007
sub returnHMMFileAccHash($){
    my $hmmFile = shift;
    my %hmmAccHash; 
    my $fh;
    if(-e $hmmFile){
	open($fh, $hmmFile) || die "Can't open HMM File $hmmFile\n";
	my $hmmpfam_parser = Bio::Tools::Hmmpfam->new(-fh => $fh);
	my @hmmpfam_feat;
	while(my $hmmpfam_feat = $hmmpfam_parser->next_result ){
	    push @hmmpfam_feat, $hmmpfam_feat;
	    #print "Found feature\n";
	}
	foreach my $feature (@hmmpfam_feat){
	    my $acc = $feature->seq_id;
	    if(not $hmmAccHash{$acc}){
		$hmmAccHash{$acc} = [];
		
	    }
	    push @{$hmmAccHash{$acc}}, $feature;
	    #print "DEBUG: pushing onto $acc ". $feature->primary_tag."\n";
	    
	}
	close($fh);
	return \%hmmAccHash;
    }
    else {

	return \%hmmAccHash; # this should be empty
    }
}

# $prunedHmmHashRef = pruneHmmHash($hmmHashRef, $pValCutoff);
# take a hash for HMM features - as in returnHMMFileAccHash($file), return a new hash with only those features that have a p-value below cutoff
# Inputs: $hmmHashRef - reference to an HMM hash
#         $pValCutoff - minimum p-value cutoff for HMM hit
# Outputs: $prunedHashRef - reference to a new hash that has only p-valueCutoff or better hits 
# Kristen Naegle
# September 15, 2007
sub pruneHMMHash($$){
    my $hmmHashRef = shift;
    my $pValCutoff = shift;

    my %prunedHash; 
    
    # key is protein 
    foreach my $key (keys %$hmmHashRef){
	if(not $prunedHash{$key}){
	    $prunedHash{$key} = [];
	}
	foreach my $feature (@{$hmmHashRef->{$key}}){
	    if($feature->score <= $pValCutoff){
		push @{$prunedHash{$key}}, $feature;
#		print "DEBUG: keeping ".$feature->primary_tag."\n";
# # wasn't missing values at prune
	    }

	}
	

    }
    return \%prunedHash;
}

# $hmmHashRec = reconcileHMMDomains($hmmHashRef);
# Now go through each key and determine if there are overlapping domain predictions.. Choose the domain with the better score
# THIS IS extremely unelegant. but it works. 
# Inputs: $hmmHashRef - reference to an HMM hash ref
# Outputs: $hmmHashRec - reconciled hmmHash where the lower scoring domain is removed when there's an overlap
# Kristen Naegle
# September 27, 2007
sub reconcileHMMDomains($){
    my $hmmHashRef = shift;
    my %hmmHashRec; 
    my $i;
    # Features are in order
    my $fprint;
    foreach my $key (keys %$hmmHashRef){
	$hmmHashRec{$key} = [];
#POOR ASSUMPTION - Need to order based on start positions for this check
	my %startHash;
	my $start; 
	print "DEBUG: Number of original domains: ".scalar(@{$hmmHashRef->{$key}})."\n";
	for($i=0; $i < scalar(@{$hmmHashRef->{$key}}); $i++){
	    $start = ($hmmHashRef->{$key}->[$i])->start;
	    if(not defined($startHash{$start})){
		$startHash{$start} = $i;
		#print "Putting $i into $start\n";
	    }
	    else{ #in the case of identical start positions just get rid of the lower p-value domain 
		my $checkIndex = $startHash{$start};
		my $newScore = ($hmmHashRef->{$key}->[$i])->score;
		my $oldScore = ($hmmHashRef->{$key}->[$checkIndex])->score;
		#print "New Score: $newScore\n";
		#print "Old Score: $oldScore\n";
		if($newScore < $oldScore){
		    
		    $startHash{$start} = $i; #overwrite only if this is a better score)
		}
	    }
	    
	}


	# Now sort keys, which are the start of sequence
	my @keysToIgnore;
	my @keysToKeep;
	my @OVERLAP;
	my @NO_OVERLAP;
	my @sortedStart = sort {$a <=> $b } keys %startHash;
	print "DEBUG: number of domains under consideration ".scalar(@sortedStart)."\n";
	for ($i= 0; $i < scalar(@sortedStart); $i++){
	    my $index = $startHash{$sortedStart[$i]}; 
	    my $indexLast;
	    my $indexNext;
	    my $OVERLAP;
	    if($i > 0){
		$indexLast = $startHash{$sortedStart[$i-1]};
		my $start = ($hmmHashRef->{$key}->[$index])->start;
		my $endLast = ($hmmHashRef->{$key}->[$indexLast])->end;
		$OVERLAP =  ($start <= $endLast); 
	    } 
	    else { $OVERLAP = 0;}
	    if($i < $#sortedStart){
		$indexNext = $startHash{$sortedStart[$i+1]};
		my $end = ($hmmHashRef->{$key}->[$index])->end;
		my $startNext = ($hmmHashRef->{$key}->[$indexNext])->start;
		my $OVERLAP_NEXT =  ($end >= $startNext); 
		$OVERLAP = $OVERLAP || $OVERLAP_NEXT;
		if ($OVERLAP) {}#print "IT's TRUE!";}
		} 
	    else { $OVERLAP = $OVERLAP || 0; }
	    if($OVERLAP){
		push @OVERLAP, $index;
  	    }
	    else {
		push @NO_OVERLAP, $index;
		print "DEBUG: Index: $index has no overlap\n";

	    }
	}
	@keysToKeep = @NO_OVERLAP;
	# # $val is reference to array of indexes that overlapped. Although 
	# #since we pushed indexes on in order of start...i can do the same thing now
# #	print "DEBUG: before runnign through overlap, number to keep ".scalar(@keysToKeep)."\n";
	# #print "DEBUG: ".scalar(@OVERLAP)." number of indexes that overlap\n";
	for($i=0; $i < $#OVERLAP; $i++){
	    my $index = $OVERLAP[$i];
	    my $indexNext = $OVERLAP[$i+1];
	    my $endLast = ($hmmHashRef->{$key}->[$index])->end;
	    my $start = ($hmmHashRef->{$key}->[$indexNext])->start;
	    if($start < $endLast){
		if(($hmmHashRef->{$key}->[$index])->score < ($hmmHashRef->{$key}->[$indexNext])->score){
		    push @keysToKeep, $index;
		    push @keysToIgnore, $indexNext;

		}
		else{
		    push @keysToKeep, $indexNext;
		    push @keysToIgnore, $index;
		}
	    }
	    else{ #if they just happen to abutt by one amino acid, keep both
		push @keysToKeep, $index;
		push @keysToKeep, $indexNext;

	    }
	
		
	}

	# #print "DEBUG: Ignoring ".scalar(@keysToIgnore)."\n";
# #	print "DEBUG: Keeping ".scalar(@keysToKeep)."\n";

	# #Have to remove those things that were pushed into the no overlap bin before tested 
	# #make sure there are no keysToIgnore in the keys to keep; 

	my %seen = ();
	my @aOnly = ();
	foreach my $item (@keysToIgnore) { $seen{$item} = 1;}
	foreach my $item (@keysToKeep) {
	    unless ($seen{$item}){
		push @aOnly, $item;
	    }
	}

	# #print "Keeping @aOnly \n";
	# #print "Removing @keysToIgnore\n";
	my $hashKeysToKeep = createUniqHash(\@aOnly);
	foreach my $featToKeep (keys %$hashKeysToKeep){
	    my $feature = $hmmHashRef->{$key}->[$featToKeep];
	    push @{$hmmHashRec{$key}}, $feature;
	  #  print "DEBUG: REconcile: keeping ".$feature->primary_tag."\n";

	}

    }
    return \%hmmHashRec;

}



# createHMMResults($fastaFile, $hmmResultsFile, $HMMFile);
# runs hmmpfam on fasta sequence file. Checks existing pfam result file for predictions generated previously (wont' overwrite). Appends these new predictions to existing file (Currently requires gi is in accession line)
# Inputs: $fastaFile - a fasta formatted file with all sequences you wish to run HMM domain predictions on 
#         $hmmResultsFile - the place where the HMM results will be stored. Alsoa running compilation of everything tested 
#         $HMMFile - the domain labels for the learning portion of HMM (e.g. Pfam_fs
# Outputs: A written output file $hmmResultsFile that can then be parsed 
# Kristen Naegle
# September 17, 2007 
sub createHMMResults{
    my($fastaFile, $hmmResultsFile, $HMMFile, $CPU);
    if(scalar(@_) == 3){
	($fastaFile, $hmmResultsFile, $HMMFile) = @_;
	
	#$CPU = 6;
	$CPU = $globalVars::CPU_DAY;
    }
    elsif(scalar(@_) == 4){
	($fastaFile, $hmmResultsFile, $HMMFile, $CPU) = @_;

    }
    
    

    # start by parsing results into hash 
    my $hmmResultsHashRef = returnHMMFileAccHash($hmmResultsFile);
    my $fastaHash = returnFastaHash($fastaFile);
    #Go through fastaFile and accept fasta sequences for creation of temporary file on which hmmpfam will be run 
    #open(FASTA, $fastaFile) || die "Can't open your Fasta sequence file $fastaFile in createHMMResults\n";

    #my $line; 
    my $tempFile = "HMMFASTA_TEMPFILE";
    my $tempHMM = "HMMRESULT_TEMPFILE";
    if(!-e $tempFile) {system("touch $tempFile");}
    if(!-e $tempFile) {system("touch $tempHMM");}
    open(TEMP, ">$tempFile") || die "Can't open temp file for writing in createHMMResults\n";
    foreach my $key (keys %$fastaHash){
# #while(defined($line = <FASTA>)){
# #	if($line =~ /gi/){
# #   my $key = parseFastaAccLine($line);
	    if($hmmResultsHashRef->{$key}){

		print "$key already looked up in pfam\n";
	    }
	    else{
		print TEMP ">".$key."\n"; ##CHECK THAT this is printing the entire sequence -- append these new fasta sequences to 
		print TEMP $fastaHash->{$key}."\n\n";

	    }
	
	    

	}
#    close(FASTA);
    close(TEMP);
    
# # NOW RUN pfam 
# #    print("Running pfam on $tempFile");
# #    system("hmmpfam $HMMFile $tempFile > $tempHMM");
    system("hmmpfam --cpu $CPU $HMMFile $tempFile > $tempHMM");
    
    #cat the two files together
    my $temp3 = "TEMP_HMM_DONE";
    if(!-e $temp3) {system("touch $temp3");}
    system("cat $hmmResultsFile $tempHMM > $temp3");
    system("cp $temp3 $hmmResultsFile");


    #rm temp
    
    system("rm $tempHMM");
    system("rm $temp3");
#    system("rm $tempFile");

}


# $key = parseFastaAccLine($line)
# Given the first line of a fasta accesion, return the accession key
# Inputs: $line - the line from fasta representing the accessions
# Outputs: $key - the accession key 
# Kristen Naegle
# September 17, 2007
sub parseFastaAccLine($){
    my $line = shift;
    my @line = split(' ', $line);
    my $key = $line[0];
    $key =~ s/>//;
    return $key;
    
}

# $featureLine = printFeatureInfo($feature)
# Returns a tab separated line of feature info 
# Inputs: $feature - feature object of HMM
# Outputs: $featureLine - primary_tag\t start \t end \t score
# Kristen Naegle
# September 17, 2007
sub printFeatureInfo($){
    my $feature= shift;
    my $line; 
    $line = $feature->primary_tag; 
    $line .= "\t".$feature->start;
    $line .= "\t".$feature->end;
    $line .= "\t".$feature->score;
    return $line;
}

# $text = printProtDomainsForFile($hmmHashRef, $accKey, $site)
# Takes an hmmHash and the key of interest and condenses the information into a column type form. Returns two columns, the first is domain info for the whole protein and the second is domain info if the site falls into a domain.
# Inputs: $hmmHashRef - ref to hmmHash parsed from an hmm results file (can use prune and reconcile functions to clean up original parse)
#         $accKey - the accesion of interest
#         $site - the site to check for domain inclusion
# Outputs: two tab separated columns of format below
# format: domain:start:end:score; domain:start:end:score; etc.
# Kristen Naegle
# September 18, 2007
sub printProtDomainsForFile($$$){
    my $hmmHash = shift;
    my $accKey = shift;
    my $site = shift;
    $site =~ s/([A-Z])//gi;
    my $dprint; 
    my $temp; 
    my $sitePrint;
    foreach my $feature (@{$hmmHash->{$accKey}}){
	$temp = $feature->primary_tag.":".$feature->start.":".$feature->end.":".$feature->score.";";
	$dprint .= $temp;
	if(($site >= $feature->start) and ($site <= $feature->end)){
	    #print "FOUND HIT";
	    if(not $sitePrint){
		$sitePrint = $temp;
	    }
	    #need to see if this new domain supercedes the old based on score
	    else{
		my @oldScore = split(":", $sitePrint);
		my $oldScore = $oldScore[$#oldScore];
		$oldScore =~ s/;//g; #strip semicolon
		print "Comparing $oldScore with ".$feature->score."\n";
		#replace if it has a better p-value
		if($oldScore > $feature->score){
		    $sitePrint = $temp;
		}

	    }
	}
    }
    if(not $dprint){
	$dprint = "~~~";
    }
    if(not $sitePrint){
	$sitePrint = "~~~";
    }
    return $dprint."\t".$sitePrint;
}

# $errorFlag = appendDomainInfo($dataFile, $domainFile, $fastaFile, $outputFile, $pvalue)
# Main function to call to run through text file and append domain predictions
###STEPS:
# 1. Load unique accessions from the input file 
# 2. Check for those accessions in the Domain file
# 3. Nonexisting accessions..check for them in teh FASTA file
# 4. Get FASTA sequences for those missing and append to $fastaFile
# 5. Now run HMM predictions on the fasta file 
# 6. Get new hash and append info
# Inputs: $dataFile - txt file of data
#         $domainFile - HMM file, this will be modified if new predictions are made
#         $fastaFile - file of fasta sequences, modified if new fasta sequences are needed for predictions
#         $outputFile - destination for datafile + appended columns of pfam predictions
#         $pvalue - minimum pvalue allowed for predictions
# Outputs: $errorFlag - returns a 1 if there are proteins in dataFile for which there are noHMM predictions due to an error in either fasta capture or domain prediction
# Kristen Naegle
# November 1, 2007 
sub appendDomainInfo($$$$$){
    my $dataFile = shift;
    my $domainFile = shift;
    my $fastaFile = shift;
    my $outputFile = shift;
    my $pvalue = shift;

    my $errorFlag = 0;



    # CHECK $dataFile accessions and HMM hash against eachother
    my $accCol = returnAccCol($dataFile);
    my $accHash = returnUniqColumnHash($dataFile, $accCol);
    my $hmmHash = returnHMMFileAccHash($domainFile);
    my $missAcc = returnAccMissingFromHash($accHash, $hmmHash);
    my $tempFastaFile = "TEMP_FASTA_FILE";

    if (scalar(@$missAcc) > 1){
	my $failed = printRequestedFastaToFile($accHash, $missAcc, $tempFastaFile, $fastaFile);
	if(scalar(@$failed) > 0){
	    print "ERROR: FASTA sequences for these accessions could not be retrieved: \n";
	    foreach my $f (@$failed){
		print "$f\n";
	    }
	}    
	# Run HMM predictions (eventually move pfam training predictions to environment variable
	my $HMMTraining = "/data/knaegle/data/pfam/Pfam_fs";
	print "Creating domain predictions for new fasta sequences \n";
	createHMMResults($tempFastaFile, $domainFile, $HMMTraining);
	
	# Recheck to make sure everything was covered:
	$hmmHash = returnHMMFileAccHash($domainFile);
	$missAcc = returnAccMissingFromHash($accHash, $hmmHash);
	if(scalar(@$missAcc) > 1){
	    print "ERROR:  Still missing accessions\n";
	    foreach my $miss (@$missAcc){
		
		print "$miss\n";
	    }
	    $errorFlag = 1;
	    return $errorFlag;
	}   
	
    }
    
	#append HMM info to file;
    # new translatehash
    my $hmmHashPruned = pruneHMMHash($hmmHash, $pvalue);
    my $hmmHashRec = reconcileHMMDomains($hmmHashPruned);
    my $giTranslateHash = returnTranslateAccHash($accHash, $hmmHash);
    open(FHX, $dataFile) || die "Can't open input data file $dataFile for reading\n";
    if(!-e $outputFile){system("touch $outputFile");}
    open(OUT, ">$outputFile") || die "Can't open output file $outputFile for writing\n";
    my $siteCol = returnColumnNumber($dataFile, "site:aligned");
    my $giCol = returnAccCol($dataFile);
#header line
    my $line = <FHX>; 
    chomp $line;
    $line .= "\tdomains:pfam\tsiteDomain:pfam\n";
    print OUT $line;
    while(defined($line = <FHX>)){
	chomp $line;
	my @line = split("\t", $line);
	my $gi = returnAcc($line, $giCol); #$line[$giCol];
	my $giKey = $giTranslateHash->{$gi};
	my $site = (split(";", $line[$siteCol]))[0];
	my $append = printProtDomainsForFile($hmmHashRec, $giKey, $site);
	if(not $append){
	    $append = " ";
	}
	print OUT $line."\t".$append."\n"; 
    }
    close(FHX);
    return $errorFlag;
}

# $missingArrRef = returnAccMissingFromHash($accHash, $hmmHash);
# Given a GI hash and an HMM hash, compare the two to see what accessions from the acc hash are missing in the HMM predictions (works for fastaHash as well)
# Inputs: $accHash - a hash of accessions for which you want domain predictions
#         $hmmHash - hash of accession for which there are domain predictions already
# Outputs: $missingArrRef - reference to array of missing accessions
# Kristen Naegle
# November 1, 2007
sub returnAccMissingFromHash($$){
    my $accHash = shift;
    my $hmmHash = shift;
    my @missingAcc;
    my $giTranslateHash = returnTranslateAccHash($accHash, $hmmHash);
    foreach my $key (keys %$giTranslateHash){
	#print "$key -> $giTranslateHash{$key}\n";
	if($giTranslateHash->{$key} eq 0){
	    push @missingAcc, $key;
# #	print "ERROR: GI number $key does not have a pfam prediction in $domainFile\n";
	}
    }
    
    return \@missingAcc;
    
}

# $translateHashRef = reutrnTranslateAccHash($$);
# Need to find the full key to HMM and Fasta hashes that describe a smaller subset of an accession
# Inputs: $accHash - reference to accession hash (min descriptor)
#         $hmmHash - hmm or fasta hash of a full key word 
# Outputs: $translateHashRef - reference to the translation, where the key is the key in $accHash and the value is the full descriptor in hmmHash
# Kristen Naegle
# November 1, 2007
sub returnTranslateAccHash($$){
    my $accHash = shift;
    my $hmmHash = shift;
    my %giTranslateHash;
    foreach my $gi (keys %$accHash){
	#my $giF = returnAccNumber($gi);
	if(not $giTranslateHash{$gi}){
	    $giTranslateHash{$gi} = 0;
	    foreach my $key (keys %$hmmHash){
		if($key =~ /\Q$gi/){
		    $giTranslateHash{$gi} = $key;
		    last;
		}
		
	    }
	}
	
    }
    return \%giTranslateHash;

}

# $fastaHash = returnFastaHash($fastaFile);
# Opens a fastaFile and stuffs all accession lines into hash, value of hash is number of times the same accession shows up (should only be once in the ideal situation)
# Inputs: $fastaFile - file of fasta sequences with accessions on a line (denoted by >) and second line is sequence
# Outputs: $fastaHash - hash of fasta sequence accessions - value is sequence
# Kristen Naegle
# November 1, 2007
# Modified December 12, 2007 - to return the sequence
sub returnFastaHash($){
    my $fastaFile = shift;
    open(FASTA, $fastaFile) || die "Can't open Fasta file $fastaFile in returnFastaHash\n";
    my $line;
    my %fastaHash;
    while(defined($line = <FASTA>)) {
	chomp $line;
	if($line =~ />/){
	    my $acc = parseFastaAccLine($line);

	    if(not $fastaHash{$acc}){
		$fastaHash{$acc} = 0;
	      
	    }
 	    #$fastaHash{$acc} += 1;
	    
	    $line=<FASTA>;
	    chomp $line;
	    my $seq = $line;
	    $line = <FASTA>;
	    chomp $line;
	    while($line){
		$seq .= $line;
		$line = <FASTA>;
		chomp $line;
	    }
	    $fastaHash{$acc} = $seq; 
	} 

    }


    close(FASTA);
    return \%fastaHash;
}


# performEntrezQueryOnArray($accArr, $fastaOutputFile);
# Given and array of accessions (example is a missingAccArr), retrieve fasta sequences and append to outputFile
# Inputs: $accArr - reference to array of accessions 
#         $outputFile - file to append fasta sequences to
# Kristen Naegle
# November 1, 2007
sub performEntrezQueryOnArray($$){
    my $accArr = shift;
    my $outputFile = shift;

    my $tempFile = "tempFastaOut";
    my $tempFile2 = "temp_895467";
    if(!-e $outputFile){ system("touch $outputFile"); }
    if(!-e $tempFile){ system("touch $tempFile"); }
    if(!-e $tempFile2){ system("touch $tempFile2"); } 
    
    foreach my $acc (@$accArr){
	performEntrezQuery($acc, $tempFile);
	system("cat $tempFile $outputFile > $tempFile2");
	system("cp $tempFile2 $outputFile");


    }
    
    system("rm $tempFile");
    system("rm $tempFile2");


}

# This will request missing fasta sequences from ncbi, and then print these and other fasata sequences to the tempFastaFile so that you can run HMMDomain analysis only on this subset.
sub printRequestedFastaToFile($$$$){
    my $accHash = shift;
    my $missAcc = shift; #only want to get fasta sequences for those domains that are missing
    my $tempFastaFile = shift;
    my $fastaFile = shift;


    my $doublyTemp = "TEMP_TEMP_FASTA_FILE";
    #Get FASTA sequences 
    my $fastaHash = returnFastaHash($fastaFile);
    my $missingFasta = returnAccMissingFromHash(arrToHash($missAcc), $fastaHash);
    if(scalar(@$missingFasta > 1)){   
	print "Retrieving ". scalar(@$missingFasta)." Fasta sequences .... \n";
	performEntrezQueryOnArray($missingFasta, $tempFastaFile);
	# append these missing to the final fastaFile
	`cat $fastaFile $tempFastaFile > $doublyTemp`;
	`mv $doublyTemp $fastaFile`; #removes doubly temp at the same time as overwriting
	
	$fastaHash = returnFastaHash($fastaFile);
	#check for missing fasta sequences.
	$missingFasta = returnAccMissingFromHash(arrToHash($missAcc), $fastaHash);
# # if(scalar(@$missingFasta) > 0){
# # 	    print "ERROR in $0, you are still mising Fasta sequences in $fastaFile even after retrieving..don't know what's wrong\n";
# # 	    foreach my $f (@$missingFasta){
# # 		print "Missing: $f\n";
# # 	    }
# # 	    exit;
# # 	}
# # 	else{
# # 	    print "Have successfully retrieved missing Fasta sequences\n";
# # 	}
	
    }
    #now create a temp file that is the fasta sequences from the accessions collected. 
    my $translateHash = returnTranslateAccHash($fastaHash, $accHash);
    open(TFASTA, ">$tempFastaFile");
    foreach my $acc (keys %$translateHash){
	print TFASTA ">".$acc."\n".$fastaHash->{$acc}."\n\n";
	
    }
    close(TFASTA);
    #my $HMMTraining = "/data/knaegle/data/pfam/Pfam_fs";
    #createHMMResults($tempFastaFile, $domainFile, $HMMTraining);
    my $fastaFail = returnAccMissingFromHash($accHash, $fastaHash);
    return $fastaFail;

}


# ($errorFlag, \%domainHash, $source, $params, $version) = returnDomainHashForDB($tempFastaFile, $acc, $pvalueCutoff, $CPU);
# Grabs and writes a single accession to a fasta file, runs domain analysis, and returns the domainHash, which has keys 'start' 'stop', 'label', and 'p_value' 
# Inputs: $tempFastaFile - location of the fasta file with sequence
#         $acc - accession number (must be able to return fasta from entrez)
#         $pvalueCutoff - minimum pvalue allowed for predictions
#         $CPU - number of processors to use
# Outputs: $errorFlag - returns a 1 if there are proteins in dataFile for which there are noHMM predictions due to an error in either fasta capture or domain prediction
#          %domainHash - hash of arrays which represent all the domains in a protein with keys described above
#          $source - returns 'pfam'
#          $params - any params used - i.e. the pfam file used and pvalue cutoff
# Kristen Naegle
# March 9, 2008
sub returnDomainHashForDB($$$$){
    my ($tempFastaFile, $acc, $pvalueCutoff, $CPU) = @_; 

    my %domainHash;
    
    my $source = 'COMPUTED PFAM';
   # my $HMMTraining = "/data/knaegle/data/pfam/Pfam_fs";
    my $HMMRoot = $globalVars::PFAM_PATH;
    my $HMMTraining = $HMMRoot.'Pfam_fs';
    my ($version, $linkFile)  = returnLinkRelease($HMMTraining);
    my $params = "pval=$pvalueCutoff";
    my $errorFlag = 0;

    my $domainFile = "domain_database_temp";
    if(!-e $domainFile){`touch $domainFile`;}
    open(D, ">$domainFile"); close(D);


     my %accHash;
     $accHash{$acc} = 1;
# #     performEntrezQueryOnArray(\@accArr, $tempFastaFile);
# #     # Run HMM predictions (eventually move pfam training predictions to environment variable
    # write fasta file
    #writeFastaFileForProteinId($proteinId, $acc, $tempFastaFile);

 #   print "Creating domain predictions for new fasta sequences \n";
    
    createHMMResults($tempFastaFile, $domainFile, $HMMTraining, $CPU);
    
    # if file is empty then set error
    #print "DEBUG: Special variable before wc $?\n";
    my $count = `wc -l < $domainFile`;
    chomp $count;
    #print "DEBUG: after: $? and count=$count\n";
    if(!-e $domainFile or !$count){
	$errorFlag = 1;
	handleError('returnDomainHashForDB', "wc operation has error in HMMToolsError: $?", \@_);
	return ($errorFlag, \%domainHash, $source, $params, $version);
    }
    print "DEBUG:  Counting $domainFile, has $count lines\n";
    



    # Recheck to make sure everything was covered:
    my $hmmHash = returnHMMFileAccHash($domainFile);
    my $missAcc = returnAccMissingFromHash(\%accHash, $hmmHash);
    print "DEBUG: Missing these accessions: @$missAcc\n";
    if(scalar(@$missAcc) > 1){ # ??? really shouldn't it be > 1
	#print "ERROR:  Did not generate pfam predictions for accession $acc\n";
	$errorFlag = 1;
	handleError('returnDomainHashForDB', 'Did not generate pfam predictions', \@_);
	return ($errorFlag, \%domainHash, $source, $params, $version);
    }   
	
    #append HMM info to file;
    # new translatehash

    ###DEBUG, not getting all the domains and don't know why

    
    my $hmmHashPruned = pruneHMMHash($hmmHash, $pvalueCutoff);
    my $hmmHashRec = reconcileHMMDomains($hmmHashPruned);

  
    my $giTranslateHash = returnTranslateAccHash(\%accHash, $hmmHash);
    my $giKey = $giTranslateHash->{$acc};
    
    my($label, $start, $stop, $p_value) = returnDomainsForAcc($hmmHashRec, $giKey);
    push @{$domainHash{'start'}}, @$start;
    push @{$domainHash{'stop'}}, @$stop;
    push @{$domainHash{'label'}}, @$label;
    push @{$domainHash{'p_value'}}, @$p_value;
    print "Number of Domains in $0: ".scalar(@$start)."\n";
    return ($errorFlag, \%domainHash, $source, $params, $version);
}

# ($errorFlag, \@start, \@stop, \@label, \@p_value, $source, $params, $version) = returnDomainsByPrediction($dbh, $proteinId, $acc, $pvalueCutoff, $CPU);
# An extra level of indercition added so computation and parsing looks the same frmo the top function handleDomainsForProteinId.  Creates and writes fasta, then calls returnDomainHashForDB.
# Inputs: $dbh - database handle
#         $proteinId - id of protein table
#         $acc - accession used for translating hashes
#         $pvalueCutoff - cutoff on pvalue to prune 
#         $CPU - number of processes to use
# Outputs: $errorFlag- if 1, there was an error somwhere. Check the ERROR Log 
#          $start - ref. to array of start positions
#          $stop - ref. to array of stop positions
#          $label - ref. to array of labels
#          $p_value - ref to array of pvalues 
#          $source - source 
#          $params - params used
#          $version - version of file used to predict
# Kristen Naegle
# June 19, 2009
sub returnDomainsByPrediction($$$$$){
    my($dbh, $proteinId, $acc, $pvalueCutoff, $CPU) = @_;
    my $tempFastaFile = "fasta_database_temp";
    if(!-e $tempFastaFile){`touch $tempFastaFile`;}
    open(T, ">$tempFastaFile"); close(T);
    #acc ..here is where acc could be problem -- also returnDomainHashForDB should also return error
    writeFastaFileForProteinId($dbh, $proteinId, $acc, $tempFastaFile);
    my ($errorFlag, $dHashRef, $source, $params, $version) = returnDomainHashForDB($tempFastaFile, $acc, $pvalueCutoff, $CPU);
    

    my @start = @{$dHashRef->{'start'}};
    my @stop  = @{$dHashRef->{'stop'}};
    my @label = @{$dHashRef->{'label'}};
    my @p_value = @{$dHashRef->{'p_value'}};

    


    my $n = scalar(@start);
    #print "Num domains: $n\n";
    if($n != scalar(@stop) || $n != scalar(@label) || $n != scalar(@p_value)){
	$errorFlag = 1; 
	handleError('handleDomainsForProteinId/returnDomainsByPrediction', "Number of lables, starts, stops, and pvalues was not the same for accession=$acc, proteinID=$proteinId", \@_);
    }
    return ($errorFlag, \@start, \@stop, \@label, \@p_value, $source, $params, $version);
}

# (\@label, \@start, \@stop, \@p_value) = returnDomainsForAcc($hmmHash, $accKey);
# Returns the information important for db entry of a domain for a given accession from an HMM Hash
# Inputs: $hmmHash - read in from a domain file  (see returnDomainHashForDB);
#         $accKey - the key to the domain hash (need to translate a direct accesion to key form)
# Outputs: Arrays for each of the domains
#          \@label - labels of the domain found (assumes hash constructed by a cutoff already
#          \@start - start positions
#          \@stop - stop positions
#          \@p_value - pvalues for that domain hit
# Kristen Naegle
# March 9, 2008
sub returnDomainsForAcc($$){
    my $hmmHash = shift;
    my $accKey = shift;
 
    
    my(@label, @start, @stop, @p_value);

    foreach my $feature (@{$hmmHash->{$accKey}}){
	push @label, $feature->primary_tag;
	push @start, $feature->start;
	push @stop, $feature->end;
	push @p_value, $feature->score;
    }
    return (\@label, \@start, \@stop, \@p_value);

}

# Given a Pfam file, parse into a domain hash object
# Inputs: $pfam_file - an already calcuated pfam file from ftp.sanger.ac.uk
# Outputs: $domainArrHash - top key is accession number and points to an array of domain hashes Domian Hash has keys: 'gene', 'start', 'stop', 'description', 'pfam_id' and 'label'
# Kristen Naegle
# June 15, 2009
sub parsePfamFile($){
    my ($pfamFile) = @_;
    
    my %hash;
    
    open(P, $pfamFile) || die "Can't open $pfamFile for reading in parsePfamFile\n";
    # gather lines that represent a protein

    while(defined(my $line = <P>)){
	my @protein_lines;
	if($line =~ m/>/){
	    push @protein_lines, $line;
	    my $newLine = <P>;
	    while($newLine !~ /^\n|\r|\r\n/){
		push @protein_lines, $newLine;
		$newLine = <P>;
	    }
	  #  print "DEBUG: Found @protein_lines for a protien\n\n";
	    my ($acc, $version, $gene, $numAA, $domainHash) = parseProteinLines(\@protein_lines);
	    if(defined($hash{$acc})){
		print "ERROR: Domains for $acc found multiple times\n";
		
	    }
	    else{
		$hash{$acc} = $domainHash;
	    }

	} # end of new protein found
	
    }

    close(P);
    return \%hash;

}

# ($error, $domainHash) = returnDomainHashFromGlobalFile()
# Uses the global pfam file and returns a domain hash
# Outputs: $domainArrHash - top key is accession number and points to an array of domain hashes Domian Hash has keys: 'gene', 'start', 'stop', 'description', 'pfam_id' and 'label'
# Kristen Naegle
# June 19, 2009
sub returnDomainHashFromGlobalFile(){
    my ($pfamFile, $release) = returnGLOBALPfamFile();
    print "DEBUG: pfam file for parsing: $pfamFile\n";
    my $domainHash = parsePfamFile($pfamFile);
    return $domainHash;
}

# Given a Pfam file, parse into a domain hash object
# Inputs: $pfam_file - an already calcuated pfam file from ftp.sanger.ac.uk
# Outputs: $domainArrHash - top key is accession number and points to an array of domain hashes Domian Hash has keys: 'gene', 'start', 'stop', 'description', 'pfam_id' and 'label'
# Kristen Naegle
# June 15, 2009
sub parsePfamFile_makeSubsetFile($$$){
    my ($pfamFile, $outputFile, $species) = @_;
    my @species = keys(%$species);
    print "DEBUG: Species @species\n";

    open(P, $pfamFile) || die "Can't open $pfamFile for reading in parsePfamFile\n";
    # gather lines that represent a protein
    open(O, ">$outputFile") || die "Can't open $outputFile for writin in parsePfamFile_makeSubsetFile\n";
    print O "\n";
    while(defined(my $line = <P>)){
	my @protein_lines;
	if($line =~ m/>/){
	    #get the species and check to see if it's right
	    my @header = split(' ', $line);
	    if (scalar(@header) != 5){
		print "ERROR! parseProteinLines_makeSubsetFile assumes the fasta style line has five fields!!, it does not! with lines: @header \n";
		exit;
	    }
	    my ($gene, $j, $acc, $numAA, $j2)= @header;
	    $gene =~ m/>.+_(.+)/;
	    my $sp = $1;
	    $sp =~ s/ //;
	    #print "DEBUG: Species of line: $sp\n";
	    if(defined ($species->{$sp})){
		#print "\tDEBUG: YES\n";
		print O $line;
		my $newLine = <P>;
		while($newLine !~ /^\n|\r|\r\n/){
		    print O $newLine;
		    $newLine = <P>;
		}
		print O "\n"; #need spacer
	    } # end indeed in speces
	    
	
	} #end found a header line
    }
    close(P);
    close(O);
}


# $domainHas = parseProteinLines($proteinLines)
# GIven an array of protein lines from parse of a pfam file, return an array of hash objects (array is in order of start positions)
# Inputs: $proteinLines - reference to an array of lines that describe a single protein (first line is header and has accession, remaining lines show domain name, positions, etc.)
# Outputs: $domainHash - top key is count for a unique identifier, this points to anohter hash with keys:
#             start, stop, label, numRepeats, gene, description, pfam_id
# Kristen Naegle
# June 15, 2009
sub parseProteinLines($){
    my ($proteinLines) = @_;
    
    my %domainHash;

    # get the accessions from the first line and the number of amino acids (just do this now and return them if ever needed)
    my $header = $proteinLines->[0];
    my @header = split(' ', $header);
    if (scalar(@header) != 5){
	print "ERROR! parseProteinLines assumes the fasta style line has five fields!!, it does not! with lines: @header \n";
	exit;
    }
    my ($gene, $j, $acc, $numAA, $j2)= @header;
    $gene =~ s/>//;
    my $version;
    if($acc =~ m/\.(\d)/){
	$version = $1;
	$acc =~ s/\.$version//;
    }
   # print "DEBUG:x Gene: $gene\tAcc:$acc\tverstion:$version\tNUM_AA: $numAA\n";
    my $count = 0;   
    # now for remaining lines get details
    for(my $i=1; $i < scalar(@$proteinLines); $i++){
	my $line = $proteinLines->[$i];
	my @line = split(' ', $line);
#	print "DEBUG: Line: @line\n";
	my $domainName = $line[0];
	my $numRepeats = $line[1];
	#print "DEBUG: name: $domainName\trepeats: $numRepeats\n";
	#now find the name and description
	my ($pfam_id, $description, $subline);
	if($line =~ m/(P[FB]\d+\.?\d)(.+)  (.+)/){
	    $pfam_id = $1;
	    $description = $2;
	    $subline = $3;
	}
	

	my @pos = split(' ', $subline);
	if(scalar(@pos) != $numRepeats){
	    print "ERROR: Number of starts and stops found does not equal expected number of repeats for domain: $domainName and accession: $acc\n";
	    exit;
	}

	foreach my $pos (@pos){
	    $pos =~ m/(\d+)-(\d+)/;
	    my $start = $1;
	    my $stop = $2;
	    my %temp;
	    $temp{'start'} = $start;
	    $temp{'stop'} = $stop;
	    $temp{'label'} = $domainName;
	    $temp{'numRepeats'} = $numRepeats;
	    $temp{'gene'} = $gene;
	    $temp{'description'} = $description;
	    $temp{'pfam_id'} = $pfam_id;
	    $domainHash{$count} = \%temp;
	    $count += 1;
	}
	#now for that acc create an array of hashes based on position
	#temp: ignore ordering

    } # end for every protein

    return ($acc, $version, $gene, $numAA, \%domainHash);
}

# ($error, \@label, \@start, \@stop, \@p_value, \@pfam_id, \@description, $source, $params, $version) = returnDomainsForProteinParsed($domainHash, $acc);
# given the parsed hash of a full pfam file, return the arrays of start, stop, label, etc. for a swissprot accession
# Inputs: $domainHash - this is the hash created by parsePfamFile
#         $acc - better be the same kind of accession that is in the parsed file (swissprot, entrez, etc.)
# Outputs: $error - returns 1 if accession wasn't found in domainhash
#   The following are all ref. to arrays of the information for a domain. Position in array is the same for everything. start, stop, label, p_value, pfam_id, description
#     p_value is -1 to indicate it was not computed, but parsed
# Kristen Naegle
# June 16, 2009
sub returnDomainsForProteinParsed($$){
    my ($domainHash, $acc) = @_;
    my (@start, @stop, @label, @p_value, @pfam_id, @description);
    my $error = 0;
    my $source = "PARSED PFAM";
    my $params = "ALL";
    my ($pfamFile, $version) = returnGLOBALPfamFile();
       
    if($acc =~ m/^gi\|\d+/){ #entrez file doesn't have gi| in it
	$acc =~ s/^gi\|//;
    }
    if(defined($domainHash->{$acc})){
	my $hash = $domainHash->{$acc};
	foreach my $d (keys %$hash){
	    my $h = $hash->{$d};
	    push @start, $h->{'start'};
	    push @label, $h->{'label'};
	    push @stop, $h->{'stop'};
	    push @pfam_id, $h->{'pfam_id'};
	    push @description, $h->{'description'};
	    push @p_value, -1;
	}
    }
    else{
	$error = 1;
	return ($error, \@label, \@start, \@stop, \@p_value, \@pfam_id, \@description);   
    }

    my $n = scalar(@start);
    #print "Num domains: $n\n";
    if($n != scalar(@stop) || $n != scalar(@label) || $n != scalar(@p_value)){
	$error = 1; 
	handleError('handleDomainsForProteinId/returnDomainsForProteinParsed', "Number of lables, starts, stops, and pvalues was not the same for accession=$acc, Accession $acc", \@_);
    }

	return ($error, \@label, \@start, \@stop, \@p_value, \@pfam_id, \@description, $source, $params, $version);   
}

# ($error, $domainHash, $source, $params) = returnDomainHashForParsed($domainHash, $acc); 
# given the parsed hash of a full pfam file, return the domain hash which has keys to arrays for start, stop, label and pvalue
# Inputs: $domainHash - this is the hash created by parsePfamFile
#         $acc - better be the same kind of accession that is in the parsed file (swissprot, entrez, etc.)
# Outputs: $error - returns 1 if accession wasn't found in domainhash
#   The following are all ref. to arrays of the information for a domain. Position in array is the same for everything. start, stop, label, p_value, pfam_id, description
#     p_value is -1 to indicate it was not computed, but parsed
# Kristen Naegle
# June 16, 2009
sub returnDomainHashForParsed($$){
    my ($domainHash, $acc) = @_;
    my $error = 0;
    my $hash;
    if($acc =~ m/^gi\|\d+/){ #entrez file doesn't have gi| in it
	$acc =~ s/^gi\|//;
    }

    if($domainHash->{$acc}){
	$hash = $domainHash->{$acc};
    }
    else{
	$error = 1;
    }

    return ($error, $hash);
}

# ($pfamFile, $release) = returnGLOBALPfamFile();
# Return the global pfam file ( a link to current release version )
# Outputs: $pfamFile - the file that is a symbolic link. Here so that I can move the file around as needed 
#          $release - the version number of the pfam release
# Kristen Naegle
# June 19, 2009
sub returnGLOBALPfamFile(){
    my $pfamRoot = $globalVars::PFAM_PATH;
#    my $pfamFile = "/data/knaegle/data/pfam/swisspfam_subset";
    my $pfamFile = $pfamRoot."swisspfam_subset";
    my ($RELEASE, $file) = returnLinkRelease($pfamFile);
    return($pfamFile, $RELEASE);

    
}

# makeSpeciesFile($pfamFile, $outputFile);
# Given a pfam file and an output to write the parsed, parse it into species based on this function
# Inputs: $pfamFile - full swisspfam file source
#         $outputFile - destination of parsed file
# Kristen Naegle
# Nov. 23, 2009
sub makeSpeciesFile($$){
    my ($pfamFile, $outputFile) = @_;


    my @species = ('ARATH', 'BOVIN', 'DROME', 'CHICK', 'HUMAN', 'MOUSE', 'RAT', 'YEAST', 'SCHPO');
# arabidopsis_thaliana		ARATH
# bos_taurus			BOVIN
# caenorhabditis_elegans	*not very many of these
# drosophila_melanogaster	DROME
# gallus_gallus			CHICK
# homo_sapiens			HUMAN
# mus_musculus			MOUSE
# rattus_norvegicus		RAT
# saccharomyces_cerevisiae	YEAST
# schizosaccharomyces_pombe	SCHPO

# create a hash from teh array 
    my %species;
    foreach my $s (@species){

	$species{$s} = 1;
    }

    parsePfamFile_makeSubsetFile($pfamFile, $outputFile, \%species);


}


# retrive Pfam files Pfam_fs and swisspfam from Pfam server
# Places them in direcotry based on global and sets up symbolic link to files
# NOT FINISHED - Found out new HMM Files built by HMMER3, not currently parsed by bioperl..will have to wait 
# Kristen Naegle
# Nov. 23, 2009
sub retrievePfamFiles(){

    my $url = 'ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/';
    `wget $url`;
    my $tempFile = "index.html";
    if(!-e $tempFile){
	print "ERROR: Could not retrieve Pfam files\n";
	return -1;
    }
    my $latestRelease = parsePfamDirIndexForLatestRelease($tempFile);
    my $baseURL = $url."Pfam".$latestRelease."/";
    my $swissPfam = "swisspfam.gz";
    my $pfam_fs = "pfam_fs.gz"; #### THIS HAS CHANGED In port to HMMER3

    #remove that index.html so we don't accidentally read it again at a later date.
    # `rm $tempFile`;
    

}
# $release = parsePfamDirIndexForLatestRelease($indexFile)
# Given the index listing of Pfam release, retrieve the latest release number
# Inputs: $indexFile - a wget local copy of a Pfam fetch from root ftp
#                     example:see retrievePfamFiles
# Outputs: $release = latest release number: e.g. 23.0
# Kristen Naegle
# November 23, 2009
sub parsePfamDirIndexForLatestRelease($){
    my ($file) = @_;
    
    open(I, $file) || die "Can't open $file for reading\n";
    
    my $line = <I>;
    my @releases = (0);
    while(defined($line = <I>)){
	if($line =~ m/Pfam(\d+\.\d)/){
	    my $release = $1;
	    push @releases, $release;
	   # print "DEBUG: found $release\n";
	}

    }
    my @sortedReleases = sort {$a <=> $b} @releases;
    my $latest = $sortedReleases[$#sortedReleases];
    close(I);
    return $latest;
}

1;
