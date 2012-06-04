use strict;
use warnings;
use DBTools::queryTools;
use fileTools;
use DBTools::resultParsing;
#use Math::Counting ':long';
use POSIX;
use Getopt::Long;
use dataIO;
use errorHandling;

# (\%hash, $label) = returnClusterHash($inputFile, $clusterSetNum);
# Returns a cluster hash of ms_ids in a file
# Inputs: $inputFile - tab separated file with AT Least an MS.id column and a cluster:setNum:label column
#         $clusterSetNum - number of cluster for which to create a hash
# Outputs: \%hash - hashsref with keys equal to cluster numbers and values is an array of MS ids
#          $label - label of that clusterSet
# June 18, 2008
# Kristen Naegle
sub returnClusterHash($$){
    my ($inputFile, $clusterSetNum) = @_;
    my $clusterCols = returnColumnNumberArr($inputFile, "cluster:$clusterSetNum");
    my $h = returnHeader($inputFile);
    my @header = split('\t', $h);
    my $clusterCol = -1;
    foreach my $col (@$clusterCols){
	if($header[$col] =~ m/cluster:$clusterSetNum(\b|\D)/){
	    $clusterCol = $col;
	}

    }
    
    if($clusterCol < 0){
	print "ERROR: Cannot find correct column\n";
	my %h;
	return \%h;
    }
    #print "ClusterCol is $clusterCol\n";
    my $MSCol = returnColumnNumber($inputFile, 'MS.id');
    my @MSIds;
    open(CF, $inputFile);
    my $header = <CF>;
    chomp $header;
    my @h = split('\t', $header);
    my $label = $h[$clusterCol];
    $label =~ s/\:cluster\:$clusterSetNum(\b|\D)//;
    my %hash;
    while(defined(my $line = <CF>)){
	chomp $line;
	my @line = split('\t', $line);
	my $MS_id = $line[$MSCol];
	push @MSIds, $MS_id;
	my $cluster = $line[$clusterCol];

	if(not defined $hash{$cluster}){
	    $hash{$cluster} = [];
	}
	push @{$hash{$cluster}}, $MS_id;
    }
    return (\%hash, $label, \@MSIds);
}

# ($proteinIdArr, $proteinCount) = convertMSIdsToProteinIds($dbh, $MSIdArr, $uniquify)
# Given an array of MS ids (i.e. foreground or background), return the proteins in that group
# Inputs: $dbh - database handle
#         $MSIdArr - reference to array of MS ids
#W        $uniquify - flag 1: reduce multiple mappings of MS to proteins to single protein 0: don't remove proteinId redundancy
# Outputs: $proteinIdArr - reference to array of protein ids
#          $proteinCount - total number of proteins returned (no redundancies)
# Kristen Naegle
# July 1, 2008
# Mod. Oct. 20, 2008 - to include uniquify
sub convertMSIdsToProteinIds($$$){
    my ($dbh, $MSIds, $uniquify) = @_;
    my $MSstring = createSTHArrayVal($MSIds);

    my $sth;
    if($uniquify){
	$sth= $dbh->prepare("SELECT DISTINCT protein_id from MS where MS.id IN $MSstring");
    }
    else{
	$sth= $dbh->prepare("SELECT protein_id from MS where MS.id IN $MSstring");
    }
    $sth->execute();
    my $proteinIdArr = returnArrayOfResultsOnCol($sth, 0);
    my $proteinCount = scalar(@$proteinIdArr);
    return($proteinIdArr, $proteinCount);
}

# ($phosphopepIdArr, $phosphopepCount) = convertMSIdsToPhosphoepIds($dbh, $MSIdArr)
# Given an array of MS ids (i.e. foreground or background), return the proteins in that group
# Inputs: $dbh - database handle
#         $MSIdArr - reference to array of MS ids
#         $uniquify - flag if 1: collapse to distinct phosphopep_ids otherwise return redundant ids if mapped redundantly from MSIdArr
# Outputs: $phosphopepIdArr - reference to array of phosphopep ids
#          $phosphopepCount - total number of phosphopeptides returned
# Kristen Naegle
# July 1, 2008
# Mod Oct. 22, 2008
sub convertMSIdsToPhosphopepIds($$$){
    my ($dbh, $MSIds, $uniquify) = @_;
    my $MSstring = createSTHArrayVal($MSIds);
    my $sth;
    if($uniquify){
	$sth = $dbh->prepare("SELECT distinct phosphopep_id from MS_phosphopep where MS_id IN $MSstring");
    }
    else{
	$sth = $dbh->prepare("SELECT phosphopep_id from MS_phosphopep where MS_id IN $MSstring");
    }
    $sth->execute();
    my $phosphopepIdArr = returnArrayOfResultsOnCol($sth, 0);
    my $phosphopepCount = scalar(@$phosphopepIdArr);
    return($phosphopepIdArr, $phosphopepCount);
}

# $string = createSTHArrayVal($arrRef);
# Convert an array to a query IN statement (i.e. returns string with '(val1, val2, .., valN)'
# Inputs: $arrRef - reference to array containing items to convert to string
# Outputs: $string - converted string, with appropriate single quotes
# Kristen Naegle
# July 1, 2008
sub createSTHArrayVal($){
    my ($array) = @_;
    my $string = '(';
    my $last = scalar(@$array);
    for(my $i=0; $i<$last-1; $i++){
	$string .= "'".$array->[$i]."'".',';
    }
    $string .= "'".$array->[$last-1]."')";
    return $string;

}

# $hasRef = returnGOHashForMSIdArr($dbh, $proteinIdArr, $uniquify);
# For an array of MS Ids, return a hash of GO hashes, where top level key is GO term type ('F', 'C', or 'P') and the next level key is the term  
# Inputs: $dbh - database handle
#         $proteinIdArr - reference to array of proteins for which to include in GO fetch
#         $uniquify - flag to return only unique protein_id counts before requesting GO terms.
# Outputs: $hashRef - reference to hash of hashes. Top key is aspect and hash{$aspect} is than a returnGOAspectHash 
# Kristen Naegle
# July 1, 2008
sub returnGOHashForMSIdArr($$$){
    my ($dbh, $MSIdArr, $uniquify) = @_;
    my %hash;
    my @aspects = ('C', 'F', 'P');
    foreach my $aspect (@aspects){
	my $hashRef = returnGOAspectHash($dbh, $MSIdArr, $aspect, $uniquify);
	$hash{$aspect} = $hashRef;
    
    }
    return \%hash;
}
# $hashRef = returnGOAspectHash($dbh, $proteinIdArr, $aspect, $uniquify)
# For a specific aspect, returns a hash with key being term and value being an array containing number of counts of that term
# Inputs: $dbh - database handle
#         $MSIdArr - reference to array of MS Ids for which to include in GO fetch (get through protein table).
#         $aspect - the particular aspect you're looking for 'C', 'F' or 'P'
#         $uniquify - flag to return only unique protein_id counts before requesting GO terms.
# Outputs: $hashRef - reference to hash of GO terms and values count
# Kristen Naegle
# July 1, 2008
# Changed July 2, 2008 to account for NULL values in fgnd and bgnd
# Mod. October 20, 2008 to allow for UNIQUIFY flag 
sub returnGOAspectHash($$$$){
    my ($dbh, $MSIdArr, $aspect, $uniquify) = @_;
    if($aspect ne 'C' && $aspect ne 'F' && $aspect ne 'P'){
	handleError('returnGOTypeHash', 'Type must be C, F or P', \@_);
    }
    my $MSIdString = createSTHArrayVal($MSIdArr);
    my $sth;
    if($uniquify){
	$sth = $dbh -> prepare("SELECT term, count(*) from protein_GO join GO on protein_GO.GO_id=GO.id where protein_id IN (select distinct protein_id from MS where MS.id IN $MSIdString) and aspect=? group by term");
    }
    else{
	$sth = $dbh -> prepare("SELECT term, count(*) from MS join protein_GO on MS.protein_id=protein_GO.protein_id join GO on protein_GO.GO_id=GO.id where MS.id IN $MSIdString and aspect=? group by term");
    }
    my $keyCol = 0;
    my @valCols = (1);
    $sth->execute($aspect);
    my $hashRef = returnHashResult($sth, 0, \@valCols);

    # Add to hash NULL Terms, the proteins that were missing that -- redo $sth execution as determined by uniquify above.

   #  my $sthNull = $dbh -> prepare("SELECT protein_id, count(*) from protein_GO join GO on protein_GO.GO_id=GO.id where protein_id IN $proteinString and aspect=? group by protein_id");
    my $sthNull = $dbh->prepare("SELECT MS.id from MS join protein_GO on MS.protein_id=protein_GO.protein_id join GO on protein_GO.GO_id=GO.id where MS.id IN $MSIdString and aspect=? group by MS.id");  #doesn't matter if uniquify or not..both will return same number of MS.ids that had annotation term
    $sthNull->execute($aspect);
    my $idsInResult = returnArrayOfResultsOnCol($sthNull,0);
    if($idsInResult->[0] == -1){
	push @{$hashRef->{'---'}}, scalar(@$MSIdArr);
    }
    else{
	my ($union, $intersection, $diff) = returnArrayStats($MSIdArr, $idsInResult);
	push @{$hashRef->{'---'}}, scalar(@$diff);
    }
    return $hashRef;
}

# Return Interesection, and diff of two arrays
# (\@union, \@intersection, \@difference) = returnArrayStats($array1, $array2)
# Inptus: $array1 - reference to one array of objects
#         $array2 - reference to second array of objects
# Outputs: \@union - the union of the two arrays
#          \@interesection - interesection of the two arrays
#          \@difference - difference of the two arrays
# Kristen Naegle
# July, 2008
sub returnArrayStats($$){
    my ($array1, $array2) = @_;
    my (@union, @intersection, @difference);
    @union = @intersection = @difference = ();
    my %count = ();
    foreach my $element (@$array1, @$array2) { 
	$count{$element}++ 
	}
    foreach my $element (keys %count) {	
	push @union, $element;	
	push @{ $count{$element} > 1 ? \@intersection : \@difference }, $element;
    }  
    return (\@union, \@intersection, \@difference);
}


# $hash = convertArrToHash($arr)
# Convert an array to hash with value equal to count and keys are objects in array
# Inputs: $arr - ref to array of objects
# Outputs: $hash - ref. to hash of keys with objects and value equal to the number of times they appear.
# Kristen Naegle
# ? 2008
sub convertArrToHash($){
    my ($arr) = @_;
    my %hash1;
    foreach my $arrE (@$arr){
	if(not defined ($hash1{$arrE})){
	    $hash1{$arrE} = 0;
	}
	$hash1{$arrE} += 1;
    }
    return \%hash1;
}
# $hashRef = returnDomainAspectHash($dbh, $MSIdArr, $aspect, $domainThreshold, $uniquify)
# For a specific source, returns a hash with key being term and value being an array containing number of counts of that term. 
# Inputs: $dbh - database handle
#         $MSIdArr - ref to array of MS.ids
#         $aspect - the particular source (e.g. pfam)
#         $domainThreshold - the threshold for the prediction required to count it
#W        $uniquify - flag: if 0 allow protein redundancy, else unquify protein ids from MS Id
# Outputs: $hashRef - reference to hash of domain labels and values count
# Kristen Naegle
# July 1, 2008
# Mod Oct. 20, 2008 - to take MS ids and uniquify option, also added a domain_threshold cutoff
sub returnDomainAspectHash($$$$$){
    my ($dbh, $MSIdArr, $aspect, $domainThreshold, $uniquify) = @_;
    my $MSIdString = createSTHArrayVal($MSIdArr);
    
    if($uniquify){
	print "ERROR!! Uniquify not enabled in returnDomainAspectHash\n";
	$uniquify = 0;
    }
    # First get unique domain.labels in the MS Id set. Then ask how many proteins have that label in the set. (Since many proteins can have domain repeats..not interested in the repeats). 
    my $sth;
    $sth = $dbh->prepare("SELECT domain.label from MS join domain on MS.protein_id=domain.protein_id where MS.id IN $MSIdString and source REGEXP ? and p_value <= ? group by label");
    $sth->execute($aspect, $domainThreshold);
    my $domainArr = returnArrayOfResultsOnCol($sth, 0);

# for each hashRef (unique domain label) find out how many distinct proteins have that label in a group
    my %hash;
    my $sth2; 
#    if($uniquify){
#	# unique means select protein_ids with a domain label
#	$sth2 = $dbh->prepare("SELECT MS.protein_id from MS join domain on MS.protein_id=domain.protein_id where MS.id IN $MSIdString  and source=? and label=? and p_value <= ? group by MS.protein_id");
#    }
#    else{
#	$sth2 = $dbh->prepare("SELECT MS.protein_id from MS join domain on MS.protein_id=domain.protein_id where MS.id IN $MSIdString  and source REGEXP ? and label=? and p_value <= ? group by MS.id");
#    }
    $sth2 = $dbh->prepare("SELECT MS.protein_id from MS join domain on MS.protein_id=domain.protein_id where MS.id IN $MSIdString  and source REGEXP ? and label=? and p_value <= ? group by MS.id");
    
# Keep running array of proteinIds that returned a domain..then ask (using SELECT not in for missing proteins).  -- Also asking here the count of proteins with that label (with min confidence of assignment)
    my @proteinIdArr;
    foreach my $label (@$domainArr){
	$sth2->execute($aspect, $label, $domainThreshold);
	my $proteinResults = returnArrayOfResultsOnCol($sth2, 0);
	my $proteinCount = scalar(@$proteinResults);
	push @proteinIdArr, @$proteinResults;
	push @{$hash{$label}}, $proteinCount;
    }
    my $uniqueProtHash = convertArrToHash(\@proteinIdArr);
    my @proteinsWLabel = keys %$uniqueProtHash;
    my $proteinIdStr = createSTHArrayVal(\@proteinsWLabel);
    my $sth3 = $dbh->prepare("SELECT DISTINCT protein_id from MS where MS.id IN $MSIdString and protein_id NOT IN $proteinIdStr");
    $sth3->execute();
    my $missingProteins = returnArrayOfResultsOnCol($sth3, 0);
    push @{$hash{'---'}}, ((scalar(@$missingProteins) == 1 && $missingProteins->[0] == -1)?0:scalar(@$missingProteins));  #have to make sure that missingProteins doesn't have a size of 1 but first value==-1 which means result returned empty
    return \%hash;
}


# $hasRef = returnDomainHashForMSIdArr($dbh, $MSIdArr, $domainThreshold, $uniquify);
# For an array of proteins, return a hash of domains
# Inputs: $dbh - database handle
#         $MSIdArr - reference to array of MS Ids for which to include in protein Domain fetch
#         $domainThreshold - p-value cutoff for counting domain of protein
#         $uniquify - flag: if 0 include redundant proteins, else reduce to unique protein set from MS ids
# Outputs: $hashRef - reference to hash of hashes. Top key is aspect (domain) and hash{$aspect} is than a returnDomainAspectHash 
# Kristen Naegle
# July 1, 2008
# Mod Oct 20, 2008 - to include unquify term
sub returnDomainHashForMSIdArr($$$$){
    my ($dbh, $MSIdArr, $domainThreshold, $uniqify) = @_;
    my %hash;

    my @aspects = ('PFAM');
    foreach my $aspect (@aspects){
	my $hashRef = returnDomainAspectHash($dbh, $MSIdArr, $aspect, $domainThreshold, $uniqify);
	$hash{$aspect} = $hashRef;
    
    }
    return \%hash;
}
# $hasRef = returnPepPredHashForMSIdArr($dbh, $MSIdArr, $scansiteCutoff, $uniquify);
# For an array of proteins, return a hash of GO hashes, where top level key is GO term type ('F', 'C', or 'P') and the next level key is the term  
# Inputs: $dbh - database handle
#         $proteinIdArr - reference to array of proteins for which to include in GO fetch
#         $uniquify - flag 1: reduce phosphopep redundancies from MS ids 0: include all phosphopeps
# Outputs: $hashRef - reference to hash of hashes. Top key is aspect and hash{$aspect} is than a returnGOAspectHash 
# Kristen Naegle
# July 1, 2008
# Modified October 15, 2008 - to include pelm Kinase annotations
# Mod. Oct. 20, 2008 - to include uniquify term.
sub returnPepPredHashForMSIdArr($$$$){
    my ($dbh, $MSIdArr, $cutoff, $uniquify) = @_;
    my %hash;
    my @aspects = ('scansite_kinase', 'scansite_bind', 'pelm_kinase');
    foreach my $aspect (@aspects){
	#print "looking for: $aspect\n";
	my $hashRef = returnPepPredAspectHash($dbh, $MSIdArr, $aspect, $cutoff, $uniquify);
	#print "FOUND ".scalar(keys %$hashRef)."predictions\n";
	$hash{$aspect} = $hashRef;
    
    }
    return \%hash;
}

# $hashRef = returnPepPredAspectHash($dbh, $phosphopepIdArr, $source, $cutoff, $uniquify)
# For a specific source, returns a hash with key being term and value being an array containing number of counts of that term
# Inputs: $dbh - database handle
#         $proteinIdArr - reference to array of proteins for which to include in GO fetch
#         $aspect - the particular aspect you're looking for 'C', 'F' or 'P'
# Outputs: $hashRef - reference to hash of GO terms and values count
# Kristen Naegle
# July 1, 2008
# Modified Sept. 22, 2008 - Now returns a --- field to count for nonreturns due to threshold cutoff. Also, didn't originally return ~~~ field since that was under scansite for both fields
# Mod. Oct. 20, 2008 - To include uniquify. Also to simplify query, changed the countLabel to ignore 
sub returnPepPredAspectHash($$$$$){
    my ($dbh, $MSIdArr, $aspect, $cutoff, $uniquify) = @_;

    my $MSIdString = createSTHArrayVal($MSIdArr);

    #handle the fact that scansite has a score cutoff but pelm does not 
    my $andScoreString;
    if($aspect =~ 'scansite'){
	$andScoreString = "and score <= $cutoff";
    }
    else{  #case is pelm instead and want a dummy score string
	$andScoreString = "and value is NOT NULL";
    }
    
    my $sth;
    my $sthNull;
    if($uniquify){
	$sth = $dbh -> prepare("SELECT value, count(*) from phosphopep_prediction where phosphopep_id IN (SELECT DISTINCT phosphopep_id from MS_phosphopep where MS_id IN $MSIdString) and source=? $andScoreString group by value");
	
    }
    else{  #not unique
	$sth = $dbh-> prepare("SELECT value, count(*) from MS_phosphopep join phosphopep_prediction on MS_phosphopep.phosphopep_id = phosphopep_prediction.phosphopep_id where MS_id IN $MSIdString and source=? $andScoreString group by value");

	
    }
    $sth->execute($aspect);
    my $keyCol = 0;
    my @valCols = (1);
    my $hashRef = returnHashResult($sth, 0, \@valCols);
    
    #handle unlabeled - Since one MS can map to many phosphopeps..need to ask how many phosphopeps did not have a label (under conditions) had no label. Now ignoring uniquify..uniquify is just wrong. 
    $sthNull = $dbh->prepare("SELECT MS_phosphopep.phosphopep_id from MS_phosphopep where MS_id IN $MSIdString and phosphopep_id NOT IN (SELECT MS_phosphopep.phosphopep_id from MS_phosphopep join phosphopep_prediction on MS_phosphopep.phosphopep_id=phosphopep_prediction.phosphopep_id where MS_id IN $MSIdString and source=? $andScoreString)");

   # $sthNull = $dbh->prepare("SELECT DISTINCT MS_id from MS_phosphopep join phosphopep_prediction on MS_phosphopep.phosphopep_id=phosphopep_prediction.phosphopep_id where MS_id IN $MSIdString and source=? $andScoreString");
    $sthNull->execute($aspect);
    my $MSIdsWithNoPred = returnArrayOfResultsOnCol($sthNull, 0);
    if($MSIdsWithNoPred->[0] == -1){
	push @{$hashRef->{'---'}}, 0;
    }
    else{ 
	push @{$hashRef->{'---'}}, scalar(@$MSIdsWithNoPred);
    }
    
    return $hashRef;
}

# ($hasRef, $countHash) = returnSiteDomainHashForMSIdArr($dbh, $MSIdArr, $uniquify);
# For an array of MS Ids find the phosphopep->pfam_sites
# Inputs: $dbh - database handle
#         $MSIdArr - reference to array of MS Ids for which to include in siteDomain Fetch
#         $uniquify - flag 1: remove phosphopep redundancy 0: don't ignore redundancy
# Outputs: $hashRef - reference to hash of hashes. Top key is aspect and hash{$aspect} is than a returnGOAspectHash 
#          $count - the number of pfam_site objects (for the most part these will be 1:1 with MS ids)
# Kristen Naegle
# July 1, 2008
# Mod Oct 2008 - to include uniquify
sub returnSiteDomainHashForMSIdArr($$$){
    my ($dbh, $MSIdArr, $uniquify) = @_;
    my %hash;

    my @aspects = ('pfam_site');
    my ($hashRef, $count);
    foreach my $aspect (@aspects){
	($hashRef, $count) = returnSiteDomainAspectHash($dbh, $MSIdArr, $uniquify);
	$hash{$aspect} = $hashRef;
    
    }
    return (\%hash, $count);
}
# ($hashRef, $count) = returnSiteDomainAspectHash($dbh, $MSIdArr, $uniquify)
# For a specific source, returns a hash with key being term and value being an array containing number of counts of that term
# Inputs: $dbh - database handle
#         $MSIdArr - reference to array of MS Ids for which to include in siteDomain Fetch
#         $aspect - pfam - this is the only one available in the format here
#         $uniquify - flag 1: remove phosphopep redundancy 0: don't ignore redundancy
# Outputs: $hashRef - reference to hash of siteDomain values and count
#          $count - the total number of objects counted for their pfam_site - since mostly this should be 1:1 correspondance with MS.id but if there is something at a boundary then will be larger
# Kristen Naegle
# July 1, 2008
# Mod Oct 2008 - to include uniquify. Also, removed aspect since this isn't used. I'm basing this on pfam_site in phosphopep
# Mod May 11, 2009 - to count MS ids uniquely (to remove pfam_site bias) unless a tryptic fragment bounds more than one domains. (UNIQUIFY option no longer works
sub returnSiteDomainAspectHash($$$){
    my ($dbh, $MSIdArr, $uniquify) = @_;
#    my $MSIdString = createSTHArrayVal($MSIdArr);
    my %countHash;
    if($uniquify){
	print "WARNING: This option doesn't work (uniquify in returnSiteDomainAspectHash)\n";
	print "Uniquify flag: $uniquify\n";
    }
    my $count = 0;
#   #  my $unlabeled;
# #     if($uniquify){
# # 	$sth = $dbh -> prepare("SELECT pfam_site, count(*) from phosphopep where phosphopep.id IN (SELECT DISTINCT phosphopep_id from MS_phosphopep where MS_id IN $MSIdString) group by pfam_site");
# #     }
# #     else{
# # 	$sth = $dbh->prepare("SELECT pfam_site, count(*) from MS_phosphopep join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id IN $MSIdString group by pfam_site");
# #     }
    my $sth = $dbh->prepare("SELECT pfam_site from phosphopep join MS_phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_phosphopep.MS_id=?");
    foreach my $MSId (@$MSIdArr){
	$sth->execute($MSId);
	my $results = returnArrayOfResultsOnCol($sth, 0);
	# uniquify the results
	my $unique = returnUniqueArray($results);
	foreach my $pfam_site (@$unique){
	    if(not defined $countHash{$pfam_site}){
		my @arr = (0);
		$countHash{$pfam_site} = \@arr;
	    }
	    $countHash{$pfam_site}->[0] += 1;
	    $count += 1;
	    
	}
	#stuff each into the hash

    }

# #     my $keyCol = 0;
# #     my @valCols = (1);
# #     $sth->execute();
# #     my $hashRef = returnHashResult($sth, 0, \@valCols);
#   #  $unlabeled = $hashRef->{'~~~'}->[0];
 #  # return $hashRef;
    return (\%countHash, $count);
}

# \@peps = returnPepAlignedForMSIdArr($dbh, $MSIdArr, $unquify);
# 
# Inputs: $dbh - database handle
#         $MSIdArr - ref. to array of MS ids 
#         $uniquify - flag to determine if phosphopeps should be non-redundant before getting aligned peptides
# Outputs: \@peps - ref to array of aligned peptides
# Kristen Naegle
# August 15, 2008
# Mod. Oct. 20, 2008 - to use MS Id array instead and offer uniquify flag
sub returnPepAlignedForMSIdArr($$$){
    my ($dbh, $MSIdArr, $uniquify) = @_;
    my @peps;
    my $sth;
    my $MSIdString = createSTHArrayVal($MSIdArr);
    if($uniquify){
	$sth = $dbh->prepare("SELECT pep_aligned from phosphopep where phosphopep.id IN (SELECT DISTINCT phosphopep_id from MS_phosphopep where MS_id IN $MSIdString)");    }
    else{
	$sth = $dbh->prepare("SELECT pep_aligned from MS_phosphopep join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id where MS_id IN $MSIdString");
    }
    $sth->execute();
    my $pepArr = returnArrayOfResultsOnCol($sth, 0);
    
    foreach my $pep (@$pepArr){
	if(length($pep) != 15){
	    
	    $pep = returnAlignedSpace($pep);
	}
	push @peps, $pep;
    }
    
    return \@peps;
}

# ($enrichmentHashRef, $numHypTestedHash) = calculateEnrichmentFromHash($dbh, $fgndSize, $bgndSize, $fgndHash, $bgndHash, $minCount, $threshold, $countLabeled)
# Given a count hash (with two top level keys, the metric and specific metric) and value equal to label followed by an array with first value equal to the count of that label, convert this to an enrichment hash where value is now a p-value.  Do it only if number of counts in foreground is >= $minCount and report only if p_value <= $threshold
#Inputs: $dbh - database handle
#        $fgndHash - count hash as described for the foreground
#        $bgndHash - count hash as described for background
#        $minCount - min. number of times label must appear in fgnd in order to calculate pval
#        $threshold - min p-value of label required to report it in enrichment hash
#        $countUnlabeled - if 1 then will include all in count of bgnd and fgnd
#Outputs: $enrichmentHash - same as input, but with pvalue
#         $numTested - need to return the number of hypotheses tested.. Same top key structure, with value equal to the total number of hypotheses tested (e.g. {'F} -> 200)
#Kristen Naegle
#October 20, 2008
# Updated November 25, 2008 - to include correction of $threshold by number of tests run.
sub calculateEnrichmentFromHash($$$$$$$$){
    my ($dbh, $fgndSize, $bgndSize, $fgndHash, $bgndHash, $minCount, $threshold, $countUnlabeled) = @_;
 #   print "Count Unlabeled: $countUnlabeled\n";
    my %pValTopHash;
    my %numTested;
    foreach my $aspect (keys %$fgndHash){
	my $fgnd_aspect = $fgndHash->{$aspect};
	my $bgnd_aspect = $bgndHash->{$aspect};
	my ($n, $N);

	if($countUnlabeled){
	    $n = $fgndSize;
	    $N = $bgndSize;
	}
	else{
	    $n = returnCountMinUnlabeled($fgnd_aspect, $fgndSize);
	    $N = returnCountMinUnlabeled($bgnd_aspect, $bgndSize);
	}
	printLineToLog("Aspect: $aspect\nForeground Size: $n\tBackground Size: $N\n");
	
	# create this hash the first time around..then remove terms below a given p-value adjusted by an input Type.  
	my %pValHash; #second level hash

	my $numHypothesesTested = 0; # zeroed for each aspect
	foreach my $term (keys %$fgnd_aspect){
	    if(!$countUnlabeled){
		if($term eq '~~~' or $term eq '---'){
		    next;
		}
	    }
	    my $fgnd_count = $fgnd_aspect->{$term}->[0];
	    my $bgnd_count = $bgnd_aspect->{$term}->[0];
	    if($fgnd_count >= $minCount){
#		print "Testing $term   k=$fgnd_count\tK=$bgnd_count\tn=$n\tN=$N\n";
		my $pval;
		if($fgnd_count > $n or $bgnd_count > $N or $n > $N or $fgnd_count > $bgnd_count){
		    handleError('calculateEnrichmentFromHash', "Your hypergeometric variables are wrong in testing ASPECT: $aspect and TERM: $term\n k=$fgnd_count\tn=$n\tK=$bgnd_count\tN=$N!!", \@_);
		    print "k=$fgnd_count\tK=$bgnd_count\tn=$n\tN=$N\n";
		    $pval = -1; 
		}
		else{
		    $numHypothesesTested += 1;
		    $pval = hyper($n, $N, $fgnd_count, $bgnd_count);
		}
		$pValHash{$term} = $pval; #stick them all in then go through and remove terms accroding to correction.
	       
	    }
	    
	}
	$pValTopHash{$aspect} = \%pValHash;
	$numTested{$aspect} = $numHypothesesTested;
    } # end the initial stuff of pvalues in. 



    return (\%pValTopHash, \%numTested);
}

# $unlabeledCount = returnCountMinUnlabeled($hash, $N);
# Removes from the count passed in the unlabeled features
# Inputs: $hash - a feature with labels equal to a metric label and values an array of count (really just the first term is used)
#         $N the number of objects you want to modify by the unlabeled count
# Outputs: $unlabeledCount - the count with unlabeled objects subtracted.
# Kristen Naegle
# Oct 21, 2008
sub returnCountMinUnlabeled($$){
    my ($hash, $N) = @_;
    my $unlabeledCount = $N;
    if($hash->{'~~~'}){
	$unlabeledCount -= $hash->{'~~~'}->[0];
    }
     if($hash->{'---'}){
	$unlabeledCount -= $hash->{'---'}->[0];
    }
    return $unlabeledCount;
}

# ($hashRef, $correctedHashRef, $numTest, $fgndHash) = calculateGOEnrichment($dbh, $fgnd, $bgnd,  $countUnlabeled, $minCount, $threshold, $CORRECTION_TYPE, $alpha);
# For a list of foreground MS Ids and a list of background MS Ids, report enrichment for terms occuring n times or more
# Returns a hash of hashes like the input foreground, except instead of count, there is a calculated probability value
# Inputs: $dbh - database handle
#         $fgnd - reference to array of MS ids of foreground
#         $bgnd - reference to array of MS ids of bacgkround
#         $countUnlabeled - flag: 1 - include unlabeled in calculations 0:exclude ~~~ and --- unlabeled 
#         $minCount - minimum count required to calculate p-value
#         $threshold - p-value threshold for reporting GO terms
#         $CORRECTION_TYPE - type of Multiple Hypotehsis testing correction
#         $alpha - target correction p-value
# REMOVED:        $uniquify - flag 1: condense protein ids from multiple MS id mappings 0: don't uniquify protein ids
# Outputs: $hashRef - reference to hash of hashses. Top level key is GO aspect, second level key is term and value is hypergeometric calculation
#          $correctedHashRef - same as hashRef with correction according to CORRECTION_TYPE and alpha
#          $numHypothesesTested - ref. to hash. Top level is GO aspect, second level is the number of hypotheses tested for that aspect
# Outputs: $hashRef - reference to hash of hashses. Top level key is domain aspect, second level key is term and value is hypergeometric calculation (top key is just to keep it in same form as other enrichment hashes
#          $correctedHashRef - reference to hash of hashes - same as hashRef by corrected for MHC
#         $numHypothesesTested - ref. to hash with top key equal domain aspect and value equal to number of hypotheses tested for that aspect.
#         $fgndHash - returns the fgnd hash with all terms and counts for those terms
# Kristen Naegle
# July 1, 2008
# Mod. Sept. 22, 2008 - to include count changes due to non-labeled proteins
# Mod. Oct. 18, 2008 - to include minCount
# Mod Oct. 22, 2008 - to allow for threshold in pval calculations and uniquify. Also moved code to another independent function for calculation 
# Mod March 18, 2009 - to include MHC Correction
sub calculateGOEnrichment($$$$$$$$){
    my ($dbh, $fgnd, $bgnd, $countUnlabeled, $minCount, $threshold, $CORRECTION_TYPE, $alpha) = @_;
# #     if($uniquify){
# # 	handleError('calculateGOEnrichment', 'TURNING OFF UNIQUIFY - not currently supported', \@_);
# #     }
    my $uniquify = 0;
    my $bgnd_GOHash = returnGOHashForMSIdArr($dbh, $bgnd, $uniquify);
    my $fgnd_GOHash = returnGOHashForMSIdArr($dbh, $fgnd, $uniquify);
    my ($proteinFgndArr, $n) = convertMSIdsToProteinIds($dbh, $fgnd, $uniquify);
    my ($proteinBgndArr, $N) = convertMSIdsToProteinIds($dbh, $bgnd, $uniquify);

    my ($enrichmentHash, $numHypothesesTested) = calculateEnrichmentFromHash($dbh, $n, $N, $fgnd_GOHash, $bgnd_GOHash, $minCount, $threshold, $countUnlabeled);
    
    my $correctedEnrichmentHash = pruneEnrichmentHash($enrichmentHash, $numHypothesesTested, $CORRECTION_TYPE, $alpha);
  
    return ($enrichmentHash, $correctedEnrichmentHash,  $numHypothesesTested, $fgnd_GOHash);
}


# ($hashRef, $correctedHashRef, $numHypothesesTested, $fgndHash) = calculateDomainEnrichment($dbh, $fgnd, $bgnd, $countUnlabeled, $minCount, $domainThreshold, $pvalThreshold, $uniquify);
# For a list of foreground protein Ids and a list of background protein Ids, report enrichment for terms occuring n times or more
# Returns a hash of hashes like the input foreground, except instead of count, there is a calculated probability value
# Inputs: $dbh - database handle
#         $fgnd - reference to array of MS ids of foreground
#         $bgnd - reference to array of MS ids of bacgkround
#         $countUnlabeled - flag: 1 - include unlabeled in calculations 0:exclude ~~~ and --- unlabeled 
#         $minCount - min. count in foreground to calculate pval
#         $domainThreshold - minimum p-value for including in count for domain of a protein
#         $pvalThreshold - min. p-value of hyge calculation required to report
#         $CORRECTION_TYPE - type of Multiple Hypotehsis testing correction
#         $alpha - target correction p-value
# REMOVED:        $uniquify - flag as to whether to uniquify protein ids from the MS id arrays before getting domain counts
# Outputs: $hashRef - reference to hash of hashses. Top level key is domain aspect, second level key is term and value is hypergeometric calculation (top key is just to keep it in same form as other enrichment hashes
#          $correctedHashRef - reference to hash of hashes - same as hashRef by corrected for MHC
#         $numHypothesesTested - ref. to hash with top key equal domain aspect and value equal to number of hypotheses tested for that aspect.
#         $fgndHash - returns the fgnd hash with all terms and counts for those terms
# Kristen Naegle
# July 1, 2008
# Mod. Oct. 18, 2008 to include minCount
# Mod Oct. 20, 2008 - to include uniquify option
# Mod Oct 22, 2008 - to use parallel code for enrichment calculations
# Mod March 18, 2009 - to include MHC Correction
sub calculateDomainEnrichment($$$$$$$$$){
    my ($dbh, $fgnd, $bgnd, $countUnlabeled, $minCount, $domainThreshold, $threshold, $CORRECTION_TYPE, $alpha) = @_;

# #     if($uniquify){
# # 	handleError('calculateDomainEnrichment', 'TURNING OFF UNIQUIFY - not currently supported', \@_);
# #     }
    my $uniquify = 0;
    
    my $bgnd_DomainHash = returnDomainHashForMSIdArr($dbh, $bgnd, $domainThreshold, $uniquify);
    my $fgnd_DomainHash = returnDomainHashForMSIdArr($dbh, $fgnd, $domainThreshold, $uniquify);
    my ($proteinFgndArr, $n) = convertMSIdsToProteinIds($dbh, $fgnd, $uniquify);
    my ($proteinBgndArr, $N) = convertMSIdsToProteinIds($dbh, $bgnd, $uniquify);

    my ($enrichmentHash, $numHypothesesTested) = calculateEnrichmentFromHash($dbh, $n, $N, $fgnd_DomainHash, $bgnd_DomainHash, $minCount, $threshold, $countUnlabeled);
    my $correctedEnrichmentHash = pruneEnrichmentHash($enrichmentHash, $numHypothesesTested, $CORRECTION_TYPE, $alpha);
    return ($enrichmentHash, $correctedEnrichmentHash, $numHypothesesTested, $fgnd_DomainHash);
}
# ($hashRef, $correctedHashRef, $numHypothesesTested, $fgndHash) = calculatePepPredEnrichment($dbh, $fgnd, $bgnd, $scansiteCutoff, $countUnlabeled, $minCount, $threshold, $CORRECTION_TYPE, $alpha);
# For a list of foreground phosphopep Ids and a list of background phosphopep Ids, report enrichment for terms occuring n times or more
# Returns a hash of hashes like the input foreground, except instead of count, there is a calculated probability value
# Inputs: $dbh - database handle
#         $fgnd - reference to array of phosphopep ids of foreground
#         $bgnd - reference to array of phosphopep ids of bacgkround
#         $scansiteCutoff - scansite cutoff value (5 is maximium - low stringency, 1 is med stringency, and 0.1 is high stringency)
#         $countUnlabeled - flag: 1 - include unlabeled in calculations 0:exclude ~~~ and --- unlabeled 
#         $minCount - minimum number of labels that must match to calculate enrichment
#         $threshold - pvalue threshold to consider for hyge report 
#         $CORRECTION_TYPE - type of Multiple Hypotehsis testing correction
#         $alpha - target correction p-value
# REMOVED:         $uniquify - flag 1: ignore redundant phosphopeps from MS assignment 0 means include them
# Outputs: $hashRef - reference to hash of hashses. Top level key is prediction type, second level key is term and value is hypergeometric calculation
#          $correctedHashRef - reference to hash of hashes - same as hashRef by corrected for MHC
#          $numHypothesesTested - ref. to hash with top level key equal to predition type and second level equal to the number of hypotheses tested
#         $fgndHash - returns the fgnd hash with all terms and counts for those terms
# Kristen Naegle
# July 1, 2008
# Modified Sept 22, 2008 - Now can vary count for fgnd and bgnd to ignore unlabeled data
# Modified October 15, 2008 - To include pelm_kinase predictions (version is hardcoded)
# Modified October 18, 2008 - To include a minCount term.
# Mod. October 20, 2008 - To include uniquify flag
# Mod. March 18, 2009 - To incorporate MHC. 
sub calculatePepPredEnrichment($$$$$$$$$){
    my ($dbh, $fgnd, $bgnd, $scansiteCutoff, $countUnlabeled, $minCount, $threshold, $CORRECTION_TYPE, $alpha) = @_;
# #     if($uniquify){
# # 	handleError('calculatePepPredEnrichment', 'TURNING OFF UNIQUIFY - not currently supported', \@_);
# #     }
    my $uniquify = 0;
    my ($phosphopepFgndArr, $n) = convertMSIdsToPhosphopepIds($dbh, $fgnd, $uniquify);
    my ($phosphopepBgndArr, $N) = convertMSIdsToPhosphopepIds($dbh, $bgnd, $uniquify);

    my $bgnd_PredHash = returnPepPredHashForMSIdArr($dbh, $bgnd, $scansiteCutoff, $uniquify);
    my $fgnd_PredHash = returnPepPredHashForMSIdArr($dbh, $fgnd, $scansiteCutoff, $uniquify);
    my ($enrichmentHash, $numHypothesesTested) = calculateEnrichmentFromHash($dbh, $n, $N, $fgnd_PredHash, $bgnd_PredHash, $minCount, $threshold, $countUnlabeled);

    #correct 
    my $correctedEnrichmentHash = pruneEnrichmentHash($enrichmentHash, $numHypothesesTested, $CORRECTION_TYPE, $alpha);
  
    return ($enrichmentHash, $correctedEnrichmentHash, $numHypothesesTested, $fgnd_PredHash);
   
}

# ($hashRef, $hashRefCorrected, $numHypothesesTested, $fgndHash) = calculateSiteDomainEnrichment($dbh, $fgnd, $bgnd, $countUnlabeled, $minCount, $threshold, $CORRECTION_TYPE, $alpha);
# For a list of foreground phosphopep Ids and a list of background phosphopep Ids, report enrichment for terms occuring n times or more
# Returns a hash of hashes like the input foreground, except instead of count, there is a calculated probability value
# Inputs: $dbh - database handle
#         $fgnd - reference to array of MS ids of foreground
#         $bgnd - reference to array of MS ids of bacgkround
#         $countUnlabeled - flag: 1 - include unlabeled in calculations 0:exclude ~~~ and --- unlabeled 
#         $minCount - min. count in fgnd required to calculate p-val
#         $threshold - min pvalue required to report significance of site Domain enrichment
#         $CORRECTION_TYPE - type of Multiple Hypotehsis testing correction
#         $alpha - target correction p-value
#  REMOVED:       $uniquify - flag 1: reduce redundancy at phosphopep level 0: don't uniquify
# Outputs: $hashRef - reference to hash of hashses. Top level key is prediction type, second level key is term and value is hypergeometric calculation
#          $correctedHashRef - reference to hash of hashes - same as hashRef by corrected for MHC
#         $numHypothesesTested - ref. to hash with top level key prediction type, and value equal to teh number of hypotheses tested for that 
#         $fgndHash - returns the fgnd hash with all terms and counts for those terms
# Kristen Naegle
# July 1, 2008
# Mod. Oct 20, 2008 - to switch to MS ids and apply uniquify option. Mod Oct 22 - to create parallel structure and add threshold.
# Mod March 11, 2009 - to count pfam_site only once for a peptide, unless a multiply phosphorylated peptide occurs in multiple domains - then handle as multiple sites.
# Mod. March 18, 2009 - To incorporate MHC. 
sub calculateSiteDomainEnrichment($$$$$$$$){
    my ($dbh, $fgnd, $bgnd, $countUnlabeled, $minCount, $threshold, $CORRECTION_TYPE, $alpha) = @_;
    my $uniquify = 0;
# #     if($uniquify){
# # 	handleError('calculateSiteDomainEnrichment', 'TURNING OFF UNIQUIFY - not currently supported', \@_);
# # 	$uniquify = 0;
# #     }

#     #my ($phosphopepFgndArr, $n) = convertMSIdsToPhosphopepIds($dbh, $fgnd, $uniquify);
#     #my ($phosphopepBgndArr, $N) = convertMSIdsToPhosphopepIds($dbh, $bgnd, $uniquify);
# #    print "DEBUG: In calculateSiteDomainEnrichment uniquify: $uniquify\n";
    my ($bgnd_siteDomainHash, $N) = returnSiteDomainHashForMSIdArr($dbh, $bgnd, $uniquify);
    my ($fgnd_siteDomainHash, $n) = returnSiteDomainHashForMSIdArr($dbh, $fgnd, $uniquify);
    my ($enrichmentHash, $numHypothesesTested) = calculateEnrichmentFromHash($dbh, $n, $N, $fgnd_siteDomainHash, $bgnd_siteDomainHash, $minCount, $threshold, $countUnlabeled);

    my $correctedEnrichmentHash = pruneEnrichmentHash($enrichmentHash, $numHypothesesTested, $CORRECTION_TYPE, $alpha);

    return ($enrichmentHash, $correctedEnrichmentHash, $numHypothesesTested, $fgnd_siteDomainHash);

}

# ($hashRef, $hashRefCorrected, $numTests, $fgndHash) = calculateMotifEnrichment($dbh, $fgnd, $bgnd, $threshold, $numAA, $CORRECTION_TYPE, $alpha);
# For a list of foreground phosphopep Ids and a list of background phosphopep Ids, report enrichment for motifs 
# Returns a hash of hashes like the input foreground, with secondary keys motifs and values equal to calculated p-value
# Inputs: $dbh - database handle
#         $fgnd - reference to array of MS ids of foreground
#         $bgnd - reference to array of MS ids of bacgkround
#         $threshold - motif threshold
#         $numAA - this is the number of amino acids on each side of the oriented phosphroylation site. $numAA must be between 1 and 7. 
#         $numPeptFraction - the fraction of the foreground that should match the parent motif before checking submotifs.
#         $CORRECTION_TYPE - type of Multiple Hypotehsis testing correction
#         $alpha - target correction p-value
#  NO Not really..but may need to       $conversionFile - the conversion file in order to process the motifs.
# REMOVED        $unquify - flag: if 0 then return all phosphopeps (redundant), else return only non-redundant phosphopep->aligned_peps
# Outputs: $hashRef - reference to hash of hashses. Top level key is prediction type, second level key is term and value is hypergeometric calculation
#          $correctedHashRefProcessed - reference to hash of hashes - same as hashRef by corrected for MHC, also processed post correction
#          $numTests - ref. to hash. Top level key is 'motifs' and second level is the count of hte log file .. wrong, but matches rest of structure
#         $fgndHash - returns the fgnd hash with all terms and counts for those terms (not returning the processed version, but could)
# Kristen Naegle
# August 15, 2008
# Mod. October 20, 2008 - to include flag for unqiue
# Mod. Feb 5, 2009 - to parse new motif output for number of tests run during enrichment
# Mod. March 18, 2009 - To incorporate MHC. NOTICE: Currently assumes conversionFile is found in local directory at 'conversion' and automatically processes motifs after enrichment.
sub calculateMotifEnrichment($$$$$$$$){
    my ($dbh, $fgnd, $bgnd, $cutoff, $numAA, $numpeptFraction, $CORRECTION_TYPE, $alpha) = @_;
    
# #     if($numAA < 1 or $numAA > 7){
# # 	print "ERROR: NUM_AA argument to calculateMotifEnrichment must be between 1 and 7 amino acids\n";
# # 	exit;

# #     }
    my $conversionFile = 'conversion';
    if($numAA != 5 and $numAA != 7){
	print "ERROR: NUM_AA argument to calculateMotifEnrichment must be either 5 or 7\n";
	exit;

    }

    my $uniquify = 0;

    my $bgnd_peps = returnPepAlignedForMSIdArr($dbh, $bgnd, $uniquify);
    my $fgnd_peps = returnPepAlignedForMSIdArr($dbh, $fgnd, $uniquify);
    # write both to a fgnd and bgnd file 
    my $fgndFile = "fgnd";
    my $bgndFile = "bgnd";
    
    if($numAA != 7){ #7 aa is the default return
	$fgnd_peps = cutAlignedPeptides($fgnd_peps, $numAA);
	$bgnd_peps = cutAlignedPeptides($bgnd_peps, $numAA);
    }

    printPepsToFile($fgndFile, $fgnd_peps);
    printPepsToFile($bgndFile, $bgnd_peps);
  #  my $numpept = int(scalar(@$fgnd_peps)/8);
    my $numpept = int(scalar(@$fgnd_peps)/$numpeptFraction);
    if($numpept < 2){
	$numpept = 2;
    }
    # run enrichment 
    my $logFile = "log.$fgndFile.$bgndFile";
    system("rm $logFile");
    #system("rm $logFile");

    my $tempFile = "TEMP_MOTIF_OUT";
    
    if($numAA ==7){
#	system("/programs/bin/c_stat_sig_peptide_refpass_depth_memory.pl --foreground $fgndFile --background $bgndFile --cutoff $cutoff --numpept $numpept > $tempFile");
	system("/programs/bin/c_stat_sig_peptide_refpass_depth_memory.phospho_acetyl.pl --foreground $fgndFile --background $bgndFile --cutoff $cutoff --numpept $numpept > $tempFile");
	
    }
    else{
	system("/programs/bin/c_stat_sig_peptide_refpass_depth_memory_5.pl --foreground $fgndFile --background $bgndFile --cutoff $cutoff --numpept $numpept > $tempFile");
    }
    # clean up and then parse file
    my ($pValHash, $fgndHash) = parseMotifLog($logFile);
    # this is wrong, but I want to keep in line with everything else
    my $numTests = parseMotifCountFromFile($tempFile);
    my %numTests;
    $numTests{'motifs'} = $numTests;
    system("rm $tempFile");

# #     print "DEBUG: IN Motif ENrichment for example looking at numTested\n";
# #     foreach my $a (keys %numTests){
# # 	print "aspect: $a\t $numTests{$a}\n";
# #     }

    my $correctedEnrichmentHash = pruneEnrichmentHash($pValHash, \%numTests, $CORRECTION_TYPE, $alpha);


    # In corrected EnrichmentHash process motifs
    my ($correctedEnrichmentHashProcessed, $fgndHashProcessed) = processMotifs($correctedEnrichmentHash, $fgndHash, $conversionFile);
    

    return ($pValHash, $correctedEnrichmentHashProcessed, \%numTests, $fgndHash);
}


# process motifs - remove those motifs from the pvalue hash that are less detailed descriptors of other motifs that belong to the exact same peptides
# Inputs: $pValHash - ref to p-value hash;
#         $countHash - count hash same keys as $p which are motifs
#         $cFile - conversion file
# Outputs: $processedP - pvalue hash that's been processed
#          $proccessedN - count hash that's been processed
# Brian Joughin
# March 18, 2009
sub processMotifs($$$) {
    # Function assumes all motifs are the same length!
    my($pValHash, $countHash, $cFile) = @_;
    my $c = readConversion($cFile);
#     #my $p = shift; # P-values
#     #my $n = shift; # Number of occurences
#     #my $c = shift; # Conversion table
    my $aspect = 'motifs';
    my $p = $pValHash->{$aspect};
    my $n = $countHash->{$aspect};

    my %newp;  my %newn; # Outputs.
    my %motifs_to_remove; # Stuff that should go away.

    my @motifs = keys %{$p}; # Array of motifs.

    foreach my $i (0..$#motifs) {   # Foreach motif i
	my @motif_i = split('', $motifs[$i]);

	foreach my $j (0..$#motifs) {
	    if ($i == $j) { next; } # Don't compare a motif to itself.
	    if ($n->{$motifs[$i]}->[0] != $n->{$motifs[$j]}->[0]) { next; } # Don't bother if the membership size is unequal.
	    
	    my @motif_j = split("", $motifs[$j]); # Foreach motif j with the same number of members as i...


	    my $match = 1; # Assume that i will be eliminated by j.
	    foreach my $aa (0..$#motif_i) { # At each amino acid position...
		# j might eliminate i if i is a ., the same as j, or a conversion-generalization of j
		if (!(($motif_i[$aa] eq $motif_j[$aa]) || ($motif_i[$aa] eq ".") || (defined $c->{$motif_i[$aa]}->{$motif_j[$aa]}))) {
		    $match = 0; # j does NOT eliminate i!
		    last; # So stop checking.
		}
	    }
	    if ($match == 1) { # if j eliminates i.
		#print STDERR "$motifs[$i] eliminated by $motifs[$j]\n";
		$motifs_to_remove{$motifs[$i]} = 1; # And remember to remove i.
		last; # Move to the next i.
	    }
	}
    }

    # Create %newc and %newp as uneliminated %c and %p.
    foreach my $motif (@motifs) {
	if (! defined $motifs_to_remove{$motif}) {
	    $newp{$motif} = $p->{$motif};
	    my @arr = []; 
	    $newn{$motif} = \@arr;    # I kept the array in the numpept hash.
	   # print "DEBUG: old count getting put in $n->{$motif}->[0]\n";
	    $newn{$motif}->[0] = $n->{$motif}->[0];
	    # print "$motif $p->{$motif} $n->{$motif}->[0]\n";
	}
    }

    my %processed; 
    $processed{$aspect} = \%newp;
    my %processedCount;
    $processedCount{$aspect} = \%newn;
	
    return (\%processed, \%processedCount); # Return results.
}

# \%multiHash = readConversion($conversionFile)
# Read a conversion file so that you can process motifs
# Inputs: $conversionFile - conversion file
# Outputs: %multiHash - hash of multis
# Brian Joughin
# March 18, 2009
sub readConversion($) {
    my ($c) = @_;
    my %multis;
    open (IN, $c) or die "Can't read conversion file $c";
    while(defined(my $line = <IN>)){
	if ($line =~ /\[/){
	    my @line = split(" ", $line);
	    $multis{$line[0]} = ();
	    my @multi = split("", $line[1]);
	    foreach my $index (1..$#multi-1){
		$multis{$line[0]}->{$multi[$index]} = 1;
	    }
	}
    }
    close IN;
    return \%multis;
}



# $newPeps = cutAlignedPeptides($peps, $numAA)
# This function takes an array of peptides and cuts them from 15-mers to numAA*2-1 mers
# This will simply re-return the same array of numAA is not valid
# Inputs: $peps - an array of peptides, assumes all are 15mers and already orietned (garbage in, garbage out)
#         $numAA - number of amino acids on each side of phosphorylation site that should be saved 
# Outptus: $newPeps - ref to array of peptides that have been cut
# Kristen Naegle
# Feb 21, 2009
sub cutAlignedPeptides($$){
    my($peps, $numAA) = @_;

    if($numAA >= 7 or $numAA < 1){
	return $peps;
    }
   
    my $diff = 7 - $numAA; # this is the number of chars to cut at each end
    
    my @newPeps;
    foreach my $p (@$peps){
	my $length = $numAA*2+1;
	my $new = substr($p, $diff, $length);
	    
	push @newPeps, $new;
    }
    
    return \@newPeps;
}


# $numTests = parseMotifCountFromFile($motifFile);
# From the piped output of a motif run, find the NUM_TESTS_FOR_KRISTEN field and return
# Inptus: $motifFile - this is /c_stat_sig_peptides.. call > motif_file
# Outputs: $numTests - the number of tests run in the motif enrichment
# Kristen Naegle
# Feb. 5, 2009
sub parseMotifCountFromFile($){
    my ($motifFile) = @_;
    
    my $motifFile2 = "TEMP_MOTIF_LINES";
    system("tail -n11 $motifFile > $motifFile2");

    open(TM, "$motifFile2") || die "Can't open $motifFile for reading\n";
    my @results = <TM>;
    chomp @results;
    close(TM);
    
    my $results = join " ", @results;
     my $numTests=0;
    if($results =~ m/NUM_TESTS\D+(\d+) \D/){
	$numTests = $1;
    }
   

    system("rm $motifFile2");
    return $numTests;

}

# (\%pValHash, \%fgndHash) = parseMotifLog($logFile)
# parses a motif log file into a hash
# Inputs: $logFile - log file of motif results (this processes it before parsing)
# Outputs: $pvalHash - ref. to hash with top key 'motifs' and value ref. to hashes with keys motif and value equal to pvalue
#          $fgndHash - ref. to foreground hash
# Kristen Naegle
# July, 2008
# Changed March 11, 2009 - now doesn't run on processed file, since we have FDR issues with processiong
# Fixed March 23, 2009 - wasn't able to parse count n from file properly if digits were too big. 
sub parseMotifLog($){
    my ($logFile) = @_;
    my $newLog = $logFile;
    # clean up the file
 #   my $newLog; # = $logFile.".processed"
    my %pValHash;
    my %fgndHash;
    #system("/programs/bin/process_motifs.pl $logFile > $newLog");
    open(LOG, $newLog) || die "Can't open $newLog for parsing\n";
    my %hash;
    my %fHash;
    while(defined(my $line = <LOG>)){
	$line =~ s/ //g;
#	my @line = split('\|', $line);
	my ($blank, $motif, $fgndStuff, $bgndStuff, $pval) = split('\|', $line);
	#print "DEBUG: motif:$motif\tfgndStuff:$fgndStuff\tbgnd$bgndStuff\tpval$pval\n";
	my @n = split('/', $fgndStuff);
	my $n = $n[0];
	$n=~ s/ //g;
	#print "$motif\t$pval\n";
	$hash{$motif} = $pval;
	push @{$fHash{$motif}}, $n;
	#print "num: $n\n";
	
	
    }
    $pValHash{'motifs'} = \%hash;
    $fgndHash{'motifs'} = \%fHash;
    close(LOG);
    return (\%pValHash, \%fgndHash);
}

# printPepsToFile($outputFile, $peps)
# Prints list of peptides to a file (i.e. create a fgnd, bgnd file)
# Input: $outputFile - name of file to write to 
#        $peps - ref to array of peptides
# Kristen Naegle
# August 15, 2008
sub printPepsToFile($$){
    my ($outputFile, $peps) = @_;
    open(FGND, ">$outputFile") || die "Can't open $outputFile for writing\n";
    foreach my $pep (@$peps){
	print FGND "$pep\n";
    }
    close(FGND);

}

# $hashEnrichment = calculateAllEnrichments($dbh, $msIds, $fgndMS, $countUnlabeled, $threshold, $motifThreshold, $numAA, $uniquify, $fractionOfForeground);
# calculates and returns all enrichment types (available to date)
# Inputs: $dbh - database handle
#         $msIds - All MS ids in background (ref to array)
#         $fgndMS - MS id array of foreground
#         $countLabel - flag to include count for all protein/peptides, or only those that are labeled. Currently only effects PhosphoPepPrediction and GOEnrichment 
#                     - if zero, include all, if 1 exclude just '~~~' and if 2 exclude '~~~' and '---' versions of unlabeled
#         $threshold - the threshold used in all runs except for motif. Motif threshold is fixed at 0.01
#         $numAA - the number of amino acids to run in motif enrichment (5 or 7 is allowed now).
#         $CORRECTION_TYPE - 'BF', 'BH', or 'NONE'
#         $fractionOfForeground - the number to divide the size of the set against (for e.g. if you want it to default to numCount=2, then this is N, if half the labels must hit than it's n/2, $numPeptFraction =2;
#         $scansiteThreshold - threshold for considering scansite cutoffs
#         $motifThreshold - the parent threshold required to dive into nodes
#         $domainThreshold - the threshold for considering a domain in enrichment
#   REMOVED      $uniquify - flag for whether to uniquify phosphopep and proteins before getting their label values _ MOD I TURN THIS FLAG ALWAYS OFF. (write error to log as well)
# Outputs: $hashEnrichment - hash of hash of hashes. Top key is description of type ('GO', 'domain', 'peptide_prediction', 'site_domain') and next hash has descriptive keys then labels and finally pvalue calcuations of those labels
#          $correctedHashEnrichment - pruned hash of enrichment based on alpha and correction type chosen
#          $numTests - has same structure as hashEnrichment, but contains the number of tests run per category
#          $fgnd - same strucutre, but contains labels 
# Kristen Naegle
# July 2, 2008
# Mod Sept. 22, 2008 - to allow for countLabel
# Mod Oct. 20, 2008 - To fix selections on MS ids only (to allow for redundancy)..Also added individual threshold and UNIQUIFY flag at this time.  Also, added check at this point to make sure that all foreground MS come from background (at not at a greater rate.
# Mod March 16, 2009 - to move correction into here (from print) and to call motif pruning since motif enrichment has changed to NOT include processing, also made scansite threshold, motifthreshold and numpeptFraction required.  Removing uniquify..setting to zero for all subfunction calls
# Mod April 1, 2009 to use calculateSingleEnrichment
# Mod July 30, 2009 to include all arguments
sub calculateAllEnrichments($$$$$$$$$$$){
    my ($dbh, $msIds, $fgndMS, $countUnlabeled, $threshold, $numAA, $CORRECTION_TYPE, $fraction, $scansiteThreshold, $motifThreshold, $domainThreshold) = @_;

    my $uniquify = 0;
    
    
  #  my $scansiteThreshold = 3;
  #  my $motifThreshold = 0.01; #eventually add this..but for now separate motif branching from threshold reporting
#    my $numPeptFraction = 8; #require 1/8th of data matches #passed in now
    #my $fraction = scalar(@$fgndMS)/$numPeptFraction; 

#    my $minCount = 2; # for now..all except motifs will require at least two
   # my $domainThreshold = 1e-6;

    my $CHECK_sample = checkFgndAgainstBgnd($fgndMS, $msIds);
    if($CHECK_sample){
	handleError('calculateAllEnrichments', "Your foreground was not sampled from your background\n", \@_);
	print "ERROR: Your foreground was not sampled from your background\n";
	print "\tCheck your foreground selection\n\n\n";
	exit;

    }
    my $CHECK = 0;
    my (%hash, %hashCorrected, %numTests, %fgnd);
    my($hash, $hashCorrected, $numTests, $fgnd); # the thing set for each enrichment type
    my @args;
    #GO enrichment
    @args = ($fraction);
    ($hash, $hashCorrected, $numTests, $fgnd) = calculateSingleEnrichment($dbh, $CHECK, 'GO', $msIds, $fgndMS, $countUnlabeled, $threshold, $CORRECTION_TYPE, @args);
    foreach my $key (keys %$hash){
	$hash{$key} = $hash->{$key};
	$hashCorrected{$key} = $hashCorrected->{$key};
	$numTests{$key} = $numTests->{$key};
	$fgnd{$key} = $fgnd->{$key};

    }
    
# motif enrichment
#    @args = ($motifThreshold, $numAA, $numPeptFraction);
    @args = ($motifThreshold, $numAA, $fraction);
    ($hash, $hashCorrected, $numTests, $fgnd) = calculateSingleEnrichment($dbh, $CHECK, 'sequence', $msIds, $fgndMS, $countUnlabeled, $threshold, $CORRECTION_TYPE, @args);
    foreach my $key (keys %$hash){
	$hash{$key} = $hash->{$key};
	$hashCorrected{$key} = $hashCorrected->{$key};
	$numTests{$key} = $numTests->{$key};
	$fgnd{$key} = $fgnd->{$key};

    }

    # peptide_prediction enrichment
    @args = ($scansiteThreshold, $fraction);
    ($hash, $hashCorrected, $numTests, $fgnd) = calculateSingleEnrichment($dbh, $CHECK, 'peptide_prediction', $msIds, $fgndMS, $countUnlabeled, $threshold, $CORRECTION_TYPE, @args);
    foreach my $key (keys %$hash){
	$hash{$key} = $hash->{$key};
	$hashCorrected{$key} = $hashCorrected->{$key};
	$numTests{$key} = $numTests->{$key};
	$fgnd{$key} = $fgnd->{$key};

    }

 # domain  enrichment
    @args = ($domainThreshold, $fraction);
    ($hash, $hashCorrected, $numTests, $fgnd) = calculateSingleEnrichment($dbh, $CHECK, 'domain', $msIds, $fgndMS, $countUnlabeled, $threshold, $CORRECTION_TYPE, @args);
    foreach my $key (keys %$hash){
	$hash{$key} = $hash->{$key};
	$hashCorrected{$key} = $hashCorrected->{$key};
	$numTests{$key} = $numTests->{$key};
	$fgnd{$key} = $fgnd->{$key};

    }

     # site_domain  enrichment
    @args = ($fraction);
    ($hash, $hashCorrected, $numTests, $fgnd) = calculateSingleEnrichment($dbh, $CHECK, 'site_domain', $msIds, $fgndMS, $countUnlabeled, $threshold, $CORRECTION_TYPE, @args);
    foreach my $key (keys %$hash){
	$hash{$key} = $hash->{$key};
	$hashCorrected{$key} = $hashCorrected->{$key};
	$numTests{$key} = $numTests->{$key};
	$fgnd{$key} = $fgnd->{$key};

    }

    #data feature enrichment
     @args = ($fraction);
    ($hash, $hashCorrected, $numTests, $fgnd) = calculateSingleEnrichment($dbh, $CHECK, 'data_feature', $msIds, $fgndMS, $countUnlabeled, $threshold, $CORRECTION_TYPE, @args);
    foreach my $key (keys %$hash){
	$hash{$key} = $hash->{$key};
	$hashCorrected{$key} = $hashCorrected->{$key};
	$numTests{$key} = $numTests->{$key};
	$fgnd{$key} = $fgnd->{$key};
	
    }
    

    return (\%hash, \%hashCorrected, \%numTests, \%fgnd);
}


# ($hashEnrichment, $correctedHashEnrichment, $numTests, $fgnd) = calculateSingleEnrichment($dbh, $CHECK, $enrichmentType, $msIds, $fgndMS, $countUnlabeled, $threshold, $CORRECTION_TYPE, @args);
# calculates and returns all enrichment types (available to date)
# Inputs: $dbh - database handle
#         $CHECK - if one, then check that the foreground comes from background. Put this bool value in so that if you call multiple times with the same fgnd (e.g. calculateAllEnrichments) you only have to perform this check before calling). 
#         $enrichmentType - type of enrichment to perform (options: GO, peptide_prediction, sequence, site_domain, domain)
#         $msIds - All MS ids in background (ref to array)
#         $fgndMS - MS id array of foreground
#         $countUnLabeled - flag to include count for all protein/peptides, or only those that are labeled. Currently only effects PhosphoPepPrediction and GOEnrichment 
#                     - if zero, include all, if 1 exclude just '~~~' and if 2 exclude '~~~' and '---' versions of unlabeled
#         $threshold - the corrected alpha threshold
#         $CORRECTION_TYPE - 'BF', 'BH', or 'NONE'
#         @args - required variable arguments dependent on enrichmentType
#                 peptide_prediction (scansiteThreshold, fraction)
#                 sequence (motifThreshold, numAA, fraction)
#                 GO (fraction)
#                 site_domain(fraction)
#                 domain (domainThreshold, fraction)
# Outputs: $hashEnrichment - hash of hash of hashes. Top key is description of type ('GO', 'domain', 'peptide_prediction', 'site_domain') and next hash has descriptive keys then labels and finally pvalue calcuations of those labels
#          $correctedHashEnrichment - pruned hash of enrichment based on alpha and correction type chosen
#          $numTests - has same structure as hashEnrichment, but contains the number of tests run per category
#          $fgnd - same strucutre, but contains labels 
# Kristen Naegle
# April 1, 2009
sub calculateSingleEnrichment{
    my ($dbh, $CHECK, $enrichmentType, $msIds, $fgndMS, $countLabel, $threshold, $CORRECTION_TYPE, @args) = @_;
    # check args
    my $alpha = $threshold;
    my $uniquify = 0;
    print "DEBUG: in Calculate Alpha=$threshold\n";
    my(%hash, %hashCorrected, %numTests, %fgnd);

    if($CHECK){
	my $CHECK_sample = checkFgndAgainstBgnd($fgndMS, $msIds);
	if($CHECK_sample){
	    handleError('calculateAllEnrichments', "Your foreground was not sampled from your background\n", \@_);
	    print "ERROR: Your foreground was not sampled from your background\n";
	    print "\tCheck your foreground selection\n\n\n";
	    exit;
	    
	}
    }
    
    if($enrichmentType eq 'sequence'){
	if(scalar(@args) != 3){
	    print "ARGS for sequence should be size 2 with values (motifTheshold, numAA, fraction)\n";
	    exit;
	}
	my ($motifThreshold, $numAA, $fraction) = @args;
	my ($pvalHashMotifs, $pvalHashMotifsCorrected, $numMotifTests, $fgndMotif) =  calculateMotifEnrichment($dbh, $fgndMS, $msIds, $motifThreshold, $numAA, $fraction, $CORRECTION_TYPE, $alpha);
	$hash{'sequence'} = $pvalHashMotifs;
	$hashCorrected{'sequence'} = $pvalHashMotifsCorrected;
	$numTests{'sequence'} = $numMotifTests;
	$fgnd{'sequence'} = $fgndMotif;
	
    }
    elsif($enrichmentType eq 'peptide_prediction'){
	# args: scansiteThreshold, minCount
	if(scalar(@args) != 2){
	    print "ARGS for peptide_prediciton should be size 2 with values (scansiteThreshold, fraction)\n";
	    exit;
	}
	my ($scansiteThreshold, $fraction) = @args;
	my $minCount = returnMinCount($fgndMS, $fraction);
	my ($pvalHashPepPred, $pvalHashPepPredCorrected, $numPepPredTests, $fgndPred) = calculatePepPredEnrichment($dbh, $fgndMS, $msIds, $scansiteThreshold, $countLabel, $minCount, $alpha, $CORRECTION_TYPE, $alpha);
	$hash{'peptide_prediction'} = $pvalHashPepPred;
	$hashCorrected{'peptide_prediction'} = $pvalHashPepPredCorrected;
	$numTests{'peptide_prediction'} = $numPepPredTests;
	$fgnd{'peptide_prediction'} = $fgndPred;
    }
    elsif($enrichmentType eq 'GO'){
	if(scalar(@args) != 1){
	    print "ARGS for GO should be size 1 with values (fraction)\n";
	    exit;
	}
	my ($fraction) = @args;
	my $minCount = returnMinCount($fgndMS, $fraction);
	    
	my ($pvalHashGO, $pvalHashGOCorrected, $numGOTests, $fgndGO) = calculateGOEnrichment($dbh, $fgndMS, $msIds, $countLabel, $minCount, $alpha, $CORRECTION_TYPE, $alpha);
	$hash{'GO'} = $pvalHashGO;
	$hashCorrected{'GO'} = $pvalHashGOCorrected;
	$numTests{'GO'} = $numGOTests;
	$fgnd{'GO'} = $fgndGO;
	
    }
    elsif($enrichmentType eq 'site_domain'){
	if(scalar(@args) != 1){
	    print "ARGS for site_domain should be size 1 with values (fraction)\n";
	}

	my ($fraction) = @args;
	my $minCount = returnMinCount($fgndMS, $fraction);
	my ($pvalHashSite, $pvalHashSiteCorrected, $numSiteDomainTests, $fgndSiteDomain) = calculateSiteDomainEnrichment($dbh, $fgndMS, $msIds, $countLabel, $minCount, $alpha, $CORRECTION_TYPE, $alpha);

	$hash{'site_domain'} = $pvalHashSite;
	$numTests{'site_domain'} = $numSiteDomainTests;
	$hashCorrected{'site_domain'} = $pvalHashSiteCorrected;
	$fgnd{'site_domain'} = $fgndSiteDomain;
	
    }
    elsif($enrichmentType eq 'domain'){
	if(scalar(@args) != 2){
	    print "ARGS for domain should be size 2 with values (domainThreshold, fraction)\n";
	    exit;
	}
	my ($domainThreshold, $fraction) = @args;
	my $minCount = returnMinCount($fgndMS, $fraction);
	my ($pvalHashDomain, $pvalHashDomainCorrected, $numDomainTests, $fgndDomain) = calculateDomainEnrichment($dbh, $fgndMS, $msIds, $countLabel, $minCount, $domainThreshold, $alpha, $CORRECTION_TYPE, $alpha);
	$hash{'domain'} = $pvalHashDomain;
	$hashCorrected{'domain'} = $pvalHashDomainCorrected;
	$numTests{'domain'} = $numDomainTests;
	$fgnd{'domain'} = $fgndDomain;

    }
    elsif($enrichmentType eq 'data_feature'){
	if(scalar(@args) != 1){
	    print "ARGS for data_feature should be size 1 with values (fraction)\n";
	}

	my ($fraction) = @args;
	my $minCount = returnMinCount($fgndMS, $fraction);
	my ($pvalHashData, $pvalHashDataCorrected, $numDataTests, $fgndData) = calculateDataFeatureEnrichment($dbh, $fgndMS, $msIds, $countLabel, $minCount, $alpha, $CORRECTION_TYPE, $alpha);
	$hash{$enrichmentType} = $pvalHashData;
	$hashCorrected{$enrichmentType} = $pvalHashDataCorrected;
	$numTests{$enrichmentType} = $numDataTests;
	$fgnd{$enrichmentType} = $fgndData;

    }


    else{

	print "ERROR\n";
	print "Expect enrichment_type to have value: sequence, domain, site_domain, peptide_prediction, or GO\n";
    }


    return (\%hash, \%hashCorrected, \%numTests, \%fgnd);

}

# # # (\%enrichmentHash, \@proteinSizeArrUnique, \@peptideSizeArr, \@peptideSizeArrUnique) = calculateRandomEnrichment($dbh, $expId, $fgndSize, $numRepeats, $threshold, $countLabel, $numAA, $uniquify, $CORRECTION_TYPE);
# # # Returns an enrichment hash (just like calculateAllEnrichment), but now $hash->{$metric}->{$label} is an array of pvalues (lose the specific label of the test. E.g. $metric='domain' and $label='pfam'
# # # Inputs: $dbh - database handle
# # #         $expId - id of experiment from which to pull randomly
# # #         $fgndSize - size of MSid array to create (in the future may change this to require specific protein and phosphosite id sizes)
# # #         $numRepeats - number of random tests to repeat
# # #         $threshold - pval threshold by which to cut-off insertion into array. Same threshold used for motif enrichment as well
# # #         $countLabel - flag for determining labels to ignore.
# # #                     - if zero, include all, if 1 exclude just '~~~' and if 2 exclude '~~~' and '---' versions of unlabeled
# # #         $numAA - number of amino acids to use in motif enrichment (5 or 7)
# # #         $uniquify - see other enrichment - dead var right now
# # #         $CORRECTION_TYPE - type of correction to use
# # # Outputs: $enrichmentHash - ref. to enrichment has as explained above
# # #          $proteinSizeArr - ref. to array of protein sizes chosen in each foreground
# # #          $peptideSizeArr - ref. to array of peptide (phosphosite) sizes chosen in each foreground. 
# # #          $CORRECTION_TYPE - type of correction to hold (will use $threshold as alpha value). 
# # # Kristen Naegle
# # # October 20, 2008         
# # sub calculateRandomEnrichment($$$$$$$$$){
# #     my ($dbh, $expId, $fgndSize, $numRepeats, $threshold, $countLabel, $numAA, $uniquify, $CORRECTION_TYPE, $fraction) = @_;
    
# #     my $msIds = returnMSIdsForExpId($dbh, $expId);

# #     my %masterHash;
# #     my %labelHash;
# #     my @proteinSizeArrUnique;
# #     my @peptideSizeArr;
# #     my @peptideSizeArrUnique;
# #     my %numTestHash;

# #     my $count = 0;
# #     for my $i (1..$numRepeats){
# # 	my $fgndMS = createRandomForeground($msIds, $fgndSize); #generate random MS fgnd
# # 	# debug lines for problems with enrichment
# # 	printLineToLog("\n");
# # 	printLineToLog("Checking MS Ids @$fgndMS\n");
# # 	my ($hashEnrichment, $numTestHash, $fgndHash) = calculateAllEnrichments($dbh, $msIds, $fgndMS, $countLabel, $threshold, $numAA, $CORRECTION_TYPE, $fraction);    
	
# # 	#print "DEBUG: Before prune have the following pvalues\n";
# # 	debug_printEnrichmentHash($hashEnrichment);
# # 	$hashEnrichment = pruneEnrichmentHashTop($hashEnrichment, $numTestHash, $CORRECTION_TYPE, $threshold);
# # 	#print "\n\n: DEBUG After prune\n";
# # 	debug_printEnrichmentHash($hashEnrichment);


# # 	my ($fgnd_proteinIds, $fgnd_count_proteinUnique) = convertMSIdsToProteinIds($dbh, $fgndMS, 1);
# #     my ($fgnd_phosphopepIds, $fgnd_count_pep) = convertMSIdsToPhosphopepIds($dbh, $fgndMS, 0);  # This count changes though drastically with ignore unlabeled
# # 	my ($fgnd_phosphopepIdsUnique, $fgnd_count_pepUnique) = convertMSIdsToPhosphopepIds($dbh, $fgndMS, 1);  # This count changes though drastically with ignore unlabeled
# # 	push @proteinSizeArrUnique, $fgnd_count_proteinUnique;
# # 	push @peptideSizeArr, $fgnd_count_pep;
# # 	push @peptideSizeArrUnique, $fgnd_count_pepUnique;
# # 	#print "Round $i found $fgnd_count_protein proteins and $fgnd_count_pep peptides\n";
	
# # 	foreach my $type (keys %$hashEnrichment){
# # 	    my $lowHash = $hashEnrichment->{$type}; #lowHash has keys equal to the metric and value equal to the p-value. want to combine these into an array of p-values
# # 	    my $lowTestHash = $numTestHash->{$type};
# # 	    foreach my $metric (keys %$lowHash){
# # 		if(not defined($masterHash{$type}->{$metric})){
# # 		    $masterHash{$type}->{$metric} = [];
# # 		}
# # 		#the number of labels in the metric hash is the number of tests run
# # 		if(not defined($numTestHash{$type}->{$metric})){
# # 		    $numTestHash{$type}->{$metric} = [];
# # 		}
# # 		if(not defined($labelHash{$type}->{$metric})){
# # 		    my %newHash;
# # 		    $labelHash{$type}->{$metric} = \%newHash;
# # 		}
# # 		my $pValHash = $lowHash->{$metric};
# # #		push @{$numTestHash{$type}->{$metric}}, scalar(keys %$pValHash);
# # 		push @{$numTestHash{$type}->{$metric}}, $lowTestHash->{$metric}; #changed Nov, 2009!! I had not separated pvalue of return with that of motif. 
# # 		foreach my $label (keys %$pValHash){
# # 		    my $pval = $pValHash->{$label};
# # 		    #if($pval < $threshold){
# # 		    push @{$masterHash{$type}->{$metric}}, $pval;
# # 		    $count += 1;
# # 		    #}
# # 		    if(not defined(($labelHash{$type}->{$metric})->{$label})){
# # 			my %anotherHash;
# # 			$anotherHash{'count'} = 0;
# # 			$anotherHash{'pval'} = [];
# # 			($labelHash{$type}->{$metric})->{$label} = \%anotherHash;
# # 		#	print "Adding $label to hash\n";
# # 		   }

# # 		    (($labelHash{$type}->{$metric})->{$label})->{'count'} += 1;
# # 		    push @{(($labelHash{$type}->{$metric})->{$label})->{'pval'}}, $pval;
# # 		}
# # 	    }
	    
# # 	}
# #     }
# #     print "FOUND $count total enrichments\n";
# #     return (\%masterHash, \%numTestHash, \%labelHash, \@proteinSizeArrUnique, \@peptideSizeArr,  \@peptideSizeArrUnique);
# # }


# printEnrichmentReport($enrichmentHash, $numTestHash, $fgndNum, $outputFile)
# Given an enrichment hash, print a report - those things that were enriched below a cutoff
# Inputs: $enrichmentHash - hash of hashes, top key is description of term and second level are labels with pValues as values
#         $numTestHash - same structure as enrichmentHash, but has values equal to the number of hypotheses tested for each aspect term
#         $fgndNum - hash with same structure, but has the count (in an array element) of the number of hits of that label in the foreground
#    REMOVED     $alpha - the threshold below which to print a result corrected by indicated test
#    REMOVED     $CORRECTION_TYPE - current options 'NONE', Benjamini and Hoetchberg 'BH'
#         $outputFile - destination of report (appends!)
# Kristen Naegle
# July 2, 2008
# November 26, 2008 - updated to include print with adjustment by Correction_TYPE 
# Modified March 19, 2009 - Now assumes that enrichmentHash is the modified hash (or not if'n you don't want)..but doesnt prune at print time. 
sub printEnrichmentReport($$$$){
    my ($enrichmentHash, $numTestHash, $fgndNum, $outputFile) = @_;
    open(OUT_R, ">>$outputFile");
    foreach my $aspect (keys %$enrichmentHash){
	my $pvalHash = $enrichmentHash->{$aspect};
	my $fgndHash = $fgndNum->{$aspect};
#	my @pValArr = values %$pvalHash;
	my $m = $numTestHash->{$aspect};
# #	my $pAdj = returnAdjPValue($CORRECTION_TYPE, \@pValArr, $m, $alpha);
# #	print OUT_R sprintf("%s,%d,%5f\n", $aspect, $m, $pAdj);
	print OUT_R sprintf("%s,%d\n", $aspect, $m);
	# #my $s = sprintf("%30s\t%3s", "term", 'pvalue');
# #	print OUT_R $s."\n";
	foreach my $term (keys %$pvalHash){
	    my $pval = $pvalHash->{$term};
	    my $n = $fgndHash->{$term}->[0];
	   # if($pval <= $pAdj){
	    my $s = sprintf("%s,%d,%10f", $term,$n, $pval);
	    print OUT_R $s."\n";
	   # }
	}
	print OUT_R "\n";
    }
    close(OUT_R);


}


# $prunedEnrichmentHash = pruneEnrichmentHash($enrichmentHash, $numTestHash, $CORRECTION_TYPE, $alpha)
# Returns a pruned enrichment hash with those terms from enrichmentHash removed that did not fall below alpha corrected according to type specified
# Inputs: $enrichmentHash - enrichment hash object, top key is aspect, second key is the label.  
#         $numTestHash - hash with the same structure, but value represents number of hypotheses tested for an aspect
#         $CORRECTION_TYPE - correction for multiple hypothesis testing. Currently supported 'NONE', 'BH' (Benjamini Hochberg FDR).
#         $alpha - target alpha value
# Outputs: $prunedEnrichmentHash - enrichment hash with fields removed not falling below alpha standard
# Kristen Naegle
# November 26, 2008 - NOT TESTED
sub pruneEnrichmentHash($$$$){
    my ($enrichmentHash, $numTestHash, $CORRECTION_TYPE, $alpha) = @_;
    my %prunedHash;
    foreach my $aspect (keys %$enrichmentHash){
	my $pValHash = $enrichmentHash->{$aspect};
	#convert all terms into a pval array
	my @pValArr = values %$pValHash;
#	print "DEBUG: PVALS: @pValArr\n";
	my $m = $numTestHash->{$aspect};
	#print "DEBUG: Calling returnAdjPValue with $CORRECTION_TYPE\t $m \t $alpha\n";
	my $pAdj = returnAdjPValue($CORRECTION_TYPE, \@pValArr, $m, $alpha); 
#	print "DEBUG: Adjusted $alpha to $pAdj\n";
	my %newHash;
	foreach my $term (keys %$pValHash){
	    if($pValHash->{$term} <= $pAdj){
		#print "DEBUG: Deleting $term with pValue $pValHash->{$term}\n";
		$newHash{$term} = $pValHash->{$term};
	    }
	}
	$prunedHash{$aspect} = \%newHash;
    }
    return \%prunedHash;
}

# $prunedHash = pruneEnrichmentHashTop(($enrichmentHash, $numTestHash, $CORRECTION_TYPE, $alpha) 
# The top function that runs on an aspect and calls pruneEnrichmentHash (which handles at the term level)
# Inputs: $enrichmentHash - enrichment hash object, top key is aspect, second key is the label.  
#         $numTestHash - hash with the same structure, but value represents number of hypotheses tested for an aspect
#         $CORRECTION_TYPE - correction for multiple hypothesis testing. Currently supported 'NONE', 'BH' (Benjamini Hochberg FDR).
#         $alpha - target alpha value
# Outputs: $prunedEnrichmentHash - enrichment hash with fields removed not falling below alpha standard
# Kristen Naegle
# November 26, 2008 - NOT TESTED
sub pruneEnrichmentHashTop($$$$){
    my ($enrichmentHash, $numTestHash, $CORRECTION_TYPE, $alpha) = @_;
    my %hash; 
    foreach my $aspect (keys %$enrichmentHash){
	my $nHash = pruneEnrichmentHash($enrichmentHash->{$aspect}, $numTestHash->{$aspect}, $CORRECTION_TYPE, $alpha);
	$hash{$aspect} = $nHash;

    }
    return \%hash;
    
}


# $pvalue = hyper($n, $N, $k, $D);
# Compute the hypergeometric of k occurences in foreground n, given D occurences in total background N
# Inputs: $n - size of foreground
#         $N - size of background
#         $k - number of hits in foreground
#         $D - number of hits in background
# Outputs: $pvalue - hypergeometric sum of probability of having $k or greater hits in foreground
# Brian Joughin - copied from /programs/bin/hyper.pl on July 1, 2008
sub hyper($$$$){
    my ($n, $N, $k, $D) = @_;

    my $sum = 0;

    if ( ($k == 0) ) { return 1; }
    if ( ($n == $k) && ($N == $D) ) { return 1; }
    my $max = $n;
    if ($D < $n) { $max = $D; }
    foreach my $i ($k..$max){
	# $sum += combination($D, $i) * combination($N-$D, $n-$i) / combination($N, $n);
	my $log = 0;
	# combination($D, $i)
	foreach my $index (0..$i-1){
	    $log += (log10($D-$index) - log10($i-$index));
#		print "$i $index 1 " . (log10($D-$index) - log10($i-$index))  . "\n";
	}
	# combination($N-$D, $n-$i)
	foreach my $index (0..$n-$i-1) {
	    $log += (log10($N-$D-$index) - log10($n-$i-$index));
#		print "$i $index 2 " . (log10($N-$D-$index) - log10($n-$i-$index))  . "\n";
	}
	# -combination($N, $n)
	foreach my $index (0..$n-1){
	    $log -= (log10($N-$index) - log10($n-$index));
#		print "$i $index 3 " . (log10($N-$index) - log10($n-$index))  . "\n";
	}
	$sum += 10**$log;
    }
    
    return $sum;
}


#  my($n) = returnCount($ids, $ids_aspect, $LABELED);
#  Returns the count in an array - i.e. if $LABELED is zero, simply returns the size of the array of ids, otherwise it returns the size - those that are unlabeled 
# Inputs: $ids - array of ids (protein, or phosphopep ids)
#         $ids_aspect - a hash with keys values and values are an array with zero position equal to the count of that..this is used for the '~~~' and '---' values
#         $LABELED - if zero, include all, if 1 exclude just '~~~' and if 2 exclude '~~~' and '---' versions of unlabeled
#  Output: $n - the count
# Kristen Naegle
# Sept. 22, 2008
sub returnCount($$$){
    my($ids, $ids_aspect, $LABELED) = @_;
    my($n);
    $n = scalar(@$ids);

    
    if ($LABELED){
	my $unlabeled = 0;
	if(defined($ids_aspect->{'~~~'})){
	    $unlabeled = $ids_aspect->{'~~~'}->[0];
	}
	if($LABELED == 2){
	    if(defined($ids_aspect->{'---'})){
		$unlabeled += $ids_aspect->{'---'}->[0];
	    }
	    
	}
	$n = $n-$unlabeled;
    }
    
    return $n;

}


#  my($n) = returnCount($ids, $ids_aspect, $LABELED);
#  Returns the count in an array - i.e. if $LABELED is zero, simply returns the size of the array of ids, otherwise it returns the size - those that are unlabeled 
# Inputs: $ids - array of ids (protein, or phosphopep ids)
#         $ids_aspect - a hash with keys values and values are an array with zero position equal to the count of that..this is used for the '~~~' and '---' values
#         $LABELED - if zero, include all, if 1 exclude just '~~~' and if 2 exclude '~~~' and '---' versions of unlabeled
#  Output: $n - the count
# Kristen Naegle
# Sept. 22, 2008
sub returnHashCount($$$){
    my($ids, $ids_aspect, $LABELED) = @_;
    my($n);
    $n = scalar(@$ids);

    
    if ($LABELED){
	my $unlabeled = 0;
	if(defined($ids_aspect->{'~~~'})){
	    $unlabeled = $ids_aspect->{'~~~'}->[0];
	}
	if($LABELED == 2){
	    if(defined($ids_aspect->{'---'})){
		$unlabeled += $ids_aspect->{'---'}->[0];
	    }
	    
	}
	$n = $n-$unlabeled;
    }
    
    return $n;

}


# $error = checkFgndAgainstBgnd($fgndArr, $bgndArr);
# Checks to make sure that a fgnd was sampled from the background given (checks the number of appearances as well that count_fgnd <= count_bgnd
# Inputs: $fgndArr - ref. to array of fgnd objects
#         $bgndArr - ref. to array of bgnd objects
# Outputs: $error - 0 if passes collection test
#                   1 if it doesn't
# Kristen Naegle
# October 20, 2008
sub checkFgndAgainstBgnd($$){
    my ($fgnd, $bgnd) = @_;
    my $bgndHash = returnCountHash($bgnd);
    my $fgndHash = returnCountHash($fgnd);
    my $error = 0;
    my $fgndCount;
    my $bgndCount;
    foreach my $f (keys %$fgndHash){
	if(not defined $bgndHash->{$f} or ($bgndHash->{$f} < $fgndHash->{$f})){
	    $error = 1;
	}
    }

    return $error;

}

# $hashRef = returnCountHash($arr)
# Given an array, return a hash with key equal to objects, and value equal to count, number of times that object appears in array
# Inputs: $arr - ref. to array of objects
# Outputs: $hash - ref. to hash described above
# Kristen Naegle
# October 20, 2008
sub returnCountHash($){
    my ($arr) = @_;
    my %hash;

    foreach my $item (@$arr){
	if(not defined $hash{$item}){
	    $hash{$item} = 0;
	}
	$hash{$item} += 1;
    }
    return \%hash;
}

# \%countHash = convertResultsHashToCoverage($hash);
# Given a hash of results (i.e. top key is a value like protein_id and it points to an array of hashes of results, convert this to a count for number predictions/annotations. 
# Inputs: $hash - as described above
#         $field - the field you want to look at specifically (i.e. 'scansite_kinase' would come from 'source'
#         $value - 
#         $nullFieldCheck - this is empty if there is not value to check for '~~~' value to remove from count. Otherewise send the name of the column (e.g. value or label) that says where this should be checked. 
#         $optHashArg - this is an optional score/pvalue argument. Assumes you want field of hash according to $optHash->{'field'} to be <= $optHashArg->{'value'}
# Outputs: $countHash - same top level structure, but now 
sub convertResultHashToCoverage{
    my ($hash, $field, $value, $optHashArg) = @_;
    my %countHash;
    my $OPT = 0;
    my $optField = $optHashArg->{'field'};
    my $optValue = $optHashArg->{'value'};
    if($optField){
	$OPT = 1;
	
    }
    foreach my $obj (keys %$hash){
	$countHash{$obj} = 0;
	my $results = $hash->{$obj};
	foreach my $result (@$results){
	    my $ADD = 1;
	    
	    if($result->{$field} eq $value or ($field eq 'pfam_site' && $result->{$field} ne '-1')){
		if($OPT){
		    if($result->{$optField} > $optValue){
			$ADD = 0;
		    }
		    
		}
		$countHash{$obj} += $ADD;
	    }

	}
    }
    return \%countHash;
}


# $metricCountHash = returnCountLabelHashForAllMetrics($dbh, $msIds, $minCount, $domainThreshold, $scansiteThreshold, $countUnlabeled);
# Needed a hash of all the metrics that keeps track of COUNTABLE labels (i.e. those that occur with at least minCount.  
# Inputs: $dbh - database handle
#         $msIds - All MS ids in background (ref to array)
#         $minCount - minimum number that label must appear
#         $domainThreshold - the threshold used to suppress bad domain predictions
#         $scansiteThreshold - the threshold used to suppress poor scansite predictions
#         $countUnLabeled - flag to include count for all protein/peptides, or only those that are labeled. Currently only effects PhosphoPepPrediction and GOEnrichment. If one, count unlabeled as well
# Outputs: $hashEnrichment - hash of hash of hashes. Top key is description of type ('GO', 'domain', 'peptide_prediction', 'site_domain') and next hash has descriptive keys then labels and finally number of times that label occurs in the given MS Id set
# Kristen Naegle
# November 7, 2008
sub returnCountLabelHashForAllMetrics($$$$$$){
    my ($dbh, $msIds, $minCount, $domainThreshold, $scansiteCutoff, $countUnlabeled) = @_;
    
    my $uniquify = 0;
    
    my $GOHash = returnGOHashForMSIdArr($dbh, $msIds, $uniquify);
    my $domainHash = returnDomainHashForMSIdArr($dbh, $msIds, $domainThreshold, $uniquify);
    my $pepPredHash = returnPepPredHashForMSIdArr($dbh, $msIds, $scansiteCutoff, $uniquify);
    my $siteDomainHash = returnSiteDomainHashForMSIdArr($dbh, $msIds, $uniquify);

    my %hash;
    $hash{'GO'} = cutCountsUnderMinFromHash($GOHash, $minCount, $countUnlabeled);
    $hash{'domain'} = cutCountsUnderMinFromHash($domainHash, $minCount, $countUnlabeled);
    $hash{'peptide_prediction'} = cutCountsUnderMinFromHash($pepPredHash, $minCount, $countUnlabeled);
    $hash{'site_domain'} = cutCountsUnderMinFromHash($siteDomainHash, $minCount, $countUnlabeled);
    return \%hash;
}

# \%cutHash = cutCountsUnderMinFromHash($hash, $minCount)
# Given a count hash for labels (with top level key equal to metric, and subsequent keys labels, with values count), cut labels from hash that are below minimum count
# Inputs: $hash - ref. to hash as described above (like returned from returnGOHashForMSIdArr)
#         $minCount - minimum number of occurences for that label required to keep
#         $countUnlabeled - this removes --- and ~~~ from the allowed labeled metrics
# Outputs: \%cutHash - ref. to new hash with labels removed that were under count
# Kristen Naegle
# November 7, 2008
sub cutCountsUnderMinFromHash($$$){
    my ($hash, $minCount, $countUnlabeled) = @_;
    my %cutHash;
    
    foreach my $metric (keys %$hash){
	my $labelHash = $hash->{$metric};
	my %newlabelHash;
	foreach my $label (keys %$labelHash){
	    my $count = ($labelHash->{$label})->[0];
	    if($count >= $minCount){
		if(!$countUnlabeled){
		    if($label eq '~~~' or $label eq '---'){
			next;
		    }
		}
		$newlabelHash{$label} = $count;
	    }
	    
	}
	$cutHash{$metric} = \%newlabelHash;

    }
    return \%cutHash;

}


# \%hash = returnProteinMetricCoverage($dbh, $MSIds, $coverageHash, $minCount, $domainThreshold, $scansiteThreshold);
# Returns a hash with keys that match other description hashes (top level e.g. GO, next is metric: F, then keys of proteins with arrays of labels). 
# Inputs: $dbh - databasehandle
#         $MSIds - arr($error, $hashEnrichment, $hashCorrected, $numTests, $fgndHash) = returnEnrichmentHashForReport($dbh, $bgnd, $fgnd, $countUnlabeled, $alpha, $numAA, $CORRECTION_TYPE, $fraction, $enrichType, @enrichArgs);ay of MS ids to map to protein ids to get coverage
#         $coverageHash - this is the list of acceptable labels to count for the various metrics (see returnCountLabelHashForAllMetrics). This hash should already be cut down for both unlabeled and minCount.  
#         $domainThrehsold - threshold of domain predictions to consider
#         $scansiteThreshold - threshold for scansite predictions to consider
# Outputs:  
sub returnObjectMetricCoverage($$$$$){
    my ($dbh, $MSIds, $coverageHash, $domainThreshold, $scansiteThreshold) = @_;
    my %hash;
    #Get Protein level info first
    my $uniquify = 1;
    my ($proteinIdArr, $proteinCount) = convertMSIdsToProteinIds($dbh, $MSIds, $uniquify);
    my ($hashPfam, $domainsAllowed, $domainsCut) = returnDomainLabelsByProteinId($dbh, $proteinIdArr, $coverageHash, $domainThreshold);
    my ($GOHash, $GOAllowed, $GOCut) = returnGOLabelsByProteinId($dbh, $proteinIdArr, $coverageHash);
    


    # Site level stuff
    my ($pepIdArr, $pepCount) = convertMSIdsToPhosphopepIds($dbh, $MSIds, $uniquify);
    my ($pepPredHash, $predAllowed, $predCut) = returnPepPredLabelsByPepId($dbh, $pepIdArr, $coverageHash, $scansiteThreshold);
    my ($pfamSiteHash, $pfamSiteAllowed, $pfamSiteCut) = returnPfamSiteLabelsByPepId($dbh, $pepIdArr, $coverageHash);
    

    $hash{'domain'} = $hashPfam;
    $hash{'peptide_prediction'} = $pepPredHash;
    $hash{'site_domain'} = $pfamSiteHash;
    $hash{'GO'} = $GOHash;

    print "DOMAINS ALLOWED: @$domainsAllowed\n\n";
    print "DOMAINS CUT: @$domainsCut\n\n";
    print "PREDICTIONS ALLOWED: @$predAllowed\n\n";
    print "PREDICTIONS CUT: @$predCut\n\n";
    print "PFAM_SITES ALLOWED: @$pfamSiteAllowed\n\n";
    print "PFAM_SITES CUT: @$pfamSiteCut\n\n";
    print "GO ALLOWED: @$GOAllowed\n\n";
    print "GO CUT: @$GOCut\n\n";

    return \%hash;
}


# (\%countPredHash, $predictionsAllowed, $predictionsCut)=returnPepPredLabelsByPepId($dbh, $pepIdArr, $allowedLabelHash, $scansiteThreshold)
# Get array of phosphopep_prediction labels that match a set of phosphopep_ids
# Inputs: $dbh - database handle - I'm getting really tired of explaining this one..am I up to a million lines of code yet?
#         $pepIdArr - array of pep Ids for which to get a prediction label count
#         $allowedLabelHash - this is the full metric coverage, within here we cut down to pfam/domain
#         $scansiteThreshold - required threshold for prediction of scansite
# Outputs: $countPredHash - ref. to hash with keys equal to the peptide ids and values equal to the number of labels of predictions that were found (and allowed)
#          $predictionsAllowed - ref to array of predictions that were allowed
#          $predictionsCut - ref. to array of predictions that were cut
# Kristen Naegle
# November 7, 2008
sub returnPepPredLabelsByPepId($$$$){
    my ($dbh, $pepIdArr, $allowedLabelHash, $scansiteThreshold) = @_;
   
    my $table = "phosphopep_prediction";
    my $field = "phosphopep_id";
    my $expression = "=";
    my ($predictionHash) = returnMultipleResultsSingleTableSearch($dbh, $table, $field, $expression, $pepIdArr);

    #cutting the hash like this with all pep_predictions, assumes no labels are the same between the metrics..this is true surely for scansite_bind and scansite_kinase.  At the time it also checked out to be true for scansite_kinase and pelm_kinase overlap
    my @predictions = ('scansite_kinase', 'scansite_bind', 'pelm_kinase');
    my $dH = $allowedLabelHash->{'peptide_prediction'};
    

    # Create a master hash for allowed keys, since I'm comparing based on all predictions
    my %masterAllowed;
    
    foreach my $prediction (@predictions){
	my $allowedLabelHashPred = $dH->{$prediction};
	foreach my $key (keys %$allowedLabelHashPred){
	    $masterAllowed{$key} = $allowedLabelHashPred->{$key};

	}
    }
	
    ($predictionHash, my $predictionsAllowed, my $predictionsCut) = cutObjectMetricLabelsByAllowable($predictionHash, "value", \%masterAllowed);
    my %siteHash;
    #count for each prediction
    $field = 'source';
    my $value = 'scansite_kinase';
    my %optArgs;
    $optArgs{'field'} = "score";
    $optArgs{'value'} = $scansiteThreshold;
    my $countScansiteKinaseHash = convertResultHashToCoverage($predictionHash, $field, $value, \%optArgs);
    $siteHash{$value} = $countScansiteKinaseHash;

    $value = 'scansite_bind';
    my $countScansiteBindHash = convertResultHashToCoverage($predictionHash, $field, $value,\%optArgs);
    $siteHash{$value} = $countScansiteBindHash;

    $value = 'pelm_kinase';
    my $countPELMKinaseHash = convertResultHashToCoverage($predictionHash, $field, $value);
    $siteHash{$value} = $countPELMKinaseHash;
    
    return (\%siteHash, $predictionsAllowed, $predictionsCut);
}


# (\%pfamSitesHash, $pfamSitesAllowed, $pfamSitesCut)=returnPfamSiteLabelsByPepId($dbh, $pepIdArr, $allowedLabelHash)
# Get array of phosphopep_prediction labels that match a set of phosphopep_ids
# Inputs: $dbh - database handle - I'm getting really tired of explaining this one..am I up to a million lines of code yet?
#         $pepIdArr - array of pep Ids for which to get a prediction label count
#         $allowedLabelHash - this is the full metric coverage, within here we cut down to pfam/domain
# Outputs: $pfamSiteHash - ref. to hash with keys equal to the peptide ids and values equal to the number of labels of pfam sites found
#          $pfamSitesAllowed - ref to array of predictions that were allowed
#          $pfamSitesCut - ref. to array of predictions that were cut
# Kristen Naegle
# November 7, 2008
sub returnPfamSiteLabelsByPepId($$$$){
    my ($dbh, $pepIdArr, $allowedLabelHash) = @_;
   
    my $table = "phosphopep";
    my $field = "id";
    my $expression = "=";
    my ($pepHash) = returnMultipleResultsSingleTableSearch($dbh, $table, $field, $expression, $pepIdArr);

    my $dH = $allowedLabelHash->{'site_domain'};
    
    ($pepHash, my $pfamSitesAllowed, my $pfamSitesCut) = cutObjectMetricLabelsByAllowable($pepHash, "pfam_site", $dH->{'pfam_site'});
    
    my %siteHash;
    #count for each prediction
    $field = 'pfam_site';
    my $value = '';
    my $countPfamSiteHash = convertResultHashToCoverage($pepHash, $field, $value);
    $siteHash{$field} = $countPfamSiteHash;

    return (\%siteHash, $pfamSitesAllowed, $pfamSitesCut);
}


# $domainPfamHash = returnDomainLablesByProteinId($dbh, $proteinIdArr, $allowedLabelHash, $domainThreshold);
# Get array of domain labels that match a set of proteins (with min domain prediction)
# Inputs: $dbh - database handle - I'm getting really tired of explaining this one..am I up to a million lines of code yet?
#         $proteinIdArr - array of protein Ids for which to get a domain label count
#         $allowedLabelHash - this is the full metric coverage, within here we cut down to pfam/domain
#         $domainThreshold - required threshold for prediction of domain to count
# Outputs: $countPfamHash - ref. to hash with keys equal to the protein ids and values equal to the number of labels of domains that were found
#          $domainsAllowed - ref to array of domains that were allowed
#          $domainsCut - ref. to array of domains that were cut
# Kristen Naegle
# November 7, 2008
sub returnDomainLabelsByProteinId($$$$){
    my ($dbh, $proteinIdArr, $allowedLabelHash, $domainThreshold) = @_;
    my($table, $field, $expression);
    $table = "domain";
    $field = "protein_id";
    $expression = "=";
    my ($domainHash) = returnMultipleResultsSingleTableSearch($dbh, $table, $field, $expression, $proteinIdArr);

    # go through here and cut domainHash labels before counting 
    #field here is label;
   # my @ks = keys %$allowedLabelHash;
   # print "Keys in LH: @ks\n";
    my $dH = $allowedLabelHash->{'domain'};
    #my @ks = keys %$dH;
    #print "Keys in DH: @ks\n";
    my $allowedLabelHashDomain = $dH->{'pfam'};
    ($domainHash, my $domainsAllowed, my $domainsCut) = cutObjectMetricLabelsByAllowable($domainHash, "label", $allowedLabelHashDomain);
#convert to count for those of pfam decent
    $field = 'source';
    my $value = 'pfam';
    my %optArgs;
    $optArgs{'field'} = "p_value";
    $optArgs{'value'} = $domainThreshold;
    my $countPfamHash = convertResultHashToCoverage($domainHash, $field, $value, \%optArgs);
    my %hashWKey;
    $hashWKey{'pfam'} = $countPfamHash;
    return (\%hashWKey, $domainsAllowed, $domainsCut);

}

# $domainPfamHash = returnDomainLablesByProteinId($dbh, $proteinIdArr, $allowedLabelHash, $domainThreshold);
# Get array of domain labels that match a set of proteins (with min domain prediction)
# Inputs: $dbh - database handle - I'm getting really tired of explaining this one..am I up to a million lines of code yet?
#         $proteinIdArr - array of protein Ids for which to get a domain label count
#         $allowedLabelHash - this is the full metric coverage, within here we cut down to pfam/domain
# Outputs: $countGOHash - ref. to hash with keys equal to the protein ids and values equal to the number of labels of GO terms found
#          $GOTermsAllowed - ref to array of domains that were allowed
#          $GOTermsCut - ref. to array of domains that were cut
# Kristen Naegle
# November 7, 2008
sub returnGOLabelsByProteinId($$$$){
    my ($dbh, $proteinIdArr, $allowedLabelHash) = @_;

    my $sth = $dbh->prepare("SELECT protein_id, GO.* from protein_GO join GO on protein_GO.GO_id=GO.id where protein_id=?");
    my %GOHash;
    foreach my $protein (@$proteinIdArr){
	$sth->execute($protein);
	$GOHash{$protein} = returnMultipleResultsForExSTH($sth);
    }
    
    my $dH = $allowedLabelHash->{'GO'};
    my @GOs = ('F', 'C', 'P');

    my %masterAllowed;
    
    foreach my $GO (@GOs){
	my $allowedLabelHashPred = $dH->{$GO};
	foreach my $key (keys %$allowedLabelHashPred){
	    $masterAllowed{$key} = $allowedLabelHashPred->{$key};

	}
    }
    
    my ($GOHash, $GOsAllowed, $GOsCut) = cutObjectMetricLabelsByAllowable(\%GOHash, "term", \%masterAllowed);

    my %hashWKey;

    my $field = 'aspect';
    foreach my $value (@GOs){
	my $GOValHash = convertResultHashToCoverage($GOHash, $field, $value);
	$hashWKey{$value} = $GOValHash;
    }
    return (\%hashWKey, $GOsAllowed, $GOsCut);

}


# $cutHash = cutObjectMetricLabelsByAllowable($objHash, $field, $allowedHash);
# Want to force removal (by setting to -1) the fields of a result that don't occurin the allowedLabelHash
# Inputs: $objHash - this is something that is returned by returnMultipleResultsSingleTableSearch. It has keys equal to an objec, that points to an array of hashes of results (where that hash has keys equal to the column field name and value equal to column value)
#         $field - the field for which to look for allowance within $allowedHash (e.g. If this is a domain cut, then the field is 'label'
#         $allowedHash - this is a subset of totalCoverageHash, has keys equal to label, and definition of that label as a key is enough to suggest allowance
# Outputs: $cutHash - same as input, but any label that did not appear in allowedHash now has values equal to -1 so it will be dismissed during conversion to count
#         \@allowed - ref. to an array of labels that were allowed
#         \@cut - ref. to an array of labels that were cut
# Kristen Naegle
# November 7, 2008
sub cutObjectMetricLabelsByAllowable($$$){
    my ($objHash, $field, $allowedHash) = @_;
 #   print "DEBUG: ACCEPTABLE labels\n";
    my @allowed = keys %$allowedHash;
 #   print "@allowed\n";
    my %cut;
    foreach my $obj (keys %$objHash){ #for example protein is the object
	my $objArr = $objHash->{$obj};
	#my @newArr = [];  #turn all the fields to -1;
	foreach my $result (@$objArr){
	    my $label = $result->{$field};
	    if(not defined $allowedHash->{$label}){
	#	print "DEBUG: Cutting $label\n";
		if(not defined $cut{$label}){
		    $cut{$label} = 1;
		}
		foreach my $f (keys %$result){
		    
		    $result->{$f} = -1; # is this really changing it in memory?
		}
	    }

	}
    }
    my @cut = keys %cut;
    return ($objHash,\@allowed, \@cut);
}

# printMatlabCoverage($dbh, $hash, $outputFile)
# Print a matlab file with the coverage structure in order to analyze.  Prints protein and peptide structs separately 
# Inputs: $dbh -database handle
#         $hash - this is the metric coverage hash 
#         $outputFile - name of file to print structs to
# Kristen Naegle
# November 9, 2008
sub printMatlabCoverage($$$){
    my ($dbh, $hash, $outputFile) = @_;

    my %proteinMetricsTop; 
    @{$proteinMetricsTop{'domain'}} = ('pfam');
    @{$proteinMetricsTop{'GO'}} = ('F', 'C', 'P');
    
 
# first print protein stuff..convert protein ids to acc_gene
    my $tKey = $hash->{'domain'};
    my $protHashEx = $tKey->{'pfam'};
    my @proteinIds = keys %$protHashEx;
   # print "DEBUG: Found protein ids: @proteinIds\n";
    my %convProtNames;
    my $sth = $dbh->prepare('SELECT * from protein where id=?');
    foreach my $protein (@proteinIds){
	$sth->execute($protein);
	my $result = returnMultipleResultsForExSTH($sth);
	my $rHash = $result->[0];
	$convProtNames{$protein} = $rHash->{'acc_gene'};
#	print "$protein translates to $rHash->{'acc_gene'}\n";
    }

    #print protein structure to file;
    printStructToFile(\%proteinMetricsTop, $hash, \%convProtNames, 'protein', $outputFile);
    

    # Now take care of site info
    my %siteMetricsTop;
    @{$siteMetricsTop{'peptide_prediction'}} = ('scansite_kinase', 'scansite_bind', 'pelm_kinase');
    @{$siteMetricsTop{'site_domain'}} = ('pfam_site');
    my $MSKey = $hash->{'peptide_prediction'};
    my $MSHashEx = $MSKey->{'scansite_kinase'};
    my @pepIds = keys %$MSHashEx;
    my %convSiteNames;
    $sth = $dbh->prepare('SELECT phosphopep.*, protein.acc_gene from phosphopep join protein on phosphopep.protein_id=protein.id where phosphopep.id=?');
    foreach my $pep (@pepIds){
	$sth->execute($pep);
	my $result = returnMultipleResultsForExSTH($sth);
	my $name = ($result->[0])->{'acc_gene'};
	foreach my $r (@$result){
	    my $site_type = $r->{'site_type'};
	    my $site_pos = $r->{'site_pos'};
	    $name .= "_".$site_type.$site_pos;
	}
	$convSiteNames{$pep} = $name;
	#print "$pep translates to $name\n";

    }
    printStructToFile(\%siteMetricsTop, $hash, \%convSiteNames, 'site', $outputFile);
     


    
}

# printStructToFle($topHash, $dataHash, $convHash, $appendName, $outputFile);
# Takes a dataHash (countHash) and prints colNames, rowNames and a data matrix of the count to a file..occording specifically to the values passed in by tophash (so you can print various aspects of a dataHash).  Uses appendName to be able to affect sructure name in Matlab, uses an append to Outputfile.
# Inputs: $topHash - hash of hashes, whose keys (top) and then next layer indicate what part of dataHash to print
#         $dataHash - this is the actual count hash. Key levels match that of the topHahs
#         $convHash - this is a hash with keys equal to the object ids in dataHash, and values equal to the string for which you would like to print in rowName (for example protein_id to acc_gene translation)
#        $appendName - name to append to structure to differentiate between different constructs
#        $outputFile - name of file for which to append these structures to
# Kristen Naegle
# November 9, 2008
sub printStructToFile($$$$$){
    my ($topHash, $dataHash, $convHash, $appendName, $outputFile) = @_;

    open(OM, ">>$outputFile") || die "Can't open $outputFile for writing\n";
    
    my @keysOfData;

    # print the order of column labels
    print OM "colNames_".$appendName." = {";
    my $first = 0;
    foreach my $top (sort(keys %$topHash)){
	my $nextArr = $topHash->{$top};
	my $next;
	foreach $next (@$nextArr){
	    print OM "'$next'\n";
	}
	if(!$first){
	    my $dataTop = $dataHash->{$top};
	    my $dataNext= $dataTop->{$nextArr->[0]};
	    @keysOfData = keys %$dataNext;
	  #  print "Found Keys: @keysOfData\n";
	    $first = 1;
	}
    }
    print OM "}\n";

    #print the row names
    print OM "rowNames_".$appendName."= {";
    foreach my $id (@keysOfData){
	print OM "'".$convHash->{$id}."'\n";
	
    }
    print OM "}\n";
    
    # print the array of data - for each id, then for each metric 
    print OM "values_".$appendName."=[\n";
    foreach my $id (@keysOfData){
	foreach my $top (sort(keys %$topHash)){
	    my $nextHash = $dataHash->{$top};
	    my $nextArr = $topHash->{$top};
	    foreach my $next (@$nextArr){
		my $bHash = $nextHash->{$next};
		print OM $bHash->{$id}."\t";
	    }
	    
	}
	print OM "\n";

    }
    print OM "]\n";


    
    close(OM);

}

# pAdj = returnFDRAdjPValue($pValArr, $m, $alpha)
# Returns p* given a pValue array where m tests were tested with target alpha
# Inputs: $pValArr - array of pvalues returned from test
#         $m - number of hypotheses tested
#         $alpha - target false positive rate
# Outputs: $pAdj - pvalue for adjustment
# Kristen Naegle
# November 25, 2008
sub returnFDRAdjPValue($$$){
    my($pValArr, $m, $alpha) = @_;
    
#     if(scalar(@$pValArr) < $m){
# 	handleError('returnFDRAdjPValue', 'You do not have a pvalue for every hypothesis tested',\@_);
# 	return -1;
#     }

    my @sort = sort {$a<=>$b} @$pValArr;
    my $pAdj = 0; 
   # print "DEBUG: Sorted array @sort\n";
   # for (my $i=$#sort; $i >=0; $i--){
    for (my $i=0; $i <= $#sort; $i++){
	
	my $pi = $sort[$i]; 
	my $compare = $alpha / $m;
	$compare = $compare * ($i + 1); # conver from zero to ones based system
	#print "DEBUG FDR: Comparing $pi to $compare\n";
	if( $pi <= $compare){
	    $pAdj = $pi;

	}
	else{
	    return $pAdj;
	}

    }

    return $pAdj; 

}

# pAdj = returnBFAdjPValue($m, $alpha)
# Returns p* given a pValue array where m tests were tested with target alpha
# Inputs: $m - number of hypotheses tested
#         $alpha - target false positive rate
# Outputs: $pAdj - pvalue for adjustment ($alpha/$m)
# Kristen Naegle
# March 18, 2009
sub returnBFAdjPValue($$){
    my($m, $alpha) = @_;
    
    my $pAdj = $alpha/$m;
    return $pAdj; 

}


# returnAdjPValue($CORRECTION_TYPE, $pValArr, $m, $alpha)
# Calls given function of correction_type to find the adjusted pValue given m hypothesis tests and a target error rate of alpha
# Inputs: $CORRECTION_TYPE - has options:
#                            BF - BonFerroni correction
#                            BH - Benjamini Hodges FDR correction
#                            NONE - no adjustment, returns p*=alpha
#        $pValArr - pValArr of all tests m
#        $m - number of hypotheses tested
#        $alpha - target p-value with multiple tests
# Outputs: $pAdj - adjusted p-value for which to allow 
# Kristen Naegle
# November 25, 2008
sub returnAdjPValue($$$$){
    my ($TYPE, $pValArr, $m, $alpha) = @_;
    my $pAdj; 
    if($TYPE eq 'BH'){
	$pAdj = returnFDRAdjPValue($pValArr, $m, $alpha);
	return $pAdj;
    }
    if($TYPE eq 'BF'){
	$pAdj = returnBFAdjPValue($m, $alpha);
	return $pAdj;
    }
    if($TYPE eq 'NONE'){
	$pAdj = $alpha;
	return $pAdj;
    }
    else{
	handleError('returnAdjPvalue', "Don't recognize correction algorithm type of $TYPE", \@_);
	return -1;
    }
    
    

}

# $minCount = returnMinCount($fgnd, $fraction)
# Given an array of forground items, return the min count based on desired fractional hit
# Inputs: $fgnd - ref. to array of objects
#         $fraction - amount of foreground that must hit
# Outputs: $minCount - resulting fractional number rounded
# Kristen Naegle
# April 1, 2009
sub returnMinCount($$){
    my($fgnd, $fraction) = @_;
    my $numpept = int(scalar(@$fgnd)/$fraction);
    if($numpept < 2){
	$numpept = 2;
    }

    return $numpept;
} 

# Prints an enrichment report for a file for all enrichments
# Inputs: $dbh - database handle
#         $inputFile - tab separated cluster file
#         $clusterSetNum - middle portion of cluster name in cluster:clusterSetNum:description
#         $expId - experiment id that cluster comes from
#         $threshold - alpha value
#         $countUnlabeled - 1 if you want to include all values regardless of labels
#         $numAA - number of amino acids
#         $CORRECTION_TYPE - type of correction, currently 'BH' and 'none'
#         $fraction - fraction that must match a test in the foreground (for example 2 if you want half the foreground to match
#         $outputFile - output file for enrichment report
# Kristen Naegle
# July 28, 2009 - moved into a function  
sub printEnrichmentReportForClusterSet($$$$$$$$$$$$){
    my ($dbh, $inputFile, $clusterSetNum, $expId, $alpha, $countUnlabeled, $numAA, $CORRECTION_TYPE, $fraction, $outputFile, $enrichType, @enrichArgs) = @_;

    my $msIds = returnMSIdsForExpId($dbh, $expId);
# Double check that file MS ids match the experiment id
    system("rm $outputFile");
    system("touch $outputFile");
  #  open(O, ">>$outputFile") || die "can't open output $outputFile for appending\n";
    print "Found ".scalar(@$msIds)." MS Ids in experiment # $expId\n";
	my ($cHash, $label) = returnClusterHash($inputFile, $clusterSetNum);
    my @clusters = keys %$cHash;
    print "FOUND ".scalar(@clusters)." clusters with label $label\n";
    my @msIdCheck;
    foreach my $c (@clusters){
	push @msIdCheck, @{$cHash->{$c}};
	
    }
    my $CHECK = 0;
    print "In those clusters found ".scalar(@msIdCheck)." ms ids\n";
    print "DEBUG: alpha=$alpha\n";
    if(@msIdCheck == @$msIds){
	print "PASS: experiment matches clusters\n";
	foreach my $c (@clusters){
	    open(O, ">>$outputFile") || die "Can't open $outputFile for writing\n";
	    my @fgndMS = @{$cHash->{$c}};
	    print O "-----------CLUSTER $c ENRICHMENT-----------\n";
	    print O "Size: ".scalar(@fgndMS)." MS phosphopeptides\n";
	    #my ($hashEnrichment, $hashCorrected, $numTests, $fgndHash) = calculateAllEnrichments($dbh, $msIds, \@fgndMS, $countUnlabeled, $threshold, $numAA, $CORRECTION_TYPE, $fraction);
	    my ($error, $hashEnrichment, $hashCorrected, $numTests, $fgndHash) = returnEnrichmentHashForReport($dbh, $msIds, \@fgndMS, $countUnlabeled, $alpha, $numAA, $CORRECTION_TYPE, $fraction, $enrichType, $CHECK, @enrichArgs);
	    if($error){
		print "Ran into error!\n";
		print "Can't print enrichment report\n";
		exit;
	    }
	    foreach my $type (keys %$hashEnrichment){
		print O "ENRICHED FOR: $type\n"; #have to just ignore error for write to a closed file, otherwise we'll lose compatability.
		close(O);
		printEnrichmentReport($hashCorrected->{$type}, $numTests->{$type}, $fgndHash->{$type}, $outputFile);
		
	    }
	  #  print "\n\n";
	}
	
    }
    else{
	print "FAIL: experiment ids do not match cluster ids\n";
	exit;
    }
    
}

# Prints an enrichment report for a file for all enrichments, for cluster #1 only, special case for supersets.
# Inputs: $dbh - database handle
#         $inputFile - tab separated cluster file
#         $clusterSetNum - middle portion of cluster name in cluster:clusterSetNum:description
#         $expId - experiment id that cluster comes from
#         $threshold - alpha value
#         $countUnlabeled - 1 if you want to include all values regardless of labels
#         $numAA - number of amino acids
#         $CORRECTION_TYPE - type of correction, currently 'BH' and 'none'
#         $fraction - fraction that must match a test in the foreground (for example 2 if you want half the foreground to match
#         $outputFile - output file for enrichment report
# Kristen Naegle
# July 28, 2009 - moved into a function  
sub printEnrichmentReportForClusterSet_clusterOneOnly($$$$$$$$$$$$){
    my ($dbh, $inputFile, $clusterSetNum, $expId, $alpha, $countUnlabeled, $numAA, $CORRECTION_TYPE, $fraction, $outputFile, $enrichType, @enrichArgs) = @_;

    my $msIds = returnMSIdsForExpId($dbh, $expId);
# Double check that file MS ids match the experiment id
    system("rm $outputFile");
    system("touch $outputFile");
  #  open(O, ">>$outputFile") || die "can't open output $outputFile for appending\n";
    print "Found ".scalar(@$msIds)." MS Ids in experiment # $expId\n";
	my ($cHash, $label) = returnClusterHash($inputFile, $clusterSetNum);
    #my @clusters = keys %$cHash;
    my @clusters = (1);
    print "FOUND ".scalar(@clusters)." clusters with label $label\n";
    my @msIdCheck;

    my $CHECK = 0;
    print "In those clusters found ".scalar(@msIdCheck)." ms ids\n";
    print "DEBUG: alpha=$alpha\n";
    #if(@msIdCheck == @$msIds){
	print "PASS: experiment matches clusters\n";
	foreach my $c (@clusters){
	    open(O, ">>$outputFile") || die "Can't open $outputFile for writing\n";
	    my @fgndMS = @{$cHash->{$c}};
	    print O "-----------CLUSTER $c ENRICHMENT-----------\n";
	    print O "Size: ".scalar(@fgndMS)." MS phosphopeptides\n";
	    #my ($hashEnrichment, $hashCorrected, $numTests, $fgndHash) = calculateAllEnrichments($dbh, $msIds, \@fgndMS, $countUnlabeled, $threshold, $numAA, $CORRECTION_TYPE, $fraction);
	    my ($error, $hashEnrichment, $hashCorrected, $numTests, $fgndHash) = returnEnrichmentHashForReport($dbh, $msIds, \@fgndMS, $countUnlabeled, $alpha, $numAA, $CORRECTION_TYPE, $fraction, $enrichType, $CHECK, @enrichArgs);
	    if($error){
		print "Ran into error!\n";
		print "Can't print enrichment report\n";
		exit;
	    }
	    foreach my $type (keys %$hashEnrichment){
		print O "ENRICHED FOR: $type\n"; #have to just ignore error for write to a closed file, otherwise we'll lose compatability.
		close(O);
		printEnrichmentReport($hashCorrected->{$type}, $numTests->{$type}, $fgndHash->{$type}, $outputFile);
		
	    }
	  #  print "\n\n";
	}
	
    #
    #else{
#	print "FAIL: experiment ids do not match cluster ids\n";
#	exit;
   # }
    
}

#  ($error, $hashEnrichment, $hashCorrected, $numTests, $fgndHash) = returnEnrichmentHashForReport($dbh, $bgnd, $fgnd, $countUnlabeled, $alpha, $numAA, $CORRECTION_TYPE, $fraction, $enrichType, @enrichArgs);
# handle enrichment based on enrichment type argument.  If only a single enrichment, then @enrichArgs is size 1 and contains a threshold
# Inputs: $dbh - database handle 
#         $bgnd - ref. to array of ms ids for full background
#         $fgnd - ref to array of ms ids for foreground
#         $countUnlabeled - 1 to include all selected
#         $alpha - corrected pvalue
#         $numAA - number of amino acids to consider in motif enrichment
#         $CORRECTION_TYPE - type of correction to perform, none, BH, etc.
#         $fraction - number of things in foreground that must match to run test
#         $enrichType - type of enrichment: 'all', 'sequence', 'peptide_prediction' 'domain' 'site_domain', 'GO'
#         @enrichArgs - if choosing 'all' this is size 3, ($scansiteThreshold, $motifThreshold, $domainThreshold) else it's size 1 and contains a threshold (for GO and domain_site this is a bogus value)
# Outptus: $error - error is 1 if there was  a problem (like number of arguments are wrong or don't recognize enrichmentType
#          $hashEnrichment - enrichment hash
#          $hashCorrected - the corrected enrichment hash according to alpha and correction type
#          $numTests - hash holding number of tests
#          $fgndHash - hash of hits in foreground
# Kristen Naegle
# July 29, 2009
sub returnEnrichmentHashForReport($$$$$$$$$$){
    my($dbh, $bgnd, $fgnd, $countUnlabeled, $alpha, $numAA, $CORRECTION_TYPE, $fraction, $enrichType, $CHECK, $enrichArgs) = @_;
    my($hashEnrichment, $hashCorrected, $numTests, $fgndHash);
    my $error = 0;

   # print "DEBUG: alpha=$alpha\n";
    if($enrichType eq 'all'){
	#if all then require @enrichArgs have stuff in this order:
	if(scalar(@$enrichArgs) != 3){
	    print "ERROR! enrichment requires ScansiteThreshold, MotifThreshold and domainThreshold\n";
	    print "You passed in @$enrichArgs\n";
	    print "DEBUG: $enrichArgs->[0]\n";
	    $error = 1;
	    return($error, $hashEnrichment, $hashCorrected, $numTests, $fgndHash);
	}
	my($scansiteThreshold, $motifThreshold, $domainThreshold) = @$enrichArgs;
	($hashEnrichment, $hashCorrected, $numTests, $fgndHash) = calculateAllEnrichments($dbh, $bgnd,$fgnd, $countUnlabeled, $alpha, $numAA, $CORRECTION_TYPE, $fraction, $scansiteThreshold, $motifThreshold, $domainThreshold)
    }
    else{
	if(scalar(@$enrichArgs) != 1){
	    $error = 1;
	    print "You chose single enrichment, enrichment arguments should be size 1 and contain a threshold value.\n";
	    return($error, $hashEnrichment, $hashCorrected, $numTests, $fgndHash);
	}
	my $threshold = $enrichArgs->[0];


	my @args;
	if($enrichType eq 'GO'){
	    #GO enrichment
	    @args = ($fraction);
	    ($hashEnrichment, $hashCorrected, $numTests, $fgndHash) = calculateSingleEnrichment($dbh, $CHECK, 'GO', $bgnd, $fgnd, $countUnlabeled, $alpha, $CORRECTION_TYPE, @args);
	    
	}
	elsif($enrichType eq 'sequence'){
    
	    @args = ($threshold, $numAA, $fraction);
	    ($hashEnrichment, $hashCorrected, $numTests, $fgndHash) = calculateSingleEnrichment($dbh, $CHECK, 'sequence', $bgnd, $fgnd, $countUnlabeled, $alpha, $CORRECTION_TYPE, @args);
	   
	    
	}
	elsif($enrichType eq 'peptide_prediction'){
	    @args = ($threshold, $fraction);
	    ($hashEnrichment, $hashCorrected, $numTests, $fgndHash) = calculateSingleEnrichment($dbh, $CHECK, 'peptide_prediction', $bgnd, $fgnd, $countUnlabeled, $alpha, $CORRECTION_TYPE, @args);
	}
	elsif($enrichType eq 'domain'){
	    # domain  enrichment
	    @args = ($threshold, $fraction);
	    ($hashEnrichment, $hashCorrected, $numTests, $fgndHash) = calculateSingleEnrichment($dbh, $CHECK, 'domain', $bgnd, $fgnd, $countUnlabeled, $alpha, $CORRECTION_TYPE, @args);
	}
	elsif($enrichType eq 'site_domain'){
	    
	    # site_domain  enrichment
	    @args = ($fraction);
	    ($hashEnrichment, $hashCorrected, $numTests, $fgndHash)= calculateSingleEnrichment($dbh, $CHECK, 'site_domain', $bgnd, $fgnd, $countUnlabeled, $alpha, $CORRECTION_TYPE, @args);
	    
	    
	}
	
	else{
	    print "ERROR: Don't recognize enrichment type $enrichType\n";
	    $error = 1;
	}
	
    }
    return($error, $hashEnrichment, $hashCorrected, $numTests, $fgndHash);

}

# (maxLabel, minLabel, maxPosChangeLabel, maxNegChangeLabel) = returnDataFeatureLabels($dataLabelArr, $dataValues)
# given data arrays of values and corresponding labels, return the labels that correspond to max value, min value, max positive change and max negative change. When no pos or no negative change occurs, return ~~~ label
# Inputs: $dataLabelArr - an array, in order of desired feature, of labels
#         $dataValues - an array of values that corresponds to the same order of the label array
# Outputs: $maxLabel - the labe corresponding to the maximum value
#          $minLabel - the label corresponding to the minimum value
#          $maxPosChangeLabel - this will be label1_label2 indicating the max positive change occurs between label1 value and label2 value
#          $maxNegChangeLabel - this will be label1_label2 indicating the max negative change occurs between label1 value and label2 value
# Kristen Naegle
# January 19, 2010
sub returnDataFeatureLabels($$){
    my($labelArr, $valueArr) = @_;
    my ($max, $min, $maxPos, $maxNeg);
    $maxPos = '~~~';
    $maxNeg = '~~~';
    #create a hash so we can look up the correct label
    my %dHash;
    for(my $i=0; $i<scalar(@$valueArr); $i++){
	my $label = $labelArr->[$i];
	my $value = $valueArr->[$i];
	$dHash{$value} = $label;
	
    }
    my @sortedData = sort {$a<=>$b} @$valueArr;
    my $maxVal = $sortedData[$#sortedData];
    my $minVal = $sortedData[0];
    $max = $dHash{$maxVal};
    $min =$dHash{$minVal};


    # do the same for maximum change
    my %cHash;
    for(my $i=1; $i<scalar(@$valueArr); $i++){
	my $label1 = $labelArr->[$i-1];
	my $label2 = $labelArr->[$i];
	my $newLabel = $label1.'_'.$label2;
	my $value = $valueArr->[$i]-$valueArr->[$i-1];
	$cHash{$value} = $newLabel;
    }
    my @changes = keys %cHash;
    #print "DEBUG: differential array is: @changes\n";
    @changes = sort {$a <=> $b} @changes;
    if($changes[$#changes] > 0){
	$maxPos = $cHash{$changes[$#changes]};
	
    }
    if($changes[0] < 0){
	$maxNeg = $cHash{$changes[0]};
    }

    return($max, $min, $maxPos, $maxNeg);

}

# \%featureHash = returnDataFeatureLabelsForMSIds($dbh, $MSIdArr);
# for an array of MS Ids get the max, min, maxPos, and maxNeg for an MSId. Returns a hash based on MS.id that points to a hash with keys equal to the quantitative data type and then values equal to the label of that type.
# Inputs: $dbh - database handle, where we get the data from for an MS id
#         $MSIdArr - array of MS ids. Ms.id
# Outputs: $featureHash - ref to hash, with keys equal to each MS id and value a hash of features
# This assumes that there is only a single run.  Must reshape runs for clustering, so this is a fine assumption.
# Kristen Naegle
# January 19, 2010
sub returnDataFeatureLabelsForMSIds($$){
    my ($dbh, $MSIdArr) = @_;
    my ($sth_data, $cHash) = returnDataForMSIdSTH($dbh);
    my %MSHash;
    # get the labels and data for each MS id and then get features
    foreach my $MSId (@$MSIdArr){
	my ($dataHash, $labels) = returnDataForMSIdNoStddev($sth_data, $MSId, $cHash);
	my @runs = keys %$dataHash;
	if(scalar(@runs) > 1){
	    print "ERROR: You have more than one run. Do not know how to handle/do not handle for returnDataFeatureLabelsForMSIds\n";
	    exit;
	}
	my $dataArr = $dataHash->{$runs[0]};
	#print "DEBUG: size of labels:".scalar(@$labels)." and data: ".scalar(@$dataArr)."\n";
	my($max, $min, $maxPos, $maxNeg) = returnDataFeatureLabels($labels, $dataArr);
	my %featureHash;
	$featureHash{'maxValue'} = $max;
	$featureHash{'minValue'} = $min;
	$featureHash{'maxPosChange'} = $maxPos;
	$featureHash{'maxNegChange'} = $maxNeg;
	$MSHash{$MSId} = \%featureHash;
    }
    return \%MSHash;
}

# (\%featureHash, $N) = returnDataFeatureHashForMSIds($dbh, $MSIdArr);
# This will create the data features for MS ids but then convert to counts so that enrichment can be calculated
# Inputs: $dbh - database handle
#         $MSIdArr -ref to array of MS ids
# Outputs: $featureHash - hash of hashes of data features.  Top key are the data features (like maxValue) and then second hash has keys equal to teh labels and values equal to the number of times in the MSIdArr that occurs
#          $N - the number of features
# Kristen Naegle
# January 19, 2010
sub returnDataFeatureHashForMSIds($$){
    my ($dbh, $MSIdArr) = @_;
    my $MSIdHash = returnDataFeatureLabelsForMSIds($dbh, $MSIdArr);
    my %newHash;
    foreach my $MSId (keys %$MSIdHash){
	my $featureHash = $MSIdHash->{$MSId};
	foreach my $feature (keys %$featureHash){
	    my $label = $featureHash->{$feature};	    
	    if(not defined($newHash{$feature})){
		my %newLabel;
		$newHash{$feature} = \%newLabel; 
	    }
	    
	    if(not defined ($newHash{$feature}->{$label})){
		my @array = (0);
		$newHash{$feature}->{$label} = \@array;
		
	    }
	    ($newHash{$feature}->{$label})->[0] += 1;
	}
    }
    my $n = scalar(@$MSIdArr);
    return (\%newHash, $n);
}


# ($hashRef, $hashRefCorrected, $numHypothesesTested, $fgndHash) = calculateDataFeatureEnrichment($dbh, $fgnd, $bgnd, $countUnlabeled, $minCount, $threshold, $CORRECTION_TYPE, $alpha);
# For a list of foreground MS Ids and a list of background MS Ids, report enrichment for terms occuring n times or more
# Returns a hash of hashes like the input foreground, except instead of count, there is a calculated probability value
# Inputs: $dbh - database handle
#         $fgnd - reference to array of MS ids of foreground
#         $bgnd - reference to array of MS ids of bacgkround
#         $countUnlabeled - flag: 1 - include unlabeled in calculations 0:exclude ~~~ and --- unlabeled 
#         $minCount - min. count in fgnd required to calculate p-val
#         $threshold - min pvalue required to report significance of site Domain enrichment
#         $CORRECTION_TYPE - type of Multiple Hypotehsis testing correction
#         $alpha - target correction p-value
# Outputs: $hashRef - reference to hash of hashses. Top level key is prediction type, second level key is term and value is hypergeometric calculation
#          $correctedHashRef - reference to hash of hashes - same as hashRef by corrected for MHC
#         $numHypothesesTested - ref. to hash with top level key prediction type, and value equal to teh number of hypotheses tested for that 
#         $fgndHash - returns the fgnd hash with all terms and counts for those terms
# Kristen Naegle
# January 19, 2010
sub calculateDataFeatureEnrichment($$$$$$$$){
    my ($dbh, $fgnd, $bgnd, $countUnlabeled, $minCount, $threshold, $CORRECTION_TYPE, $alpha) = @_;
    my $uniquify = 0;
    
    my ($bgnd_hash, $N) = returnDataFeatureHashForMSIds($dbh, $bgnd);
    my ($fgnd_hash, $n) = returnDataFeatureHashForMSIds($dbh, $fgnd);
    my ($enrichmentHash, $numHypothesesTested) = calculateEnrichmentFromHash($dbh, $n, $N, $fgnd_hash, $bgnd_hash, $minCount, $threshold, $countUnlabeled);

    my $correctedEnrichmentHash = pruneEnrichmentHash($enrichmentHash, $numHypothesesTested, $CORRECTION_TYPE, $alpha);

    return ($enrichmentHash, $correctedEnrichmentHash, $numHypothesesTested, $fgnd_hash);

}


1;
