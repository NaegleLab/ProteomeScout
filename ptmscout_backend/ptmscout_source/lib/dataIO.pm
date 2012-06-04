use strict;
use warnings;
use DBTools::insertFunctions;
use globalVars;
use entrezTools;
use List::Util qw(sum);

#(\%hash, $labelArr) = returnDataHashForExpId($dbh, $expId)
# returns an hash with keys=ms_id and arrays of data and a string description of the order of the data. 
# Inputs: $dbh - database handle
#         $expId - experiment Id 
# Outputs: \%hash - reference to hash, keys are MS_ids and values are an array of data
#          $labelArr - an array of labels describing the order of data in values of hash
# Kristen Naegle
# May ?, 2008
sub returnDataHashForExpId($$){
    my ($dbh, $expId) = @_;
    my (%hash, @names); 

# #     # first return the MS_ids in an experiment, foreach of those get their data
# #     # find out what data column is (time, condition?). Get the labels 
# #     # handle NA (if value is zero and NA is 1..return NaN)
# #     #combine type:label 

    my $MS_ids = returnMSIdsForExpId($dbh, $expId);
    my ($sth_data, $cHash) = returnDataForMSIdSTH($dbh); #includes run
    my $labelArr;
    foreach my $MS_id (@$MS_ids){
	#if(not defined($hash{$MS_id})){
	 #   $hash{$MS_id} = []; #handling runs?
	#}
	my ($dataHash, $labels) = returnDataForMSId($sth_data, $MS_id, $cHash); #create this funciton and test (especially for NaN's)
	#push @{$hash{$MS_id}}, @$data;
	$hash{$MS_id} = $dataHash;
	$labelArr = $labels;
    }
    return(\%hash, $labelArr);



}


# printMinDataDescForExpId($dbh, $expId, $outputFile)
# Prints a text file to outputFile that has first column of MS_ids and remaining columns are data (header describes)
# Inputs: $dbh - database handle
#         $expId - experiment id to write
#         $outputFile - destination of data
#         $AVG - if this value is one, the same data points in each run will be averaged
# Kristen Naegle
# May ?, 2008
# November 14, 2008 - changing to output gene_site1_site2.. as another data header
sub printMinDataDescForExpId($$$$){
    my ($dbh, $expId, $outputFile, $AVG) = @_;
    my($expHash, $labels) = returnDataHashForExpId($dbh, $expId);
    open(FH_D, ">$outputFile") || die "Can't open $outputFile for writing\n";
    #print the header
    print FH_D "MS_id\tgene_site";
    ###TEMP
    print FH_D "\tpep\trun";
    ###END TEMP
    my $runArray = returnRunArray($dbh, $expId);
    ($labels) = returnDataLabelsReshaped($labels, $runArray, $AVG);


    for(my $i=0; $i < scalar(@$labels); $i++){
	print FH_D "\tdata_".$labels->[$i];
    }
   
    print FH_D "\n";

    # get a translation hash from MS_id to a gene_site description
    my @MSIdArr = keys %$expHash;
    my $transHash = translateMSIdsToDescriptor($dbh, \@MSIdArr);
    # print the data lines
    foreach my $MS_id (keys %$expHash){
	my $dataHash = $expHash->{$MS_id};	
	$dataHash  = returnDataHashReshaped($dataHash, $runArray, $AVG);

	foreach my $run (keys %$dataHash){ 
	    print FH_D "$MS_id\t".$transHash->{$MS_id};
	    my $r = returnMSInfoById($dbh, $MS_id);
	    print FH_D "\t".$r->{'phosphopep'};
	    print FH_D "\t".$run;
	    my $data = $dataHash->{$run};
	    for(my $j=0; $j < scalar(@$data);  $j++){
		print FH_D "\t".$data->[$j];
	
	    }

	    print FH_D "\n";
	}
    }
    
    
    close(FH_D);
}

# ($newLabels) = returnDataLabelsReshaped($labels, $runArr, $AVG)
# When averaging, need to add stdev labels to label array
# Inputs: $labels - ref. to array of labels 
#         $runArr - the order of the runs to append to the label array
#         $AVG - boolean value - if 
# Outputs: $newLabels - ref to array of labels, has incoming with appended stdev labels, or has an array of run types appended to the data type
# Kristen Naegle
# Jan. 10, 2009
sub returnDataLabelsReshaped($$$){
    my ($labels, $runArr, $AVG) = @_;
    my @stdevLabels; 
    my @newLabels;

    if($AVG){
	for(my $i=0; $i<scalar(@$labels); $i++){
	    #going to append stddev to the middle field
	    my $label = $labels->[$i];

	    my $newLabel = "stddev-".$label; #"stddev-".$label[0]."_".$label[1];
	    push @stdevLabels, $newLabel;
	    
	}

	@newLabels = (@$labels, @stdevLabels);
    }
    else{
	# if not averaging, then reshape the label array to append the run type
	for(my $i=0; $i<scalar(@$runArr); $i++){
	    my $run = $runArr->[$i];
	    for(my $j=0; $j<scalar(@$labels); $j++){
		my $newLabel = $run."-".$labels->[$j];
		push @newLabels, $newLabel;
	    }

	}
	
    }
    return \@newLabels;

}




# $newDataHash = returnDataHashReshaped($dataHash, $runs, $AVG)
# when averaged required for multiple runs, create a new hash that has run type average and includes stddevs.  Stddevs are listed in the same order as primary data, but appended to the last half. If not averaging, simply reshape according to order indicated by the run array.
# Inputs: $dataHash - hash with keys equal to runs and values equal to array of @data in desired order
#         $runs - ref. to array of runs in the order desired.  When runs are mssing 'NA' is put in place.
#         $AVG - boolean value that says whether ot average or reshape data
# Outputs: $dataHash - same structure, but all runs are replaced or reshaped
# Kristen Naegle
# Jan 10, 2010
sub returnDataHashReshaped($$$){
    my ($dataHash, $runArr, $AVG) = @_;
    my %newHash;
    

    if($AVG){
	#get how many data points we'll need to walk through
	my @runs = keys(%$dataHash);
	my $dataEx = $dataHash->{$runs[0]};
	my @newData;
	my @newStdDev;
	for(my $i=0; $i<scalar(@$dataEx); $i++){
	    my @dataPt; 
	    foreach my $run (@runs){
		my $d = $dataHash->{$run};
		print "DEBUG: data point: ".$d->[$i]."\n";
		push @dataPt, $d->[$i]; ##add handle here to only push on if $d->[$i] is a value
	    }
	    my($avg, $stddev) = returnAvgStddevDataArray(\@dataPt);
	    push @newData, $avg; #pushing on average;
	    push @newStdDev, $stddev;
#	print "DEBUG: size of data now is".@newData."\n";
	}
	my @dataAll = (@newData, @newStdDev);
	$newHash{'average'} = \@dataAll;
    }
    else{
	# first find the number of data points
	my @newData;
	my $numPoints = 0;
	my $runCount = 0;
	while(!$numPoints){
	    my $run = $runArr->[$runCount];
	    if(defined($dataHash->{$run})){
		my $dataArr = $dataHash->{$run};
		$numPoints = scalar(@$dataArr);
	    }
	    else{
		$runCount+=1;
	    }
	} #should now have $numPoints
	print "DEBUG: number of data points $numPoints\n";
	for(my $i=0; $i<scalar(@$runArr); $i++){
	    my $run = $runArr->[$i];
	    if(defined($dataHash->{$run})){
		
		push @newData, @{$dataHash->{$run}}; 
	    }
	    else{
		for(my $j=0; $j<$numPoints; $j++){
		    push @newData, 'NaN';
		}
	    }

	}
	
	$newHash{'reshaped'} = \@newData;
    }

    return \%newHash;

}

# ($avg, $stddev) = returnAvgStddevDataArray($dataArray);
# returns the arithmetic mean and stddeviation of a data array
# Inputs: $dataArray - reference to an array, of variable length
# Outputs: $avg - arithmetic average of the data array
#          $stddev - standard deviation of array
# Kristen Naegle
# Jan 10, 2010
sub returnAvgStddevDataArray($){
    my ($dataArr) = @_;
    
    my $avg = sum(@$dataArr)/@$dataArr;
    
    my $sum = 0;
    my $stddev; 
    if(scalar(@$dataArr) == 1){
	$avg = $dataArr->[0];
	$stddev = 'NA';
	return ($avg, $stddev);
    }
    foreach my $pt (@$dataArr){
	my $m = $pt-$avg;
	
	$sum += $m*$m;
	#print "DEBUG: sum at $pt $sum\n";

    }
    $sum = $sum/@$dataArr;
#    print "DEBUG: dividing by ".@$dataArr."\n";
    $stddev = sqrt($sum);
    #force them to 4 floating point sig digs
    $avg = sprintf("%0.4f",$avg);
    $stddev = sprintf("%0.4f", $stddev);
    return($avg, $stddev);

}

# $MSTableHash = returnMSInfoById($dbh, $MSId)
# A small function I needed for debug
# Given an MS id return MS info hash (in order to print phosphopep)
# Inputs: $dbh- database handle
# Outputs: $MSId - id to MS table
# Kristen Naegle
# July 8, 2009
sub returnMSInfoById($$){
    my ($dbh, $MSId) = @_;
    my $sth = $dbh->prepare('SELECT * from MS where id=?');
    $sth->execute($MSId);
    my $result = returnMultipleResultsForExSTH($sth);
    my $r = $result->[0];
    return $r;
    

}

# \%convHash = translateMSIdsToDescriptor($dbh, $MSIdArr)
# Assumes MSIdArr is an array of unique MS ids for which to translate.  Returns a hash with keys equal to the MSIds and values a descriptor of accGene_site1_site2..etc.
# Inputs: $dbh -database handle
#         $MSIdArr - ref. to array of MSIds to translate
# Outputs: \%convHahs - ref. to hash of descriptors for keys MSIds
# Kristen Naegle
# November 14, 2008
sub translateMSIdsToDescriptor($$){
    my ($dbh, $MSIdArr) = @_;
    my %convHash;
    my $sth = $dbh->prepare('SELECT phosphopep.*, protein.acc_gene from MS_phosphopep join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id join protein on phosphopep.protein_id=protein.id where MS_phosphopep.MS_id=?');
    foreach my $MSId (@$MSIdArr){
	$sth->execute($MSId);
	my $result = returnMultipleResultsForExSTH($sth);
	my $gene = ($result->[0])->{'acc_gene'};
	my $descStr = $gene;
	foreach my $r (@$result){
	    my $site_type = $r->{'site_type'};
	    my $site_pos = $r->{'site_pos'};
	    $descStr .= "_".$site_type.$site_pos;
	}
	$convHash{$MSId} = $descStr;
    }

    return \%convHash;
}


# printExperimentOutput($dbh, $outputFile, $expId)
# Prints an output file of experiment (just the aligned peptides and data), no predictions or GO terms 
# Inputs: $dbh - databae handle
#         $outputFile - name of file to write to
#         $expId - id to experiment table
# Kristen Naegle
# June 18, 2008
# Modified January 26, 2009 - to include print for scansite predictions, GO annotations and pfam domains
# Modified March 16, 2009 - to print only one line per MS Id
sub printExperimentOutput($$$){
    my ($dbh, $outputFile, $expId) = @_;
    #get the MS_ids
    my($dataHash, $labels) = returnDataHashForExpId($dbh, $expId);
    my %MSSeen;
    open(OUT_EXP, ">$outputFile") || die "Can't open $outputFile for writing\n";
    my $hashRef = returnMSInfoByExpId($dbh, $expId);
    my $header = "MS.id\tpep:tryps\tacc:gene\tacc\tname:long\tsite\tpep:aligned";
    for(my $i=0; $i < scalar(@$labels); $i++){
	$header .= "\tdata:".$labels->[$i];
    }
    $header .= "\tscansite:kinase\tscansite:bind\n";
    print OUT_EXP $header;

    foreach my $key (keys %$hashRef){
	my $sHash = $hashRef->{$key};
	my $MS_id = $sHash->{'MS.id'};
	if(not defined $MSSeen{$MS_id}){
	    $MSSeen{$MS_id} = 1;
	    print OUT_EXP "$sHash->{'MS.id'}\t$sHash->{'MS.phosphopep'}\t$sHash->{'protein.acc_gene'}\t$sHash->{'acc.value'}\t$sHash->{'protein.name'}\t";
	    # get aligned peptides
	    # get data
	    my $aHash = returnAlignedPeptidesForMSId($dbh, $MS_id);
	    my $siteType;
	    my $aligned;
	    foreach my $k (keys %$aHash){
		my $aRef = $aHash->{$k};
		$siteType .= $aRef->{'phosphopep.type'}.$aRef->{'phosphopep.site'}.';';
		$aligned .= $aRef->{'phosphopep.aligned'}.';';
	    }
	    print OUT_EXP "$siteType\t$aligned";
	    my $data = $dataHash->{$MS_id};
	    for(my $j=0; $j < scalar(@$data);  $j++){
		print OUT_EXP "\t".$data->[$j];
	    }
	    # get scansite predictions
	    my ($kinase, $bind) = returnScansitePredictions($dbh, $MS_id);
	    print OUT_EXP "\t".$kinase."\t".$bind;
	    # get pfam domains
	    
	    # get site domain
	    
	    print OUT_EXP "\n";
	}


    }
 

}

# ($kinase, $bind) = returnScansitePredictionsInTxt($dbh, $MSId)
# Given an MS.id, find all scansite predicitions for the aligned phosphopeptides, assumes order by site number in return
# Inputs: $dbh - database handle
#         $MSId - MS.id 
# Outputs: $kinase - scansite kinase predictions of the type kinase:score, kinase:score; kinase:score (peptides seperated by semicolon
#          $bind - same as above, only binding predictions 
# Kristen Naegle
# January 27, 2009
sub returnScansitePredictions($$){
    my ($dbh, $MSId) = @_;
    
    my $sth = $dbh->prepare("SELECT phosphopep.id, value, score from MS join MS_phosphopep on MS.id=MS_phosphopep.MS_id join phosphopep on phosphopep.id=MS_phosphopep.phosphopep_id join phosphopep_prediction on phosphopep.id=phosphopep_prediction.phosphopep_id where MS.id=? and source='scansite_kinase' order by site_type, score");
    $sth->execute($MSId);
    my ($resultsKinase) = returnMultipleResultsForExSTH($sth);
    my $kinase = scansiteLineMaker($resultsKinase);
    $sth = $dbh->prepare("SELECT phosphopep.id, value, score from MS join MS_phosphopep on MS.id=MS_phosphopep.MS_id join phosphopep on phosphopep.id=MS_phosphopep.phosphopep_id join phosphopep_prediction on phosphopep.id=phosphopep_prediction.phosphopep_id where MS.id=? and source='scansite_bind' order by site_type, score");
    $sth->execute($MSId);
    my ($resultsBind) = returnMultipleResultsForExSTH($sth);
    my $bind = scansiteLineMaker($resultsBind);
    return($kinase, $bind);

}
# ($text) = scansiteLineMaker($resultArr)
# Given an array of hash results from a scansite return, append everything into a single tab line entry
# Inputs: $resultArr - ref. to array of result hash - see returnScansitePredictions, has fields id (phosphopep), value and score
# Outputs: $text - diff. phosphopep predictions separated by ; and same are separated by , with value:score for each
# Kristen Naegle
# January 27, 2009
sub scansiteLineMaker($){
    my ($results) = @_;
    my $r1 = $results->[0];
    my $id;
    my $line;
    if(defined($r1->{'id'})){
	$id = $r1->{'id'};
    }
    foreach my $r (@$results){
	if($r->{'value'} ne '-1'){
	    $line .= $r->{'value'}.":".$r->{'score'};
	
	    if($r->{'id'} == $id){
		$line .= ",";
	    }
	    else{
		$line .=";";
	    }
	    $id = $r->{'id'};
	}
	else{
	    $line .= '~~~'; 
	}
    }
    return $line;
}

# printFastFile($dbh, $pepHash, $outputFile)
# Print the FASTA format file from a hash, where keys are the description and values are the sequence
# Inputs: $dbh - database handle
#         $pepHash - hash with keys equal to desired header and values equal to sequence
#         $outputFile - file to write fasta formatted output to
# Kristen Naegle
# July 19, 2008
sub printFastaFile($$$){
    my($dbh, $pepHash, $outputFile) = @_;  

    open(PEP_FASTA, ">$outputFile") || die "Can't open $outputFile for writing FASTA sequences\n";
    foreach my $id (keys %$pepHash){
	print PEP_FASTA ">$id\n";
	my $pep = returnAlignedSpace($pepHash->{$id});

	$pep =~ s/ /_/g;
	print PEP_FASTA $pep."\n";
    }
    close(PEP_FASTA);
}

# printPhosphoPepsToFasta($dbh, $MSidArr, $outputFile)
# Given an array of MS ids, print aligned peptides to a fasta file
# Inputs: $dbh - database handle
#         $MSidArr - array of MS ids 
#         $outputFile - file to write fasta output to
# Kristen Naegle
# July 19, 2008
sub printPhosphoPepsToFasta($$$){
    my ($dbh, $MSidArr, $outputFile) = @_;
    my %hash;
    foreach my $msId (@$MSidArr){
	my ($pepIds, $peps) = returnPhosphopepTableForMSId($dbh, $msId);
	for (my $i=0; $i < scalar(@$pepIds); $i++){
	    my $pepId = $pepIds->[$i];
	    my $pep = $peps->[$i];
	    if(not defined $hash{$pepId}){
		$hash{$pepId} = $pep;
	    }
	}
    }
    printFastaFile($dbh, \%hash, $outputFile);
 }

# printPhosphoPepsSitesSpeciesToFasta($dbh, $MSidArr, $siteType, $species, $outputFile)
# Given an array of MS ids, print aligned peptides to a fasta file, only for a particular site type and species.
# Inputs: $dbh - database handle
#         $MSidArr - array of MS ids 
#         $siteType - type of site 'Y', 'S' or 'T'
#         $species - two word name of species from protein table
#         $outputFile - file to write fasta output to
# Kristen Naegle
# July 23, 2008
sub printPhosphoPepsSiteSpeciesToFasta($$$$$){
    my ($dbh, $MSidArr, $siteType, $species, $outputFile) = @_;
    my %hash;
    foreach my $msId (@$MSidArr){
	my ($pepIds, $peps) = returnPhosphopepTableSiteSpeciesForMSId($dbh, $msId, $siteType, $species);
	for (my $i=0; $i < scalar(@$pepIds); $i++){
	    my $pepId = $pepIds->[$i];
	    my $pep = $peps->[$i];
	    if(not defined $hash{$pepId}){
		$hash{$pepId} = $pep;
	    }
	}
    }
    printFastaFile($dbh, \%hash, $outputFile);
 }

# $aligned_spaced = returnAlignedSpace($peptide);
# Given a pep_aligned, return a value with 15 characters, spaces filling in for N and C terminals
# Inputs: $peptide - phosphopep.pep_aligned (i.e. something already aligned, but missing spaces for c/n-term cutoff
# Ouputs: $aligned_spaced - 15mer
# Kristen Naegle
# July 19, 2008
sub returnAlignedSpace($){
    my ($pep) = @_;
    my $length = length($pep);
    my $spaced_pep;
    if($length==15){
	$spaced_pep = $pep;
    } 
    else{
	if($pep =~ /([s|y|t])/){
	    my $match = $1;
	    #print "MATCHED: $match\n";
	    my $pos = rindex($pep, $match);
	    my $diff = 7 - $pos;
	    if($diff > 0){
		my $spaces;
		for (my $i=1; $i <= $diff; $i++){
		    $spaces .= " ";
		}
		$spaced_pep = $spaces.$pep;
	    }
	    else{
		my $spaces;
		for(my $i=1; $i<= (15-$length); $i++){
		    $spaces .= " ";
		} 
		$spaced_pep = $pep.$spaces;
	    }
	}
	else{
	    $spaced_pep = $pep;
	    print "ERROR: No aligned position to be found\n";
	    
	}
	
	
    }
    return $spaced_pep;
}

# printProteinFile($dbh, $expId, $outputFile)
# Prints protein information to file for an experiment. Currently handles just the GO terms
# Inputs: $dbh - database handle
#         $expId - eperiment id
#         $outputFile - place to write .txt file
# Kristen Naegle
# July 17, 2008
sub printProteinFile($$$){
    my($dbh, $expId, $outputFile) = @_;
    open(PROT_OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
    print PROT_OUT "acc\tname\tGO:F\tGO:C\tGO:P\n";
    # for each unique protein Id print a F, C and P GO term column (GO_term:lable; 
    my $proteinIds = returnProteinIdsForExpId($dbh, $expId);
    my @aspects = ('F', 'C', 'P');
    foreach my $proteinId (@$proteinIds){
	my $accArr = returnAccValuesByProteinId($dbh, $proteinId);
	my $accHash = returnAccHashByTypeFromArr($accArr);
	my $accGI = $accHash->{'gi'};
	my $acc = $accGI->[0];
	my $protHash = returnProteinDescForProteinId($dbh, $proteinId);
	print PROT_OUT "$acc\t$protHash->{'name'}";
	
	my $hashRef = returnGOHashForProteinId($dbh, $proteinId);
	foreach my $aspect (@aspects){
	    if(not defined($hashRef->{$aspect})){
		print PROT_OUT "\t~~~";
		
	    }
	    else{
		my $GOString = returnGOString($hashRef->{$aspect});
		print PROT_OUT "\t$GOString";
	    }
	}
	print PROT_OUT "\n";
	
    }

    close(PROT_OUT);
}

# $GoString = returnGOString($GOArr)
# returns a string for of GO terms from array of hashes
# Inputs: $GOArr - ref. to array of GO hashes (keys 'GO' and 'term');
# Outputs: $string - string concatenation of GO terms
# Kristen Naegle
# July 17. 2008
sub returnGOString($){
    my ($GOArr) = @_;
    my $string;
    foreach my $GO (@$GOArr){
	$string .= $GO->{'GO'}.":".$GO->{'term'}.";";
    }
    return $string;
}

# printTableFromHash($fileName, $hash)
# Given a hash that describes the field names and values, make a tab separated table file based on Name.
# Inputs: $fileName - destination of table print
#         $cols - ref. to array of columns, in order of table insert
#         $hashArr - array of hashes with keys equal to field names and values equal to table value
# Kristen Naegle
# Dec. 29, 2009
sub printTableFromHash($$$){
    my ($fileName, $cols, $hashArr) = @_;

    open(OUT, ">$fileName") || die "can't open $fileName for writing\n";
    my $first = $hashArr->[0];
    my @fields = @$cols;
    
    foreach my $hash (@$hashArr){
	#check for totally null fields
	my $id = $hash->{'id'};
	if(!$id or $id eq 'NULL'){
	    #print "FOUND NULL\n";
	    next;
	}

	for(my $i=0; $i<$#fields; $i++){
	    my $value = $hash->{$fields[$i]};
	#    if(not defined $value){ next;} #handle left outer join empties
	    chomp $value;
	    print OUT $value."\t";
	}
	my $value = $hash->{$fields[$#fields]};
	chomp $value;
	print OUT $value."\n";
	
    }
    
    close(OUT);
    
}

# (\@tableCols, $hashArr) = returnDBResultsForTable($dbh, $tableName)
# This will return hash results for given table name, selecting for a particular table distinctly with all manner of possible joins (or those that are for HTC project, not all of ptmscout schema)
# Inputs: $dbh - database handle
#         $tableName - name of table
#         $experiments - ref to array of experiments to include
# Kristen Naegle
# Dec 30, 2009
sub returnDBResultsForTable($$$){
    my($dbh, $tableName, $experiments) = @_;
 
    my $expStr = makeSelectINString($experiments, 0);

   #check for table name error ?
    my $sth = $dbh->prepare("SELECT distinct $tableName.* from MS join experiment on MS.experiment_id=experiment.id join MS_phosphopep on MS.id=MS_phosphopep.MS_id left outer join data on MS.id=data.MS_id join phosphopep on MS_phosphopep.phosphopep_id=phosphopep.id left outer join phosphopep_prediction on phosphopep.id=phosphopep_prediction.phosphopep_id join protein on MS.protein_id=protein.id join domain on protein.id=domain.protein_id join acc on acc.protein_id=protein.id left outer join protein_GO on protein_GO.protein_id=protein.id left outer join GO on protein_GO.GO_id=GO.id where MS.experiment_id IN $expStr");
    $sth->execute();
    my $results = returnMultipleResultsForExSTH($sth);

    # need to get the column names in the proper order
    my $sthCols = $dbh->prepare("SELECT * from $tableName where 1=0;");
    $sthCols->execute();
    my @cols = @{$sth->{NAME}}; 
    $sthCols->finish;

    return (\@cols, $results);


}


# createDBFiles($dbh, $dir);
# Create all db files in a target directory.  
# DB files created for tables: see returnTableOrder for all tables included.
# inputs: $dbh - database handle
#         $dir - target directory where tableName.txt will be written
#         $experimentIds - ref to array of experiments
# Kristen Naegle
# Dec 30, 2009
sub createDBFiles($$$){
    my($dbh, $dir, $experimentIds) = @_;
    
    if($dir !~ m/^\//){
	$dir .= "/";
    }

    my $tables = returnTableOrder();
    foreach my $table (@$tables){
	my ($cols, $results) = returnDBResultsForTable($dbh, $table, $experimentIds);
	my $fileName = $dir.$table.".txt";
	printTableFromHash($fileName, $cols, $results);

    }


}

# @tables = returnTableOrder();
# Returns the appropriate table order for loading for ptmscout database slice necessary for enrichment (i.e. no expression information and no ambiguity)
# Outputs: $table - ref to array of tables in order of load. Meets FK constraints
# Kristen Naegle
# Dec 30, 2009
sub returnTableOrder(){
    my @tables = ('experiment', 'protein', 'MS', 'data', 'phosphopep', 'MS_phosphopep', 'phosphopep_prediction', 'domain', 'acc', 'GO', 'protein_GO');
    return \@tables;
}

# loadDBFiles($dbh, $dir);
#Given a directory and database handle, load all files into appropriate tables
# Inputs: $dbh - database handle
#         $dir - directory location of files, must be labeled tableName.txt, see createDBFiles
# Kristen Naegle
# Dec 31, 2009 
sub loadDBFiles($$){
    my ($dbh, $dir) = @_;
   
    my $tables = returnTableOrder();
    if($dir !~ m/^\//){
	$dir .= "/";
    }    

    for(my $i=0; $i<scalar(@$tables); $i++){

	my $tableName = $tables->[$i];
	print "DEBUG: Table $tableName\n";
	my $fileName = $dir.$tableName.".txt";
	if(!-e $fileName){
	    print "ERROR. Could not find $fileName for load\n";
	    print "Continuing with remaining files\n";
	}
	else{
	    my $sth = $dbh->prepare("LOAD DATA LOCAL INFILE '$fileName' INTO TABLE $tableName FIELDS TERMINATED BY '\t' LINES TERMINATED by '\n'");
	    $sth->execute();
	    $dbh->commit();
	}


    }


}

# insertResultIntoTable($tableName, $hash); 
# Inserts a single result hash into the given table based on field names and values of hash
# Inputs: $tableName - name of table to insert
#         $hash - hash ref where keys are table field names and values are corresponding values



# $hashArr = returnResultsHashFromDBFile($fileName);
# Takes a DBFile such as what's written in createDBFiles and reverses process back into array of hashes
# Inputs: $file
# Outputs: $hashArr - array of hashes with keys equal to the table field and value equal to table value
# Kristen Naegle
# Dec 31, 2009
sub returnResultsHashFromDBFile($){
    my ($fileName) = @_;
    my @results;

    open(F, $fileName) || die "Can't open file $fileName for reading\n";
    #first line is header with order of table fields
    my $header = <F>;
    chomp $header;
    my @header = split('\t', $header);
    while(defined (my $line = <F>)){
	chomp $line;
	my @line = split('\t', $line);
	my %hash;
	for(my $i=0; $i<=$#line; $i++){
	    $hash{$header[$i]} = $line[$i];
	}
	push @results, \%hash;

    }

    close(F);
    return \@results;

}


# exprotSchema($db, $user, $password, $outputFile);
# Calls mysqldump for schema only on a database and writes to the given file
# Inputs: $db - name of database
#         $user - user
#         $password - password
#         $outputFile - destination of mysqldump
# Kristen Naegle
# Dec. 29,2009
# Never mind. Leaving this here for posterity. But the -p flag prompts for a password, so instead I'm going to create a one time schema that will be packaged up each time. 
sub exportSchema($$$$){
    my ($db, $user, $password, $outputFile) = @_;
    
    `mysqldump -u $user -p $password $db --no-data > $outputFile`;

}

# printRefFromPMID($PMID, $outputFile)
# Prints a file in form Field\tValue for fields and values of a reference. Prints an empty file if it was a bad PMID fetch
# Inputs: $PMID - pubmed id
#         $outputFile - destination of results
# Kristen Naegle
# January 9, 2010
sub printRefFromPMID($$){
    my($PMID, $outputFile) = @_;
    

    my $refHash = returnRefFromPMID($PMID);
    open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
    foreach my $key (keys %$refHash){
	#convert fields for experiment table
	my $newKey = $key;
	if($key eq 'authors'){
	    $newKey = 'author';
	}
	elsif($key eq 'title'){
	    $newKey = 'name';
	}
	elsif($key eq 'date'){
	    $newKey = 'pub_date';
	}
	elsif($key eq 'medline_page'){
	    $newKey = 'pages';
	}
	
	print OUT "$newKey\t".$refHash->{$key}."\n";

    }

    close(OUT);

}

1;
