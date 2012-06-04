use strict;
use warnings;
use DBTools::dbIO;
use commonTools;
use DBTools::insertFunctions;
use DBTools::queryTools;
use globalVars;

# getNCBIFile($number)
# Given a part of the ncbi database based on the number, get it and unzip it to /data/knaegle/data/NCBI_DB/File
# Inputs: $number - the number of the refseq to get
# Outputs: writes file to system
# Kristen Naegle
# April 30, 2008
sub getNCBIFile($){
    my ($number) = @_;
    my $file = returnFileNameMV($number);
    my $url = makeNCBIMVurl($number);
    `wget $url`;
    my $path = $globalVars::REFSEQ_PATH;
    `mv $file $path.`;
    `gunzip $path$file`;

}

# $url = makeNCBIMVurl($number)
# Create the url to retreive the vertebrate mammalian file according to the number
# Inputs: $number - number portion of file url
# Outputs: $url - url that can be retrived
# Kristen Naegle
# April 30, 2008
sub makeNCBIMVurl($){
    my ($number) = @_;

    my $url = 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/';
    $url .= returnFileNameMV($number);
    return $url;

}

# $fileName = returnFileNameMV($number) 
# Returns the file name of vertebrate mammalian refseq file according to number
# Inputs: $number - number portion of name field
# Outputs: $fileName - name of file
# Kristen Naegle
# April 30, 2008
sub returnFileNameMV($){
    my ($number) = @_;
    my $file = 'vertebrate_mammalian'.$number.'.protein.faa.gz';
    return $file;

}


# addFileToDB($dbh, $fileName, $date, $CHECK)
# Adds a retrieved refseq file to refseq database. 
# Inputs: $dbh - database handle to refseq database
#         $fileName - name of local file to load
#         $date - Month Year that you are loading this (ideally, corresponds to release date of refseq files
#         $CHECK - BOOLEAN value, set to 0 if this is the first time this has been loaded, if partially loaded set to 1 to prevent redundancies. SERIOUS Slowdown to set to 1
# Outputs: writes file to database
# Kristen Naegle
# April 30, 2008
sub addFileToDB($$$$){
    my ($dbh, $fileName, $date, $CHECK) = @_;
    my $path = $globalVars::REFSEQ_PATH;
    my $file = $path.$fileName;
    $file =~ s/.gz//;
    open(FH_RF, $file) || die "Can't open refseq fasta file $file for reading\n";
    my $line = <FH_RF>;
    print "Testing parser\n";
    my $sth = returnInsertStatementRefSeq($dbh);


    my %args;
    $args{'gi'} = '?';
    my $statement = createSelectString('id', 'refseq', \%args);
    my $sthCheckExist = $dbh->prepare($statement);
    while(defined($line) && $line =~ /^>/){
	my($gi, $xp, $name, $species) = parseRefSeqFasta($line);
	$line = <FH_RF>;
	chomp $line;
	my $sequence = $line;
	while(defined($line = <FH_RF>) && $line !~ /^>/){
	    #$line = <FH_RF>;
	    chomp $line;
	    $sequence .= $line; 
	}
# # 	#print "$gi\t $xp\t $name\t $species\n";
# # 	#$line = <FH_RF>;
# # 	#chomp $line;
# # 	#print $sequence."\n";
# # 	#check for existence
# # #	my %args;
# # #	$args{'gi'} = $gi;
# # 	#my $exist = checkForExistence($dbh, 'id', 'refseq', \%args);
	my $exist = -1;
	if($CHECK){
	    $sthCheckExist->execute($gi);
	    $exist = returnSingleResultOnCol($sthCheckExist, 0);
	}
	if($exist == -1){
	    if($gi && $xp && $name && $species && $sequence){
		insertRefSeqEntry($sth, $gi, $xp, $name, $species, $sequence, $date);
	    }
	}
#	}
	

    }
    close(FH_RF);
    

}

# ($gi, $xp, $name, $species) = parseRefSeqFasta($line)
# Given a fasta header line of refseq file, parse and return relevant items
# Inputs: $line - fasta header line
# Outputs:$gi - gi accession value (includes gi|)
#         $xp - refseq accession
#         $name - protein name
#         $species - two word species name
# Kristen Naegle
# April 30, 2008
sub parseRefSeqFasta($){
    my ($line) = @_;
    chomp $line;
    my ($gi, $xp, $name, $species);
    if($line !~ /^>/){
	print "ERROR: Not a fasta header line in parseRefSeqFasta\n";
    }
    else{
	my @line = split('\|', $line);
	
	$gi = "gi|".$line[1];
	$xp = $line[3];
	my $rest = $line[4];
#	$rest =~ m/\[(.+)\]$|\[.+\]\[(.+)\]$/;
	$rest =~ m/\[([\w\s]+)\]$/;
	$species = lc($1);
	
	$name = $rest;
	$name =~ s/PREDICTED: //;
        $name =~ s/^ //;
	#$name =~ s/\[$species\]$//;
	}
    return ($gi, $xp, $name, $species);
}

# clearRefSeqDB($dbh);
# Given refseq database handle clear the refseq table (does not commit)
# Inputs: $dbh - database handle to refseq database
# Kristen Naegle
# Nov. 24, 2009
sub clearRefSeqDB($){
    my ($dbh) = @_;
    
    my $sth = $dbh->prepare('delete from refseq');
    $sth->execute();
    
}

# $fileNumbersLoaded = addAllFilesToDB($dbh, $CHECK);
# automatically parse all current protein files and add them to the database. Assumes you have clared the database
# Inputs: $dbh - database handle to refseq
#         $CHECK - boolean as to whether to check first for proteins before loading - only set to 1 if you have a partial load
# Outputs: $fileNumbersLoaded - ref. to an array of file numbers loaded from current load
# Kristen Naegle
# Nov. 24, 2009
sub addAllFilesToDB($$){
    my($dbh, $CHECK) = @_;

    my $fileNumbers = returnCurrentRefSeqFileNumbers();
    my $date = returnMonthYear();
    foreach my $number (@$fileNumbers){
	getNCBIFile($number);
	my $fileName = returnFileNameMV($number);
	addFileToDB($dbh, $fileName, $date, $CHECK);
	$dbh->commit();

    }
    return $fileNumbers;
}


# $fileNumberArr = returnCurrentRefSeqFileNumbers();
# Hits the ncbi refseq ftp current release folder and parses it to find all file numbers for appropriate protein material
# Outputs: $fileNumberArr - reference to array of file numbers
# Kristen Naegle
# Nov. 24, 2009
sub returnCurrentRefSeqFileNumbers(){
    
    my $url = 'ftp://ftp.ncbi.nih.gov/refseq/release/vertebrate_mammalian/';
    `wget $url`;
    my $tempFile = "index.html";
    open(TF, $tempFile) || die "Can't open $tempFile for reading\n";
    my @numbers; 
    while(defined(my $line=<TF>)){
	if($line =~ m/vertebrate_mammalian(\d+)\.protein\.faa\.gz/){
	    push @numbers, $1;
	}

    }
    close(TF);
    `rm $tempFile`;
    return \@numbers
}

1;
